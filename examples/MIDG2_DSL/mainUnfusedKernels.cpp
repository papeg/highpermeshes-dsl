// remove all DUNE dependencies from Jakob Schenk's gridIteration.hh implementation
// performance results are no longer hurt
// merge version created from midg_cpp_modified and gridIteration.hh by Ayesha Afzal

#include <array>
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <chrono>
#include <numeric>
#include <map>
#include <iterator>
#include <mutex>

#include <HighPerMeshes.hpp>

#include "data3dN03.hpp" 					                  //!< application-dependent discontinuous Galerkin's cubic order node information
#include "RKCoeff.hpp" 						                  //!< application-dependent Runge-Kutta coefficients

#define PML true
#define SAMPLING true

using namespace HPM;

using Vec3D = HPM::dataType::Vec3D;
using Mat3D = HPM::dataType::Mat3D;

#ifdef PML

class Pml
{
    Vec3D coeff_field;
    Vec3D coeff_aux;
    
    public:
    Vec3D sigma;

    Pml(Vec3D sigma): sigma(sigma)
    {
        coeff_field = {
            sigma.x - sigma.y - sigma.z,
            sigma.y - sigma.z - sigma.x,
            sigma.z - sigma.x - sigma.y
        };
        
        coeff_aux = {
            sigma.x * sigma.x - sigma.x * sigma.y - sigma.x * sigma.z + sigma.y * sigma.z,
            sigma.y * sigma.y - sigma.y * sigma.z - sigma.y * sigma.x + sigma.z * sigma.x,
            sigma.z * sigma.z - sigma.z * sigma.x - sigma.z * sigma.y + sigma.x * sigma.y
        };
    }
    
    void apply(Vec3D& rhs, const Vec3D& field, const Vec3D& aux, Vec3D& rhs_aux, double epsilon)
    {
        rhs.x += coeff_field.x * field.x - aux.x / epsilon;
        rhs.y += coeff_field.y * field.y - aux.y / epsilon;
        rhs.z += coeff_field.z * field.z - aux.z / epsilon;
        
        rhs_aux.x = coeff_aux.x * epsilon * field.x - sigma.x * aux.x;
        rhs_aux.y = coeff_aux.y * epsilon * field.y - sigma.y * aux.y;
        rhs_aux.z = coeff_aux.z * epsilon * field.z - sigma.z * aux.z;
    }
};

#endif

class Material
{
public:
    double epsilon_local;
    double epsilon_neighbour;
    double mu_local;
    double mu_neighbour;
    double impedance_local;
    double impedance_neighbour;
    double impedance;
    double conductance_local;
    double conductance_neighbour;
    double conductance;

    Material(double epsilon_local, double epsilon_neighbour, double mu_local, double mu_neighbour): epsilon_local(epsilon_local), epsilon_neighbour(epsilon_neighbour), mu_local(mu_local), mu_neighbour(mu_neighbour)
    {
        impedance_local = sqrt(mu_local / epsilon_local);
        impedance_neighbour = sqrt(mu_neighbour / mu_local);
        impedance = impedance_local + impedance_neighbour;
        
        conductance_local = 1.0 / impedance_local;
        conductance_neighbour = 1.0 / impedance_neighbour;
        conductance = conductance_local + conductance_neighbour;
    }
};

HPM::dataType::Real SourceFunc(double t, double t0, double f, double w, double A)
{
    return A * sin(2 * M_PI * (t-t0) * f) * exp(-((t-t0) * (t-t0) / (w*w*2)));
}

HPM::dataType::Real SourceFuncDerivative(double t, double t0, double f, double w, double A)
{
    return A * exp(-(t-t0)*(t-t0)/ (2 * w * w)) * (2 * M_PI * f * cos(2 * M_PI * f * (t - t0)) - ((t - t0)/(w*w)) * sin(2 * M_PI * f * (t - t0)));
}

class SourceExcitation
{
public:
    HPM::dataType::Vec3D E;
    HPM::dataType::Vec3D H;

    SourceExcitation(HPM::dataType::Vec3D E, HPM::dataType::Vec3D H): E(E), H(H) {}

    SourceExcitation(): E({0.0, 0.0, 0.0}), H({0, 0, 0}) {}

    SourceExcitation operator-() const
    {
        SourceExcitation negative;
        negative.E = -E;
        negative.H = -H;
        return negative;
    }
};

SourceExcitation TFSF_Source(double t, double dz, double sourceFreq, double sourceWidth, double sourceAmp, HPM::dataType::Vec3D normals, double incoming_angle, double polarisation_angle)
{
    if (t - dz < 0) {
        return SourceExcitation();
    }
    double sourceValue = SourceFunc(t-dz, 5 * sourceWidth, sourceFreq, sourceWidth, sourceAmp);

    double sin_polarisation_angle = sin(polarisation_angle);
    double cos_polarisation_angle = cos(polarisation_angle);
    double sin_incoming_angle = sin(incoming_angle);
    double cos_incoming_angle = cos(incoming_angle);

    HPM::dataType::Vec3D E0 = {
        sourceValue * cos_polarisation_angle * cos_incoming_angle,
        sourceValue * sin_polarisation_angle,
        sourceValue * sin_incoming_angle * cos_polarisation_angle
    };

    HPM::dataType::Vec3D H0 = {
        - sourceValue * sin_polarisation_angle * cos_incoming_angle,
        sourceValue * cos_polarisation_angle,
        - sourceValue * sin_incoming_angle * sin_polarisation_angle
    };

    HPM::dataType::Vec3D FxE = CrossProduct(E0, normals);
    HPM::dataType::Vec3D excitationValueE = CrossProduct(normals, FxE);

    HPM::dataType::Vec3D FxH = CrossProduct(H0, normals);
    HPM::dataType::Vec3D excitationValueH = CrossProduct(normals, FxH);
    
    return SourceExcitation(excitationValueE, excitationValueH);
}

int main(int argc, char** argv)
{
    if (argc != 2) {
        std::cerr << "invalid call! pass config file path as argument" << std::endl;
        exit(1);
    }

    HPM::drts::Runtime hpm { HPM::GetBuffer{} };

    using CoordinateType = HPM::dataType::Coord3D;
    using RealType = HPM::dataType::Real;
    using Vec3D = HPM::dataType::Vec3D;
    using Mat3D = HPM::dataType::Mat3D;
    using TFSF = HPM::mesh::TFSF;

    //////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                         Reading Configuration, DG and Mesh Files                                     //
    //////////////////////////////////////////////////////////////////////////////////////////////////////////
    auto maint1 = std::chrono::high_resolution_clock::now();

    /** \brief read configuration file */
    HPM::auxiliary::ConfigParser CFG(argv[1]);
    const RealType startTime = CFG.GetValue<RealType>("StartTime"); 		//!< get the value of a user-specific starting time
    const RealType finalTime = CFG.GetValue<RealType>("FinalTime"); 		//!< get the value of a user-specific stop time


    std::cout << "simulating from " << startTime << " to " << finalTime << std::endl;

    /** \brief read mesh file */
    const std::string meshFile = CFG.GetValue<std::string>("MeshFile"); 	//!< get the name of a user-specific mesh file
    using Mesh = HPM::mesh::Mesh<CoordinateType, HPM::entity::Simplex>;
    // TODO: double to realtype
    const Mesh mesh = Mesh::template CreateFromFileWithGroups<HPM::auxiliary::GambitMeshFileReader>(meshFile);

    /** \brief read application-dependent discontinuous Galerkin's stuff */
    constexpr std::size_t order = 3;
    using DG = DgNodes<RealType, Vec3D, order>;
    HPM::DG::DgNodesMap<DG, Mesh> DgNodeMap(mesh);

    #ifdef PML
    /** \brief read PML parameters from config file */
    Vec3D sigma =  {CFG.GetValue<RealType>("SigmaX"), CFG.GetValue<RealType>("SigmaY"), CFG.GetValue<RealType>("SigmaZ")};
    for (size_t i = 0; i < 3; i++) {
        if (sigma[i] == -1.0) {
            sigma[i] = 10.0 / (2 * mesh.domain.pmlThickness);
        }
    }
    Pml pml(sigma);

    std::cout << "PML: sigma set to " << pml.sigma << std::endl;
    #endif
    
    // magnetic permeability
    const RealType mu = 1.0;

    auto upwind = CFG.GetValue<RealType>("Upwind");
    std::cout << "upwind paramter set to " << upwind << std::endl;

    // //////////////////////////////////////////////////////////////////////////////////////////////////////////
    // //    All three kernels (Maxwell's Surface Kernel, Maxwell's Volume Kernel, Runge-Kutta kernel)         //
    // //////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    auto AllCells { mesh.GetEntityRange<Mesh::CellDimension>() } ;

    constexpr auto Dofs = ::HPM::dof::MakeDofs<0, 0, 0, DG::numVolNodes, 0>();

#ifdef SAMPLING
    auto sampling = CFG.GetValue<int>("Sampling");
    auto planeXYLocation = CFG.GetValue<int>("SamplingXYPlaneLocation");
    if (planeXYLocation == -1) {
        planeXYLocation = mesh.domain.getCenter().z;
    }
    auto planeXZLocation = CFG.GetValue<int>("SamplingXZPlaneLocation");
    if (planeXZLocation == -1) {
        planeXZLocation = mesh.domain.getCenter().y;
    }
    auto resolution = CFG.GetValue<int>("SamplingResolution");
    auto samplePointsXY = HPM::DG::getSamplePointsXY<DG>(mesh.domain, planeXYLocation, resolution);
    auto samplePointsXZ = HPM::DG::getSamplePointsXZ<DG>(mesh.domain, planeXZLocation, resolution);
    if (sampling) {
        std::cout << "Sampling XY plane at z = " << planeXYLocation << " with resolution " << resolution << "²" << std::endl;
        std::cout << "Sampling XZ plane at y = " << planeXZLocation << " with resolution " << resolution << "²" << std::endl;
    } else {
        std::cout << "Sampling disabled" << std::endl;
    }
    
    auto detector = CFG.GetValue<int>("Detector");
    auto detectorSamplePoints = HPM::DG::getDetectorSamplePoints<DG>(CFG.GetValue<std::string>("DetectorCoordinates"));
    if (detector) {
        std::cout << "Detectors enabled at " << detectorSamplePoints; 
    } else {
        std::cout << "Detectors disabled" << std::endl;
    }
#endif

    /** \brief load initial conditions for fields */
    auto fieldH { hpm.GetBuffer<CoordinateType>(mesh, Dofs) }; 
    auto fieldE { hpm.GetBuffer<CoordinateType>(mesh, Dofs) };

    #ifdef PML
    // TODO: make smaller
    auto fieldQ { hpm.GetBuffer<CoordinateType>(mesh, Dofs) };
    auto fieldP { hpm.GetBuffer<CoordinateType>(mesh, Dofs) };
    #endif

    auto initialCondition = CFG.GetValue<int>("InitialCondition");
    if (initialCondition) {
        std::cout << "writing initial condition" << std::endl;
    } else {
        std::cout << "no initial condition" << std::endl;
    }
    auto sourceFreq = CFG.GetValue<RealType>("SourceFreq");
    auto sourceAmp = CFG.GetValue<RealType>("SourceAmp");
    auto sourceWidth = CFG.GetValue<RealType>("SourceWidth");
    
    Vec3D sourceSphereCenter{CFG.GetValue<RealType>("SourceSphereCenterX"), CFG.GetValue<RealType>("SourceSphereCenterY"), CFG.GetValue<RealType>("SourceSphereCenterZ")};
    if (sourceSphereCenter.x == -1) sourceSphereCenter.x = mesh.domain.getCenter().x;
    if (sourceSphereCenter.y == -1) sourceSphereCenter.y = mesh.domain.getCenter().y;
    if (sourceSphereCenter.z == -1) sourceSphereCenter.z = mesh.domain.getCenter().z;
    auto sourceSphereRadius = CFG.GetValue<RealType>("SourceSphereRadius");

    auto sourceSphere = CFG.GetValue<int>("SourceSphere");
    std::vector<std::tuple<size_t, std::vector<size_t>>> sourceNodes;
    
    HPM::SequentialDispatcher body;

    if (initialCondition || sourceSphere 
#ifdef SAMPLING
    || sampling || detector
#endif 
    ) {
        body.Execute(
            HPM::ForEachEntity(
            AllCells,
            std::tuple(Write(Cell(fieldE))),
            [&] (const auto& cell, auto &&, auto lvs)
            {
                if (initialCondition) {
                    HPM::ForEach(DG::numVolNodes, [&](const auto& n) {
                        auto& fieldE = std::get<0>(lvs);
                        const auto& nodeCoords = DG::LocalToGlobal(DG::referenceCoords[n], cell.GetTopology().GetNodes());
                        fieldE[n].y = sin(M_PI * nodeCoords.x) * sin(M_PI * nodeCoords.z); 	//!< initial conditions for y component of electric field
                    });
                }
                #ifdef SAMPLING
                if (sampling) {
                    // fill sample points
                    samplePointsXY.checkCellForSamplePoints(cell.GetTopology().GetIndex(), cell.GetTopology().GetVertices());
                    samplePointsXZ.checkCellForSamplePoints(cell.GetTopology().GetIndex(), cell.GetTopology().GetVertices());
                }
                if (detector) {
                    detectorSamplePoints.checkCellForSamplePoints(cell.GetTopology().GetIndex(), cell.GetTopology().GetVertices());
                }
                #endif
                if (sourceSphere) {
                    std::vector<size_t> volNodes;
                    HPM::ForEach(DG::numVolNodes, [&](const std::size_t n) {
                        const auto& nodeCoords = DG::LocalToGlobal(DG::referenceCoords[n], cell.GetTopology().GetNodes()); 		                         //!< reference to global nodal coordinates
                        if (mesh.domain.CoordinateIsInRadiusOfOther(nodeCoords, sourceSphereRadius, sourceSphereCenter)) {
                            volNodes.push_back(n);
                        }
                    });
                    if (volNodes.size() > 0) {
                        sourceNodes.push_back(std::make_tuple(cell.GetTopology().GetIndex(), volNodes));
                    }
                }
            })
        );
    }

    if (sourceSphere) {
        std::cout << "sourcing sphere at " << sourceSphereCenter << " with radius " << sourceSphereRadius << std::endl;
        for (size_t i = 0; i < sourceNodes.size(); i++) {
            std::cout << "including " << std::get<1>(sourceNodes[i]).size() << " node(s) in element " << std::get<0>(sourceNodes[i]) << std::endl;
        }
        std::cout << "frequency " << sourceFreq << ", width " << sourceWidth << ",and amplitude " << sourceAmp << std::endl;
    } else {
        std::cout << "no pulsing" << std::endl;
    }


    /** \brief create storage for intermediate fields*/
    auto resH { hpm.GetBuffer<CoordinateType>(mesh, Dofs) };
    auto resE { hpm.GetBuffer<CoordinateType>(mesh, Dofs) };
    auto rhsH { hpm.GetBuffer<CoordinateType>(mesh, Dofs) };
    auto rhsE { hpm.GetBuffer<CoordinateType>(mesh, Dofs) };

    #ifdef PML
    /** \brief create storage for intermediate auxiliary fields*/
    auto resQ { hpm.GetBuffer<CoordinateType>(mesh, Dofs) };
    auto resP { hpm.GetBuffer<CoordinateType>(mesh, Dofs) };
    auto rhsQ { hpm.GetBuffer<CoordinateType>(mesh, Dofs) };
    auto rhsP { hpm.GetBuffer<CoordinateType>(mesh, Dofs) };
    #endif

    /** \brief determine time step size (polynomial order-based and algorithmic-specific) */
    RealType timeStep = 1.0e6;

    HPM::ForEachEntity(
        AllCells, 
        [&] (const auto& cell)
            {
                const RealType face_normal_scaling_factor = 2.0 / cell.GetGeometry().GetAbsJacobianDeterminant();

                HPM::ForEachSubEntity(cell, [&](const auto& face){
                    timeStep = std::min(timeStep, 1.0 / (face.GetGeometry().GetNormal() * face_normal_scaling_factor).Norm());
                }); 
            });

    timeStep = finalTime / floor(finalTime * (order + 1) * (order + 1)/(.5 * timeStep) );
    std::cout << "time step: " << timeStep << std::endl;

    // check if timestep is reasonable
    assert(timeStep > 5.0e-4);

    // write out source function and derivative
    std::ofstream FuncValuesFile;
    FuncValuesFile.open("func.values", std::fstream::trunc);
    FuncValuesFile << "t,f,d\n";
    double i = 0.0;
    std::cout.precision(17);
    while(true) {
        double t = i * timeStep; // s
        double f = SourceFunc(t, 5 * sourceWidth, sourceFreq, sourceWidth, sourceAmp);
        double d = SourceFuncDerivative(t, 5 * sourceWidth, sourceFreq, sourceWidth, sourceAmp);
        FuncValuesFile << std::fixed << t << "," << f << "," << d << '\n';
        i += 1.0;
        if (t > finalTime) {
            std::cout << "ended at t = " << t << ", iter = " << i << std::endl;
            break;
        }
    }
    FuncValuesFile.close();
    
    auto sourcePlaneXY = CFG.GetValue<int>("SourcePlaneXY");
    if (sourcePlaneXY) {
        std::cout << "plane wave with amplitude " << sourceAmp << ", frequency " << sourceFreq << " and width " << sourceWidth << std::endl;
    } else {
        std::cout << "no plane wave" << std::endl;
    }
    auto materialEnabled = CFG.GetValue<int>("MaterialEnabled");
    if (materialEnabled) {
        std::cout << "material parameter enabled" << std::endl;
    } else {
        std::cout << "ignoring material parameter" << std::endl;
    }

    auto maint2 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> setup_duration = maint2 - maint1;
    std::cout << "Setup time in seconds: " << setup_duration.count() << std::endl;
    double aggregate_time1 = 0.0, aggregate_time2 = 0.0, aggregate_time3 = 0.0;

    {
        auto t1 = std::chrono::high_resolution_clock::now();
        
        /** \brief Maxwell's surface kernel */
        auto surfaceKernelLoop = HPM::ForEachIncidence<2>(
            AllCells,
            std::tuple(
                Read(ContainingMeshElement(fieldH)),
                Read(ContainingMeshElement(fieldE)),
                Read(NeighboringMeshElementOrSelf(fieldH)),
                Read(NeighboringMeshElementOrSelf(fieldE)),
                Write(ContainingMeshElement(rhsH)),
                Write(ContainingMeshElement(rhsE))
                #ifdef PML
                ,Read(ContainingMeshElement(fieldQ)),
                Read(ContainingMeshElement(fieldP)),
                Read(NeighboringMeshElementOrSelf(fieldQ)),
                Read(NeighboringMeshElementOrSelf(fieldP)),
                Write(ContainingMeshElement(rhsQ)),
                Write(ContainingMeshElement(rhsP))
                #endif
            ),
            [&](const auto &element, const auto &face, const auto &iter, auto &lvs) {
                const std::size_t face_index = face.GetTopology().GetLocalIndex();
                const RealType face_normal_scaling_factor = 2.0 / element.GetGeometry().GetAbsJacobianDeterminant();

                const Vec3D &face_normal = face.GetGeometry().GetNormal() * face_normal_scaling_factor; //!< get all normal coordinates for each face of an element
                const RealType Edg = face_normal.Norm() * 0.5;                                          //!< get edge length for each face
                const Vec3D &face_unit_normal = face.GetGeometry().GetUnitNormal();
                const auto &localMap{DgNodeMap.Get(element, face)};
                

                RealType epsilon_local = element.GetTopology().GetMaterial();
                const auto neighbour = face.GetTopology().GetNeighboringCell().GetTopology();
                RealType epsilon_neighbour = neighbour.GetMaterial();

                TFSF tfsf = mesh.tfsf[element.GetTopology().GetIndex()][face_index];

                bool element_is_pml = mesh.pml[element.GetTopology().GetIndex()];

                Material material = Material(epsilon_local, epsilon_neighbour, mu, mu);

                HPM::ForEach(DG::NumSurfaceNodes, [&](const std::size_t m) {
                    SourceExcitation excitationValue = SourceExcitation();
                    if (sourcePlaneXY) {
                        const auto& coords = DG::LocalToGlobal(DG::referenceCoords[localMap[m][0]], element.GetTopology().GetNodes());
                        const auto t = timeStep * ((iter / 5) + rk4c<RealType>[iter % 5]);
                        if (tfsf == TFSF::FromScatteredToTotal) {
                            RealType dz = coords.z - mesh.domain.tfsfBegin.z;
                            excitationValue = TFSF_Source(t, dz, sourceFreq, sourceWidth, sourceAmp, face_unit_normal, 0, 0);
                        } else if (tfsf == TFSF::FromTotalToScattered) {
                            RealType dz = coords.z - mesh.domain.tfsfBegin.z;
                            excitationValue = - TFSF_Source(t, dz, sourceFreq, sourceWidth, sourceAmp, face_unit_normal, 0, 0);
                        }
                    }
                    const auto &fieldH = std::get<0>(lvs);
                    const auto &fieldE = std::get<1>(lvs);

                    auto &NeighboringFieldH = std::get<2>(lvs);
                    auto &NeighboringFieldE = std::get<3>(lvs);

                    const Vec3D &dH = Edg * HPM::DG::Delta(fieldH, NeighboringFieldH, m, localMap) + excitationValue.H; //!< fields differences
                    const Vec3D &dE = Edg * HPM::DG::DirectionalDelta(fieldE, NeighboringFieldE, face, m, localMap) + excitationValue.E;

                    const Vec3D &flux_H = (upwind*((dH) + ((dH)*face_unit_normal) * face_unit_normal) - material.impedance_local * CrossProduct(face_unit_normal, (dE))) / material.impedance * 2.0; //!< fields fluxes
                    const Vec3D &flux_E = (upwind*((dE) - ((dE)*face_unit_normal) * face_unit_normal) + material.impedance_local * CrossProduct(face_unit_normal, (dH))) / material.conductance * 2.0;

                    auto &rhsH = std::get<4>(lvs);
                    auto &rhsE = std::get<5>(lvs);

                    HPM::ForEach(DG::numVolNodes, [&](const std::size_t n) {
                        rhsH[n] += DG::LIFT[face_index][m][n] * flux_H;
                        rhsE[n] += DG::LIFT[face_index][m][n] * flux_E;
                    });

                    #ifdef PML
                    if (element_is_pml) {
                        const auto &fieldQ = std::get<6>(lvs);
                        const auto &fieldP = std::get<7>(lvs);

                        auto &NeighboringFieldQ = std::get<8>(lvs);
                        auto &NeighboringFieldP = std::get<9>(lvs);
                        
                        const Vec3D &dQ = Edg * HPM::DG::Delta(fieldQ, NeighboringFieldQ, m, localMap) + excitationValue.H; //!< fields differences
                        const Vec3D &dP = Edg * HPM::DG::DirectionalDelta(fieldP, NeighboringFieldP, face, m, localMap) + excitationValue.E;
                        
                        const Vec3D &flux_Q = (upwind*((dQ) - ((dQ)*face_unit_normal) * face_unit_normal) - material.impedance_local * CrossProduct(face_unit_normal, (dQ))) / material.impedance * 2.0;
                        const Vec3D &flux_P = (upwind*((dP) - ((dP)*face_unit_normal) * face_unit_normal) + material.impedance_local * CrossProduct(face_unit_normal, (dP))) / material.conductance * 2.0; //!< fields fluxes

                        auto &rhsQ = std::get<10>(lvs);
                        auto &rhsP = std::get<11>(lvs);
                        
                        HPM::ForEach(DG::numVolNodes, [&](const std::size_t n) {
                            rhsQ[n] += DG::LIFT[face_index][m][n] * flux_Q;
                            rhsP[n] += DG::LIFT[face_index][m][n] * flux_P;
                        });
                    }
                    #endif
                });
            },
            HPM::internal::OpenMP_ForEachIncidence<Mesh::CellDimension, 2>{}
        );

        auto t2 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> duration = t2 - t1;
        aggregate_time1 += duration.count();
        t1 = std::chrono::high_resolution_clock::now();

        /** \brief Maxwell's volume kernel */
        auto volumeKernelLoop = HPM::ForEachEntity(
            AllCells,
            std::tuple(
                Read(Cell(fieldH)),
                Read(Cell(fieldE)),
                Cell(rhsH),
                Cell(rhsE)
                #ifdef PML
                ,Read(Cell(fieldQ)),
                Read(Cell(fieldP)),
                Cell(rhsQ),
                Cell(rhsP)
                #endif
            ),
            [&](const auto &element, const auto &iter, auto &lvs) {
                const Mat3D &D = element.GetGeometry().GetInverseJacobian() * 2.0;
                RealType epsilon = element.GetTopology().GetMaterial();
                bool element_is_pml = mesh.pml[element.GetTopology().GetIndex()];
                if (!materialEnabled) {
                    epsilon = 1.0;
                }

                HPM::ForEach(DG::numVolNodes, [&](const std::size_t n) {
                    Mat3D derivative_E, derivative_H; //!< derivative of fields w.r.t reference oordinates
                    Mat3D derivative_P, derivative_Q; //!< derivative of fields w.r.t reference coordinates

                    const auto &fieldH = std::get<0>(lvs);
                    const auto &fieldE = std::get<1>(lvs);
                    #ifdef PML
                    const auto &fieldQ = std::get<4>(lvs);
                    const auto &fieldP = std::get<5>(lvs);
                    #endif

                    HPM::ForEach(DG::numVolNodes, [&](const std::size_t m) {
                        derivative_H += DyadicProduct(DG::derivative[n][m], fieldH[m]);
                        derivative_E += DyadicProduct(DG::derivative[n][m], fieldE[m]);
                        #ifdef PML
                        if (element_is_pml) {
                            derivative_Q += DyadicProduct(DG::derivative[n][m], fieldQ[m]);
                            derivative_P += DyadicProduct(DG::derivative[n][m], fieldP[m]);
                        }
                        #endif
                    });

                    auto &rhsH = std::get<2>(lvs);
                    auto &rhsE = std::get<3>(lvs);
                    #ifdef PML
                    auto &rhsQ = std::get<6>(lvs);
                    auto &rhsP = std::get<7>(lvs);
                    #endif

                    rhsH[n] += -Curl(D, derivative_E) / mu; //!< first half of right-hand-side of fields
                    rhsE[n] += Curl(D, derivative_H) / epsilon;

                    #ifdef PML
                    if (element_is_pml) {
                        rhsQ[n] += -Curl(D, derivative_P) / mu;
                        rhsP[n] += Curl(D, derivative_Q) / epsilon;
                        
                        pml.apply(rhsH[n], fieldH[n], fieldQ[n], rhsQ[n], epsilon);
                        pml.apply(rhsE[n], fieldE[n], fieldP[n], rhsP[n], epsilon);
                    }
                    #endif
                });
                if (sourceSphere) {
                    // dipole source
                    auto &rhsE = std::get<3>(lvs);

                    for (const auto& node: sourceNodes) {
                        if (element.GetTopology().GetIndex() == std::get<0>(node)) {
                            const auto t = timeStep * ((iter / 5) + rk4c<RealType>[iter % 5]);
                            for (const auto& n: std::get<1>(node)) {
                                rhsE[n].x += SourceFunc(t, 5 * sourceWidth, sourceFreq, sourceWidth, sourceAmp);
                            }
                        }
                    }
                }
            },
            HPM::internal::OpenMP_ForEachEntity<Mesh::CellDimension>{} 
        );

        t2 = std::chrono::high_resolution_clock::now();
        duration = t2 - t1;
        aggregate_time2 += duration.count();
        t1 = std::chrono::high_resolution_clock::now();

        /** \brief Runge-Kutta integrtion kernel */
        auto rungeKuttaLoop =
            HPM::ForEachEntity(
                AllCells,
                std::tuple(
                    Write(Cell(fieldH)),
                    Write(Cell(fieldE)),
                    Cell(rhsH),
                    Cell(rhsE),
                    Cell(resH),
                    Cell(resE)
                    #ifdef PML
                    ,Write(Cell(fieldQ)),
                    Write(Cell(fieldP)),
                    Cell(rhsQ),
                    Cell(rhsP),
                    Cell(resQ),
                    Cell(resP)
                    #endif
                ),
                [&](const auto &element, const auto &iter, auto &lvs) {
                    
                    bool element_is_pml = mesh.pml[element.GetTopology().GetIndex()];

                    const auto &RKstage = RungeKuttaCoeff<RealType>::rk4[iter % 5];

                    auto &fieldH = std::get<0>(lvs);
                    auto &fieldE = std::get<1>(lvs);
                    auto &rhsH = std::get<2>(lvs);
                    auto &rhsE = std::get<3>(lvs);
                    auto &resH = std::get<4>(lvs);
                    auto &resE = std::get<5>(lvs);
                    #ifdef PML
                    if (element_is_pml) {
                        auto &fieldQ = std::get<6>(lvs);
                        auto &fieldP = std::get<7>(lvs);
                        auto &rhsQ = std::get<8>(lvs);
                        auto &rhsP = std::get<9>(lvs);
                        auto &resQ = std::get<10>(lvs);
                        auto &resP = std::get<11>(lvs);

                        HPM::ForEach(DG::numVolNodes, [&](const std::size_t n) {
                            resQ[n] = RKstage[0] * resQ[n] + timeStep * rhsQ[n];
                            resP[n] = RKstage[0] * resP[n] + timeStep * rhsP[n]; //!< residual fields

                            fieldQ[n] += RKstage[1] * resQ[n];
                            fieldP[n] += RKstage[1] * resP[n]; //!< updated fields

                            assign_to_entries(rhsQ[n], 0.0);
                            assign_to_entries(rhsP[n], 0.0); //TODO
                        });
                    }
                    #endif

                    HPM::ForEach(DG::numVolNodes, [&](const std::size_t n) {
                        resH[n] = RKstage[0] * resH[n] + timeStep * rhsH[n]; //!< residual fields
                        resE[n] = RKstage[0] * resE[n] + timeStep * rhsE[n];

                        fieldH[n] += RKstage[1] * resH[n]; //!< updated fields
                        fieldE[n] += RKstage[1] * resE[n];

                        assign_to_entries(rhsH[n], 0.0); //TODO
                        assign_to_entries(rhsE[n], 0.0);
                    });
                },
                HPM::internal::OpenMP_ForEachEntity<Mesh::CellDimension>{} 
            );


        #ifdef SAMPLING
        std::ofstream XYValuesFile;
        std::ofstream XZValuesFile;
        std::ofstream DetectorValuesFile;
        if (sampling) {
            std::string XYValuesFilePath = CFG.GetValue<std::string>("SamplingXYValuesFile");
            XYValuesFile.open(XYValuesFilePath.c_str());
            XYValuesFile << "step,i,eps,x,y,z,Ex,Hy\n";

            std::string XZValuesFilePath = CFG.GetValue<std::string>("SamplingXZValuesFile");
            XZValuesFile.open(XZValuesFilePath.c_str());
            XZValuesFile << "step,i,eps,x,y,z,Ex,Hy\n";
        }
        if (detector) {
            DetectorValuesFile.open(CFG.GetValue<std::string>("DetectorValuesFile").c_str());
            DetectorValuesFile << "step,t,i,x,y,z,Ex,Ey,Ez,Hx,Hy,Hz\n";
        }
        std::mutex XYValuesFileLock;
        std::mutex XZValuesFileLock;
        std::mutex DetectorValuesFileLock;

        //using namespace HPM::DG::util;

        auto writeValuesLoop =
            HPM::ForEachEntity(
                AllCells,
                std::tuple(
                    Read(Cell(fieldE)),
                    Read(Cell(fieldH))
                ),
                [&](const auto &cell, const auto &iter, auto &lvs) {
                    if (iter % 10 == 0) {
                        if (sampling && iter % 200 == 0) {
                            auto &fieldE = std::get<0>(lvs);
                            const auto cellIndex = cell.GetTopology().GetIndex();
                            auto pointIterator = samplePointsXY.indexMap.equal_range(cellIndex);

                            for (auto it = pointIterator.first; it != pointIterator.second; it++) {
                                const auto pointIndex = it->second;
                                const auto point = samplePointsXY.points[pointIndex];
                                const auto weights = samplePointsXY.interpolationWeights[pointIndex];

                                Vec3D weightedE = {0.0, 0.0, 0.0};
                                Vec3D weightedH = {0.0, 0.0, 0.0};

                                HPM::ForEach(
                                    DG::numVolNodes,
                                    [&](const size_t n) {
                                        weightedE += fieldE[n] * weights[n];
                                        weightedH += fieldE[n] * weights[n];
                                    });

                                std::lock_guard<std::mutex> lock(XYValuesFileLock);
                                XYValuesFile << iter << "," << pointIndex << "," << cell.GetTopology().GetMaterial() << "," << point.x << "," << point.y << "," << point.z << "," << weightedE.x << "," << weightedH.y << "\n";
                            }
                            pointIterator = samplePointsXZ.indexMap.equal_range(cellIndex);
                            for (auto it = pointIterator.first; it != pointIterator.second; it++) {
                                const auto pointIndex = it->second;
                                const auto point = samplePointsXZ.points[pointIndex];
                                const auto weights = samplePointsXZ.interpolationWeights[pointIndex];

                                Vec3D weightedE = {0.0, 0.0, 0.0};
                                Vec3D weightedH = {0.0, 0.0, 0.0};

                                HPM::ForEach(
                                    DG::numVolNodes,
                                    [&](const size_t n) {
                                        weightedE += fieldE[n] * weights[n];
                                        weightedH += fieldE[n] * weights[n];
                                    });

                                std::lock_guard<std::mutex> lock(XZValuesFileLock);
                                XZValuesFile << iter << "," << pointIndex << "," << cell.GetTopology().GetMaterial() << "," << point.x << "," << point.y << "," << point.z << "," << weightedE.x << "," << weightedH.y << "\n";
                            }
                        }
                        if (detector) {
                            auto &fieldE = std::get<0>(lvs);
                            auto &fieldH = std::get<1>(lvs);
                            const auto cellIndex = cell.GetTopology().GetIndex();
                            auto pointIterator = detectorSamplePoints.indexMap.equal_range(cellIndex);
                            for (auto it = pointIterator.first; it != pointIterator.second; it++) {
                                const auto pointIndex = it->second;
                                const auto point = detectorSamplePoints.points[pointIndex];
                                const auto weights = detectorSamplePoints.interpolationWeights[pointIndex];

                                Vec3D weightedE = {0.0, 0.0, 0.0};
                                Vec3D weightedH = {0.0, 0.0, 0.0};

                                HPM::ForEach(
                                    DG::numVolNodes,
                                    [&](const size_t n) {
                                        weightedE += fieldE[n] * weights[n];
                                        weightedH += fieldE[n] * weights[n];
                                    });

                                const auto t = timeStep * (iter / 5);
                                std::lock_guard<std::mutex> lock(DetectorValuesFileLock);
                                DetectorValuesFile << iter << "," << t << "," << pointIndex << "," << point.x << "," << point.y << "," << point.z << "," << weightedE.x << "," << weightedE.y << "," << weightedE.z << "," << weightedH.x << "," << weightedH.y << "," << weightedH.z << "\n";
                            }
                        }
                    }
                },
                HPM::internal::OpenMP_ForEachEntity<Mesh::CellDimension>{}
            );

        body.Execute( HPM::iterator::Range<size_t> { static_cast<std::size_t>(((finalTime - startTime) / timeStep) * 5) },
                    surfaceKernelLoop, volumeKernelLoop, rungeKuttaLoop, writeValuesLoop);
        
        if (sampling) {
            XYValuesFile.close();
            XZValuesFile.close();
        }
        if (detector) {     
            DetectorValuesFile.close();
        }
        #else
        body.Execute( HPM::iterator::Range<size_t> { static_cast<std::size_t>(((finalTime - startTime) / timeStep) * 5) },
                    surfaceKernelLoop, volumeKernelLoop, rungeKuttaLoop);
        #endif


        t2 = std::chrono::high_resolution_clock::now();
        duration = t2 - t1;
        aggregate_time3 += duration.count();
    }
    std::cout << "Aggregate execution time for Surface kernel       = " << aggregate_time1 * 1000 << " ms" << std::endl;
    std::cout << "Aggregate execution time for Volume kernel        = " << aggregate_time2 * 1000 << " ms" << std::endl;
    std::cout << "Aggregate execution time for RK kernel            = " << aggregate_time3 * 1000 << " ms" << std::endl;
    std::cout << "Aggregate all kernel execution time               = " << (aggregate_time1 + aggregate_time2 + aggregate_time3) * 1000 << " ms" << std::endl;
    std::cout << "Individual Execution time of Surface kernel       = " << (aggregate_time1 * 1000) / (finalTime / timeStep * 5) << " ms" << std::endl;
    std::cout << "Individual Execution time of Volume kernel        = " << (aggregate_time2 * 1000) / (finalTime / timeStep * 5) << " ms" << std::endl;
    std::cout << "Individual Execution time of RK kernel            = " << (aggregate_time3 * 1000) / (finalTime / timeStep * 5) << " ms" << std::endl;
    std::cout << "Individual all kernel execution time              = " << ((aggregate_time1 + aggregate_time2 + aggregate_time3) * 1000) / (finalTime / timeStep * 5) << " ms" << std::endl;

    /** \brief find maximum & minimum values for Ey*/
    double maxErrorEy = 0;
    double minEy = std::numeric_limits<RealType>::max();
    double maxEy = std::numeric_limits<RealType>::lowest();

    body.Execute(
        HPM::ForEachEntity(
        AllCells,
        std::tuple( Read(Cell(fieldE)) ),
        [&](const auto& element, auto&&, auto lvs)
        {    
            const auto& fieldE = std::get<0>(lvs);

            HPM::ForEach(DG::numVolNodes, [&] (const std::size_t n)
            {
                const auto& nodeCoords = DG::LocalToGlobal(DG::referenceCoords[n], element.GetTopology().GetNodes()); 		                         //!< reference to global nodal coordinates
                const RealType exactEy = sin(M_PI * nodeCoords.x) * sin(M_PI * nodeCoords.z) * cos(sqrt(2.) * M_PI * finalTime); 	                            //!< exact analytical electrical field value in y direction
                if (initialCondition) {
                    maxErrorEy = std::max(maxErrorEy, std::abs(exactEy - fieldE[n].y));   //!< maximum error in electrical field value in y direction
                }
                minEy = std::min(minEy, fieldE[n].y);                                 //!< minimum electric field value in y direction
                maxEy = std::max(maxEy, fieldE[n].y);  
            });
        }
        )
    );

    if (initialCondition) {
        std::cout << "\nt=" << finalTime
            << " Ey in [ " << minEy
            << ", " << maxEy
            << " ] with max nodal error " << maxErrorEy
            << std::endl;
    } else {
        std::cout << "\nt=" << finalTime
            << " Ey in [ " << minEy
            << ", " << maxEy
            << " ]" << std::endl;
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                         Shutdown of the runtime system                                               //
    //////////////////////////////////////////////////////////////////////////////////////////////////////////

    return EXIT_SUCCESS;
}
