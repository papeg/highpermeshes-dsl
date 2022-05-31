// Copyright (c) 2017-2020
//
// Distributed under the MIT Software License
// (See accompanying file LICENSE)

#ifndef MISC_DG_HPP
#define MISC_DG_HPP

#include <HighPerMeshes/dsl/meshes/Mesh.hpp>
#include <string>


namespace HPM::DG
{
    template <std::size_t NumSurfaceNodes>
    using SurfaceMap = std::array<std::array<std::size_t, 2>, NumSurfaceNodes>;

    template <typename BufferT_1, typename BufferT_2, typename T>
    auto Delta(const BufferT_1& buffer, const BufferT_2& neighboring_buffer, const int index, const T& face_node_to_cell_mapping)
    {
        return neighboring_buffer[face_node_to_cell_mapping[index][1]] - buffer[face_node_to_cell_mapping[index][0]];
    }

    template <typename FaceT>
    auto Direction(const FaceT& face) -> ::HPM::dataType::Real
    {
        return (face.GetTopology().HasNeighboringCell() ? 1.0 : -1.0);
    }

    template <typename BufferT_1, typename BufferT_2, typename FaceT, typename T>
    auto DirectionalDelta(const BufferT_1& buffer, const BufferT_2& neighboring_buffer, const FaceT& face, const int index, const T& face_node_to_cell_mapping)
    {
        return Direction(face) * neighboring_buffer[face_node_to_cell_mapping[index][1]] - buffer[face_node_to_cell_mapping[index][0]];
    }

    // \todo { Not sure if threshold is the right name - Stefan G. 23.07.2019 }
    template <typename DgInfo, typename EntityT, typename FaceT>
    static auto ComputeForOneFace(const EntityT& element, const FaceT& face, double threshold = 1.0E-4)
    {
        const std::size_t face_index = face.GetTopology().GetLocalIndex();
        const auto& element_nodes = element.GetTopology().GetNodes();
        const auto& neighboring_element_nodes = face.GetTopology().GetNeighboringCell().GetTopology().GetNodes();
        const std::size_t neighboring_face_index = face.GetTopology().GetLocalIndexOfNeighboringFace();
        SurfaceMap<DgInfo::NumSurfaceNodes> result;

        for (std::size_t n1 = 0; n1 < DgInfo::NumSurfaceNodes; ++n1)
        {
            result[n1][0] = DgInfo::localMask[face_index][n1];
            
            // Go through all degrees of freedom for a surface of the neighboring face and find the one that has a distance below a certain threshhold. If it can't be found it throws an exception.
            for (std::size_t n2 = 0; n2 < DgInfo::NumSurfaceNodes; ++n2)
            {
                // Find normalized distance between these nodes and check if it is (almost) zero
                const auto& d = DgInfo::LocalToGlobal(DgInfo::referenceCoords[DgInfo::localMask[face_index][n1]], element_nodes) -
                                DgInfo::LocalToGlobal(DgInfo::referenceCoords[DgInfo::localMask[neighboring_face_index][n2]], neighboring_element_nodes);

                if (d.Norm() < threshold)
                {
                    result[n1][1] = DgInfo::localMask[neighboring_face_index][n2];
                    break;
                }

                if (n2 == DgInfo::NumSurfaceNodes - 1)
                {
                    throw std::runtime_error("Couldn't find matching neighbor node!");
                }
            }
        }

        return result;
    }

    template <typename DgInfo, typename MeshT>
    struct DgNodesMap
    {
        // \todo { 4 seems like a magic number to me and should probably be inferred by the Mesh topology - Stefan G. 23.07.2019 }
        // \todo { What is the meaning behind this map for the last index? Map[element index][local face index][DOF of face][??] - Stefan G. 23.07.2019 }
        using Map = std::vector<std::array<SurfaceMap<DgInfo::NumSurfaceNodes>, 4>>;

        DgNodesMap(const MeshT& mesh) : mesh(mesh), map(mesh.GetNumEntities())
        {
            for (auto const& element : mesh.GetEntities())
            {
                const std::size_t element_index = element.GetTopology().GetIndex();
                for (auto const& face : element.GetTopology().GetSubEntities())
                {
                    const std::size_t face_index = face.GetTopology().GetLocalIndex();
                    map[element_index][face_index] = ComputeForOneFace<DgInfo>(element, face);
                }
            }
        }

        template <typename FaceT>
        auto Get(const typename MeshT::CellT& element, const FaceT& face) const -> const SurfaceMap<DgInfo::NumSurfaceNodes>&
        {
            static_assert(std::is_same_v<typename FaceT::ParentEntityT, typename MeshT::CellT>, "face must be of face type in the given mesh");

            return map[element.GetTopology().GetIndex()][face.GetTopology().GetLocalIndex()];
        }

        const MeshT& mesh;
        Map map;
    };
    
    using Vec3D = HPM::dataType::Vec3D;

    template <typename DgInfo>
    auto calcBarycentricVector(HPM::dataType::Vec3D samplePoint, std::array<HPM::dataType::Vec3D, 4> tetVertices)
    {
        using Vec3D = HPM::dataType::Vec3D;

        Vec3D coordDifference = samplePoint - tetVertices[3];
        HPM::dataType::Mat3D matrixT;
        int numOfOtherVertices = 3;
        for (int row = 0; row < numOfOtherVertices; ++row)
        {
            matrixT[row] = tetVertices[row] - tetVertices[3];
        }

        double determinantT = matrixT.Determinant();

        auto transposedCoFactorMatrixT = matrixT.Cofactor();

        for (int row = 0; row < transposedCoFactorMatrixT.NumColumns(); row++)
            std::transform(transposedCoFactorMatrixT[row].data.begin(), transposedCoFactorMatrixT[row].data.end(), transposedCoFactorMatrixT[row].data.begin(), [&determinantT](auto &c) { return c / determinantT; });

        std::vector<double> barycentricVector{0.0, 0.0, 0.0, 0.0};

        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                barycentricVector[i] += transposedCoFactorMatrixT[i][j] * coordDifference[j];
            }
        }
        barycentricVector[3] = 1 - barycentricVector[0] - barycentricVector[1] - barycentricVector[2];

        return barycentricVector;
    }

    template <typename DgInfo>
    bool isContainedIn(const std::vector<double> &barycentricVector)
    {
        constexpr double tolerance = -1.0e-5;
        std::vector<double>::const_iterator minValue = std::min_element(barycentricVector.begin(), barycentricVector.end());

        return (*minValue) > tolerance;
    }

    template <typename DgInfo>
    double calcJacobiPolynom(double x, double alpha, double beta, size_t jacobiOrder)
    {
        double alphaOld = 0.0, alphaNew = 0.0, betaNew = 0.0, h1 = 0.0;
        double gamma0 = 0.0, gamma1 = 0.0;
        double alphaPlusBeta = alpha+beta, alphaPlusBetaPlusOne = alpha+beta+1.0, alphaPlusOne = alpha+1.0, betaPlusOne = beta+1.0;

        std::vector<double> PL;
        double prow;
        double xMinusBetaNew;

        gamma0 = std::pow(2.0,alphaPlusBetaPlusOne)/(alphaPlusBetaPlusOne)*HPM::math::Gamma(alphaPlusOne)*HPM::math::Gamma(betaPlusOne)/HPM::math::Gamma(alphaPlusBetaPlusOne);

        if (0 == jacobiOrder) {
            return 1.0/sqrt(gamma0);

        } else {
            PL.push_back(1.0/std::sqrt(gamma0));
        }

        gamma1 = (alphaPlusOne)*(betaPlusOne)/(alphaPlusBeta+3.0)*gamma0;

        prow = ((alphaPlusBeta+2.0)*x/2.0 + (alpha-beta)/2.0) / sqrt(gamma1);

        if (1 == jacobiOrder) {
            return prow;
        } else {
            PL.push_back(prow);
        }

        alphaOld = 2.0/(2.0+alphaPlusBeta)*std::sqrt((alphaPlusOne)*(betaPlusOne)/(alphaPlusBeta+3.0));

        for (size_t i = 1; i <= jacobiOrder - 1; ++i) {  //i = 1, 2
            h1 = 2.0*i+alphaPlusBeta;
            alphaNew = 2.0/(h1+2.0)*std::sqrt((i+1)*(i+alphaPlusBetaPlusOne)*(i+alphaPlusOne)*(i+betaPlusOne)/(h1+1.0)/(h1+3.0));
            betaNew = - ((alpha*alpha)-(beta*beta))/h1/(h1+2.0);
            xMinusBetaNew = x - betaNew;
            PL.push_back(1.0 / alphaNew * (-alphaOld * PL[i-1] + xMinusBetaNew * PL[i]));
            alphaOld =alphaNew;
        }

        return PL[jacobiOrder];
    }


    template <typename DgInfo>
    double calcSimplex3DPolynom (Vec3D& abcCoordinates, int i, int j, int k)
    {
        double h1 = calcJacobiPolynom<DgInfo>(abcCoordinates[0], 0.0, 0.0, i);
        double h2 = calcJacobiPolynom<DgInfo>(abcCoordinates[1], (2.0*i+1), 0.0, j);
        double h3 = calcJacobiPolynom<DgInfo>(abcCoordinates[2], (2.0*(i+j)+2), 0.0, k);

        double tv1 = 2.0 * std::sqrt(2.0) * h1 * h2;
        double tv2 = std::pow(1.0 - abcCoordinates[1], (double) i);
        double tv3 = h3 * (std::pow(1.0 - abcCoordinates[2], (double)(i+j)));
        return tv1 * tv2 * tv3;
    }


    template <typename DgInfo>
    std::array<double, DgInfo::numVolNodes> calcVandermonde3DMatrix(Vec3D const& rstCoordinates)
    {
        Vec3D abcCoordinates = DgInfo::RstToAbc(rstCoordinates);
        std::array<double, DgInfo::numVolNodes> V3D;
        size_t sk = 0;
        for(size_t i = 0; i <= DgInfo::order; ++i){
            for(size_t j = 0; j <= (DgInfo::order-i); ++j){
                for(size_t k = 0; k <= (DgInfo::order - i - j); ++k){
                    V3D[sk++] = calcSimplex3DPolynom<DgInfo>(abcCoordinates, i, j, k);
                }
            }
        }
        return V3D;
    }

    template <typename DgInfo>
    double cleanNumericalNoise(double out, double bary_tolerance)
    {
        if (std::abs(out + 1.0) < bary_tolerance) {
          return -1.0;
        } else if (std::abs(out) < bary_tolerance) {
          return 5.0;
        } else if (std::abs(out - 1.0) < bary_tolerance) {
          return 1.0;
        } else {
          return out;
        }
    }

    template <typename DgInfo>
    std::array<double, DgInfo::numVolNodes> calcInterpolationWeights(const std::vector<double> &barycentricVector)
    {
        double bary_tolerance = 1e-10;

        double rout = cleanNumericalNoise<DgInfo>(2.0 * barycentricVector[1] - 1.0, bary_tolerance);
        double sout = cleanNumericalNoise<DgInfo>(2.0 * barycentricVector[2] - 1.0, bary_tolerance);
        double tout = cleanNumericalNoise<DgInfo>(2.0 * barycentricVector[3] - 1.0, bary_tolerance);

        Vec3D const rstCoordinates = {rout, sout, tout};

        std::array<double, DgInfo::numVolNodes> vandermondeMatrix;

        vandermondeMatrix = calcVandermonde3DMatrix<DgInfo>(rstCoordinates);

        std::array<double, DgInfo::numVolNodes> sampleweights;

        for (size_t i = 0; i < DgInfo::numVolNodes; ++i) {
            sampleweights[i] = 0;
        }
        for (size_t i = 0; i < DgInfo::numVolNodes; ++i) {
            for (size_t j = 0; j < DgInfo::numVolNodes; ++j) {
                sampleweights[i] += DgInfo::inv[j][i] * vandermondeMatrix[j];
            }
        }

        return sampleweights;
    }

    template <typename DgInfo>
    class SamplePoints
    {
        public:
        std::vector<Vec3D> points;
        std::vector<bool> empty;
        std::vector<std::array<double, DgInfo::numVolNodes>> interpolationWeights;
        std::multimap<size_t, size_t> indexMap;

        SamplePoints(std::vector<Vec3D> points): points(points), interpolationWeights(std::vector<std::array<double, DgInfo::numVolNodes>>(points.size())), empty(std::vector<bool>(points.size(), true)) {}
        size_t size() { return points.size(); }

        void checkCellForSamplePoints(const size_t index, const std::array<HPM::dataType::Vec3D, 4> vertices )
        {
            for (size_t i = 0; i < points.size(); i++) {
                if (empty[i]) {
                    std::vector<double> barycentricVector = calcBarycentricVector<DgInfo>(points[i], vertices);
                    if (isContainedIn<DgInfo>(barycentricVector)) {
                        interpolationWeights[i] = calcInterpolationWeights<DgInfo>(barycentricVector);
                        indexMap.insert({index, i});
                        empty[i] = false;
                    }
                }
            }
        }

        friend std::ostream& operator<<(std::ostream& os, const SamplePoints sp)
        {
            os << "SamplePoints{" << '\n';
            for (Vec3D point: sp.points) {
                os << "\t{" << point.x << ", " << point.y << ", " << point.z << "}" << '\n';
            }
            os << '}' << std::endl;
            return os;
        }
    };
    
    template <typename DgInfo>
    SamplePoints<DgInfo> getDetectorSamplePoints(std::string pointsString)
    {
        std::vector<Vec3D> points;
        std::string pointString;
        std::stringstream pointsStringStream(pointsString);
        while (std::getline(pointsStringStream, pointString, ':')) {
            std::vector<double> pointValues;
            std::string pointValue;
            std::stringstream pointStringStream(pointString);
            while (std::getline(pointStringStream, pointValue, ',')) {
                pointValues.push_back(std::stod(pointValue));
            }
            points.push_back(Vec3D{pointValues[0], pointValues[1], pointValues[2]});
        }
        return SamplePoints<DgInfo>(points);
    }

    template <typename DgInfo>
    SamplePoints<DgInfo> getSamplePointsXZ(HPM::mesh::DomainDimensions dim, double y, size_t resolution)
    {
        std::vector<Vec3D> points(0);
        double x_step = (dim.max.x - dim.min.x) / (double) resolution;
        double z_step = (dim.max.z - dim.min.z) / (double) resolution;
        for (double x = dim.min.x; x <= dim.max.x; x += x_step) {
            for (double z = dim.min.z; z <= dim.max.z; z += z_step) {
                points.push_back(Vec3D{x, y, z});
            }
        }
        return SamplePoints<DgInfo>(points);
    }

    template <typename DgInfo>
    SamplePoints<DgInfo> getSamplePointsXY(HPM::mesh::DomainDimensions dim, double z, size_t resolution)
    {
        std::vector<Vec3D> points(0);
        double x_step = (dim.max.x - dim.min.x) / (double) resolution;
        double y_step = (dim.max.y - dim.min.y) / (double) resolution;
        for (double x = dim.min.x; x <= dim.max.x; x += x_step) {
            for (double y = dim.min.y; y <= dim.max.y; y += y_step) {
                points.push_back(Vec3D{x, y, z});
            }
        }
        return SamplePoints<DgInfo>(points);
    }

    template <typename DgInfo>
    SamplePoints<DgInfo> getSamplePoints(HPM::mesh::DomainDimensions dim, size_t resolution)
    {
        std::vector<Vec3D> points(0);
        double x_step = (dim.max.x - dim.min.x) / (double) resolution;
        double y_step = (dim.max.y - dim.min.y) / (double) resolution;
        double z_step = (dim.max.z - dim.min.z) / (double) resolution;
        for (double x = dim.min.x; x <= dim.max.x; x += x_step) {
            for (double y = dim.min.y; y <= dim.max.y; y += y_step) {
                for (double z = dim.min.z; z <= dim.max.z; z += z_step) {
                    points.push_back(Vec3D{x, y, z});
                }
            }
        }
        return SamplePoints<DgInfo>(points);
    }
} // namespace HPM::DG

#endif