#include <HighPerMeshes.hpp>
#include <HighPerMeshes/drts/UsingOpenCL.hpp>

using namespace HPM;

int main()
{
    const auxiliary::ConfigParser CFG("config.cfg");

    const std::string oclPlatformName = CFG.GetValue<std::string>("oclPlatformName"); 
    const std::string oclDeviceName = CFG.GetValue<std::string>("oclDeviceName"); 

    OpenCLHandler _ocl(oclPlatformName, oclDeviceName);
    _ocl.LoadKernelsFromBinary("memcopy.aocx", {"memcopy"});


    //The runtime determines the configuration of HighPerMeshes. 
    //The GetBuffer class determines that we use a normal buffer to allocate space 
    drts::Runtime hpm{
        GetBuffer<OpenCLHandler::SVMAllocator>{}
    };

    // We also provide a parser to read from a config file if needed. In this example we get the path to a GAMBIT neutral file.
    const std::string meshFile = CFG.GetValue<std::string>("MeshFile"); 


    // The next step initializes a mesh
    // For this purpose we define a mesh class that needs two types as information: A CoordinateType that tells us which dimensionality the mesh has and how to store the coordinates and a topology class that can be used to define the mesh topology, i.e. how nodes are connected to each other. 
    // For this purpose we currently provide he Simplex class to define meshes for simplexes of any dimensionality

    // The CoordinateType tells us which data type and which dimensionality to use for a given mesh. 
    using CoordinateType = dataType::Vec<double, 3>;

    using Mesh = mesh::Mesh<CoordinateType, HPM::entity::Simplex>;
    const Mesh mesh = Mesh::template CreateFromFile<HPM::auxiliary::GambitMeshFileReader>(meshFile);

    //We store the CellDimension of the mesh because we're goint to use it more often later in the example.
    static constexpr auto CellDimension = Mesh::CellDimension;

    // We can determine what entities we want to iterate over by using the member functions of the mesh
    // `mesh.GetEntityRange` allows us to iterate over all entities of a certain dimension, in this case we want to iterate over each cell
    const auto AllCells{
        mesh.GetEntityRange<CellDimension>()
    };

    // Degrees of freedom, which are referred to `dofs` most of the time, allow us to associate entities of the mesh with space in a buffer.
    // In this example, we define just one degree of freedom for each face and each cell in the mesh 
    constexpr auto dofs= dof::MakeDofs<0, 0, 0, 1, 0>();

    // Here we allocate a buffer and see the benefits of leaving the buffer generation to the runtime system.
    // Independent of technology used, we can call the runtime's GetBuffer function to allocate data of type int for each entity in the mesh corresponding to the specified degrees of freedom.

    auto buffer_in{
        hpm.GetBuffer<int>(mesh, dofs, _ocl.GetSVMAllocator<int>())
    };

    auto buffer_out{
        hpm.GetBuffer<int>(mesh, dofs, _ocl.GetSVMAllocator<int>())
    };

    _ocl.SetKernelArg("memcopy", 0, buffer_in);
    _ocl.SetKernelArg("memcopy", 1, buffer_out);
    _ocl.SetKernelArg("memcopy", 2, int(buffer_in.GetSize()));

    // The dispatcher is used to dispatch kernels to some parallelization strategy.
    // In this case, we use a SequentialDispatcher to just execute the specified kernels
    SequentialDispatcher dispatcher;

    _ocl.MapSVMBuffer(buffer_out);
    _ocl.MapSVMBuffer(buffer_in);

    dispatcher.Execute(
        iterator::Range{1}, 
        ForEachEntity(
            AllCells,
            std::tuple(Write(Cell(buffer_in))),
            [&](const auto &, auto, auto local_view) {

                auto& bufferAccess_out = dof::GetDofs<CellDimension>(std::get<0>(local_view));

                const auto dof = 0;

                bufferAccess_out[dof] = 43; 

            })
    );

#if 1
    _ocl.UnmapSVMBuffer(buffer_out);
    _ocl.UnmapSVMBuffer(buffer_in);

    _ocl.EnqueueKernel("memcopy");

    _ocl.MapSVMBuffer(buffer_out);
    _ocl.MapSVMBuffer(buffer_in);
#else
    dispatcher.Execute(
        iterator::Range{1}, 
        ForEachEntity(
            AllCells,
            std::tuple(Write(Cell(buffer_out)),Read(Cell(buffer_in))),
            [&](const auto &, auto, auto local_view) {

                auto& bufferAccess_out = dof::GetDofs<CellDimension>(std::get<0>(local_view));
                auto& bufferAccess_in = dof::GetDofs<CellDimension>(std::get<1>(local_view));

                const auto dof = 0;

                bufferAccess_out[dof] = bufferAccess_in[dof];

            })
    );
#endif

    dispatcher.Execute(
        iterator::Range{1}, 
        ForEachEntity(
            AllCells,
            std::tuple(Read(Cell(buffer_out))),
            [&](const auto &, auto, auto local_view) {

                auto& bufferAccess_in = dof::GetDofs<CellDimension>(std::get<0>(local_view));

                const auto dof = 0;

                if(bufferAccess_in[dof] != 43) std::cout<<bufferAccess_in[dof]<<' ';

            })
    );

    std::cout<<std::endl;

}
