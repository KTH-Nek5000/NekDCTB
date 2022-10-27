#include <adios2.h>
#include <mpi.h>
#include <string>
#include <iostream>
#include <ctime>
#include <vector>

int lx1, ly1, lz1, lx2, ly2, lz2;
int lxyz1, lxyz2; 
int nelv, nelt;
int nelgv, nelgt;
int nelbv, nelbt;

std::size_t total1, start1, count1;
std::size_t total2, start2, count2;
std::size_t total3, start3, count3;
std::size_t init_total, init_start, init_count;
adios2::ADIOS adios;
adios2::IO io;
adios2::Engine writer;
adios2::Variable<int> init_int_const;
adios2::Variable<double> init_double_const;
std::vector<int> vINT;
std::vector<double> vDOUBLE;
adios2::Variable<double> x;
adios2::Variable<double> y;
adios2::Variable<double> z;
adios2::Variable<double> p;
adios2::Variable<double> vx;
adios2::Variable<double> vy;
adios2::Variable<double> vz;
adios2::Variable<double> t;
adios2::Variable<double> bm1;
double dataTime = 0.0;
std::clock_t startT;
std::clock_t startTotal;
int rank, size, i;
bool if3d;
bool firstStep;
extern "C" void async_catalyst_setup_(
    const int *lx1_in,
    const int *ly1_in,
    const int *lz1_in,
    const int *lx2_in,
    const int *ly2_in,
    const int *lz2_in,
    const int *nelv_in,
    const int *nelt_in,
    const int *nelbv_in,
    const int *nelbt_in,
    const int *nelgv_in,
    const int *nelgt_in,
    const int *iostep_in,
    const int *if3d_in,
    const double *t_start_in,
    const double *dt_in,
    const double *xml_in,
    const double *yml_in,
    const double *zml_in,
    const double *pr,
    const double *v,
    const double *u,
    const double *w,
    const double *tem,
    const double *bm,
    const int *comm_int
){
    /*
    
    
    */
    startTotal = std::clock();
    std::string configFile="config/config.xml";
    MPI_Comm comm = MPI_Comm_f2c(*comm_int);
    adios = adios2::ADIOS(configFile, comm);
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);
    io = adios.DeclareIO("writer");
    writer = io.Open("globalArray", adios2::Mode::Write, comm);
    /*
    adios2::IO io_xyz = adios.DeclareIO("writer0");
    std::string file = "geo" + std::to_string(size) + ".bp";
    adios2::Engine writer0 = io_xyz.Open(file, adios2::Mode::Write,comm);
    */
    if (!io.InConfigFile())
    {
        // if not defined by user, we can change the default settings
        // BPFile is the default engine
        io.SetEngine("BPFile");
        io.SetParameters({{"num_threads", "1"}});

        // ISO-POSIX file output is the default transport (called "File")
        // Passing parameters to the transport
#ifdef _WIN32
        io.AddTransport("File", {{"Library", "stdio"}});
#else
        io.AddTransport("File", {{"Library", "posix"}});
#endif
    }
    //writer0.BeginStep();
    writer.BeginStep();

    lx1 = *lx1_in;
    ly1 = *ly1_in;
    lz1 = *lz1_in;
    lx2 = *lx2_in;
    ly2 = *ly2_in;
    lz2 = *lz2_in;
    nelgt = *nelgt_in;
    nelgv = *nelgv_in;
    if3d = true;
    if(*if3d_in==0) if3d = false;
    nelt = *nelt_in;
    nelv = *nelv_in;
    /* To make sure the in-situ core with rank i recieves the data from the simulation core with rank i*(size/insitu_size) to i*(size/insitu_size+1)-1 */ 
    //nelbv = *nelbv_in;
    //nelbt = *nelbt_in;
    
    std::vector<int>temp(size);
    std::vector<int>temp1(size);
    MPI_Allgather(&nelv, 1, MPI_INT, temp.data(), 1,MPI_INT, comm);
    MPI_Allgather(&nelt, 1, MPI_INT, temp1.data(), 1,MPI_INT, comm);
    nelbt = 0;
    nelbv = 0;
    for(i=0;i<rank;++i){
        nelbv += temp[i];
        nelbt += temp1[i];
    }
    
    init_count = 12;
    init_total = init_count*static_cast<std::size_t>(size);
    init_start = init_count*static_cast<std::size_t>(rank);
    vINT.resize(init_count);
    vINT[0] = lx1;
    vINT[1] = ly1;
    vINT[2] = lz1;
    vINT[3] = lx2;
    vINT[4] = ly2;
    vINT[5] = lz2;
    vINT[6] = nelv;
    vINT[7] = nelt;
    vINT[8] = nelgv;
    vINT[9] = nelgt;
    vINT[10] = *iostep_in;
    vINT[11] = *if3d_in;
    init_int_const = io.DefineVariable<int>("INT_CONST", {init_total}, {init_start}, {init_count});
    
    init_count = 2;
    init_total = init_count*static_cast<std::size_t>(size);
    init_start = init_count*static_cast<std::size_t>(rank);
    vDOUBLE.resize(init_count);
    vDOUBLE[0] = *t_start_in;
    vDOUBLE[1] = *dt_in;
    init_double_const = io.DefineVariable<double>("DOUBLE_CONST", {init_total}, {init_start}, {init_count});
    lxyz1 = lx1 * ly1 * lz1;
    lxyz2 = lx2 * ly2 * lz2;
    total1 = static_cast<std::size_t>(lxyz1 * nelgv);
    start1 = static_cast<std::size_t>(lxyz1 * nelbv);
    count1 = static_cast<std::size_t>(lxyz1 * nelv);
    total2 = static_cast<std::size_t>(lxyz2 * nelgv);
    start2 = static_cast<std::size_t>(lxyz2 * nelbv);
    count2 = static_cast<std::size_t>(lxyz2 * nelv);
    total3 = static_cast<std::size_t>(lxyz1 * nelgt);
    start3 = static_cast<std::size_t>(lxyz1 * nelbt);
    count3 = static_cast<std::size_t>(lxyz1 * nelt);
    x = io.DefineVariable<double>("X", {total3}, {start3}, {count3});
    y = io.DefineVariable<double>("Y", {total3}, {start3}, {count3});
    z = io.DefineVariable<double>("Z", {total3}, {start3}, {count3});
    p = io.DefineVariable<double>("P", {total2}, {start2}, {count2});
    vx = io.DefineVariable<double>("VX", {total1}, {start1}, {count1});
    vy = io.DefineVariable<double>("VY", {total1}, {start1}, {count1});
    vz = io.DefineVariable<double>("VZ", {total1}, {start1}, {count1});
    t = io.DefineVariable<double>("T", {total3}, {start3}, {count3});
    bm1 = io.DefineVariable<double>("BM1", {total3}, {start3}, {count3});
    
    //writer = io.Open("globalArray", adios2::Mode::Write);
    //writer.BeginStep();
    writer.Put<int>(init_int_const, vINT.data());
    writer.Put<double>(init_double_const, vDOUBLE.data());
    writer.EndStep();
    writer.BeginStep();
    writer.Put<double>(x, xml_in);
    writer.Put<double>(y, yml_in);
    writer.Put<double>(z, zml_in);
    writer.Put<double>(p, pr);
    writer.Put<double>(vx, v);
    writer.Put<double>(vy, u);
    writer.Put<double>(vz, w);
    writer.Put<double>(t, tem);
    writer.Put<double>(bm1, bm);
    writer.EndStep();
    /*
    adios2::Variable<double> x_xyz = io_xyz.DefineVariable<double>("X", {total1}, {start1}, {count1});
    adios2::Variable<double> y_xyz = io_xyz.DefineVariable<double>("Y", {total1}, {start1}, {count1});
    adios2::Variable<double> z_xyz = io_xyz.DefineVariable<double>("Z", {total1}, {start1}, {count1});
    writer0.Put<double>(x_xyz, xml_in);
    writer0.Put<double>(y_xyz, yml_in);
    writer0.Put<double>(z_xyz, zml_in);
    writer0.EndStep();
    writer0.Close();
    */
    firstStep=true;
    if(!rank) std::cout << "In-Situ setting done" << std::endl;
    std::cout << "Nek rank: " << rank << " count: " << nelt << " , start: " << nelbt << " , total: " << nelgt << " , nelbv: " << *nelbt_in << std::endl;
}

extern "C" void async_catalyst_update_(
    const double *xml_in,
    const double *yml_in,
    const double *zml_in,
    const double *pr,
    const double *v,
    const double *u,
    const double *w,
    const double *temp,
    const double *bm
){
    startT = std::clock();
    writer.BeginStep();
    /*
    if(firstStep){
    	writer.Put<double>(x, xml_in,adios2::Mode::Sync);
    	writer.Put<double>(y, yml_in,adios2::Mode::Sync);
    	writer.Put<double>(z, zml_in,adios2::Mode::Sync);
	firstStep=false;
    }
    */
    writer.Put<double>(p, pr);
    writer.Put<double>(vx, v);
    writer.Put<double>(vy, u);
    writer.Put<double>(vz, w);
    writer.Put<double>(t, temp);
    writer.Put<double>(bm1, bm);
    writer.EndStep();
    dataTime += (std::clock() - startT) / (double) CLOCKS_PER_SEC;
}

extern "C" void async_catalyst_finalize_(){
    writer.Close();
    
    std::cout <<  "rank: " << rank << " sin-situ time: " << dataTime << "s, total time: " << (std::clock() - startTotal) / (double) CLOCKS_PER_SEC << "s. " << std::endl;
}


