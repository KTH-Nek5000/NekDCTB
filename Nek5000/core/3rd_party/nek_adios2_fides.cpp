#include <adios2.h>
#include <mpi.h>
#include <string>
#include <iostream>
#include <ctime>
#include <vector>

int lx1, ly1, lz1, lx2, ly2, lz2;
int lxyz1, lxyz2, lxyz3;
int nelv, nelt;
int nelgv, nelgt;
int nelbv, nelbt;

std::size_t total1, start1, count1;
std::size_t total3, start3, count3;
std::size_t init_total, init_start, init_count;
adios2::ADIOS adios;
adios2::IO io;
adios2::Engine writer;
adios2::Variable<int> init_int_const;
adios2::Variable<double> init_double_const;
std::vector<int> vINT;
std::vector<double> vDOUBLE;
std::vector<int> connectivity;
std::vector<float> points;
adios2::Variable<double> x;
adios2::Variable<double> y;
adios2::Variable<double> z;
adios2::Variable<double> p;
adios2::Variable<double> vx;
adios2::Variable<double> vy;
adios2::Variable<double> vz;
adios2::Variable<double> t;
adios2::Variable<int> vstep;
adios2::Variable<int> ADIOS_connectivity;
adios2::Variable<float> ADIOS_points;
double dataTime = 0.0;
std::clock_t startT;
std::clock_t startTotal;
int rank, size, i;
int lxy1, lyz3, lx3, ly3, lz3;
bool if3d;
bool firstStep;
int step = 0;

void convert_points_connectivity(
    const double *x_in,
    const double *y_in,
    const double *z_in,
    float *points_out,
    int *connectivity_out,
    const int lx1_in,
    const int ly1_in,
    const int lz1_in,
    const int nelt_in,
    const int nelbt_in,
    const bool if3d_in)
{
    /**
     * @param x_in, y_in, z_in the position of points in x, y, z directions
     * @param lx_in, ly_in, lz_in the number of points in one element in x, y, z  directions
     * @param nelt_in the number of elements assigned to this rank
     * @param nelbt_in the index of the first elements assigned to this rank
     * @param if3d_in the bool indicate if the simulation case is 3d 
     * @param points_out all the positions of the points
     * @param connectivity_out the index of points each cell connected to 
     */
    int idx, ii, jj, kk;
    int points_offset, start_index, index;
    int lxyz1_ = lx1_in * ly1_in * lz1_in;
    int lxy1_ = lx1_in * ly1_in;
    int lx3_ = lx1_in - 1;
    int ly3_ = ly1_in - 1;
    int lz3_ = lz1_in - 1;
    if (!if3d_in)
        lz3_ = lz1_in;
    int lyz3_ = ly3_ * lz3_;
    int lxyz3_ = lx3_ * lyz3_;
    if (if3d_in)
    {
        for (idx = 0; idx < nelt_in; ++idx)
        {
            points_offset = lxyz1_ * (idx + nelbt_in);
            start_index = lxyz3_ * idx * 8;
            for (ii = 0; ii < lx3_; ++ii)
            {
                for (jj = 0; jj < ly3_; ++jj)
                {
                    for (kk = 0; kk < lz3_; ++kk)
                    {
                        index = start_index + (ii * lyz3_ + jj * lz3_ + kk) * 8;
                        connectivity_out[index] = kk * lxy1_ + jj * lx1_in + ii + points_offset;
                        connectivity_out[index + 1] = kk * lxy1_ + jj * lx1_in + (ii + 1) + points_offset;
                        connectivity_out[index + 2] = kk * lxy1_ + (jj + 1) * lx1_in + (ii + 1) + points_offset;
                        connectivity_out[index + 3] = kk * lxy1_ + (jj + 1) * lx1_in + ii + points_offset;
                        connectivity_out[index + 4] = (kk + 1) * lxy1_ + jj * lx1_in + ii + points_offset;
                        connectivity_out[index + 5] = (kk + 1) * lxy1_ + jj * lx1_in + (ii + 1) + points_offset;
                        connectivity_out[index + 6] = (kk + 1) * lxy1_ + (jj + 1) * lx1_in + (ii + 1) + points_offset;
                        connectivity_out[index + 7] = (kk + 1) * lxy1_ + (jj + 1) * lx1_in + ii + points_offset;
                    }
                }
            }
        }
    }
    else
    {
        for (idx = 0; idx < nelt_in; ++idx)
        {
            points_offset = lxyz1_ * (idx + nelbt_in);
            start_index = lxyz3_ * idx * 4;
            for (ii = 0; ii < lx3_; ++ii)
            {
                for (jj = 0; jj < ly3_; ++jj)
                {
                    index = start_index + (ii * ly3_ + jj) * 4;
                    connectivity_out[index] = jj * lx1_in + ii + points_offset;
                    connectivity_out[index + 1] = jj * lx1_in + (ii + 1) + points_offset;
                    connectivity_out[index + 2] = (jj + 1) * lx1_in + (ii + 1) + points_offset;
                    connectivity_out[index + 3] = (jj + 1) * lx1_in + ii + points_offset;
                }
            }
        }
    }
    for (idx = 0; idx < nelt_in * lxyz1_; ++idx)
    {
        points_out[i * 3] = static_cast<float>(x_in[idx]);
        points_out[i * 3 + 1] = static_cast<float>(y_in[idx]);
        points_out[i * 3 + 2] = static_cast<float>(z_in[idx]);
    }
}

extern "C" void adios2_fides_setup_(
    const int *lx1_in,
    const int *ly1_in,
    const int *lz1_in,
    const int *nelv_in,
    const int *nelbv_in,
    const int *nelgv_in,
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
    const int *comm_int)
{
    startTotal = std::clock();
    std::string configFile = "config/config.xml";
    MPI_Comm comm = MPI_Comm_f2c(*comm_int);
    adios = adios2::ADIOS(configFile, comm);
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);
    io = adios.DeclareIO("writer");

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

    /*
    Assume the simulated 2d points are as following:
    0 - 1 - 2 - 3 - 4 - 5 - 6 - 7
    |   |   |   |   |   |   |   |
    8 - 9 - 10- 11- 12- 13- 14- 15 
    |   |   |   |   |   |   |   |
    16- 17- 18- 19- 20- 21- 22- 23
    |   |   |   |   |   |   |   |
    24- 25- 26- 27- 28- 29- 30- 31
    |   |   |   |   |   |   |   |
    32- 33- 34- 35- 36- 37- 38- 39
    |   |   |   |   |   |   |   |
    40- 41- 42- 43- 44- 45- 46- 47
    |   |   |   |   |   |   |   |
    48- 49- 50- 51- 52- 53- 54- 55
    |   |   |   |   |   |   |   |
    56- 57- 58- 59- 60- 61- 62- 63
    Nek5000 would first devide them into following 16 elements, if lx=2 and ly=2. Each element has one mesh with 2x2 points and one cell for the vtk.

    |-------|-------|-------|-------|
    | ele 0 | ele 1 | ele 2 | ele 3 |
    | 0 - 1 | 2 - 3 | 4 - 5 | 6 - 7 |
    | |   | | |   | | |   | | |   | |
    | 8 - 9 | 10- 11| 12- 13| 14- 15|
    |-------|-------|-------|-------|
    | ele 4 | ele 5 | ele 6 | ele 7 |
    | 16- 17| 18- 19| 20- 21| 22- 23|
    | |   | | |   | | |   | | |   | |
    | 24- 25| 26- 27| 28- 29| 30- 31|
    |-------|-------|-------|-------|
    | ele 8 | ele 9 | ele 10| ele 11|
    | 32- 33| 34- 35| 36- 37| 38- 39|
    | |   | | |   | | |   | | |   | |
    | 40- 41| 42- 43| 44- 45| 46- 47|
    |-------|-------|-------|-------|
    | ele 12| ele 13| ele 14| ele 15|
    | 48- 49| 50- 51| 52- 53| 54- 55|
    | |   | | |   | | |   | | |   | |
    | 56- 57| 58- 59| 60- 61| 62- 63|
    |-------|-------|-------|-------|
    Work load distribution of Nek5000 would then based on the element. For instance, if four MPI ranks are used, ele 0,1,4,and 5 would be assigned to rank 0; ele 2,3,6,and 7 be would assigned to rank 1; ele 8,9,12,and 13 would be assigned to rank 2; ele 10,11,14,and 15 would be assigned to rank 3. In this case, nelv = 4, nelbv = 0,4,8,12. The indexes of all the points would be:
    |-------|-------|-------|-------|
    | num 0 | num 1 | num 4 | num 5 |
    | 0 - 1 | 4 - 5 | 16- 17| 20- 21|
    | |   | | |   | | |   | | |   | |
    | 2 - 3 | 6 - 7 | 18- 19| 22- 23|
    |-------|-------|-------|-------|
    | num 2 | num 3 | num 6 | num 7 |
    | 8 - 9 | 12- 13| 24- 25| 28- 29|
    | |   | | |   | | |   | | |   | |
    | 10- 11| 14- 15| 26- 27| 30- 31|
    |-------|-------|-------|-------|
    | num 8 | num 9 | num 12| num 13|
    | 32- 33| 36- 37| 48- 49| 52- 53|
    | |   | | |   | | |   | | |   | |
    | 34- 35| 38- 39| 50- 51| 54- 55|
    |-------|-------|-------|-------|
    | num 10| num 11| num 14| num 15|
    | 40- 41| 44- 45| 56- 57| 60- 61|
    | |   | | |   | | |   | | |   | |
    | 42- 43| 46- 47| 58- 59| 62- 63|
    |-------|-------|-------|-------|
    If 16 MPI ranks are used, ele 0 would be assigned to rank0; ele 1 would be assigned ro rank 1, etc. The indexes of all the points would be:
    |-------|-------|-------|-------|
    | num 0 | num 1 | num 2 | num 3 |
    | 0 - 1 | 4 - 5 | 8 - 9 | 12- 13|
    | |   | | |   | | |   | | |   | |
    | 2 - 3 | 6 - 7 | 10- 11| 14- 15|
    |-------|-------|-------|-------|
    | num 4 | num 5 | num 6 | num 7 |
    | 16- 17| 20- 21| 24- 25| 28- 29|
    | |   | | |   | | |   | | |   | |
    | 18- 19| 22- 23| 26- 27| 30- 31|
    |-------|-------|-------|-------|
    | num 8 | num 9 | num 10| num 11|
    | 32- 33| 36- 37| 40- 41| 44- 45|
    | |   | | |   | | |   | | |   | |
    | 34- 35| 38- 39| 42- 43| 46- 47|
    |-------|-------|-------|-------|
    | num 12| num 13| num 14| num 15|
    | 48- 49| 52- 53| 56- 57| 60- 61|
    | |   | | |   | | |   | | |   | |
    | 50- 51| 54- 55| 58- 59| 62- 63|
    |-------|-------|-------|-------|
    We don't build cells between elements, as long as number of reader ranks is the factor of the total number of elements or all the points within one element would be read by one reader rank, visualization and data analysis doesn't require these cells would be correct, regardless of the writer ranks, because Nek5000 balances the workload and covers all the elements. 
    */

    lx1 = *lx1_in; // number of points in one element in x direction
    ly1 = *ly1_in; // number of points in one element in y direction
    lz1 = *lz1_in; // number of points in one element in z direction

    nelgv = *nelgv_in; // total number of elements
    if3d = true;       //
    if (*if3d_in == 0)
        if3d = false; // bool variable to show if the case is 3d
    nelv = *nelv_in;  // number of elements calculated by this rank

    std::vector<int> temp(size);
    MPI_Allgather(&nelv, 1, MPI_INT, temp.data(), 1, MPI_INT, comm);
    nelbv = 0; // the index of the first element calculated by this rank
    for (i = 0; i < rank; ++i)
    {
        nelbv += temp[i];
    }

    init_count = 7;
    init_total = init_count * static_cast<std::size_t>(size);
    init_start = init_count * static_cast<std::size_t>(rank);
    vINT.resize(init_count);
    vINT[0] = lx1;
    vINT[1] = ly1;
    vINT[2] = lz1;
    vINT[3] = nelv;
    vINT[4] = nelgv;
    vINT[5] = *iostep_in; // in-situ step is executed every *iostep_in simulation step
    vINT[6] = *if3d_in;
    init_int_const = io.DefineVariable<int>("INT_CONST", {init_total}, {init_start}, {init_count});

    init_count = 2;
    init_total = init_count * static_cast<std::size_t>(size);
    init_start = init_count * static_cast<std::size_t>(rank);
    vDOUBLE.resize(init_count);
    vDOUBLE[0] = *t_start_in; // the physical initial time Nek50000 would simulate
    vDOUBLE[1] = *dt_in;      // the physical time changed between two simulation step
    init_double_const = io.DefineVariable<double>("DOUBLE_CONST", {init_total}, {init_start}, {init_count});
    lxyz1 = lx1 * ly1 * lz1;
    lxy1 = lx1 * ly1;
    lx3 = lx1 - 1;
    ly3 = ly1 - 1;
    lz3 = lz1 - 1;
    lyz3 = ly3 * lz3;
    if (if3d)
        lxyz3 = lx3 * ly3 * lz3;
    else
        lxyz3 = lx3 * ly3;
    total1 = static_cast<std::size_t>(lxyz1 * nelgv);
    start1 = static_cast<std::size_t>(lxyz1 * nelbv);
    count1 = static_cast<std::size_t>(lxyz1 * nelv);
    total3 = static_cast<std::size_t>(lxyz3 * nelgv);
    start3 = static_cast<std::size_t>(lxyz3 * nelbv);
    count3 = static_cast<std::size_t>(lxyz3 * nelv);

    /*The pressure, velocity and temperature of each point*/
    p = io.DefineVariable<double>("P", {total1}, {start1}, {count1});
    vx = io.DefineVariable<double>("VX", {total1}, {start1}, {count1});
    vy = io.DefineVariable<double>("VY", {total1}, {start1}, {count1});
    vz = io.DefineVariable<double>("VZ", {total1}, {start1}, {count1});
    t = io.DefineVariable<double>("T", {total1}, {start1}, {count1});
    vstep = io.DefineVariable<int>("step");

    /* Fides schema */
    io.DefineAttribute<std::string>("Fides_Data_Model", "unstructured_single");
    if (if3d)
    {
        io.DefineAttribute<std::string>("Fides_Cell_Type", "hexahedron");
        ADIOS_connectivity = io.DefineVariable<int>("connectivity", {total3, 8}, {start3, 0}, {count3, 8});
        connectivity.resize(count3 * 8);
    }
    else
    {
        io.DefineAttribute<std::string>("Fides_Cell_Type", "quad");
        ADIOS_connectivity = io.DefineVariable<int>("connectivity", {total3, 4}, {start3, 0}, {count3, 4});
        connectivity.resize(count3 * 4);
    }

    ADIOS_points = io.DefineVariable<float>("points", {total1, 3}, {start1, 0}, {count1, 3});
    points.resize(count1 * 3);
    convert_points_connectivity(xml_in, yml_in, zml_in, points.data(), connectivity.data(), lx1, ly1, lz1, nelv, nelbv, if3d);

    io.DefineAttribute<std::string>("Fides_Coordinates_Variable", "points");
    io.DefineAttribute<std::string>("Fides_Connectivity_Variable", "connectivity");
    //io.DefineAttribute<std::string>("Fides_Time_Variable", "step");

    std::vector<std::string> varList = {"P", "T", "VX", "VY", "VZ"};
    std::vector<std::string> assocList = {"points", "points", "points", "points", "points"};
    io.DefineAttribute<std::string>("Fides_Variable_List", varList.data(), varList.size());
    io.DefineAttribute<std::string>("Fides_Variable_Associations", assocList.data(), assocList.size());

    writer = io.Open("globalArray", adios2::Mode::Write);
    /* 
    Here is the first step in adios2 
    Some Nek5000 related data are provided in this step
    */
    writer.BeginStep();
    writer.Put<int>(init_int_const, vINT.data());
    writer.Put<double>(init_double_const, vDOUBLE.data());
    /*
    But vtk related information such as points and connectivity are communicated in this step 
    Also the initail values of p (pressure), vx (velocity in x direction), vy (velocity in y direction), vz (velocity in z direction), and t (temperature) are communicated. 
    */
    writer.Put<int>(ADIOS_connectivity, connectivity.data());
    writer.Put<float>(ADIOS_points, points.data());
    writer.Put<double>(p, pr);
    writer.Put<double>(vx, v);
    writer.Put<double>(vy, u);
    writer.Put<double>(vz, w);
    writer.Put<double>(t, tem);
    writer.Put<int>(vstep, step);
    writer.EndStep();
    if (!rank)
        std::cout << "In-Situ setting done" << std::endl;
    std::cout << "Nek rank: " << rank << " count: " << nelv << " , start: " << nelbv << " , total: " << nelgv << std::endl;
}

extern "C" void adios2_fides_update_(
    const double *pr,
    const double *v,
    const double *u,
    const double *w,
    const double *temp)
{
    startT = std::clock();
    writer.BeginStep();
    writer.Put<double>(p, pr);
    writer.Put<double>(vx, v);
    writer.Put<double>(vy, u);
    writer.Put<double>(vz, w);
    writer.Put<double>(t, temp);
    ++step;
    writer.Put<int>(vstep, step);
    writer.EndStep();
    dataTime += (std::clock() - startT) / (double)CLOCKS_PER_SEC;
}

extern "C" void adios2_fides_finalize_()
{
    writer.Close();

    std::cout << "rank: " << rank << " sin-situ time: " << dataTime << "s, total time: " << (std::clock() - startTotal) / (double)CLOCKS_PER_SEC << "s. " << std::endl;
}
