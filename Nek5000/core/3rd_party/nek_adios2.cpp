#include <adios2.h>
#include <string>
#include <iostream>
#include <ctime>

adios2::ADIOS adios;
adios2::IO io;
adios2::IO ior;
adios2::IO io_head;
adios2::Engine writer;
adios2::Engine writer_head;
adios2::Engine readr; 
adios2::Variable<double> x;
adios2::Variable<double> y;
adios2::Variable<double> z;
adios2::Variable<int> lglelw;
adios2::Variable<float> p;
adios2::Variable<float> vx;
adios2::Variable<float> vy;
adios2::Variable<float> vz;
adios2::Variable<double> bm1;
adios2::Variable<float> vxr;
adios2::Variable<float> vyr;
adios2::Variable<float> vzr;
adios2::Variable<float> prr;
adios2::Variable<int> lglelr;
std::vector<float> vVXr;
std::vector<float> vVYr;
std::vector<float> vVZr;
std::vector<float> vVPrr;
std::vector<float> vVXw;
std::vector<float> vVYw;
std::vector<float> vVZw;
std::vector<float> vVPrw;
std::vector<int> vVLGLELr;
// adios2::Variable<double> t;
double dataTime = 0.0;
std::clock_t startT;
std::clock_t elapsedT;
int rank, size;
int ifile;
int ifilew;
int nxyze;
bool firstWrite;


extern "C" void adios2_setup_(
    const int *nval,
    const int *nelvin,
    const int *nelb,
    const int *nelgv,
    const int *nelgt,
    const double *xml,
    const double *yml,
    const double *zml,
    const int *comm_int
){

    // Set Up Adios2 engines
    std::string configFile="config/config.xml";
    MPI_Comm comm = MPI_Comm_f2c(*comm_int);
    adios = adios2::ADIOS(configFile, comm);
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);
    io = adios.DeclareIO("writer");
    io_head = adios.DeclareIO("writer0");
    if (!io.InConfigFile())
    {
        io.SetEngine("BPFile");
        io.SetParameters({{"num_threads", "1"}});
    }
    ior = adios.DeclareIO("inputReader");
    if (!ior.InConfigFile())
    {
        ior.SetEngine("BPFile");
        ior.SetParameters({{"num_threads", "1"}});
    }
    
    // Partition data residing in the velocity mesh. //
    // Number of elements
    unsigned int nelv = static_cast<unsigned int>((*nelvin));
    // First element in my rank
    unsigned int start = static_cast<unsigned int>(*nelb);
    // Start entry in global array
    start *= static_cast<unsigned int>(*nval);
    // n is count, i.e number of entries in the array in my rank
    unsigned int n = static_cast<unsigned int> (*nval) * nelv;
    // gn is the total size of the arrays, not per io rank 
    unsigned int gn = static_cast<unsigned int>((*nelgv)*(*nval));
    std::cout << rank << ": " << gn << ", " << start << "," << n << std::endl;
    // Define arrays
    x = io_head.DefineVariable<double>("X", {gn}, {start}, {n});
    y = io_head.DefineVariable<double>("Y", {gn}, {start}, {n});
    z = io_head.DefineVariable<double>("Z", {gn}, {start}, {n});
    p = io.DefineVariable<float>("P_OUT", {gn}, {start}, {n});
    vx = io.DefineVariable<float>("VX_OUT", {gn}, {start}, {n});
    vy = io.DefineVariable<float>("VY_OUT", {gn}, {start}, {n});
    vz = io.DefineVariable<float>("VZ_OUT", {gn}, {start}, {n});
    bm1 = io.DefineVariable<double>("BM1", {gn}, {start}, {n});
    nxyze = n;
    firstWrite = true;

    // Partition data residing in the temperature mesh. //
    unsigned int nelt = static_cast<unsigned int>((*nelvin));
    start = static_cast<unsigned int>(*nelb);
    ifile = 0 ;
    ifilew = 0 ;
    start = start * static_cast<unsigned int>(*nval);
    n = static_cast<unsigned int> (*nval) * nelt;
    gn = static_cast<unsigned int>((*nelgt)*(*nval));

    // Partition element data //
    nelt = static_cast<unsigned int>((*nelvin));
    start = static_cast<unsigned int>(*nelb);
    ifile = 0 ;
    ifilew = 0 ;
    start = start;
    n = static_cast<unsigned int> (nelt);
    gn = static_cast<unsigned int>((*nelgt));
    lglelw = io.DefineVariable<int>("LGLEL_OUT", {gn}, {start}, {n});

    // Write the geometry of the case //
    writer_head = io_head.Open("geo.bp", adios2::Mode::Write);
    writer_head.Put<double>(x, xml);
    writer_head.Put<double>(y, yml);
    writer_head.Put<double>(z, yml);
    writer_head.Close();
    if(!rank)
	std::cout << "geo.bp written" << std::endl;
}

extern "C" void adios2_update_(
    const int *lglel,
    const double *pr,
    const double *v,
    const double *u,
    const double *w,
    const double *temp
){
    startT = std::clock();
    ifilew=ifilew+1;
    std:: string fileName= "out.f0000"+ std::to_string(ifilew) +".bp";

    if (firstWrite)
    {      
	vVXw.resize(nxyze);
	vVYw.resize(nxyze);
	vVZw.resize(nxyze);
	vVPrw.resize(nxyze);

	firstWrite=false;
    }

    for (int i=0; i<nxyze; i++){
    	vVPrw[i] = static_cast<float>(pr[i]);
    	vVXw[i]  = static_cast<float>(v[i]);
    	vVYw[i]  = static_cast<float>(u[i]);
    	vVZw[i]  = static_cast<float>(w[i]);
    }

    writer = io.Open(fileName, adios2::Mode::Write);
    writer.BeginStep();
    writer.Put<float>(p,vVPrw.data());
    writer.Put<float>(vx,vVXw.data());
    writer.Put<float>(vy,vVYw.data());
    writer.Put<float>(vz,vVZw.data());
    writer.Put<int>(lglelw,lglel);
    writer.EndStep();
    writer.Close();
    dataTime += (std::clock() - startT) / (double) CLOCKS_PER_SEC;
}

extern "C" void adios2_stream_(
    const int *lglel,
    const float *pr,
    const float *v,
    const float *u,
    const float *w,
    const double *mass1,
    const double *temp
){
    startT = std::clock();
    writer.BeginStep();
    writer.Put<float>(p, pr);
    writer.Put<float>(vx, v);
    writer.Put<float>(vy, u);
    writer.Put<float>(vz, w);
    writer.Put<double>(bm1, mass1);
    writer.Put<int>(lglelw, lglel);
    writer.EndStep();
    dataTime += (std::clock() - startT) / (double) CLOCKS_PER_SEC;
}

extern "C" void adios2_read_(
    int *lglelrr,
    double *pr,
    double *v,
    double *u,
    double *w,
    const int *nval,
    const int *nelvin,
    const int *nelb,
    const int *nelgv,
    const int *nelgt,
    const int *comm_int,
    char *fname
){

    startT = std::clock();
    // Partition velocity mesh data //
    unsigned int nelv = static_cast<unsigned int>((*nelvin));
    unsigned int start = static_cast<unsigned int>(*nelb);
    start *= static_cast<unsigned int>(*nval);
    unsigned int n = static_cast<unsigned int> (*nval) * nelv;
    unsigned int gn = static_cast<unsigned int>((*nelgv)*(*nval));

    // Partition element data //
    unsigned int nelv3 = static_cast<unsigned int>((*nelvin));
    unsigned int start3 = static_cast<unsigned int>(*nelb);
    unsigned int n3 = static_cast<unsigned int> (nelv);
    unsigned int gn3 = static_cast<unsigned int>((*nelgv));
    
    // Advance the reader 
    ifile=ifile+1;
    // Make a new name every time the reader is called
    std:: string fileName= "out.f0000"+ std::to_string(ifile) +".bp";
    
    // Read the new file
    startT = std::clock();
    readr = ior.Open(fileName,adios2::Mode::Read); 
    int step = 1;
    bool firstStep = true;
    // Inquire the variables
    vxr = ior.InquireVariable<float>("VX_OUT");
    vyr = ior.InquireVariable<float>("VY_OUT");
    vzr = ior.InquireVariable<float>("VZ_OUT");
    prr = ior.InquireVariable<float>("P_OUT");
    lglelr = ior.InquireVariable<int>("LGLEL_OUT");
    elapsedT = (std::clock() - startT) / (double) CLOCKS_PER_SEC;
    if (rank==0){
    std::cout <<  "rank: " << rank << " open and inquire: " << elapsedT << "s." << std::endl;
    }

    // Resize buffers    
    startT = std::clock();
    if (firstStep)
    {
        readr.LockReaderSelections();      
	vVXr.resize(n);
	vVYr.resize(n);
	vVZr.resize(n);
	vVPrr.resize(n);
	vVLGLELr.resize(n3);
	firstStep=false;
    }
    elapsedT = (std::clock() - startT) / (double) CLOCKS_PER_SEC;
    if (rank==0){
    std::cout <<  "rank: " << rank << " resize vector: " << elapsedT << "s." << std::endl;
    }

    //Select data per rank
    startT = std::clock();
    vxr.SetSelection({{start}, {n}});
    vyr.SetSelection({{start}, {n}});
    vzr.SetSelection({{start}, {n}});
    prr.SetSelection({{start}, {n}});
    lglelr.SetSelection({{start3}, {n3}});
    elapsedT = (std::clock() - startT) / (double) CLOCKS_PER_SEC;
    if (rank==0){
    std::cout <<  "rank: " << rank << " set selection: " << elapsedT << "s." << std::endl;
    }

    // Read the data
    startT = std::clock();
    readr.Get<float>(vxr,vVXr.data());
    readr.Get<float>(vyr,vVYr.data());
    readr.Get<float>(vzr,vVZr.data());
    readr.Get<float>(prr,vVPrr.data());
    readr.Get<int>(lglelr,vVLGLELr.data());
    elapsedT = (std::clock() - startT) / (double) CLOCKS_PER_SEC;
    if (rank==0){ 
    std::cout <<  "rank: " << rank << " get data: " << elapsedT << "s." << std::endl;
    }
    readr.Close();	
    
    // Copy the data into nek array
    startT = std::clock();
    for (int i=0; i<vVZr.size(); i++){
    	pr[i]=static_cast<double>(vVPrr[i]);
    	v[i]=static_cast<double>(vVXr[i]);
    	u[i]=static_cast<double>(vVYr[i]);
    	w[i]=static_cast<double>(vVZr[i]);
    }
    for (int i=0; i<vVLGLELr.size(); i++){
    	lglelrr[i]=vVLGLELr[i];
    }
    elapsedT = (std::clock() - startT) / (double) CLOCKS_PER_SEC;
    if (rank==0){  
    std::cout <<  "rank: " << rank << " copy into nek: " << elapsedT << "s." << std::endl;
    }
    startT = std::clock();

    dataTime += (std::clock() - startT) / (double) CLOCKS_PER_SEC;
}

extern "C" void adios2_finalize_(){
    //writer.Close();
    std::cout <<  "rank: " << rank << " in-situ time: " << dataTime << "s." << std::endl;
}
