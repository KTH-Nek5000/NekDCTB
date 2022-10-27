#include <adios2.h>
#include <mpi.h>
#include <string>
#include <vector>
#include <iostream>
#include <time.h>
#include <algorithm>
#include "nek_catalyst.h"
#include "myCPPythonAdaptorAPI.h"

int world_rank, world_size;
int rank, size;
int nek_size;
int lx1,ly1,lz1,lx2,ly2,lz2, lxyz1,lxyz2;
int nelv, nelt, nelgt, nelgv;
int nelbv, nelbt;
int dim;
int i, li, lj, lk, idx1, idx2;
int dx12, dy12, dz12;
int iostep;
int step;
int flag;
double simulation_start, simulation_dt, simulation_time;
char nameP[] = "pressure";
char nameV[] = "velocity";
char nameT[] = "temperature";
std::clock_t startT;
std::clock_t startTotal;
double setupTime, insituTime, totalTime;

MPI_Comm newcomm;
std::size_t total, count, start;
std::size_t total2, count2, start2;
std::size_t total3, count3, start3;
std::size_t total4, count4, start4;

adios2::Engine reader;

std::vector<int> vINT_CONST;
std::vector<double> vDOUBLE_CONST;
std::vector<double> vXM1;
std::vector<double> vYM1;
std::vector<double> vZM1;

std::vector<double> vVX;
std::vector<double> vVY;
std::vector<double> vVZ;
std::vector<double> vP;
std::vector<double> vT;
std::vector<double> vBM1;


// These are the variables read from global array
adios2::Variable<double> p;
adios2::Variable<double> vx;
adios2::Variable<double> vy;
adios2::Variable<double> vz;
adios2::Variable<double> t;
adios2::Variable<double> bm1;

adios2::Variable<int> init_int_vec1;
adios2::Variable<double> init_double_vec1;
adios2::Variable<double> xm1;
adios2::Variable<double> ym1;
adios2::Variable<double> zm1;

adios2::IO inIO;
bool firstStep;
int main(int argc, char *argv[]){
    MPI_Init(&argc, &argv);
    startTotal = std::clock();
    startT = std::clock();
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    int colour = 1;
    MPI_Comm_split(MPI_COMM_WORLD, colour, world_rank, &newcomm);
    MPI_Comm_rank(newcomm, &rank);
    MPI_Comm_size(newcomm, &size);
    nek_size=world_size-size;
    if(nek_size % size != 0){
    	if(!rank) std::cout << "The number of cores to run in-situ analysis should be a factor of the number of cores to run simulation!" << std::endl;
	MPI_Finalize();
        return 9;
    }
    std::string configFile="config/config.xml";
    adios2::ADIOS adios(configFile, newcomm);
    // Reader opens insituMPI engine, for data transit
    inIO = adios.DeclareIO("reader");
    // The reader opens the global array where the writer
    // in the nek executable writes.
    reader = inIO.Open("globalArray", adios2::Mode::Read, newcomm);

    // Define variables to know starts and counts, etc.
    adios2::StepStatus status = reader.BeginStep();
    if (status != adios2::StepStatus::OK){
        reader.EndStep();
        reader.Close();
	std::cout << "Reader has not done begin step." << std::endl;
        MPI_Finalize();
        return 3;
    }
    std::size_t int_vec1_size, int_vec1_start, int_vec1_count;
    std::size_t double_vec1_size, double_vec1_start, double_vec1_count;
    int len_nek = 12;
    init_int_vec1 = inIO.InquireVariable<int>("INT_CONST");
    if (!init_int_vec1){
        std::cout << "Error: NO variable INT_CONST found. Unable to proceed. \n Exiting. " << std::endl;
        return 4;
    }
    int_vec1_size = init_int_vec1.Shape()[0];
    if(int_vec1_size != static_cast<std::size_t>(len_nek * nek_size)){
        std::cout << int_vec1_size << " vs " << len_nek * nek_size << " \n Error: INT_CONST length not match" << std::endl;
        return 4;
    }
    /*
    if(!rank){
        int_vec1_count = static_cast<std::size_t>(len_nek);
	int_vec1_start = 0;
    }else{
        int_vec1_count = static_cast<std::size_t>((static_cast<int>(int_vec1_size) - len_nek) / (size-1));
        int_vec1_start = static_cast<std::size_t>(rank - 1) * int_vec1_count + static_cast<std::size_t>(len_nek);
    }
    */
    int_vec1_count = static_cast<std::size_t>(len_nek * nek_size / size);
    int_vec1_start = static_cast<std::size_t>(rank) * int_vec1_count;
    init_int_vec1.SetSelection({{int_vec1_start}, {int_vec1_count}});
    vINT_CONST.resize(int_vec1_count);
    reader.Get<int>(init_int_vec1,vINT_CONST.data());
    
    int len_double_const_nek = 2;
    init_double_vec1 = inIO.InquireVariable<double>("DOUBLE_CONST");
    if (!init_double_vec1){
        std::cout << "Error: NO variable DOUBLE_CONST found. Unable to proceed. \n Exiting. " << std::endl;
        return 4;
    }
    double_vec1_size = init_double_vec1.Shape()[0];
    if(double_vec1_size != static_cast<std::size_t>(len_double_const_nek * nek_size)){
        std::cout << double_vec1_size << " vs " << len_double_const_nek * nek_size << " \n Error: DOUBLE_CONST length not match" << std::endl;
        return 4;
    }
    /*
    if(!rank){
        double_vec1_count = static_cast<std::size_t>(len_double_const_nek);
        double_vec1_start = 0;
    }else{
        double_vec1_count = static_cast<std::size_t>((static_cast<int>(double_vec1_size) - len_double_const_nek) / (size-1));
        double_vec1_start = static_cast<std::size_t>(rank-1) * double_vec1_count + static_cast<std::size_t>(len_double_const_nek);
    }
    */
    double_vec1_count = static_cast<std::size_t>(len_double_const_nek * nek_size / size);
    double_vec1_start = static_cast<std::size_t>(rank) * double_vec1_count;
    vDOUBLE_CONST.resize(double_vec1_count);
    init_double_vec1.SetSelection({{double_vec1_start}, {double_vec1_count}});
    reader.Get<double>(init_double_vec1,vDOUBLE_CONST.data());
    reader.EndStep();
    reader.BeginStep();
    //reader.PerformGets();
    //reader.EndStep();
    
    /*
    std::vector<int>vINT_CONST2(len_nek);
    if(!rank) for(i=0;i<len_nek;++i){vINT_CONST2[i]=vINT_CONST[i];}
    std::vector<double>vDOUBLE_CONST2(len_double_const_nek);
    if(!rank) for(i=0;i<len_double_const_nek;++i){vDOUBLE_CONST2[i]=vDOUBLE_CONST[i];}
    //if(!rank)std::cout << "Before bcast" << std::endl;
    MPI_Bcast(vINT_CONST2.data(), len_nek, MPI_INT, 0, newcomm);
    MPI_Bcast(vDOUBLE_CONST2.data(), len_double_const_nek, MPI_DOUBLE, 0, newcomm);
    //if(!rank)std::cout << "After bcast" << std::endl;

    lx1 = vINT_CONST2[0];
    ly1 = vINT_CONST2[1];
    lz1 = vINT_CONST2[2];
    lx2 = vINT_CONST2[3];
    ly2 = vINT_CONST2[4];
    lz2 = vINT_CONST2[5];
    nelgv = vINT_CONST2[6];
    nelgt = vINT_CONST2[7];
    iostep = vINT_CONST2[8];
    dim = 3;
    if(!vINT_CONST2[9]) dim = 2;
    */
    lx1 = vINT_CONST[0];
    ly1 = vINT_CONST[1];
    lz1 = vINT_CONST[2];
    lx2 = vINT_CONST[3];
    ly2 = vINT_CONST[4];
    lz2 = vINT_CONST[5];
    nelgv = vINT_CONST[8];
    nelgt = vINT_CONST[9];
    iostep = vINT_CONST[10];
    dim = 3;
    if(!vINT_CONST[11]) dim = 2;
    nelv = 0;
    nelt = 0;
    for(i = 0; i * size < nek_size; ++i){
    	nelv += vINT_CONST[i * len_nek + 6];
	nelt += vINT_CONST[i * len_nek + 7];
    }
    lxyz1=lx1*ly1*lz1;
    lxyz2=lx2*ly2*lz2;
    dx12 = static_cast<int>((lx1 - lx2)/2);
    dy12 = static_cast<int>((ly1 - ly2)/2);
    dz12 = static_cast<int>((lz1 - lz2)/2);
    if(!rank) std::cout << "dim: " << dim << std::endl;
    simulation_start = vDOUBLE_CONST[0];
    simulation_dt = vDOUBLE_CONST[1];
    /*
    nelv = static_cast<int>(nelgv / size);
    if(rank < nelgv % size) ++nelv;
    nelt = static_cast<int>(nelgt / size);
    if(rank < nelgt % size) ++nelt;
    */
    step = 0;
    std::vector<int>temp(size);
    std::vector<int>temp1(size);
    
    MPI_Allgather(&nelv, 1, MPI_INT, temp.data(), 1, MPI_INT, newcomm);
    MPI_Allgather(&nelt, 1, MPI_INT, temp1.data(), 1, MPI_INT, newcomm);
    int check_nelgv = 0;
    int check_nelgt = 0;
    nelbt = 0;
    nelbv = 0;
    for(i=0;i<rank;++i){
        nelbv += temp[i];
        nelbt += temp1[i];
	check_nelgv += temp[i];
	check_nelgt += temp1[i];
    }
    for(i=rank;i<size;++i){
	check_nelgv += temp[i];
        check_nelgt += temp1[i];
    }
    if(check_nelgv != nelgv || check_nelgt != nelgt){
    	if(!rank) std::cout << "Elements are missing!" << std::endl;
	MPI_Finalize();
        return 8;
    }
    std::cout << rank << " : " << nelgv << " , " << nelbv << " , " << nelv << std::endl;
    /*
    status = reader.BeginStep();
    if (status != adios2::StepStatus::OK){
        reader.EndStep();
        reader.Close();
        MPI_Finalize();
        return 3;
    }
    */
    //total3 = static_cast<std::size_t>(lxyz1 * nelgt);
    xm1 = inIO.InquireVariable<double>("X");
    if(!xm1){
        std::cout << "Error: NO variable X found. Unable to proceed. \n Exiting. " << std::endl;
        return 4;
    }
    total3 = xm1.Shape()[0];
    ym1 = inIO.InquireVariable<double>("Y");
    if(!ym1){
        std::cout << "Error: NO variable Y found. Unable to proceed. \n Exiting. " << std::endl;
        return 4;
    }
    zm1 = inIO.InquireVariable<double>("Z");
    if(!zm1){
        std::cout << "Error: NO variable Z found. Unable to proceed. \n Exiting. " << std::endl;
        return 4;
    }
    if(total3 != static_cast<std::size_t>(lxyz1 * nelgt)){
        std::cout << total3  << " vs " << static_cast<std::size_t>(lxyz1 * nelgt) << " \n Error: INT_CONST length not match" << std::endl;
        return 4;
    } 
    start3 = static_cast<std::size_t>(lxyz1 * nelbt);
    count3 = static_cast<std::size_t>(lxyz1 * nelt);
    std::cout << rank << " : " << total3 << " , " << start3 << " , " << count3 << std::endl;
    xm1.SetSelection({{start3},{count3}});
    ym1.SetSelection({{start3},{count3}});
    zm1.SetSelection({{start3},{count3}});
    vXM1.resize(count3);
    vYM1.resize(count3);
    vZM1.resize(count3);
    reader.Get<double>(xm1,vXM1.data());
    reader.Get<double>(ym1,vYM1.data());
    reader.Get<double>(zm1,vZM1.data());
    firstStep=true;
    if(!rank)std::cout << "1st reader step" << std::endl;
    /*
    status = reader.BeginStep();
    if (status != adios2::StepStatus::OK){
        reader.EndStep();
        reader.Close();
        MPI_Finalize();
        return 3;
    }
    */
    

    myCPPythonAdaptorAPI::CoProcessorInitialize(&newcomm,NULL);
    catalyst_usrpipe();
   
    p = inIO.InquireVariable<double>("P");
    if (!p){
        std::cout << "Error: NO variable P found. Unable to proceed. \n Exiting. " << std::endl;
        return 4;
    }
    vx = inIO.InquireVariable<double>("VX");
    if (!vx){
        std::cout << "Error: NO variable VX found. Unable to proceed. \n Exiting. " << std::endl;
        return 4;
    }
    vy = inIO.InquireVariable<double>("VY");
    if (!vy){
        std::cout << "Error: NO variable VY found. Unable to proceed. \n Exiting. " << std::endl;
        return 4;
    }
    vz = inIO.InquireVariable<double>("VZ");
    if (!vz){
        std::cout << "Error: NO variable VZ found. Unable to proceed. \n Exiting. " << std::endl;
        return 4;
    }
    t = inIO.InquireVariable<double>("T");
    if (!t){
        std::cout << "Error: NO variable T found. Unable to proceed. \n Exiting. " << std::endl;
        return 4;
    }
    bm1 = inIO.InquireVariable<double>("BM1");
    if (!bm1){
        std::cout << "Error: NO variable BM1 found. Unable to proceed. \n Exiting. " << std::endl;
        return 4;
    }
    total = vx.Shape()[0];
    if(total != static_cast<std::size_t>(lxyz1*nelgv)){
         std::cout << total << " vs " << lxyz1*nelgv << " \n Error: V length not match" << std::endl;
        return 4;
    }
    start = static_cast<std::size_t>(lxyz1 * nelbv);
    count = static_cast<std::size_t>(lxyz1 * nelv);

    total2 = p.Shape()[0];
    if(total2 != static_cast<std::size_t>(lxyz2*nelgv)){
         std::cout << total2 << " vs " << lxyz2*nelgv << " \n Error: P length not match" << std::endl;
        return 4;
    }
    start2 = static_cast<std::size_t>(lxyz2 * nelbv);
    count2 = static_cast<std::size_t>(lxyz2 * nelv);

    vP.resize(count2);
    vVX.resize(count);
    vVY.resize(count);
    vVZ.resize(count);
    vT.resize(count3);
    vBM1.resize(count3);

    p.SetSelection({{start2}, {count2}});
    vx.SetSelection({{start}, {count}});
    vy.SetSelection({{start}, {count}});
    vz.SetSelection({{start}, {count}});
    t.SetSelection({{start3}, {count3}});
    bm1.SetSelection({{start3}, {count3}});

    reader.Get<double>(p, vP.data());
    reader.Get<double>(vx, vVX.data());
    reader.Get<double>(vy, vVY.data());
    reader.Get<double>(vz, vVZ.data());
    reader.Get<double>(t, vT.data());
    reader.Get<double>(bm1, vBM1.data());
    if(!rank)std::cout << "2nd reader step" << std::endl;
    reader.EndStep();
    //myCPPythonAdaptorAPI::CoProcessorInitialize(&newcomm,NULL);
    //catalyst_usrpipe();
    setupTime = (std::clock() - startT) / (double) CLOCKS_PER_SEC;
    insituTime = 0.0;
    step = iostep;
    if(!rank) std::cout << "Initialization done." << std::endl;
    std::vector<double>vP2(count, 0.0);
    try{
    while(true){
        status = reader.BeginStep();
        if(status == adios2::StepStatus::OK){
	    p = inIO.InquireVariable<double>("P");
    	    if (!p){
        	std::cout << "Error: NO variable P found. Unable to proceed. \n Exiting. " << std::endl;
        	return 4;
    	    }
    	    vx = inIO.InquireVariable<double>("VX");
    	    if (!vx){
        	std::cout << "Error: NO variable VX found. Unable to proceed. \n Exiting. " << std::endl;
        	return 4;
     	    }
    	    vy = inIO.InquireVariable<double>("VY");
    	    if (!vy){
        	std::cout << "Error: NO variable VY found. Unable to proceed. \n Exiting. " << std::endl;
        	return 4;
    	    }
    	    vz = inIO.InquireVariable<double>("VZ");
    	    if (!vz){
        	std::cout << "Error: NO variable VZ found. Unable to proceed. \n Exiting. " << std::endl;
        	return 4;
    	    }
    	    t = inIO.InquireVariable<double>("T");
    	    if (!t){
        	std::cout << "Error: NO variable T found. Unable to proceed. \n Exiting. " << std::endl;
        	return 4;
    	    }
    	    bm1 = inIO.InquireVariable<double>("BM1");
    	    if (!bm1){
        	std::cout << "Error: NO variable BM1 found. Unable to proceed. \n Exiting. " << std::endl;
        	return 4;
     	    }
	    /*
	    if(firstStep){
	    	xm1 = inIO.InquireVariable<double>("X");
    		if(!xm1){
        	    std::cout << "Error: NO variable X found. Unable to proceed. \n Exiting. " << std::endl;
        	    return 4;
    		}
    		ym1 = inIO.InquireVariable<double>("Y");
    		if(!ym1){
        	     std::cout << "Error: NO variable Y found. Unable to proceed. \n Exiting. " << std::endl;
        	     return 4;
    		}
    		zm1 = inIO.InquireVariable<double>("Z");
    		if(!zm1){
        	    std::cout << "Error: NO variable Z found. Unable to proceed. \n Exiting. " << std::endl;
      		    return 4;
   		}
	        	
    	    	xm1.SetSelection({{start3},{count3}});
    		ym1.SetSelection({{start3},{count3}});
    		zm1.SetSelection({{start3},{count3}});
		
    		reader.Get<double>(xm1,vXM1.data(),adios2::Mode::Sync);
    		reader.Get<double>(ym1,vYM1.data(),adios2::Mode::Sync);
    		reader.Get<double>(zm1,vZM1.data(),adios2::Mode::Sync);
		firstStep=false;
	    }
	    */
	    
            p.SetSelection({{start2}, {count2}});
            vx.SetSelection({{start}, {count}});
            vy.SetSelection({{start}, {count}});
            vz.SetSelection({{start}, {count}});
            t.SetSelection({{start3}, {count3}});
            bm1.SetSelection({{start3}, {count3}});
	    
            reader.Get<double>(p, vP.data());
            reader.Get<double>(vx, vVX.data());
            reader.Get<double>(vy, vVY.data());
            reader.Get<double>(vz, vVZ.data());
            reader.Get<double>(t, vT.data());
            reader.Get<double>(bm1, vBM1.data());
            reader.EndStep();
            startT = std::clock();
	    for(i=0;i<nelv;++i){
	    	for(li=0;li<lx2;++li){
		    for(lj=0;lj<ly2;++lj){
		    	for(lk=0;lk<lz2;++lk){
			    idx1=i*lxyz2+lz2*ly2*li+lz2*lj+lk;
			    idx2=i*lxyz1+lz1*ly1*(li+dx12)+lz1*(lj+dy12)+lk+dz12;
			    vP2[idx2]=vP[idx1];
			}
		    }
		}
	    }
            simulation_time = simulation_start + static_cast<double>(step) * simulation_dt;
            requestdatadescription(&step, &simulation_time, &flag);
            if(flag){
                needtocreategrid(&flag);
                creategrid(vXM1.data(), vYM1.data(), vZM1.data(), &lx1, &ly1, &lz1, &nelt, &dim);
                add_scalar_field(vP2.data(), nameP);
                add_vector_field(vVX.data(), vVY.data(), vVZ.data(), &dim, nameV);
                add_scalar_field(vT.data(), nameT);
                coprocess();
            }
	    step += iostep;
            insituTime += (std::clock() - startT) / (double) CLOCKS_PER_SEC;
        }else if(status == adios2::StepStatus::EndOfStream){
            if(!rank) std::cout << "End of stream" << std::endl;
            reader.Close();
            coprocessorfinalize();
            break;
        }
    }}
    catch (std::invalid_argument &e)
    {
        std::cout << "Invalid argument exception, STOPPING PROGRAM from rank "
                  << rank << "\n";
        std::cout << e.what() << "\n";
    }
    catch (std::ios_base::failure &e)
    {
        std::cout << "IO System base failure exception, STOPPING PROGRAM "
                     "from rank "
                  << rank << "\n";
        std::cout << e.what() << "\n";
    }
    catch (std::exception &e)
    {
        std::cout << "Exception, STOPPING PROGRAM from rank " << rank << "\n";
        std::cout << e.what() << "\n";
    }
    totalTime = (std::clock() - startTotal) / (double) CLOCKS_PER_SEC;
    std::cout << "catalystSpace rank: " << rank << " in-situ time: " << insituTime << " s, total time: " << totalTime << " s." << std::endl;   
    MPI_Finalize();
    return 0;    
}

