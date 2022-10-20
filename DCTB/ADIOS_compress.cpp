#include <adios2.h>
#include <mpi.h>
#include <string>
#include <vector>
#include <iostream>
#include <algorithm>

int main(int argc, char *argv[]){
    int world_rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    int colour = 1;
    MPI_Comm newcomm;
    MPI_Comm_split(MPI_COMM_WORLD, colour, world_rank, &newcomm);
    int rank, size;
    MPI_Comm_rank(newcomm, &rank);
    MPI_Comm_size(newcomm, &size);
    unsigned int urank, usize;
    urank = static_cast<unsigned int>(rank);
    usize = static_cast<unsigned int>(size);

    std::string configFile="config/config.xml";
    adios2::ADIOS adios(configFile, newcomm);
    // Reader opens insituMPI engine, for data transit
    adios2::IO inIO = adios.DeclareIO("reader");
    // The writer uses the BP engine for compression
    adios2::IO outIO = adios.DeclareIO("writer2");
    // The reader opens the global array where the writer
    // in the nek executable writes.
    adios2::Engine reader = inIO.Open("globalArray", adios2::Mode::Read, newcomm);
    adios2::Engine writer;
    // These are the variables read from global array
    adios2::Variable<double> p;
    adios2::Variable<double> vx;
    adios2::Variable<double> vy;
    adios2::Variable<double> vz;
    adios2::Variable<int> lglelw;
    // adios2::Variable<double> t;
    
    // These are the variables written
    adios2::Variable<double> p_out;
    adios2::Variable<double> vx_out;
    adios2::Variable<double> vy_out;
    adios2::Variable<double> vz_out;
    adios2::Variable<int> lglelw_out;
    // adios2::Variable<double> t_out;
    
    // These are some temporal variables (maybe?)
    std::vector<double> vP;
    std::vector<double> vVX;
    std::vector<double> vVY;
    std::vector<double> vVZ;
    std::vector<int> vVLGLELW;
    // std::vector<double> vT;
    
    // Define variables to know starts and counts, etc.
    std::size_t total_size, my_start, my_count;
    std::size_t total_size2, my_start2, my_count2;
    std::size_t total_size3, my_start3, my_count3;
    bool firststep = true;

    unsigned int step = 1;
    try
    {
    while (true)
    {
        adios2::StepStatus status = reader.BeginStep();
        if (status != adios2::StepStatus::OK)
        {
            break;
        }
	// Check that the variables are in global array

        p = inIO.InquireVariable<double>("P");
        if (!p)
        {
            std::cout << "Error: NO variable P found. Unable to proceed. "
                            "Exiting. "
                        << std::endl;
            break;
        }
        vx = inIO.InquireVariable<double>("VX");
        if (!vx)
        {
            std::cout << "Error: NO variable VX found. Unable to proceed. "
                            "Exiting. "
                        << std::endl;
            break;
        }
        vy = inIO.InquireVariable<double>("VY");
        if (!vy)
        {
            std::cout << "Error: NO variable VY found. Unable to proceed. "
                            "Exiting. "
                        << std::endl;
            break;
        }
        vz = inIO.InquireVariable<double>("VZ");
        if (!vz)
        {
            std::cout << "Error: NO variable VZ found. Unable to proceed. "
                            "Exiting. "
                        << std::endl;
            break;
        }
        lglelw = inIO.InquireVariable<int>("LGLEL");
        if (!lglelw)
        {
            std::cout << "Error: NO variable LGLEL found. Unable to proceed. "
                            "Exiting. "
                        << std::endl;
            break;
        }
        // t = inIO.InquireVariable<double>("T");
        // if (!t)
        // {
        //     std::cout << "Error: NO variable T found. Unable to proceed. "
        //                     "Exiting. "
        //                 << std::endl;
        //     break;
        // }
	
	// Set all the vector sizes the first time you write
        if(firststep){
            total_size = vz.Shape()[0];
            my_count = total_size / usize;
            my_start = my_count * urank;
            if(total_size % usize != 0){
                if(urank < (total_size%usize)){
                    ++my_count;
                    my_start += urank;
                }else{
                    my_start += (total_size%usize);
                }
            }
            vP.resize(my_count);
            vVX.resize(my_count);
            vVY.resize(my_count);
            vVZ.resize(my_count);
            // total_size2 = t.Shape()[0];
            // my_start2 = (total_size2/size) * rank;
            // my_count2 = (total_size2/size);
            // if(total_size2 % usize != 0){
            //     if(urank < (total_size2%usize)){
            //         ++my_count2;
            //         my_start2 += urank;
            //     }else{
            //         my_start2 += (total_size2%usize);
            //     }
            // }
            // vT.resize(my_count2);
            total_size3 = lglelw.Shape()[0];
            my_start3 = (total_size3/size) * rank;
            my_count3 = (total_size3/size);
            if(total_size3 % usize != 0){
                if(urank < (total_size3%usize)){
                    ++my_count3;
                    my_start3 += urank;
                }else{
                    my_start3 += (total_size3%usize);
                }
            }
            vVLGLELW.resize(my_count3);
            firststep = false;
	    //std::cout << rank << ": " << total_size << ", " << my_start << ", " << my_count << std::endl;
            
	    // Define the variable where the writing will be done. 
	    // The total size, where in the vector to start writing
	    // and how many entries in the vector will be written.
	    p_out = outIO.DefineVariable<double>("P_OUT", {total_size}, {my_start}, {my_count});
            vx_out = outIO.DefineVariable<double>("VX_OUT", {total_size}, {my_start}, {my_count});
            vy_out = outIO.DefineVariable<double>("VY_OUT", {total_size}, {my_start}, {my_count});
            vz_out = outIO.DefineVariable<double>("VZ_OUT", {total_size}, {my_start}, {my_count});
            lglelw_out = outIO.DefineVariable<int>("LGLEL_OUT", {total_size3}, {my_start3}, {my_count3});
            // t_out = outIO.DefineVariable<double>("T_OUT", {total_size2}, {my_start2}, {my_count2});
        }

	// Set how much of the vector in the global array will be read
        p.SetSelection({{my_start}, {my_count}});
        vx.SetSelection({{my_start}, {my_count}});
        vy.SetSelection({{my_start}, {my_count}});
        vz.SetSelection({{my_start}, {my_count}});
        // t.SetSelection({{my_start2}, {my_count2}});
        lglelw.SetSelection({{my_start3}, {my_count3}});
	
	// Read from vz into vVZ
        reader.Get<double>(p, vP.data());
        reader.Get<double>(vx, vVX.data());
        reader.Get<double>(vy, vVY.data());
        reader.Get<double>(vz, vVZ.data());
        // reader.Get<double>(t, vT.data());
        reader.Get<int>(lglelw, vVLGLELW.data());
        reader.EndStep();

	// Write the outputs
	std:: string fileName= "out.f0000"+ std::to_string(step) +".bp";
        //writer = outIO.Open("out.bp", adios2::Mode::Write, newcomm);
        writer = outIO.Open(fileName, adios2::Mode::Write, newcomm);
        writer.BeginStep();
        writer.Put<double>(p_out, vP.data());
        writer.Put<double>(vx_out, vVX.data());
        writer.Put<double>(vy_out, vVY.data());
        writer.Put<double>(vz_out,  vVZ.data());
        //writer.Put<double>(t_out, vT.data());
        writer.Put<int>(lglelw_out,  vVLGLELW.data());
        writer.EndStep();
        writer.Close();
        ++step;
    }
    reader.Close();
    }
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
    MPI_Finalize();
    return 0;
}
