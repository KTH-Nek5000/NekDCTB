# Nek5000 with in-situ adios fides? 

File changed compared to the original Nek5000:
1. core/drive.f <- the main function with split MPI communicator
2. core/drive1.f <- deactivated prepost function
3. core/3rd_party/nek_in_situ.f <- additional in-situ functions
4. core/3rd_party/nek_adios2_fides.cpp + adios2_fides.f <- additional files for adios2_fides. 
5. core/makefile.template + makenek.inc <- modified compile files 