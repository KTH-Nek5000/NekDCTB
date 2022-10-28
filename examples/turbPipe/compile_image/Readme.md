# How to try Nek5000 with in-situ image generation? 
1. Install Paraview 5.9 with Catalyst with the instructions (3.1-3.3) on: 
https://github.com/KTH-Nek5000/InSituPackage

2. Compile the synchronous in-situ image generation:  
2.1 Change the `PARAVIEW` an `OSMESA` path in `synch_compile_script`   
2.2 Compile with `synch_compile_script`   


3. Compile the asynchronous in-situ image generation:   
3.1 Change the `PARAVIEW` an `OSMESA` path in `insitu/make.sh`   
3.2 Compile with `asynch_compile_script`

File changed compared to the original Nek5000:
1. core/drive.f <- the main function with split MPI communicator
2. core/drive1.f <- deactivated prepost function
3. core/3rd_party/nek_in_situ.f <- additional in-situ functions
4. core/3rd_party/nek_catalyst.cxx + catalyst.f <- additional files for synchronous in-situ image generation. 
5. core/3rd_party/nek_async_catalyst.cpp + async_catalyst.f <- additional files for asynchronous in-situ image generation. 
6. core/makefile.template + makenek.inc <- modified compile files 