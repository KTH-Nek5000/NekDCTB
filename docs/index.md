# NekDCTB
NekDCTB is a data compression toolbox for in-situ and static data compression of files generated with the CFD solver Nek5000 [Link](https://nek5000.mcs.anl.gov/). 
The method used, however, is functional for any other spectral element solver. NekDCTB supports lossy and lossless compression, and uses the library ADIOS2
[Link](https://github.com/ornladios) for the pusrpose of encoding the data and for data transint during asyncrhonous operation. The target of the toolbox are 3D CFD
flow fields with billions of gridpoints that are not only a major computational challenge but also a concern for storage. 
