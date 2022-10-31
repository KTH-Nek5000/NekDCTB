# NekDCTB
Data Compression Tool Box for Nek5000

NekDCTB is a data compression toolbox for in-situ and static data compression of files generated with the CFD solver Nek5000. The method used, however, is functional for any other spectral element solver. NekDCTB supports lossy and lossless compression, and uses the library ADIOS2 for the pusrpose of encoding the data and for data transint during asyncrhonous operation. The target of the toolbox are 3D CFD flow fields with billions of gridpoints that are not only a major computational challenge but also a concern for storage. The project has been built on top of other tools developed at KTH - Royal Institute of Technology in Sweden, specifically the Nek5000 KTH Framework, however core functions have been shipped with this distribution in order to provide a stand-alone toolbox.


## Documentation
The documentation and installation procedures can be found at: https://kth-nek5000.github.io/NekDCTB/index.html
