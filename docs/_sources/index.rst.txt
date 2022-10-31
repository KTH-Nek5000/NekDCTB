NekDCTB
=======

NekDCTB is a data compression toolbox for in-situ and static data
compression of files generated with the CFD solver
`Nek5000 <https://nek5000.mcs.anl.gov/>`__. The method used, however, is
functional for any other spectral element solver. NekDCTB supports lossy
and lossless compression, and uses the library
`ADIOS2 <https://github.com/ornladios>`__ for the pusrpose of encoding
the data and for data transint during asyncrhonous operation. The target
of the toolbox are 3D CFD flow fields with billions of gridpoints that
are not only a major computational challenge but also a concern for
storage. The project has been built on top of other tools developed at
KTH - Royal Institute of Technology in Sweden, specifically the `Nek5000
KTH Framework <https://github.com/KTH-Nek5000/KTH_Framework>`__, however
core functions have been shipped with this distribution in order to
provide a stand-alone toolbox.

.. toctree::
   :maxdepth: 3
   :caption: Contents:

   ./pages/installation
   ./pages/list_of_codes
   ./pages/background
   ./pages/results
   ./pages/casename_par
   ./pages/casename_usr
   ./pages/compile_script
   ./pages/makefile_usr_inc
   ./pages/truncation
   ./pages/encoding

