.. _casename_par:

<casename>.par
==============

The following runtime parameters must be added to the .par file

.. code:: sh

   [_IOTRUNC]
   FILETOCOMP           = "<casename>0.f"
   SCOMP                = yes
   SDECOMP              = yes
   INSITUCOMP           = no
   NUMFILE              = 5
   TRUNCSTEP            = 5
   READSTEP             = 7
   TARGETERR            = 1e-3

The parameters control the following characteristics:

-  ``FILETOCOMP`` [string] specifies the name of the files in case that
   the targets are files that have already been written to disk.

-  ``SCOMP`` [yes/no] Activates/Deactivates one of the modes available.
   if yes, already written xxx0.fxxx files from nek will be compressed
   into out.xxx.bp files.

-  ``SDECOMP`` [yes/no] Activates/Deactivates one of the modes
   available. if yes, out.xxx.bp files will be decompressed into
   standard nek xxxx0.fxxx files.

-  ``INSITUCOMP`` [yes/no] Activates/Deactivates one of the modes
   available. if yes, nek5000 will execute normally and will produced
   compressed outputs out.xxx.bp with a frequency given by ``TRUNCSTEP``

-  ``NUMFILE`` [integer] Indicates the number of files to be compressed
   or decompressed. This only applies if ``SCOMP`` and/or ``SDECOMP``
   are set to yes.

-  ``TRUNCSTEP`` [integer] Determines the frequency (in number of time
   steps) in which files will be compressed and written when
   ``INSITUCOMP`` is set to yes.

-  ``READSTEP`` [integer] Determines the frequency (in number of time
   steps) in which files will be decompressed when ``INSITUCOMP`` is set
   to yes. This feature has no real use in a normal run but has been
   kept for easily user debugging. It is recommended that this is set to
   a number higher than the total number of time steps in a run.

-  ``TARGETERR`` [real] Determines the maximun value of the root mean
   squared error (Weighted L2 norm) between the compressed data and the
   original one. This parameter controls how much compression will be
   obtained by regulating the accepted error by the user.
