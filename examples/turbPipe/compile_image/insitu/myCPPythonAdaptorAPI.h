/*=========================================================================

  Program:   ParaView
  Module:    myCPPythonAdaptorAPI.h

  Copyright (c) Kitware, Inc.
  All rights reserved.
  See Copyright.txt or http://www.paraview.org/HTML/Copyright.html for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#ifndef myCPPythonAdaptorAPI_h
#define myCPPythonAdaptorAPI_h
#include "mpi.h"
#include <vtkMPI.h>
#include "vtkCPPythonAdaptorAPI.h"
 
/// Similar to vtkCPAdaptorAPI provides the implementation for API exposed to
/// typical adaptor, such as C, Fortran, except that is adds the ability to
/// initialize the coprocessor with Python capabilities.
class myCPPythonAdaptorAPI : public vtkCPPythonAdaptorAPI
{
public:
  vtkTypeMacro(myCPPythonAdaptorAPI, vtkCPPythonAdaptorAPI);

  /// Call at the start of the simulation. Users can still call
  /// CoProcessorInitialize() without arguments, in which case Python
  /// interpretor will not be initialized and hence unavailable.
  static void CoProcessorInitialize(MPI_Comm* handle, const char* pythonFileName);
  static void CoProcessorFinalize();

protected:
  myCPPythonAdaptorAPI();
  ~myCPPythonAdaptorAPI();
  static vtkMPICommunicatorOpaqueComm* Comm;
  
private:
  myCPPythonAdaptorAPI(const myCPPythonAdaptorAPI&) = delete;
  void operator=(const myCPPythonAdaptorAPI&) = delete;
};

#endif
// HeaderTest-Exclude: myCPPythonAdaptorAPI.h

