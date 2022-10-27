/*=========================================================================

  Program:   ParaView
  Module:    $RCSfile$

  Copyright (c) Kitware, Inc.
  All rights reserved.
  See Copyright.txt or http://www.paraview.org/HTML/Copyright.html for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "myCPPythonAdaptorAPI.h"

#include "vtkCPProcessor.h"
#include "vtkCPPythonScriptPipeline.h"

//----------------------------------------------------------------------------
vtkMPICommunicatorOpaqueComm* myCPPythonAdaptorAPI::Comm = NULL;
void myCPPythonAdaptorAPI::CoProcessorInitialize(MPI_Comm* handle, const char* pythonFileName)
{
  if (!Superclass::CoProcessor)
  {
    Superclass::CoProcessor = vtkCPProcessor::New();
    myCPPythonAdaptorAPI::Comm = new vtkMPICommunicatorOpaqueComm(handle);
    Superclass::CoProcessor->Initialize(*myCPPythonAdaptorAPI::Comm);
  }
  // needed to initialize vtkCPDataDescription.
  Superclass::CoProcessorInitialize(pythonFileName);

  if (pythonFileName)
  {
    vtkCPPythonScriptPipeline* pipeline = vtkCPPythonScriptPipeline::New();
    pipeline->Initialize(pythonFileName);
    Superclass::CoProcessor->AddPipeline(pipeline);
    pipeline->FastDelete();
  }
}

void myCPPythonAdaptorAPI::CoProcessorFinalize()
{
  if(myCPPythonAdaptorAPI::Comm){
    delete myCPPythonAdaptorAPI::Comm;
    myCPPythonAdaptorAPI::Comm = NULL;
  }
  Superclass::CoProcessorFinalize();
}
