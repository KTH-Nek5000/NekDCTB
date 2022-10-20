#include <vtkIdTypeArray.h>
#include <vtkUnsignedCharArray.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkUnstructuredGrid.h>
#include <vtkCellArray.h>
#include "vtkDoubleArray.h"
#include "vtkCPPythonAdaptorAPI.h"
#include "vtkCPDataDescription.h"
#include "vtkCPInputDataDescription.h"
#include "vtkCPProcessor.h"


vtkUnstructuredGrid *grid = NULL;

extern "C" void testcoproc_(int *flag)
{
  if (vtkCPPythonAdaptorAPI::GetCoProcessorData()) {
    *flag = 1;
  }
  else {
    *flag = 0;
  }
}

extern "C" void creategrid_(const double *x, const double *y, const double *z,
			    const int *lx1, const int *ly1, const int *lz1, 
			    const int *lelt, const int *dim, const int *flag)
{

  if (!vtkCPPythonAdaptorAPI::GetCoProcessorData()) {
    vtkGenericWarningMacro("Unable to access CoProcessorData.");
    return;
  }


  const vtkIdType points_per_el = (vtkIdType) *lx1 *
    (vtkIdType) * ly1 * (vtkIdType) *lz1;

#ifdef AMR
  if (grid != NULL) {
    if ((*flag) == 0) {
      return;
    }
    else {
      vtkCPPythonAdaptorAPI::GetCoProcessorData()->GetInputDescriptionByName("input")->Reset();
      grid->Delete();
    }
  }
#else
  if (grid != NULL)
    return;
#endif

  grid = vtkUnstructuredGrid::New();
  vtkPoints *points = vtkPoints::New();
  
  for (uint i = 0; i < points_per_el * (*lelt); i++) 
    points->InsertNextPoint((float) x[i], (float) y[i], (float) z[i]);
  grid->SetPoints(points);
  points->Delete();

  const uint h_p_lelt = (((*lx1) - 1) * 
			 ((*ly1) - 1) * ((*dim) == 3 ? ((*lz1 - 1)) : 1));
  const vtkIdType h_lelt = h_p_lelt * (*lelt);  
  const vtkIdType total = ((*dim) == 3 ? (9 * h_lelt) : (5 * h_lelt));
  
  vtkIdTypeArray *id_list = vtkIdTypeArray::New();
  id_list->SetNumberOfValues(total);
  vtkIdType *idp  = id_list->GetPointer(0);

  vtkUnsignedCharArray *cell_type = vtkUnsignedCharArray::New();
  cell_type->SetNumberOfValues(h_lelt);
  unsigned char *ctp = cell_type->GetPointer(0);

  vtkIdTypeArray *cell_location = vtkIdTypeArray::New();
  cell_location->SetNumberOfValues(h_lelt);
  vtkIdType *clp = cell_location->GetPointer(0);

  vtkIdType hex_idx = 0;
  vtkIdType el_idx = 0;
  for (uint i = 0; i < (*lelt); i++, el_idx++) {
    vtkIdType points_offset = points_per_el * el_idx;
    for (uint ii = 0; ii < ((*lx1) - 1); ii++) {
      for (uint jj = 0; jj < ((*ly1) - 1); jj++) {
	if ((*dim) == 2) {
	  *idp++ = 4;
	  *idp++ = jj * (*lx1) + ii + points_offset;
	  *idp++ = jj * (*lx1) + (ii + 1) + points_offset;
	  *idp++ = (jj + 1) * (*lx1) + (ii + 1) + points_offset;
	  *idp++ = (jj + 1) * (*lx1) + (ii) + points_offset;
	  *ctp++ = VTK_QUAD;
	  *clp++ = 5 * (hex_idx++);
	}
	else {
	  for (uint kk = 0; kk < (*lz1 - 1); kk++) {
	    *idp++ = 8;
	    *idp++ = kk * ((*ly1) * (*lx1)) + jj * (*lx1) + ii + points_offset;
	    *idp++ = kk * ((*ly1) * (*lx1)) + jj * (*lx1) + (ii + 1) + points_offset;
	    *idp++ = kk * ((*ly1) * (*lx1)) + (jj + 1) * (*lx1) + (ii + 1) + points_offset;
	    *idp++ = kk * ((*ly1) * (*lx1)) + (jj + 1) * (*lx1) + (ii) + points_offset;
	    *idp++ = (kk + 1) * ((*ly1) * (*lx1)) + jj * (*lx1) + ii + points_offset;
	    *idp++ = (kk + 1) * ((*ly1) * (*lx1)) + jj * (*lx1) + (ii + 1) + points_offset;
	    *idp++ = (kk + 1) * ((*ly1) * (*lx1)) + (jj + 1) * (*lx1) + (ii + 1) + points_offset;
	    *idp++ = (kk + 1) * ((*ly1) * (*lx1)) + (jj + 1) * (*lx1) + (ii) + points_offset;
	    *ctp++ = VTK_HEXAHEDRON;
	    *clp++ = 9 * (hex_idx++);
	  }
	}	  
      }
    }
  }

  vtkCellArray *cells = vtkCellArray::New();
  cells->SetCells(h_lelt, id_list);
  id_list->Delete();

  grid->SetCells(cell_type, cell_location, cells);
  cell_type->Delete();
  cell_location->Delete();
  cells->Delete();

  vtkCPPythonAdaptorAPI::GetCoProcessorData()->GetInputDescriptionByName("input")->SetGrid(grid);

}


extern "C" void add_scalar_field_(double *data, char *name) 
{
  vtkCPInputDataDescription* id =
    vtkCPPythonAdaptorAPI::GetCoProcessorData()->GetInputDescriptionByName("input");

  vtkUnstructuredGrid *ugrid = vtkUnstructuredGrid::SafeDownCast(id->GetGrid());
  
  if (!ugrid) {
    vtkGenericWarningMacro("No grid defined");
    return;
  }  

  if (id->IsFieldNeeded(name, vtkDataObject::POINT)) {
    vtkDoubleArray *field = vtkDoubleArray::New();
    field->SetName(name);
    field->SetArray(data, ugrid->GetNumberOfPoints(), 1);
    ugrid->GetPointData()->AddArray(field);
    field->Delete();
  }
}


extern "C" void add_vector_field_(double *xdata, double *ydata, double *zdata, 
				  int *dim, char *name) 
{
  vtkCPInputDataDescription* id =
    vtkCPPythonAdaptorAPI::GetCoProcessorData()->GetInputDescriptionByName("input");

  vtkUnstructuredGrid *ugrid = vtkUnstructuredGrid::SafeDownCast(id->GetGrid());
  
  if (!ugrid) {
    vtkGenericWarningMacro("No grid defined");
    return;
  }  

  if (id->IsFieldNeeded(name, vtkDataObject::POINT)) {
    vtkDoubleArray *field = vtkDoubleArray::New();
    field->SetName(name);
    field->SetNumberOfComponents(3);
    const vtkIdType n = ugrid->GetNumberOfPoints();
    field->SetNumberOfTuples(n);
    // TODO avoid a deep copy here
    for (vtkIdType i = 0; i < n; i++) {
      field->SetTuple3(i, xdata[i], ydata[i], ((*dim) == 3 ? zdata[i] : 0.0));
    }
    ugrid->GetPointData()->AddArray(field);
    field->Delete();
  }
}
