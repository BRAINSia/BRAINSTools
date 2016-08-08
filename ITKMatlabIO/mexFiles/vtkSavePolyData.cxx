#include "mex.h"
#include "matrix.h"
#include "vtkCommon.h"
#include "vtkPolyData.h"
#include "vtkPoints.h"
#include "vtkPointData.h"
#include "vtkCellArray.h"
#include "vtkDoubleArray.h"
#include "vtkDataSetAttributes.h"
#include "vtkSmartPointer.h"
#include "vtkLine.h"
#include "vtkXMLPolyDataWriter.h"
#include "vtkPolyDataWriter.h"
#include "vtksys/SystemTools.hxx"
#include "vtkErrorObserver.h"
#include "vtkVersion.h"
namespace
{
void
WritePolyData(const vtkPolyData *pd, const char *filename, bool writeBinary = false, bool writeCompressed = false )
{
  const std::string ext(vtksys::SystemTools::GetFilenameExtension(filename) );
  // all the casting is because vtk SetInput isn't const-correct.
  vtkDataObject *dataObject =
    const_cast<vtkDataObject *>(static_cast<const vtkDataObject *>(pd) );

  vtkSmartPointer<ErrorObserver> errorObserver =
    vtkSmartPointer<ErrorObserver>::New();

  if( ext == ".vtp" )
    {
    vtkSmartPointer<vtkXMLPolyDataWriter> writer =
      vtkSmartPointer<vtkXMLPolyDataWriter>::New();

    writer->AddObserver(vtkCommand::ErrorEvent, errorObserver);
    writer->AddObserver(vtkCommand::WarningEvent, errorObserver);

    if( writeBinary )
      {
      writer->SetDataModeToBinary();
      }
    else
      {
      writer->SetDataModeToAscii();
      }
    if( writeCompressed )
      {
      writer->SetCompressorTypeToZLib();
      }
    else
      {
      writer->SetCompressorTypeToNone();
      }
#if VTK_MAJOR_VERSION < 6
    writer->SetInput(dataObject);
#else
    writer->SetInputData(dataObject);
#endif
    writer->SetFileName(filename);
    writer->Write();
    if( errorObserver->GetWarning() )
      {
      std::string error("Caught warning! ");
      error += errorObserver->GetWarningMessage();
      mexWarnMsgTxt(error.c_str() );
      }

    if( errorObserver->GetError() )
      {
      std::string error("Caught error! ");
      error += errorObserver->GetErrorMessage();
      mexWarnMsgTxt(error.c_str() );
      }

    }
  else
    {
    vtkSmartPointer<vtkPolyDataWriter> writer =
      vtkSmartPointer<vtkPolyDataWriter>::New();
    writer->AddObserver(vtkCommand::ErrorEvent, errorObserver);
    writer->AddObserver(vtkCommand::WarningEvent, errorObserver);
    if( writeBinary )
      {
      writer->SetFileTypeToBinary();
      }
#if VTK_MAJOR_VERSION < 6
    writer->SetInput(dataObject);
#else
    writer->SetInputData(dataObject);
#endif
    writer->SetFileName(filename);
    writer->Write();
    if( errorObserver->GetWarning() )
      {
      std::string error("Caught warning! ");
      error += errorObserver->GetWarningMessage();
      mexWarnMsgTxt(error.c_str() );
      }

    if( errorObserver->GetError() )
      {
      std::string error("Caught error! ");
      error += errorObserver->GetErrorMessage();
      mexWarnMsgTxt(error.c_str() );
      }
    }
}

}

void vtkSavePolyData(int, mxArray *[],
                     int nrhs, const mxArray *prhs[])
{
  if( nrhs < 2 || !mxIsChar(prhs[0]) )
    {
    mexErrMsgTxt("vtkSavePolyData: requires two args: one string, one struct");
    return;
    }
  const mxArray * const filenameMx(prhs[0]);
  const int             filenameLen(mxGetM(filenameMx) * mxGetN(filenameMx) + 1);
  char *                filename = static_cast<char *>(mxCalloc(filenameLen, sizeof(mxChar) ) );
  mxGetString(filenameMx, filename, filenameLen);

  const mxArray *const structMx = prhs[1];
  const mwSize         numberOfFields(mxGetNumberOfFields(structMx) );
  if( numberOfFields == 0 )
    {
    mexErrMsgTxt("No data in structure array");
    return;
    }
  const mxArray *fibers = mxGetFieldByNumber(structMx, 0, FIELDNAME_INDEX_fibers);
  if( fibers == 0 )
    {
    mexErrMsgTxt("No fiber data in struct");
    }

  // create and populate the fibers
  vtkSmartPointer<vtkPolyData>  polyData = vtkSmartPointer<vtkPolyData>::New();
  vtkSmartPointer<vtkPoints>    points = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();

  const mwIndex numFibers(mxGetNumberOfElements(fibers) );
  size_t        numPoints = 0;

  // populate points and lines
  vtkIdType counter = 0;
  for( mwIndex i = 0; i < numFibers; ++i )
    {
    mxArray *     curFiber = mxGetCell(fibers, i);
    const double *curFiberPoints = static_cast<const double *>(mxGetData(curFiber) );
    const mwSize *dims = mxGetDimensions(curFiber);
    numPoints += dims[1];

    // handle ids in parallel
    vtkIdType *ids = new vtkIdType[dims[1]];

    double curPt[3];
    for( mwIndex j = 0; j < dims[1]; ++j, ++counter )
      {
      for( mwIndex k = 0; k < 3; ++k )
        {
        curPt[k] = curFiberPoints[(j * dims[0]) + k];
        }
      // add point
      points->InsertNextPoint(curPt[0], curPt[1], curPt[2]);
      // add point id into points
      ids[j] = counter;
      }
    // make a new line, add the ids, insert into lines
    vtkSmartPointer<vtkLine> curLine = vtkSmartPointer<vtkLine>::New();
    curLine->Initialize(dims[1], ids, points);
    lines->InsertNextCell(curLine);
    delete [] ids;
    }
  polyData->SetPoints(points);
  polyData->SetLines(lines);

  vtkPointData *pointData = polyData->GetPointData();
  // handle decorations
  for( mwSize i = 1; i < numberOfFields; ++i )
    {
    mxArray *     currentMxDecoration = mxGetFieldByNumber(structMx, 0, i);
    const double *currentData = static_cast<const double *>(mxGetData(currentMxDecoration) );
    const mwSize *mxDimensions = mxGetDimensions(currentMxDecoration);
    const char *  curName = mxGetFieldNameByNumber(structMx, i);

    if( mxDimensions[0] != numPoints )
      {
      mexErrMsgTxt("mismatch between number of points and number of values in decoratin");
      return;
      }

    vtkSmartPointer<vtkDoubleArray> curDecoration = vtkSmartPointer<vtkDoubleArray>::New();
    curDecoration->SetNumberOfComponents(mxDimensions[1]);
    curDecoration->SetName(curName);
    double *tmp = new double[mxDimensions[1]];
    for( mwSize j = 0; j < mxDimensions[0]; ++j )
      {
      for( mwSize k = 0; k < mxDimensions[1]; ++k )
        {
        tmp[k] = currentData[(mxDimensions[0] * k) + j];
        }
      curDecoration->InsertNextTuple(tmp);
      }
    delete [] tmp;
    int idx = pointData->AddArray(curDecoration);

    switch( mxDimensions[1] )
      {
      case 1:
        pointData->SetActiveAttribute(idx, vtkDataSetAttributes::SCALARS);
        break;
      case 9:
      default:
        // TODO -- try and track what's coming out of matlab --
        // normals? tensors?
        pointData->SetActiveAttribute(idx, vtkDataSetAttributes::TENSORS);
      }
    }
  WritePolyData(polyData, filename, true, true);
}

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
  try
    {
    vtkSavePolyData(nlhs, plhs, nrhs, prhs);
    }
  catch( ... )
    {
    printf("Exception in vtkSavePolyData\n");
    mexErrMsgTxt("Exception in vtkSavePolyData");
    }

}
