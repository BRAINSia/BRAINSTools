#include "mex.h"
#include <sstream>
#include <vector>
#include <sstream>
#include "vtkPointData.h"
#include "vtkDoubleArray.h"
#include "vtkFloatArray.h"
#include "vtkCellArray.h"
#include "vtkSmartPointer.h"
#include "vtkPolyData.h"
#include "vtkPolyDataReader.h"
#include "vtkXMLPolyDataReader.h"
#include "vtksys/SystemTools.hxx"
#include "vtkErrorObserver.h"

#ifdef __DEBUG_THIS
#define DBTRACE(s) mexPrintf("%d: %s\n", __LINE__, s)
#else
#define DBTRACE(s) /* */
#endif
namespace // anon namespace gives all these file scope; only
          // mexFunction is a public symbol...
{
vtkSmartPointer<vtkPolyData>
ReadPolyData(const char *filename)
{
  const std::string ext(vtksys::SystemTools::GetFilenameExtension(filename) );

  vtkSmartPointer<ErrorObserver> errorObserver =
    vtkSmartPointer<ErrorObserver>::New();

  if( ext == ".vtp" )
    {
    vtkSmartPointer<vtkXMLPolyDataReader> reader =
      vtkSmartPointer<vtkXMLPolyDataReader>::New();
    reader->AddObserver(vtkCommand::ErrorEvent, errorObserver);
    reader->AddObserver(vtkCommand::WarningEvent, errorObserver);

    reader->SetFileName(filename);
    reader->Update();

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
      return NULL;
      }
    return vtkSmartPointer<vtkPolyData>(reader->GetOutput() );
    }
  else
    {
    vtkSmartPointer<vtkPolyDataReader> reader =
      vtkSmartPointer<vtkPolyDataReader>::New();

    reader->AddObserver(vtkCommand::ErrorEvent, errorObserver);
    reader->AddObserver(vtkCommand::WarningEvent, errorObserver);

    reader->SetFileName(filename);
    reader->Update();

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
      return NULL;
      }

    return vtkSmartPointer<vtkPolyData>(reader->GetOutput() );
    }
}

inline
char * mxStrdup(const char *s)
{
  char *rval = static_cast<char *>(mxCalloc(strlen(s) + 1, sizeof(char) ) );

  strcpy(rval, s);
  return rval;
}

template <typename TvtkArray, typename TScalar>
mxArray * _makeMxArray(vtkDataArray *currentData)
{
  TvtkArray *curArray = TvtkArray::SafeDownCast(currentData);

  mwSize   numTuples = curArray->GetNumberOfTuples();
  mwSize   numComponents = curArray->GetNumberOfComponents();
  mxArray *curMx = mxCreateDoubleMatrix(numTuples, numComponents, mxREAL);
  double * mxData = static_cast<double *>(mxGetData(curMx) );

  for( vtkIdType j = 0; j < static_cast<vtkIdType>(numTuples); ++j )
    {
    double *curTuple = curArray->GetTuple(j);
    for( mwSize k = 0; k < numComponents; ++k )
      {
      mxData[(numTuples * k) + j] = curTuple[k];
      }
    }
  return curMx;
}

mxArray * makeMxArray(vtkDataArray *currentData)
{
  vtkDoubleArray *curDoubleArray = vtkDoubleArray::SafeDownCast(currentData);

  if( curDoubleArray == 0 )
    {
    vtkFloatArray *curFloatArray = vtkFloatArray::SafeDownCast(currentData);
    if( curFloatArray == 0 )
      {
      return 0;
      }
    return _makeMxArray<vtkFloatArray, float>(currentData);
    }
  return _makeMxArray<vtkDoubleArray, double>(currentData);
}

} // anon namespace
void
vtkLoadPolyData(int /* nhls */, mxArray *plhs[],
                int nrhs, const mxArray *prhs[])
{
  if( nrhs < 1 || !mxIsChar(prhs[0]) )
    {
    mexErrMsgTxt("vtkLoadPolyData requires 1 argument, the name of the VTK or VTP file to load");
    return;
    }
  const size_t filenamelen = mxGetM(prhs[0]) * mxGetN(prhs[0]) + 1;
  char * const filename = static_cast<char *>(mxCalloc(filenamelen, sizeof(mxChar) ) );
  mxGetString(prhs[0], filename, filenamelen);
  // go read the file
  vtkSmartPointer<vtkPolyData> polyData;
  polyData = ReadPolyData(filename);
  DBTRACE("after ReadPolyData");
  if( polyData == NULL )
    {
    return;
    }

  //
  // the matlab Struct will have a named component for the fibers,
  // and a named component for each of the point data decorator arrays.
  // the name of the decorator arrays is the same as what is in the
  // VTK/VTP file
  const char * *fieldnames = static_cast<const char * *>(mxCalloc(1, sizeof(char *) ) );
  mwSize        numElements(0);

  // get the lines
  vtkCellArray *lines = polyData->GetLines();
  const mwSize  numLines(polyData->GetNumberOfLines() );

  if( numLines > 0 )
    {
    ++numElements;
    fieldnames = static_cast<const char * *>(mxRealloc(fieldnames, sizeof(char *) * numElements) );
    fieldnames[0] = mxStrdup("fibers");
    }

  // do a census to get all the decorators
  vtkPointData *pd = polyData->GetPointData();
  mwSize        pointDataStart(numElements);

  if( pd != 0 )
    {
    numElements += pd->GetNumberOfArrays();
    fieldnames = static_cast<const char * *>(mxRealloc(fieldnames, sizeof(char *) * numElements) );
    for( int i = 0; i < pd->GetNumberOfArrays(); ++i )
      {
      const char *arrayName = pd->GetArrayName(i);
      if( arrayName )
        {
        fieldnames[pointDataStart + i] = mxStrdup(arrayName);
        }
      else
        {
        // I guess they can be unnamed, but normally this shouldn't happen.
        std::stringstream ss; ss << "unnamed" << std::setw(4) << i;
        fieldnames[pointDataStart + i] = mxStrdup(ss.str().c_str() );
        }
      }
    }
  DBTRACE("after checking point data");

  // create a structure to hold the lines
  mxArray *structMx = mxCreateStructMatrix(1, 1, numElements, fieldnames);

  vtkPoints *points = polyData->GetPoints();
  mwSize     linesDims[2];
  linesDims[0] = 1;
  linesDims[1] = numLines;
  mxArray *mxLines = mxCreateCellArray(2, linesDims);

  mwSize linedims[2];
  linedims[0] = 3;

  vtkIdType *indices;
  vtkIdType  numberOfPoints;

  lines->InitTraversal();
  for( mwIndex lineIndex = 0; lines->GetNextCell(numberOfPoints, indices) != 0; ++lineIndex )
    {
    linedims[1] = numberOfPoints;
    mxArray *curLine = mxCreateNumericArray(2, linedims, mxDOUBLE_CLASS, mxREAL);
    double * curLinePoints = static_cast<double *>(mxGetData(curLine) );
    for( unsigned i = 0; i < numberOfPoints; ++i )
      {
      double point[3];
      points->GetPoint(indices[i], point);
      for( unsigned j = 0; j < 3; ++j )
        {
        curLinePoints[(i * 3) + j] = point[j];
        }
      }
    mxSetCell(mxLines, lineIndex, curLine);
    }
  mxSetFieldByNumber(structMx, 0, 0, mxLines);
  DBTRACE("after getting lines");
  //
  // point data
  for( mwSize i = pointDataStart; i < numElements; ++i )
    {
    vtkDataArray *currentData = pd->GetArray(fieldnames[i]);
    mxArray *     curMx = makeMxArray(currentData);;
    if( curMx != 0 )
      {
      mxSetFieldByNumber(structMx, 0, i, curMx);
      }
    }
  DBTRACE("after getting decorations");
  plhs[0] = structMx;
}

void mexFunction(int nhls, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
  vtkLoadPolyData(nhls, plhs, nrhs, prhs);
}
