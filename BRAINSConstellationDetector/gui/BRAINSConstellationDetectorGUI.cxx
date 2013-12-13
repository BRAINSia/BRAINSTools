/*
 * Author: Wei Lu
 * at Psychiatry Imaging Lab,
 * University of Iowa Health Care 2010
 */

#include "BRAINSConstellationDetectorGUICLP.h"

#include "QVTKInteractionCallback.h"
#include "QSliceViewer.h"
#include "QLabelList.h"
#include "QDelLabelDialogs.h"
#include "QFileDialogs.h"
#include "QHelpDialog.h"

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkStatisticsImageFilter.h"
#include "itkImageToVTKImageFilter.h"

#include "vtkCommand.h"
#include "vtkImageData.h"
#include "vtkMatrix4x4.h"
#include "vtkImageReslice.h"
#include "vtkLookupTable.h"
#include "vtkImageMapToColors.h"
#include "vtkCamera.h"
#include "vtkImageActor.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkInteractorStyleImage.h"

#include <QApplication>
#include <QDir>
#include <QString>
#include <QWidget>
#include <QSlider>
#include <QLabel>
#include <QVBoxLayout>
#include <QPushButton>
#include <QListWidget>

#include <cstdio>
#include <iostream>
#include <fstream>
#include <cstring>
#include <map>
#include <vector>

#include "BRAINSCommonLib.h"

// Image, filter typedef
namespace
{ // put in anon namespace so they don't collide with other files.
const unsigned int LocalImageDimension = 3;
typedef short LocalPixelType;
}

typedef itk::Image<LocalPixelType, LocalImageDimension> ImageType;
typedef ImageType::Pointer                         ImagePointerType;
typedef ImageType::PointType                       ImagePointType;
typedef ImageType::SpacingType                     ImageSpacingType;
typedef ImageType::SizeType                        ImageSizeType;
typedef ImageType::DirectionType                   ImageDirectionType;
typedef ImageType::IndexType                       ImageIndexType;
typedef ImagePointType::CoordRepType               ImageCoordType;

typedef itk::ImageFileReader<ImageType>       ReaderType;
typedef itk::ImageFileWriter<ImageType>       WriterType;
typedef itk::StatisticsImageFilter<ImageType> StatisticsFilterType;
typedef itk::ImageToVTKImageFilter<ImageType> ConnectorType;

int main(int argc, char *argv[])
{
  PARSE_ARGS;
  BRAINSRegisterAlternateIO();

  if( inputVolume.compare("") == 0 )
    {
    std::cerr << "To run the program please specify the input volume filename." << std::endl;
    std::cerr << "Type " << argv[0] << " -h for more help." << std::endl;
    return EXIT_FAILURE;
    }

  // ------------------------------------
  // Instantiate widgets

  int          WIN_WIDTH = 1200;  // application width
  int          WIN_HEIGHT = 500;  // application height
  double       PRECISION = 1000.; // means .001 precision for slice number
  QApplication app(argc, argv);

  QWidget base;
  base.setWindowTitle("BRAINS Constellation Viewer");
  QGridLayout layout(&base);

  // image info widet
  QLabel infoLabel("- Image Volume Info -");
  layout.addWidget(&infoLabel, 1, 1, 1, 2);
  infoLabel.setAlignment(Qt::AlignHCenter);

  QLabel info;
  layout.addWidget(&info, 2, 1, 8, 2);
  info.setAlignment(Qt::AlignHCenter);
  info.setMaximumWidth(WIN_WIDTH / 5);
  info.setIndent(5);
  info.setWordWrap(1);
  info.setFrameStyle(QFrame::Panel | QFrame::Raised);
  info.setLineWidth(2);

  // slice viewer widget
  QLabel label1("Axial");
  QLabel label2("Sagittal");
  QLabel label3("Coronal");

  label1.setAlignment(Qt::AlignHCenter);
  label2.setAlignment(Qt::AlignHCenter);
  label3.setAlignment(Qt::AlignHCenter);

  layout.addWidget(&label1, 1, 3, 1, 3);
  layout.addWidget(&label2, 1, 6, 1, 3);
  layout.addWidget(&label3, 1, 9, 1, 3);

  QSliceViewer sliceViewer1(1);
  QSliceViewer sliceViewer2(2);
  QSliceViewer sliceViewer3(3);

  layout.addWidget(&sliceViewer1, 2, 3, 9, 3);
  layout.addWidget(&sliceViewer2, 2, 6, 9, 3);
  layout.addWidget(&sliceViewer3, 2, 9, 9, 3);

  // slider bar widget
  QSlider slider1(Qt::Horizontal);
  QSlider slider2(Qt::Horizontal);
  QSlider slider3(Qt::Horizontal);

  layout.addWidget(&slider1, 11, 3, 1, 3);
  layout.addWidget(&slider2, 11, 6, 1, 3);
  layout.addWidget(&slider3, 11, 9, 1, 3);

  // status label widget
  QLabel status1("");
  QLabel status2("");
  QLabel status3("");

  layout.addWidget(&status1, 12, 3, 1, 3);
  layout.addWidget(&status2, 12, 6, 1, 3);
  layout.addWidget(&status3, 12, 9, 1, 3);

  status1.setAlignment(Qt::AlignHCenter);
  status2.setAlignment(Qt::AlignHCenter);
  status3.setAlignment(Qt::AlignHCenter);

  // label list widget
  QLabelList labelList;
  layout.addWidget(&labelList, 10, 1, 4, 2);
  labelList.SetInputVolume(inputVolume);
  labelList.SetInputLandmarks(inputLandmarks);
  labelList.SetOutputLandmarks(outputLandmarks);

  // add label dialog widget
  QPushButton addLabelButton("Add");
  layout.addWidget(&addLabelButton, 13, 3);

  // remove label dialog widget
  QPushButton removeLabelButton("Remove");
  layout.addWidget(&removeLabelButton, 13, 4);

  // remove all dialog widget
  QPushButton removeAllButton("Remove All");
  layout.addWidget(&removeAllButton, 13, 5);

  // save landmarks widget
  QPushButton saveButton("Save");
  layout.addWidget(&saveButton, 13, 6);

  // save as landmarks widget
  QPushButton saveAsButton("Save As...");
  layout.addWidget(&saveAsButton, 13, 7);

  // help widget
  QPushButton helpButton("Help");
  layout.addWidget(&helpButton, 13, 10);
  QHelpDialog helpDialog;

  // quit widget
  QPushButton quitButton("Quit");
  layout.addWidget(&quitButton, 13, 11);
  QObject::connect( &quitButton, SIGNAL( clicked() ),
                    &app, SLOT( quit() ) );

  // label deletion dialog
  QDelLabelDialogs delAllButtonDialog("Are you sure to remove all the labels?");
  QDelLabelDialogs delCurrButtonDialog("Are you sure to remove the current label?");
  QDelLabelDialogs delCurrKeyDialog("Are you sure to remove the current label?");
  QDelLabelDialogs doubleClickLabelDialog("Are you sure to remove the current label?");

  // ------------------------------------
  // Read input image from ITK
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName(inputVolume);
  try
    {
    reader->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    std::cerr << " Error while reading image file(s) with ITK:\n "
              << err << std::endl;
    }

  // ------------------------------------
  // Prepare for rendering image with VTK

  // get some basic information of the image using ITK
  ImagePointerType   image      = reader->GetOutput();
  ImagePointType     origin     = image->GetOrigin();
  ImageSpacingType   spacing    = image->GetSpacing();
  ImageSizeType      size       = image->GetLargestPossibleRegion().GetSize();
  ImageDirectionType direction  = image->GetDirection();

  // get the image extent in index space ( IJK coordinate system )
  int            indexExtent[6];
  ImageIndexType startIndex = image->GetLargestPossibleRegion().GetIndex();
  ImageIndexType stopIndex;
  stopIndex[0]    = startIndex[0] + size[0] - 1;
  stopIndex[1]    = startIndex[1] + size[1] - 1;
  stopIndex[2]    = startIndex[2] + size[2] - 1;
  indexExtent[0]  = startIndex[0];
  indexExtent[2]  = startIndex[1];
  indexExtent[4]  = startIndex[2];
  indexExtent[1]  = stopIndex[0];
  indexExtent[3]  = stopIndex[1];
  indexExtent[5]  = stopIndex[2];

  // calculate the physical image extent of input image with identity
  // direction
  ImageDirectionType direction2;
  direction2.Fill(0);
  direction2(0, 0) = 1.0;
  direction2(1, 1) = 1.0;
  direction2(2, 2) = 1.0;
  image->SetDirection(direction2);

  double         physicalExtentIdentity[6];
  ImagePointType physicalStartLocation;
  ImagePointType physicalStopLocation;
  image->TransformIndexToPhysicalPoint(startIndex, physicalStartLocation);
  image->TransformIndexToPhysicalPoint(stopIndex, physicalStopLocation);
  physicalExtentIdentity[0] = physicalStartLocation[0];
  physicalExtentIdentity[2] = physicalStartLocation[1];
  physicalExtentIdentity[4] = physicalStartLocation[2];
  physicalExtentIdentity[1] = physicalStopLocation[0];
  physicalExtentIdentity[3] = physicalStopLocation[1];
  physicalExtentIdentity[5] = physicalStopLocation[2];

  // calculate the physical center of the image with identity direction
  double center[4];
  center[0] = ( physicalExtentIdentity[1] + physicalExtentIdentity[0] ) / 2.0;
  center[1] = ( physicalExtentIdentity[3] + physicalExtentIdentity[2] ) / 2.0;
  center[2] = ( physicalExtentIdentity[5] + physicalExtentIdentity[4] ) / 2.0;
  center[3] = 1.0;

  // recalculate the image extent of a directional image in physical space
  double physicalExtent[6];
  image->SetDirection(direction);
  image->TransformIndexToPhysicalPoint(startIndex, physicalStartLocation);
  image->TransformIndexToPhysicalPoint(stopIndex, physicalStopLocation);
  physicalExtent[0] = physicalStartLocation[0];
  physicalExtent[2] = physicalStartLocation[1];
  physicalExtent[4] = physicalStartLocation[2];
  physicalExtent[1] = physicalStopLocation[0];
  physicalExtent[3] = physicalStopLocation[1];
  physicalExtent[5] = physicalStopLocation[2];

  // calculate the physical center of the image with original direction
  double center2[3];
  center2[0] = ( physicalExtent[1] + physicalExtent[0] ) / 2.0;
  center2[1] = ( physicalExtent[3] + physicalExtent[2] ) / 2.0;
  center2[2] = ( physicalExtent[5] + physicalExtent[4] ) / 2.0;

  // compute the bound of image
  double bound[6];
  bound[0] =
    physicalExtent[0] < physicalExtent[1] ?
    physicalExtent[0] : physicalExtent[1];
  bound[1] =
    physicalExtent[0] >= physicalExtent[1] ?
    physicalExtent[0] : physicalExtent[1];
  bound[2] =
    physicalExtent[2] < physicalExtent[3] ?
    physicalExtent[2] : physicalExtent[3];
  bound[3] =
    physicalExtent[2] >= physicalExtent[3] ?
    physicalExtent[2] : physicalExtent[3];
  bound[4] =
    physicalExtent[4] < physicalExtent[5] ?
    physicalExtent[4] : physicalExtent[5];
  bound[5] =
    physicalExtent[4] >= physicalExtent[5] ?
    physicalExtent[4] : physicalExtent[5];

  // set the range of sliders
  slider1.setRange(bound[4] * PRECISION, bound[5] * PRECISION);
  slider2.setRange(bound[0] * PRECISION, bound[1] * PRECISION);
  slider3.setRange(bound[2] * PRECISION, bound[3] * PRECISION);

  // set the initial positon of the sliders
  slider1.setSliderPosition( ( bound[4] + bound[5] ) / 2 * PRECISION );
  slider2.setSliderPosition( ( bound[0] + bound[1] ) / 2 * PRECISION );
  slider3.setSliderPosition( ( bound[2] + bound[3] ) / 2 * PRECISION );

  // pass the bound info to slice viewer
  sliceViewer1.SetBound(bound);
  sliceViewer2.SetBound(bound);
  sliceViewer3.SetBound(bound);

  // show welcome info and some image information
  QString infoText;
  infoText =
    "<p><b>Physical Spacing:</b><br /><font size=2>("
    + infoText.number(spacing[0]) + " "
    + infoText.number(spacing[1]) + " "
    + infoText.number(spacing[2]) + ")</font></p>"
    + "<p><b>Physical Origin:</b><br /><font size=2>("
    + infoText.number(origin[0]) + " "
    + infoText.number(origin[1]) + " "
    + infoText.number(origin[2]) + ")</font></p>"
    + "<p><b>Physical Center:</b><br /><font size=2>("
    + infoText.number(center2[0]) + " "
    + infoText.number(center2[1]) + " "
    + infoText.number(center2[2]) + ")</font></p>"
    + "<p><b>Physical Direction:</b><br /><font size=2>"
    + infoText.number( direction(0, 0) ) + " "
    + infoText.number( direction(0, 1) )
    + " " + infoText.number( direction(0, 2) )
    + "<br />"
    + infoText.number( direction(1, 0) ) + " "
    + infoText.number( direction(1, 1) )
    + " " + infoText.number( direction(1, 2) )
    + "<br />"
    + infoText.number( direction(2, 0) ) + " "
    + infoText.number( direction(2, 1) )
    + " " + infoText.number( direction(2, 2) )
    + "</font></p>"
    + "<p><b>Index Extent:</b><br /><font size=2>["
    + infoText.number(indexExtent[0]) + " "
    + infoText.number(indexExtent[1]) + "; "
    + infoText.number(indexExtent[2]) + " "
    + infoText.number(indexExtent[3]) + "; "
    + infoText.number(indexExtent[4]) + " "
    + infoText.number(indexExtent[5]) + "]</font></p>"
    + "<p><b>Physical Extent:</b><br /><font size=2>("
    + infoText.number(bound[0]) + " "
    + infoText.number(bound[1]) + "; "
    + infoText.number(bound[2]) + " "
    + infoText.number(bound[3]) + "; "
    + infoText.number(bound[4]) + " "
    + infoText.number(bound[5]) + ")</font></p>";

  info.setText(infoText);

  // ------------------------------------
  // ITK to VTK

  ConnectorType::Pointer connector = ConnectorType::New();
  connector->SetInput(image);

  // ------------------------------------
  // get image slices with VTK Image Reslice

  // Matrices for axial, coronal, and sagittal view orientations
  // in LAS coordinate
  static double axialElements[16] =
    {
    1, 0, 0, 0,
    0, -1, 0, 0,
    0, 0, 1, 0,
    0, 0, 0, 1
    };

  // in ASL coordinate
  static double sagittalElements[16] =
    {
    0, 0, 1, 0,
    1, 0, 0, 0, // later we need to change it from -1 to 1
    0, 1, 0, 0,
    0, 0, 0, 1
    };

  // in LSA coordinate
  static double coronalElements[16] =
    {
    1, 0, 0, 0,
    0, 0, -1, 0,
    0, 1, 0, 0,
    0, 0, 0, 1
    };

  // ITK -> VTK pipeline does not transfer direction infomation!
  double rotation[16] =
    {
    0, 0, 0, 0,
    0, 0, 0, 0,
    0, 0, 0, 0,
    0, 0, 0, 1
    };

  unsigned int i, j;
  for( i = 0; i < LocalImageDimension; ++i )
    {
    for( j = 0; j < LocalImageDimension; ++j )
      {
      rotation[4 * i + j] = direction[i][j];
      }
    }

  // compute the reslice axes location
  double        matrix[16];
  double        invRotation[16];
  vtkMatrix4x4 *resliceAxes1 = vtkMatrix4x4::New();
  vtkMatrix4x4 *resliceAxes2 = vtkMatrix4x4::New();
  vtkMatrix4x4 *resliceAxes3 = vtkMatrix4x4::New();
  resliceAxes1->Invert(rotation, invRotation);

  resliceAxes1->Multiply4x4(invRotation, axialElements, matrix);
  resliceAxes1->DeepCopy(matrix);
  resliceAxes1->SetElement(0, 3, center[0]);
  resliceAxes1->SetElement(1, 3, center[1]);
  resliceAxes1->SetElement(2, 3, center[2]);

  resliceAxes2->Multiply4x4(invRotation, sagittalElements, matrix);
  resliceAxes2->DeepCopy(matrix);
  resliceAxes2->SetElement(0, 3, center[0]);
  resliceAxes2->SetElement(1, 3, center[1]);
  resliceAxes2->SetElement(2, 3, center[2]);

  resliceAxes3->Multiply4x4(invRotation, coronalElements, matrix);
  resliceAxes3->DeepCopy(matrix);
  resliceAxes3->SetElement(0, 3, center[0]);
  resliceAxes3->SetElement(1, 3, center[1]);
  resliceAxes3->SetElement(2, 3, center[2]);

  // extract a slice in the desired orientation
  vtkImageReslice *reslice1 = vtkImageReslice::New();
  reslice1->SetInputConnection( connector->GetImporter()->GetOutputPort() );
  reslice1->SetOutputDimensionality(2);
  reslice1->SetResliceAxes(resliceAxes1);
  reslice1->SetInterpolationModeToLinear();

  vtkImageReslice *reslice2 = vtkImageReslice::New();
  reslice2->SetInputConnection( connector->GetImporter()->GetOutputPort() );
  reslice2->SetOutputDimensionality(2);
  reslice2->SetResliceAxes(resliceAxes2);
  reslice2->SetInterpolationModeToLinear();

  vtkImageReslice *reslice3 = vtkImageReslice::New();
  reslice3->SetInputConnection( connector->GetImporter()->GetOutputPort() );
  reslice3->SetOutputDimensionality(2);
  reslice3->SetResliceAxes(resliceAxes3);
  reslice3->SetInterpolationModeToLinear();

  // find the slice size
  reslice1->Update();
  reslice2->Update();
  reslice3->Update();

  // ------------------------------------
  // Render the image with VTK

  // create a colormap lookup table
  StatisticsFilterType::Pointer statisticsFilter = StatisticsFilterType::New();
  statisticsFilter->SetInput( reader->GetOutput() );
  statisticsFilter->Update();
  LocalPixelType minPixelValue = statisticsFilter->GetMinimum();
  LocalPixelType maxPixelValue = statisticsFilter->GetMaximum();

  vtkLookupTable *table = vtkLookupTable::New();
  table->SetRange(minPixelValue, maxPixelValue);   // image intensity range
  table->SetValueRange(0.0, 1.0);
  table->SetSaturationRange(0.0, 0.0);   // monochromatic
  table->SetRampToLinear();
  table->Build();

  // map the image through the lookup table
  vtkImageMapToColors *color1 = vtkImageMapToColors::New();
  color1->SetLookupTable(table);
  color1->SetInputConnection( reslice1->GetOutputPort() );

  vtkImageMapToColors *color2 = vtkImageMapToColors::New();
  color2->SetLookupTable(table);
  color2->SetInputConnection( reslice2->GetOutputPort() );

  vtkImageMapToColors *color3 = vtkImageMapToColors::New();
  color3->SetLookupTable(table);
  color3->SetInputConnection( reslice3->GetOutputPort() );

  // set up the actor
  vtkImageActor *actor1 = vtkImageActor::New();
  actor1->SetInput( color1->GetOutput() );

  vtkImageActor *actor2 = vtkImageActor::New();
  actor2->SetInput( color2->GetOutput() );

  vtkImageActor *actor3 = vtkImageActor::New();
  actor3->SetInput( color3->GetOutput() );

  // set up the renderer
  vtkRenderer *renderer1 = vtkRenderer::New();
  renderer1->AddActor(actor1);
  renderer1->SetBackground(0.45, 0.65, 0.45);

  vtkRenderer *renderer2 = vtkRenderer::New();
  renderer2->AddActor(actor2);
  renderer2->SetBackground(0.45, 0.65, 0.45);

  vtkRenderer *renderer3 = vtkRenderer::New();
  renderer3->AddActor(actor3);
  renderer3->SetBackground(0.45, 0.65, 0.45);

  // set up the render window
  vtkRenderWindow *window1 = vtkRenderWindow::New();
  window1->AddRenderer(renderer1);
  sliceViewer1.SetRenderWindow(window1);

  vtkRenderWindow *window2 = vtkRenderWindow::New();
  window2->AddRenderer(renderer2);
  sliceViewer2.SetRenderWindow(window2);

  vtkRenderWindow *window3 = vtkRenderWindow::New();
  window3->AddRenderer(renderer3);
  sliceViewer3.SetRenderWindow(window3);

  // set up the interaction
  vtkInteractorStyleImage *imageStyle1 = vtkInteractorStyleImage::New();
  QVTKInteractor *         interactor1 = QVTKInteractor::New();
  interactor1->SetInteractorStyle(imageStyle1);
  window1->SetInteractor(interactor1);

  window1->Render();

  QVTKInteractionCallback callback1("L:\nP:\nS:", 1, bound);
  callback1.SetImageReslice(reslice1);
  callback1.SetInteractor(interactor1);
  callback1.SetIndexExtent(indexExtent);
  callback1.SetPhysicalExtentIdentity(physicalExtentIdentity);
  callback1.SetPhysicalExtent(physicalExtent);
  callback1.SetPrecision(PRECISION);

  vtkInteractorStyleImage *imageStyle2 = vtkInteractorStyleImage::New();
  QVTKInteractor *         interactor2 = QVTKInteractor::New();
  interactor2->SetInteractorStyle(imageStyle2);
  window2->SetInteractor(interactor2);

  window2->Render();

  QVTKInteractionCallback callback2("L:\nP:\nS:", 2, bound);
  callback2.SetImageReslice(reslice2);
  callback2.SetInteractor(interactor2);
  callback2.SetIndexExtent(indexExtent);
  callback2.SetPhysicalExtentIdentity(physicalExtentIdentity);
  callback2.SetPhysicalExtent(physicalExtent);
  callback2.SetType(2);
  callback2.SetPrecision(PRECISION);

  vtkInteractorStyleImage *imageStyle3 = vtkInteractorStyleImage::New();
  QVTKInteractor *         interactor3 = QVTKInteractor::New();
  interactor3->SetInteractorStyle(imageStyle3);
  window3->SetInteractor(interactor3);

  window3->Render();

  QVTKInteractionCallback callback3("L:\nP:\nS:", 3, bound);
  callback3.SetImageReslice(reslice3);
  callback3.SetInteractor(interactor3);
  callback3.SetIndexExtent(indexExtent);
  callback3.SetPhysicalExtentIdentity(physicalExtentIdentity);
  callback3.SetPhysicalExtent(physicalExtent);
  callback3.SetType(3);
  callback3.SetPrecision(PRECISION);

  callback1.SetDirection(rotation);
  callback2.SetDirection(rotation);
  callback3.SetDirection(rotation);

  double origin16[16] =
    {
    1, 0, 0, 0,
    0, 1, 0, 0,
    0, 0, 1, 0,
    0, 0, 0, 1
    };
  origin16[3] = origin[0];
  origin16[7] = origin[1];
  origin16[11] = origin[2];
  callback1.SetOrigin(origin16);
  callback2.SetOrigin(origin16);
  callback3.SetOrigin(origin16);

  double spacing16[16] =
    {
    1, 0, 0, 0,
    0, 1, 0, 0,
    0, 0, 1, 0,
    0, 0, 0, 1
    };
  spacing16[0] = spacing[0];
  spacing16[5] = spacing[1];
  spacing16[10] = spacing[2];
  callback1.SetSpacing(spacing16);
  callback2.SetSpacing(spacing16);
  callback3.SetSpacing(spacing16);

  imageStyle1->AddObserver(vtkCommand::MouseMoveEvent, &callback1);
  imageStyle1->AddObserver(vtkCommand::LeftButtonPressEvent, &callback1);
  imageStyle1->AddObserver(vtkCommand::LeftButtonReleaseEvent, &callback1);
  imageStyle1->AddObserver(vtkCommand::KeyPressEvent, &callback1);

  imageStyle2->AddObserver(vtkCommand::MouseMoveEvent, &callback2);
  imageStyle2->AddObserver(vtkCommand::LeftButtonPressEvent, &callback2);
  imageStyle2->AddObserver(vtkCommand::LeftButtonReleaseEvent, &callback2);
  imageStyle2->AddObserver(vtkCommand::KeyPressEvent, &callback2);

  imageStyle3->AddObserver(vtkCommand::MouseMoveEvent, &callback3);
  imageStyle3->AddObserver(vtkCommand::LeftButtonPressEvent, &callback3);
  imageStyle3->AddObserver(vtkCommand::LeftButtonReleaseEvent, &callback3);
  imageStyle3->AddObserver(vtkCommand::KeyPressEvent, &callback3);

  // ------------------------------------
  // Load & display landmarks from file
  labelList.loadLandmarks();
  std::map<QString, std::vector<double> >           landmarks = labelList.GetLandmarks();
  std::map<QString, std::vector<double> >::iterator it;
  for( it = landmarks.begin(); it != landmarks.end(); ++it )
    {
    if( it->first.compare("") != 0 )
      {
      // add points to slice viewer
      double *point = new double[3];
      for( i = 0; i < LocalImageDimension; ++i )
        {
        point[i] = ( it->second )[i];
        }
      sliceViewer1.createLabel(point);
      sliceViewer2.createLabel(point);
      sliceViewer3.createLabel(point);

      // add points to label list
      QString pointString =
        pointString.number(point[0]) + " "
        + pointString.number(point[1]) + " "
        + pointString.number(point[2]);
      labelList.createListItem(pointString, it->first);
      }
    }

  // ------------------------------------
  // Setup core communication

  // connect slice viewer, slider bar, and status label
  QObject::connect( &callback1, SIGNAL( textChanged(const QString &) ),
                    &status1, SLOT( setText(const QString &) ) );
  QObject::connect( &callback1, SIGNAL( valueChangedZ(int) ),
                    &slider1, SLOT( setValue(int) ) );
  QObject::connect( &slider1, SIGNAL( valueChanged(int) ),
                    &callback1, SLOT( setValueZLocation(int) ) );

  QObject::connect( &callback2, SIGNAL( textChanged(const QString &) ),
                    &status2, SLOT( setText(const QString &) ) );
  QObject::connect( &callback2, SIGNAL( valueChangedZ(int) ),
                    &slider2, SLOT( setValue(int) ) );
  QObject::connect( &slider2, SIGNAL( valueChanged(int) ),
                    &callback2, SLOT( setValueZLocation(int) ) );

  QObject::connect( &callback3, SIGNAL( textChanged(const QString &) ),
                    &status3, SLOT( setText(const QString &) ) );
  QObject::connect( &callback3, SIGNAL( valueChangedZ(int) ),
                    &slider3, SLOT( setValue(int) ) );
  QObject::connect( &slider3, SIGNAL( valueChanged(int) ),
                    &callback3, SLOT( setValueZLocation(int) ) );

  // connect slice viewers
  QObject::connect( &callback2, SIGNAL( valueChangedY(double) ),
                    &callback1, SLOT( setValueYLocation(double) ) );
  QObject::connect( &callback3, SIGNAL( valueChangedY(double) ),
                    &callback1, SLOT( setValueYLocation(double) ) );

  QObject::connect( &callback1, SIGNAL( valueChangedX(double) ),
                    &callback2, SLOT( setValueXLocation(double) ) );
  QObject::connect( &callback3, SIGNAL( valueChangedX(double) ),
                    &callback2, SLOT( setValueXLocation(double) ) );

  QObject::connect( &callback1, SIGNAL( valueChangedY(double) ),
                    &callback3, SLOT( setValueYLocation(double) ) );
  QObject::connect( &callback2, SIGNAL( valueChangedX(double) ),
                    &callback3, SLOT( setValueXLocation(double) ) );

  // connect mouse/keyboard callback with fiducial labels
  QObject::connect( &callback1, SIGNAL( moveActor(double *) ),
                    &sliceViewer1, SLOT( moveLabelSlot(double *) ) );
  QObject::connect( &callback1, SIGNAL( moveActor(double *) ),
                    &sliceViewer2, SLOT( moveLabelSlot(double *) ) );
  QObject::connect( &callback1, SIGNAL( moveActor(double *) ),
                    &sliceViewer3, SLOT( moveLabelSlot(double *) ) );

  QObject::connect( &callback1, SIGNAL( createActor() ),
                    &sliceViewer1, SLOT( createLabelSlot() ) );
  QObject::connect( &callback1, SIGNAL( createActor() ),
                    &sliceViewer2, SLOT( createLabelSlot() ) );
  QObject::connect( &callback1, SIGNAL( createActor() ),
                    &sliceViewer3, SLOT( createLabelSlot() ) );

  QObject::connect( &callback1, SIGNAL( switchActor() ),
                    &sliceViewer1, SLOT( switchLabelSlot() ) );
  QObject::connect( &callback1, SIGNAL( switchActor() ),
                    &sliceViewer2, SLOT( switchLabelSlot() ) );
  QObject::connect( &callback1, SIGNAL( switchActor() ),
                    &sliceViewer3, SLOT( switchLabelSlot() ) );

  QObject::connect( &callback2, SIGNAL( moveActor(double *) ),
                    &sliceViewer2, SLOT( moveLabelSlot(double *) ) );
  QObject::connect( &callback2, SIGNAL( moveActor(double *) ),
                    &sliceViewer3, SLOT( moveLabelSlot(double *) ) );
  QObject::connect( &callback2, SIGNAL( moveActor(double *) ),
                    &sliceViewer1, SLOT( moveLabelSlot(double *) ) );

  QObject::connect( &callback2, SIGNAL( createActor() ),
                    &sliceViewer2, SLOT( createLabelSlot() ) );
  QObject::connect( &callback2, SIGNAL( createActor() ),
                    &sliceViewer3, SLOT( createLabelSlot() ) );
  QObject::connect( &callback2, SIGNAL( createActor() ),
                    &sliceViewer1, SLOT( createLabelSlot() ) );

  QObject::connect( &callback2, SIGNAL( switchActor() ),
                    &sliceViewer2, SLOT( switchLabelSlot() ) );
  QObject::connect( &callback2, SIGNAL( switchActor() ),
                    &sliceViewer3, SLOT( switchLabelSlot() ) );
  QObject::connect( &callback2, SIGNAL( switchActor() ),
                    &sliceViewer1, SLOT( switchLabelSlot() ) );

  QObject::connect( &callback3, SIGNAL( moveActor(double *) ),
                    &sliceViewer3, SLOT( moveLabelSlot(double *) ) );
  QObject::connect( &callback3, SIGNAL( moveActor(double *) ),
                    &sliceViewer1, SLOT( moveLabelSlot(double *) ) );
  QObject::connect( &callback3, SIGNAL( moveActor(double *) ),
                    &sliceViewer2, SLOT( moveLabelSlot(double *) ) );

  QObject::connect( &callback3, SIGNAL( createActor() ),
                    &sliceViewer3, SLOT( createLabelSlot() ) );
  QObject::connect( &callback3, SIGNAL( createActor() ),
                    &sliceViewer1, SLOT( createLabelSlot() ) );
  QObject::connect( &callback3, SIGNAL( createActor() ),
                    &sliceViewer2, SLOT( createLabelSlot() ) );

  QObject::connect( &callback3, SIGNAL( switchActor() ),
                    &sliceViewer3, SLOT( switchLabelSlot() ) );
  QObject::connect( &callback3, SIGNAL( switchActor() ),
                    &sliceViewer1, SLOT( switchLabelSlot() ) );
  QObject::connect( &callback3, SIGNAL( switchActor() ),
                    &sliceViewer2, SLOT( switchLabelSlot() ) );

  QObject::connect( &callback1, SIGNAL( deleteActor() ),
                    &delCurrKeyDialog, SLOT( exec() ) );
  QObject::connect( &callback2, SIGNAL( deleteActor() ),
                    &delCurrKeyDialog, SLOT( exec() ) );
  QObject::connect( &callback3, SIGNAL( deleteActor() ),
                    &delCurrKeyDialog, SLOT( exec() ) );
  QObject::connect( &delCurrKeyDialog, SIGNAL( accepted() ),
                    &sliceViewer1, SLOT( deleteLabelSlot() ) );
  QObject::connect( &delCurrKeyDialog, SIGNAL( accepted() ),
                    &sliceViewer2, SLOT( deleteLabelSlot() ) );
  QObject::connect( &delCurrKeyDialog, SIGNAL( accepted() ),
                    &sliceViewer3, SLOT( deleteLabelSlot() ) );
  QObject::connect( &delCurrKeyDialog, SIGNAL( accepted() ),
                    &labelList, SLOT( deleteListItemSlot() ) );

  QObject::connect( &labelList, SIGNAL( itemClicked(QListWidgetItem *) ),
                    &sliceViewer1, SLOT( pickLabelSlot(QListWidgetItem *) ) );
  QObject::connect( &labelList, SIGNAL( itemClicked(QListWidgetItem *) ),
                    &sliceViewer2, SLOT( pickLabelSlot(QListWidgetItem *) ) );
  QObject::connect( &labelList, SIGNAL( itemClicked(QListWidgetItem *) ),
                    &sliceViewer3, SLOT( pickLabelSlot(QListWidgetItem *) ) );

  QObject::connect( &labelList, SIGNAL( itemClicked(QListWidgetItem *) ),
                    &labelList, SLOT( cancelHighlight(QListWidgetItem *) ) );
  QObject::connect( &labelList, SIGNAL( itemClicked(QListWidgetItem *) ),
                    &labelList, SLOT( cancelHighlight(QListWidgetItem *) ) );
  QObject::connect( &labelList, SIGNAL( itemClicked(QListWidgetItem *) ),
                    &labelList, SLOT( cancelHighlight(QListWidgetItem *) ) );

  // connect mouse/keyboard callback with fiducial label list
  QObject::connect( &callback1, SIGNAL( createListItem(const QString &) ),
                    &labelList, SLOT( createListItemSlot(const QString &) ) );
  QObject::connect( &callback1, SIGNAL( editListItem(const QString &) ),
                    &labelList, SLOT( editListItemSlot(const QString &) ) );
  QObject::connect( &callback1, SIGNAL( switchListItem() ),
                    &labelList, SLOT( switchListItemSlot() ) );

  QObject::connect( &callback2, SIGNAL( createListItem(const QString &) ),
                    &labelList, SLOT( createListItemSlot(const QString &) ) );
  QObject::connect( &callback2, SIGNAL( editListItem(const QString &) ),
                    &labelList, SLOT( editListItemSlot(const QString &) ) );
  QObject::connect( &callback2, SIGNAL( switchListItem() ),
                    &labelList, SLOT( switchListItemSlot() ) );

  QObject::connect( &callback3, SIGNAL( createListItem(const QString &) ),
                    &labelList, SLOT( createListItemSlot(const QString &) ) );
  QObject::connect( &callback3, SIGNAL( editListItem(const QString &) ),
                    &labelList, SLOT( editListItemSlot(const QString &) ) );
  QObject::connect( &callback3, SIGNAL( switchListItem() ),
                    &labelList, SLOT( switchListItemSlot() ) );

  QObject::connect( &labelList, SIGNAL( itemDoubleClicked(QListWidgetItem *) ),
                    &doubleClickLabelDialog, SLOT( exec(QListWidgetItem *) ) );
  QObject::connect( &doubleClickLabelDialog, SIGNAL( accepted(QListWidgetItem *) ),
                    &labelList, SLOT( deleteListItemMouseSlot(QListWidgetItem *) ) );
  QObject::connect( &doubleClickLabelDialog, SIGNAL( accepted(QListWidgetItem *) ),
                    &sliceViewer1, SLOT( deleteLabelMouseSlot(QListWidgetItem *) ) );
  QObject::connect( &doubleClickLabelDialog, SIGNAL( accepted(QListWidgetItem *) ),
                    &sliceViewer2, SLOT( deleteLabelMouseSlot(QListWidgetItem *) ) );
  QObject::connect( &doubleClickLabelDialog, SIGNAL( accepted(QListWidgetItem *) ),
                    &sliceViewer3, SLOT( deleteLabelMouseSlot(QListWidgetItem *) ) );

  // callback or list (send slice update request) >>
  // list (send label position) >> callback (update label position)
  QObject::connect( &labelList, SIGNAL( sliceChangedList() ),
                    &labelList, SLOT( sliceChangedSlot() ) );

  QObject::connect( &labelList, SIGNAL( sendLabelPosition(double *) ),
                    &callback1, SLOT( receiveLabelPos(double *) ) );
  QObject::connect( &labelList, SIGNAL( sendLabelPosition(double *) ),
                    &callback2, SLOT( receiveLabelPos(double *) ) );
  QObject::connect( &labelList, SIGNAL( sendLabelPosition(double *) ),
                    &callback3, SLOT( receiveLabelPos(double *) ) );

  QObject::connect( &labelList, SIGNAL( sendLabelPosition(double *) ),
                    &callback1, SLOT( receiveLabelPos(double *) ) );
  QObject::connect( &labelList, SIGNAL( sendLabelPosition(double *) ),
                    &callback2, SLOT( receiveLabelPos(double *) ) );
  QObject::connect( &labelList, SIGNAL( sendLabelPosition(double *) ),
                    &callback3, SLOT( receiveLabelPos(double *) ) );

  QObject::connect( &labelList, SIGNAL( sendLabelPosition(double *) ),
                    &callback1, SLOT( receiveLabelPos(double *) ) );
  QObject::connect( &labelList, SIGNAL( sendLabelPosition(double *) ),
                    &callback2, SLOT( receiveLabelPos(double *) ) );
  QObject::connect( &labelList, SIGNAL( sendLabelPosition(double *) ),
                    &callback3, SLOT( receiveLabelPos(double *) ) );

  // add label widget to slice viewer, status bar and label list
  QObject::connect( &addLabelButton, SIGNAL( clicked() ),
                    &sliceViewer1, SLOT( createLabelSlot() ) );
  QObject::connect( &addLabelButton, SIGNAL( clicked() ),
                    &sliceViewer2, SLOT( createLabelSlot() ) );
  QObject::connect( &addLabelButton, SIGNAL( clicked() ),
                    &sliceViewer3, SLOT( createLabelSlot() ) );

  QObject::connect( &addLabelButton, SIGNAL( clicked() ),
                    &callback1, SLOT( createListAddButtonSlot() ) );

  // remove label widget to slice viewer, status bar and label list
  QObject::connect( &removeLabelButton, SIGNAL( clicked() ),
                    &delCurrButtonDialog, SLOT( exec() ) );

  QObject::connect( &delCurrButtonDialog, SIGNAL( accepted() ),
                    &sliceViewer1, SLOT( deleteLabelSlot() ) );
  QObject::connect( &delCurrButtonDialog, SIGNAL( accepted() ),
                    &sliceViewer2, SLOT( deleteLabelSlot() ) );
  QObject::connect( &delCurrButtonDialog, SIGNAL( accepted() ),
                    &sliceViewer3, SLOT( deleteLabelSlot() ) );
  QObject::connect( &delCurrButtonDialog, SIGNAL( accepted() ),
                    &labelList, SLOT( deleteListItemSlot() ) );

  // remove all widget
  QObject::connect( &removeAllButton, SIGNAL( clicked() ),
                    &delAllButtonDialog, SLOT( exec() ) );

  QObject::connect( &delAllButtonDialog, SIGNAL( accepted() ),
                    &sliceViewer1, SLOT( deleteAllLabelSlot() ) );
  QObject::connect( &delAllButtonDialog, SIGNAL( accepted() ),
                    &sliceViewer2, SLOT( deleteAllLabelSlot() ) );
  QObject::connect( &delAllButtonDialog, SIGNAL( accepted() ),
                    &sliceViewer3, SLOT( deleteAllLabelSlot() ) );
  QObject::connect( &delAllButtonDialog, SIGNAL( accepted() ),
                    &labelList, SLOT( deleteListSlot() ) );

  // check visibility
  QObject::connect( &callback1, SIGNAL( checkVisibility() ),
                    &labelList, SLOT( checkVisibilitySlot() ) );
  QObject::connect( &callback2, SIGNAL( checkVisibility() ),
                    &labelList, SLOT( checkVisibilitySlot() ) );
  QObject::connect( &callback3, SIGNAL( checkVisibility() ),
                    &labelList, SLOT( checkVisibilitySlot() ) );

  QObject::connect( &callback1, SIGNAL( visibilityCallback(double *) ),
                    &labelList, SLOT( checkVisibilitySlot(double *) ) );
  QObject::connect( &callback2, SIGNAL( visibilityCallback(double *) ),
                    &labelList, SLOT( checkVisibilitySlot(double *) ) );
  QObject::connect( &callback3, SIGNAL( visibilityCallback(double *) ),
                    &labelList, SLOT( checkVisibilitySlot(double *) ) );

  QObject::connect( &labelList, SIGNAL( itemDoubleClicked(QListWidgetItem *) ),
                    &labelList, SLOT( checkVisibilitySlot(QListWidgetItem *) ) );

  QObject::connect( &labelList, SIGNAL( visibilityTable(int *) ),
                    &sliceViewer1, SLOT( visibilityUpdate(int *) ) );
  QObject::connect( &labelList, SIGNAL( visibilityTable(int *) ),
                    &sliceViewer2, SLOT( visibilityUpdate(int *) ) );
  QObject::connect( &labelList, SIGNAL( visibilityTable(int *) ),
                    &sliceViewer3, SLOT( visibilityUpdate(int *) ) );

  // handle mouse wheel event
  QObject::connect( &callback1, SIGNAL( wheelChanged() ),
                    &labelList, SLOT( ackWheelChanged() ) );
  QObject::connect( &callback2, SIGNAL( wheelChanged() ),
                    &labelList, SLOT( ackWheelChanged() ) );
  QObject::connect( &callback3, SIGNAL( wheelChanged() ),
                    &labelList, SLOT( ackWheelChanged() ) );

  QObject::connect( &labelList, SIGNAL( sendLabelPositions(double *) ),
                    &sliceViewer1, SLOT( wheelSlot(double *) ) );
  QObject::connect( &labelList, SIGNAL( sendLabelPositions(double *) ),
                    &sliceViewer2, SLOT( wheelSlot(double *) ) );
  QObject::connect( &labelList, SIGNAL( sendLabelPositions(double *) ),
                    &sliceViewer3, SLOT( wheelSlot(double *) ) );

  QObject::connect( &saveButton, SIGNAL( clicked() ),
                    &labelList, SLOT( saveLandmarks() ) );

  QObject::connect( &saveAsButton, SIGNAL( clicked() ),
                    &labelList, SLOT( saveAsLandmarks() ) );

  QObject::connect( &helpButton, SIGNAL( clicked() ),
                    &helpDialog, SLOT( open() ) );

  // ------------------------------------
  // Qt initialize

  base.resize(WIN_WIDTH, WIN_HEIGHT);
  base.show();
  app.exec();

  // ------------------------------------
  // clean up all VTK components

  table->Delete();

  interactor1->Delete();
  imageStyle1->Delete();
  window1->Delete();
  renderer1->Delete();
  actor1->Delete();
  color1->Delete();
  reslice1->Delete();
  resliceAxes1->Delete();

  interactor2->Delete();
  imageStyle2->Delete();
  window2->Delete();
  renderer2->Delete();
  actor2->Delete();
  color2->Delete();
  reslice2->Delete();
  resliceAxes2->Delete();

  interactor3->Delete();
  imageStyle3->Delete();
  window3->Delete();
  renderer3->Delete();
  actor3->Delete();
  color3->Delete();
  reslice3->Delete();
  resliceAxes3->Delete();

  return 0;
}
