/*=========================================================================
 *
 *  Copyright SINAPSE: Scalable Informatics for Neuroscience, Processing and Software Engineering
 *            The University of Iowa
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/
/*
 * Author: Wei Lu
 * at Psychiatry Imaging Lab,
 * University of Iowa Health Care 2010
 */

#include "QVTKInteractionCallback.h"

#include "vtkImageData.h"
#include "vtkMatrix4x4.h"
#include "vtkCamera.h"
#include "vtkRenderer.h"
#include "vtkRendererCollection.h"
#include "vtkRenderWindow.h"
#include "vtkInteractorStyleImage.h"
#include "vtkVersion.h"

// The mouse motion callback, to turn "Slicing" on and off
void
QVTKInteractionCallback::Execute(vtkObject *, unsigned long myEvent, void *)
{
  int      type = m_type;
  double * direction = m_direction;
  double * physicalExtent = m_physicalExtent;
  double * physicalExtentIdentity = m_physicalExtentIdentity;
  double * bound = m_bound;

  /*
  double spacing[3];
  spacing[0] = m_spacing[0];
  spacing[1] = m_spacing[5];
  spacing[2] = m_spacing[10];

  // bad spacing
  if ( spacing[0] * spacing[1] * spacing[2] == 0 )
    {
    std::cerr << "Bad spacing. Halt." << std::endl;
    exit( -1 );
    }
  */

  // int *indexExtent = m_indexExtent;

  vtkRenderWindowInteractor * interactor = m_interactor;
  int                         currPos[2];

  interactor->GetEventPosition(currPos);
  int eventSize[2];
  interactor->GetSize(eventSize);
  double cPos[3]; // camera position
  interactor->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->GetActiveCamera()->GetPosition(cPos);

  // compute the screen size of image slice
  double R; // physical width/height ratio
  if (type == 2)
  { // sagittal view mode ( ASL )
    R = (bound[3] - bound[2]) / (bound[5] - bound[4]);
  }
  else if (type == 3)
  { // coronal view mode ( LSA )
    R = (bound[1] - bound[0]) / (bound[5] - bound[4]);
  }
  else
  { // axial view mode ( LAS )
    R = (bound[1] - bound[0]) / (bound[3] - bound[2]);
  }
  int h = 0.96 * eventSize[1] / sqrt(1.0 + R * R); // sreen height of image
                                                   // slice
  int w = R * h;                                   // screen width of image
                                                   // slice

  // wheeling effect
  if (m_cPos < 0.0)
  {
    m_cPos = cPos[2];
  }
  else if (m_cPos != cPos[2])
  {
    double r = m_cPos / cPos[2];
    h *= r;
    w *= r;
    if (m_r != r)
    {
      m_r = r;
      emit wheelChanged();
    }
  }

  // center mapping due to image direction compensation for the reslice axes
  int    centerMapped;
  double centerMappedDouble;
  if (type == 2)
  { // sagittal view mode ( ASL )
    centerMappedDouble = (direction[0] + 2 * direction[4] + 4 * direction[8]) / 2;
  }
  else if (type == 3)
  { // coronal view mode ( LSA )
    centerMappedDouble = (direction[1] + 2 * direction[5] + 4 * direction[9]) / 2;
  }
  else
  { // axial view mode ( LAS )
    centerMappedDouble = (direction[2] + 2 * direction[6] + 4 * direction[10]) / 2;
  }
  if (centerMappedDouble < 0)
  {
    centerMappedDouble *= -1;
  }
  centerMapped = centerMappedDouble + EPS;

  // move the center point that we are slicing through
  vtkImageReslice * reslice = m_imageReslice;
#if (VTK_MAJOR_VERSION < 6)
  reslice->GetOutput()->UpdateInformation();
#else
  reslice->UpdateInformation();
#endif
  vtkMatrix4x4 * matrix = reslice->GetResliceAxes();
  double         zPhysicalLocation; // pyhysical position along normal vector
  double         center[4] = { (*matrix)[0][3], (*matrix)[1][3], (*matrix)[2][3], 1.0 };
  QString        labelPosString;   // label postion for viewer-list communication
  double         labelPosRatio[3]; // label postion for inter-viewer
                                   // communication

  // handling keyboard myEvent especially for fiducial labeling
  if (myEvent == vtkCommand::KeyPressEvent)
  {
    char key = m_interactor->GetKeyCode();
    if (key == 'C' || key == 'c') // create a new actor
    {
      QString initPos = initPos.number((physicalExtent[0] + physicalExtent[1]) / 2.0) + " " +
                        initPos.number((physicalExtent[2] + physicalExtent[3]) / 2.0) + " " +
                        initPos.number((physicalExtent[4] + physicalExtent[5]) / 2.0);
      emit createListItem(initPos);
      emit createActor();
    }
    if (key == 'S' || key == 's') // switch from actors
    {
      emit switchListItem();
      emit switchActor();
    }
    if (key == 'D' || key == 'd') // delete current actor
    {
      emit deleteListItem();
      emit deleteActor();
    }
    emit checkVisibility();
  }

  // handling mouse myEvent
  else if (myEvent == vtkCommand::LeftButtonPressEvent)
  {
    m_action = 1;

    // Sharing last mouse position of slice viewer for slider bar control
    m_lastPos[0] = currPos[0];
    m_lastPos[1] = currPos[1];

    // linear mapping
    double zFactor = (center[centerMapped] - physicalExtentIdentity[2 * centerMapped]) /
                     (physicalExtentIdentity[2 * centerMapped + 1] - physicalExtentIdentity[2 * centerMapped]);

    // show current physical position
    if (type == 2)
    { // sagittal view mode ( ASL )
      zPhysicalLocation = zFactor * (physicalExtent[1] - physicalExtent[0]) + physicalExtent[0];
      if (zPhysicalLocation < bound[0])
      {
        zPhysicalLocation = bound[0];
      }
      if (zPhysicalLocation > bound[1])
      {
        zPhysicalLocation = bound[1];
      }

      // store the physical location
      m_valueSendX = (bound[3] - bound[2]) * (currPos[0] - (eventSize[0] - w) / 2.0) / w + bound[2];
      m_valueSendY = (bound[5] - bound[4]) * (currPos[1] - (eventSize[1] - h) / 2.0) / h + bound[4];
      m_valueReceiveZ = zPhysicalLocation;

      /*
      // regulation due to resolution
      m_valueSendX = spacing[1-m_type%2] * round( m_valueSendX / spacing[1-m_type%2] );
      m_valueSendY = spacing[m_type/2+1] * round( m_valueSendY / spacing[m_type/2+1] );
      m_valueReceiveZ = spacing[(m_type+1)%3] * round( m_valueReceiveZ / spacing[(m_type+1)%3] );
      */

      m_text_physicalLocation = m_text_physicalLocation.number(m_valueReceiveZ);

      // label postion for list-viewer communication
      labelPosString = m_text_physicalLocation.number(m_valueReceiveZ) + " " +
                       m_text_physicalLocation.number(m_valueSendX) + " " +
                       m_text_physicalLocation.number(m_valueSendY);

      // label postion ratio for inter-viewer communication
      labelPosRatio[0] = (m_valueReceiveZ - bound[0]) / (bound[1] - bound[0]);
      labelPosRatio[1] = (m_valueSendX - bound[2]) / (bound[3] - bound[2]);
      labelPosRatio[2] = (m_valueSendY - bound[4]) / (bound[5] - bound[4]);
    }
    else if (type == 3)
    { // coronal view mode ( LSA )
      zPhysicalLocation = (1 - zFactor) * (physicalExtent[2] - physicalExtent[3]) + physicalExtent[3];
      if (zPhysicalLocation < bound[2])
      {
        zPhysicalLocation = bound[2];
      }
      if (zPhysicalLocation > bound[3])
      {
        zPhysicalLocation = bound[3];
      }

      // store the physical location
      m_valueSendX = (bound[1] - bound[0]) * (currPos[0] - (eventSize[0] - w) / 2.0) / w + bound[0];
      m_valueSendY = (bound[5] - bound[4]) * (currPos[1] - (eventSize[1] - h) / 2.0) / h + bound[4];
      m_valueReceiveZ = zPhysicalLocation;

      /*
      // regulation due to resolution
      m_valueSendX = spacing[1-m_type%2] * round( m_valueSendX / spacing[1-m_type%2] );
      m_valueSendY = spacing[m_type/2+1] * round( m_valueSendY / spacing[m_type/2+1] );
      m_valueReceiveZ = spacing[(m_type+1)%3] * round( m_valueReceiveZ / spacing[(m_type+1)%3] );
      */

      m_text_physicalLocation = m_text_physicalLocation.number(m_valueReceiveZ);

      // label postion for list-viewer communication
      labelPosString = m_text_physicalLocation.number(m_valueSendX) + " " +
                       m_text_physicalLocation.number(m_valueReceiveZ) + " " +
                       m_text_physicalLocation.number(m_valueSendY);

      // label postion ratio for inter-viewer communication
      labelPosRatio[0] = (m_valueSendX - bound[0]) / (bound[1] - bound[0]);
      labelPosRatio[1] = (m_valueReceiveZ - bound[2]) / (bound[3] - bound[2]);
      labelPosRatio[2] = (m_valueSendY - bound[4]) / (bound[5] - bound[4]);
    }
    else
    { // axial view mode ( LAS )
      zPhysicalLocation = zFactor * (physicalExtent[5] - physicalExtent[4]) + physicalExtent[4];
      if (zPhysicalLocation < bound[4])
      {
        zPhysicalLocation = bound[4];
      }
      if (zPhysicalLocation > bound[5])
      {
        zPhysicalLocation = bound[5];
      }

      // store the physical location
      m_valueSendX = (bound[1] - bound[0]) * (currPos[0] - (eventSize[0] - w) / 2.0) / w + bound[0];
      m_valueSendY = (bound[2] - bound[3]) * (currPos[1] - (eventSize[1] - h) / 2.0) / h + bound[3];
      m_valueReceiveZ = zPhysicalLocation;

      /*
      // regulation due to resolution
      m_valueSendX = spacing[1-m_type%2] * round( m_valueSendX / spacing[1-m_type%2] );
      m_valueSendY = spacing[m_type/2+1] * round( m_valueSendY / spacing[m_type/2+1] );
      m_valueReceiveZ = spacing[(m_type+1)%3] * round( m_valueReceiveZ / spacing[(m_type+1)%3] );
      */

      m_text_physicalLocation = m_text_physicalLocation.number(m_valueReceiveZ);

      // label postion for list-viewer communication
      labelPosString = m_text_physicalLocation.number(m_valueSendX) + " " +
                       m_text_physicalLocation.number(m_valueSendY) + " " +
                       m_text_physicalLocation.number(m_valueReceiveZ);

      // label postion ratio for inter-viewer communication
      labelPosRatio[0] = (m_valueSendX - bound[0]) / (bound[1] - bound[0]);
      labelPosRatio[1] = (m_valueSendY - bound[2]) / (bound[3] - bound[2]);
      labelPosRatio[2] = (m_valueReceiveZ - bound[4]) / (bound[5] - bound[4]);
    }
    // signal to fiducial label
    emit moveActor(labelPosRatio);

    // signal to neighbor slice viewers
    emit valueChangedX(m_valueSendX);
    emit valueChangedY(m_valueSendY);

    // signal to slider bar
    emit valueChangedZ(m_valueReceiveZ * m_precision);

    /*
    // signal to neighbor slice viewers
    emit valueChangedX( round( m_valueSendX / spacing[1-m_type%2] ) );
    emit valueChangedY( round( m_valueSendY / spacing[ m_type / 2 + 1 ] ) );

    // signal to slider bar
    emit valueChangedZ( round( m_valueReceiveZ / spacing[ (m_type + 1)  % 3 ] ) );
    */

    // signal to status label
    emit textChanged(m_text_physicalLocation);

    // signal to label list
    emit editListItem(labelPosString);

    emit checkVisibility();
  }
  else if (myEvent == vtkCommand::LeftButtonReleaseEvent)
  {
    m_action = 0;
  }
  else if (myEvent == vtkCommand::MouseMoveEvent)
  {
    if (m_action == 1)
    { // change a slice to show
      // Sharing last mouse position of slice viewer for slider bar control
      m_lastPos[0] = currPos[0];
      m_lastPos[1] = currPos[1];

      // linear mapping
      double zFactor = (center[centerMapped] - physicalExtentIdentity[2 * centerMapped]) /
                       (physicalExtentIdentity[2 * centerMapped + 1] - physicalExtentIdentity[2 * centerMapped]);

      // show current physical position
      if (type == 2)
      { // sagittal view mode ( ASL )
        zPhysicalLocation = zFactor * (physicalExtent[1] - physicalExtent[0]) + physicalExtent[0];
        if (zPhysicalLocation < bound[0])
        {
          zPhysicalLocation = bound[0];
        }
        if (zPhysicalLocation > bound[1])
        {
          zPhysicalLocation = bound[1];
        }

        // store the physical location
        m_valueSendX = (bound[3] - bound[2]) * (currPos[0] - (eventSize[0] - w) / 2.0) / w + bound[2];
        m_valueSendY = (bound[5] - bound[4]) * (currPos[1] - (eventSize[1] - h) / 2.0) / h + bound[4];
        m_valueReceiveZ = zPhysicalLocation;

        /*
        // regulation due to resolution
        m_valueSendX = spacing[1-m_type%2] * round( m_valueSendX / spacing[1-m_type%2] );
        m_valueSendY = spacing[m_type/2+1] * round( m_valueSendY / spacing[m_type/2+1] );
        m_valueReceiveZ = spacing[(m_type+1)%3] * round( m_valueReceiveZ / spacing[(m_type+1)%3] );
        */

        m_text_physicalLocation = m_text_physicalLocation.number(m_valueReceiveZ);

        // label postion for list-viewer communication
        labelPosString = m_text_physicalLocation.number(m_valueReceiveZ) + " " +
                         m_text_physicalLocation.number(m_valueSendX) + " " +
                         m_text_physicalLocation.number(m_valueSendY);

        // label postion ratio for inter-viewer communication
        labelPosRatio[0] = (m_valueReceiveZ - bound[0]) / (bound[1] - bound[0]);
        labelPosRatio[1] = (m_valueSendX - bound[2]) / (bound[3] - bound[2]);
        labelPosRatio[2] = (m_valueSendY - bound[4]) / (bound[5] - bound[4]);
      }
      else if (type == 3)
      { // coronal view mode ( LSA )
        zPhysicalLocation = (1 - zFactor) * (physicalExtent[2] - physicalExtent[3]) + physicalExtent[3];
        if (zPhysicalLocation < bound[2])
        {
          zPhysicalLocation = bound[2];
        }
        if (zPhysicalLocation > bound[3])
        {
          zPhysicalLocation = bound[3];
        }

        // store the physical location
        m_valueSendX = (bound[1] - bound[0]) * (currPos[0] - (eventSize[0] - w) / 2.0) / w + bound[0];
        m_valueSendY = (bound[5] - bound[4]) * (currPos[1] - (eventSize[1] - h) / 2.0) / h + bound[4];
        m_valueReceiveZ = zPhysicalLocation;

        /*
        // regulation due to resolution
        m_valueSendX = spacing[1-m_type%2] * round( m_valueSendX / spacing[1-m_type%2] );
        m_valueSendY = spacing[m_type/2+1] * round( m_valueSendY / spacing[m_type/2+1] );
        m_valueReceiveZ = spacing[(m_type+1)%3] * round( m_valueReceiveZ / spacing[(m_type+1)%3] );
        */

        m_text_physicalLocation = m_text_physicalLocation.number(m_valueReceiveZ);

        // label postion for list-viewer communication
        labelPosString = m_text_physicalLocation.number(m_valueSendX) + " " +
                         m_text_physicalLocation.number(m_valueReceiveZ) + " " +
                         m_text_physicalLocation.number(m_valueSendY);

        // label postion ratio for inter-viewer communication
        labelPosRatio[0] = (m_valueSendX - bound[0]) / (bound[1] - bound[0]);
        labelPosRatio[1] = (m_valueReceiveZ - bound[2]) / (bound[3] - bound[2]);
        labelPosRatio[2] = (m_valueSendY - bound[4]) / (bound[5] - bound[4]);
      }
      else
      { // axial view mode ( LAS )
        zPhysicalLocation = zFactor * (physicalExtent[5] - physicalExtent[4]) + physicalExtent[4];
        if (zPhysicalLocation < bound[4])
        {
          zPhysicalLocation = bound[4];
        }
        if (zPhysicalLocation > bound[5])
        {
          zPhysicalLocation = bound[5];
        }

        // store the physical location
        m_valueSendX = (bound[1] - bound[0]) * (currPos[0] - (eventSize[0] - w) / 2.0) / w + bound[0];
        m_valueSendY = (bound[2] - bound[3]) * (currPos[1] - (eventSize[1] - h) / 2.0) / h + bound[3];
        m_valueReceiveZ = zPhysicalLocation;

        /*
        // regulation due to resolution
        m_valueSendX = spacing[1-m_type%2] * round( m_valueSendX / spacing[1-m_type%2] );
        m_valueSendY = spacing[m_type/2+1] * round( m_valueSendY / spacing[m_type/2+1] );
        m_valueReceiveZ = spacing[(m_type+1)%3] * round( m_valueReceiveZ / spacing[(m_type+1)%3] );
        */

        m_text_physicalLocation = m_text_physicalLocation.number(m_valueReceiveZ);

        // label postion for list-viewer communication
        labelPosString = m_text_physicalLocation.number(m_valueSendX) + " " +
                         m_text_physicalLocation.number(m_valueSendY) + " " +
                         m_text_physicalLocation.number(m_valueReceiveZ);

        // label postion ratio for inter-viewer communication
        labelPosRatio[0] = (m_valueSendX - bound[0]) / (bound[1] - bound[0]);
        labelPosRatio[1] = (m_valueSendY - bound[2]) / (bound[3] - bound[2]);
        labelPosRatio[2] = (m_valueReceiveZ - bound[4]) / (bound[5] - bound[4]);
      }

      // signal to fiducial label
      emit moveActor(labelPosRatio);

      // signal to neighbor slice viewers
      emit valueChangedX(m_valueSendX);
      emit valueChangedY(m_valueSendY);

      // signal to slider bar
      emit valueChangedZ(m_valueReceiveZ * m_precision);

      /*
       // signal to neighbor slice viewers
       emit valueChangedX( round( m_valueSendX / spacing[1-m_type%2] ) );
       emit valueChangedY( round( m_valueSendY / spacing[ m_type / 2 + 1 ] ) );

       // signal to slider bar
       emit valueChangedZ( round( m_valueReceiveZ / spacing[ (m_type + 1)  % 3 ] ) );
       */

      // signal to status label
      emit textChanged(m_text_physicalLocation);

      // signal to label list
      emit editListItem(labelPosString);
      emit checkVisibility();
    }
    else
    {
      vtkInteractorStyle * style = vtkInteractorStyle::SafeDownCast(interactor->GetInteractorStyle());
      if (style)
      {
        style->OnMouseMove();
      }
    }
  }
}

QVTKInteractionCallback::QVTKInteractionCallback(const QString & text,
                                                 const int       type,
                                                 double *        bound,
                                                 QObject *       myParent)
  : QObject(myParent)
  , m_origin(nullptr)
  , m_valueSendX(0)
  , m_valueSendY(0)
  , m_valueSendZ(0)
{
  EPS = .3;
  m_action = 0;
  m_imageReslice = nullptr;
  m_interactor = nullptr;
  m_direction = nullptr;
  m_physicalExtentIdentity = nullptr;
  m_physicalExtent = nullptr;
  m_bound = bound;
  m_spacing = nullptr;
  m_indexExtent = nullptr;
  m_text_physicalLocation = text;
  m_type = type;
  m_lastPos[0] = 0;
  m_lastPos[1] = 0;
  m_cPos = -1.0;
  m_r = 0;
  m_valueReceiveX = -1;
  m_valueReceiveY = -1;
  m_precision = 1000.;

  if (m_type == 2)
  { // sagittal view mode ( ASL )
    m_valueReceiveZ = (m_bound[0] + m_bound[1]) / 2;
  }
  else if (m_type == 3)
  { // coronal view mode ( LSA )
    m_valueReceiveZ = (m_bound[2] + m_bound[3]) / 2;
  }
  else
  { // axial view mode ( LAS )
    m_valueReceiveZ = (m_bound[4] + m_bound[5]) / 2;
  }
}

// slot function for recieve signal from slider bar
void
QVTKInteractionCallback::setValueZLocation(const int zLocation)
{
  if ((m_valueReceiveZ >= zLocation / m_precision - 0.5) && (m_valueReceiveZ < zLocation / m_precision + 0.5))
  {
    return;
  }
  m_valueReceiveZ = zLocation / m_precision;

  /*
  double spacing[3];
  spacing[0] = m_spacing[0];
  spacing[1] = m_spacing[5];
  spacing[2] = m_spacing[10];

  // avoid infinite loop and update internal z-axis location
  if ( m_valueReceiveZ == zLocation * spacing[ (m_type + 1)  % 3 ] )
    {
    return;
    }
  m_valueReceiveZ = zLocation * spacing[ (m_type + 1)  % 3 ];
  */

#if (VTK_MAJOR_VERSION < 6)
  m_imageReslice->GetOutput()->UpdateInformation();
#else
  m_imageReslice->UpdateInformation();
#endif
  vtkMatrix4x4 * matrix = m_imageReslice->GetResliceAxes();

  // center mapping due to image direction compensation for the reslice axes
  int      centerMapped;
  double   centerMappedDouble;
  int      type = m_type;
  double * direction = m_direction;
  double * physicalExtent = m_physicalExtent;
  double * physicalExtentIdentity = m_physicalExtentIdentity;

  if (type == 2)
  { // sagittal view mode ( ASL )
    centerMappedDouble = (direction[0] + 2 * direction[4] + 4 * direction[8]) / 2;
  }
  else if (type == 3)
  { // coronal view mode ( LSA )
    centerMappedDouble = (direction[1] + 2 * direction[5] + 4 * direction[9]) / 2;
  }
  else
  { // axial view mode ( LAS )
    centerMappedDouble = (direction[2] + 2 * direction[6] + 4 * direction[10]) / 2;
  }
  if (centerMappedDouble < 0)
  {
    centerMappedDouble *= -1;
  }
  centerMapped = centerMappedDouble + EPS;

  double * bound = m_bound;
  double   center;
  int      eventSize[2];
  m_interactor->GetSize(eventSize);

  // show current physical position
  if (type == 2)
  { // sagittal view mode ( ASL )
    if (m_valueReceiveZ < bound[0])
    {
      m_valueReceiveZ = bound[0];
    }
    if (m_valueReceiveZ > bound[1])
    {
      m_valueReceiveZ = bound[1];
    }

    // calculate the center value from m_valueReceiveZ
    center = (m_valueReceiveZ - physicalExtent[0]) / (physicalExtent[1] - physicalExtent[0]) *
               (physicalExtentIdentity[2 * centerMapped + 1] - physicalExtentIdentity[2 * centerMapped]) +
             physicalExtentIdentity[2 * centerMapped];

    // store the physical location and the z-axis location
    m_text_physicalLocation = m_text_physicalLocation.number(m_valueReceiveZ);
  }
  else if (type == 3)
  { // coronal view mode ( LSA )
    if (m_valueReceiveZ < bound[2])
    {
      m_valueReceiveZ = bound[2];
    }
    if (m_valueReceiveZ > bound[3])
    {
      m_valueReceiveZ = bound[3];
    }

    // calculate the center value from m_valueReceiveZ
    center = physicalExtentIdentity[2 * centerMapped + 1] -
             (m_valueReceiveZ - physicalExtent[3]) / (physicalExtent[2] - physicalExtent[3]) *
               (physicalExtentIdentity[2 * centerMapped + 1] - physicalExtentIdentity[2 * centerMapped]);

    // store the physical location and the z-axis location
    m_text_physicalLocation = m_text_physicalLocation.number(m_valueReceiveZ);
  }
  else
  { // axial view mode ( LAS )
    if (m_valueReceiveZ < bound[4])
    {
      m_valueReceiveZ = bound[4];
    }
    if (m_valueReceiveZ > bound[5])
    {
      m_valueReceiveZ = bound[5];
    }

    // calculate the center value from m_valueReceiveZ
    center = (m_valueReceiveZ - physicalExtent[4]) / (physicalExtent[5] - physicalExtent[4]) *
               (physicalExtentIdentity[2 * centerMapped + 1] - physicalExtentIdentity[2 * centerMapped]) +
             physicalExtentIdentity[2 * centerMapped];

    // store the physical location and the z-axis location
    m_text_physicalLocation = m_text_physicalLocation.number(m_valueReceiveZ);
  }

  matrix->SetElement(centerMapped, 3, center);

  m_interactor->Render();

  // signal to status label
  emit textChanged(m_text_physicalLocation);

  // transmit current slice number to list for visiblity check
  // designed for slider bar effect
  double tag[2];
  tag[0] = m_valueReceiveZ;
  tag[1] = double(m_type);

  emit visibilityCallback(tag);
}

// this slot function will change the slice according to signal from x-channel
// slice viewer
void
QVTKInteractionCallback::setValueXLocation(const double xLocation)
{
  // avoid infinite loop and update internal z-axis location
  if (m_valueReceiveX == xLocation)
  {
    return;
  }
  m_valueReceiveX = xLocation;

  /*
  double spacing[3];
  spacing[0] = m_spacing[0];
  spacing[1] = m_spacing[5];
  spacing[2] = m_spacing[10];

  // avoid infinite loop and update internal z-axis location
  if ( m_valueReceiveX == xLocation * spacing[1-m_type%2] )
    {
    return;
    }
  m_valueReceiveX = xLocation * spacing[1-m_type%2];
  */

#if (VTK_MAJOR_VERSION < 6)
  m_imageReslice->GetOutput()->UpdateInformation();
#else
  m_imageReslice->UpdateInformation();
#endif
  vtkMatrix4x4 * matrix = m_imageReslice->GetResliceAxes();

  // center mapping due to image direction compensation for the reslice axes
  int      centerMapped;
  double   centerMappedDouble;
  int      type = m_type;
  double * direction = m_direction;
  double * physicalExtent = m_physicalExtent;
  double * physicalExtentIdentity = m_physicalExtentIdentity;

  if (type == 2)
  { // sagittal view mode ( ASL )
    centerMappedDouble = (direction[0] + 2 * direction[4] + 4 * direction[8]) / 2;
  }
  else if (type == 3)
  { // coronal view mode ( LSA )
    centerMappedDouble = (direction[1] + 2 * direction[5] + 4 * direction[9]) / 2;
  }
  else
  { // axial view mode ( LAS )
    centerMappedDouble = (direction[2] + 2 * direction[6] + 4 * direction[10]) / 2;
  }
  if (centerMappedDouble < 0)
  {
    centerMappedDouble *= -1;
  }
  centerMapped = centerMappedDouble + EPS;

  double * bound = m_bound;
  double   center;
  int      eventSize[2];
  m_interactor->GetSize(eventSize);

  // show current physical position
  if (type == 2)
  { // sagittal view mode ( ASL )
    if (m_valueReceiveX < bound[0])
    {
      m_valueReceiveX = bound[0];
    }
    if (m_valueReceiveX > bound[1])
    {
      m_valueReceiveX = bound[1];
    }

    // calculate the center value from m_valueReceiveX
    center = (m_valueReceiveX - physicalExtent[0]) / (physicalExtent[1] - physicalExtent[0]) *
               (physicalExtentIdentity[2 * centerMapped + 1] - physicalExtentIdentity[2 * centerMapped]) +
             physicalExtentIdentity[2 * centerMapped];

    // store the physical location and the z-axis location
    m_text_physicalLocation = m_text_physicalLocation.number(m_valueReceiveX);
  }
  else if (type == 3)
  { // coronal view mode ( LSA )
    if (m_valueReceiveX < bound[2])
    {
      m_valueReceiveX = bound[2];
    }
    if (m_valueReceiveX > bound[3])
    {
      m_valueReceiveX = bound[3];
    }

    // calculate the center value from m_valueReceiveX
    center = physicalExtentIdentity[2 * centerMapped + 1] -
             (m_valueReceiveX - physicalExtent[3]) / (physicalExtent[2] - physicalExtent[3]) *
               (physicalExtentIdentity[2 * centerMapped + 1] - physicalExtentIdentity[2 * centerMapped]);

    // store the physical location and the z-axis location
    m_text_physicalLocation = m_text_physicalLocation.number(m_valueReceiveX);
  }
  else
  { // axial view mode ( LAS )
    if (m_valueReceiveX < bound[4])
    {
      m_valueReceiveX = bound[4];
    }
    if (m_valueReceiveX > bound[5])
    {
      m_valueReceiveX = bound[5];
    }

    // calculate the center value from m_valueReceiveX
    center = (m_valueReceiveX - physicalExtent[4]) / (physicalExtent[5] - physicalExtent[4]) *
               (physicalExtentIdentity[2 * centerMapped + 1] - physicalExtentIdentity[2 * centerMapped]) +
             physicalExtentIdentity[2 * centerMapped];

    // store the physical location and the z-axis location
    m_text_physicalLocation = m_text_physicalLocation.number(m_valueReceiveX);
  }

  matrix->SetElement(centerMapped, 3, center);

  m_interactor->Render();

  // signal to status label and slider bar
  emit textChanged(m_text_physicalLocation);

  emit valueChangedZ(m_valueReceiveX * m_precision);

  // emit valueChangedZ( round( m_valueReceiveX / spacing[1-m_type%2] ) );
}

// this slot function will change the slice according to signal from y-channel
// slice viewer
void
QVTKInteractionCallback::setValueYLocation(const double yLocation)
{
  // avoid infinite loop and update internal z-axis location
  if (m_valueReceiveY == yLocation)
  {
    return;
  }
  m_valueReceiveY = yLocation;

  /*
  double spacing[3];
  spacing[0] = m_spacing[0];
  spacing[1] = m_spacing[5];
  spacing[2] = m_spacing[10];

  // avoid infinite loop and update internal z-axis location
  if ( m_valueReceiveY == yLocation * spacing[ m_type / 2 + 1 ] )
    {
    return;
    }
  m_valueReceiveY = yLocation * spacing[ m_type / 2 + 1 ];
  */

#if (VTK_MAJOR_VERSION < 6)
  m_imageReslice->GetOutput()->UpdateInformation();
#else
  m_imageReslice->UpdateInformation();
#endif
  vtkMatrix4x4 * matrix = m_imageReslice->GetResliceAxes();

  // center mapping due to image direction compensation for the reslice axes
  int      centerMapped;
  double   centerMappedDouble;
  int      type = m_type;
  double * direction = m_direction;
  double * physicalExtent = m_physicalExtent;
  double * physicalExtentIdentity = m_physicalExtentIdentity;

  if (type == 2)
  { // sagittal view mode ( ASL )
    centerMappedDouble = (direction[0] + 2 * direction[4] + 4 * direction[8]) / 2;
  }
  else if (type == 3)
  { // coronal view mode ( LSA )
    centerMappedDouble = (direction[1] + 2 * direction[5] + 4 * direction[9]) / 2;
  }
  else
  { // axial view mode ( LAS )
    centerMappedDouble = (direction[2] + 2 * direction[6] + 4 * direction[10]) / 2;
  }
  if (centerMappedDouble < 0)
  {
    centerMappedDouble *= -1;
  }
  centerMapped = centerMappedDouble + EPS;

  double * bound = m_bound;
  double   center;
  int      eventSize[2];
  m_interactor->GetSize(eventSize);

  // show current physical position
  if (type == 2)
  { // sagittal view mode ( ASL )
    if (m_valueReceiveY < bound[0])
    {
      m_valueReceiveY = bound[0];
    }
    if (m_valueReceiveY > bound[1])
    {
      m_valueReceiveY = bound[1];
    }

    // calculate the center value from m_valueReceiveY
    center = (m_valueReceiveY - physicalExtent[0]) / (physicalExtent[1] - physicalExtent[0]) *
               (physicalExtentIdentity[2 * centerMapped + 1] - physicalExtentIdentity[2 * centerMapped]) +
             physicalExtentIdentity[2 * centerMapped];

    // store the physical location and the z-axis location
    m_text_physicalLocation = m_text_physicalLocation.number(m_valueReceiveY);
  }
  else if (type == 3)
  { // coronal view mode ( LSA )
    if (m_valueReceiveY < bound[2])
    {
      m_valueReceiveY = bound[2];
    }
    if (m_valueReceiveY > bound[3])
    {
      m_valueReceiveY = bound[3];
    }

    // calculate the center value from m_valueReceiveY
    center = physicalExtentIdentity[2 * centerMapped + 1] -
             (m_valueReceiveY - physicalExtent[3]) / (physicalExtent[2] - physicalExtent[3]) *
               (physicalExtentIdentity[2 * centerMapped + 1] - physicalExtentIdentity[2 * centerMapped]);

    // store the physical location and the z-axis location
    m_text_physicalLocation = m_text_physicalLocation.number(m_valueReceiveY);
  }
  else
  { // axial view mode ( LAS )
    if (m_valueReceiveY < bound[4])
    {
      m_valueReceiveY = bound[4];
    }
    if (m_valueReceiveY > bound[5])
    {
      m_valueReceiveY = bound[5];
    }

    // calculate the center value from m_valueReceiveY
    center = (m_valueReceiveY - physicalExtent[4]) / (physicalExtent[5] - physicalExtent[4]) *
               (physicalExtentIdentity[2 * centerMapped + 1] - physicalExtentIdentity[2 * centerMapped]) +
             physicalExtentIdentity[2 * centerMapped];

    // store the physical location and the z-axis location
    m_text_physicalLocation = m_text_physicalLocation.number(m_valueReceiveY);
  }

  matrix->SetElement(centerMapped, 3, center);

  m_interactor->Render();

  // signal to status label and slider bar
  emit textChanged(m_text_physicalLocation);
  emit valueChangedZ(m_valueReceiveY * m_precision);

  // emit valueChangedZ( round( m_valueReceiveY / spacing[ m_type / 2 + 1 ] ) );
}

void
QVTKInteractionCallback::receiveLabelPos(double * pos)
{
#if (VTK_MAJOR_VERSION < 6)
  m_imageReslice->GetOutput()->UpdateInformation();
#else
  m_imageReslice->UpdateInformation();
#endif
  vtkMatrix4x4 * matrix = m_imageReslice->GetResliceAxes();

  // center mapping due to image direction compensation for the reslice axes
  int      centerMapped;
  double   centerMappedDouble;
  int      type = m_type;
  double * direction = m_direction;
  double * physicalExtent = m_physicalExtent;
  double * physicalExtentIdentity = m_physicalExtentIdentity;

  if (type == 2)
  { // sagittal view mode ( ASL )
    centerMappedDouble = (direction[0] + 2 * direction[4] + 4 * direction[8]) / 2;
  }
  else if (type == 3)
  { // coronal view mode ( LSA )
    centerMappedDouble = (direction[1] + 2 * direction[5] + 4 * direction[9]) / 2;
  }
  else
  { // axial view mode ( LAS )
    centerMappedDouble = (direction[2] + 2 * direction[6] + 4 * direction[10]) / 2;
  }
  if (centerMappedDouble < 0)
  {
    centerMappedDouble *= -1;
  }
  centerMapped = centerMappedDouble + EPS;

  double * bound = m_bound;
  double   center;
  double   posLocal;
  int      eventSize[2];
  m_interactor->GetSize(eventSize);

  // show current physical position
  if (type == 2)
  { // sagittal view mode ( ASL )
    posLocal = pos[0];
    if (pos[0] < bound[0])
    {
      pos[0] = bound[0];
    }
    if (pos[0] > bound[1])
    {
      pos[0] = bound[1];
    }

    // calculate the center value from pos[0]
    center = (pos[0] - physicalExtent[0]) / (physicalExtent[1] - physicalExtent[0]) *
               (physicalExtentIdentity[2 * centerMapped + 1] - physicalExtentIdentity[2 * centerMapped]) +
             physicalExtentIdentity[2 * centerMapped];

    // store the physical location and the z-axis location
    m_text_physicalLocation = m_text_physicalLocation.number(pos[0]);
  }
  else if (type == 3)
  { // coronal view mode ( LSA )
    posLocal = pos[1];
    if (pos[1] < bound[2])
    {
      pos[1] = bound[2];
    }
    if (pos[1] > bound[3])
    {
      pos[1] = bound[3];
    }

    // calculate the center value from pos[1]
    center = physicalExtentIdentity[2 * centerMapped + 1] -
             (pos[1] - physicalExtent[3]) / (physicalExtent[2] - physicalExtent[3]) *
               (physicalExtentIdentity[2 * centerMapped + 1] - physicalExtentIdentity[2 * centerMapped]);

    // store the physical location and the z-axis location
    m_text_physicalLocation = m_text_physicalLocation.number(pos[1]);
  }
  else
  { // axial view mode ( LAS )
    posLocal = pos[2];
    if (pos[2] < bound[4])
    {
      pos[2] = bound[4];
    }
    if (pos[2] > bound[5])
    {
      pos[2] = bound[5];
    }

    // calculate the center value from pos[2]
    center = (pos[2] - physicalExtent[4]) / (physicalExtent[5] - physicalExtent[4]) *
               (physicalExtentIdentity[2 * centerMapped + 1] - physicalExtentIdentity[2 * centerMapped]) +
             physicalExtentIdentity[2 * centerMapped];

    // store the physical location and the z-axis location
    m_text_physicalLocation = m_text_physicalLocation.number(pos[2]);
  }

  matrix->SetElement(centerMapped, 3, center);

  m_interactor->Render();

  // signal to status label and slider bar
  emit textChanged(m_text_physicalLocation);
  emit valueChangedZ(posLocal * m_precision);
}

void
QVTKInteractionCallback::createListAddButtonSlot()
{
  double * physicalExtent = m_physicalExtent;
  QString  initPos = initPos.number((physicalExtent[0] + physicalExtent[1]) / 2.0) + " " +
                    initPos.number((physicalExtent[2] + physicalExtent[3]) / 2.0) + " " +
                    initPos.number((physicalExtent[4] + physicalExtent[5]) / 2.0);
  emit createListItem(initPos);
}

// void QVTKInteractionCallback::deleteListButtonSlot()
// {
// }
