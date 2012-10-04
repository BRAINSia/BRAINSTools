/*
 * Author: Wei Lu
 * at Psychiatry Imaging Lab,
 * University of Iowa Health Care 2010
 */

#ifndef __QVTKInteractionCallback_H
#define __QVTKInteractionCallback_H

#include "vtkCommand.h"
#include "vtkImageReslice.h"
#include "vtkRenderWindowInteractor.h"

#include <QObject>
#include <QString>

#include <cmath>

#include <QDebug>

// The mouse motion callback, to turn "Slicing" on and off
class QVTKInteractionCallback : public QObject, public vtkCommand
{
  Q_OBJECT;
public:

  QVTKInteractionCallback( const QString & text, const int type, double *bound, QObject *myParent = 0 );

  void SetImageReslice( vtkImageReslice *reslice )
  {
    m_imageReslice = reslice;
  }

  vtkImageReslice * GetImageReslice()
  {
    return m_imageReslice;
  }

  void SetInteractor( vtkRenderWindowInteractor *interactor )
  {
    m_interactor = interactor;
  }

  vtkRenderWindowInteractor * GetInteractor()
  {
    return m_interactor;
  }

  void SetDirection( double *direction )
  {
    m_direction = direction;
  }

  double * GetDirection()
  {
    return m_direction;
  }

  void SetOrigin( double *origin )
  {
    m_origin = origin;
  }

  double * GetOrigin()
  {
    return m_origin;
  }

  void SetSpacing( double *spacing )
  {
    m_spacing = spacing;
  }

  double * GetSpacing()
  {
    return m_spacing;
  }

  void SetPhysicalExtentIdentity( double *physicalExtentIdentity )
  {
    m_physicalExtentIdentity = physicalExtentIdentity;
  }

  double * GetPhysicalExtentIdentity()
  {
    return m_physicalExtentIdentity;
  }

  void SetPhysicalExtent( double *physicalExtent )
  {
    m_physicalExtent = physicalExtent;
  }

  double * GetPhysicalExtent()
  {
    return m_physicalExtent;
  }

  void SetBound( double *bound )
  {
    m_bound = bound;
  }

  double * GetBound()
  {
    return m_bound;
  }

  void SetIndexExtent( int *indexExtent )
  {
    m_indexExtent = indexExtent;
  }

  int * GetIndexExtent()
  {
    return m_indexExtent;
  }

  void SetType( int type )
  {
    m_type = type;
  }

  int GetType()
  {
    return m_type;
  }

  void SetPrecision( double precision )
  {
    m_precision = precision;
  }

  double GetPrecision()
  {
    return m_precision;
  }

  void Execute( vtkObject *, unsigned long event, void * );

  const QString & textPhysicalLocation() const
  {
    return m_text_physicalLocation;
  }

public slots:

  /*
   * This slot responds to the signal from slider bar, including update and invoke rerendering
   * Due to int ui of QSlider, we pass int zLocation = 1000 x sliceValue
   * Later we have something like double m_zLocation = zLocation / 1000
   */
  void setValueZLocation( const int zLocation );

  // this slot corresponds to commnication among slice viewers
  void setValueXLocation( const double xLocation );

  void setValueYLocation( const double yLocation );

  void receiveLabelPos( double *pos ); // set the label positon from pos info

  // from label list

  void createListAddButtonSlot(); // receive the request from add button

  void deleteListButtonSlot(); // receive the request from remove button

signals:

  void textChanged( const QString & ); // update status label

  void valueChangedZ( int ); // update slider bar, and due to int ui of QSlider

  // update neighbor slice viewer
  void valueChangedX( double );

  void valueChangedY( double );

  // update fiducial label
  void createActor();

  void switchActor();

  void deleteActor();

  void moveActor( double * );

  // update label list
  void createListItem( const QString & );

  void editListItem( const QString & );

  void switchListItem();

  void deleteListItem();

  void checkVisibility(); // check visibility request to list

  // transmit current slice number to list for visiblity check
  // designed for slider bar effect
  void visibilityCallback( double *tag );

  void wheelChanged(); // indicate the move of mouse wheel

protected:

  // Round error
  double EPS;

  // Actions:
  // Idle = 0;
  // Get location = 1;
  int m_action;

  // Pointer to vtkImageReslice
  vtkImageReslice *m_imageReslice;

  // Pointer to the interactor
  vtkRenderWindowInteractor *m_interactor;

  // Pointer to the direction
  double *m_direction;

  // Pointer to the origin
  double *m_origin;

  // Pointer to the spacing
  double *m_spacing;

  // Pointer to the physical extent of input image with identity direction
  double *m_physicalExtentIdentity;

  // Pointer to the physical extent
  double *m_physicalExtent;

  // Pointer to the physical bound ( = ordered extent )
  double *m_bound;

  // Pointer to the index extent
  int *m_indexExtent;

  // Callback type: axial, sagittal, or coronal
  int m_type;

  // Storing physical postion
  QString m_text_physicalLocation;

  // Storing z-axis postion for x channel communication
  double m_valueSendX;

  double m_valueReceiveX;

  // Storing z-axis postion for y channel communication
  double m_valueSendY;

  double m_valueReceiveY;

  // Storing z-axis postion for slice viewer / slider bar communication
  double m_valueSendZ;

  double m_valueReceiveZ;

  // Sharing last mouse position of slice viewer for slider bar control
  double m_lastPos[2];

  // initial camera postion, helping to account wheeling effect
  double m_cPos;

  // internal ratio for screen size change
  double m_r;

  // slice number precision due to int ui of QSlider
  double m_precision;
};

#endif // _QVTKInteractionCallback_H
