/*
 * Author: Wei Lu
 * at Psychiatry Imaging Lab,
 * University of Iowa Health Care 2010
 */

#include "QLabelList.h"
#include "itkNumberToString.h"

void QLabelList::createListItemSlot(const QString & label)
{
  // copy label to m_label
  m_label.clear();
  m_label = m_label.append(label);

  QString newLabel = QString("newLabel_%1: ( %2 )").arg( QString::number(m_color + 1) ).arg(label);
  this->addItem(newLabel);
  this->setCurrentRow(this->count() - 1);
  this->setCurrentItem( this->item( this->currentRow() ) );
  QColor color( ( m_color % 3 / 2 / 2.0 + ( m_color * m_color ) % 3 / 4.0 ) * 255,
                ( ( m_color + 1 ) % 3 / 2 / 2.0 + m_color % 7 / 12.0 ) * 255,
                ( ( m_color + 2 ) % 3 / 2 / 2.0 + m_color % 4 / 6.0 ) * 255,
                127 );
  QBrush brush(color);
  this->item( this->currentRow() )->setBackground(brush);

  cancelHighlight( this->currentItem() );

  m_color++;
}

void QLabelList::createListItem(const QString & label, const QString & name)
{
  // copy label to m_label
  m_label.clear();
  m_label = m_label.append(label);

  QString newLabel = QString("%1: ( %2 )").arg(name).arg(label);
  this->addItem(newLabel);
  this->setCurrentRow(this->count() - 1);
  this->setCurrentItem( this->item( this->currentRow() ) );
  QColor color( ( m_color % 3 / 2 / 2.0 + ( m_color * m_color ) % 3 / 4.0 ) * 255,
                ( ( m_color + 1 ) % 3 / 2 / 2.0 + m_color % 7 / 12.0 ) * 255,
                ( ( m_color + 2 ) % 3 / 2 / 2.0 + m_color % 4 / 6.0 ) * 255,
                127 );
  QBrush brush(color);
  this->item( this->currentRow() )->setBackground(brush);

  cancelHighlight( this->currentItem() );

  m_color++;
}

void QLabelList::createListItemAddButtonSlot()
{
  createListItemSlot(m_label);
}

void QLabelList::editListItemSlot(const QString & label)
{
  if( this->currentItem() != NULL )
    {
    QString lastLabel = this->currentItem()->text();
    QString newLabel = QString("%1: ( %2 )").arg( lastLabel.section(':', 0, 0) ).arg(label);
    this->currentItem()->setText(newLabel);
    }
}

void QLabelList::switchListItemSlot()
{
  if( this->currentItem() != NULL )
    {
    if( this->currentRow() == this->count() - 1 )
      {
      this->setCurrentRow(0);
      this->setCurrentItem( this->item( this->currentRow() ) );
      }
    else
      {
      this->setCurrentRow(this->currentRow() + 1);
      this->setCurrentItem( this->item( this->currentRow() ) );
      }

    cancelHighlight( this->currentItem() );
    }
}

void QLabelList::deleteListItemSlot()
{
  if( this->currentItem() != NULL )
    {
    delete this->currentItem();
    if( this->count() > 0 )
      {
      this->setCurrentRow(0);
      }
    else
      {
      this->setCurrentRow(-1);
      }
    this->setCurrentItem( this->item( this->currentRow() ) );
    cancelHighlight( this->currentItem() );
    }
}

void QLabelList::deleteListSlot()
{
  if( this->currentItem() != NULL )
    {
    this->setCurrentItem( this->item(0) );
    while( this->currentItem() != NULL )
      {
      delete this->currentItem();
      }

    this->setCurrentRow(-1);
    // reset the color seed
    m_color = 0;
    }
}

void QLabelList::deleteListItemMouseSlot(QListWidgetItem *)
{
  deleteListItemSlot();
}

void QLabelList::cancelHighlight(QListWidgetItem *)
{
  if( this->currentItem() != NULL )
    {
    // disable the highlight color
    const QColor color = this->item( this->currentRow() )->background().color();
    QPalette     pal = this->palette();
    pal.setColor(QPalette::Highlight, color);
    this->setPalette(pal);

    // signal to sliceChangeSlot indicating the update of slice viewer
    emit sliceChangedList();
    }
}

void QLabelList::sliceChangedSlot()
{
  if( this->currentItem() != NULL )
    {
    QString textLabel = this->item( this->currentRow() )->text();
    double  myPos[3];
    myPos[0] = textLabel.section(' ', 2, 2).toDouble();
    myPos[1] = textLabel.section(' ', 3, 3).toDouble();
    myPos[2] = textLabel.section(' ', 4, 4).toDouble();
    emit sendLabelPosition(myPos);
    }
}

#define MAX_LABEL_NUM 100

void QLabelList::checkVisibilitySlot()
{
  if( this->currentItem() != NULL )
    {
    int table[3 * MAX_LABEL_NUM];

    int     currRow = this->currentRow();
    QString textLabel = this->item(currRow)->text();
    double  currPos[3];
    currPos[0] = textLabel.section(' ', 2, 2).toDouble();
    currPos[1] = textLabel.section(' ', 3, 3).toDouble();
    currPos[2] = textLabel.section(' ', 4, 4).toDouble();
    double slice;

    int i;
    for( i = 0; i < 3; ++i )
      { // for different viewers
      this->setCurrentRow(0);
      while( this->item( this->currentRow() ) != NULL )
        {
        textLabel = this->item( this->currentRow() )->text();
        slice = textLabel.section(' ', 2 + i, 2 + i).toDouble();
        if( ( slice >= currPos[i] - .5 )
            && ( slice <  currPos[i] + .5 ) )
          {
          table[i * this->count() + this->currentRow()] = 1;
          }
        else
          {
          table[i * this->count() + this->currentRow()] = 0;
          }
        this->setCurrentRow(this->currentRow() + 1);
        }
      }
    this->setCurrentRow(currRow);
    emit visibilityTable(table);
    }
}

void QLabelList::checkVisibilitySlot(double *tag)
{
  if( this->currentItem() != NULL )
    {
    int     table[3 * MAX_LABEL_NUM];
    int     currRow = this->currentRow();
    QString textLabel = this->item(currRow)->text();
    double  currPos[3];
    currPos[0] = textLabel.section(' ', 2, 2).toDouble();
    currPos[1] = textLabel.section(' ', 3, 3).toDouble();
    currPos[2] = textLabel.section(' ', 4, 4).toDouble();
    double slice;

    currPos[( int(tag[1]) + 1 ) % 3] = tag[0];

    int i;
    for( i = 0; i < 3; ++i )
      { // for different viewers
      this->setCurrentRow(0);
      while( this->item( this->currentRow() ) != NULL )
        {
        textLabel = this->item( this->currentRow() )->text();
        slice = textLabel.section(' ', 2 + i, 2 + i).toDouble();
        if( ( slice >= currPos[i] - .5 )
            && ( slice <  currPos[i] + .5 ) )
          {
          table[i * this->count() + this->currentRow()] = 1;
          }
        else
          {
          table[i * this->count() + this->currentRow()] = 0;
          }
        this->setCurrentRow(this->currentRow() + 1);
        }
      }
    this->setCurrentRow(currRow);
    emit visibilityTable(table);
    }
}

void QLabelList::checkVisibilitySlot(QListWidgetItem *)
{
  checkVisibilitySlot();
}

void QLabelList::ackWheelChanged()
{
  if( this->currentItem() != NULL )
    {
    double  labelPos[3 * MAX_LABEL_NUM];
    int     currRow = this->currentRow();
    double  currPos[3];
    QString textLabel;

    this->setCurrentRow(0);
    while( this->item( this->currentRow() ) != NULL )
      {
      textLabel = this->item( this->currentRow() )->text();
      currPos[0] = textLabel.section(' ', 2, 2).toDouble();
      currPos[1] = textLabel.section(' ', 3, 3).toDouble();
      currPos[2] = textLabel.section(' ', 4, 4).toDouble();

      int i;
      for( i = 0; i < 3; ++i )
        {
        labelPos[3 * this->currentRow() + i] = currPos[i];
        }

      this->setCurrentRow(this->currentRow() + 1);
      }

    this->setCurrentRow(currRow);
    emit sendLabelPositions(labelPos);
    }
}

void QLabelList::readLandmarks()
{
  int currRow = this->currentRow();

  this->setCurrentRow(0);
  m_landmarks.clear();
  while( this->item( this->currentRow() ) != NULL )
    {
    std::vector<double> labelPos;
    QString             textLabel = this->item( this->currentRow() )->text();

    QString name      = textLabel.section(':', 0, 0);
    labelPos.push_back( textLabel.section(' ', 2, 2).toDouble() );
    labelPos.push_back( textLabel.section(' ', 3, 3).toDouble() );
    labelPos.push_back( textLabel.section(' ', 4, 4).toDouble() );
    m_landmarks.insert( std::make_pair<QString, std::vector<double> >(name, labelPos) );
    this->setCurrentRow(this->currentRow() + 1);
    }

  this->setCurrentRow(currRow);
}

void QLabelList::loadLandmarks()
{
  if( m_inputLandmarks.compare("") != 0 )
    {
    m_landmarks.clear();
    QFile input(m_inputLandmarks);
    if( !input.open(QIODevice::ReadOnly) )
      {
      std::cerr << "Cannot load landmark file!" << std::endl;
      exit(-1);
      }
    QTextStream myfile(&input);
    QString     line = myfile.readLine();
    while( !line.isNull() )
      {
      if( !line.startsWith('#') )
        {
        QString name = line.section(',', 0, 0);
        if( name.compare("") != 0 )
          {
          std::vector<double> labelPos;
          labelPos.push_back( -line.section(',', 1, 1).toDouble() );
          labelPos.push_back( -line.section(',', 2, 2).toDouble() );
          labelPos.push_back( line.section(',', 3, 3).toDouble() );
          m_landmarks.insert( std::make_pair<QString, std::vector<double> >(name, labelPos) );
          }
        }
      line = myfile.readLine();
      }

    input.close();
    }
}

void QLabelList::saveLandmarks()
{
  // if user specifies output landmarks filename
  if( m_outputLandmarks.compare("") != 0 )
    {
    writeLandmarks();
    }
  // if user specifies input landmarks filename but not the output one
  else if( m_inputLandmarks.compare("") != 0 )
    {
    m_outputLandmarks = m_inputLandmarks;
    writeLandmarks();
    }
  // if user specifies neither the input landmarks filename nor the output one
  else
    {
    saveAsLandmarks();
    }
}

void QLabelList::saveAsLandmarks()
{
  QFileDialogs fileDialog;

  m_outputLandmarks = fileDialog.saveLandmarksFile();
  if( m_outputLandmarks.compare("") != 0 )
    {
    writeLandmarks();
    }
}

void QLabelList::writeLandmarks()
{
  itk::NumberToString<double> doubleConvert;

  assert(m_outputLandmarks.compare("") != 0);
  assert(m_inputVolume.compare("") != 0);

  // get proper filename for fcsv and mrml file
  QFileInfo landmarksFileInfo(m_outputLandmarks);
  QFileInfo imageFileInfo(m_inputVolume);
  m_outputLandmarks = landmarksFileInfo.absoluteFilePath();
  QString landmarksFullFilenameWithoutExtension =
    landmarksFileInfo.absolutePath() + "/" + landmarksFileInfo.completeBaseName();
  QString imageFullFilename = imageFileInfo.absoluteFilePath();
  QString imageFilenameWithoutPath = imageFileInfo.fileName();

  // read in latest landmarks from label list class
  this->readLandmarks();

    {
    QFile output(m_outputLandmarks);
    if( !output.open(QIODevice::WriteOnly | QIODevice::Truncate) )
      {
      std::cerr << "Cannot save landmark file!" << std::endl;
      exit(-1);
      }

    QTextStream myfile(&output);
    myfile << "#Fiducial List file " << m_outputLandmarks << "\n";
    myfile << "#numPoints = " << m_landmarks.size() << "\n";
    myfile << "#symbolScale = 5\n";
    myfile << "#visibility = 1\n";
    myfile << "#textScale = 4.5\n";
    myfile << "#color = 0.4,1,1\n";
    myfile << "#selectedColor = 1,0.5,0.5\n";
    myfile << "#label,x,y,z,sel,vis\n";

    // note no LPS -> RAS is needed
    LandmarksMapType::iterator it;
    for( it = m_landmarks.begin(); it != m_landmarks.end(); ++it )
      {
      if( ( it->first ).compare("") != 0 )
        {
        myfile << it->first << ","
               << doubleConvert(-( it->second )[0]).c_str() << ","
               << doubleConvert(-( it->second )[1]).c_str() << ","
               << doubleConvert( ( it->second )[2]).c_str() << ",1,1\n";
        }
      }

    output.close();
    }

  // write scene to mrml file
  QString mrmlFullFilename = landmarksFullFilenameWithoutExtension + ".mrml";
    {
    QFile output(mrmlFullFilename);
    if( !output.open(QIODevice::WriteOnly | QIODevice::Truncate) )
      {
      std::cerr << "Cannot write mrml file!" << std::endl;
      exit(-1);
      }
    QTextStream myfile(&output);

    myfile
      <<
      "<MRML userTags=\"\">\n<Selection\n id=\"vtkMRMLSelectionNode1\" name=\"vtkMRMLSelectionNode1\" hideFromEditors=\"true\" selectable=\"true\" selected=\"false\" activeVolumeID=\"vtkMRMLScalarVolumeNode1\" secondaryVolumeID=\"NULL\" activeLabelVolumeID=\"NULL\" activeFiducialListID=\"NULL\" activeROIListID=\"NULL\" activeCameraID=\"NULL\" activeViewID=\"NULL\" activeLayoutID=\"vtkMRMLLayoutNode1\"></Selection>\n<Interaction\n id=\"vtkMRMLInteractionNode1\" name=\"vtkMRMLInteractionNode1\" hideFromEditors=\"true\" selectable=\"true\" selected=\"false\" currentInteractionMode=\"ViewTransform\" lastInteractionMode=\"ViewTransform\"></Interaction>\n<Layout\n id=\"vtkMRMLLayoutNode1\" name=\"vtkMRMLLayoutNode1\" hideFromEditors=\"true\" selectable=\"true\" selected=\"false\" currentViewArrangement=\"3\" guiPanelVisibility=\"1\" bottomPanelVisibility =\"0\" guiPanelLR=\"0\" numberOfCompareViewRows=\"0\" numberOfCompareViewColumns=\"0\" numberOfLightboxRows=\"1\" numberOfLightboxColumns=\"1\" mainPanelSize=\"400\" secondaryPanelSize=\"400\"></Layout>\n<View\n id=\"vtkMRMLViewNode1\" name=\"vtkMRMLViewNode1\" hideFromEditors=\"true\" selectable=\"true\" selected=\"false\" active=\"false\" fieldOfView=\"200\" letterSize=\"0.05\" boxVisible=\"true\" fiducialsVisible=\"true\" fiducialLabelsVisible=\"true\" axisLabelsVisible=\"true\" backgroundColor=\"0.70196 0.70196 0.90588\" animationMode=\"Off\" viewAxisMode=\"LookFrom\" spinDegrees=\"2\" spinMs=\"5\" spinDirection=\"YawLeft\" rotateDegrees=\"5\" rockLength=\"200\" rockCount=\"0\" stereoType=\"NoStereo\" renderMode=\"Perspective\"></View>\n<Camera\n id=\"vtkMRMLCameraNode1\" name=\"vtkMRMLCameraNode1\" hideFromEditors=\"true\" selectable=\"true\" selected=\"false\" position=\"0 500 0\" focalPoint=\"0 0 0\" viewUp=\"0 0 1\" parallelProjection=\"false\" parallelScale=\"1\" active=\"false\"></Camera>\n<TGParameters\n id=\"vtkMRMLChangeTrackerNode1\" name=\"vtkMRMLChangeTrackerNode1\" hideFromEditors=\"true\" selectable=\"true\" selected=\"false\" ROIMin=\"-1 -1 -1\" ROIMax=\"-1 -1 -1\" SegmentThresholdMin=\"-1\" SegmentThresholdMax=\"-1\" Analysis_Intensity_Flag=\"0\" Analysis_Deformable_Flag=\"0\" UseITK=\"1\"></TGParameters>\n<VolumeRenderingSelection\n id=\"vtkMRMLVolumeRenderingSelectionNode1\" name=\"vtkMRMLVolumeRenderingSelectionNode1\" hideFromEditors=\"true\" selectable=\"true\" selected=\"false\" activeVolumeID=\"NULL\" activeVolumeRenderingID=\"NULL\"></VolumeRenderingSelection>\n<Slice\n id=\"vtkMRMLSliceNode1\" name=\"Green\" hideFromEditors=\"true\" selectable=\"true\" selected=\"false\" fieldOfView=\"376.154 255 1\" dimensions=\"326 221 1\" activeSlice=\"0\" layoutGridRows=\"1\" layoutGridColumns=\"1\" sliceToRAS=\"-1 -0 0 -127.5 -0 -0 1 -148 0 1 0 127.5 -0 0 -0 1\" layoutName=\"Green\" orientation=\"Coronal\" jumpMode=\"1\" sliceVisibility=\"false\" widgetVisibility=\"false\" useLabelOutline=\"false\" sliceSpacingMode=\"0\" prescribedSliceSpacing=\"1 1 1\"></Slice>\n<SliceComposite\n id=\"vtkMRMLSliceCompositeNode1\" name=\"vtkMRMLSliceCompositeNode1\" hideFromEditors=\"true\" selectable=\"true\" selected=\"false\" backgroundVolumeID=\"vtkMRMLScalarVolumeNode1\" foregroundVolumeID=\"\" labelVolumeID=\"\" compositing=\"0\" labelOpacity=\"1\" linkedControl=\"1\" foregroundGrid=\"0\" backgroundGrid=\"0\" labelGrid=\"1\" fiducialVisibility=\"1\" fiducialLabelVisibility=\"1\" sliceIntersectionVisibility=\"0\" layoutName=\"Green\" annotationMode=\"All\" doPropagateVolumeSelection=\"1\"></SliceComposite>\n<Slice\n id=\"vtkMRMLSliceNode2\" name=\"Red\" hideFromEditors=\"true\" selectable=\"true\" selected=\"false\" fieldOfView=\"376.705 255 1\" dimensions=\"325 220 1\" activeSlice=\"0\" layoutGridRows=\"1\" layoutGridColumns=\"1\" sliceToRAS=\"-1 0 -0 -127.5 0 1 0 -127.5 -0 0 1 114 0 -0 0 1\" layoutName=\"Red\" orientation=\"Axial\" jumpMode=\"1\" sliceVisibility=\"false\" widgetVisibility=\"false\" useLabelOutline=\"false\" sliceSpacingMode=\"0\" prescribedSliceSpacing=\"1 1 1\"></Slice>\n<SliceComposite\n id=\"vtkMRMLSliceCompositeNode2\" name=\"vtkMRMLSliceCompositeNode2\" hideFromEditors=\"true\" selectable=\"true\" selected=\"false\" backgroundVolumeID=\"vtkMRMLScalarVolumeNode1\" foregroundVolumeID=\"\" labelVolumeID=\"\" compositing=\"0\" labelOpacity=\"1\" linkedControl=\"1\" foregroundGrid=\"0\" backgroundGrid=\"0\" labelGrid=\"1\" fiducialVisibility=\"1\" fiducialLabelVisibility=\"1\" sliceIntersectionVisibility=\"0\" layoutName=\"Red\" annotationMode=\"All\" doPropagateVolumeSelection=\"1\"></SliceComposite>\n<Slice\n id=\"vtkMRMLSliceNode3\" name=\"Yellow\" hideFromEditors=\"true\" selectable=\"true\" selected=\"false\" fieldOfView=\"375 255 1\" dimensions=\"325 221 1\" activeSlice=\"0\" layoutGridRows=\"1\" layoutGridColumns=\"1\" sliceToRAS=\"-0 0 1 -127 -1 -0 0 -127.5 -0 1 -0 127.5 0 -0 0 1\" layoutName=\"Yellow\" orientation=\"Sagittal\" jumpMode=\"1\" sliceVisibility=\"false\" widgetVisibility=\"false\" useLabelOutline=\"false\" sliceSpacingMode=\"0\" prescribedSliceSpacing=\"1 1 1\"></Slice>\n<SliceComposite\n id=\"vtkMRMLSliceCompositeNode3\" name=\"vtkMRMLSliceCompositeNode3\" hideFromEditors=\"true\" selectable=\"true\" selected=\"false\" backgroundVolumeID=\"vtkMRMLScalarVolumeNode1\" foregroundVolumeID=\"\" labelVolumeID=\"\" compositing=\"0\" labelOpacity=\"1\" linkedControl=\"1\" foregroundGrid=\"0\" backgroundGrid=\"0\" labelGrid=\"1\" fiducialVisibility=\"1\" fiducialLabelVisibility=\"1\" sliceIntersectionVisibility=\"0\" layoutName=\"Yellow\" annotationMode=\"All\" doPropagateVolumeSelection=\"1\"></SliceComposite>\n<Crosshair\n id=\"vtkMRMLCrosshairNode1\" name=\"vtkMRMLCrosshairNode1\" hideFromEditors=\"true\" selectable=\"true\" selected=\"false\" crosshairMode=\"NoCrosshair\" crosshairBehavior=\"Normal\" crosshairThickness=\"Fine\" crosshairRAS=\"-127 -148 114\"></Crosshair>\n<ScriptedModule\n id=\"vtkMRMLScriptedModuleNode1\" name=\"vtkMRMLScriptedModuleNode1\" hideFromEditors=\"true\" selectable=\"true\" selected=\"false\" ModuleName =\"Editor\" parameter0= \"label 1\"></ScriptedModule>\n";

    myfile << "<FiducialList\n id=\"vtkMRMLFiducialListNode1\" name=\"" << landmarksFullFilenameWithoutExtension
           <<
      "\" hideFromEditors=\"false\" selectable=\"true\" selected=\"false\" storageNodeRef=\"vtkMRMLFiducialListStorageNode1\" userTags=\"\" symbolScale=\"5\" symbolType=\"11\" textScale=\"4.5\" visibility=\"1\" color=\"0.4 1 1\" selectedcolor=\"1 0.5 0.5\" ambient=\"0\" diffuse=\"1\" specular=\"0\" power=\"1\" locked=\"0\" opacity=\"1\" fiducials=\"\n";

    LandmarksMapType::iterator it;
    unsigned int               index = 0;
    for( it = m_landmarks.begin(); it != m_landmarks.end(); ++it )
      {
      myfile << "id " << it->first << " labeltext " << it->first << " xyz "
             << doubleConvert( ( it->second )[0]).c_str() << " "
             << doubleConvert( ( it->second )[1]).c_str() << " "
             << doubleConvert( ( it->second )[2]).c_str();
                         if( ++index < m_landmarks.size() )
                           {
                           myfile << " orientationwxyz 0 0 0 1 selected 1 visibility 1\n";
                           }
                         else
                           {
                           myfile << " orientationwxyz 0 0 0 1 selected 1 visibility 1\"></FiducialList>\n";
                           }
                         }

                         myfile << "<FiducialListStorage\n id=\"vtkMRMLFiducialListStorageNode1\" "
                         << "name=\"vtkMRMLFiducialListStorageNode1\" hideFromEditors=\"true\" "
                         << "selectable=\"true\" selected=\"false\" fileName=\""
                         << m_outputLandmarks
                         << "\" useCompression=\"1\" readState=\"0\" writeState=\"0\"></FiducialListStorage>\n"
                         << "<VolumeArchetypeStorage\n id=\"vtkMRMLVolumeArchetypeStorageNode1\" "
                         << "name=\"vtkMRMLVolumeArchetypeStorageNode1\" hideFromEditors=\"true\""
                         << " selectable=\"true\" selected=\"false\" fileName=\""
                         << imageFullFilename
                         << "\" useCompression=\"1\" readState=\"0\" writeState=\"0\" "
                         <<
                         "centerImage=\"0\" singleFile=\"0\" UseOrientationFromFile=\"1\"></VolumeArchetypeStorage>\n"
                         << "<Volume\n id=\"vtkMRMLScalarVolumeNode1\" name=\"" << imageFilenameWithoutPath
                         << "\" hideFromEditors=\"false\" selectable=\"true\" selected=\"false\" "
                         << "storageNodeRef=\"vtkMRMLVolumeArchetypeStorageNode1\" userTags=\"\" "
                         << "displayNodeRef=\"vtkMRMLScalarVolumeDisplayNode1\" ijkToRASDirections=\"-1 "
                         << "  -0   0 -0   -0   1 0 1 0 \" spacing=\"1 1 1\" origin=\"-0 -255 -0\" "
                         << "labelMap=\"0\"></Volume>\n<VolumeDisplay\n id=\"vtkMRMLScalarVolumeDisplayNode1\" "
                         << "name=\"vtkMRMLScalarVolumeDisplayNode1\" hideFromEditors=\"true\" "
                         << "selectable=\"true\" selected=\"false\" color=\"0.5 0.5 0.5\" selectedColor=\"1 "
                         << "0 0\" selectedAmbient=\"0.4\" ambient=\"0\" diffuse=\"1\" selectedSpecular=\"0.5\" "
                         << "specular=\"0\" power=\"1\" opacity=\"1\" visibility=\"true\" clipping=\"false\" "
                         << "sliceIntersectionVisibility=\"false\" backfaceCulling=\"true\" scalarVisibility=\"false\" "
                         << "vectorVisibility=\"false\" tensorVisibility=\"false\" autoScalarRange=\"true\" "
                         << "scalarRange=\"0 100\" colorNodeRef=\"vtkMRMLColorTableNodeGrey\"  "
                         << "window=\"204\" level=\"153\" upperThreshold=\"32767\" lowerThreshold=\"-32768\" "
                         <<
                         "interpolate=\"1\" autoWindowLevel=\"1\" applyThreshold=\"0\" autoThreshold=\"0\"></VolumeDisplay>\n</MRML>\n";

                         output.close();
                         }
                         }
