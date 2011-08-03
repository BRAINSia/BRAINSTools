/*
 * Author: Wei Lu
 * at Psychiatry Imaging Lab,
 * University of Iowa Health Care 2010
 */

#include "QHelpDialog.h"

QHelpDialog::QHelpDialog(QWidget *myParent) : QDialog(myParent)
{
  QString text(
    "This help message provides the basic usage of the GUI tool. Please visit our wiki page for more information:<br /><b>http://www.nitrc.org/plugins/mwiki/index.php/brainscdetector:MainPage</b><br /><br /><b>Select</b><br />To select a landmark:<br />* Click on the landmark in the label list<br /><br /><b>Move</b><br />To change the location of a landmark:<br />* Select the landmark<br />* Move either one of the slider bar to find a good slice<br />* Click at a point of the corresponding viewer to move the landmark to the current place<br />* As the other viewers change slices accordingly, adjust the landmark location in 3D<br /><br /><b>Zoom in/out</b><br />To zoom in/out of a viewer:<br />* Hover the mouse on that viewer<br />* Right click then drag up/down to zoom in/out about the center<br /><br /><b>Add</b><br />To add a new landmark:<br />* Click the Add button<br /><br /><b>Remove</b><br />To remove a landmark:<br />* Double click the landmark in the label list, then choose Accept in the dialog window, or<br />* Select the landmark, then click on the Remove button in the GUI<br /><br />Other operations such as <b>Remove All</b> landmarks, <b>Save</b> landmarks, etc can be achieved by clicking a proper button on the GUI.");

  this->setWindowTitle("GUI Help Info");
  m_label.setText(text);
  m_label.setTextInteractionFlags(Qt::TextBrowserInteraction);
  m_layout = new QGridLayout(this);
  m_layout->addWidget(&m_label, 1, 1, Qt::AlignLeft);
}
