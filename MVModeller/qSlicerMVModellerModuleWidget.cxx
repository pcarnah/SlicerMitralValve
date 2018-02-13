/*==============================================================================

  Program: 3D Slicer

  Portions (c) Copyright Brigham and Women's Hospital (BWH) All Rights Reserved.

  See COPYRIGHT.txt
  or http://www.slicer.org/copyright/copyright.txt for details.

  Unless required by applicable law or agreed to in writing, software
  distributed under the License is distributed on an "AS IS" BASIS,
  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  See the License for the specific language governing permissions and
  limitations under the License.

==============================================================================*/

// Qt includes
#include <QDebug>

// SlicerQt includes
#include "qSlicerMVModellerModuleWidget.h"
#include "ui_qSlicerMVModellerModuleWidget.h"

#include "vtkMRMLMarkupsFiducialNode.h"
#include "vtkMRMLScene.h"

#include "vtkSlicerMVModellerLogic.h"

//-----------------------------------------------------------------------------
/// \ingroup Slicer_QtModules_ExtensionTemplate
class qSlicerMVModellerModuleWidgetPrivate: public Ui_qSlicerMVModellerModuleWidget
{
public:
  qSlicerMVModellerModuleWidgetPrivate();
};

//-----------------------------------------------------------------------------
// qSlicerMVModellerModuleWidgetPrivate methods

//-----------------------------------------------------------------------------
qSlicerMVModellerModuleWidgetPrivate::qSlicerMVModellerModuleWidgetPrivate()
{
}

//-----------------------------------------------------------------------------
// qSlicerMVModellerModuleWidget methods

//-----------------------------------------------------------------------------
qSlicerMVModellerModuleWidget::qSlicerMVModellerModuleWidget(QWidget* _parent)
  : Superclass( _parent )
  , d_ptr( new qSlicerMVModellerModuleWidgetPrivate )
{
}

//-----------------------------------------------------------------------------
qSlicerMVModellerModuleWidget::~qSlicerMVModellerModuleWidget()
{
}

//-----------------------------------------------------------------------------
void qSlicerMVModellerModuleWidget::setup()
{
  Q_D(qSlicerMVModellerModuleWidget);
  d->setupUi(this);
  this->Superclass::setup();

  d->MVOpening->setMRMLScene(this->mrmlScene());
  connect(this, SIGNAL(mrmlSceneChanged(vtkMRMLScene*)), d->MVOpening, SLOT(setMRMLScene(vtkMRMLScene*)));

  d->PlaneStepWidget->setMRMLScene(this->mrmlScene());
  connect(this, SIGNAL(mrmlSceneChanged(vtkMRMLScene*)), d->PlaneStepWidget, SLOT(setMRMLScene(vtkMRMLScene*)));

  // Connections
  connect(d->MVOpening, SIGNAL(closeMVOpening(vtkMRMLNode*)), this, SLOT(closeMVOpening(vtkMRMLNode*)));
  connect(d->MVOpening, SIGNAL(drawMVOpening(vtkMRMLNode*)), this, SLOT(drawMVOpening(vtkMRMLNode*)));
  connect(d->PlaneStepWidget, SIGNAL(selectMVPlane(int)), this, SLOT(selectMVPlane(int)));
  connect(d->PlaneStepWidget, SIGNAL(beginDrawPlane()), this, SLOT(drawPlaneProfile()));
  connect(d->PlaneStepWidget, SIGNAL(endDrawPlane(int)), this, SLOT(endPlaneProfile(int)));
}

//-----------------------------------------------------------------------------
void qSlicerMVModellerModuleWidget::drawMVOpening(vtkMRMLNode * node)
{
    vtkSlicerMVModellerLogic::SafeDownCast(this->logic())->drawMVOpening(node);
}

//-----------------------------------------------------------------------------
void qSlicerMVModellerModuleWidget::closeMVOpening(vtkMRMLNode * node)
{
    vtkSlicerMVModellerLogic::SafeDownCast(this->logic())->closeMVOpening(node);
}

//-----------------------------------------------------------------------------
void qSlicerMVModellerModuleWidget::selectMVPlane(const int & i)
{
    vtkSlicerMVModellerLogic::SafeDownCast(this->logic())->selectMVPlane(i);
}

//-----------------------------------------------------------------------------
void qSlicerMVModellerModuleWidget::drawPlaneProfile()
{
    vtkSlicerMVModellerLogic::SafeDownCast(this->logic())->beginDrawPlaneSpline();
}

//-----------------------------------------------------------------------------
void qSlicerMVModellerModuleWidget::endPlaneProfile(const int &i)
{
    vtkSlicerMVModellerLogic::SafeDownCast(this->logic())->closePlaneSpline(i);
}
