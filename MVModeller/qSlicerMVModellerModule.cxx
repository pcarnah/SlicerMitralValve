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
#include <QtPlugin>

// MVModeller Logic includes
#include <vtkSlicerMVModellerLogic.h>

// MVModeller includes
#include "qSlicerMVModellerModule.h"
#include "qSlicerMVModellerModuleWidget.h"

//-----------------------------------------------------------------------------
Q_EXPORT_PLUGIN2(qSlicerMVModellerModule, qSlicerMVModellerModule);

//-----------------------------------------------------------------------------
/// \ingroup Slicer_QtModules_ExtensionTemplate
class qSlicerMVModellerModulePrivate
{
public:
  qSlicerMVModellerModulePrivate();
};

//-----------------------------------------------------------------------------
// qSlicerMVModellerModulePrivate methods

//-----------------------------------------------------------------------------
qSlicerMVModellerModulePrivate::qSlicerMVModellerModulePrivate()
{
}

//-----------------------------------------------------------------------------
// qSlicerMVModellerModule methods

//-----------------------------------------------------------------------------
qSlicerMVModellerModule::qSlicerMVModellerModule(QObject* _parent)
  : Superclass(_parent)
  , d_ptr(new qSlicerMVModellerModulePrivate)
{
}

//-----------------------------------------------------------------------------
qSlicerMVModellerModule::~qSlicerMVModellerModule()
{
}

//-----------------------------------------------------------------------------
QString qSlicerMVModellerModule::helpText() const
{
  return "This is a loadable module that can be bundled in an extension";
}

//-----------------------------------------------------------------------------
QString qSlicerMVModellerModule::acknowledgementText() const
{
  return "This work was partially funded by NIH grant NXNNXXNNNNNN-NNXN";
}

//-----------------------------------------------------------------------------
QStringList qSlicerMVModellerModule::contributors() const
{
  QStringList moduleContributors;
  moduleContributors << QString("John Doe (AnyWare Corp.)");
  return moduleContributors;
}

//-----------------------------------------------------------------------------
QIcon qSlicerMVModellerModule::icon() const
{
  return QIcon(":/Icons/MVModeller.png");
}

//-----------------------------------------------------------------------------
QStringList qSlicerMVModellerModule::categories() const
{
  return QStringList() << "Examples";
}

//-----------------------------------------------------------------------------
QStringList qSlicerMVModellerModule::dependencies() const
{
  return QStringList();
}

//-----------------------------------------------------------------------------
void qSlicerMVModellerModule::setup()
{
  this->Superclass::setup();
}

//-----------------------------------------------------------------------------
qSlicerAbstractModuleRepresentation* qSlicerMVModellerModule
::createWidgetRepresentation()
{
  return new qSlicerMVModellerModuleWidget;
}

//-----------------------------------------------------------------------------
vtkMRMLAbstractLogic* qSlicerMVModellerModule::createLogic()
{
  return vtkSlicerMVModellerLogic::New();
}
