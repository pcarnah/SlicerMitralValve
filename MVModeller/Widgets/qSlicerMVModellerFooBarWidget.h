/*==============================================================================

  Program: 3D Slicer

  Copyright (c) Kitware Inc.

  See COPYRIGHT.txt
  or http://www.slicer.org/copyright/copyright.txt for details.

  Unless required by applicable law or agreed to in writing, software
  distributed under the License is distributed on an "AS IS" BASIS,
  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  See the License for the specific language governing permissions and
  limitations under the License.

  This file was originally developed by Jean-Christophe Fillion-Robin, Kitware Inc.
  and was partially funded by NIH grant 3P41RR013218-12S1

==============================================================================*/

#ifndef __qSlicerMVModellerFooBarWidget_h
#define __qSlicerMVModellerFooBarWidget_h

// Qt includes
#include <QWidget>

// FooBar Widgets includes
#include "qSlicerMVModellerModuleWidgetsExport.h"

class qSlicerMVModellerFooBarWidgetPrivate;

/// \ingroup Slicer_QtModules_MVModeller
class Q_SLICER_MODULE_MVMODELLER_WIDGETS_EXPORT qSlicerMVModellerFooBarWidget
  : public QWidget
{
  Q_OBJECT
public:
  typedef QWidget Superclass;
  qSlicerMVModellerFooBarWidget(QWidget *parent=0);
  virtual ~qSlicerMVModellerFooBarWidget();

protected slots:

protected:
  QScopedPointer<qSlicerMVModellerFooBarWidgetPrivate> d_ptr;

private:
  Q_DECLARE_PRIVATE(qSlicerMVModellerFooBarWidget);
  Q_DISABLE_COPY(qSlicerMVModellerFooBarWidget);
};

#endif
