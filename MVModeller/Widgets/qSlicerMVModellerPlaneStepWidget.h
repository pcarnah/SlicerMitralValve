#ifndef QSLICERMVMODELLERPLANESTEPWIDGET_H
#define QSLICERMVMODELLERPLANESTEPWIDGET_H

#include <QWidget>

// PlaneStep Widgets includes
#include "qSlicerMVModellerModuleWidgetsExport.h"

// MRML Includes
#include "qMRMLWidget.h"

// CTK includes
#include <ctkVTKObject.h>

class qSlicerMVModellerPlaneStepWidgetPrivate;

/// \ingroup Slicer_QtModules_MVModeller
class Q_SLICER_MODULE_MVMODELLER_WIDGETS_EXPORT qSlicerMVModellerPlaneStepWidget
        : public qMRMLWidget
{
    Q_OBJECT
    QVTK_OBJECT
    public:
        typedef qMRMLWidget Superclass;
    qSlicerMVModellerPlaneStepWidget(QWidget *parent = 0);
    virtual ~qSlicerMVModellerPlaneStepWidget();

public slots:
    virtual void setMRMLScene(vtkMRMLScene* newScene) override;

    void on_sliderSelectPlane_valueChanged(double);
    void on_beginDrawPlaneButton_clicked();
    void on_endDrawPlaneButton_clicked();

signals:
    void selectMVPlane(const int& i);
    void beginDrawPlane();
    void endDrawPlane(const int& i);

protected slots:


protected:
    QScopedPointer<qSlicerMVModellerPlaneStepWidgetPrivate> d_ptr;

private:
    Q_DECLARE_PRIVATE(qSlicerMVModellerPlaneStepWidget);
    Q_DISABLE_COPY(qSlicerMVModellerPlaneStepWidget);
};

#endif // QSLICERMVMODELLERPLANESTEPWIDGET_H
