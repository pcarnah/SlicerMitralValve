#ifndef QSLICERMVMODELLERMVOPENINGWIDGET_H
#define QSLICERMVMODELLERMVOPENINGWIDGET_H

#include <QWidget>

// MVOpening Widgets includes
#include "qSlicerMVModellerModuleWidgetsExport.h"

// MRML Includes
#include "qMRMLWidget.h"

// CTK includes
#include <ctkVTKObject.h>

class qSlicerMVModellerMVOpeningWidgetPrivate;
class vtkMRMLNode;

/// \ingroup Slicer_QtModules_MVModeller
class Q_SLICER_MODULE_MVMODELLER_WIDGETS_EXPORT qSlicerMVModellerMVOpeningWidget
        : public qMRMLWidget
{
    Q_OBJECT
    QVTK_OBJECT
    public:
        typedef qMRMLWidget Superclass;
    qSlicerMVModellerMVOpeningWidget(QWidget *parent = 0);
    virtual ~qSlicerMVModellerMVOpeningWidget();

public slots:
    virtual void setMRMLScene(vtkMRMLScene* newScene) override;
    void on_buttonBeginOpening_clicked();
    void on_buttonConfirmOpening_clicked();

signals:
    void drawMVOpening(vtkMRMLNode*);
    void closeMVOpening(vtkMRMLNode*);

protected slots:


protected:
    QScopedPointer<qSlicerMVModellerMVOpeningWidgetPrivate> d_ptr;

private:
    Q_DECLARE_PRIVATE(qSlicerMVModellerMVOpeningWidget);
    Q_DISABLE_COPY(qSlicerMVModellerMVOpeningWidget);
};


#endif // QSLICERMVMODELLERMVOPENINGWIDGET_H
