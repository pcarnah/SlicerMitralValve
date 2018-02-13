#include "qSlicerMVModellerPlaneStepWidget.h"
#include "ui_qSlicerMVModellerPlaneStepWidget.h"

//-----------------------------------------------------------------------------
/// \ingroup Slicer_QtModules_MVModeller
class qSlicerMVModellerPlaneStepWidgetPrivate
        : public Ui_qSlicerMVModellerPlaneStepWidget
{
    Q_DECLARE_PUBLIC(qSlicerMVModellerPlaneStepWidget);
protected:
    qSlicerMVModellerPlaneStepWidget* const q_ptr;

public:
    qSlicerMVModellerPlaneStepWidgetPrivate(
            qSlicerMVModellerPlaneStepWidget& object);
    virtual void setupUi(qSlicerMVModellerPlaneStepWidget*);
};

// --------------------------------------------------------------------------
qSlicerMVModellerPlaneStepWidgetPrivate
::qSlicerMVModellerPlaneStepWidgetPrivate(
        qSlicerMVModellerPlaneStepWidget& object)
    : q_ptr(&object)
{
}

// --------------------------------------------------------------------------
void qSlicerMVModellerPlaneStepWidgetPrivate
::setupUi(qSlicerMVModellerPlaneStepWidget* widget)
{
    this->Ui_qSlicerMVModellerPlaneStepWidget::setupUi(widget);
}

//-----------------------------------------------------------------------------
// qSlicerMVModellerPlaneStepWidget methods

//-----------------------------------------------------------------------------
qSlicerMVModellerPlaneStepWidget
::qSlicerMVModellerPlaneStepWidget(QWidget* parentWidget)
    : Superclass(parentWidget)
    , d_ptr(new qSlicerMVModellerPlaneStepWidgetPrivate(*this))
{
    Q_D(qSlicerMVModellerPlaneStepWidget);
    d->setupUi(this);
}

//-----------------------------------------------------------------------------
qSlicerMVModellerPlaneStepWidget
::~qSlicerMVModellerPlaneStepWidget()
{
}

//-----------------------------------------------------------------------------
void qSlicerMVModellerPlaneStepWidget::setMRMLScene(vtkMRMLScene* newScene)
{
    Q_D(qSlicerMVModellerPlaneStepWidget);
    if(newScene == this->mrmlScene())
        return;

    Superclass::setMRMLScene(newScene);
}

// --------------------------------------------------------------------------
void qSlicerMVModellerPlaneStepWidget::on_sliderSelectPlane_valueChanged(double i)
{
    int index = static_cast<int>(i) - 1;

    emit selectMVPlane(index);
}

// --------------------------------------------------------------------------
void qSlicerMVModellerPlaneStepWidget::on_beginDrawPlaneButton_clicked()
{
    Q_D(qSlicerMVModellerPlaneStepWidget);

    //d->sliderSelectPlane->setEnabled(false);
    emit beginDrawPlane();
}

// --------------------------------------------------------------------------
void qSlicerMVModellerPlaneStepWidget::on_endDrawPlaneButton_clicked()
{
    Q_D(qSlicerMVModellerPlaneStepWidget);

    //d->sliderSelectPlane->setEnabled(true);

    emit endDrawPlane(static_cast<int>(d->sliderSelectPlane->value()));
}
