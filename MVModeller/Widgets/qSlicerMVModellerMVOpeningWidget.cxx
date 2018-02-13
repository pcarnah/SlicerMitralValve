#include "qSlicerMVModellerMVOpeningWidget.h"
#include "ui_qSlicerMVModellerMVOpeningWidget.h"

#include <QDebug>

#include "vtkMRML.h"
#include "vtkMRMLNode.h"
#include "vtkMRMLMarkupsFiducialNode.h"
#include "vtkMRMLScene.h"
#include "vtkObjectBase.h"
#include "vtkMRMLSliceNode.h"

//-----------------------------------------------------------------------------
/// \ingroup Slicer_QtModules_MVModeller
class qSlicerMVModellerMVOpeningWidgetPrivate
        : public Ui_qSlicerMVModellerMVOpeningWidget
{
    Q_DECLARE_PUBLIC(qSlicerMVModellerMVOpeningWidget);
protected:
    qSlicerMVModellerMVOpeningWidget* const q_ptr;

public:
    qSlicerMVModellerMVOpeningWidgetPrivate(
            qSlicerMVModellerMVOpeningWidget& object);
    virtual void setupUi(qSlicerMVModellerMVOpeningWidget*);
};

// --------------------------------------------------------------------------
qSlicerMVModellerMVOpeningWidgetPrivate
::qSlicerMVModellerMVOpeningWidgetPrivate(
        qSlicerMVModellerMVOpeningWidget& object)
    : q_ptr(&object)
{
}

// --------------------------------------------------------------------------
void qSlicerMVModellerMVOpeningWidgetPrivate
::setupUi(qSlicerMVModellerMVOpeningWidget* widget)
{
    this->Ui_qSlicerMVModellerMVOpeningWidget::setupUi(widget);
}

//-----------------------------------------------------------------------------
// qSlicerMVModellerMVOpeningWidget methods

//-----------------------------------------------------------------------------
qSlicerMVModellerMVOpeningWidget
::qSlicerMVModellerMVOpeningWidget(QWidget* parentWidget)
    : Superclass(parentWidget)
    , d_ptr(new qSlicerMVModellerMVOpeningWidgetPrivate(*this))
{
    Q_D(qSlicerMVModellerMVOpeningWidget);
    d->setupUi(this);

    //hide currently unused controls (to be implemented)
    d->buttonSetOpening->hide();
    d->sliderOpeningPlane->hide();

    d->comboSourceSelector->setMRMLScene(this->mrmlScene());
    connect(this, SIGNAL(mrmlSceneChanged(vtkMRMLScene*)), d->comboSourceSelector, SLOT(setMRMLScene(vtkMRMLScene*)));
}

//-----------------------------------------------------------------------------
qSlicerMVModellerMVOpeningWidget
::~qSlicerMVModellerMVOpeningWidget()
{
}

//-----------------------------------------------------------------------------
void qSlicerMVModellerMVOpeningWidget::setMRMLScene(vtkMRMLScene* newScene)
{
    Q_D(qSlicerMVModellerMVOpeningWidget);
    if(newScene == this->mrmlScene())
        return;

    Superclass::setMRMLScene(newScene);
}

// --------------------------------------------------------------------------
void qSlicerMVModellerMVOpeningWidget::on_buttonBeginOpening_clicked()
{
    Q_D(qSlicerMVModellerMVOpeningWidget);

    emit drawMVOpening(d->comboSourceSelector->currentNode());
}

// --------------------------------------------------------------------------
void qSlicerMVModellerMVOpeningWidget::on_buttonConfirmOpening_clicked()
{
    Q_D(qSlicerMVModellerMVOpeningWidget);

	qWarning() << "emit confirm opening";

    emit closeMVOpening(d->comboSourceSelector->currentNode());
}

