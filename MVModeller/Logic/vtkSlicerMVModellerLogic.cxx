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

// MVModeller Logic includes
#include "vtkSlicerMVModellerLogic.h"

// MRML includes
#include <vtkMRMLScene.h>
#include <vtkMRMLNode.h>
#include <vtkMRMLMarkupsFiducialNode.h>
#include <vtkMRMLInteractionNode.h>
#include <vtkMRMLSelectionNode.h>
#include <vtkMRMLModelNode.h>
#include <vtkMRMLModelDisplayNode.h>
#include <vtkMRMLTransformNode.h>
#include <vtkMRMLSliceNode.h>

// VTK includes
#include <vtkIntArray.h>
#include <vtkNew.h>
#include <vtkObjectFactory.h>
#include <vtkPlane.h>
#include <vtkPlaneSource.h>
#include <vtkCardinalSpline.h>
#include <vtkPoints.h>
#include <vtkCellArray.h>
#include <vtkPolyData.h>
#include <vtkTubeFilter.h>
#include <vtkPolyDataAlgorithm.h>
#include <vtkPlaneSource.h>
#include <vtkMath.h>
#include <vtkTransform.h>
#include <vtkIdList.h>
#include <vtkMatrix4x4.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkCollection.h>

// STD includes
#include <cassert>
#include <cmath>


//----------------------------------------------------------------------------
vtkStandardNewMacro(vtkSlicerMVModellerLogic);

//----------------------------------------------------------------------------
vtkSlicerMVModellerLogic::vtkSlicerMVModellerLogic()
{
}

//----------------------------------------------------------------------------
vtkSlicerMVModellerLogic::~vtkSlicerMVModellerLogic()
{
}

//----------------------------------------------------------------------------
void vtkSlicerMVModellerLogic::PrintSelf(ostream& os, vtkIndent indent)
{
    this->Superclass::PrintSelf(os, indent);
}

//---------------------------------------------------------------------------
void vtkSlicerMVModellerLogic::SetMRMLSceneInternal(vtkMRMLScene * newScene)
{
    vtkNew<vtkIntArray> events;
    events->InsertNextValue(vtkMRMLScene::NodeAddedEvent);
    events->InsertNextValue(vtkMRMLScene::NodeRemovedEvent);
    events->InsertNextValue(vtkMRMLScene::EndBatchProcessEvent);
    this->SetAndObserveMRMLSceneEventsInternal(newScene, events.GetPointer());
}

//-----------------------------------------------------------------------------
void vtkSlicerMVModellerLogic::RegisterNodes()
{
    assert(this->GetMRMLScene() != 0);
}

//---------------------------------------------------------------------------
void vtkSlicerMVModellerLogic::UpdateFromMRMLScene()
{
    assert(this->GetMRMLScene() != 0);
}

//---------------------------------------------------------------------------
void vtkSlicerMVModellerLogic
::OnMRMLSceneNodeAdded(vtkMRMLNode* vtkNotUsed(node))
{
}

//---------------------------------------------------------------------------
void vtkSlicerMVModellerLogic
::OnMRMLSceneNodeRemoved(vtkMRMLNode* vtkNotUsed(node))
{
}

//---------------------------------------------------------------------------
void vtkSlicerMVModellerLogic::drawMVOpening(vtkMRMLNode * node)
{
    vtkWarningMacro("Entering draw opening vtk logic\n");
    vtkMRMLApplicationLogic *mrmlAppLogic = this->GetMRMLApplicationLogic();
    vtkMRMLInteractionNode *inode = mrmlAppLogic->GetInteractionNode();
    vtkMRMLSelectionNode *snode = mrmlAppLogic->GetSelectionNode();

    snode->SetReferenceActivePlaceNodeClassName("vtkMRMLMarkupsFiducialNode");

    inode->SwitchToPersistentPlaceMode();
    inode->SetCurrentInteractionMode(vtkMRMLInteractionNode::Place);
}

//---------------------------------------------------------------------------
void vtkSlicerMVModellerLogic::closeMVOpening(vtkMRMLNode * node)
{
    vtkWarningMacro("Entering close opening vtk logic\n");
    vtkMRMLApplicationLogic *mrmlAppLogic = this->GetMRMLApplicationLogic();
    vtkMRMLInteractionNode *inode = mrmlAppLogic->GetInteractionNode();
    vtkMRMLSelectionNode *snode = mrmlAppLogic->GetSelectionNode();

    inode->SwitchToViewTransformMode();
    inode->SwitchToSinglePlaceMode();
    inode->SetCurrentInteractionMode(vtkMRMLInteractionNode::ViewTransform);

    if( !node )
    {
        vtkWarningMacro("Bad node pointer" << endl);
    }

    nodeToPolyCardinalSpline( vtkMRMLMarkupsFiducialNode::SafeDownCast(node) );
}

//---------------------------------------------------------------------------
void vtkSlicerMVModellerLogic::nodeToPolyCardinalSpline(vtkMRMLMarkupsFiducialNode* sourceNode)
{
    if (!sourceNode)
    {
        vtkWarningMacro("Bad node pointer" << endl);
        return;
    }

    //Get number of control points and initialize vectors
    int nCtrlPoints = sourceNode->GetNumberOfFiducials();
    vtkWarningMacro("nCtrlPoints" << nCtrlPoints << endl);
    double pos[3] = { 0, 0, 0 };
    double centroid[3] = { 0, 0, 0 };

    vtkNew<vtkCardinalSpline> aSplineX;
    vtkNew<vtkCardinalSpline> aSplineY;
    vtkNew<vtkCardinalSpline> aSplineZ;

    aSplineX->ClosedOn();
    aSplineY->ClosedOn();
    aSplineZ->ClosedOn();

    //Add control points and calculate centroid
    for(int i = 0; i < nCtrlPoints; i++)
    {
        sourceNode->GetNthFiducialPosition(i, pos);

        aSplineX->AddPoint(i, pos[0]);
        aSplineY->AddPoint(i, pos[1]);
        aSplineZ->AddPoint(i, pos[2]);

        centroid[0] += pos[0];
        centroid[1] += pos[1];
        centroid[2] += pos[2];
    }

    centroid[0] = centroid[0] / nCtrlPoints;
    centroid[1] = centroid[1] / nCtrlPoints;
    centroid[2] = centroid[2] / nCtrlPoints;

    //    Interpolate x, y and z by using the three spline filters and create new points
    int nInterpolatedPoints = 52 * (nCtrlPoints - 1);
    vtkWarningMacro("nInterpolatedPoints" << nInterpolatedPoints << endl);
    vtkNew<vtkPoints> points;

    double r[2] = {0,0};
    aSplineX->GetParametricRange(r);

    double t = r[0];
    int p = 0;
    double tStep = (static_cast<double>(nCtrlPoints) - 1) / (static_cast<double>(nInterpolatedPoints) - 1);
    int nOutputPoints = 0;

    while(t < r[1] + 1.0)
    {
        points->InsertPoint(p, aSplineX->Evaluate(t), aSplineY->Evaluate(t), aSplineZ->Evaluate(t));
        t = t + tStep;
        p++;
    }

    nOutputPoints = p;


    //Create model defining mv opening
    vtkNew<vtkCellArray> lines;
    lines->InsertNextCell(nOutputPoints);
    for(int i = 0; i < nOutputPoints; i++)
    {
        lines->InsertCellPoint(i);
    }

    vtkNew<vtkPolyData> outputPoly;
    outputPoly->SetPoints(points.GetPointer());
    outputPoly->SetLines(lines.GetPointer());

    vtkNew<vtkTubeFilter> tubeFilter;
    tubeFilter->SetInputData(outputPoly.GetPointer());
    tubeFilter->SetRadius(0.5);
    tubeFilter->SetNumberOfSides(20);
    tubeFilter->CappingOn();
    tubeFilter->Update();

    vtkNew<vtkMRMLModelNode> model;
    this->GetMRMLScene()->AddNode(model.GetPointer());
    model->SetScene(this->GetMRMLScene());
    vtkDebugMacro("Added model node");

    model->CreateDefaultDisplayNodes();

    model->SetAndObservePolyData(tubeFilter->GetOutput());
    model->Modified();

    vtkMRMLDisplayNode* displayNode = model->GetDisplayNode();
    if (displayNode)
    {
        displayNode->SetActiveScalarName("Curvature");
        displayNode->SliceIntersectionVisibilityOn();
    }
    else
    {
        vtkWarningMacro("Couldn't get display node");
    }

    vtkMRMLModelNode* lastPlane;
    double lastNorm[3], lastCenter[3];

    ///Temp logic for creating a plane from polydata points
    int planeStep = nOutputPoints / 11;
    for (int i = 0; i < nOutputPoints - 1; i += planeStep)
    {
        double pos1[3], pos2[3];
        points->GetPoint(i, pos1);
        points->GetPoint(i + 1, pos2);

        double norm[3] = { pos2[0] - pos1[0], pos2[1] - pos1[1], pos2[2] - pos1[2] };


        vtkNew<vtkPlaneSource> planeSource;
        planeSource->SetOrigin(centroid);

        double v1[3], v2[3];

        v1[0] = pos1[0] - centroid[0];
        v1[1] = pos1[1] - centroid[1];
        v1[2] = pos1[2] - centroid[2];

        vtkMath::Cross(norm, v1, v2);
        vtkMath::Normalize(v2);

        //Set plane scale where x is centroid  -> interesction point and y is vertical through the valve
        double yScale = 50;
        double xScale = 3.5;

        planeSource->SetPoint1(centroid[0] + xScale * v1[0], centroid[1] + xScale * v1[1], centroid[2] + xScale * v1[2]);
        planeSource->SetPoint2(centroid[0] + yScale * v2[0], centroid[1] + yScale * v2[1], centroid[2] + yScale * v2[2]);

        planeSource->SetNormal(norm);
        planeSource->SetCenter(pos1);

        planeSource->Update();

        vtkNew<vtkMRMLModelNode> planeModel;
        planeModel->SetAndObservePolyData(planeSource->GetOutput());
        planeModel->SetName("PlaneModel");
        planeModel->Modified();


        planeModel->SetScene(this->GetMRMLScene());

        this->GetMRMLScene()->AddNode(planeModel.GetPointer());
        vtkDebugMacro("Added model node");

        if (!planeModel->GetDisplayNode())
        {
            planeModel->CreateDefaultDisplayNodes();
        }

        vtkMRMLDisplayNode* displayNodePlane = planeModel->GetDisplayNode();
        if (displayNode)
        {
            displayNodePlane->SetActiveScalarName("Plane");
            displayNodePlane->BackfaceCullingOff();
            displayNodePlane->FrontfaceCullingOff();
            displayNodePlane->SliceIntersectionVisibilityOn();
            displayNodePlane->SetOpacity(0.6);
        }

        lastPlane = planeModel.GetPointer();
        lastNorm[0] = norm[0];
        lastNorm[1] = norm[1];
        lastNorm[2] = norm[2];
        planeSource->GetCenter(lastCenter);
    }

    lastPlane->GetDisplayNode()->SetColor(0,0,200);



    vtkSmartPointer<vtkCollection> planes = vtkSmartPointer<vtkCollection>::Take(this->GetMRMLScene()->GetNodesByClassByName("vtkMRMLModelNode", "PlaneModel"));
    vtkSmartPointer<vtkMRMLSliceNode> yellowSlice = vtkMRMLSliceNode::SafeDownCast(this->GetMRMLScene()->GetNodeByID("vtkMRMLSliceNodeYellow"));
    if(yellowSlice)
    {
        vtkSmartPointer<vtkMRMLModelNode> model = vtkMRMLModelNode::SafeDownCast(planes->GetItemAsObject(0));
        vtkSmartPointer<vtkPoints> points = model->GetPolyData()->GetPoints();

        vtkNew<vtkPlaneSource> plane;
        plane->SetOrigin(points->GetPoint(1));
        plane->SetPoint1(points->GetPoint(0));
        plane->SetPoint2(points->GetPoint(2));

        double N[3], T[3], P[3];
        //N = normal
        plane->GetNormal(N);
        //T = [0,1,0]
        T[0] = 0;
        T[1] = 1;
        T[2] = 0;
        //P = plane center
        plane->GetCenter(P);
        yellowSlice->SetSliceToRASByNTP(N[0], N[1], N[2], T[0] , T[1] , T[2], P[0], P[1], P[2], 0);
    }
    ///End temp logic
    ///Move into full calcs of all planes when complete
}

//---------------------------------------------------------------------------
void vtkSlicerMVModellerLogic::generateOpeningPlanes()
{

}
