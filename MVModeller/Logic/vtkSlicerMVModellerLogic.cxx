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
#include <vtkMRMLSliceLogic.h>
#include <vtkMRMLScalarVolumeNode.h>

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
#include <vtkImageReslice.h>
#include <vtkImageData.h>
#include <vtkImageMapper.h>

// STD includes
#include <cassert>
#include <cmath>


//----------------------------------------------------------------------------
vtkStandardNewMacro(vtkSlicerMVModellerLogic);

//----------------------------------------------------------------------------
vtkSlicerMVModellerLogic::vtkSlicerMVModellerLogic()
{
    profile = 0;
    planes = vtkSmartPointer<vtkCollection>::New();
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


    vtkSmartPointer<vtkMRMLMarkupsFiducialNode> fidNode = vtkMRMLMarkupsFiducialNode::SafeDownCast(node);
    vtkSmartPointer<vtkPolyData> poly = nodeToPolyCardinalSpline( fidNode, true );

    vtkNew<vtkMRMLModelNode> model;
    model->SetName("MVProfileSpline");

    this->GetMRMLScene()->AddNode(model.GetPointer());
    model->SetScene(this->GetMRMLScene());
    vtkDebugMacro("Added model node");

    model->CreateDefaultDisplayNodes();

    model->SetAndObservePolyData(poly);
    model->Modified();

    profile = poly;


    vtkMRMLDisplayNode* displayNode = model->GetDisplayNode();
    if (displayNode)
    {
        displayNode->SetActiveScalarName("Curvature");
        displayNode->SliceIntersectionVisibilityOn();
        displayNode->SetLineWidth(2);
    }
    else
    {
        vtkWarningMacro("Couldn't get display node");
    }

    generateOpeningPlanes();

    fidNode->GetDisplayNode()->VisibilityOff();
}

//---------------------------------------------------------------------------
vtkSmartPointer<vtkPolyData> vtkSlicerMVModellerLogic::nodeToPolyCardinalSpline(vtkMRMLMarkupsFiducialNode* sourceNode, bool closed)
{
    vtkSmartPointer<vtkPolyData> outputPoly = vtkSmartPointer<vtkPolyData>::New();
    if (!sourceNode)
    {
        vtkWarningMacro("Bad node pointer" << endl);
        return outputPoly;
    }

    //Get number of control points and initialize vectors
    int nCtrlPoints = sourceNode->GetNumberOfFiducials();
    double pos[3] = { 0, 0, 0 };

    vtkNew<vtkCardinalSpline> aSplineX;
    vtkNew<vtkCardinalSpline> aSplineY;
    vtkNew<vtkCardinalSpline> aSplineZ;

    if (closed)
    {
        aSplineX->ClosedOn();
        aSplineY->ClosedOn();
        aSplineZ->ClosedOn();
    }
    else
    {
        aSplineX->ClosedOff();
        aSplineY->ClosedOff();
        aSplineZ->ClosedOff();
    }


    //Add control points and calculate centroid
    for (int i = 0; i < nCtrlPoints; i++)
    {
        sourceNode->GetNthFiducialPosition(i, pos);

        aSplineX->AddPoint(i, pos[0]);
        aSplineY->AddPoint(i, pos[1]);
        aSplineZ->AddPoint(i, pos[2]);
    }

    //    Interpolate x, y and z by using the three spline filters and create new points
    int nInterpolatedPoints = 52 * (nCtrlPoints - 1);
    vtkNew<vtkPoints> points;

    double r[2] = {0,0};
    aSplineX->GetParametricRange(r);

    double t = r[0];
    int p = 0;
    double tStep = (static_cast<double>(nCtrlPoints) - 1) / (static_cast<double>(nInterpolatedPoints) - 1);
    int nOutputPoints = 0;

    if (closed)
    {
        while (t < r[1] + 1.0)
        {
            points->InsertPoint(p, aSplineX->Evaluate(t), aSplineY->Evaluate(t), aSplineZ->Evaluate(t));
            t = t + tStep;
            p++;
        }
        points->InsertPoint(p, aSplineX->Evaluate(r[0]), aSplineY->Evaluate(r[0]), aSplineZ->Evaluate(r[0]));
        p = p + 1;
        points->InsertPoint(p, aSplineX->Evaluate(r[0] + tStep), aSplineY->Evaluate(r[0] + tStep), aSplineZ->Evaluate(r[0] + tStep));
    }
    else
    {
        while (t < r[1])
        {
            points->InsertPoint(p, aSplineX->Evaluate(t), aSplineY->Evaluate(t), aSplineZ->Evaluate(t));
            t = t + tStep;
            p++;
        }
    }

    nOutputPoints = p;


    //Create model defining mv opening
    vtkNew<vtkCellArray> lines;
    lines->InsertNextCell(nOutputPoints);
    for(int i = 0; i < nOutputPoints; i++)
    {
        lines->InsertCellPoint(i);
    }

    outputPoly->SetPoints(points.GetPointer());
    outputPoly->SetLines(lines.GetPointer());

    return outputPoly;
}

//---------------------------------------------------------------------------
void vtkSlicerMVModellerLogic::generateOpeningPlanes()
{
    if (!profile)
    {
        vtkWarningMacro("Couldn't find model defining MV opening profile.");
        return;
    }

    if( !planes )
    {
        planes = vtkSmartPointer<vtkCollection>::New();
    }

    planes->RemoveAllItems();

    vtkSmartPointer<vtkPoints> points = profile->GetPoints();

    int nPoints = points->GetNumberOfPoints();
    int nCtrlPoints = nPoints / 10;

    double centroid[3] = { 0, 0, 0 };
    double pos[3];

    for (int i = 0; i < nCtrlPoints; i++)
    {
        points->GetPoint(i * 10, pos);

        centroid[0] += pos[0];
        centroid[1] += pos[1];
        centroid[2] += pos[2];
    }

    centroid[0] = centroid[0] / nCtrlPoints;
    centroid[1] = centroid[1] / nCtrlPoints;
    centroid[2] = centroid[2] / nCtrlPoints;

    int planeStep = nPoints / 11;
    for (int i = 0; i < 11; i++)
    {
        double pos1[3], pos2[3];
        points->GetPoint(i * planeStep, pos1);
        points->GetPoint((i * planeStep) + 1, pos2);

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
        if (displayNodePlane)
        {
            displayNodePlane->SetActiveScalarName("Plane");
            displayNodePlane->BackfaceCullingOff();
            displayNodePlane->FrontfaceCullingOff();
            displayNodePlane->SliceIntersectionVisibilityOn();
            displayNodePlane->SetOpacity(0.6);

            if( i != 0 && i != 5 && i && 10 )
            {
                displayNodePlane->VisibilityOff();
            }
        }

        planes->AddItem(planeSource.GetPointer());
    }
}

//---------------------------------------------------------------------------
void vtkSlicerMVModellerLogic::selectMVPlane(const int & i)
{
    if( !planes )
    {
        return;
    }

    vtkSmartPointer<vtkMRMLSliceNode> yellowSlice = vtkMRMLSliceNode::SafeDownCast(this->GetMRMLScene()->GetNodeByID("vtkMRMLSliceNodeYellow"));
    if(yellowSlice && i < planes->GetNumberOfItems())
    {
        vtkSmartPointer<vtkPlaneSource> plane = vtkPlaneSource::SafeDownCast( planes->GetItemAsObject(i) );

        double N[3], T[3], P[3];

        //N = normal
        plane->GetNormal(N);

        //T = [0,1,0]
        T[0] = 0;
        T[1] = N[0] >= 0 ? 1 : -1;      //Adjust direction so slice is oreiented up
        T[2] = 0;

        //P = plane center
        plane->GetCenter(P);

        vtkDebugMacro("Plane#: " << i << " N:" << N[0] << "," << N[1] << "," << N[2] << " P:" << P[0] << "," << P[1] << "," << P[2]);
        yellowSlice->SetSliceToRASByNTP(N[0], N[1], N[2], T[0] , T[1] , T[2], P[0], P[1], P[2], 0);
    }
}
