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
#include <vtkNew.h>
#include <vtkPlane.h>
#include <vtkPlaneSource.h>
#include <vtkCardinalSpline.h>
#include <vtkPoints.h>
#include <vtkCellArray.h>
#include <vtkPolyData.h>
#include <vtkMath.h>
#include <vtkTransform.h>
#include <vtkIdList.h>
#include <vtkMatrix4x4.h>
#include <vtkCollection.h>
#include <vtkAppendPolyData.h>
#include <vtkDelaunay3D.h>
#include <vtkDelaunay2D.h>
#include <vtkUnstructuredGrid.h>
#include <vtkLineSource.h>
#include <vtkSplineFilter.h>
#include <vtkIdTypeArray.h>
#include <vtkPolyDataNormals.h>
#include <vtkDataSetSurfaceFilter.h>
#include <vtkCleanPolyData.h>
#include <vtkTriangleFilter.h>
#include <vtkWedge.h>
#include <vtkPolyLine.h>
#include <vtkKochanekSpline.h>
#include <vtkFillHolesFilter.h>
#include <vtkDecimatePro.h>
#include <vtkQuad.h>
#include <vtkPolygon.h>
#include <vtkCellIterator.h>

// STD includes
#include <cassert>
#include <cmath>
#include <array>
#include <limits>

using namespace std;


//----------------------------------------------------------------------------
vtkStandardNewMacro(vtkSlicerMVModellerLogic);

//----------------------------------------------------------------------------
vtkSlicerMVModellerLogic::vtkSlicerMVModellerLogic()
{
    profile = 0;
    fidNode = 0;
    planes = vtkSmartPointer<vtkCollection>::New();
    leafletSplines = vtkSmartPointer<vtkCollection>::New();
    for(int i = 0; i < 11; i++)
    {
        leafletSplines->AddItem(vtkSmartPointer<vtkPolyData>::New());
    }
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
vtkSmartPointer<vtkPolyData> vtkSlicerMVModellerLogic::nodeToPolyCardinalSpline(vtkMRMLMarkupsFiducialNode* sourceNode, bool closed, int nSubs)
{
    vtkSmartPointer<vtkPolyData> outputPoly = vtkSmartPointer<vtkPolyData>::New();
    if (!sourceNode)
    {
        vtkWarningMacro("Bad node pointer" << endl);
        return outputPoly;
    }

    if( closed )
    {
        //Get number of control points and initialize vectors
        int nCtrlPoints = sourceNode->GetNumberOfFiducials();
        double pos[3] = { 0, 0, 0 };

        vtkNew<vtkCardinalSpline> aSplineX;
        vtkNew<vtkCardinalSpline> aSplineY;
        vtkNew<vtkCardinalSpline> aSplineZ;

        aSplineX->ClosedOn();
        aSplineY->ClosedOn();
        aSplineZ->ClosedOn();

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

        double r[2] = { 0, 0 };
        aSplineX->GetParametricRange(r);

        double t = r[0];
        int p = 0;
        double tStep = (static_cast<double>(nCtrlPoints)-1) / (static_cast<double>(nInterpolatedPoints)-1);
        int nOutputPoints = 0;


        while (t < r[1] + 1.0)
        {
            points->InsertPoint(p, aSplineX->Evaluate(t), aSplineY->Evaluate(t), aSplineZ->Evaluate(t));
            t = t + tStep;
            p++;
        }
        points->InsertPoint(p, aSplineX->Evaluate(r[0]), aSplineY->Evaluate(r[0]), aSplineZ->Evaluate(r[0]));
        p = p + 1;
        points->InsertPoint(p, aSplineX->Evaluate(r[0] + tStep), aSplineY->Evaluate(r[0] + tStep), aSplineZ->Evaluate(r[0] + tStep));



        nOutputPoints = p;


        //Create model defining mv opening
        vtkNew<vtkCellArray> lines;
        lines->InsertNextCell(nOutputPoints);
        for (int i = 0; i < nOutputPoints; i++)
        {
            lines->InsertCellPoint(i);
        }

        outputPoly->SetPoints(points.GetPointer());
        outputPoly->SetLines(lines.GetPointer());
    }
    else
    {
        vtkSmartPointer<vtkPolyData> inputPoly = vtkSmartPointer<vtkPolyData>::New();

        //Get number of control points and initialize vectors
        int nCtrlPoints = sourceNode->GetNumberOfFiducials();
        double pos[3] = { 0, 0, 0 };

        vtkNew<vtkPoints> points;
        vtkNew<vtkCellArray> lines;
        lines->InsertNextCell(nCtrlPoints);
        for (int i = 0; i < nCtrlPoints; i++)
        {
            sourceNode->GetNthFiducialPosition(i, pos);

            points->InsertPoint(i, pos);
            lines->InsertCellPoint(i);
        }

        inputPoly->SetPoints(points.GetPointer());
        inputPoly->SetLines(lines.GetPointer());

        vtkNew<vtkSplineFilter> spline;
        spline->SetInputData(inputPoly);
        spline->SetSubdivideToSpecified();
        spline->SetNumberOfSubdivisions(nSubs);
        spline->Update();

        outputPoly->DeepCopy(spline->GetOutput());
    }

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
            displayNodePlane->VisibilityOff();
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

//---------------------------------------------------------------------------
void vtkSlicerMVModellerLogic::beginDrawPlaneSpline()
{
    vtkMRMLApplicationLogic *mrmlAppLogic = this->GetMRMLApplicationLogic();
    vtkMRMLInteractionNode *inode = mrmlAppLogic->GetInteractionNode();
    vtkMRMLSelectionNode *snode = mrmlAppLogic->GetSelectionNode();

    snode->SetReferenceActivePlaceNodeClassName("vtkMRMLMarkupsFiducialNode");

    if( !fidNode )
    {
        fidNode = vtkSmartPointer<vtkMRMLMarkupsFiducialNode>::New();
        fidNode->SetName("P");
        fidNode->SetScene(this->GetMRMLScene());

        this->GetMRMLScene()->AddNode(fidNode);
    }
    else
    {
        fidNode->RemoveAllMarkups();

        if( fidNode->GetScene() == NULL )
        {
            fidNode->SetScene(this->GetMRMLScene());

            this->GetMRMLScene()->AddNode(fidNode);
        }
    }

    snode->SetReferenceActivePlaceNodeID(fidNode->GetID());

    inode->SwitchToPersistentPlaceMode();
    inode->SetCurrentInteractionMode(vtkMRMLInteractionNode::Place);
}

//---------------------------------------------------------------------------
void vtkSlicerMVModellerLogic::closePlaneSpline(int planeNum)
{
    if( !fidNode )
    {
        return;
    }

    vtkMRMLApplicationLogic *mrmlAppLogic = this->GetMRMLApplicationLogic();
    vtkMRMLInteractionNode *inode = mrmlAppLogic->GetInteractionNode();
    vtkMRMLSelectionNode *snode = mrmlAppLogic->GetSelectionNode();

    inode->SwitchToViewTransformMode();
    inode->SwitchToSinglePlaceMode();
    inode->SetCurrentInteractionMode(vtkMRMLInteractionNode::ViewTransform);

    vtkSmartPointer<vtkPolyData> poly = nodeToPolyCardinalSpline( fidNode, false, LEAFLET_SPLINE_SUBDIVISIONS );

    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    points->DeepCopy(poly->GetPoints());
    vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();

    double p1[3], p2[3], p3[3];
    points->GetPoint(0, p1);
    points->GetPoint(points->GetNumberOfPoints() - 1, p2);

    double x,y,z;
    if( p1[2] < p2[2] )
    {
        swap(p1, p2);
    }

    p3[0] = p1[0];
    p3[1] = p1[1];
    p3[2] = p2[2];

    vtkNew<vtkLineSource> line1, line2;
    line1->SetPoint1(p1);
    line1->SetPoint2(p3);
    line1->SetResolution(25);

    line2->SetPoint1(p3);
    line2->SetPoint2(p2);
    line2->SetResolution(25);

    vtkNew<vtkAppendPolyData> append;
    append->AddInputData(poly);
    append->AddInputConnection(line1->GetOutputPort());
    append->AddInputConnection(line2->GetOutputPort());
    append->Update();


    vtkNew<vtkMRMLModelNode> model;
    model->SetName("MVLeafletSpline");

    this->GetMRMLScene()->AddNode(model.GetPointer());
    model->SetScene(this->GetMRMLScene());
    vtkDebugMacro("Added model node");

    model->CreateDefaultDisplayNodes();

    model->SetAndObservePolyData(append->GetOutput());
    model->Modified();

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

    if( planeNum <= leafletSplines->GetNumberOfItems() )
    {
        leafletSplines->ReplaceItem(planeNum - 1, model->GetPolyData());
    }

    fidNode->RemoveAllMarkups();
}

//---------------------------------------------------------------------------
void vtkSlicerMVModellerLogic::generateSurface()
{
    vtkNew<vtkAppendPolyData> append;
    int pointsPerSpline = 0, splineCount = 0;

    for(int i = 0; i < leafletSplines->GetNumberOfItems(); i++)
    {
        vtkSmartPointer<vtkPolyData> poly = vtkPolyData::SafeDownCast(leafletSplines->GetItemAsObject(i));
        if(poly->GetNumberOfPoints() != 0)
        {
            append->AddInputData(poly);
            pointsPerSpline = poly->GetNumberOfPoints();
            splineCount++;
        }
    }
    vtkSmartPointer<vtkPolyData> poly = vtkPolyData::SafeDownCast(leafletSplines->GetItemAsObject(0));
    if(poly->GetNumberOfPoints() != 0)
    {
        append->AddInputData(poly);
    }
    append->Update();

    vtkSmartPointer<vtkPolyData> inputPoly = vtkSmartPointer<vtkPolyData>::New();
    inputPoly->DeepCopy(append->GetOutput());

    vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();

    for(int i = 0; i < pointsPerSpline; i+=2)
    {
        lines->InsertNextCell(splineCount + 1);
        for(int j = 0; j < splineCount; j++)
        {
            lines->InsertCellPoint(i + (j * pointsPerSpline));
        }
        lines->InsertCellPoint(i);
    }


    inputPoly->SetLines(lines);

    vtkNew<vtkSplineFilter> spline;
    spline->SetInputData(inputPoly);
    spline->SetSubdivideToSpecified();
    spline->SetNumberOfSubdivisions(70);
    spline->GetSpline()->ClosedOn();
    spline->Update();

    vtkSmartPointer<vtkPolyData> linesPoly = spline->GetOutput();

    vtkSmartPointer<vtkPolyData> polys = vtkSmartPointer<vtkPolyData>::New();
    polys->Allocate();

    vtkSmartPointer<vtkPolygon> face1 = vtkSmartPointer<vtkPolygon>::New();
    vtkSmartPointer<vtkPolygon> face2 = vtkSmartPointer<vtkPolygon>::New();

    vtkSmartPointer<vtkCellIterator> lineIter1 = linesPoly->NewCellIterator();
    lineIter1->InitTraversal();
    vtkSmartPointer<vtkCellIterator> lineIter2 = linesPoly->NewCellIterator();
    lineIter2->InitTraversal();
    lineIter2->GoToNextCell();

    for(int i = 0; i < linesPoly->GetNumberOfCells() - 1; i++)
    {
        vtkNew<vtkIdList> outerLine1;
        outerLine1->DeepCopy(linesPoly->GetCell(i)->GetPointIds());
        vtkNew<vtkIdList> outerLine2;
        outerLine2->DeepCopy(linesPoly->GetCell(i + 1)->GetPointIds());

        for(int j = 0; j < outerLine1->GetNumberOfIds() - 1; j++)
        {
            vtkSmartPointer<vtkQuad> wedge = vtkSmartPointer<vtkQuad>::New();
            wedge->GetPointIds()->SetId(0, outerLine1->GetId(j));
            wedge->GetPointIds()->SetId(1, outerLine2->GetId(j));
            wedge->GetPointIds()->SetId(2, outerLine2->GetId(j + 1));
            wedge->GetPointIds()->SetId(3, outerLine1->GetId(j + 1));
            polys->InsertNextCell(wedge->GetCellType(), wedge->GetPointIds());
        }

        face1->GetPointIds()->InsertId(i, outerLine1->GetId(0));
        face2->GetPointIds()->InsertId(i, outerLine1->GetId(outerLine1->GetNumberOfIds() - 1));

        if( i == linesPoly->GetNumberOfCells() - 2)
        {
            face1->GetPointIds()->InsertId(i+1, outerLine2->GetId(0));
            face2->GetPointIds()->InsertId(i+1, outerLine2->GetId(outerLine1->GetNumberOfIds() - 1));
        }

    }

    polys->InsertNextCell(face1->GetCellType(), face1->GetPointIds());
    int cell = polys->InsertNextCell(face2->GetCellType(), face2->GetPointIds());
	polys->ReverseCell(cell);

    polys->SetPoints(linesPoly->GetPoints());

    vtkNew<vtkTriangleFilter> tri;
    tri->SetInputData(polys);
    tri->Update();

    vtkNew<vtkFillHolesFilter> fill;
    fill->SetInputConnection(tri->GetOutputPort());
    fill->Update();

    vtkNew<vtkPolyDataNormals> norm;
    norm->SetInputConnection(tri->GetOutputPort());
    norm->AutoOrientNormalsOn();
    norm->SetFeatureAngle(60);
    norm->Update();

    vtkNew<vtkMRMLModelNode> model;
    model->SetName("MVSurface");

    this->GetMRMLScene()->AddNode(model.GetPointer());
    model->SetScene(this->GetMRMLScene());

    model->CreateDefaultDisplayNodes();

    model->SetAndObservePolyData(norm->GetOutput());
    model->Modified();

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
}

vtkCollection *vtkSlicerMVModellerLogic::getLeafletSplines() const
{
    return leafletSplines;
}

//---------------------------------------------------------------------------
vtkPolyData *vtkSlicerMVModellerLogic::getProfile() const
{
    return profile;
}

//---------------------------------------------------------------------------
vtkPolyData *vtkSlicerMVModellerLogic::getMergedLeafletSplines()
{
    vtkNew<vtkAppendPolyData> append;
    int pointsPerSpline = 0, splineCount = 0;

    for(int i = 0; i < leafletSplines->GetNumberOfItems(); i++)
    {
        vtkSmartPointer<vtkPolyData> poly = vtkPolyData::SafeDownCast(leafletSplines->GetItemAsObject(i));
        if(poly && poly->GetNumberOfPoints() > 0)
        {
            append->AddInputData(poly);
            pointsPerSpline = poly->GetNumberOfPoints();
            splineCount++;
        }
    }
    append->Update();

    vtkSmartPointer<vtkPolyData> inputPoly = vtkSmartPointer<vtkPolyData>::New();
    inputPoly->DeepCopy(append->GetOutput());

    vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();

    for(int i = 0; i < pointsPerSpline; i++)
    {
        lines->InsertNextCell(splineCount);
        for(int j = 0; j < splineCount; j++)
        {
            lines->InsertCellPoint(i + (j * pointsPerSpline));
        }

    }

    inputPoly->SetLines(lines);

    vtkNew<vtkSplineFilter> spline;
    spline->SetInputData(inputPoly);
    spline->SetSubdivideToSpecified();
    spline->SetNumberOfSubdivisions(70);
    spline->GetSpline()->ClosedOff();
    spline->Update();

    append->AddInputData(spline->GetOutput());
    append->Update();

    vtkPolyData *poly = vtkPolyData::New();
    poly->DeepCopy(append->GetOutput());
    return poly;
}
