import os
import unittest
import vtk, qt, ctk, slicer
import SimpleITK as sitk
import sitkUtils
from slicer.ScriptedLoadableModule import *
import HeartValveLib
import numpy as np
import math
import logging


#
# MVSegmenter
#

class MVSegmenter(ScriptedLoadableModule):
    """Uses ScriptedLoadableModule base class, available at:
    https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
    """

    def __init__(self, parent):
        ScriptedLoadableModule.__init__(self, parent)
        self.parent.title = "Mitral Valve Segmenter"
        self.parent.categories = ["Cardiac"]
        self.parent.dependencies = []
        self.parent.contributors = [
            "Patrick Carnahan (Robarts Research Institute)"]
        self.parent.helpText = """
This module implements an algorithm for automatic mitral valve segmentation using ITK.
"""
        self.parent.helpText += self.getDefaultModuleDocumentationLink()
        self.parent.acknowledgementText = """
This file was originally developed by Jean-Christophe Fillion-Robin, Kitware Inc.
and Steve Pieper, Isomics, Inc. and was partially funded by NIH grant 3P41RR013218-12S1.
"""  # replace with organization, grant and thanks.


#
# MVSegmenterWidget
#

class MVSegmenterWidget(ScriptedLoadableModuleWidget):
    """Uses ScriptedLoadableModuleWidget base class, available at:
    https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
    """

    def __init__(self, parent):
        ScriptedLoadableModuleWidget.__init__(self, parent)

        self.logic = MVSegmenterLogic()

    def setup(self):
        ScriptedLoadableModuleWidget.setup(self)

        # Vertical spacing between sections
        vSpace = 10

        # Instantiate and connect widgets ...

        #
        # Parameters Area
        #
        parametersCollapsibleButton = ctk.ctkCollapsibleButton()
        parametersCollapsibleButton.text = "Parameters"
        self.layout.addWidget(parametersCollapsibleButton)

        # Layout within the dummy collapsible button
        parametersFormLayout = qt.QFormLayout(parametersCollapsibleButton)

        #
        #   Heart valve node selector
        #
        self.heartValveSelector = slicer.qMRMLNodeComboBox()
        self.heartValveSelector.nodeTypes = ["vtkMRMLScriptedModuleNode"]
        self.heartValveSelector.setNodeTypeLabel("HeartValve", "vtkMRMLScriptedModuleNode")
        self.heartValveSelector.baseName = "HeartValve"
        self.heartValveSelector.addAttribute("vtkMRMLScriptedModuleNode", "ModuleName", "HeartValve")
        self.heartValveSelector.addEnabled = False
        self.heartValveSelector.removeEnabled = True
        self.heartValveSelector.noneEnabled = False
        self.heartValveSelector.showHidden = True  # scripted module nodes are hidden by default
        self.heartValveSelector.renameEnabled = True
        self.heartValveSelector.setMRMLScene(slicer.mrmlScene)
        self.heartValveSelector.setToolTip("Select heart valve node where annulus was defined")
        parametersFormLayout.addRow("Heart valve: ", self.heartValveSelector)

        #
        # input volume selector
        #
        self.inputSelector = slicer.qMRMLNodeComboBox()
        self.inputSelector.nodeTypes = ["vtkMRMLScalarVolumeNode"]
        self.inputSelector.selectNodeUponCreation = True
        self.inputSelector.addEnabled = False
        self.inputSelector.removeEnabled = False
        self.inputSelector.noneEnabled = False
        self.inputSelector.showHidden = False
        self.inputSelector.showChildNodeTypes = False
        self.inputSelector.setMRMLScene(slicer.mrmlScene)
        self.inputSelector.setToolTip("Pick the input to the algorithm.")
        parametersFormLayout.addRow("Input Volume", self.inputSelector)

        #
        # output segmentation selector
        #
        self.outputSegmentationSelector = slicer.qMRMLNodeComboBox()
        self.outputSegmentationSelector.nodeTypes = ["vtkMRMLSegmentationNode"]
        self.outputSegmentationSelector.selectNodeUponCreation = True
        self.outputSegmentationSelector.addEnabled = True
        self.outputSegmentationSelector.removeEnabled = False
        self.outputSegmentationSelector.noneEnabled = True
        self.outputSegmentationSelector.showHidden = False
        self.outputSegmentationSelector.showChildNodeTypes = False
        self.outputSegmentationSelector.setMRMLScene(slicer.mrmlScene)
        self.outputSegmentationSelector.setToolTip("Pick the output to the algorithm.")
        parametersFormLayout.addRow("Output Segmentation", self.outputSegmentationSelector)

        #
        # check box to trigger taking screen shots for later use in tutorials
        #
        self.enableScreenshotsFlagCheckBox = qt.QCheckBox()
        self.enableScreenshotsFlagCheckBox.checked = 0
        self.enableScreenshotsFlagCheckBox.setToolTip(
            "If checked, take screen shots for tutorials. Use Save Data to write them to disk.")
        parametersFormLayout.addRow("Enable Screenshots", self.enableScreenshotsFlagCheckBox)

        # Add vertical spacer
        self.layout.addSpacing(vSpace)

        #
        #  First Phase Segmentation
        #
        firstPassCollapsibleButton = ctk.ctkCollapsibleButton()
        firstPassCollapsibleButton.text = "Blood Pool Segmentation"
        self.layout.addWidget(firstPassCollapsibleButton)

        # Layout within the dummy collapsible button
        firstPassFormLayout = qt.QFormLayout(firstPassCollapsibleButton)

        #
        # Initialize BP Button
        #
        self.initBPButton = qt.QPushButton("Initialize Segmentation")
        self.initBPButton.toolTip = "Run the initial segmentation pass."
        self.initBPButton.enabled = False
        firstPassFormLayout.addRow(self.initBPButton)

        #
        #  Increment First Pass Buttons
        #
        incrementFirstHBox = qt.QHBoxLayout()
        incrementFirstHBox.addStretch(5)

        self.incrementFirstButton50 = qt.QPushButton("+50")
        self.incrementFirstButton50.toolTip = "Run the algorithm for 50 more iterations"
        self.incrementFirstButton50.enabled = False
        incrementFirstHBox.addWidget(self.incrementFirstButton50)

        self.incrementFirstButton100 = qt.QPushButton("+100")
        self.incrementFirstButton100.toolTip = "Run the algorithm for 100 more iterations"
        self.incrementFirstButton100.enabled = False
        incrementFirstHBox.addWidget(self.incrementFirstButton100)

        self.incrementFirstButton500 = qt.QPushButton("+500")
        self.incrementFirstButton500.toolTip = "Run the algorithm for 500 more iterations"
        self.incrementFirstButton500.enabled = False
        incrementFirstHBox.addWidget(self.incrementFirstButton500)

        self.undoButtonBP = qt.QPushButton("Undo")
        self.undoButtonBP.toolTip = "Undo previous step"
        self.undoButtonBP.enabled = False
        incrementFirstHBox.addWidget(self.undoButtonBP)

        self.redoButtonBP = qt.QPushButton("Redo")
        self.redoButtonBP.toolTip = "Redo previous step"
        self.redoButtonBP.enabled = False
        incrementFirstHBox.addWidget(self.redoButtonBP)

        firstPassFormLayout.addRow("Increment Segmentation", incrementFirstHBox)

        # Add vertical spacer
        self.layout.addSpacing(vSpace)

        #
        #  Second Phase Segmentation
        #
        secondPassCollapsibleButton = ctk.ctkCollapsibleButton()
        secondPassCollapsibleButton.text = "Leaflet Segmentation"
        self.layout.addWidget(secondPassCollapsibleButton)

        # Layout within the dummy collapsible button
        secondPassFormLayout = qt.QFormLayout(secondPassCollapsibleButton)

        # Initialize
        self.initLeafletButton = qt.QPushButton("Initialize Segmentation")
        self.initLeafletButton.toolTip = "Run the initial leaflet segmentation pass."
        self.initLeafletButton.enabled = False
        secondPassFormLayout.addRow(self.initLeafletButton)

        #
        # Increment Buttons
        #

        incrementHBox = qt.QHBoxLayout()
        incrementHBox.addStretch(5)

        self.incrementButton10 = qt.QPushButton("+10")
        self.incrementButton10.toolTip = "Run the algorithm for 10 more iterations."
        self.incrementButton10.enabled = False
        incrementHBox.addWidget(self.incrementButton10)

        self.incrementButton50 = qt.QPushButton("+50")
        self.incrementButton50.toolTip = "Run the algorithm for 50 more iterations."
        self.incrementButton50.enabled = False
        incrementHBox.addWidget(self.incrementButton50)

        self.incrementButton200 = qt.QPushButton("+200")
        self.incrementButton200.toolTip = "Run the algorithm for 200 more iterations."
        self.incrementButton200.enabled = False
        incrementHBox.addWidget(self.incrementButton200)

        self.undoButtonLeaflet = qt.QPushButton("Undo")
        self.undoButtonLeaflet.toolTip = "Undo previous step"
        self.undoButtonLeaflet.enabled = False
        incrementHBox.addWidget(self.undoButtonLeaflet)

        self.redoButtonLeaflet = qt.QPushButton("Redo")
        self.redoButtonLeaflet.toolTip = "Redo previous step"
        self.redoButtonLeaflet.enabled = False
        incrementHBox.addWidget(self.redoButtonLeaflet)

        secondPassFormLayout.addRow("Increment Segmentation", incrementHBox)

        # Add vertical spacer
        self.layout.addSpacing(vSpace)

        #
        #  Manual Adjustment
        #
        manualAdjCollapsibleButton = ctk.ctkCollapsibleButton()
        manualAdjCollapsibleButton.text = "Manual Adjustment"
        manualAdjCollapsibleButton.collapsed = True
        # self.layout.addWidget(manualAdjCollapsibleButton)

        manualAdjFormLayout = qt.QFormLayout(manualAdjCollapsibleButton)

        # Markups node selector
        self.markupsSelector = slicer.qMRMLNodeComboBox()
        self.markupsSelector.nodeTypes = ["vtkMRMLMarkupsFiducialNode"]
        self.markupsSelector.selectNodeUponCreation = True
        self.markupsSelector.addEnabled = True
        self.markupsSelector.removeEnabled = False
        self.markupsSelector.renameEnabled = True
        self.markupsSelector.noneEnabled = False
        self.markupsSelector.showHidden = False
        self.markupsSelector.showChildNodeTypes = False
        self.markupsSelector.setMRMLScene(slicer.mrmlScene)
        self.markupsSelector.setToolTip("Pick the markups node used for adjusting model.")
        manualAdjFormLayout.addRow("Markups Node", self.markupsSelector)

        # Activate button
        self.generateSurfaceMarkups = qt.QPushButton("Generate Surface Points")
        self.generateSurfaceMarkups.toolTip = "Generates points on the surface of the model that will allow the model to be deformed"
        self.generateSurfaceMarkups.enabled = False
        manualAdjFormLayout.addRow(self.generateSurfaceMarkups)

        # Add vertical spacer
        self.layout.addSpacing(vSpace)

        #
        #  Export Inner Surface Model
        #
        generateSurfaceMoldCollapsibleButton = ctk.ctkCollapsibleButton()
        generateSurfaceMoldCollapsibleButton.text = "Generate Mold"
        generateSurfaceMoldCollapsibleButton.collapsed = True
        self.layout.addWidget(generateSurfaceMoldCollapsibleButton)

        exportModelFormLayout = qt.QFormLayout(generateSurfaceMoldCollapsibleButton)

        self.baseDepthSlider = ctk.ctkSliderWidget()
        self.baseDepthSlider.singleStep = 0.1
        self.baseDepthSlider.pageStep = 1
        self.baseDepthSlider.minimum = -20
        self.baseDepthSlider.maximum = 0
        self.baseDepthSlider.value = -10
        self.baseDepthSlider.decimals = 1
        self.baseDepthSlider.tracking = False
        self.baseDepthSlider.setToolTip("Adjust the depth of the base surface for the mold")
        exportModelFormLayout.addRow("Base Clipping Depth", self.baseDepthSlider)

        # Export button
        self.generateMoldButton = qt.QPushButton("Generate Mold")
        self.generateMoldButton.toolTip = "Generate the inner surface mold in the segmentation node."
        self.generateMoldButton.enabled = False
        exportModelFormLayout.addRow(self.generateMoldButton)

        self.exportMoldButton = qt.QPushButton("Export Mold to Model")
        self.exportMoldButton.toolTip = "Export mold segmentation node to model."
        self.exportMoldButton.enabled = False
        exportModelFormLayout.addRow(self.exportMoldButton)

        # Add vertical spacer
        self.layout.addSpacing(vSpace)

        # connections
        self.initBPButton.connect('clicked(bool)', self.onInitBPButton)
        self.incrementFirstButton50.connect('clicked(bool)', self.onIncrementFirst50Button)
        self.incrementFirstButton100.connect('clicked(bool)', self.onIncrementFirst100Button)
        self.incrementFirstButton500.connect('clicked(bool)', self.onIncrementFirst500Button)
        self.undoButtonBP.connect('clicked(bool)', self.onUndoButtonBP)
        self.redoButtonBP.connect('clicked(bool)', self.onRedoButtonBP)

        self.initLeafletButton.connect('clicked(bool)', self.onInitLeafletButton)
        self.incrementButton10.connect('clicked(bool)', self.onIncrement10Button)
        self.incrementButton50.connect('clicked(bool)', self.onIncrement50Button)
        self.incrementButton200.connect('clicked(bool)', self.onIncrement200Button)
        self.undoButtonLeaflet.connect('clicked(bool)', self.onUndoButtonLeaflet)
        self.redoButtonLeaflet.connect('clicked(bool)', self.onRedoButtonLeaflet)

        self.heartValveSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.onSelect)
        self.inputSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.onSelect)
        self.outputSegmentationSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.onSelect)

        self.markupsSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.onSelect)
        self.generateSurfaceMarkups.connect('clicked(bool)', self.onGenerateSurfaceMarkups)

        self.generateMoldButton.connect('clicked(bool)', self.onGenerateModelButton)
        self.exportMoldButton.connect('clicked(bool)', self.onExportModelButton)

        # Add vertical spacer
        self.layout.addStretch(1)

        # Refresh Apply button state
        self.onSelect()

    def cleanup(self):
        pass

    def onSelect(self):
        self.initBPButton.enabled = self.heartValveSelector.currentNode() and self.inputSelector.currentNode() and self.outputSegmentationSelector.currentNode()
        # self.generateSurfaceMarkups.enabled = self.heartValveSelector.currentNode() and self.outputSegmentationSelector.currentNode() and self.markupsSelector.currentNode()
        self.generateMoldButton.enabled = self.heartValveSelector.currentNode() and self.outputSegmentationSelector.currentNode() and self.inputSelector.currentNode()
        self.exportMoldButton.enabled = self.heartValveSelector.currentNode() and self.outputSegmentationSelector.currentNode() \
                                        and self.outputSegmentationSelector.currentNode().GetClosedSurfaceRepresentation('Final_Mold')

    def onInitBPButton(self):
        try:
            # This can be a long operation - indicate it to the user
            qt.QApplication.setOverrideCursor(qt.Qt.WaitCursor)

            enableScreenshotsFlag = self.enableScreenshotsFlagCheckBox.checked
            self.logic.initBPSeg(self.inputSelector.currentNode(), self.heartValveSelector.currentNode(),
                                 self.outputSegmentationSelector.currentNode())
            self.incrementFirstButton50.enabled = True
            self.incrementFirstButton100.enabled = True
            self.incrementFirstButton500.enabled = True
            self.initLeafletButton.enabled = True
        finally:
            qt.QApplication.restoreOverrideCursor()

    def onIncrementFirst50Button(self):
        try:
            # This can be a long operation - indicate it to the user
            qt.QApplication.setOverrideCursor(qt.Qt.WaitCursor)

            self.logic.iterateFirstPass(50, self.outputSegmentationSelector.currentNode())
            self.undoButtonBP.enabled = True

        finally:
            qt.QApplication.restoreOverrideCursor()

    def onIncrementFirst100Button(self):
        try:
            # This can be a long operation - indicate it to the user
            qt.QApplication.setOverrideCursor(qt.Qt.WaitCursor)

            self.logic.iterateFirstPass(100, self.outputSegmentationSelector.currentNode())
            self.undoButtonBP.enabled = True

        finally:
            qt.QApplication.restoreOverrideCursor()

    def onIncrementFirst500Button(self):
        try:
            # This can be a long operation - indicate it to the user
            qt.QApplication.setOverrideCursor(qt.Qt.WaitCursor)

            self.logic.iterateFirstPass(500, self.outputSegmentationSelector.currentNode())
            self.undoButtonBP.enabled = True

        finally:
            qt.QApplication.restoreOverrideCursor()

    def onInitLeafletButton(self):
        try:
            # This can be a long operation - indicate it to the user
            qt.QApplication.setOverrideCursor(qt.Qt.WaitCursor)

            self.logic.initLeafletSeg(self.outputSegmentationSelector.currentNode())
            self.incrementButton10.enabled = True
            self.incrementButton50.enabled = True
            self.incrementButton200.enabled = True

        finally:
            qt.QApplication.restoreOverrideCursor()

    def onIncrement10Button(self):
        self.logic.iterateSecondPass(10, self.outputSegmentationSelector.currentNode())
        self.undoButtonLeaflet.enabled = True

    def onIncrement50Button(self):
        try:
            # This can be a long operation - indicate it to the user
            qt.QApplication.setOverrideCursor(qt.Qt.WaitCursor)

            self.logic.iterateSecondPass(50, self.outputSegmentationSelector.currentNode())
            self.undoButtonLeaflet.enabled = True

        finally:
            qt.QApplication.restoreOverrideCursor()

    def onIncrement200Button(self):
        try:
            # This can be a long operation - indicate it to the user
            qt.QApplication.setOverrideCursor(qt.Qt.WaitCursor)

            self.logic.iterateSecondPass(200, self.outputSegmentationSelector.currentNode())
            self.undoButtonLeaflet.enabled = True

        finally:
            qt.QApplication.restoreOverrideCursor()

    def onUndoButtonBP(self):
        self.logic.undoBPIteration(self.outputSegmentationSelector.currentNode())
        self.redoButtonBP.enabled = True
        self.undoButtonBP.enabled = False

    def onRedoButtonBP(self):
        self.logic.redoBPIteration(self.outputSegmentationSelector.currentNode())
        self.redoButtonBP.enabled = False
        self.undoButtonBP.enabled = True

    def onUndoButtonLeaflet(self):
        self.logic.undoLeafletIteration(self.outputSegmentationSelector.currentNode())
        self.redoButtonLeaflet.enabled = True
        self.undoButtonLeaflet.enabled = False

    def onRedoButtonLeaflet(self):
        self.logic.redoLeafletIteration(self.outputSegmentationSelector.currentNode())
        self.redoButtonLeaflet.enabled = False
        self.undoButtonLeaflet.enabled = True

    def onGenerateSurfaceMarkups(self):
        success = self.logic.generateSurfaceMarkups(self.outputSegmentationSelector.currentNode(),
                                                    self.heartValveSelector.currentNode(),
                                                    self.markupsSelector.currentNode())

    def onGenerateModelButton(self):
        try:
            # This can be a long operation - indicate it to the user
            qt.QApplication.setOverrideCursor(qt.Qt.WaitCursor)

            self.logic.generateSurfaceMold(self.outputSegmentationSelector.currentNode(),
                                           self.heartValveSelector.currentNode(),
                                           float(self.baseDepthSlider.value),
                                           self.inputSelector.currentNode())
            self.onSelect()

        finally:
            qt.QApplication.restoreOverrideCursor()

    def onExportModelButton(self):
        try:
            # This can be a long operation - indicate it to the user
            qt.QApplication.setOverrideCursor(qt.Qt.WaitCursor)

            self.logic.exportSurfaceMold(self.outputSegmentationSelector.currentNode())

        finally:
            qt.QApplication.restoreOverrideCursor()

    def onGenerateBasePlate(self):
        return

#
# MVSegmenterLogic
#

class MVSegmenterLogic(ScriptedLoadableModuleLogic):
    """This class should implement all the actual
    computation done by your module.  The interface
    should be such that other python code can import
    this class and make use of the functionality without
    requiring an instance of the Widget.
    Uses ScriptedLoadableModuleLogic base class, available at:
    https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
    """

    def __init__(self):
        ScriptedLoadableModuleLogic.__init__(self)
        self._speedImg = None
        self._levelSet = None
        self._bpLevelSet = None

        self._nextBpLevelSet = None
        self._prevBpLevelSet = None
        self._prevLeafletLevelSet = None
        self._nextLeafletLevelSet = None

        self.moldBasePlate = None

    def hasImageData(self, volumeNode):
        """This is an example logic method that
        returns true if the passed in volume
        node has valid image data
        """
        if not volumeNode:
            logging.debug('hasImageData failed: no volume node')
            return False
        if volumeNode.GetImageData() is None:
            logging.debug('hasImageData failed: no image data in volume node')
            return False
        return True

    def isValidInputOutputData(self, inputVolumeNode, outputVolumeNode):
        """Validates if the output is not the same as input
        """
        if not inputVolumeNode:
            logging.debug('isValidInputOutputData failed: no input volume node defined')
            return False
        if not outputVolumeNode:
            logging.debug('isValidInputOutputData failed: no output volume node defined')
            return False
        if inputVolumeNode.GetID() == outputVolumeNode.GetID():
            logging.debug(
                'isValidInputOutputData failed: input and output volume is the same. Create a new volume for output to avoid this error.')
            return False
        return True

    def takeScreenshot(self, name, description, type=-1):
        # show the message even if not taking a screen shot
        slicer.util.delayDisplay(
            'Take screenshot: ' + description + '.\nResult is available in the Annotations module.', 3000)

        lm = slicer.app.layoutManager()
        # switch on the type to get the requested window
        widget = 0
        if type == slicer.qMRMLScreenShotDialog.FullLayout:
            # full layout
            widget = lm.viewport()
        elif type == slicer.qMRMLScreenShotDialog.ThreeD:
            # just the 3D window
            widget = lm.threeDWidget(0).threeDView()
        elif type == slicer.qMRMLScreenShotDialog.Red:
            # red slice window
            widget = lm.sliceWidget("Red")
        elif type == slicer.qMRMLScreenShotDialog.Yellow:
            # yellow slice window
            widget = lm.sliceWidget("Yellow")
        elif type == slicer.qMRMLScreenShotDialog.Green:
            # green slice window
            widget = lm.sliceWidget("Green")
        else:
            # default to using the full window
            widget = slicer.util.mainWindow()
            # reset the type so that the node is set correctly
            type = slicer.qMRMLScreenShotDialog.FullLayout

        # grab and convert to vtk image data
        qimage = ctk.ctkWidgetsUtils.grabWidget(widget)
        imageData = vtk.vtkImageData()
        slicer.qMRMLUtils().qImageToVtkImageData(qimage, imageData)

        annotationLogic = slicer.modules.annotations.logic()
        annotationLogic.CreateSnapShot(name, description, type, 1, imageData)

    def initBPSeg(self, inputVolume, heartValveNode, outputSeg):
        """
        Run the actual algorithm
        :param heartValveNode:
        """

        if not inputVolume or not heartValveNode or not outputSeg:
            logging.error("Missing parameter")
            return

        valveModel = HeartValveLib.getValveModel(heartValveNode)
        if valveModel.getAnnulusContourMarkupNode().GetNumberOfFiducials() == 0:
            logging.error("Annulus contour not defined")
            return

        if valveModel.getProbeToRasTransformNode():
            outputSeg.SetAndObserveTransformNodeID(valveModel.getProbeToRasTransformNode().GetID())

        # calculate speed image from input volume
        # Uses DiscreteGaussian -> GradientMagnitude -> Sigmoid filters
        speedImg = sitkUtils.PullVolumeFromSlicer(inputVolume)

        blurFilter = sitk.DiscreteGaussianImageFilter()
        blurFilter.SetMaximumError(0.25)
        blurFilter.SetMaximumKernelWidth(32)
        blurFilter.SetUseImageSpacing(True)
        blurFilter.SetVariance(1.5)
        speedImg = blurFilter.Execute(speedImg)

        gradMag = sitk.GradientMagnitudeImageFilter()
        gradMag.SetUseImageSpacing(True)
        speedImg = gradMag.Execute(speedImg)

        sigmoid = sitk.SigmoidImageFilter()
        sigmoid.SetOutputMinimum(0)
        sigmoid.SetOutputMaximum(1.0)
        sigmoid.SetAlpha(-5.0)
        sigmoid.SetBeta(10.0)
        speedImg = sigmoid.Execute(speedImg)

        self._speedImg = speedImg

        # compute initial level set

        # Find annulus center
        markups = valveModel.getAnnulusContourMarkupNode()
        pos = np.zeros(3)
        centroid = np.zeros(3)
        for i in range(markups.GetNumberOfFiducials()):
            markups.GetNthFiducialPosition(i, pos)
            centroid += pos

        centroid = centroid / markups.GetNumberOfFiducials()
        centroid = centroid + valveModel.getAnnulusContourPlane()[1] * -10
        centroid = np.array(self.rasToIJK(centroid, inputVolume))
        if centroid[2] <= 0:
            centroid[2] = 1

        # Run fast marching based on annulus center
        fastMarching = sitk.FastMarchingImageFilter()
        fastMarching.SetTrialPoints([centroid.astype('uint32').tolist()])
        fmarch = fastMarching.Execute(speedImg)

        thresh = sitk.BinaryThresholdImageFilter()
        thresh.SetLowerThreshold(0)
        thresh.SetUpperThreshold(10)
        thresh.SetInsideValue(1)
        thresh.SetOutsideValue(0)
        mask = thresh.Execute(fmarch)

        signedDis = sitk.SignedDanielssonDistanceMapImageFilter()
        levelSet = signedDis.Execute(mask)

        # Run first pass of geodesic active contour
        geodesicActiveContour = sitk.GeodesicActiveContourLevelSetImageFilter()
        geodesicActiveContour.SetCurvatureScaling(0.8)
        geodesicActiveContour.SetAdvectionScaling(1.2)
        geodesicActiveContour.SetPropagationScaling(1.0)
        geodesicActiveContour.SetMaximumRMSError(0.0001)
        geodesicActiveContour.SetNumberOfIterations(500)

        out_mask = geodesicActiveContour.Execute(levelSet, speedImg)

        self._bpLevelSet = out_mask

        threshold = sitk.BinaryThresholdImageFilter()
        threshold.SetInsideValue(1)
        threshold.SetLowerThreshold(-1000.0)
        threshold.SetOutsideValue(0)
        threshold.SetUpperThreshold(0.0)
        out_mask = threshold.Execute(out_mask)

        self.pushITKImageToSegmentation(out_mask, outputSeg, 'BP Segmentation')

    def iterateFirstPass(self, nIter, outputSeg):
        geodesicActiveContour = sitk.GeodesicActiveContourLevelSetImageFilter()
        geodesicActiveContour.SetCurvatureScaling(1.2)
        geodesicActiveContour.SetAdvectionScaling(1.0)
        geodesicActiveContour.SetPropagationScaling(0.9)
        geodesicActiveContour.SetMaximumRMSError(0.00001)
        geodesicActiveContour.SetNumberOfIterations(nIter)

        out_mask = geodesicActiveContour.Execute(self._bpLevelSet, self._speedImg)

        self._prevBpLevelSet = self._bpLevelSet
        self._bpLevelSet = out_mask

        threshold = sitk.BinaryThresholdImageFilter()
        threshold.SetInsideValue(1)
        threshold.SetLowerThreshold(-1000.0)
        threshold.SetOutsideValue(0)
        threshold.SetUpperThreshold(0.0)
        out_mask = threshold.Execute(out_mask)

        self.pushITKImageToSegmentation(out_mask, outputSeg, 'BP Segmentation')

    def initLeafletSeg(self, outputSeg):

        # Get region bordering initial blood-pool segmentation
        threshold = sitk.BinaryThresholdImageFilter()
        threshold.SetInsideValue(1)
        threshold.SetLowerThreshold(-1000.0)
        threshold.SetOutsideValue(0)
        threshold.SetUpperThreshold(0.0)
        distMap = threshold.Execute(self._bpLevelSet)

        signedDis = sitk.SignedDanielssonDistanceMapImageFilter()
        distMap = signedDis.Execute(distMap)

        distThreshold = sitk.BinaryThresholdImageFilter()
        distThreshold.SetInsideValue(1)
        distThreshold.SetLowerThreshold(1)
        distThreshold.SetOutsideValue(0)
        distThreshold.SetUpperThreshold(11)
        leafletMask = distThreshold.Execute(distMap)

        # Run second pass to get final leaflet segmentation

        levelSet = signedDis.Execute(leafletMask)

        geodesicActiveContour2 = sitk.GeodesicActiveContourLevelSetImageFilter()
        geodesicActiveContour2.SetCurvatureScaling(1.0)
        geodesicActiveContour2.SetAdvectionScaling(0.1)
        geodesicActiveContour2.SetPropagationScaling(-0.6)
        geodesicActiveContour2.SetMaximumRMSError(0.0001)
        geodesicActiveContour2.SetNumberOfIterations(300)
        out_mask = geodesicActiveContour2.Execute(levelSet, self._speedImg)

        self._levelSet = out_mask

        threshold = sitk.BinaryThresholdImageFilter()
        threshold.SetInsideValue(1)
        threshold.SetLowerThreshold(-1000.0)
        threshold.SetOutsideValue(0)
        threshold.SetUpperThreshold(0.0)
        out_mask = threshold.Execute(out_mask)

        self.pushITKImageToSegmentation(out_mask, outputSeg, 'Leaflet Segmentation')

        return out_mask

    def iterateSecondPass(self, nIter, outputSeg):

        geodesicActiveContour2 = sitk.GeodesicActiveContourLevelSetImageFilter()
        geodesicActiveContour2.SetCurvatureScaling(0.9)
        geodesicActiveContour2.SetAdvectionScaling(0.1)
        geodesicActiveContour2.SetPropagationScaling(-0.4)
        geodesicActiveContour2.SetMaximumRMSError(0.0001)
        geodesicActiveContour2.SetNumberOfIterations(nIter)
        out_mask = geodesicActiveContour2.Execute(self._levelSet, self._speedImg)

        self._prevLeafletLevelSet = self._levelSet
        self._levelSet = out_mask

        threshold = sitk.BinaryThresholdImageFilter()
        threshold.SetInsideValue(1)
        threshold.SetLowerThreshold(-1000.0)
        threshold.SetOutsideValue(0)
        threshold.SetUpperThreshold(0.0)
        out_mask = threshold.Execute(out_mask)

        self.pushITKImageToSegmentation(out_mask, outputSeg, 'Leaflet Segmentation')

        return out_mask

    def undoBPIteration(self, outputSeg):
        self._nextBpLevelSet = self._bpLevelSet
        self._bpLevelSet = self._prevBpLevelSet

        threshold = sitk.BinaryThresholdImageFilter()
        threshold.SetInsideValue(1)
        threshold.SetLowerThreshold(-1000.0)
        threshold.SetOutsideValue(0)
        threshold.SetUpperThreshold(0.0)

        self.pushITKImageToSegmentation(threshold.Execute(self._bpLevelSet), outputSeg, 'BP Segmentation')

    def redoBPIteration(self, outputSeg):
        self._prevBpLevelSet = self._bpLevelSet
        self._bpLevelSet = self._nextBpLevelSet

        threshold = sitk.BinaryThresholdImageFilter()
        threshold.SetInsideValue(1)
        threshold.SetLowerThreshold(-1000.0)
        threshold.SetOutsideValue(0)
        threshold.SetUpperThreshold(0.0)

        self.pushITKImageToSegmentation(threshold.Execute(self._bpLevelSet), outputSeg, 'BP Segmentation')

    def undoLeafletIteration(self, outputSeg):
        self._nextLeafletLevelSet = self._levelSet
        self._levelSet = self._prevLeafletLevelSet

        threshold = sitk.BinaryThresholdImageFilter()
        threshold.SetInsideValue(1)
        threshold.SetLowerThreshold(-1000.0)
        threshold.SetOutsideValue(0)
        threshold.SetUpperThreshold(0.0)

        self.pushITKImageToSegmentation(threshold.Execute(self._levelSet), outputSeg, 'Leaflet Segmentation')

        return self._levelSet

    def redoLeafletIteration(self, outputSeg):
        self._prevLeafletLevelSet = self._levelSet
        self._levelSet = self._nextLeafletLevelSet

        threshold = sitk.BinaryThresholdImageFilter()
        threshold.SetInsideValue(1)
        threshold.SetLowerThreshold(-1000.0)
        threshold.SetOutsideValue(0)
        threshold.SetUpperThreshold(0.0)

        self.pushITKImageToSegmentation(threshold.Execute(self._levelSet), outputSeg, 'Leaflet Segmentation')

        return self._levelSet

    def pushITKImageToSegmentation(self, img, segmentationNode, segmentId='Leaflet Segmentation'):
        if segmentationNode.GetSegmentation().GetSegmentIndex(segmentId) == -1:
            segmentationNode.GetSegmentation().AddEmptySegment(segmentId)

        if segmentationNode.GetSegmentation().GetConversionParameter('Decimation factor') != '0.0':
            segmentationNode.GetSegmentation().SetConversionParameter('Decimation factor', '0.0')
            segmentationNode.RemoveClosedSurfaceRepresentation()
            segmentationNode.CreateClosedSurfaceRepresentation()

        if segmentationNode.GetSegmentation().GetConversionParameter('Smoothing factor') != '0.5':
            segmentationNode.GetSegmentation().SetConversionParameter('Smoothing factor', '0.5')
            segmentationNode.RemoveClosedSurfaceRepresentation()
            segmentationNode.CreateClosedSurfaceRepresentation()

        # Create temporary label map node
        tempNode = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLLabelMapVolumeNode', 'temp_labelmap')
        tempNode.SetAndObserveTransformNodeID(segmentationNode.GetTransformNodeID())
        sitkUtils.PushVolumeToSlicer(img, tempNode)

        segmentationIds = vtk.vtkStringArray()
        segmentationIds.InsertNextValue(segmentId)

        slicer.modules.segmentations.logic().ImportLabelmapToSegmentationNode(tempNode, segmentationNode,
                                                                              segmentationIds)

        slicer.mrmlScene.RemoveNode(tempNode)

    def rasToIJK(self, point, volume):
        matrix = vtk.vtkMatrix4x4()
        volume.GetRASToIJKMatrix(matrix)

        homogeneousPoint = [point[0], point[1], point[2], 1]
        outPoint = matrix.MultiplyPoint(homogeneousPoint)

        return outPoint[0:3]

    def generateSurfaceMarkups(self, segNode, heartValveNode, markupsNode):
        if not segNode or not heartValveNode or not markupsNode:
            logging.error("Missing parameter")
            return False

        valveModel = HeartValveLib.getValveModel(heartValveNode)
        if valveModel.getAnnulusContourMarkupNode().GetNumberOfFiducials() == 0:
            logging.error("Annulus contour not defined")
            return False

        leafletModel = segNode.GetClosedSurfaceRepresentation('Leaflet Segmentation')
        if leafletModel is None:
            logging.error("Missing segmentation")
            return False

        obb = vtk.vtkOBBTree()
        obb.SetDataSet(leafletModel)
        obb.BuildLocator()

        markupsNode.RemoveAllMarkups()
        markupsNode.SetAndObserveTransformNodeID(segNode.GetTransformNodeID())

        fixedMarkupsNode = markupsNode.GetNodeReference('fixedNodeRef')
        if not fixedMarkupsNode:
            fixedMarkupsNode = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLMarkupsFiducialNode',
                                                                  markupsNode.GetName() + '-fixed')
            markupsNode.AddNodeReferenceID('fixedNodeRef', fixedMarkupsNode.GetID())

        fixedMarkupsNode.RemoveAllMarkups()

        # Get points from leaflet model where line from point to annulus center does not self intersect
        a0 = np.zeros(3)
        contourPlane = valveModel.getAnnulusContourPlane()
        points = vtk.vtkPoints()
        ids = vtk.vtkIdTypeArray()
        for i in range(leafletModel.GetNumberOfPoints()):
            leafletModel.GetPoint(i, a0)
            r = obb.IntersectWithLine(a0, contourPlane[0], points, None)
            # If only 1 intersection point, line does not cross through leaflet model as the line always intersects at a0
            if points.GetNumberOfPoints() == 1:
                index = ids.InsertNextValue(i)

        # Extract the points lying on the inside of the model
        selectionNode = vtk.vtkSelectionNode()
        selectionNode.SetFieldType(vtk.vtkSelectionNode.POINT)
        selectionNode.SetContentType(vtk.vtkSelectionNode.INDICES)
        selectionNode.SetSelectionList(ids)

        selection = vtk.vtkSelection()
        selection.AddNode(selectionNode)

        extractSelection = vtk.vtkExtractSelection()
        extractSelection.SetInputData(0, leafletModel)
        extractSelection.SetInputData(1, selection)
        extractSelection.Update()

        geom = vtk.vtkGeometryFilter()
        geom.SetInputConnection(extractSelection.GetOutputPort())
        geom.Update()

        # Cleaner ensures points are spaced out evenly
        cleaner = vtk.vtkCleanPolyData()
        cleaner.SetTolerance(0.07)
        cleaner.SetInputConnection(geom.GetOutputPort())
        cleaner.Update()

        # Populate markups node with model points
        mNodeModifyState = markupsNode.StartModify()
        fNodeModifyState = fixedMarkupsNode.StartModify()
        points = cleaner.GetOutput().GetPoints()
        for i in range(points.GetNumberOfPoints()):
            markupsNode.AddFiducialFromArray(points.GetPoint(i))
            fixedMarkupsNode.AddFiducialFromArray(points.GetPoint(i))

        markupsNode.EndModify(mNodeModifyState)
        fixedMarkupsNode.EndModify(fNodeModifyState)
        
        markupsNode.GetMarkupsDisplayNode().SetTextScale(0)

        fixedMarkupsNode.GetMarkupsDisplayNode().SetVisibility(0)
        fixedMarkupsNode.SetLocked(1)

        return True

    def pushModelToSegmentation(self, segNode, model, name):
        import vtkSegmentationCorePython as vtkSegmentationCore

        if segNode.GetSegmentation().GetSegmentIndex(name) != -1:
            segNode.GetSegmentation().RemoveSegment(name)
        # Add segment from polydata
        segment = vtkSegmentationCore.vtkSegment()
        segment.SetName(name)
        segment.AddRepresentation(
            vtkSegmentationCore.vtkSegmentationConverter.GetSegmentationClosedSurfaceRepresentationName(),
            model)
        segNode.GetSegmentation().AddSegment(segment, name)



    def generateSurfaceMold(self, segNode, heartValveNode, depth, volume):
        # Check that parameters exist
        if not segNode or not heartValveNode or not depth or not volume:
            logging.error("Missing parameter")
            return None

        valveModel = HeartValveLib.getValveModel(heartValveNode)

        # Get the proximal surface of the valve
        extractedSurface = self.extractInnerSurfaceModel(segNode, valveModel)
        if not extractedSurface:
            return None

        # Get the annulus projected onto the proximal surface
        projectedAnnulus = self.generateProjectedAnnulus(extractedSurface, valveModel)
        if not projectedAnnulus:
            return None

        self.pushModelToSegmentation(segNode, projectedAnnulus, 'Projected_Annulus')

        contourPlane = valveModel.getAnnulusContourPlane()

        # Create clipped leaflet mold across middle
        baseClippingPlane = vtk.vtkPlane()
        baseClippingPlane.SetNormal(contourPlane[1])
        baseClippingPlane.SetOrigin(contourPlane[0] + contourPlane[1] * depth)

        midClippingPlane = vtk.vtkPlane()
        midClippingPlane.SetNormal(contourPlane[1])
        midClippingPlane.SetOrigin(contourPlane[0] + (contourPlane[1] * 2))

        topMold, bottomMold = self.buildMoldHalves(extractedSurface, midClippingPlane, baseClippingPlane)

        # Put top, bottom and base of mold together

        append = vtk.vtkAppendPolyData()
        append.AddInputData(topMold)
        append.AddInputData(bottomMold)
        append.Update()

        normAuto = vtk.vtkPolyDataNormals()
        normAuto.ConsistencyOn()
        normAuto.SetInputConnection(append.GetOutputPort())
        normAuto.Update()

        mold = vtk.vtkPolyData()
        mold.DeepCopy(normAuto.GetOutput())

        self.pushModelToSegmentation(segNode, mold, 'Final_Mold')

        # Remake closed surface representation after adding mold (makes it generated model from labelmap)
        segNode.RemoveClosedSurfaceRepresentation()
        segNode.CreateClosedSurfaceRepresentation()

        # Create segment editor to get access to effects
        segmentEditorWidget = slicer.qMRMLSegmentEditorWidget()
        segmentEditorWidget.setMRMLScene(slicer.mrmlScene)
        segmentEditorNode = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLSegmentEditorNode")
        segmentEditorNode.SetOverwriteMode(segmentEditorNode.OverwriteNone)
        segmentEditorWidget.setMRMLSegmentEditorNode(segmentEditorNode)
        segmentEditorWidget.setSegmentationNode(segNode)
        segmentEditorWidget.setMasterVolumeNode(volume)
        segmentEditorWidget.setCurrentSegmentID('Final_Mold')

        # Smoothing
        segmentEditorWidget.setActiveEffectByName("Smoothing")
        effect = segmentEditorWidget.activeEffect()
        effect.setParameter("SmoothingMethod", "GAUSSIAN")
        effect.setParameter("GaussianStandardDeviationMm", 0.8)
        effect.self().onApply()

        # Clean up
        segmentEditorWidget = None
        slicer.mrmlScene.RemoveNode(segmentEditorNode)

    def exportSurfaceMold(self, segNode):

        # Get segmentation closed surface representations
        segMold = segNode.GetClosedSurfaceRepresentation('Final_Mold')
        annulusMold = segNode.GetClosedSurfaceRepresentation('Projected_Annulus')
        if not segMold or not annulusMold:
            logging.error("Missing mold segmentation")
            return

        # Fill holes (remove boudnary edges) and clean
        holeFill = vtk.vtkFillHolesFilter()
        holeFill.SetInputData(segMold)
        holeFill.SetHoleSize(holeFill.GetHoleSizeMaxValue())
        holeFill.Update()

        clean = vtk.vtkCleanPolyData()
        clean.SetInputConnection(holeFill.GetOutputPort())
        clean.Update()

        # Remove disconnected fragments
        conn = vtk.vtkConnectivityFilter()
        conn.SetInputConnection(clean.GetOutputPort())
        conn.SetExtractionModeToLargestRegion()
        conn.Update()

        normAuto = vtk.vtkPolyDataNormals()
        normAuto.ConsistencyOff()
        normAuto.SetInputConnection(conn.GetOutputPort())
        normAuto.Update()

        # Decimate to 80%
        decimate = vtk.vtkDecimatePro()
        decimate.SetTargetReduction(0.8)
        decimate.SetInputConnection(normAuto.GetOutputPort())
        decimate.Update()

        holeFill = vtk.vtkFillHolesFilter()
        holeFill.SetInputConnection(decimate.GetOutputPort())
        holeFill.SetHoleSize(holeFill.GetHoleSizeMaxValue())
        holeFill.Update()

        normAuto = vtk.vtkPolyDataNormals()
        normAuto.ConsistencyOn()
        normAuto.AutoOrientNormalsOn()
        normAuto.SetInputConnection(holeFill.GetOutputPort())
        normAuto.Update()


        # Perform boolean operations to get the base mold with the annulus projection
        boolFilter = vtk.vtkLoopBooleanPolyDataFilter()
        boolFilter.SetInputConnection(0, normAuto.GetOutputPort())
        boolFilter.SetInputData(1, annulusMold)

        boolFilter.SetOperationToDifference()
        boolFilter.Update()

        moldModel = vtk.vtkPolyData()
        moldModel.DeepCopy(boolFilter.GetOutput())

        self.addOrUpdateModel(moldModel, 'Final_Mold_Model', segNode.GetTransformNodeID())

        boolFilter.SetOperationToIntersection()
        boolFilter.Update()

        annulusModel = vtk.vtkPolyData()
        annulusModel.DeepCopy(boolFilter.GetOutput())
        self.addOrUpdateModel(annulusModel, 'Projected_Annulus_Model', segNode.GetTransformNodeID())


    def extractInnerSurfaceModel(self, segNode, valveModel, segName = 'Leaflet Segmentation'):
        if not segNode or not valveModel:
            logging.error("Missing parameter")
            return None

        if valveModel.getAnnulusContourMarkupNode().GetNumberOfFiducials() == 0:
            logging.error("Annulus contour not defined")
            return None

        annulusPlane = valveModel.getAnnulusContourPlane()


        leafletModel = segNode.GetClosedSurfaceRepresentation(segName)
        bpModel = segNode.GetClosedSurfaceRepresentation('BP Segmentation')
        if leafletModel is None or bpModel is None:
            logging.error("Missing segmentation")
            return None

        # Decimate leaflet and BP polydata for efficiency
        decimate = vtk.vtkDecimatePro()
        decimate.SetTargetReduction(0.5)
        decimate.PreserveTopologyOn()
        decimate.BoundaryVertexDeletionOff()

        decimate.SetInputData(leafletModel)
        decimate.Update()
        leafletModel = vtk.vtkPolyData()
        leafletModel.DeepCopy(decimate.GetOutput())

        decimate.SetInputData(bpModel)
        decimate.Update()
        bpModel = vtk.vtkPolyData()
        bpModel.DeepCopy(decimate.GetOutput())

        # Extract surface within a distance of 4 from bp segmentation
        # Done for performance reasons as it quickly eliminates polydata that we do not need and makes the loop phase faster
        imp = vtk.vtkImplicitPolyDataDistance()
        imp.SetInput(bpModel)

        clip = vtk.vtkClipPolyData()
        clip.SetClipFunction(imp)
        clip.GenerateClipScalarsOn()
        clip.InsideOutOn()
        clip.SetValue(4.0)
        clip.SetInputData(leafletModel)
        clip.Update()

        clipped = vtk.vtkPolyData()
        clipped.DeepCopy(clip.GetOutput())

        # OBBTree for determining self intersection of rays
        obb = vtk.vtkOBBTree()
        obb.SetDataSet(leafletModel)
        obb.BuildLocator()

        locator = vtk.vtkPointLocator()
        locator.SetDataSet(bpModel)
        locator.BuildLocator()

        # Loop over remaining points and build scalar array using different techniques for top and bottom half
        a0 = np.zeros(3)
        p = annulusPlane[0] + annulusPlane[1] * 2
        points = vtk.vtkPoints()
        normals = clipped.GetPointData().GetNormals()
        bpNormals = bpModel.GetPointData().GetNormals()
        scalars = vtk.vtkFloatArray()
        scalars.SetNumberOfValues(clipped.GetNumberOfPoints())
        for i in range(clipped.GetNumberOfPoints()):
            clipped.GetPoint(i, a0)
            if np.dot(annulusPlane[1], a0 - p) > 0:
                # Point is above annulus plane, set scalar based on self intersection (scalar value will determine clipping)

                r = obb.IntersectWithLine(a0, annulusPlane[0], points, None)
                # If only 1 intersection point, line does not cross through leaflet model as the line always intersects at a0
                if points.GetNumberOfPoints() == 1:
                    scalars.SetValue(i, 10)     # Set scalar to large value so point will be kept
                else:
                    scalars.SetValue(i, -10)    # Set scalar to small value so point will be discarded
            else:
                # Point is below annulus plane
                # Get the closest point on bp surface, find angle between 2 normals in radians
                closestPoint = locator.FindClosestPoint(a0)
                v = np.array(bpNormals.GetTuple(closestPoint))
                n = np.array(normals.GetTuple(i))
                angle = math.acos(np.dot(n, v) / np.linalg.norm(n) / np.linalg.norm(v))
                scalars.SetValue(i, angle)      # Set scalar to angle

        # Scalars now angles in radians that we can threshold
        clipped.GetPointData().SetScalars(scalars)

        # Clip based on scalar values (keep scalars bigger than value)
        clip2 = vtk.vtkClipPolyData()
        clip2.GenerateClipScalarsOff()
        clip2.SetValue(1.74533)  # 100 degrees threshold in radians
        clip2.SetInputData(clipped)

        conn = vtk.vtkConnectivityFilter()
        conn.SetInputConnection(clip2.GetOutputPort())
        conn.SetExtractionModeToLargestRegion()
        conn.Update()

        # Fill small holes resulting from extraction and clean poly data
        fill = vtk.vtkFillHolesFilter()
        fill.SetHoleSize(3)
        fill.SetInputConnection(conn.GetOutputPort())
        fill.Update()

        clean = vtk.vtkCleanPolyData()
        clean.SetInputConnection(fill.GetOutputPort())
        clean.Update()

        # Fix normals
        normClean = vtk.vtkPolyDataNormals()
        normClean.ConsistencyOn()
        normClean.FlipNormalsOn()
        normClean.SetInputConnection(clean.GetOutputPort())
        normClean.Update()

        innerModel = vtk.vtkPolyData()
        innerModel.DeepCopy(normClean.GetOutput())
        return innerModel


    def generateProjectedAnnulus(self, extractedLeaflet, valveModel):
        if not extractedLeaflet or not valveModel:
            logging.error("Missing parameter")
            return None

        if valveModel.getAnnulusContourMarkupNode().GetNumberOfFiducials() == 0:
            logging.error("Annulus contour not defined")
            return None

        # Project defined annulus onto inner surface

        # Use OBBTree to find intersection with model
        obb = vtk.vtkOBBTree()
        obb.SetDataSet(extractedLeaflet)
        obb.BuildLocator()

        annulusMarkups = valveModel.getAnnulusContourMarkupNode()
        contourPlane = valveModel.getAnnulusContourPlane()
        pos = np.zeros(3)
        points = vtk.vtkPoints()
        projPoints = vtk.vtkPoints()

        # Project annulus inwards towards center
        center = contourPlane[0]
        for i in range(annulusMarkups.GetNumberOfFiducials()):
            annulusMarkups.GetNthFiducialPosition(i, pos)
            r = obb.IntersectWithLine(pos + pos - center, center, points, None)
            if r != 0:
                projPoints.InsertNextPoint(points.GetPoint(0))

        lines = vtk.vtkCellArray()
        lines.InsertNextCell(projPoints.GetNumberOfPoints())
        for i in range(projPoints.GetNumberOfPoints()):
            lines.InsertCellPoint(i)

        # Create spline polydata
        projContour = vtk.vtkPolyData()
        projContour.SetPoints(projPoints)
        projContour.SetLines(lines)

        splineFilter = vtk.vtkSplineFilter()
        splineFilter.SetNumberOfSubdivisions(500)
        splineFilter.GetSpline().ClosedOn()
        splineFilter.SetInputData(projContour)
        splineFilter.Update()

        # Close spline
        strip = vtk.vtkStripper()
        strip.SetInputConnection(splineFilter.GetOutputPort())
        strip.Update()

        # Create tube from spline fitted projected annulus
        tubeFilter = vtk.vtkTubeFilter()
        tubeFilter.SetRadius(1) # Radius of 1 determined through trial and error on printed models
        tubeFilter.SetNumberOfSides(20)
        tubeFilter.CappingOff()
        tubeFilter.SetInputConnection(strip.GetOutputPort())
        tubeFilter.Update()

        cleanTube = vtk.vtkCleanPolyData()
        cleanTube.SetInputConnection(tubeFilter.GetOutputPort())
        cleanTube.Update()

        # Push fitted annulus onto segmentation node
        annulusFittedModel = vtk.vtkPolyData()
        annulusFittedModel.DeepCopy(cleanTube.GetOutput())
        return annulusFittedModel

    def buildMoldHalves(self, extractedSurface, midClippingPlane, baseClippingPlane):

        # Split top and bottom halves of surface
        clipMid = vtk.vtkClipPolyData()
        clipMid.SetClipFunction(midClippingPlane)
        clipMid.SetInputData(extractedSurface)
        clipMid.GenerateClippedOutputOn()
        clipMid.Update()

        # Fill across mid clip plane
        cutter = vtk.vtkCutter()
        cutter.SetCutFunction(midClippingPlane)
        cutter.SetInputData(extractedSurface)
        cutter.Update()

        loop = vtk.vtkContourLoopExtraction()
        loop.SetNormal(midClippingPlane.GetNormal())
        loop.SetLoopClosureToAll()
        loop.SetInputConnection(cutter.GetOutputPort())
        loop.Update()

        tri = vtk.vtkTriangleFilter()
        tri.SetInputConnection(loop.GetOutputPort())
        tri.Update()

        # Extrusion towards annulus centroid to thicken leaflet walls inwards

        clean = vtk.vtkCleanPolyData()
        clean.SetInputConnection(clipMid.GetOutputPort())
        clean.Update()

        extrudeIn = vtk.vtkLinearExtrusionFilter()
        extrudeIn.CappingOn()
        extrudeIn.SetExtrusionTypeToPointExtrusion()
        extrudeIn.SetScaleFactor(-0.6)
        extrudeIn.SetExtrusionPoint(midClippingPlane.GetOrigin())
        extrudeIn.SetInputConnection(clean.GetOutputPort())
        extrudeIn.Update()

        normAuto = vtk.vtkPolyDataNormals()
        normAuto.ConsistencyOn()
        normAuto.AutoOrientNormalsOn()
        normAuto.SetInputConnection(extrudeIn.GetOutputPort())
        normAuto.Update()

        # Need clean then fill to close extruded model
        clean = vtk.vtkCleanPolyData()
        clean.SetInputConnection(normAuto.GetOutputPort())
        clean.Update()

        # Make normals point outwards for final model
        normAuto = vtk.vtkPolyDataNormals()
        normAuto.ConsistencyOn()
        normAuto.AutoOrientNormalsOn()
        normAuto.SetInputConnection(clean.GetOutputPort())
        normAuto.Update()

        topMold = vtk.vtkPolyData()
        topMold.DeepCopy(normAuto.GetOutput())

        # Add fill back on to clipped bottom mold and clean
        append = vtk.vtkAppendPolyData()
        append.AddInputConnection(clipMid.GetClippedOutputPort())
        append.AddInputConnection(tri.GetOutputPort())
        append.Update()

        clean = vtk.vtkCleanPolyData()
        clean.SetInputConnection(append.GetOutputPort())
        clean.Update()

        # Extrude bottom part of mold down with filled clip plane
        extrudeDown = vtk.vtkLinearExtrusionFilter()
        extrudeDown.SetExtrusionTypeToVectorExtrusion()
        extrudeDown.SetVector(np.array(midClippingPlane.GetNormal()) * -1)
        extrudeDown.SetScaleFactor(40)
        extrudeDown.SetInputConnection(clean.GetOutputPort())
        extrudeDown.CappingOff()
        extrudeDown.Update()

        normAuto = vtk.vtkPolyDataNormals()
        normAuto.ConsistencyOn()
        normAuto.AutoOrientNormalsOn()
        normAuto.SetInputConnection(extrudeDown.GetOutputPort())
        normAuto.Update()

        append = vtk.vtkAppendPolyData()
        append.AddInputConnection(normAuto.GetOutputPort())
        append.AddInputConnection(clean.GetOutputPort())
        append.Update()

        clean = vtk.vtkCleanPolyData()
        clean.SetInputConnection(append.GetOutputPort())
        clean.Update()

        # Perform the bottom clipping at the specified depth
        clipBase = vtk.vtkClipPolyData()
        clipBase.SetClipFunction(baseClippingPlane)
        clipBase.SetInputConnection(clean.GetOutputPort())
        clipBase.Update()

        # Fill bottom clip plane
        cutter = vtk.vtkCutter()
        cutter.SetCutFunction(baseClippingPlane)
        cutter.SetInputConnection(clean.GetOutputPort())
        cutter.Update()

        loop = vtk.vtkContourLoopExtraction()
        loop.SetNormal(baseClippingPlane.GetNormal())
        loop.SetLoopClosureToAll()
        loop.SetInputConnection(cutter.GetOutputPort())
        loop.Update()

        tri = vtk.vtkTriangleFilter()
        tri.SetInputConnection(loop.GetOutputPort())
        tri.Update()

        # Flip bottom surface so normal points out
        reverse = vtk.vtkReverseSense()
        reverse.SetInputConnection(tri.GetOutputPort())
        reverse.Update()

        appendBottom = vtk.vtkAppendPolyData()
        appendBottom.AddInputConnection(clipBase.GetOutputPort())
        appendBottom.AddInputConnection(reverse.GetOutputPort())
        appendBottom.Update()

        clean = vtk.vtkCleanPolyData()
        clean.SetInputConnection(appendBottom.GetOutputPort())
        clean.Update()

        normAuto = vtk.vtkPolyDataNormals()
        normAuto.ConsistencyOn()
        normAuto.AutoOrientNormalsOn()
        normAuto.SetInputConnection(clean.GetOutputPort())
        normAuto.Update()

        bottomMold = vtk.vtkPolyData()
        bottomMold.DeepCopy(normAuto.GetOutput())

        return topMold, bottomMold

    def addOrUpdateModel(self, model, name, tformId = None):
        node = slicer.util.getFirstNodeByName(name)
        if not node:
            node = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLModelNode', name)
            node.CreateDefaultDisplayNodes()

        node.SetAndObservePolyData(model)

        if tformId:
            node.SetAndObserveTransformNodeID(tformId)





class MVSegmenterTest(ScriptedLoadableModuleTest):
    """
    This is the test case for your scripted module.
    Uses ScriptedLoadableModuleTest base class, available at:
    https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
    """

    def setUp(self):
        """ Do whatever is needed to reset the state - typically a scene clear will be enough.
        """
        slicer.mrmlScene.Clear(0)

    def runTest(self):
        """Run as few or as many tests as needed here.
        """
        self.setUp()
        self.test_MVSegmenter1()

    def test_MVSegmenter1(self):
        """ Ideally you should have several levels of tests.  At the lowest level
        tests should exercise the functionality of the logic with different inputs
        (both valid and invalid).  At higher levels your tests should emulate the
        way the user would interact with your code and confirm that it still works
        the way you intended.
        One of the most important features of the tests is that it should alert other
        developers when their changes will have an impact on the behavior of your
        module.  For example, if a developer removes a feature that you depend on,
        your test should break so they know that the feature is needed.
        """

        self.delayDisplay("No tests are implemented")
