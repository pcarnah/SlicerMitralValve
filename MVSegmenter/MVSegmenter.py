import os
import sys
import unittest
import importlib
from pathlib import Path
import vtk, qt, ctk, slicer
import SimpleITK as sitk
import sitkUtils
from slicer.ScriptedLoadableModule import *
import HeartValveLib
import numpy as np
import math
import logging
import platform


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
        self.papillaryMarkupsNode = None
        self.papillaryMarkupNodeObserver = None

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

        # Add vertical spacer
        self.layout.addSpacing(vSpace)

        self.runDeepMVButton = qt.QPushButton("Run DeepMV")
        self.runDeepMVButton.toolTip = "Run DeepMV Segmentation"
        self.runDeepMVButton.enabled = False
        self.layout.addWidget(self.runDeepMVButton)

        # Add vertical spacer
        self.layout.addSpacing(vSpace)

        # Semi-Automatic Segmentation
        semiAutoCollapsibleButton = ctk.ctkCollapsibleButton()
        semiAutoCollapsibleButton.text = "Semi-Automatic Segmentation"
        semiAutoCollapsibleButton.collapsed = True
        self.layout.addWidget(semiAutoCollapsibleButton)
        semiAutoFormLayout = qt.QFormLayout(semiAutoCollapsibleButton)

        #
        #  First Phase Segmentation
        #
        firstPassCollapsibleButton = ctk.ctkCollapsibleButton()
        firstPassCollapsibleButton.text = "Blood Pool Segmentation"
        semiAutoFormLayout.addWidget(firstPassCollapsibleButton)

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
        semiAutoFormLayout.addRow(secondPassCollapsibleButton)

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
        #  Export Inner Surface Model
        #
        generateSurfaceMoldCollapsibleButton = ctk.ctkCollapsibleButton()
        generateSurfaceMoldCollapsibleButton.text = "Generate Mold"
        generateSurfaceMoldCollapsibleButton.collapsed = False
        self.layout.addWidget(generateSurfaceMoldCollapsibleButton)

        exportModelFormLayout = qt.QFormLayout(generateSurfaceMoldCollapsibleButton)

        self.baseDepthSlider = ctk.ctkSliderWidget()
        self.baseDepthSlider.singleStep = 0.1
        self.baseDepthSlider.pageStep = 1
        self.baseDepthSlider.minimum = -20
        self.baseDepthSlider.maximum = 0
        self.baseDepthSlider.value = -12.5
        self.baseDepthSlider.decimals = 1
        self.baseDepthSlider.tracking = False
        self.baseDepthSlider.setToolTip("Adjust the depth of the base surface for the mold")
        exportModelFormLayout.addRow("Base Clipping Depth", self.baseDepthSlider)

        # Export button
        self.generateMoldButton = qt.QPushButton("Generate Mold")
        self.generateMoldButton.toolTip = "Generate the inner surface mold in the segmentation node."
        self.generateMoldButton.enabled = False
        exportModelFormLayout.addRow(self.generateMoldButton)

        # Annulus Projection Offset
        self.annulusOffsetSlider = ctk.ctkSliderWidget()
        self.annulusOffsetSlider.singleStep = 0.05
        self.annulusOffsetSlider.pageStep = 1
        self.annulusOffsetSlider.minimum = -5
        self.annulusOffsetSlider.maximum = 1
        self.annulusOffsetSlider.value = -1
        self.annulusOffsetSlider.decimals = 2
        self.annulusOffsetSlider.tracking = False
        self.annulusOffsetSlider.setToolTip("Adjust the offset of the annulus projection")
        exportModelFormLayout.addRow("Annulus Projection Offset", self.annulusOffsetSlider)

        self.projectAnnulusButton = qt.QPushButton("Project Annulus")
        self.projectAnnulusButton.toolTip = "Project annulus onto mold"
        self.projectAnnulusButton.enabled = False
        exportModelFormLayout.addRow(self.projectAnnulusButton)

        self.subtractAnnulusButton = qt.QPushButton("Subtract Annulus")
        self.subtractAnnulusButton.toolTip = "Subtract annulus from mold"
        self.subtractAnnulusButton.enabled = False
        exportModelFormLayout.addRow(self.subtractAnnulusButton)

        papillaryHBox = qt.QHBoxLayout()
        # papillaryHBox.addStretch(20)

        self.addPapillary1Button = qt.QPushButton("Place Papillary Tips")
        self.addPapillary1Button.toolTip = "Place the papillary muscle tips"
        papillaryHBox.addWidget(self.addPapillary1Button)

        self.deleteLastPapillaryButton = qt.QPushButton("Delete Last Tip Markup")
        self.deleteLastPapillaryButton.toolTip = "Delete the last papillary muscle tip"
        papillaryHBox.addWidget(self.deleteLastPapillaryButton)

        self.deleteAllPapillarryButton = qt.QPushButton("Delete All Tip Markups")
        self.deleteAllPapillarryButton.toolTip = "Delete all papillary muscle tips"
        papillaryHBox.addWidget(self.deleteAllPapillarryButton)

        exportModelFormLayout.addRow("Papillary Placement", papillaryHBox)

        self.exportMoldButton = qt.QPushButton("Export Mold to Model")
        self.exportMoldButton.toolTip = "Export mold segmentation node to model."
        self.exportMoldButton.enabled = False
        exportModelFormLayout.addRow(self.exportMoldButton)

        # Add vertical spacer
        self.layout.addSpacing(vSpace)

        # connections
        self.runDeepMVButton.connect('clicked(bool)', self.onRunDeepMV)
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

        # self.markupsSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.onSelect)
        # self.generateSurfaceMarkups.connect('clicked(bool)', self.onGenerateSurfaceMarkups)

        self.generateMoldButton.connect('clicked(bool)', self.onGenerateModelButton)
        self.projectAnnulusButton.connect('clicked(bool)', self.onProjectAnnulusButton)
        self.subtractAnnulusButton.connect('clicked(bool)', self.onSubtractAnnulusButton)
        self.addPapillary1Button.connect('clicked(bool)', self.onAddPapillaryButton)
        self.deleteLastPapillaryButton.connect('clicked(bool)', self.onDeleteLastPapillaryButton)
        self.deleteAllPapillarryButton.connect('clicked(bool)', self.onDeleteAllPapillaryButton)
        self.exportMoldButton.connect('clicked(bool)', self.onExportModelButton)

        # Add vertical spacer
        self.layout.addStretch(1)

        # Refresh Apply button state
        self.onSelect()

    def cleanup(self):
        self.papillaryMarkupsNode.RemoveObserver(self.papillaryMarkupNodeObserver)
        self.papillaryMarkupsNode = None
        self.papillaryMarkupNodeObserver = None

    def setupPapillaryMarkups(self):
        if self.papillaryMarkupsNode:
            if self.papillaryMarkupsNode.GetScene():
                return
            else:
                # Cleanup if node removed from scene or scene cleared
                self.papillaryMarkupsNode.RemoveObserver(self.papillaryMarkupNodeObserver)
                self.papillaryMarkupsNode = None
                self.papillaryMarkupNodeObserver = None

        self.papillaryMarkupsNode = slicer.util.getFirstNodeByClassByName('vtkMRMLMarkupsFiducialNode',
                                                                          'Papillary Tips Markup')
        if not self.papillaryMarkupsNode:
            self.papillaryMarkupsNode = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLMarkupsFiducialNode',
                                                                           'Papillary Tips Markup')

        if not self.papillaryMarkupsNode.GetDisplayNode():
            self.papillaryMarkupsNode.CreateDefaultDisplayNodes()

        self.papillaryMarkupsNode.GetDisplayNode().SetTextScale(2.0)
        self.papillaryMarkupsNode.GetDisplayNode().SetGlyphSize(3.0)
        self.papillaryMarkupsNode.GetDisplayNode().SetUseGlyphScale(False)
        self.papillaryMarkupsNode.GetDisplayNode().SetOpacity(0.6)

        self.papillaryMarkupNodeObserver = \
            self.papillaryMarkupsNode.AddObserver(slicer.vtkMRMLMarkupsNode.PointModifiedEvent,
                                                  self.onPapillaryMarkupNodeModified)

    def onSelect(self):
        self.setupPapillaryMarkups()

        self.runDeepMVButton.enabled = self.heartValveSelector.currentNode() and self.inputSelector.currentNode() and self.outputSegmentationSelector.currentNode()
        self.initBPButton.enabled = self.heartValveSelector.currentNode() and self.inputSelector.currentNode() and self.outputSegmentationSelector.currentNode()
        # self.generateSurfaceMarkups.enabled = self.heartValveSelector.currentNode() and self.outputSegmentationSelector.currentNode() and self.markupsSelector.currentNode()
        self.generateMoldButton.enabled = self.heartValveSelector.currentNode() \
                                          and self.outputSegmentationSelector.currentNode() \
                                          and self.inputSelector.currentNode() \
                                          and self.outputSegmentationSelector.currentNode().GetSegmentation().GetSegment('Leaflet Segmentation')

        self.projectAnnulusButton.enabled = self.heartValveSelector.currentNode() and self.outputSegmentationSelector.currentNode() \
                                        and self.outputSegmentationSelector.currentNode().GetSegmentation().GetSegment('Final_Mold')

        numberOfPoints = self.papillaryMarkupsNode.GetNumberOfDefinedControlPoints()
        self.exportMoldButton.enabled = self.projectAnnulusButton.enabled and self.outputSegmentationSelector.currentNode().GetSegmentation().GetSegment('Projected_Annulus') \
                                        and numberOfPoints >= 2
        self.subtractAnnulusButton.enabled = self.projectAnnulusButton.enabled and self.outputSegmentationSelector.currentNode().GetSegmentation().GetSegment('Projected_Annulus')

        self.addPapillary1Button.enabled = self.generateMoldButton.enabled and numberOfPoints < 2
        self.deleteLastPapillaryButton.enabled = numberOfPoints > 0
        self.deleteAllPapillarryButton.enabled = numberOfPoints > 0

    def onRunDeepMV(self):
        try:
            # This can be a long operation - indicate it to the user
            qt.QApplication.setOverrideCursor(qt.Qt.WaitCursor)

            self.logic.runDeepMV(self.heartValveSelector.currentNode(), self.inputSelector.currentNode(), self.outputSegmentationSelector.currentNode())
            self.onSelect()
        finally:
            qt.QApplication.restoreOverrideCursor()

    def onInitBPButton(self):
        try:
            # This can be a long operation - indicate it to the user
            qt.QApplication.setOverrideCursor(qt.Qt.WaitCursor)

            self.logic.initBPSeg(self.inputSelector.currentNode(), self.heartValveSelector.currentNode(),
                                 self.outputSegmentationSelector.currentNode())
            self.incrementFirstButton50.enabled = True
            self.incrementFirstButton100.enabled = True
            self.incrementFirstButton500.enabled = True
            self.initLeafletButton.enabled = True
            self.onSelect()
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
            self.onSelect()

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
        # Will disable button when undo stack is empty
        self.undoButtonBP.enabled = self.logic.undoBPIteration(self.outputSegmentationSelector.currentNode())

        # Enable redo button on undo
        self.redoButtonBP.enabled = True


    def onRedoButtonBP(self):
        # Will disable button when redo stack is empty
        self.redoButtonBP.enabled = self.logic.redoBPIteration(self.outputSegmentationSelector.currentNode())

        # Enable undo button on redo
        self.undoButtonBP.enabled = True

    def onUndoButtonLeaflet(self):
        # Will disable button when undo stack is empty
        self.undoButtonLeaflet.enabled = self.logic.undoLeafletIteration(self.outputSegmentationSelector.currentNode())

        # Enable redo button on undo
        self.redoButtonLeaflet.enabled = True

    def onRedoButtonLeaflet(self):
        # Will disable button when redo stack is empty
        self.redoButtonLeaflet.enabled = self.logic.redoLeafletIteration(self.outputSegmentationSelector.currentNode())

        # Enable undo button on redo
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

    def onProjectAnnulusButton(self):
        try:
            # This can be a long operation - indicate it to the user
            qt.QApplication.setOverrideCursor(qt.Qt.WaitCursor)

            self.logic.projectAnnulus(self.outputSegmentationSelector.currentNode(),
                                           self.heartValveSelector.currentNode(),
                                           float(self.annulusOffsetSlider.value))
            self.onSelect()

        finally:
            qt.QApplication.restoreOverrideCursor()

    def onSubtractAnnulusButton(self):
        try:
            # This can be a long operation - indicate it to the user
            qt.QApplication.setOverrideCursor(qt.Qt.WaitCursor)

            self.logic.subtractAnnulusSegmentation(self.outputSegmentationSelector.currentNode(),
                                                   self.inputSelector.currentNode())
            self.onSelect()

        finally:
            qt.QApplication.restoreOverrideCursor()


    def onAddPapillaryButton(self):
        self.papillaryMarkupsNode.SetAndObserveTransformNodeID(self.outputSegmentationSelector.currentNode().GetTransformNodeID())
        slicer.modules.markups.logic().SetActiveListID(self.papillaryMarkupsNode)
        slicer.app.applicationLogic().GetInteractionNode().SetPlaceModePersistence(1)
        slicer.app.applicationLogic().GetInteractionNode().SetCurrentInteractionMode(1)


    def onDeleteLastPapillaryButton(self):
        numberOfPoints = self.papillaryMarkupsNode.GetNumberOfDefinedControlPoints()
        if numberOfPoints > 0:
            self.papillaryMarkupsNode.RemoveMarkup(numberOfPoints - 1)
        self.onSelect()

    def onDeleteAllPapillaryButton(self):
        while self.papillaryMarkupsNode.GetNumberOfDefinedControlPoints() > 0:
            self.papillaryMarkupsNode.RemoveMarkup(self.papillaryMarkupsNode.GetNumberOfDefinedControlPoints() - 1)
        self.onSelect()

    def onPapillaryMarkupNodeModified(self, unusedArg1=None, unusedArg2=None, unusedArg3=None):
        numberOfPoints = self.papillaryMarkupsNode.GetNumberOfDefinedControlPoints()

        if numberOfPoints >= 2:
            slicer.app.applicationLogic().GetInteractionNode().SetCurrentInteractionMode(2)

        self.onSelect()

    def onExportModelButton(self):
        try:
            # This can be a long operation - indicate it to the user
            qt.QApplication.setOverrideCursor(qt.Qt.WaitCursor)

            self.logic.exportSurfaceMold(self.outputSegmentationSelector.currentNode(), self.papillaryMarkupsNode)

            self.onSelect()

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
        self._speedImgRefNode = None
        self._leafletLevelSet = None
        self._bpLevelSet = None

        self._undoBPLevelSetStack = []
        self._redoBPLevelSetStack = []
        self._undoLeafletLevelSetStack = []
        self._redoLeafletLevelSetStack = []

        self.moldBasePlate = None


    def initBPSeg(self, inputVolume, heartValveNode, outputSeg):
        """
        Initialize the blood pool segmentation
        :param inputVolume: The reference image volume to performa segmentation on
        :param heartValveNode: The SlicerHeart MRML node containing annulus definition
        :param outputSeg: Segmentation node to save output to
        :return None
        """

        if not inputVolume or not heartValveNode or not outputSeg:
            logging.debug("initBPSeg failed: Missing parameter")
            return

        valveModel = HeartValveLib.getValveModel(heartValveNode)
        if valveModel.getAnnulusContourMarkupNode().GetNumberOfFiducials() == 0:
            logging.debug("initBPSeg failed: Annulus contour not defined")
            return

        if valveModel.getProbeToRasTransformNode():
            outputSeg.SetAndObserveTransformNodeID(valveModel.getProbeToRasTransformNode().GetID())

        if not outputSeg.GetNodeReference(outputSeg.GetReferenceImageGeometryReferenceRole()):
            outputSeg.SetReferenceImageGeometryParameterFromVolumeNode(inputVolume)

        self._speedImgRefNode = inputVolume
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

        self.updateBPLevelSet(out_mask)

        threshold = sitk.BinaryThresholdImageFilter()
        threshold.SetInsideValue(1)
        threshold.SetLowerThreshold(-1000.0)
        threshold.SetOutsideValue(0)
        threshold.SetUpperThreshold(0.0)
        out_mask = threshold.Execute(out_mask)

        self.pushITKImageToSegmentation(out_mask, outputSeg, 'BP Segmentation')

    def iterateFirstPass(self, nIter, outputSeg):
        """
        Iterates the blood pool segmentation by nIter amount
        :param nIter: Number of iterations to run the active contour algorithm for
        :param outputSeg: Segmentation node to save output to
        :return: None
        """
        self.updateBPLevelSetFromSegmentation(outputSeg)

        geodesicActiveContour = sitk.GeodesicActiveContourLevelSetImageFilter()
        geodesicActiveContour.SetCurvatureScaling(1.2)
        geodesicActiveContour.SetAdvectionScaling(1.0)
        geodesicActiveContour.SetPropagationScaling(0.9)
        geodesicActiveContour.SetMaximumRMSError(0.00001)
        geodesicActiveContour.SetNumberOfIterations(nIter)

        out_mask = geodesicActiveContour.Execute(self._bpLevelSet, self._speedImg)

        self.updateBPLevelSet(out_mask)

        threshold = sitk.BinaryThresholdImageFilter()
        threshold.SetInsideValue(1)
        threshold.SetLowerThreshold(-1000.0)
        threshold.SetOutsideValue(0)
        threshold.SetUpperThreshold(0.0)
        out_mask = threshold.Execute(out_mask)

        self.pushITKImageToSegmentation(out_mask, outputSeg, 'BP Segmentation')

    def initLeafletSeg(self, outputSeg):
        """
        Initializes the leaflet segmentation based on the boundary of the blood pool segmentation
        :param outputSeg: Segmentation node to save output to
        :return: None
        """
        self.updateBPLevelSetFromSegmentation(outputSeg)

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

        self.updateLeafletLevelSet(out_mask)

        threshold = sitk.BinaryThresholdImageFilter()
        threshold.SetInsideValue(1)
        threshold.SetLowerThreshold(-1000.0)
        threshold.SetOutsideValue(0)
        threshold.SetUpperThreshold(0.0)
        out_mask = threshold.Execute(out_mask)

        self.pushITKImageToSegmentation(out_mask, outputSeg, 'Leaflet Segmentation')

        return out_mask

    def iterateSecondPass(self, nIter, outputSeg):
        """
        Iterates the leaflet segmentation by nIter amount
        :param nIter: Number of iterations to run the active contour algorithm for
        :param outputSeg: Segmentation node to save the output to
        :return: None
        """
        self.updateLeafletLevelSetFromSegmentation(outputSeg)

        geodesicActiveContour2 = sitk.GeodesicActiveContourLevelSetImageFilter()
        geodesicActiveContour2.SetCurvatureScaling(0.9)
        geodesicActiveContour2.SetAdvectionScaling(0.1)
        geodesicActiveContour2.SetPropagationScaling(-0.4)
        geodesicActiveContour2.SetMaximumRMSError(0.0001)
        geodesicActiveContour2.SetNumberOfIterations(nIter)
        out_mask = geodesicActiveContour2.Execute(self._leafletLevelSet, self._speedImg)

        self.updateLeafletLevelSet(out_mask)

        threshold = sitk.BinaryThresholdImageFilter()
        threshold.SetInsideValue(1)
        threshold.SetLowerThreshold(-1000.0)
        threshold.SetOutsideValue(0)
        threshold.SetUpperThreshold(0.0)
        out_mask = threshold.Execute(out_mask)

        self.pushITKImageToSegmentation(out_mask, outputSeg, 'Leaflet Segmentation')

        return out_mask

    def updateBPLevelSet(self, levelSet):
        """
        Update the blood pool level set instance variable. Maintains the undo stack.
        :param levelSet: New level set
        :return: None
        """
        if self._bpLevelSet:
            # Only update if there has been a change to preserve undo pool
            if not self.levelSetsEqual(levelSet, self._bpLevelSet):
                self._undoBPLevelSetStack.append(self._bpLevelSet)
                self._bpLevelSet = levelSet
        else:
            self._bpLevelSet = levelSet

    def undoBPIteration(self, outputSeg):
        """
        Performs the undo operation on the blood pool level set
        :param outputSeg: Segmentation node used to save output
        :return: True if the stack has additional items, False if it is empty
        """
        if not self._undoBPLevelSetStack:
            logging.debug("undoBPIteration failed: stack empty")
            return False

        self._redoBPLevelSetStack.append(self._bpLevelSet)
        self._bpLevelSet = self._undoBPLevelSetStack.pop()

        threshold = sitk.BinaryThresholdImageFilter()
        threshold.SetInsideValue(1)
        threshold.SetLowerThreshold(-1000.0)
        threshold.SetOutsideValue(0)
        threshold.SetUpperThreshold(0.0)

        self.pushITKImageToSegmentation(threshold.Execute(self._bpLevelSet), outputSeg, 'BP Segmentation')

        # Return true if stack is not empty
        if self._undoBPLevelSetStack:
            return True
        else:
            return False

    def redoBPIteration(self, outputSeg):
        """
        Performs the redo operation on the blood pool level set
        :param outputSeg: Segmentation node used to save output
        :return: True if the stack has additional items, False if it is empty
        """
        if not self._redoBPLevelSetStack:
            logging.debug("redoBPIteration failed: stack empty")
            return False

        self._undoBPLevelSetStack.append(self._bpLevelSet)
        self._bpLevelSet = self._redoBPLevelSetStack.pop()

        threshold = sitk.BinaryThresholdImageFilter()
        threshold.SetInsideValue(1)
        threshold.SetLowerThreshold(-1000.0)
        threshold.SetOutsideValue(0)
        threshold.SetUpperThreshold(0.0)

        self.pushITKImageToSegmentation(threshold.Execute(self._bpLevelSet), outputSeg, 'BP Segmentation')

        # Return true if stack is not empty
        if self._redoBPLevelSetStack:
            return True
        else:
            return False

    def updateLeafletLevelSet(self, levelSet):
        """
        Update the leaflet level set instance variable. Maintains the undo stack.
        :param levelSet: New level set
        :return: None
        """
        if self._leafletLevelSet:
            # Only update if there has been a change to preserve undo pool
            if not self.levelSetsEqual(levelSet, self._leafletLevelSet):
                self._undoLeafletLevelSetStack.append(self._leafletLevelSet)
                self._leafletLevelSet = levelSet
        else:
            self._leafletLevelSet = levelSet

    def undoLeafletIteration(self, outputSeg):
        """
        Performs the undo operation on the leaflet level set
        :param outputSeg: Segmentation node used to save output
        :return: True if the stack has additional items, False if it is empty
        """
        if not self._undoLeafletLevelSetStack:
            logging.debug("undoLeafletIteration failed: stack empty")
            return False

        self._redoLeafletLevelSetStack.append(self._leafletLevelSet)
        self._leafletLevelSet = self._undoLeafletLevelSetStack.pop()

        threshold = sitk.BinaryThresholdImageFilter()
        threshold.SetInsideValue(1)
        threshold.SetLowerThreshold(-1000.0)
        threshold.SetOutsideValue(0)
        threshold.SetUpperThreshold(0.0)

        self.pushITKImageToSegmentation(threshold.Execute(self._leafletLevelSet), outputSeg, 'Leaflet Segmentation')

        # Return true if stack is not empty
        if self._undoLeafletLevelSetStack:
            return True
        else:
            return False


    def redoLeafletIteration(self, outputSeg):
        """
        Performs the redo operation on the leaflet level set
        :param outputSeg: Segmentation node used to save output
        :return: True if the stack has additional items, False if it is empty
        """

        if not self._redoLeafletLevelSetStack:
            logging.debug("redoLeafletIteration failed: stack empty")
            return False

        self._undoLeafletLevelSetStack.append(self._leafletLevelSet)
        self._leafletLevelSet = self._redoLeafletLevelSetStack.pop()

        threshold = sitk.BinaryThresholdImageFilter()
        threshold.SetInsideValue(1)
        threshold.SetLowerThreshold(-1000.0)
        threshold.SetOutsideValue(0)
        threshold.SetUpperThreshold(0.0)

        self.pushITKImageToSegmentation(threshold.Execute(self._leafletLevelSet), outputSeg, 'Leaflet Segmentation')

        # Return true if stack is not empty
        if self._redoLeafletLevelSetStack:
            return True
        else:
            return False

    def updateBPLevelSetFromSegmentation(self, segNode, segmentId = 'BP Segmentation'):
        """
        Retrieves the new level set from the Segmentation node and update the blood pool level set for subsequent operations.
        :param segNode: Segmentation node used to store output
        :param segmentId: The segment ID to access
        :return: None
        """
        # Get new binary mask from segmentation node
        mask = self.pullITKImageFromSegmentation(segNode, segmentId, self._speedImgRefNode)

        # Get level set from mask
        signedDis = sitk.SignedDanielssonDistanceMapImageFilter()
        levelSet = signedDis.Execute(mask)

        self.updateBPLevelSet(levelSet)

    def updateLeafletLevelSetFromSegmentation(self, segNode, segmentId = 'Leaflet Segmentation'):
        """
        Retrieves the new level set from the Segmentation node and update the leaflet level set for subsequent operations.
        :param segNode: Segmentation node used to store output
        :param segmentId: The segment ID to access
        :return: None
        """
        # Get new binary mask from segmentation node
        mask = self.pullITKImageFromSegmentation(segNode, segmentId, self._speedImgRefNode)

        # Get level set from mask
        signedDis = sitk.SignedDanielssonDistanceMapImageFilter()
        levelSet = signedDis.Execute(mask)

        self.updateLeafletLevelSet(levelSet)

    def levelSetsEqual(self, lvlset1, lvlset2):
        """
        Compares if 2 level sets cover the same area to determine equality.
        :param lvlset1: First level set for comparison
        :param lvlset2: Second level set for comparison
        :return: True if equal, False otherwise
        """

        # Convert level sets to label maps
        threshold = sitk.BinaryThresholdImageFilter()
        threshold.SetInsideValue(1)
        threshold.SetLowerThreshold(-1000.0)
        threshold.SetOutsideValue(0)
        threshold.SetUpperThreshold(0.0)
        im = threshold.Execute(lvlset1)
        im2 = threshold.Execute(lvlset2)

        # Check if binary label maps match
        stat = sitk.StatisticsImageFilter()
        stat.Execute(sitk.NotEqual(im, im2))
        max = stat.GetMaximum()

        # Level sets are equal if max of pixelwise not is 0
        return max == 0

    def pushITKImageToSegmentation(self, img, segmentationNode, segmentId='Leaflet Segmentation'):
        """
        Pushes an itk image to a Segmentation MRML node
        :param img: The itk image
        :param segmentationNode: The output segmentation ndoe
        :param segmentId: The segment ID to store in
        :return: None
        """
        if segmentationNode.GetSegmentation().GetSegmentIndex(segmentId) == -1:
            segmentationNode.GetSegmentation().AddEmptySegment(segmentId)

        if segmentationNode.GetSegmentation().GetConversionParameter('Decimation factor') != '0.5':
            segmentationNode.GetSegmentation().SetConversionParameter('Decimation factor', '0.5')

        if segmentationNode.GetSegmentation().GetConversionParameter('Smoothing factor') != '0.5':
            segmentationNode.GetSegmentation().SetConversionParameter('Smoothing factor', '0.5')

        slicer.modules.SegmentEditorWidget.editor.mrmlSegmentEditorNode().SetOverwriteMode(
            slicer.modules.SegmentEditorWidget.editor.mrmlSegmentEditorNode().OverwriteNone)


        # Create temporary label map node
        tempNode = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLLabelMapVolumeNode', 'temp_labelmap')
        tempNode.SetAndObserveTransformNodeID(segmentationNode.GetTransformNodeID())
        sitkUtils.PushVolumeToSlicer(img, tempNode)

        segmentationIds = vtk.vtkStringArray()
        segmentationIds.InsertNextValue(segmentId)

        slicer.modules.segmentations.logic().ImportLabelmapToSegmentationNode(tempNode, segmentationNode,
                                                                              segmentationIds)

        segmentationNode.RemoveClosedSurfaceRepresentation()
        segmentationNode.CreateClosedSurfaceRepresentation()

        slicer.mrmlScene.RemoveNode(tempNode)

    def pullITKImageFromSegmentation(self, segmentationNode, segmentId, refNode = None):
        """
        Retrieves an itk image from a Segmentation MRML node
        :param segmentationNode: The output segmentation node
        :param segmentId: The segment ID to access
        :return: The itk image
        """
        if not segmentationNode or not segmentId:
            logging.debug("pullITKImageFromSegmentation failed: Missing parameter")
            return

        if segmentationNode.GetSegmentation().GetSegmentIndex(segmentId) == -1:
            logging.debug('pullITKImageFromSegmentation failed: Segment not found - ' + segmentId)
            return

        # Create temporary label map node to get itk image from slicer volume
        tempNode = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLLabelMapVolumeNode', 'temp_labelmap')
        tempNode.SetAndObserveTransformNodeID(segmentationNode.GetTransformNodeID())

        segmentationIds = vtk.vtkStringArray()
        segmentationIds.InsertNextValue(segmentId)

        # Export segmentation labelmap to temporary node
        if refNode:
            slicer.modules.segmentations.logic().ExportSegmentsToLabelmapNode(segmentationNode, segmentationIds,
                                                                              tempNode, refNode)
        else:
            slicer.modules.segmentations.logic().ExportSegmentsToLabelmapNode(segmentationNode, segmentationIds,
                                                                              tempNode)
        # Get the itk image
        itkImg = sitkUtils.PullVolumeFromSlicer(tempNode)

        # Remove the temporary node
        slicer.mrmlScene.RemoveNode(tempNode)

        return itkImg


    def rasToIJK(self, point, volume):
        """
        Convert a RAS point to an IJK point for a reference volume. Uses the volume node's RAStoIJK matrix.
        :param point: The RAS input point.
        :param volume: The reference volume
        :return: The output point as a list of 3 points
        """
        matrix = vtk.vtkMatrix4x4()
        volume.GetRASToIJKMatrix(matrix)

        homogeneousPoint = [point[0], point[1], point[2], 1]
        outPoint = matrix.MultiplyPoint(homogeneousPoint)

        return outPoint[0:3]

    def generateSurfaceMarkups(self, segNode, heartValveNode, markupsNode):
        """
        Generate evenly spaced markups on the surface of the model.
        :param segNode: Segmentation node storing model
        :param heartValveNode: SlicerHeart HeartValve MRML node contatining annulus definition
        :param markupsNode: Ouput markups node
        :return: True if success, False otherwise
        """
        if not segNode or not heartValveNode or not markupsNode:
            logging.debug("generateSurfaceMarkups failed: Missing parameter")
            return False

        valveModel = HeartValveLib.getValveModel(heartValveNode)
        if valveModel.getAnnulusContourMarkupNode().GetNumberOfFiducials() == 0:
            logging.debug("generateSurfaceMarkups failed: Annulus contour not defined")
            return False

        leafletModel = segNode.GetClosedSurfaceInternalRepresentation('Leaflet Segmentation')
        if leafletModel is None:
            logging.debug("generateSurfaceMarkups failed: Missing segmentation")
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
        """
        Method to push a vtkPolyData model to a segmentation node
        :param segNode: The Segmentation node
        :param model: vtkPolyData model
        :param name: Segment name to push to
        :return: None
        """
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
        """
        Generate the complete surface mold from the segmentation. Clips the bottom of the mold to a specified depth.
        Will output mold to Segmentation node.
        :param segNode: The Segmentation node
        :param heartValveNode: SlicerHeart HeartValve MRML node contatining annulus definition
        :param depth: Clipping depth
        :param volume: Reference volume
        :return: None
        """
        # Check that parameters exist
        if not segNode or not heartValveNode or not depth or not volume:
            logging.debug("generateSurfaceMold failed: Missing parameter")
            return None

        valveModel = HeartValveLib.getValveModel(heartValveNode)
        segNode.CreateClosedSurfaceRepresentation()

        # Get the proximal surface of the valve
        extractedSurface = self.extractInnerSurfaceModel(segNode, valveModel)
        if not extractedSurface:
            return None

        contourPlane = valveModel.getAnnulusContourPlane()

        # Create clipped leaflet mold across middle
        baseClippingPlane = vtk.vtkPlane()
        baseClippingPlane.SetNormal(contourPlane[1])
        baseClippingPlane.SetOrigin(contourPlane[0] + contourPlane[1] * depth)

        midClippingPlane = vtk.vtkPlane()
        midClippingPlane.SetNormal(contourPlane[1])
        midClippingPlane.SetOrigin(contourPlane[0])
        # midClippingPlane.SetOrigin(contourPlane[0] + (contourPlane[1] * 2))

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

    def projectAnnulus(self, segNode, heartValveNode, offset = 0):
        """
        Project the annulus onto the mold model using the optional offset
        :param segNode: Segmentation node
        :param heartValveNode: Heart valve node containing annulus
        :param offset: Optional offset value for projection
        :return: None
        """
        # Check that parameters exist
        if not segNode or not heartValveNode:
            logging.debug("projectAnnulus failed: Missing parameter")
            return None

        valveModel = HeartValveLib.getValveModel(heartValveNode)

        # Get segmentation closed surface representations
        segMold = segNode.GetClosedSurfaceInternalRepresentation('Final_Mold')
        if not segMold:
            logging.debug("projectAnnulus: Missing mold segmentation")
            return None


        # Get the annulus projected onto the proximal surface
        projectedAnnulus, stiffener = self.generateProjectedAnnulus(segMold, valveModel, offset)
        if not projectedAnnulus:
            return None

        self.pushModelToSegmentation(segNode, projectedAnnulus, 'Projected_Annulus')
        self.pushModelToSegmentation(segNode, stiffener, 'Stiffener_Surface')

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
        segmentEditorWidget.setMasterVolumeNode(valveModel.getValveVolumeNode())
        segmentEditorWidget.setCurrentSegmentID('Stiffener_Surface')

        # Smoothing
        segmentEditorWidget.setActiveEffectByName("Smoothing")
        effect = segmentEditorWidget.activeEffect()
        effect.setParameter("SmoothingMethod", "GAUSSIAN")
        effect.setParameter("GaussianStandardDeviationMm", 1.0)
        effect.self().onApply()

        # Clean up
        segmentEditorWidget = None
        slicer.mrmlScene.RemoveNode(segmentEditorNode)

    def subtractAnnulusSegmentation(self, segNode, volume):
        """
        Performs binary subtraction in the segmentation node
        :param segNode: Segmentation node
        :param volume: Master volume
        :return: None
        """
        # Check that parameters exist
        if not segNode or not volume:
            logging.debug("subtractAnnulusSegmentation failed: Missing parameter")
            return None

        if not segNode.GetSegmentation().GetSegment('Projected_Annulus') or not segNode.GetSegmentation().GetSegment('Final_Mold'):
            logging.debug("subtractAnnulusSegmentation failed: Missing segment")

        # Create segment editor to get access to effects
        segmentEditorWidget = slicer.qMRMLSegmentEditorWidget()
        segmentEditorWidget.setMRMLScene(slicer.mrmlScene)
        segmentEditorNode = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLSegmentEditorNode")
        segmentEditorNode.SetOverwriteMode(segmentEditorNode.OverwriteNone)
        segmentEditorWidget.setMRMLSegmentEditorNode(segmentEditorNode)
        segmentEditorWidget.setSegmentationNode(segNode)
        segmentEditorWidget.setMasterVolumeNode(volume)
        segmentEditorWidget.setCurrentSegmentID('Final_Mold')

        # Subtraction
        segmentEditorWidget.setActiveEffectByName("Logical operators")
        effect = segmentEditorWidget.activeEffect()
        effect.setParameter("Operation", "SUBTRACT")
        effect.setParameter("Bypass masking", 1)
        effect.setParameter("ModifierSegmentID", 'Projected_Annulus')
        effect.self().onApply()

        # Smoothing
        segmentEditorWidget.setActiveEffectByName("Smoothing")
        effect = segmentEditorWidget.activeEffect()
        effect.setParameter("SmoothingMethod", "GAUSSIAN")
        effect.setParameter("GaussianStandardDeviationMm", 0.8)
        effect.self().onApply()

        # Clean up
        segmentEditorWidget = None
        slicer.mrmlScene.RemoveNode(segmentEditorNode)

    def exportSurfaceMold(self, segNode, papillaryMarkupsNode):
        """
        Export the surface mold from the Segmentation node to Models.
        :param segNode: Segmentation node containing mold
        :return: None
        """

        # Get segmentation closed surface representations
        segMold = segNode.GetClosedSurfaceInternalRepresentation('Final_Mold')
        annulusMold = segNode.GetClosedSurfaceInternalRepresentation('Projected_Annulus')
        stiffener = segNode.GetClosedSurfaceInternalRepresentation('Stiffener_Surface')
        if not segMold or not annulusMold or not stiffener:
            logging.debug("exportSurfaceMold failed: Missing mold segmentation")
            return

        # Set default polydata extension to stl
        defaultModelStorageNode = slicer.vtkMRMLModelStorageNode()
        defaultModelStorageNode.SetDefaultWriteFileExtension('stl')
        slicer.mrmlScene.AddDefaultNode(defaultModelStorageNode)

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


        # Export closed surface meshes to model nodes
        moldModel = vtk.vtkPolyData()
        moldModel.DeepCopy(normAuto.GetOutput())

        self.addOrUpdateModel(moldModel, 'Final_Mold_Model', segNode.GetTransformNodeID(), segNode.GetSegmentation().GetSegment('Final_Mold').GetColor())

        annulusModel = vtk.vtkPolyData()
        annulusModel.DeepCopy(annulusMold)
        self.addOrUpdateModel(annulusModel, 'Projected_Annulus_Model', segNode.GetTransformNodeID(), segNode.GetSegmentation().GetSegment('Projected_Annulus').GetColor())

        # Decimate stiffener ring to 98%
        decimate = vtk.vtkDecimatePro()
        decimate.SetTargetReduction(0.98)
        decimate.SetInputData(stiffener)
        decimate.Update()

        stiffenerModel = vtk.vtkPolyData()
        stiffenerModel.DeepCopy(decimate.GetOutput())
        self.addOrUpdateModel(stiffenerModel, 'Stiffener_Model', segNode.GetTransformNodeID(),
                              segNode.GetSegmentation().GetSegment('Stiffener_Surface').GetColor())

        papillaryModel = vtk.vtkPolyData()
        papAppend = vtk.vtkAppendPolyData()

        for i in range(papillaryMarkupsNode.GetNumberOfDefinedControlPoints()):
            p = np.array([0,0,0])
            papillaryMarkupsNode.GetNthControlPointPosition(i, p)

            sphereSource = vtk.vtkSphereSource()
            sphereSource.SetCenter(p)
            sphereSource.SetRadius(2)
            sphereSource.Update()

            papAppend.AddInputConnection(sphereSource.GetOutputPort())
        papAppend.Update()
        papillaryModel.DeepCopy(papAppend.GetOutput())

        self.addOrUpdateModel(papillaryModel, 'Papillary_Model', segNode.GetTransformNodeID())


    def extractInnerSurfaceModel(self, segNode, valveModel, segName = 'Leaflet Segmentation'):
        """
        Extracts the inner surface (proximal to image probe) from the segmentation. Uses surface normals of leaflet
        segmentation along with the annulus normal and blood pool definitions to determine points on inside of segmentaion.
        :param segNode: The segmentation node
        :param valveModel: SlicerHeart HeartValve MRML node containing annulus definition
        :param segName: Name of segment containing leaflet segmentation
        :return: vtkPolyData model of inner surface
        """
        if not segNode or not valveModel:
            logging.debug("extractInnerSurfaceModel failed: Missing parameter")
            return None

        if valveModel.getAnnulusContourMarkupNode().GetNumberOfFiducials() == 0:
            logging.debug("extractInnerSurfaceModel failed: Annulus contour not defined")
            return None

        annulusPlane = valveModel.getAnnulusContourPlane()

        segNode.CreateClosedSurfaceRepresentation()
        leafletModel = segNode.GetClosedSurfaceInternalRepresentation(segNode.GetSegmentation().GetSegmentIdBySegmentName(segName))
        if leafletModel is None:
            logging.debug("extractInnerSurfaceModel failed: Missing segmentation")
            return None

        # Decimate leaflet polydata for efficiency
        decimate = vtk.vtkDecimatePro()
        decimate.SetTargetReduction(0.6)
        decimate.PreserveTopologyOn()
        decimate.BoundaryVertexDeletionOff()

        decimate.SetInputData(leafletModel)
        decimate.Update()
        leafletModel = vtk.vtkPolyData()
        leafletModel.DeepCopy(decimate.GetOutput())

        clipped = vtk.vtkPolyData()
        clipped.DeepCopy(leafletModel)

        # Get Annulus minimum radius
        annulusPoints = valveModel.getAnnulusContourModelNode().GetPolyData().GetPoints()
        minDis = float('inf')
        for i in range(annulusPoints.GetNumberOfPoints()):
            p = annulusPoints.GetPoint(i)
            dis = np.linalg.norm(np.subtract(annulusPlane[0], p))
            minDis = min(minDis, dis)

        # Create tube along annulus plane normal to use for bottom half extraction
        lineSource = vtk.vtkLineSource()
        lineSource.SetPoint1(annulusPlane[0] + annulusPlane[1] * 20)
        lineSource.SetPoint2(annulusPlane[0] - annulusPlane[1] * 20)
        lineSource.Update()

        tubeFilter = vtk.vtkTubeFilter()
        tubeFilter.SetRadius(0.5 * minDis)
        tubeFilter.CappingOff()
        tubeFilter.SetNumberOfSides(50)
        tubeFilter.SetInputConnection(lineSource.GetOutputPort())
        tubeFilter.Update()

        # OBBTree for determining self intersection of rays
        obb = vtk.vtkOBBTree()
        obb.SetDataSet(leafletModel)
        obb.BuildLocator()

        locator = vtk.vtkPointLocator()
        locator.SetDataSet(tubeFilter.GetOutput())
        locator.BuildLocator()


        # Loop over remaining points and build scalar array using different techniques for top and bottom half
        a0 = np.zeros(3)
        p = annulusPlane[0] + annulusPlane[1] * 2
        points = vtk.vtkPoints()
        normals = clipped.GetPointData().GetNormals()
        tubeNormals = tubeFilter.GetOutput().GetPointData().GetNormals()
        scalars = vtk.vtkFloatArray()
        scalars.SetNumberOfValues(clipped.GetNumberOfPoints())
        for i in range(clipped.GetNumberOfPoints()):
            clipped.GetPoint(i, a0)
            if np.dot(annulusPlane[1], a0 - p) > 0:
                # Point is above annulus plane, set scalar based on self intersection (scalar value will determine clipping)
                r = obb.IntersectWithLine(a0, annulusPlane[0] + annulusPlane[1] * 5, points, None)

                # If not match try with point below annulus plane
                if points.GetNumberOfPoints() != 1:
                    r = obb.IntersectWithLine(a0, annulusPlane[0] - annulusPlane[1] * 5, points, None)

                # If only 1 intersection point, line does not cross through leaflet model as the line always intersects at a0
                if points.GetNumberOfPoints() == 1:
                    scalars.SetValue(i, 10)     # Set scalar to large value so point will be kept
                else:
                    scalars.SetValue(i, -10)    # Set scalar to small value so point will be discarded
            else:
                # Point is below annulus plane
                # Get the closest point on tube surface, find angle between 2 normals in radians
                closestPoint = locator.FindClosestPoint(a0)
                v = np.array(tubeNormals.GetTuple(closestPoint))
                n = np.array(normals.GetTuple(i))
                angle = math.acos(np.dot(n, v) / np.linalg.norm(n) / np.linalg.norm(v))
                scalars.SetValue(i, angle)      # Set scalar to angle

        # Scalars now angles in radians that we can threshold
        clipped.GetPointData().SetScalars(scalars)

        # Clip based on scalar values (keep scalars bigger than value)
        clip2 = vtk.vtkClipPolyData()
        clip2.GenerateClipScalarsOff()
        clip2.SetValue(1.5)  # 86 degrees threshold in radians
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


    def generateProjectedAnnulus(self, extractedLeaflet, valveModel, offset = 0):
        """
        Projects the annulus definition inwards onto the surface model for mold
        :param extractedLeaflet: vtkPolyData surface model
        :param valveModel: SlicerHeart HeartValve MRML node containing annulus definition
        :return: vtkPolyData model of projected annulus as a tube
        """
        if not extractedLeaflet or not valveModel:
            logging.debug("generateProjectedAnnulus failed: Missing parameter")
            return None

        if valveModel.getAnnulusContourMarkupNode().GetNumberOfFiducials() == 0:
            logging.debug("generateProjectedAnnulus failed: Annulus contour not defined")
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
        normals = vtk.vtkFloatArray()
        normals.SetNumberOfComponents(3)
        projPoints = vtk.vtkPoints()

        # Project annulus inwards towards center
        center = contourPlane[0] + offset * 5 * contourPlane[1]
        for i in range(annulusMarkups.GetNumberOfFiducials()):
            annulusMarkups.GetNthFiducialPosition(i, pos)
            # Get point above annulus to project towards
            stiffenerPos = pos + contourPlane[1] * 12
            r = obb.IntersectWithLine(pos + pos - center, center, points, None)
            if r != 0:
                projPoints.InsertNextPoint(points.GetPoint(0))
                normals.InsertNextTuple3(*((stiffenerPos - center) / np.linalg.norm(stiffenerPos - center)))

        projPoints.InsertNextPoint(projPoints.GetPoint(0))
        normals.InsertNextTuple3(*normals.GetTuple3(0))

        lines = vtk.vtkCellArray()
        lines.InsertNextCell(projPoints.GetNumberOfPoints())
        for i in range(projPoints.GetNumberOfPoints()):
            lines.InsertCellPoint(i)

        # Create spline polydata
        projContour = vtk.vtkPolyData()
        projContour.SetPoints(projPoints)
        projContour.SetLines(lines)
        projContour.GetPointData().SetNormals(normals)

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

        # Generate stiffener surface from mold outwards
        ext = vtk.vtkLinearExtrusionFilter()
        ext.SetExtrusionTypeToNormalExtrusion()
        ext.SetInputConnection(strip.GetOutputPort())
        ext.SetScaleFactor(100)
        ext.CappingOn()
        ext.Update()

        # Clean up stiffener surface
        norm = vtk.vtkPolyDataNormals()
        norm.ConsistencyOn()
        norm.FlipNormalsOn()
        norm.SplittingOn()
        norm.SetInputConnection(ext.GetOutputPort())
        norm.Update()

        ext2 = vtk.vtkLinearExtrusionFilter()
        ext2.SetExtrusionTypeToNormalExtrusion()
        ext2.SetInputConnection(norm.GetOutputPort())
        ext2.SetScaleFactor(1.75)
        ext2.CappingOn()
        ext2.Update()

        norm = vtk.vtkPolyDataNormals()
        norm.ConsistencyOn()
        norm.SplittingOn()
        norm.AutoOrientNormalsOn()
        norm.SetInputConnection(ext2.GetOutputPort())
        norm.Update()

        cleanStiffener = vtk.vtkCleanPolyData()
        cleanStiffener.SetInputConnection(norm.GetOutputPort())
        cleanStiffener.Update()

        # Push fitted annulus onto segmentation node
        annulusFittedModel = vtk.vtkPolyData()
        annulusFittedModel.DeepCopy(cleanTube.GetOutput())

        stiffener = vtk.vtkPolyData()
        stiffener.DeepCopy(cleanStiffener.GetOutput())

        return annulusFittedModel, stiffener

    def buildMoldHalves(self, extractedSurface, midClippingPlane, baseClippingPlane):
        """
        Constructs the top half of mold by extruding inwards towards annulus center, and bottom half of mold by extruding downwards and clipping.
        :param extractedSurface: Extarcted inner surface vtkPolyData model
        :param midClippingPlane: vtkPlane definition of middle surface plane
        :param baseClippingPlane: vtkPlane definition of bottom clipping plane
        :return: vtkPolyData models (topHalf, bottomHalf)
        """

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

    def addOrUpdateModel(self, model, name, tformId = None, color = None):
        """
        Add or update a Model node
        :param model: vtkPolyData model
        :param name: Name of model node
        :param tformId: Transform ID to apply to model
        :param color: Optional color of model to set
        :return: None
        """
        node = slicer.util.getFirstNodeByName(name)
        if not node:
            node = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLModelNode', name)
            node.CreateDefaultDisplayNodes()

        node.SetAndObservePolyData(model)

        if tformId:
            node.SetAndObserveTransformNodeID(tformId)

        if color:
            node.GetDisplayNode().SetColor(color)


    def runDeepMV(self, heartValveNode,  volumeNode, outputSeg):
        try:
            import monai
            import torch
            from monai.data import Dataset, DataLoader
            from monai.transforms import (Compose, LoadNiftid, Orientationd, ScaleIntensityd,
                                          AddChanneld, ToTensord, CropForegroundd, Spacingd, RandSpatialCropSamplesd,
                                          RandAffined, RandCropByPosNegLabeld, AsDiscreted, Rand3DElasticd,
                                          LabelToContour)
            from monai.handlers import (CheckpointLoader, SegmentationSaver)
            from monai.networks.layers import Norm
            from monai.networks import predict_segmentation
            from monai.networks.nets import UNet
            from monai.inferers import SlidingWindowInferer
            from monai.engines import SupervisedEvaluator
            import shutil
        except ImportError as error:
            if platform.system() == 'Darwin':
                slicer.util.pip_install('torch==1.8.1 torchvision==0.9.1')
            else:
                slicer.util.pip_install('torch==1.8.1+cpu torchvision==0.9.1+cpu -f https://download.pytorch.org/whl/torch_stable.html')
            slicer.util.pip_install('monai[nibabel,skimage,pillow,gdown,ignite,torchvision,tqdm,lmdb,psutil,tensorboard]==0.3')
            importlib.invalidate_caches()
            import monai
            import torch
            from monai.data import Dataset, DataLoader
            from monai.transforms import (Compose, LoadNiftid, Orientationd, ScaleIntensityd,
                                          AddChanneld, ToTensord, CropForegroundd, Spacingd, RandSpatialCropSamplesd,
                                          RandAffined, RandCropByPosNegLabeld, AsDiscreted, Rand3DElasticd,
                                          LabelToContour)
            from monai.handlers import (CheckpointLoader, SegmentationSaver)
            from monai.networks.layers import Norm
            from monai.networks import predict_segmentation
            from monai.networks.nets import UNet
            from monai.inferers import SlidingWindowInferer
            from monai.engines import SupervisedEvaluator
            import shutil

        if monai.__version__ != '0.3.0':
            if platform.system() == 'Darwin':
                slicer.util.pip_install('torch==1.8.1 torchvision==0.9.1')
            else:
                slicer.util.pip_install(
                    'torch==1.8.1+cpu torchvision==0.9.1+cpu -f https://download.pytorch.org/whl/torch_stable.html')
            slicer.util.pip_install('monai[nibabel,skimage,pillow,gdown,ignite,torchvision,tqdm,lmdb,psutil,tensorboard]==0.3')
            logging.error("Requires MONAI version 0.3. Please restart Slicer.")
            return

        valveModel = HeartValveLib.getValveModel(heartValveNode)

        if valveModel.getProbeToRasTransformNode():
            outputSeg.SetAndObserveTransformNodeID(valveModel.getProbeToRasTransformNode().GetID())

        if not outputSeg.GetNodeReference(outputSeg.GetReferenceImageGeometryReferenceRole()):
            outputSeg.SetReferenceImageGeometryParameterFromVolumeNode(volumeNode)

        # Get tempfile path
        path = Path(slicer.app.temporaryPath).joinpath('deepmv')
        path.mkdir(parents=True, exist_ok=True)


        inputFilePath = path.joinpath('{}.nii'.format(volumeNode.GetName()))
        img = sitkUtils.PullVolumeFromSlicer(volumeNode)
        sitk.WriteImage(img, str(inputFilePath), True)

        monai.config.print_config()
        print(str(path))

        # Need to first write image here as nifti

        # Load and segment images
        images = [str(p.absolute()) for p in path.glob("*.nii")]
        d = [{"image": im} for im in images]
        keys = ("image")

        # Define transforms for image and segmentation
        xform = Compose([
            LoadNiftid(keys),
            AddChanneld(keys),
            Spacingd(keys, 0.3, diagonal=True, mode='bilinear'),
            Orientationd(keys, axcodes='RAS'),
            ScaleIntensityd("image"),
            CropForegroundd(keys, source_key="image"),
            ToTensord(keys)
        ])

        # ds = CacheDataset(d, xform)
        ds = Dataset(d, xform)
        loader = DataLoader(ds, batch_size=1, shuffle=False, num_workers=0)

        monai.utils.first(loader)

        device = torch.device('cpu')
        # net = UNet(dimensions=3, in_channels=1, out_channels=1, channels=(16, 32, 64, 128, 256),
        #            strides=(2, 2, 2, 2), num_res_units=2, norm=Norm.BATCH).to(device)

        modelPath = Path(__file__).parent.joinpath('Resources/model.md')
        net = torch.load(str(modelPath))

        evaluator = SupervisedEvaluator(
            device=device,
            val_data_loader=loader,
            network=net,
            inferer=SlidingWindowInferer((96, 96, 96), sw_batch_size=6),
        )

        # modelPath = Path(__file__).parent.joinpath('Resources/model.pt')
        # checkpoint = torch.load(str(modelPath), map_location=device)
        # net.load_state_dict(checkpoint['net'])
        #
        # modelPath2 = modelPath.parent.joinpath('model.md')
        # torch.save(net, str(modelPath2))
        #
        #
        # checkpoint_loader = CheckpointLoader(str(modelPath), {'net': net}, map_location=device)
        # checkpoint_loader.attach(evaluator)

        prediction_saver = SegmentationSaver(
            output_dir=path,
            name="evaluator",
            dtype=np.dtype('float64'),
            batch_transform=lambda batch: batch["image_meta_dict"],
            output_transform=lambda output: predict_segmentation(output['pred'])
        )
        prediction_saver.attach(evaluator)

        # Running segmentation
        evaluator.run()

        # Retrieve segmentations afterwards and include them into segmentation node
        segPath = path.joinpath(volumeNode.GetName()).joinpath('{}_seg.nii.gz'.format(volumeNode.GetName()))
        seg = sitk.ReadImage(str(segPath))
        self.pushITKImageToSegmentation(seg, outputSeg, 'Leaflet Segmentation')

        # Cleanup temporary files
        shutil.rmtree(str(path))




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
