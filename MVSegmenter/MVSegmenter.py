import os
import unittest
import vtk, qt, ctk, slicer
import SimpleITK as sitk
import sitkUtils
from slicer.ScriptedLoadableModule import *
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
        self.parent.categories = ["Examples"]
        self.parent.dependencies = []
        self.parent.contributors = ["Patrick Carnahan (Robarts Research Institute)"] # replace with "Firstname Lastname (Organization)"
        self.parent.helpText = """
This module implements an algorithm for automatic mitral valve segmentation using ITK.
"""
        self.parent.helpText += self.getDefaultModuleDocumentationLink()
        self.parent.acknowledgementText = """
This file was originally developed by Jean-Christophe Fillion-Robin, Kitware Inc.
and Steve Pieper, Isomics, Inc. and was partially funded by NIH grant 3P41RR013218-12S1.
""" # replace with organization, grant and thanks.

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
        self.inputSelector.setMRMLScene( slicer.mrmlScene )
        self.inputSelector.setToolTip( "Pick the input to the algorithm." )
        parametersFormLayout.addRow("Input Volume", self.inputSelector)

        #
        # input mask selector
        #
        self.inputMask = slicer.qMRMLNodeComboBox()
        self.inputMask.nodeTypes = ["vtkMRMLLabelMapVolumeNode"]
        self.inputMask.selectNodeUponCreation = True
        self.inputMask.addEnabled = True
        self.inputMask.removeEnabled = False
        self.inputMask.noneEnabled = False
        self.inputMask.showHidden = False
        self.inputMask.showChildNodeTypes = False
        self.inputMask.setMRMLScene(slicer.mrmlScene)
        self.inputMask.setToolTip("Pick the input to the algorithm.")
        parametersFormLayout.addRow("Input Mask", self.inputMask)

        #
        # output volume selector
        #
        self.outputSelector = slicer.qMRMLNodeComboBox()
        self.outputSelector.nodeTypes = ["vtkMRMLLabelMapVolumeNode"]
        self.outputSelector.selectNodeUponCreation = True
        self.outputSelector.addEnabled = True
        self.outputSelector.removeEnabled = True
        self.outputSelector.noneEnabled = True
        self.outputSelector.showHidden = False
        self.outputSelector.showChildNodeTypes = False
        self.outputSelector.setMRMLScene( slicer.mrmlScene )
        self.outputSelector.setToolTip( "Pick the output to the algorithm." )
        parametersFormLayout.addRow("Output Volume", self.outputSelector)


        #
        # check box to trigger taking screen shots for later use in tutorials
        #
        self.enableScreenshotsFlagCheckBox = qt.QCheckBox()
        self.enableScreenshotsFlagCheckBox.checked = 0
        self.enableScreenshotsFlagCheckBox.setToolTip("If checked, take screen shots for tutorials. Use Save Data to write them to disk.")
        parametersFormLayout.addRow("Enable Screenshots", self.enableScreenshotsFlagCheckBox)


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

        firstPassFormLayout.addRow("Increment Segmentation", incrementFirstHBox)

        #
        #  Second Phase Segmentation
        #
        secondPassCollapsibleButton = ctk.ctkCollapsibleButton()
        secondPassCollapsibleButton.text = "Leaflet Segmentation"
        self.layout.addWidget(secondPassCollapsibleButton)

        # Layout within the dummy collapsible button
        secondPassFormLayout = qt.QFormLayout(secondPassCollapsibleButton)

        #Initialize
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

        self.undoButton = qt.QPushButton("Undo")
        self.undoButton.toolTip = "Undo previous step"
        self.undoButton.enabled = False
        incrementHBox.addWidget(self.undoButton)

        self.redoButton = qt.QPushButton("Redo")
        self.redoButton.toolTip = "Redo previous step"
        self.redoButton.enabled = False
        incrementHBox.addWidget(self.redoButton)

        secondPassFormLayout.addRow("Increment Segmentation", incrementHBox)

        # connections
        self.initBPButton.connect('clicked(bool)', self.onInitBPButton)
        self.incrementFirstButton50.connect('clicked(bool)', self.onIncrementFirst50Button)
        self.incrementFirstButton100.connect('clicked(bool)', self.onIncrementFirst100Button)
        self.incrementFirstButton500.connect('clicked(bool)', self.onIncrementFirst500Button)

        self.initLeafletButton.connect('clicked(bool)', self.onInitLeafletButton)
        self.incrementButton10.connect('clicked(bool)', self.onIncrement10Button)
        self.incrementButton50.connect('clicked(bool)', self.onIncrement50Button)
        self.undoButton.connect('clicked(bool)', self.onUndoButton)
        self.redoButton.connect('clicked(bool)', self.onRedoButton)

        self.inputSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.onSelect)
        self.inputMask.connect("currentNodeChanged(vtkMRMLNode*)", self.onSelect)
        self.outputSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.onSelect)


        # Add vertical spacer
        self.layout.addStretch(1)

        # Refresh Apply button state
        self.onSelect()

    def cleanup(self):
        pass

    def onSelect(self):
        self.initBPButton.enabled = self.inputSelector.currentNode() and self.inputMask.currentNode() and self.outputSelector.currentNode()

    def onInitBPButton(self):
        try:
            # This can be a long operation - indicate it to the user
            qt.QApplication.setOverrideCursor(qt.Qt.WaitCursor)

            enableScreenshotsFlag = self.enableScreenshotsFlagCheckBox.checked
            self.logic.initBPSeg(self.inputSelector.currentNode(), self.inputMask.currentNode(),
                                 self.outputSelector.currentNode())
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

            self.logic.iterateFirstPass(50, self.outputSelector.currentNode())
        finally:
            qt.QApplication.restoreOverrideCursor()

    def onIncrementFirst100Button(self):
        try:
            # This can be a long operation - indicate it to the user
            qt.QApplication.setOverrideCursor(qt.Qt.WaitCursor)

            self.logic.iterateFirstPass(100, self.outputSelector.currentNode())
        finally:
            qt.QApplication.restoreOverrideCursor()

    def onIncrementFirst500Button(self):
        try:
            # This can be a long operation - indicate it to the user
            qt.QApplication.setOverrideCursor(qt.Qt.WaitCursor)

            self.logic.iterateFirstPass(500, self.outputSelector.currentNode())

        finally:
            qt.QApplication.restoreOverrideCursor()

    def onInitLeafletButton(self):
        try:
            # This can be a long operation - indicate it to the user
            qt.QApplication.setOverrideCursor(qt.Qt.WaitCursor)

            self.logic.initLeafletSeg(self.outputSelector.currentNode())
            self.incrementButton10.enabled = True
            self.incrementButton50.enabled = True

        finally:
            qt.QApplication.restoreOverrideCursor()

    def onIncrement10Button(self):
        self.logic.iterateSecondPass(10, self.outputSelector.currentNode())
        self.undoButton.enabled = True

    def onIncrement50Button(self):
        try:
            # This can be a long operation - indicate it to the user
            qt.QApplication.setOverrideCursor(qt.Qt.WaitCursor)

            self.logic.iterateSecondPass(50, self.outputSelector.currentNode())
            self.undoButton.enabled = True

        finally:
            qt.QApplication.restoreOverrideCursor()

    def onUndoButton(self):
        self.logic.undoIteration(self.outputSelector.currentNode())
        self.redoButton.enabled = True
        self.undoButton.enabled = False

    def onRedoButton(self):
        self.logic.redoIteration(self.outputSelector.currentNode())
        self.redoButton.enabled = False
        self.undoButton.enabled = True

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

        self._prevLevelSet = None
        self._nextLevelSet = None

    def hasImageData(self,volumeNode):
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
        if inputVolumeNode.GetID()==outputVolumeNode.GetID():
            logging.debug('isValidInputOutputData failed: input and output volume is the same. Create a new volume for output to avoid this error.')
            return False
        return True

    def takeScreenshot(self,name,description,type=-1):
        # show the message even if not taking a screen shot
        slicer.util.delayDisplay('Take screenshot: '+description+'.\nResult is available in the Annotations module.', 3000)

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
        slicer.qMRMLUtils().qImageToVtkImageData(qimage,imageData)

        annotationLogic = slicer.modules.annotations.logic()
        annotationLogic.CreateSnapShot(name, description, type, 1, imageData)

    def initBPSeg(self, inputVolume, inputMask, outputVolume):
        """
        Run the actual algorithm
        """

        # calculate speed image from input volume
        # Uses DiscreteGaussian -> GradientMagnitude -> Sigmoid filters
        speedImg = sitkUtils.PullVolumeFromSlicer(inputVolume)

        blurFilter = sitk.DiscreteGaussianImageFilter()
        blurFilter.SetMaximumError(0.25)
        blurFilter.SetMaximumKernelWidth(64)
        blurFilter.SetUseImageSpacing(False)
        blurFilter.SetVariance(5)
        speedImg = blurFilter.Execute(speedImg)

        gradMag = sitk.GradientMagnitudeImageFilter()
        gradMag.SetUseImageSpacing(False)
        speedImg = gradMag.Execute(speedImg)

        sigmoid = sitk.SigmoidImageFilter()
        sigmoid.SetOutputMinimum(0)
        sigmoid.SetOutputMaximum(1.0)
        sigmoid.SetAlpha(-5.0)
        sigmoid.SetBeta(10.0)
        speedImg = sigmoid.Execute(speedImg)

        self._speedImg = speedImg

        # compute initial level set
        levelSet = sitkUtils.PullVolumeFromSlicer(inputMask)

        signedDis = sitk.SignedDanielssonDistanceMapImageFilter()
        levelSet = signedDis.Execute(levelSet)


        # Run first pass of geodesic active contour
        # TODO adjust active contour parameters / make user entered as different image data may need different values to behave
        # goal is to find middle ground values that work for most cases
        # could try simplified parameters aimed at fixing leaks/too small and adjust real parameters as needed here
        geodesicActiveContour = sitk.GeodesicActiveContourLevelSetImageFilter()
        geodesicActiveContour.SetCurvatureScaling(2.7)
        geodesicActiveContour.SetAdvectionScaling(0.8)
        geodesicActiveContour.SetPropagationScaling(3)
        geodesicActiveContour.SetMaximumRMSError(0.0001)
        geodesicActiveContour.SetNumberOfIterations(1000)

        out_mask = geodesicActiveContour.Execute(levelSet, speedImg)
        print(geodesicActiveContour.GetElapsedIterations())

        self._bpLevelSet = out_mask

        threshold = sitk.BinaryThresholdImageFilter()
        threshold.SetInsideValue(1)
        threshold.SetLowerThreshold(-1000.0)
        threshold.SetOutsideValue(0)
        threshold.SetUpperThreshold(0.0)
        out_mask = threshold.Execute(out_mask)

        sitkUtils.PushVolumeToSlicer(out_mask, outputVolume)

    def iterateFirstPass(self, nIter, outputVolume):
        geodesicActiveContour = sitk.GeodesicActiveContourLevelSetImageFilter()
        geodesicActiveContour.SetCurvatureScaling(2.7)
        geodesicActiveContour.SetAdvectionScaling(0.9)
        geodesicActiveContour.SetPropagationScaling(1.8)
        geodesicActiveContour.SetMaximumRMSError(0.00001)
        geodesicActiveContour.SetNumberOfIterations(nIter)

        out_mask = geodesicActiveContour.Execute(self._bpLevelSet, self._speedImg)
        print(geodesicActiveContour.GetElapsedIterations())

        self._bpLevelSet = out_mask

        threshold = sitk.BinaryThresholdImageFilter()
        threshold.SetInsideValue(1)
        threshold.SetLowerThreshold(-1000.0)
        threshold.SetOutsideValue(0)
        threshold.SetUpperThreshold(0.0)
        out_mask = threshold.Execute(out_mask)

        sitkUtils.PushVolumeToSlicer(out_mask, outputVolume)

    def initLeafletSeg(self, outputVolume):

        # Get region bordering initial blood-pool segmentation
        threshold = sitk.BinaryThresholdImageFilter()
        threshold.SetInsideValue(1)
        threshold.SetLowerThreshold(-1000.0)
        threshold.SetOutsideValue(0)
        threshold.SetUpperThreshold(0.0)
        distMap = threshold.Execute(self._bpLevelSet)

        distFilter = sitk.DanielssonDistanceMapImageFilter()
        distFilter.SetInputIsBinary(True)
        distFilter.SetSquaredDistance(False)
        distFilter.SetUseImageSpacing(False)
        distMap = distFilter.Execute(distMap)

        distThreshold = sitk.BinaryThresholdImageFilter()
        distThreshold.SetInsideValue(1)
        distThreshold.SetLowerThreshold(0.5)
        distThreshold.SetOutsideValue(0)
        distThreshold.SetUpperThreshold(15)
        leafletMask = distThreshold.Execute(distMap)

        # Run second pass to get final leaflet segmentation
        signedDis = sitk.SignedDanielssonDistanceMapImageFilter()
        levelSet = signedDis.Execute(leafletMask)

        geodesicActiveContour2 = sitk.GeodesicActiveContourLevelSetImageFilter()
        geodesicActiveContour2.SetCurvatureScaling(0.8)
        geodesicActiveContour2.SetAdvectionScaling(0.1)
        geodesicActiveContour2.SetPropagationScaling(-0.4)
        geodesicActiveContour2.SetMaximumRMSError(0.0001)
        geodesicActiveContour2.SetNumberOfIterations(400)
        out_mask = geodesicActiveContour2.Execute(levelSet, self._speedImg)

        self._levelSet = out_mask

        threshold = sitk.BinaryThresholdImageFilter()
        threshold.SetInsideValue(1)
        threshold.SetLowerThreshold(-1000.0)
        threshold.SetOutsideValue(0)
        threshold.SetUpperThreshold(0.0)
        out_mask = threshold.Execute(out_mask)

        print(geodesicActiveContour2.GetElapsedIterations())

        sitkUtils.PushVolumeToSlicer(out_mask, outputVolume)

        return out_mask

    def iterateSecondPass(self, nIter, outputVolumeNode):

        geodesicActiveContour2 = sitk.GeodesicActiveContourLevelSetImageFilter()
        geodesicActiveContour2.SetCurvatureScaling(0.8)
        geodesicActiveContour2.SetAdvectionScaling(0.1)
        geodesicActiveContour2.SetPropagationScaling(-0.4)
        geodesicActiveContour2.SetMaximumRMSError(0.0001)
        geodesicActiveContour2.SetNumberOfIterations(nIter)
        out_mask = geodesicActiveContour2.Execute(self._levelSet, self._speedImg)

        self._prevLevelSet = self._levelSet
        self._levelSet = out_mask

        print(geodesicActiveContour2.GetElapsedIterations())

        threshold = sitk.BinaryThresholdImageFilter()
        threshold.SetInsideValue(1)
        threshold.SetLowerThreshold(-1000.0)
        threshold.SetOutsideValue(0)
        threshold.SetUpperThreshold(0.0)
        out_mask = threshold.Execute(out_mask)

        sitkUtils.PushVolumeToSlicer(out_mask, outputVolumeNode)

        return out_mask

    def undoIteration(self, outputVolumeNode):

        self._nextLevelSet = self._levelSet
        self._levelSet = self._prevLevelSet

        threshold = sitk.BinaryThresholdImageFilter()
        threshold.SetInsideValue(1)
        threshold.SetLowerThreshold(-1000.0)
        threshold.SetOutsideValue(0)
        threshold.SetUpperThreshold(0.0)

        sitkUtils.PushVolumeToSlicer(threshold.Execute(self._levelSet), outputVolumeNode)

        return self._levelSet

    def redoIteration(self, outputVolumeNode):

        self._prevLevelSet = self._levelSet
        self._levelSet = self._nextLevelSet

        threshold = sitk.BinaryThresholdImageFilter()
        threshold.SetInsideValue(1)
        threshold.SetLowerThreshold(-1000.0)
        threshold.SetOutsideValue(0)
        threshold.SetUpperThreshold(0.0)

        sitkUtils.PushVolumeToSlicer(threshold.Execute(self._levelSet), outputVolumeNode)

        return self._levelSet



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

        self.delayDisplay("Starting the test")
        #
        # first, get some data
        #
        import urllib
        downloads = (
            ('http://slicer.kitware.com/midas3/download?items=5767', 'FA.nrrd', slicer.util.loadVolume),
        )

        for url,name,loader in downloads:
            filePath = slicer.app.temporaryPath + '/' + name
            if not os.path.exists(filePath) or os.stat(filePath).st_size == 0:
                logging.info('Requesting download %s from %s...\n' % (name, url))
                urllib.urlretrieve(url, filePath)
            if loader:
                logging.info('Loading %s...' % (name,))
                loader(filePath)
        self.delayDisplay('Finished with download and loading')

        volumeNode = slicer.util.getNode(pattern="FA")
        logic = MVSegmenterLogic()
        self.assertIsNotNone( logic.hasImageData(volumeNode) )
        self.delayDisplay('Test passed!')
