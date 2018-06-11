import os
import unittest
import vtk, qt, ctk, slicer
from slicer.ScriptedLoadableModule import *
import slicer.util
import logging
import numpy as np
from slicer.util import VTKObservationMixin

#
# BiplaneRegistration
#

class BiplaneRegistration(ScriptedLoadableModule):
    """Uses ScriptedLoadableModule base class, available at:
    https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
    """

    def __init__(self, parent):
        ScriptedLoadableModule.__init__(self, parent)
        self.parent.title = "Biplane Registration" #
        self.parent.categories = ["Examples"]
        self.parent.dependencies = []
        self.parent.contributors = ["Patrick Carnahan"] # replace with "Firstname Lastname (Organization)"
        self.parent.helpText = """
This is a module that allows a philips bi-plane image exported as jpg to be converted and registered in 3D space.
"""
        self.parent.helpText += self.getDefaultModuleDocumentationLink()
        self.parent.acknowledgementText = """
This file was originally developed by Jean-Christophe Fillion-Robin, Kitware Inc.
and Steve Pieper, Isomics, Inc. and was partially funded by NIH grant 3P41RR013218-12S1.
""" # replace with organization, grant and thanks.

#
# BiplaneRegistrationWidget
#

class BiplaneRegistrationWidget(ScriptedLoadableModuleWidget, VTKObservationMixin):
    """Uses ScriptedLoadableModuleWidget base class, available at:
    https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
    """

    def __init__(self, parent):
        ScriptedLoadableModuleWidget.__init__(self, parent)
        VTKObservationMixin.__init__(self)

        self.logic = BiplaneRegistrationLogic()
        self.xTransformObserver = None
        self.yTransformObserver = None

    def setup(self):
        ScriptedLoadableModuleWidget.setup(self)

        # Instantiate and connect widgets ...

        #
        # Parameters Area
        #
        parametersCollapsibleButton = ctk.ctkCollapsibleButton()
        parametersCollapsibleButton.text = "Initialize"
        self.layout.addWidget(parametersCollapsibleButton)

        # Layout within the dummy collapsible button
        parametersFormLayoutIO = qt.QFormLayout(parametersCollapsibleButton)

        #
        # input volume selectors
        #
        self.inputSelectorBiPlane = slicer.qMRMLNodeComboBox()
        self.inputSelectorBiPlane.nodeTypes = ["vtkMRMLScalarVolumeNode", "vtkMRMLVectorVolumeNode"]
        self.inputSelectorBiPlane.selectNodeUponCreation = True
        self.inputSelectorBiPlane.addEnabled = False
        self.inputSelectorBiPlane.removeEnabled = False
        self.inputSelectorBiPlane.noneEnabled = False
        self.inputSelectorBiPlane.showHidden = False
        self.inputSelectorBiPlane.showChildNodeTypes = False
        self.inputSelectorBiPlane.setMRMLScene(slicer.mrmlScene)
        self.inputSelectorBiPlane.setToolTip("Pick the input to the algorithm.")
        parametersFormLayoutIO.addRow("BiPlane Volume: ", self.inputSelectorBiPlane)


        self.inputSelectorFixed = slicer.qMRMLNodeComboBox()
        self.inputSelectorFixed.nodeTypes = ["vtkMRMLScalarVolumeNode"]
        self.inputSelectorFixed.selectNodeUponCreation = True
        self.inputSelectorFixed.addEnabled = False
        self.inputSelectorFixed.removeEnabled = False
        self.inputSelectorFixed.noneEnabled = False
        self.inputSelectorFixed.showHidden = False
        self.inputSelectorFixed.showChildNodeTypes = False
        self.inputSelectorFixed.setMRMLScene(slicer.mrmlScene)
        self.inputSelectorFixed.setToolTip("Pick the input to the algorithm.")
        parametersFormLayoutIO.addRow("Fixed Volume: ", self.inputSelectorFixed)

        #
        # output volume selector
        #
        self.outputSelector = slicer.qMRMLNodeComboBox()
        self.outputSelector.nodeTypes = ["vtkMRMLVectorVolumeNode"]
        self.outputSelector.selectNodeUponCreation = True
        self.outputSelector.addEnabled = True
        self.outputSelector.removeEnabled = True
        self.outputSelector.noneEnabled = True
        self.outputSelector.showHidden = False
        self.outputSelector.showChildNodeTypes = False
        self.outputSelector.setMRMLScene( slicer.mrmlScene )
        self.outputSelector.setToolTip( "Pick the output to the algorithm." )
        parametersFormLayoutIO.addRow("Output Volume: ", self.outputSelector)

        self.outputSelector2 = slicer.qMRMLNodeComboBox()
        self.outputSelector2.nodeTypes = ["vtkMRMLVectorVolumeNode"]
        self.outputSelector2.selectNodeUponCreation = True
        self.outputSelector2.addEnabled = True
        self.outputSelector2.removeEnabled = True
        self.outputSelector2.noneEnabled = True
        self.outputSelector2.showHidden = False
        self.outputSelector2.showChildNodeTypes = False
        self.outputSelector2.setMRMLScene(slicer.mrmlScene)
        self.outputSelector2.setToolTip("Pick the output to the algorithm.")
        parametersFormLayoutIO.addRow("Output Volume: ", self.outputSelector2)


        planeAnglesHBoxLayout = qt.QHBoxLayout()
        planeAnglesHBoxLayout.addStretch(5)

        self.onlyInt = qt.QIntValidator()
        self.onlyInt.setRange(0,90)

        self.sagittalRotationInput = qt.QLineEdit()
        self.sagittalRotationInput.setValidator(self.onlyInt)
        self.sagittalRotationInput.setFixedWidth(60)
        planeAnglesHBoxLayout.addWidget(self.sagittalRotationInput)

        planeAnglesHBoxLayout.addStretch(1)

        self.coronalRotationInput = qt.QLineEdit()
        self.coronalRotationInput.setValidator(self.onlyInt)
        self.coronalRotationInput.setFixedWidth(60)
        planeAnglesHBoxLayout.addWidget(self.coronalRotationInput)
        planeAnglesHBoxLayout.addStretch(5)

        parametersFormLayoutIO.addRow("Plane Angles", planeAnglesHBoxLayout)


        #
        # check box to trigger taking screen shots for later use in tutorials
        #
        self.enableSliceLockFlagCheckBox = qt.QCheckBox()
        self.enableSliceLockFlagCheckBox.checked = 1
        self.enableSliceLockFlagCheckBox.setToolTip("If checked, lock yellow and green slices to the 2 moving images.")
        parametersFormLayoutIO.addRow("Lock Slices to Volumes", self.enableSliceLockFlagCheckBox)

        #
        # Apply Button
        #
        self.initTransformButton = qt.QPushButton("Apply")
        self.initTransformButton.toolTip = "Run the algorithm."
        self.initTransformButton.enabled = False
        parametersFormLayoutIO.addRow(self.initTransformButton)


        ###############################################################################################################
        #
        # Auto Registration Section
        #
        ###############################################################################################################
        parametersCollapsibleButton = ctk.ctkCollapsibleButton()
        parametersCollapsibleButton.text = "Auto Register"
        self.layout.addWidget(parametersCollapsibleButton)

        # Layout within the dummy collapsible button
        parametersFormLayoutReg = qt.QFormLayout(parametersCollapsibleButton)

        self.outputGreyVolume = slicer.qMRMLNodeComboBox()
        self.outputGreyVolume.nodeTypes = ["vtkMRMLScalarVolumeNode"]
        self.outputGreyVolume.selectNodeUponCreation = True
        self.outputGreyVolume.addEnabled = True
        self.outputGreyVolume.removeEnabled = True
        self.outputGreyVolume.noneEnabled = True
        self.outputGreyVolume.showHidden = False
        self.outputGreyVolume.showChildNodeTypes = False
        self.outputGreyVolume.setMRMLScene(slicer.mrmlScene)
        self.outputGreyVolume.setToolTip("Pick the output to the algorithm.")
        parametersFormLayoutReg.addRow("Output Volume: ", self.outputGreyVolume)

        self.convertGreyscaleButton = qt.QPushButton("Convert")
        self.convertGreyscaleButton.toolTip = "Extract greyscale image."
        self.convertGreyscaleButton.enabled = True
        parametersFormLayoutReg.addRow(self.convertGreyscaleButton)
        parametersCollapsibleButton.hide()

        ###############################################################################################################
        #
        # Edit Transform Section
        #
        ###############################################################################################################
        editCollapsibleButton = ctk.ctkCollapsibleButton()
        editCollapsibleButton.text = "Edit Transform"
        self.layout.addWidget(editCollapsibleButton)

        editFormLayout = qt.QFormLayout(editCollapsibleButton)

        self.editTransform = slicer.qMRMLNodeComboBox()
        self.editTransform.nodeTypes = ["vtkMRMLTransformNode"]
        self.editTransform.selectNodeUponCreation = True
        self.editTransform.addEnabled = False
        self.editTransform.removeEnabled = False
        self.editTransform.noneEnabled = False
        self.editTransform.showHidden = False
        self.editTransform.showChildNodeTypes = False
        self.editTransform.setMRMLScene(slicer.mrmlScene)
        self.editTransform.setToolTip("Pick the output to the algorithm.")
        editFormLayout.addRow("Transform: ", self.editTransform)

        self.translationSliders = slicer.qMRMLTransformSliders()
        self.translationSliders.SingleStep = 0.5
        editFormLayout.addRow(self.translationSliders)

        self.rotationSliders = slicer.qMRMLTransformSliders()
        self.rotationSliders.TypeOfTransform = self.rotationSliders.ROTATION
        self.rotationSliders.Title = "Rotation"
        self.rotationSliders.SingleStep = 0.5
        editFormLayout.addRow(self.rotationSliders)

        # connections
        self.initTransformButton.connect('clicked(bool)', self.onInitTransformButton)
        self.inputSelectorFixed.connect("currentNodeChanged(vtkMRMLNode*)", self.onSelect)
        self.inputSelectorBiPlane.connect("currentNodeChanged(vtkMRMLNode*)", self.onSelect)
        self.outputSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.onSelect)
        self.outputSelector2.connect("currentNodeChanged(vtkMRMLNode*)", self.onSelect)
        self.sagittalRotationInput.connect("textEdited(QString)", self.onSelect)
        self.coronalRotationInput.connect("textEdited(QString)", self.onSelect)
        self.enableSliceLockFlagCheckBox.connect("stateChanged(int)", self.onEnableSliceLockFlagChanged)

        self.editTransform.connect("currentNodeChanged(vtkMRMLNode*)", self.onTransformSelect)

        self.convertGreyscaleButton.connect('clicked(bool)', self.onConvertGreyscaleButtonClicked)

        # Add vertical spacer
        self.layout.addStretch(1)

        # Refresh Apply button state
        self.onSelect()

    def cleanup(self):
        if self.xTransformObserver:
            self.RemoveObserver(self.xTransformObserver)
        if self.yTransformObserver:
            self.RemoveObserver(self.yTransformObserver)

    def onSelect(self):
        self.initTransformButton.enabled = self.inputSelectorBiPlane.currentNode() and self.inputSelectorFixed.currentNode() \
                                           and self.outputSelector.currentNode() and self.outputSelector2.currentNode() \
                                            and self.sagittalRotationInput.text.isnumeric() and self.coronalRotationInput.text.isnumeric()

    def onTransformSelect(self):
        self.translationSliders.enabled = self.editTransform.currentNode()
        self.rotationSliders.enabled = self.translationSliders.enabled

        node = self.editTransform.currentNode()
        self.translationSliders.setMRMLTransformNode(node)
        self.rotationSliders.setMRMLTransformNode(node)
        self.translationSliders.update()
        self.rotationSliders.update()

    def onInitTransformButton(self):
        self.logic.initVolumes(self.inputSelectorBiPlane.currentNode(), None, self.outputSelector2.currentNode(),
                               self.outputSelector.currentNode())

        self.logic.initTransform(self.inputSelectorFixed.currentNode(), self.outputSelector.currentNode(),
                                 self.outputSelector2.currentNode(), int(self.sagittalRotationInput.text),
                                 int(self.coronalRotationInput.text))

        # add observer for xplane transform
        if self.outputSelector.currentNode().GetTransformNodeID() and slicer.mrmlScene.GetNodeByID(
                self.outputSelector.currentNode().GetTransformNodeID()):
            xTransform = slicer.mrmlScene.GetNodeByID(self.outputSelector.currentNode().GetTransformNodeID())

            # remove existing observer
            if self.xTransformObserver:
                self.RemoveObserver(self.xTransformObserver)

            self.xTransformObserver = self.addObserver(xTransform, slicer.vtkMRMLTransformableNode.TransformModifiedEvent,
                             self.onVolumeXTransformNodeModified)

        # add observer for yplane transform
        if self.outputSelector2.currentNode().GetTransformNodeID() and slicer.mrmlScene.GetNodeByID(
                self.outputSelector2.currentNode().GetTransformNodeID()):
            yTransform = slicer.mrmlScene.GetNodeByID(self.outputSelector2.currentNode().GetTransformNodeID())

            # remove existing observer
            if self.yTransformObserver:
                self.RemoveObserver(self.yTransformObserver)

            self.yTransformObserver = self.addObserver(yTransform, slicer.vtkMRMLTransformableNode.TransformModifiedEvent,
                             self.onVolumeYTransformNodeModified)

    def onVolumeXTransformNodeModified(self, observer, eventid):
        if self.enableSliceLockFlagCheckBox.checked:
            self.logic.alignYellowSlice(self.outputSelector.currentNode())

    def onVolumeYTransformNodeModified(self, observer, eventid):
        if self.enableSliceLockFlagCheckBox.checked:
            self.logic.alignGreenSlice(self.outputSelector2.currentNode())

    def onEnableSliceLockFlagChanged(self):
        if self.enableSliceLockFlagCheckBox.checked:
            self.logic.alignYellowSliceSlice(self.outputSelector.currentNode())
            self.logic.alignGreenSlice(self.outputSelector2.currentNode())

    def onConvertGreyscaleButtonClicked(self):
        self.logic.extractGreyImage(self.outputSelector.currentNode(), self.outputGreyVolume.currentNode())

#
# BiplaneRegistrationLogic
#

class BiplaneRegistrationLogic(ScriptedLoadableModuleLogic):
    """This class should implement all the actual
    computation done by your module.  The interface
    should be such that other python code can import
    this class and make use of the functionality without
    requiring an instance of the Widget.
    Uses ScriptedLoadableModuleLogic base class, available at:
    https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
    """

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

    def initVolumes(self, inputVolumeBiPlane, inputVolumeFixed, outVolumeY, outVolumeX):
        # extract 2 image planes form side-by-side biplane
        # using fixed values for 1024x768 philips biplane
        biplane = slicer.util.arrayFromVolume(inputVolumeBiPlane)
        im = biplane[0, :, :, :]

        # xplane extraction
        imx = np.zeros((355, 350, 1, 3))
        imx[:, :, 0, :] = im[235:590, 140:490]
        idx = (abs(imx[:, :, :, 1] - imx[:, :, :, 2]) < 15) * (
                abs(imx[:, :, :, 1] - imx[:, :, :, 0]) < 15)
        imxgrey = np.zeros(imx.shape)
        imxgrey[idx] = imx[idx]

        # yplane extraction
        imy = np.zeros((355, 350, 1, 3))
        imy[:, :, 0, :] = im[235:590, 606:956]

        # set volume data
        outVolumeX.SetSpacing([0.29462782549439480183368515087702] * 3)     # spacing computed from markers in dicom file which have 10mm spacing
        slicer.util.updateVolumeFromArray(outVolumeX, imx)
        outVolumeY.SetSpacing([0.29462782549439480183368515087702] * 3)
        slicer.util.updateVolumeFromArray(outVolumeY, imy)

        # set slice views to show created volumes
        yellow_logic = slicer.app.layoutManager().sliceWidget("Yellow").sliceLogic()
        yellow_logic.GetSliceCompositeNode().SetBackgroundVolumeID(outVolumeX.GetID())
        yellow_logic.GetSliceCompositeNode().SetForegroundVolumeID(inputVolumeFixed.GetID())
        yellow_logic.GetSliceCompositeNode().SetForegroundOpacity(0.5)

        green_logic = slicer.app.layoutManager().sliceWidget("Green").sliceLogic()
        green_logic.GetSliceCompositeNode().SetBackgroundVolumeID(outVolumeY.GetID())
        green_logic.GetSliceCompositeNode().SetForegroundVolumeID(inputVolumeFixed.GetID())
        green_logic.GetSliceCompositeNode().SetForegroundOpacity(0.5)

    def initTransform(self, inputVolumeFixed, outVolumeX, outVolumeY, xAngle, yAngle):
        # initialize transform nodes if they don't already exist
        # need individual transforms for 2 planar images plus a master transform to move both together
        xTransformNode = None
        if outVolumeX.GetTransformNodeID():
            xTransformNode = slicer.mrmlScene.GetNodeByID(outVolumeX.GetTransformNodeID())
        if not xTransformNode or xTransformNode.GetName() != 'ASPlaneTransform':
            xTransformNode = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLTransformNode')
            xTransformNode.SetName('ASPlaneTransform')
            outVolumeX.SetAndObserveTransformNodeID(xTransformNode.GetID())

        yTransformNode = None
        if outVolumeY.GetTransformNodeID():
            yTransformNode = slicer.mrmlScene.GetNodeByID(outVolumeY.GetTransformNodeID())
        if not yTransformNode or yTransformNode.GetName() != 'RotPlaneTransform':
            yTransformNode = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLTransformNode')
            yTransformNode.SetName('RotPlaneTransform')
            outVolumeY.SetAndObserveTransformNodeID(yTransformNode.GetID())

        parentTransformNode = None
        if xTransformNode.GetTransformNodeID() and xTransformNode.GetTransformNodeID() == yTransformNode.GetTransformNodeID():
            parentTransformNode = slicer.mrmlScene.GetNodeByID(xTransformNode.GetTransformNodeID())
        if not parentTransformNode or parentTransformNode.GetName() != 'MasterTransform':
            parentTransformNode = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLTransformNode')
            parentTransformNode.SetName('MasterTransform')
            xTransformNode.SetAndObserveTransformNodeID(parentTransformNode.GetID())
            yTransformNode.SetAndObserveTransformNodeID(parentTransformNode.GetID())


        # compute approximate initial translation based on volume centers
        fixedBounds = np.zeros(6, 'double')
        inputVolumeFixed.GetBounds(fixedBounds)
        fixedCenter = np.array([(fixedBounds[1] + fixedBounds[0]) / 2, (fixedBounds[3] + fixedBounds[2]) / 2,
                       (fixedBounds[5] + fixedBounds[4]) / 2], 'double')

        xBounds = np.zeros(6, 'double')
        outVolumeX.GetBounds(xBounds)
        xCenter = np.array([(xBounds[1] + xBounds[0]) / 2, (xBounds[3] + xBounds[2]) / 2, (xBounds[5] + xBounds[4]) / 2] , 'double')

        yBounds = np.zeros(6, 'double')
        outVolumeY.GetBounds(yBounds)
        yCenter = np.array([(yBounds[1] + yBounds[0]) / 2, (yBounds[3] + yBounds[2]) / 2, (yBounds[5] + yBounds[4]) / 2], 'double')


        # apply rotation to sagittal plane
        xTransform = vtk.vtkTransform()
        xTransform.PostMultiply()

        translationTransform = vtk.vtkTransform()
        translationTransform.Translate(0, -51.4, -7.5)      # move center marked point to origin
        xTransform.Concatenate(translationTransform)

        rotTransform = vtk.vtkTransform()
        rotTransform.RotateZ(xAngle)
        xTransform.Concatenate(rotTransform)

        xTransform.Translate(0, 51.4, 7.5)          # move center back
        xTransformNode.SetMatrixTransformToParent(xTransform.GetMatrix())

        # apply rotation to coronal plane
        yTransform = vtk.vtkTransform()
        yTransform.PostMultiply()

        translationTransform = vtk.vtkTransform()
        translationTransform.Translate(0, -51.7, -7.7) # move center marked point to origin
        yTransform.Concatenate(translationTransform)

        rotTransform = vtk.vtkTransform()
        rotTransform.RotateZ(yAngle)
        yTransform.Concatenate(rotTransform)

        yTransform.Translate(0, 51.4, 7.7) # move center back
        yTransformNode.SetMatrixTransformToParent(yTransform.GetMatrix())

        # initialize master transform

        tr = vtk.vtkTransform()
        tr.Translate(fixedCenter - xCenter)

        parentTransformNode.SetMatrixTransformToParent(tr.GetMatrix())

        self.alignYellowSlice(outVolumeX)
        self.alignGreenSlice(outVolumeY)

    def alignYellowSlice(self, volume):
        yellow = slicer.mrmlScene.GetNodeByID('vtkMRMLSliceNodeYellow')
        yellow.RotateToVolumePlane(volume)
        slicer.app.layoutManager().sliceWidget('Yellow').sliceLogic().FitSliceToAll()

    def alignGreenSlice(self, volume):
        green = slicer.mrmlScene.GetNodeByID('vtkMRMLSliceNodeGreen')
        green.RotateToVolumePlane(volume)
        slicer.app.layoutManager().sliceWidget('Green').sliceLogic().FitSliceToAll()

    def extractGreyImage(self, volumeIn, volumeOut):
        vol_out = slicer.util.arrayFromVolume(volumeIn)
        idx = (abs(vol_out[:, :, :, 1] - vol_out[:, :, :, 2]) < 15) * (
                    abs(vol_out[:, :, :, 1] - vol_out[:, :, :, 0]) < 15)
        vol_out_grey = np.zeros(vol_out.shape, 'uint8')
        vol_out_grey[idx] = vol_out[idx]
        vol_out_grey = vol_out_grey[:, :, :, 0]
        slicer.util.updateVolumeFromArray(volumeOut, vol_out_grey)
        volumeOut.SetAndObserveTransformNodeID(volumeIn.GetTransformNodeID())
        volumeOut.SetSpacing(volumeIn.GetSpacing())



class BiplaneRegistrationTest(ScriptedLoadableModuleTest):
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
        self.test_BiplaneRegistration1()

    def test_BiplaneRegistration1(self):
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
        logic = BiplaneRegistrationLogic()
        self.assertIsNotNone( logic.hasImageData(volumeNode) )
        self.delayDisplay('Test passed!')
