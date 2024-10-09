# This example script demonstrates how a grow-cut
# operation can be performed without graphical user interface.

# The first half of the script downloads input data and creates seed segments.
# The second half of the script converts segments to merged labelmap (that's the required
# input format for grow-cut filter), computes the complete segmentation, and writes
# results into new segments.

# Generate input data
################################################

import vtkSegmentationCorePython as vtkSegmentationCore 
import vtkSlicerSegmentationsModuleLogicPython as vtkSlicerSegmentationsModuleLogic
import SampleData

# Load master volume
sampleDataLogic = SampleData.SampleDataLogic()
masterVolumeNode = sampleDataLogic.downloadMRBrainTumor1()

# Create segmentation
segmentationNode = slicer.vtkMRMLSegmentationNode()
slicer.mrmlScene.AddNode(segmentationNode)
segmentationNode.CreateDefaultDisplayNodes() # only needed for display
segmentationNode.SetReferenceImageGeometryParameterFromVolumeNode(masterVolumeNode)

# Create seed segment inside tumor
tumorSeed = vtk.vtkSphereSource()
tumorSeed.SetCenter(-6, 30, 28)
tumorSeed.SetRadius(10)
tumorSeed.Update()
segmentationNode.AddSegmentFromClosedSurfaceRepresentation(tumorSeed.GetOutput(), "Diamond", [1.0,0.0,0.0])

# Create seed segment outside tumor
backgroundSeedPositions = [[0,65,32], [1, -14, 30], [0, 28, -7], [0,30,64], [31, 33, 27], [-42, 30, 27]]
append = vtk.vtkAppendPolyData()
for backgroundSeedPosition in backgroundSeedPositions:
  backgroundSeed = vtk.vtkSphereSource()
  backgroundSeed.SetCenter(backgroundSeedPosition)
  backgroundSeed.SetRadius(10)
  backgroundSeed.Update()
  append.AddInputData(backgroundSeed.GetOutput())

append.Update()
backgroundSegmentId = segmentationNode.AddSegmentFromClosedSurfaceRepresentation(append.GetOutput(), "Background", [0.0,1.0,0.0])

# Perform grow-cut
################################################

# Get segment IDs to be processed (we will process all)
selectedSegmentIds = vtk.vtkStringArray()
segmentationNode.GetSegmentation().GetSegmentIDs(selectedSegmentIds)

# Get merged labelmap extent
mergedLabelmapGeometryImage = vtkSegmentationCore.vtkOrientedImageData()
segmentationNode.GetSegmentation().SetImageGeometryFromCommonLabelmapGeometry(mergedLabelmapGeometryImage, selectedSegmentIds)
labelsEffectiveExtent = mergedLabelmapGeometryImage.GetExtent()

# Compute extent that will be passed to the algorithm (labels effective extent slightly expanded).
# We assume that the master volume has the same origin, spacing, etc. as the segmentation.
masterImageData =  vtkSlicerSegmentationsModuleLogic.vtkSlicerSegmentationsModuleLogic.CreateOrientedImageDataFromVolumeNode(masterVolumeNode)
masterImageData.UnRegister(None)
masterImageExtent = masterImageData.GetExtent()
margin = [17, 17, 17]
labelsExpandedExtent = [
  max(masterImageExtent[0], labelsEffectiveExtent[0]-margin[0]),
  min(masterImageExtent[1], labelsEffectiveExtent[1]+margin[0]),
  max(masterImageExtent[2], labelsEffectiveExtent[2]-margin[1]),
  min(masterImageExtent[3], labelsEffectiveExtent[3]+margin[1]),
  max(masterImageExtent[4], labelsEffectiveExtent[4]-margin[2]),
  min(masterImageExtent[5], labelsEffectiveExtent[5]+margin[2]) ]
mergedLabelmapGeometryImage.SetExtent(labelsExpandedExtent)

# Clip master image data to a smaller extent to reduce computation time
masterImageClipper = vtk.vtkImageConstantPad()
masterImageClipper.SetInputData(masterImageData)
masterImageClipper.SetOutputWholeExtent(mergedLabelmapGeometryImage.GetExtent())
masterImageClipper.Update()

# Create merged labelmap input for growcut
mergedImage = vtkSegmentationCore.vtkOrientedImageData()
segmentationNode.GenerateMergedLabelmapForAllSegments(mergedImage,
  vtkSegmentationCore.vtkSegmentation.EXTENT_UNION_OF_EFFECTIVE_SEGMENTS, mergedLabelmapGeometryImage, selectedSegmentIds)

# Perform grow-cut segmentation
growCutFilter = vtkSlicerSegmentationsModuleLogic.vtkImageGrowCutSegment()
growCutFilter.SetIntensityVolume(masterImageClipper.GetOutput())
growCutFilter.SetSeedLabelVolume(mergedImage)
growCutFilter.Update()
# Convert to oriented image data
resultImage = vtkSegmentationCore.vtkOrientedImageData()
resultImage.ShallowCopy(growCutFilter.GetOutput())
resultImage.CopyDirections(mergedLabelmapGeometryImage)

# Update segmentation from grow-cut result
slicer.vtkSlicerSegmentationsModuleLogic.ImportLabelmapToSegmentationNode(resultImage, segmentationNode, selectedSegmentIds)

# Delete the background segment to make the tumor visible
# segmentationNode.RemoveSegment(backgroundSegmentId)