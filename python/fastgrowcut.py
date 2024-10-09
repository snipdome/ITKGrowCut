

import vtkSegmentationCorePython as vtkSegmentationCore 
import vtkSlicerSegmentationsModuleLogicPython as vtkSlicerSegmentationsModuleLogic
import vtkITK


# Create master volume
masterVolumeNode = slicer.util.loadVolume("D:\datasets\diamcad\RSD0711-1B.tif")

extent = masterVolumeNode.GetImageData().GetExtent()
# set origin
origin = masterVolumeNode.GetImageData().GetOrigin()
origin = [0, 0, 0]
masterVolumeNode.GetImageData().SetOrigin(origin)
# set direction
direction = masterVolumeNode.GetImageData().GetDirectionMatrix()
print(direction)
# new_direction = [[-1, 0, 0], [0, -1, 0], [0, 0, 1]]
# for i in range(3):
#     for j in range(3):
#         direction.SetElement(i, j, new_direction[i][j])
# masterVolumeNode.GetImageData().SetDirectionMatrix(direction)

# Create segmentation
segmentationNode = slicer.vtkMRMLSegmentationNode()
slicer.mrmlScene.AddNode(segmentationNode)
segmentationNode.CreateDefaultDisplayNodes() # only needed for display
segmentationNode.SetReferenceImageGeometryParameterFromVolumeNode(masterVolumeNode)


# Create seed segment inside tumor
diamond_seed = vtk.vtkSphereSource()
diamond_seed.SetCenter(-extent[1]/2, -extent[3]/2, extent[5]/2)
diamond_seed.SetRadius(10)
diamond_seed.Update()
segmentationNode.AddSegmentFromClosedSurfaceRepresentation(diamond_seed.GetOutput(), "Diamond", [1.0,0.0,0.0])

# Create seed segment outside tumor
# corners
backgroundSeedPositions = [
        [-extent[0], -extent[2], extent[4]],
        [-extent[1], -extent[2], extent[4]],
        [-extent[0], -extent[3], extent[4]],
        [-extent[1], -extent[3], extent[4]],
        [-extent[0], -extent[2], extent[5]],
        [-extent[1], -extent[2], extent[5]],
        [-extent[0], -extent[3], extent[5]],
        [-extent[1], -extent[3], extent[5]]
    ]
append = vtk.vtkAppendPolyData()
for backgroundSeedPosition in backgroundSeedPositions:
  backgroundSeed = vtk.vtkSphereSource()
  backgroundSeed.SetCenter(backgroundSeedPosition)
  backgroundSeed.SetRadius(10)
  backgroundSeed.Update()
  append.AddInputData(backgroundSeed.GetOutput())

append.Update()
segmentationNode.AddSegmentFromClosedSurfaceRepresentation(append.GetOutput(), "Background", [0.0,1.0,0.0])

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
# growCutFilter = vtkSlicerSegmentationsModuleLogic.vtkImageGrowCutSegment()
# growCutFilter.SetIntensityVolume(masterImageClipper.GetOutput())
# growCutFilter.SetSeedLabelVolume(mergedImage)
# growCutFilter.Update()

fastgrowcut = vtkITK.vtkITKGrowCut()
fastgrowcut.SetIntensityVolume(masterImageClipper.GetOutput())
fastgrowcut.SetSeedLabelVolume(mergedImage)
fastgrowcut.Update()
growCutFilter = fastgrowcut

# Convert to oriented image data
resultImage = vtkSegmentationCore.vtkOrientedImageData()
#resultImage.ShallowCopy(mergedLabelmapGeometryImage.GetImageData())
resultImage.ShallowCopy(growCutFilter.GetOutput())
resultImage.CopyDirections(mergedLabelmapGeometryImage)

# Update segmentation from grow-cut result
slicer.vtkSlicerSegmentationsModuleLogic.ImportLabelmapToSegmentationNode(resultImage, segmentationNode, selectedSegmentIds)

