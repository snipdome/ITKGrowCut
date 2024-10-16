

import vtkSegmentationCorePython as vtkSegmentationCore 
import vtkSlicerSegmentationsModuleLogicPython as vtkSlicerSegmentationsModuleLogic
import vtkITK
import numpy as np
import time

# Parameters
z_diamond_init_factor = 0.5


def kernel(volume):
  volume = np.swapaxes(volume, 0, 2)
  print(f"Volume shape: {volume.shape}")

  # print voxel size
  voxel_size = masterVolumeNode.GetImageData().GetSpacing()
  print(f"Voxel size: {voxel_size}")

  sub_vols = 8
  sub_vols_size = [int(x/sub_vols) for x in volume.shape]
  sub_vols_nums = [sub_vols//2-1, sub_vols//2+1]
  z_diamond_center = np.array([volume.shape[0]//2, volume.shape[1]//2, int(volume.shape[2]*z_diamond_init_factor)], dtype=int)
  z_diamond_center_lims = [z_diamond_center[2]-sub_vols_size[2], z_diamond_center[2]+sub_vols_size[2]]
  z_diamond_center_lims = [max(0, z_diamond_center_lims[0]), min(volume.shape[2], z_diamond_center_lims[1])]
  print(f"Diamond center subvolume is in:")
  print(f"X: {sub_vols_nums[0]*sub_vols_size[0]}:{sub_vols_nums[1]*sub_vols_size[0]}")
  print(f"Y: {sub_vols_nums[0]*sub_vols_size[1]}:{sub_vols_nums[1]*sub_vols_size[1]}")
  print(f"Z: {z_diamond_center_lims[0]}:{z_diamond_center_lims[1]}")
  diamond_central_subvolume = volume[
      sub_vols_nums[0]*sub_vols_size[0]:sub_vols_nums[1]*sub_vols_size[0],
      sub_vols_nums[0]*sub_vols_size[1]:sub_vols_nums[1]*sub_vols_size[1],
      z_diamond_center_lims[0]:z_diamond_center_lims[1]
  ]
  diamond_value = np.median( diamond_central_subvolume )
  print(f"Diamond value: {diamond_value}")
  del diamond_central_subvolume


  # get the bounding box of the diamond 
  # by projecting the maximum along the three axes
  # volume between 0.8 and 1.2  of the diamond value is considered diamond
  bool_volume = (volume > 0.8*diamond_value) & (volume < 1.2*diamond_value)
  x_proj = np.max(bool_volume, axis=0).astype(int)
  y_proj = np.max(bool_volume, axis=1).astype(int)
  z_proj = np.max(bool_volume, axis=2).astype(int)

  from skimage.segmentation import flood, flood_fill
  from matplotlib import pyplot as plt
  def detect_bounding_box(proj, seedpoint,savepath=None): # proj is a 2D numpy array
    if savepath is not None:
      plt.imshow(proj)
      plt.savefig(savepath)
    plt.close()
    #proj2 = flood_fill(proj, seedpoint, 2)
    proj2 = flood(proj, seedpoint, connectivity=1)
    # if savepath is not None:
    #   plt.imshow(proj2)
    #   plt.savefig(savepath.replace(".png", "_2.png"))
    # get the four corners of the bounding box where the projection is 2
    #bool_proj2 = proj2-proj
    bool_proj2= proj2
    x_min = np.where(bool_proj2.sum(axis=0))[0][0]
    x_max = np.where(bool_proj2.sum(axis=0))[0][-1]
    y_min = np.where(bool_proj2.sum(axis=1))[0][0]
    y_max = np.where(bool_proj2.sum(axis=1))[0][-1]
    if savepath is not None:
      plt.imshow(proj2)
      # draw the bounding box
      plt.plot([x_min, x_max, x_max, x_min, x_min], [y_min, y_min, y_max, y_max, y_min], 'r')
      plt.savefig(savepath.replace(".png", "_2.png"))
    plt.close()
    return x_min, x_max, y_min, y_max # in matplotlib, the first two coordinates goes along the horizontal axis, the second two along the vertical axis
    
  seedpoint_x = (z_diamond_center[1], z_diamond_center[2])
  z_min, z_max, y_min, y_max = detect_bounding_box(x_proj, seedpoint_x, "proj_x.png")
  print(f'Seedpoint: {seedpoint_x}')
  print(f'Diamond bounding box (x): {y_min}:{y_max}, {z_min}:{z_max}')
  seedpoint_y = (z_diamond_center[0], z_diamond_center[2])
  z_min, z_max, x_min, x_max = detect_bounding_box(y_proj, seedpoint_y, "proj_y.png") 
  print(f'Seedpoint: {seedpoint_y}')
  print(f'Diamond bounding box (y): {x_min}:{x_max}, {z_min}:{z_max}')
  seedpoint_z = (z_diamond_center[0], z_diamond_center[1])
  y_min, y_max, x_min, x_max = detect_bounding_box(z_proj, seedpoint_z, "proj_z.png")
  print(f'Seedpoint: {seedpoint_z}')
  print(f'Diamond bounding box (z): {x_min}:{x_max}, {y_min}:{y_max}')


  print(f"Diamond bounding box:")
  print(f"X: {x_min}:{x_max}")
  print(f"Y: {y_min}:{y_max}")
  print(f"Z: {z_min}:{z_max}")

  del x_proj, y_proj, z_proj, bool_volume

  new_shape = [x_max-x_min, y_max-y_min, z_max-z_min]


  extent = masterVolumeNode.GetImageData().GetExtent()
  # print(f"Extent:")
  # print(f"X: {extent[0]}:{extent[1]}")
  # print(f"Y: {extent[2]}:{extent[3]}")
  # print(f"Z: {extent[4]}:{extent[5]}")
  # set origin
  origin = masterVolumeNode.GetImageData().GetOrigin()
  origin = [0, 0, 0]
  masterVolumeNode.GetImageData().SetOrigin(origin)

  # Create segmentation
  segmentationNode = slicer.vtkMRMLSegmentationNode()
  slicer.mrmlScene.AddNode(segmentationNode)
  segmentationNode.CreateDefaultDisplayNodes() # only needed for display
  segmentationNode.SetReferenceImageGeometryParameterFromVolumeNode(masterVolumeNode)


  # Create seed segment inside diamond
  diamond_seed = vtk.vtkSphereSource()
  diamond_seed.SetCenter(-extent[1]/2, -extent[3]/2, extent[5]*z_diamond_init_factor)
  diamond_seed.SetRadius(10)
  diamond_seed.Update()
  segmentationNode.AddSegmentFromClosedSurfaceRepresentation(diamond_seed.GetOutput(), "Diamond", [1.0,0.0,0.0])

  # Create seed segment outside diamond
  print(f"Background seed positions:")
  print(f"X: {x_min}:{x_max}")
  print(f"Y: {y_min}:{y_max}")
  print(f"Z: {z_min}:{z_max}")
  backgroundSeedPositions = [
          [-x_min, -y_min, z_min],
          [-x_max, -y_min, z_min],
          [-x_min, -y_max, z_min],
          [-x_max, -y_max, z_min],
          [-x_min, -y_min, z_max],
          [-x_max, -y_min, z_max],
          [-x_min, -y_max, z_max],
          [-x_max, -y_max, z_max]
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
  margin = [6, 6, 6]
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
  fastgrowcut = vtkITK.vtkITKGrowCut()
  fastgrowcut.SetIntensityVolume(masterImageClipper.GetOutput())
  fastgrowcut.SetSeedLabelVolume(mergedImage)
  time_fastgrowcut = time.time()
  fastgrowcut.Update()
  time_fastgrowcut = time.time() - time_fastgrowcut

  # Convert to oriented image data
  resultImage = vtkSegmentationCore.vtkOrientedImageData()
  ###resultImage.ShallowCopy(mergedLabelmapGeometryImage.GetOutput())
  resultImage.ShallowCopy(fastgrowcut.GetOutput())
  resultImage.CopyDirections(mergedLabelmapGeometryImage)

  # Update segmentation from grow-cut result
  slicer.vtkSlicerSegmentationsModuleLogic.ImportLabelmapToSegmentationNode(resultImage, segmentationNode, selectedSegmentIds)

  # print timing
  print(f"FastGrowCut time: {time_fastgrowcut:.2f} s")



if __name__ == "__main__":

  # Create master volume
  #masterVolumeNode = slicer.util.loadVolume("D:\datasets\diamcad\RSD0711-1B.tif")
  masterVolumeNode = slicer.util.loadVolume("D:\datasets\diamcad\SFM0822-1.tif")

  # get numpy array from volume
  volume = slicer.util.arrayFromVolume(masterVolumeNode)

  kernel(volume)