cmake_minimum_required(VERSION 3.16.3)
# project(GrowCut)

# set(GrowCut_LIBRARIES GrowCut)

# include /home/snip/Source/ITK/CMake for ITKModuleExternal
# /home/snip/Source/ITK/CMake


# include_directories(/home/snip/Source/ITK/CMake)

# if(NOT ITK_SOURCE_DIR)
#   find_package(ITK REQUIRED)
#   list(APPEND CMAKE_MODULE_PATH ${ITK_CMAKE_DIR})
#   include(ITKModuleExternal)
# else()
#   set(ITK_DIR ${CMAKE_BINARY_DIR})
#   itk_module_impl()
# endif()



project(FastGrowCut LANGUAGES CXX)

FIND_PACKAGE(ITK REQUIRED)

IF(ITK_FOUND)
  INCLUDE(${ITK_USE_FILE})
ELSE(ITK_FOUND)
  MESSAGE(FATAL_ERROR "ITK not found. Please set ITK_DIR.")
ENDIF(ITK_FOUND)

include(${ITK_USE_FILE})

SET(SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)
SET(INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include)

# INCLUDE_DIRECTORIES(/home/snip/Source/ITK/build-5.4.0/Modules/Core/Common)
# INCLUDE_DIRECTORIES(/home/snip/Source/ITK/Modules/Core/Common/include)
# INCLUDE_DIRECTORIES(/home/snip/Source/ITK/Modules/Core/SpatialObjects/include)
# INCLUDE_DIRECTORIES(/home/snip/Source/ITK/Modules/Filtering/ImageAdaptors/include)
# INCLUDE_DIRECTORIES(/home/snip/Source/ITK/Modules/Filtering/ImageFunction/include)
# INCLUDE_DIRECTORIES(/home/snip/Source/ITK/Modules/Filtering/ImageGrid/include)
# INCLUDE_DIRECTORIES(/home/snip/Source/ITK/Modules/Filtering/ImageIntensity/include)
# INCLUDE_DIRECTORIES(/home/snip/Source/ITK/Modules/Filtering/ImageStatistics/include)
# INCLUDE_DIRECTORIES(/home/snip/Source/ITK/Modules/IO/ImageBase/include)
# INCLUDE_DIRECTORIES(/home/snip/Source/ITK/Modules/IO/ImageFunction/include)
# INCLUDE_DIRECTORIES(/home/snip/Source/ITK/Modules/IO/ImageGrid/include)
# INCLUDE_DIRECTORIES(/home/snip/Source/ITK/Modules/IO/ImageIntensity/include)
# INCLUDE_DIRECTORIES(/home/snip/Source/ITK/Modules/IO/ImageStatistics/include)
# INCLUDE_DIRECTORIES(/home/snip/Source/ITK/Modules/IO/ImageAdaptor/include)
# vnl
INCLUDE_DIRECTORIES(/home/snip/Source/ITK/Modules/ThirdParty/VNL/src/vxl/v3p/vnl)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/include)

ADD_LIBRARY(FibHeap SHARED ${SOURCE_DIR}/FibHeap.cxx ${INCLUDE_DIR}/FibHeap.h)
# link
TARGET_LINK_LIBRARIES(FibHeap ${ITK_LIBRARIES})

ADD_LIBRARY(FastGrowCut SHARED ${SOURCE_DIR}/itkFastGrowCut.cxx ${INCLUDE_DIR}/itkFastGrowCut.h)
TARGET_LINK_LIBRARIES(FastGrowCut ${ITK_LIBRARIES})

ADD_EXECUTABLE(test_basic ${SOURCE_DIR}/test_basic.cxx)
TARGET_LINK_LIBRARIES(test_basic FibHeap FastGrowCut ${ITK_LIBRARIES})


