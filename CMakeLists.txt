
cmake_minimum_required(VERSION 2.8)

project(FourierImage)

find_package(ITK REQUIRED)
include(${ITK_USE_FILE})
if (ITKVtkGlue_LOADED)
  find_package(VTK REQUIRED)
  include(${VTK_USE_FILE})
endif()


## Executables



ADD_EXECUTABLE(ImageRegistration1.cxx )
TARGET_LINK_LIBRARIES(ImageRegistration1 ${ITK_LIBRARIES})
  
ADD_EXECUTABLE(ImageRegistration2.cxx )
TARGET_LINK_LIBRARIES(ImageRegistration2 ${ITK_LIBRARIES})
  
ADD_EXECUTABLE(ImageRegistration3.cxx )
TARGET_LINK_LIBRARIES(ImageRegistration3 ${ITK_LIBRARIES})





















