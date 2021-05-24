
cmake_minimum_required(VERSION 2.8)

project(FourierImage)

find_package(ITK REQUIRED)
include(${ITK_USE_FILE})
if (ITKVtkGlue_LOADED)
  find_package(VTK REQUIRED)
  include(${VTK_USE_FILE})
endif()

add_executable(Versor3dResgistration Versor3dResgistration.cxx)
target_link_libraries(Versor3dResgistration ${ITK_LIBRARIES})

add_executable(VersorScale3dResgistration VersorScale3dResgistration.cxx)
target_link_libraries(VersorScale3dResgistration ${ITK_LIBRARIES})

add_executable(AffineRegistration AffineRegistration.cxx)
target_link_libraries(AffineRegistration ${ITK_LIBRARIES})