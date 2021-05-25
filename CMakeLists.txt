
cmake_minimum_required(VERSION 2.8)

project(FourierImage)

find_package(ITK REQUIRED)
include(${ITK_USE_FILE})
if (ITKVtkGlue_LOADED)
  find_package(VTK REQUIRED)
  include(${VTK_USE_FILE})
endif()

add_executable(Versor3dRegistration Versor3dRegistration.cxx)
target_link_libraries(Versor3dRegistration ${ITK_LIBRARIES})

add_executable(VersorScale3dRegistration VersorScale3dRegistration.cxx)
target_link_libraries(VersorScale3dRegistration ${ITK_LIBRARIES})

add_executable(AffineRegistration AffineRegistration.cxx)
target_link_libraries(AffineRegistration ${ITK_LIBRARIES})