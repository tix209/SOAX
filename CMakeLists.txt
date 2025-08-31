cmake_minimum_required(VERSION 3.16)

project(SOAX)

# Ensure all targets use position independent code for Qt compatibility
add_compile_options(-fPIC)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Compiler-specific optimizations
if(CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /O2 /W3")
    # Disable some verbose warnings
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /wd4996")
elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2 -Wall -fPIC")
endif()

# --- Static linking for ITK, Boost, VTK ---
set(BUILD_SHARED_LIBS OFF)
set(ITK_STATIC ON)
set(Boost_USE_STATIC_LIBS ON)
# VTK static linking (if available)
set(VTK_BUILD_SHARED_LIBS OFF)

# Find packages
# ITK with VTK integration
find_package(ITK REQUIRED COMPONENTS
  ITKCommon
  ITKIOImageBase
  ITKIOPNG
  ITKIOMeta
  ITKMetaIO
  ITKIOJPEG
  ITKIOTIFF
  ITKIOBMP
  ITKImageIntensity
  ITKImageFunction
  ITKVtkGlue
  ITKFEM
)
include(${ITK_USE_FILE})

# ITK version check - allow ITK 4.x for legacy compatibility
if(ITK_VERSION_MAJOR LESS 4)
    message(FATAL_ERROR "ITK 4.0 or later is required")
endif()

# Find VTK first (basic find to get version)
find_package(VTK REQUIRED)

# Check VTK version requirement
if(VTK_VERSION VERSION_LESS "8.2")
    message(FATAL_ERROR "VTK 8.2 or later is required. Found VTK ${VTK_VERSION}")
endif()

# Now find VTK with specific components based on version
if(VTK_VERSION VERSION_GREATER_EQUAL "9.0")
    # VTK 9.x component names (no 'vtk' prefix)
    set(VTK_COMPONENTS
      CommonCore
      CommonDataModel
      CommonExecutionModel
      FiltersCore
      FiltersGeneral
      FiltersGeometry
      FiltersSources
      ImagingCore
      InteractionStyle
      InteractionWidgets
      IOImage
      RenderingCore
      RenderingOpenGL2
      RenderingVolume
      RenderingVolumeOpenGL2
      GUISupportQt
    )
else()
    # VTK 8.x component names (with 'vtk' prefix)
    set(VTK_COMPONENTS
      vtkCommonCore
      vtkCommonDataModel
      vtkCommonExecutionModel
      vtkFiltersCore
      vtkFiltersGeneral
      vtkFiltersGeometry
      vtkFiltersSources
      vtkImagingCore
      vtkInteractionStyle
      vtkInteractionWidgets
      vtkIOImage
      vtkRenderingCore
      vtkRenderingOpenGL2
      vtkRenderingVolume
      vtkRenderingVolumeOpenGL2
      vtkGUISupportQt
    )
endif()

# Find VTK with the determined components
find_package(VTK REQUIRED COMPONENTS ${VTK_COMPONENTS})

# Find Qt5 or Qt6
find_package(Qt6 QUIET COMPONENTS Core Widgets OpenGL)
if(NOT Qt6_FOUND)
    find_package(Qt5 REQUIRED COMPONENTS Core Widgets OpenGL X11Extras)
    set(Qt_VERSION_MAJOR 5)
else()
    set(Qt_VERSION_MAJOR 6)
endif()
# Qt will be linked dynamically (shared)

# Enable automatic MOC, UIC, RCC
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)
set(CMAKE_AUTOUIC ON)

# Add Qt include directories
if(Qt_VERSION_MAJOR EQUAL 6)
    include_directories(${Qt6Widgets_INCLUDE_DIRS} ${Qt6Core_INCLUDE_DIRS} ${Qt6OpenGL_INCLUDE_DIRS})
else()
    include_directories(${Qt5Widgets_INCLUDE_DIRS} ${Qt5Core_INCLUDE_DIRS} ${Qt5OpenGL_INCLUDE_DIRS} ${Qt5X11Extras_INCLUDE_DIRS})
endif()

# Add VTK include directories
include_directories(${VTK_INCLUDE_DIRS})

## Enable these options for Windows build
if (WIN32)
  set(Boost_USE_STATIC_LIBS        ON)
  set(Boost_USE_MULTITHREADED      ON)
  set(Boost_USE_STATIC_RUNTIME    OFF)
endif()

find_package(Boost COMPONENTS
  program_options filesystem system REQUIRED)

if(Boost_FOUND)
  include_directories(${Boost_INCLUDE_DIRS})
endif(Boost_FOUND)

set(gui_srcs
  viewer.cc
  main_window.cc
  parameters_dialog.cc
  view_options_dialog.cc
  analysis_options_dialog.cc
  )

set(common_srcs
  global.h
  snake.h
  snake.cc
  solver_bank.h
  solver_bank.cc
  junctions.h
  junctions.cc
  snake_tip.h
  snake_tip.cc
  snake_tip_set.cc
  snake_tip_set.h
  utility.h
  utility.cc
  )

set(moc_gui_srcs
  main_window.h
  viewer.h
  parameters_dialog.h
  view_options_dialog.h
  analysis_options_dialog.h
  )

set(rcc_srcs icons.qrc)

# Qt resource compilation and MOC files are handled automatically with AUTOMOC/AUTORCC

if (APPLE)
  set(ICON_FOLDER icon)
  set(ICON_NAME letter-x.icns)
  set(ICON_FILE ${ICON_FOLDER}/${ICON_NAME})
  set_source_files_properties(${ICON_FILE} PROPERTIES
    MACOSX_PACKAGE_LOCATION Resources)

  set(COPYRIGHT "Copyright (c) 2015 Lehigh University. All rights reserved.")
  set(IDENTIFIER "edu.lehigh.soax")

  set(MACOSX_BUNDLE_BUNDLE_NAME ${PROJECT_NAME})
  set(MACOSX_BUNDLE_BUNDLE_VERSION ${PROJECT_VERSION})
  set(MACOSX_BUNDLE_LONG_VERSION_STRING ${PROJECT_VERSION})
  set(MACOSX_BUNDLE_COPYRIGHT ${COPYRIGHT})
  set(MACOSX_BUNDLE_GUI_IDENTIFIER ${IDENTIFIER})
  set(MACOSX_BUNDLE_ICON_FILE ${ICON_NAME})
  add_executable(soax MACOSX_BUNDLE soax.cc ${common_srcs} ${gui_srcs} 
    ${rcc_srcs} multisnake.cc ${ICON_FILE})
else ()
  add_executable(soax soax.cc ${common_srcs} ${gui_srcs} 
    ${rcc_srcs} multisnake.cc)
endif ()

add_executable(batch_soax batch_soax.cc ${common_srcs} multisnake.cc)
set_target_properties(batch_soax PROPERTIES AUTOMOC ON)
add_executable(best_snake best_snake.cc ${common_srcs} multisnake.cc)
set_target_properties(best_snake PROPERTIES AUTOMOC ON)
add_executable(batch_length batch_length.cc ${common_srcs} multisnake.cc)
set_target_properties(batch_length PROPERTIES AUTOMOC ON)
add_executable(batch_resample batch_resample.cc)

# Disable AUTOMOC for batch targets that don't use Qt
set_target_properties(batch_resample PROPERTIES AUTOMOC OFF)

# Link libraries
if(Qt_VERSION_MAJOR EQUAL 6)
    set(QT_LIBRARIES Qt6::Core Qt6::Widgets Qt6::OpenGL)
    include_directories(${Qt6Widgets_INCLUDE_DIRS} ${Qt6Core_INCLUDE_DIRS} ${Qt6OpenGL_INCLUDE_DIRS})
else()
    set(QT_LIBRARIES Qt5::Core Qt5::Widgets Qt5::OpenGL Qt5::X11Extras)
    include_directories(${Qt5Widgets_INCLUDE_DIRS} ${Qt5Core_INCLUDE_DIRS} ${Qt5OpenGL_INCLUDE_DIRS} ${Qt5X11Extras_INCLUDE_DIRS})
endif()

# Ensure ITK IO modules are linked for static builds
if(ITK_FOUND)
  set(ITK_IO_LIBS
    ITKIOImageBase
    ITKIOPNG
    ITKIOMeta
    ITKMetaIO
    ITKIOJPEG
    ITKIOTIFF
    ITKIOBMP
    # Add more IO modules here if needed
  )
  foreach(lib ${ITK_IO_LIBS})
    list(APPEND ITK_LIBRARIES ${${lib}_LIBRARY})
  endforeach()
endif()

# VTK static linking (if available)
if(VTK_FOUND)
  set(VTK_STATIC_LIBS ${VTK_LIBRARIES})
endif()

target_link_libraries(soax
  ${QT_LIBRARIES}
  ${VTK_STATIC_LIBS}
  ${ITK_LIBRARIES}
  )

target_link_libraries(batch_soax
  Qt5::Core
  ${ITK_LIBRARIES}
  ${Boost_PROGRAM_OPTIONS_LIBRARY}
  ${Boost_FILESYSTEM_LIBRARY}
  ${Boost_SYSTEM_LIBRARY}
  )

target_link_libraries(best_snake
  Qt5::Core
  ${ITK_LIBRARIES}
  ${Boost_PROGRAM_OPTIONS_LIBRARY}
  ${Boost_FILESYSTEM_LIBRARY}
  ${Boost_SYSTEM_LIBRARY}
  )

target_link_libraries(batch_length
  Qt5::Core
  ${ITK_LIBRARIES}
  ${Boost_PROGRAM_OPTIONS_LIBRARY}
  ${Boost_FILESYSTEM_LIBRARY}
  ${Boost_SYSTEM_LIBRARY}
  )

target_link_libraries(batch_resample
  ${ITK_LIBRARIES}
  ${Boost_PROGRAM_OPTIONS_LIBRARY}
  ${Boost_FILESYSTEM_LIBRARY}
  ${Boost_SYSTEM_LIBRARY}
  )
