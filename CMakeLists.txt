cmake_minimum_required(VERSION 2.8)

project(SOAX3D)

find_package(ITK REQUIRED)
if(ITK_FOUND)
  include(${ITK_USE_FILE})
endif(ITK_FOUND)

find_package(VTK 6.0 REQUIRED NO_MODULE)
if(VTK_FOUND)
  include(${VTK_USE_FILE})
endif(VTK_FOUND)

find_package(Qt4 REQUIRED)
if(QT_FOUND)
  include(${QT_USE_FILE})
endif(QT_FOUND)

set(Boost_USE_STATIC_LIBS        ON)
set(Boost_USE_MULTITHREADED      ON)
set(Boost_USE_STATIC_RUNTIME    OFF)

find_package(Boost 1.54.0 COMPONENTS program_options REQUIRED)

if(Boost_FOUND)
  include_directories(${Boost_INCLUDE_DIRS})
endif(Boost_FOUND)

set(gui_srcs
  gui.cc
  # viewer.cc
  main_window.cc
  # parameters_dialog.cc
  # view_options_dialog.cc
  # analysis_options_dialog.cc
  )

set(common_srcs
  global.h
  )

set(moc_gui_srcs
  main_window.h
  # viewer.h
  # parameters_dialog.h
  # view_options_dialog.h
  # analysis_options_dialog.h
  )

set(rcc_srcs icons.qrc)

add_definitions(-Wall)

qt4_wrap_cpp(mocs ${moc_gui_srcs})
qt4_add_resources(rccs ${rcc_srcs})

add_executable(soax3D ${common_srcs} ${gui_srcs} ${mocs} ${rccs})

target_link_libraries(soax3D
  ${QT_LIBRARIES}
  ${VTK_LIBRARIES}
  ${ITK_LIBRARIES}
  ${Boost_PROGRAM_OPTIONS_LIBRARY}
  )
