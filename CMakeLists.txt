cmake_minimum_required(VERSION 2.8)

project(SOAX)

find_package(ITK 4 REQUIRED)
if(ITK_FOUND)
  include(${ITK_USE_FILE})
endif(ITK_FOUND)

find_package(VTK REQUIRED NO_MODULE)
if(VTK_FOUND)
  include(${VTK_USE_FILE})
endif(VTK_FOUND)

find_package(Qt4 REQUIRED)
if(QT_FOUND)
  include(${QT_USE_FILE})
endif(QT_FOUND)

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
#  multisnake.h
#  multisnake.cc
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

# add_definitions(-Wall)

qt4_wrap_cpp(multisnake_moc multisnake.h)
qt4_wrap_cpp(mocs ${moc_gui_srcs})
qt4_add_resources(rccs ${rcc_srcs})

add_executable(soax soax.cc ${common_srcs} ${gui_srcs} ${mocs}
  ${rccs} ${multisnake_moc} multisnake.cc multisnake.h)

add_executable(batch_soax batch_soax.cc ${common_srcs}
  ${multisnake_moc} multisnake.cc multisnake.h)
# add_executable(eval_gt evaluate_gt.cc ${common_srcs})
# add_executable(eval_f evaluate_f.cc ${common_srcs})
add_executable(best_snake best_snake.cc ${common_srcs}
  ${multisnake_moc} multisnake.cc multisnake.h)

# add_executable(snr image_snr.cc ${common_srcs})
add_executable(syn2 synthetic_tamara.cc ${common_srcs}
  ${multisnake_moc} multisnake.cc multisnake.h)
# add_executable(snapshots batch_snapshots.cc)
add_executable(batch_resample batch_resample.cc)

target_link_libraries(soax
  ${QT_LIBRARIES}
  ${VTK_LIBRARIES}
  ${ITK_LIBRARIES}
  ${Boost_PROGRAM_OPTIONS_LIBRARY}
  )

target_link_libraries(batch_soax
  ${QT_LIBRARIES}
  ${ITK_LIBRARIES}
  ${Boost_PROGRAM_OPTIONS_LIBRARY}
  ${Boost_FILESYSTEM_LIBRARY}
  ${Boost_SYSTEM_LIBRARY}
  )

# target_link_libraries(eval_gt
#   ${ITK_LIBRARIES}
#   ${Boost_PROGRAM_OPTIONS_LIBRARY}
#   ${Boost_FILESYSTEM_LIBRARY}
#   ${Boost_SYSTEM_LIBRARY}
#   )

# target_link_libraries(eval_f
#   ${ITK_LIBRARIES}
#   ${Boost_PROGRAM_OPTIONS_LIBRARY}
#   ${Boost_FILESYSTEM_LIBRARY}
#   ${Boost_SYSTEM_LIBRARY}
#   )

target_link_libraries(best_snake
  ${QT_LIBRARIES}
  ${ITK_LIBRARIES}
  ${Boost_PROGRAM_OPTIONS_LIBRARY}
  ${Boost_FILESYSTEM_LIBRARY}
  ${Boost_SYSTEM_LIBRARY}
  )

# target_link_libraries(bestt
#   ${ITK_LIBRARIES}
#   ${Boost_PROGRAM_OPTIONS_LIBRARY}
#   ${Boost_FILESYSTEM_LIBRARY}
#   ${Boost_SYSTEM_LIBRARY}
#   )

target_link_libraries(syn2 ${ITK_LIBRARIES})

target_link_libraries(batch_resample
  ${ITK_LIBRARIES}
  ${Boost_PROGRAM_OPTIONS_LIBRARY}
  ${Boost_FILESYSTEM_LIBRARY}
  ${Boost_SYSTEM_LIBRARY}
  )
