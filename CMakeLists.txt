cmake_minimum_required(VERSION 2.8)
 
if(POLICY CMP0020)
  cmake_policy(SET CMP0020 NEW)
endif()
 
PROJECT(MSCS)
 
find_package(VTK REQUIRED)
include(${VTK_USE_FILE})
 
if(${VTK_VERSION} VERSION_GREATER "6" AND VTK_QT_VERSION VERSION_GREATER "4")
  # Instruct CMake to run moc automatically when needed.
  set(CMAKE_AUTOMOC ON)
  find_package(Qt5Widgets REQUIRED QUIET)
else()
  find_package(Qt4 REQUIRED)
  include(${QT_USE_FILE})
endif()
 
include_directories(${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR} src/pmcell/ src/gui/)
 
file(GLOB UI_FILES forms/*.ui)
file(GLOB QT_WRAP src/gui/*.h /src/pmcell/*.h)
file(GLOB CXX_FILES src/gui/*.cpp src/*.cpp src/pmcell/*.cpp)
file(GLOB RCC_FILES resources/*.qrc)
 
if(${VTK_VERSION} VERSION_GREATER "6" AND VTK_QT_VERSION VERSION_GREATER "4")
  qt5_wrap_ui(UISrcs ${UI_FILES} )
  # CMAKE_AUTOMOC in ON so the MocHdrs will be automatically wrapped.
  add_executable(MSCS MACOSX_BUNDLE
		${CXX_FILES} ${UISrcs} ${QT_WRAP} ${RCC_FILES})
  qt5_use_modules(MSCS Core Gui)
	qt5_add_resources(qrc_FILES ${RCC_FILES})
  target_link_libraries(MSCS ${VTK_LIBRARIES} ${qrc_FILES})
else()
  QT4_WRAP_UI(UISrcs ${UI_FILES})
  QT4_WRAP_CPP(MOCSrcs ${QT_WRAP})
	qt4_add_resources(qrc_FILES ${RCC_FILE})
	add_executable(MSCS MACOSX_BUNDLE ${CXX_FILES} ${UISrcs} ${MOCSrcs} ${qrc_FILES})
 
  if(VTK_LIBRARIES)
    if(${VTK_VERSION} VERSION_LESS "6")
      target_link_libraries(MSCS ${VTK_LIBRARIES} QVTK)
    else()
      target_link_libraries(MSCS ${VTK_LIBRARIES})
    endif()
  else()
    target_link_libraries(MSCS vtkHybrid QVTK vtkViews ${QT_LIBRARIES})
  endif()
endif()
