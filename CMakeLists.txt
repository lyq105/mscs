cmake_minimum_required(VERSION 2.8)

if(POLICY CMP0020)
	cmake_policy(SET CMP0020 NEW)
endif()

PROJECT(MSCS)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

find_package(Qt4 REQUIRED)
include(${QT_USE_FILE})

#set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)
#set(CMAKE_AUTOUIC ON)

include_directories(${CMAKE_CURRENT_SOURCE_DIR} 
										${CMAKE_CURRENT_BINARY_DIR} 
										src/elastic/ 
										src/gui/
										)

file(GLOB UI_FILES forms/*.ui)
file(GLOB QT_WRAP src/gui/*.h)
file(GLOB CXX_FILES src/*.cpp
										src/gui/*.cpp 
										src/elastic/*.cpp)
#file(GLOB RCC_FILES resources/*.qrc)
set(RCC_FILES Mscs.qrc)

qt4_wrap_ui(UISrcs ${UI_FILES})
qt4_wrap_cpp(MOCSrcs ${QT_WRAP})
#qt4_add_resources(QRCFILES ${RCC_FILE})

add_executable(MSCS 
	#WIN32
 	${CXX_FILES} 
	${UISrcs} 
	${MOCSrcs} 
	Mscs.qrc)

target_link_libraries ( MSCS
  ${VTK_LIBRARIES}
	${QT_LIBRARIES}
  QVTK 
  )

IF(WIN32) # Check if we are on Windows
	SET(QT_USE_QTMAIN TRUE)
  if(MSVC) # Check if we are using the Visual Studio compiler
		set_target_properties(MSCS PROPERTIES LINK_FLAGS_RELEASE "/SUBSYSTEM:WINDOWS")
  elseif(CMAKE_COMPILER_IS_GNUCXX)
    # SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mwindows") # Not tested
  else()
    message(SEND_ERROR "You are using an unsupported Windows compiler! (Not MSVC or GCC)")
  endif()
elseif(UNIX)
  # Nothing special required
else()
  message(SEND_ERROR "You are on an unsupported platform! (Not Win32 or Unix)")
ENDIF()
