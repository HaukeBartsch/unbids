cmake_minimum_required(VERSION 3.20.1)

if(COMMAND CMAKE_POLICY)
  cmake_policy(SET CMP0003 NEW)
endif()

project(PR2MASK)

if (EXISTS "/Users/hauke/src/itk/InsightToolkit-5.3.0/bin")
  set(ITK_DIR /Users/hauke/src/itk/InsightToolkit-5.3.0/bin)
else()
  set(ITK_DIR /opt/InsightToolkit-5.3.0/bin)
endif()

#set(CMAKE_CXX_STANDARD 20)
#set(Qt5_DIR /usr/local/Cellar/qt@5/5.15.10/lib/cmake/Qt5)
#SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -framework CoreFoundation")

# find the correct path for isysroot with `xcrun --show-sdk-path`
IF("${CMAKE_SYSTEM}" MATCHES "Linux")
   if (EXISTS "/opt/InsightToolkit-5.3.0/bin")
      set(ITK_DIR /opt/InsightToolkit-5.3.0/bin)
   endif()

   set (CMAKE_CXX_STANDARD 17)
   SET(CMAKE_CXX_FLAGS "-std=c++17 -I/usr/local/include -I/usr/include")
ELSEIF("${CMAKE_HOST_SYSTEM_PROCESSOR}" MATCHES "arm64")
   # InsightToolkit-5.3.0			
   if (EXISTS "../InsightToolkit-5.3.0/bin")
      set(ITK_DIR ..//InsightToolkit-5.3.0/bin)
   endif()
   message(STATUS CMAKE_HOST_SYSTEM_PROCESSOR = ${CMAKE_HOST_SYSTEM_PROCESSOR})
   set (CMAKE_CXX_STANDARD 17)
ELSE()
   SET(CMAKE_CXX_FLAGS "-std=c++17 -isysroot /Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX.sdk -I/usr/local/include")
ENDIF()

IF(UNIX AND NOT APPLE)
  set(LINUX TRUE)
  message(STATUS ">>> Linux")
ELSE()
  message(STATUS ">>> MacOS")
ENDIF()


FIND_PACKAGE(ITK)

IF(ITK_FOUND)
  INCLUDE(${ITK_USE_FILE})
ELSE(ITK_FOUND)
  MESSAGE(FATAL_ERROR "Cannot build ITKApps without ITK.  Please set ITK_DIR.")
ENDIF(ITK_FOUND)

# Boost
set(boost_min_ver 1.50.0)
set(boost_libs system filesystem timer date_time program_options)
find_package(Boost ${boost_min_ver})

if(Boost_FOUND)
  find_package(Boost ${boost_min_ver} COMPONENTS ${boost_libs})
endif()

find_package(TBB REQUIRED)

if(TBB_FOUND)
  MESSAGE(STATUS ">>> Found TBB")
ENDIF()

# Freetype
find_package(Freetype REQUIRED)

set( GDCM_LIBRARIES
    gdcmCommon gdcmDICT gdcmDSED gdcmIOD gdcmMEXD gdcmMSFF gdcmcharls
    gdcmjpeg12 gdcmjpeg16 gdcmjpeg8 ITKIOGDCM
)

if (EXISTS "../dcmtk")
  set(DCMTK_DIR ../dcmtk)
  MESSAGE(STATUS ">>> Found local DCMTK")
endif()

find_package(DCMTK REQUIRED)

find_package(LIBZIP REQUIRED)
if (LIBZIP_FOUND)
   MESSAGE(STATUS ">>> Found LIBZIP")
ELSE(LIBZIP_FOUND)
  MESSAGE(FATAL_ERROR "Cannot build app without libzip.")
ENDIF()


# target_link_libraries(your_target ${Boost_LIBRARIES})
add_executable(uglify uglify.cxx)
target_include_directories(uglify PUBLIC ${FREETYPE_INCLUDE_DIRS})
target_link_libraries(uglify ${ITK_LIBRARIES} ${Boost_LIBRARIES} TBB::tbb ${FREETYPE_LIBRARIES})
