cmake_minimum_required(VERSION 3.20.1)

if(COMMAND CMAKE_POLICY)
  cmake_policy(SET CMP0003 NEW)
endif()

project(UNBIDS)

set (CMAKE_CXX_STANDARD 17)

set(itk_min_ver 5.3.0)
find_package(ITK ${itk_min_ver})
FIND_PACKAGE(ITK REQUIRED)
INCLUDE(${ITK_USE_FILE})
if(ITK_FOUND)
  MESSAGE(STATUS ">>> Found ITK")
ENDIF()


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

#set(ZLIB_USE_STATIC_LIBS "ON")
SET(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
find_package(ZLIB REQUIRED)
if(ZLIB_FOUND)
  MESSAGE(STATUS ">>> Found ZLIB")
ENDIF()

add_executable(uglify uglify.cxx)
target_include_directories(uglify PUBLIC)
target_link_libraries(uglify ${ITK_LIBRARIES} ${Boost_LIBRARIES} TBB::tbb)
