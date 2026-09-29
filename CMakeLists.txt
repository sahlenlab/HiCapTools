cmake_minimum_required(VERSION 3.10)

project(HiCapTools CXX)

option(HICAPTOOLS_USE_SYSTEM_BAMTOOLS "Prefer an installed BamTools library" ON)
option(HICAPTOOLS_USE_SYSTEM_BOOST "Prefer installed Boost headers" ON)

#find_package(OpenMP REQUIRED)
#include(FindOpenMP)
#if(OPENMP_FOUND)
#   set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
#    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
#    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
#endif()
 set(MAC FALSE)
 if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
	set(MAC TRUE)
 endif()
  
 if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    # require at least gcc 4.9
    if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 4.9)
        message(FATAL_ERROR "GCC version must be at least 4.9!")
    endif()
 endif()

file(GLOB_RECURSE SOURCES "${CMAKE_SOURCE_DIR}/external/*.cpp")
file(GLOB_RECURSE SOURCES "${CMAKE_SOURCE_DIR}/ExterAlg/*.cpp")

# set our library and executable destination dirs
  set( EXECUTABLE_OUTPUT_PATH "${CMAKE_SOURCE_DIR}/bin" )
  
  file (GLOB CPPFILES "${CMAKE_SOURCE_DIR}/src/*.cpp")
  file (GLOB HFILES "${CMAKE_SOURCE_DIR}/include/*.h") 
  file (GLOB HPPFILES "${CMAKE_SOURCE_DIR}/include/*.hpp")

  # Ignore AppleDouble metadata files that can appear after macOS transfers.
  list(FILTER CPPFILES EXCLUDE REGEX "/\\._[^/]+$")
  list(FILTER HFILES EXCLUDE REGEX "/\\._[^/]+$")
  list(FILTER HPPFILES EXCLUDE REGEX "/\\._[^/]+$")
  list(FILTER SOURCES EXCLUDE REGEX "/\\._[^/]+$")

  add_executable(HiCapTools ${CPPFILES} ${HFILES} ${HPPFILES} ${SOURCES})
   
  set_property(TARGET HiCapTools PROPERTY CXX_STANDARD 11)
  set_property(TARGET HiCapTools PROPERTY CXX_STANDARD_REQUIRED ON)

  find_package(ZLIB REQUIRED)
  add_library(libbamtools UNKNOWN IMPORTED)

  if(HICAPTOOLS_USE_SYSTEM_BAMTOOLS)
	find_path(BAMTOOLS_INCLUDE_DIR api/BamReader.h
	  HINTS ENV BAMTOOLS_ROOT
	  PATHS /opt/homebrew/include/bamtools /usr/local/include/bamtools
	  PATH_SUFFIXES include/bamtools include)
	find_library(BAMTOOLS_LIBRARY NAMES bamtools
	  HINTS ENV BAMTOOLS_ROOT
	  PATHS /opt/homebrew/lib /usr/local/lib
	  PATH_SUFFIXES lib lib64)
  endif()

  if(BAMTOOLS_INCLUDE_DIR AND BAMTOOLS_LIBRARY)
	message(STATUS "Using system BamTools: ${BAMTOOLS_LIBRARY}")
	set_target_properties(libbamtools PROPERTIES IMPORTED_LOCATION "${BAMTOOLS_LIBRARY}")
	include_directories(BEFORE "${BAMTOOLS_INCLUDE_DIR}")
  else()
	if(APPLE)
	  set(BUNDLED_BAMTOOLS_LIBRARY "${CMAKE_SOURCE_DIR}/bamtools/libbamtools.2.4.1.dylib")
	else()
	  set(BUNDLED_BAMTOOLS_LIBRARY "${CMAKE_SOURCE_DIR}/bamtools/libbamtools.so.2.4.1")
	endif()
	if(NOT EXISTS "${BUNDLED_BAMTOOLS_LIBRARY}")
	  message(FATAL_ERROR "BamTools was not found and bundled library is missing: ${BUNDLED_BAMTOOLS_LIBRARY}")
	endif()
	message(STATUS "Using bundled BamTools: ${BUNDLED_BAMTOOLS_LIBRARY}")
	set_target_properties(libbamtools PROPERTIES IMPORTED_LOCATION "${BUNDLED_BAMTOOLS_LIBRARY}")
	if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU" AND CMAKE_CXX_COMPILER_VERSION VERSION_GREATER 5.1)
	  target_compile_definitions(HiCapTools PRIVATE _GLIBCXX_USE_CXX11_ABI=0)
	endif()
	if(UNIX AND NOT APPLE)
	  set_target_properties(HiCapTools PROPERTIES BUILD_RPATH "\$ORIGIN/../bamtools")
	endif()
  endif()
  if(HICAPTOOLS_USE_SYSTEM_BOOST)
	find_path(BOOST_INCLUDE_DIR boost/accumulators/accumulators.hpp
	  HINTS ENV BOOST_HOME ENV BOOST_ROOT ENV EBROOTBOOST
	  PATHS /opt/homebrew/include /usr/local/include
	  PATH_SUFFIXES include)
  endif()

  if(BOOST_INCLUDE_DIR)
	message(STATUS "Using system Boost headers: ${BOOST_INCLUDE_DIR}")
	include_directories(BEFORE "${BOOST_INCLUDE_DIR}")
  elseif(EXISTS "${CMAKE_SOURCE_DIR}/boost/boost/accumulators/accumulators.hpp")
	message(STATUS "Using bundled Boost headers")
	include_directories("${CMAKE_SOURCE_DIR}/boost")
  else()
	message(FATAL_ERROR
	  "Boost headers were not found. On Dardel, run 'ml PDC', 'ml PrgEnv-gnu', and 'ml boost' before building.")
  endif()
  include_directories( "${CMAKE_SOURCE_DIR}/ExterAlg")
  include_directories("${CMAKE_SOURCE_DIR}/include")
  include_directories("${CMAKE_SOURCE_DIR}/external")
  target_link_libraries( HiCapTools libbamtools ${ZLIB_LIBRARIES})
