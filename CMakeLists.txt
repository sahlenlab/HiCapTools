cmake_minimum_required (VERSION 3.5)

project(HiCapTools CXX)

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
    if (CMAKE_CXX_COMPILER_VERSION VERSION_GREATER 5.1)
		set (CMAKE_CXX_FLAGS "-D_GLIBCXX_USE_CXX11_ABI=0")
	endif()
 endif()

file(GLOB_RECURSE SOURCES "${CMAKE_SOURCE_DIR}/external/*.cpp")
file(GLOB_RECURSE SOURCES "${CMAKE_SOURCE_DIR}/ExterAlg/*.cpp")

# set our library and executable destination dirs
  set( EXECUTABLE_OUTPUT_PATH "${CMAKE_SOURCE_DIR}/bin" )
  
  file (GLOB CPPFILES "${CMAKE_SOURCE_DIR}/src/*.cpp")
  file (GLOB HFILES "${CMAKE_SOURCE_DIR}/include/*.h") 
  file (GLOB HPPFILES "${CMAKE_SOURCE_DIR}/include/*.hpp")

  add_executable(HiCapTools ${CPPFILES} ${HFILES} ${HPPFILES} ${SOURCES})
   
  set_property(TARGET HiCapTools PROPERTY CXX_STANDARD 11)
  set_property(TARGET HiCapTools PROPERTY CXX_STANDARD_REQUIRED ON)

  find_package(ZLIB REQUIRED)
  find_package(Boost REQUIRED)
  add_library(libbamtools UNKNOWN IMPORTED)
  if(MAC)
    set(BAMTOOLS_ROOT "" CACHE PATH "Root directory of a BamTools installation")
    find_path(BAMTOOLS_INCLUDE_DIR
      NAMES api/BamReader.h
      HINTS
        "${BAMTOOLS_ROOT}/include/bamtools"
        "$ENV{BAMTOOLS_ROOT}/include/bamtools"
      PATHS
        /opt/homebrew/include/bamtools
        /usr/local/include/bamtools
    )
    find_library(BAMTOOLS_LIBRARY
      NAMES bamtools
      HINTS
        "${BAMTOOLS_ROOT}/lib"
        "$ENV{BAMTOOLS_ROOT}/lib"
      PATHS
        /opt/homebrew/lib
        /usr/local/lib
    )
    if(NOT BAMTOOLS_INCLUDE_DIR OR NOT BAMTOOLS_LIBRARY)
      message(FATAL_ERROR
        "A native BamTools installation is required on macOS. "
        "Install it with 'brew install bamtools', or set BAMTOOLS_ROOT."
      )
    endif()
    set_target_properties(libbamtools PROPERTIES IMPORTED_LOCATION "${BAMTOOLS_LIBRARY}")
    target_include_directories(HiCapTools BEFORE PRIVATE "${BAMTOOLS_INCLUDE_DIR}")
    message(STATUS "Using system BamTools: ${BAMTOOLS_LIBRARY}")
  elseif(UNIX)
		set_target_properties(libbamtools PROPERTIES IMPORTED_LOCATION "${CMAKE_SOURCE_DIR}/bamtools/libbamtools.so.2.4.1")
  endif()
  include_directories(SYSTEM ${Boost_INCLUDE_DIRS})
  message(STATUS "Using system Boost headers: ${Boost_INCLUDE_DIRS}")
  include_directories( "${CMAKE_SOURCE_DIR}/ExterAlg")
  include_directories("${CMAKE_SOURCE_DIR}/include")
  include_directories("${CMAKE_SOURCE_DIR}/external")
  target_link_libraries( HiCapTools libbamtools ${ZLIB_LIBRARIES})
