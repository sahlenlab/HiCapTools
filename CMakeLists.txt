cmake_minimum_required (VERSION 3.1)

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
  add_library(libbamtools STATIC IMPORTED)
  if(MAC)
	set_target_properties(libbamtools PROPERTIES IMPORTED_LOCATION "${CMAKE_SOURCE_DIR}/bamtools/libbamtools.2.4.1.dylib")	
  endif()
  if(UNIX AND NOT MAC)
	set_target_properties(libbamtools PROPERTIES IMPORTED_LOCATION "${CMAKE_SOURCE_DIR}/bamtools/libbamtools.so.2.4.1")
  endif()
  include_directories( "${CMAKE_SOURCE_DIR}/boost")
  include_directories( "${CMAKE_SOURCE_DIR}/ExterAlg")
  include_directories("${CMAKE_SOURCE_DIR}/include")
  include_directories("${CMAKE_SOURCE_DIR}/external")
  target_link_libraries( HiCapTools libbamtools ${ZLIB_LIBRARIES})
