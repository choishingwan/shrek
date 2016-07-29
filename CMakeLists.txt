project(shrek)
cmake_minimum_required(VERSION 3.1 FATAL_ERROR)
include(ExternalProject)

# Set the output folder where your program will be created
set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/bin)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})
set(LIBRARY_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/lib)
include_directories(inc)

# ==========
# Check compiler version
# ==========
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
	if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 4.8)
		message(FATAL_ERROR "GCC version must be at least 4.8!")
	endif()

elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")

	if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 3.2)
		message(FATAL_ERROR "Clang version must be at least 3.2!")
	endif()

else()
	message(WARNING "You are using an unsupported compiler!
	Compilation has only been tested with Clang and GCC.")
endif()

# ================
# External projects
# ================

# Armadillo
#============
ExternalProject_Add(ArmadilloDownload
	PREFIX ${LIBRARY_OUTPUT_PATH}/armadillo-7.100.3
	DOWNLOAD_COMMAND wget http://sourceforge.net/projects/arma/files/armadillo-7.100.3.tar.xz
	CONFIGURE_COMMAND tar -xf ../armadillo-7.100.3.tar.xz -C ../ArmadilloDownload  --strip-components=1
	BUILD_IN_SOURCE 1
	BUILD_COMMAND cmake .
	INSTALL_COMMAND make
	UPDATE_COMMAND ""
	)

set(Armadillo_DIR ${LIBRARY_OUTPUT_PATH}/armadillo-7.100.3/src/ArmadilloDownload)
set(ARMADILLO_LIBRARY ${Armadillo_DIR}/libarmadillo.so)
set(ARMADILLO_INCLUDE_DIR ${Armadillo_DIR}/include/)
include_directories(${ARMADILLO_INCLUDE_DIR})

# OpenBlas
#============
ExternalProject_Add(OpenBlasDownload
	PREFIX ${LIBRARY_OUTPUT_PATH}/OpenBLAS-0.2.18
	GIT_REPOSITORY https://github.com/xianyi/OpenBLAS.git
	GIT_TAG "v0.2.18"
	CONFIGURE_COMMAND ""
	BUILD_IN_SOURCE 1
	BUILD_COMMAND make
	INSTALL_COMMAND ""
	UPDATE_COMMAND ""
	)

set(OpenBlas_DIR ${LIBRARY_OUTPUT_PATH}/OpenBLAS-0.2.18/src/OpenBlasDownload)
set(OPENBLAS_LIBRARY ${OpenBlas_DIR}/libopenblas.so)

FIND_PACKAGE(OpenMP REQUIRED)
if(OPENMP_FOUND)
	message("OPENMP FOUND")
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
	set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

# Boost
#========
set(BOOST_DIR ${LIBRARY_OUTPUT_PATH}/boost)
include_directories(BOOST_DIR)

# Custom Lib
#============

# usefulltools
add_library(usefultools SHARED "src/usefulTools.cpp") 
add_dependencies(usefultools ArmadilloDownload)
target_compile_features(usefultools PRIVATE cxx_range_for)

# genotypefilehandler
add_library(genotypefilehandler SHARED "src/genotypefilehandler.cpp") 
add_dependencies(genotypefilehandler ArmadilloDownload)
target_compile_features(genotypefilehandler PRIVATE cxx_range_for)

# command
add_library(commander SHARED "src/command.cpp") 
add_dependencies(commander ArmadilloDownload)
target_compile_features(commander PRIVATE cxx_range_for)

# snp
add_library(snp SHARED "src/snp.cpp") 
add_dependencies(snp ArmadilloDownload)
target_compile_features(snp PRIVATE cxx_range_for)

# region
add_library(region SHARED "src/region.cpp") 
add_dependencies(region ArmadilloDownload)
target_compile_features(region PRIVATE cxx_range_for)

# snpestimation
add_library(snpestimation SHARED "src/snpestimation.cpp") 
add_dependencies(snpestimation ArmadilloDownload)
target_compile_features(snpestimation PRIVATE cxx_range_for)

# genotype
add_library(genotype SHARED "src/genotype.cpp") 
add_dependencies(genotype ArmadilloDownload)
target_compile_features(genotype PRIVATE cxx_range_for)

# linkage
add_library(linkage SHARED "src/linkage.cpp") 
add_dependencies(linkage ArmadilloDownload)
target_compile_features(linkage PRIVATE cxx_range_for)

# decomposition
add_library(decomposition SHARED "src/decomposition.cpp") 
add_dependencies(decomposition ArmadilloDownload)
target_compile_features(decomposition PRIVATE cxx_range_for)

add_executable(shrek ${CMAKE_SOURCE_DIR}/src/main.cpp)
add_dependencies(shrek ArmadilloDownload OpenBlasDownload)
target_link_libraries(shrek
	genotypefilehandler
	usefultools
	commander
	snp
	region
	snpestimation
	genotype
	linkage
	decomposition
	${OPENBLAS_LIBRARY}
	${ARMADILLO_LIBRARIES}
	-O2 lapack -pthread -fopenmp
	)
target_compile_features(shrek PRIVATE cxx_range_for)
