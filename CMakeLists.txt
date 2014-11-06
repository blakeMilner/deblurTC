cmake_minimum_required(VERSION 2.4)

PROJECT(Kisip)

MARK_AS_ADVANCED(CLEAR CMAKE_VERBOSE_MAKEFILE)
OPTION(CMAKE_VERBOSE_MAKEFILE "if all commands will be echoed to the console during the make" ON)

OPTION(ENABLE_DEBUG "enable debug support" OFF)
OPTION(ENABLE_X86_64 "enable x86_64 specific optimizations" OFF)

IF (ENABLE_DEBUG)
    SET (OPT_FLAGS "-DDEBUG -O0 -g -Wall -W -fno-inline")
ELSEIF (ENABLE_X86_64)
    SET (OPT_FLAGS "-O5 -march=k8 -mfpmath=sse -m64 -ffast-math -pipe")
ELSE()
    SET (OPT_FLAGS "-O5 -ffast-math")		# This is the default
ENDIF()

SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OPT_FLAGS}")

SET(MATH_LIBRARY "-lm")

MARK_AS_ADVANCED(CMAKE_C_FLAGS)

MACRO(CHECK_REQUIRED_LIB2 upper lower header lower2 header2)
    FIND_PATH(${upper}_INCLUDE_PATH ${header} $ENV{${upper}DIR}/include ${INCLUDE_DIRS})
    IF(${upper}_INCLUDE_PATH)
    ELSE (${upper}_INCLUDE_PATH)
                FIND_PATH(${upper}_INCLUDE_PATH ${header2} $ENV{${upper}DIR}/include ${INCLUDE_DIRS})
    ENDIF(${upper}_INCLUDE_PATH)
    
    IF(${upper}_INCLUDE_PATH)
                FIND_LIBRARY(${upper}_LIBRARY NAMES ${lower} ${lower2} PATHS $ENV{${upper}DIR}/lib ${LIB_DIRS})         
                IF(${upper}_LIBRARY)
                ELSE(${upper}_LIBRARY)
                        MESSAGE(SEND_ERROR "ERROR: ${upper} not found. please install ${upper} first!")
                ENDIF(${upper}_LIBRARY)
    ELSE(${upper}_INCLUDE_PATH)
                MESSAGE(SEND_ERROR "ERROR: ${upper} not found. please install ${upper} first!")
    ENDIF(${upper}_INCLUDE_PATH)
ENDMACRO(CHECK_REQUIRED_LIB2)

MACRO(CHECK_REQUIRED_LIB upper lower header)
    CHECK_REQUIRED_LIB2(${upper} ${lower} ${header} "" "")
ENDMACRO(CHECK_REQUIRED_LIB)

#
# this looks for GSL support
# it will define the following values
#
# GSL_INCLUDE_PATH = where gsl/gsl_linalg.h can be found
# GSL_LIBRARY      = the library to link against (gsl)
# GSL_CBLAS_LIBRARY  = the library to link against (gslcblas)

FIND_PATH(GSL_INCLUDE_PATH gsl/gsl_linalg.h
        $ENV{GSLDIR}/include $ENV{HOME}/include /usr/local/include /usr/include)

FIND_LIBRARY(GSL_LIBRARY
                NAMES gsl
                PATHS $ENV{GSLDIR}/lib $ENV{HOME}/lib /usr/lib /usr/lib64 /usr/local/lib)

FIND_LIBRARY(GSL_CBLAS_LIBRARY
                NAMES gslcblas
                PATHS $ENV{GSLDIR}/lib $ENV{HOME}/lib /usr/lib /usr/lib64 /usr/local/lib)

IF(GSL_LIBRARY)
ELSE(GSL_LIBRARY)
        MESSAGE(SEND_ERROR "ERROR: GSL not found. GSL is now required to use Kisip")
ENDIF(GSL_LIBRARY)

#
# this looks for FFTW3 support
# it will define the following values
#
# FFTW_INCLUDE_PATH = where fftw3.h can be found
# FFTW_LIBRARY      = the library to link against (fftw3)

FIND_PATH(FFTW_INCLUDE_PATH fftw3.h
        $ENV{FFTWDIR}/include $ENV{HOME}/include /usr/local/include /usr/include)

FIND_LIBRARY(FFTW_LIBRARY
                NAMES fftw3
                PATHS $ENV{FFTWDIR}/lib $ENV{HOME}/lib /usr/lib /usr/lib64 /usr/local/lib)

IF(FFTW_LIBRARY)
ELSE(FFTW_LIBRARY)
        MESSAGE(SEND_ERROR "ERROR: FFTW not found. FFTW is now required to use Kisip")
ENDIF(FFTW_LIBRARY)

# find mpi
INCLUDE(${CMAKE_ROOT}/Modules/FindMPI.cmake)

set(kisip_source_files				
						apodize.c assemble.c bs_ave.c bs_init.c 
						bs_red.c bs_snrt.c calc_subfields.c 
						chkphase.c ctracker.c entry.c 
						free_vect.c getinfo.c  
						init_matrix.c init_shift.c iwlspr.c 
						mpi_master_rec.c 
						mpi_mastering.c mpi_setmoddata.c mpi_setslavetype.c 
						mpi_shutdown.c mpi_slaving.c phs_init.c 
						putslaveinfo.c radial.c 
						reconstruct.c rpr.c 
						stats.c subfielding.c surfit.c 
						w_func.c 
)

set(kisip_source_files ${kisip_source_files} fileops.c)

ADD_EXECUTABLE(kisip ${kisip_source_files})
SET_TARGET_PROPERTIES(kisip PROPERTIES LINKER_LANGUAGE C)

INCLUDE_DIRECTORIES(${Kisip_SOURCE_DIR} ${FFTW_INCLUDE_PATH} ${GSL_INCLUDE_PATH} ${MPI_INCLUDE_PATH})

TARGET_LINK_LIBRARIES(kisip ${MATH_LIBRARY} ${FFTW_LIBRARY} ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY} ${MPI_LIBRARY} ${MPI_EXTRA_LIBRARY})

INSTALL_TARGETS(/bin kisip)



