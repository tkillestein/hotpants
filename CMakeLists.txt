cmake_minimum_required(VERSION 3.20)
project(hotpants VERSION 5.1.11 LANGUAGES C)

set(CMAKE_C_STANDARD 17)
set(CMAKE_C_STANDARD_REQUIRED ON)

add_compile_options(-funroll-loops -O3 -Wall -D_GNU_SOURCE)

# ---------------------------------------------------------------------------
# Dependencies
# ---------------------------------------------------------------------------

find_package(PkgConfig REQUIRED)
pkg_check_modules(CFITSIO REQUIRED cfitsio)

# Prefer OpenBLAS over Netlib BLAS for better single-threaded performance
# and support for runtime threading control via openblas_set_num_threads().
find_library(OPENBLAS_LIBRARY NAMES openblas)
if(OPENBLAS_LIBRARY)
    message(STATUS "Found OpenBLAS: ${OPENBLAS_LIBRARY}")
    set(BLAS_LIBRARIES ${OPENBLAS_LIBRARY})
endif()

find_package(LAPACK REQUIRED)
find_library(LAPACKE_LIBRARY NAMES lapacke REQUIRED
    DOC "LAPACKE C interface library")
find_path(LAPACKE_INCLUDE_DIR NAMES lapacke.h REQUIRED
    DOC "Directory containing lapacke.h")

# OpenMP is optional — wired up for parallelism work (FFTW3 and region loops).
find_package(OpenMP)

# FFTW3 is required for FFT-accelerated convolution.
# Also find fftw3_omp for multi-threaded FFT execution.
pkg_check_modules(FFTW3 REQUIRED fftw3)
find_library(FFTW3_OMP_LIBRARY NAMES fftw3_omp fftw3_omp_LIBRARY
    DOC "FFTW3 OpenMP support library")
if(FFTW3_OMP_LIBRARY)
    message(STATUS "Found FFTW3 OpenMP library: ${FFTW3_OMP_LIBRARY}")
else()
    message(WARNING "FFTW3 OpenMP library (libfftw3_omp) not found; FFT will be single-threaded. "
                    "On Debian/Ubuntu: apt-get install libfftw3-dev. "
                    "On macOS: brew install fftw --with-openmp")
endif()
message(STATUS "FFTW3 found: ${FFTW3_LIBRARIES}")

# ---------------------------------------------------------------------------
# hotpants (shared library for Python and C API)
# ---------------------------------------------------------------------------

add_library(hotpants SHARED
    src/globals.c
    src/allocate.c
    src/basis.c
    src/basis_gaussian.c
    src/basis_delta.c
    src/alard.c
    src/functions.c
    src/hotpants_wrapper.c
)

target_include_directories(hotpants PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/src
    ${CFITSIO_INCLUDE_DIRS}
    ${LAPACKE_INCLUDE_DIR}
    ${FFTW3_INCLUDE_DIRS}
)

target_link_libraries(hotpants PUBLIC
    ${CFITSIO_LIBRARIES}
    LAPACK::LAPACK
    ${LAPACKE_LIBRARY}
    ${FFTW3_LIBRARIES}
    m
)

if(OpenMP_C_FOUND)
    target_link_libraries(hotpants PUBLIC OpenMP::OpenMP_C)
endif()

if(FFTW3_OMP_LIBRARY)
    target_link_libraries(hotpants PUBLIC ${FFTW3_OMP_LIBRARY})
endif()

# ---------------------------------------------------------------------------
# hotpants CLI (main executable, uses the shared library)
# ---------------------------------------------------------------------------

add_executable(hotpants_cli
    src/main.c
    src/vargs.c
)

target_include_directories(hotpants_cli PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src
    ${CFITSIO_INCLUDE_DIRS}
    ${LAPACKE_INCLUDE_DIR}
)

target_link_libraries(hotpants_cli PRIVATE
    hotpants
)

# Rename to 'hotpants' in installation
set_target_properties(hotpants_cli PROPERTIES OUTPUT_NAME hotpants)

# ---------------------------------------------------------------------------
# maskim (standalone mask application utility)
# ---------------------------------------------------------------------------

add_executable(maskim src/maskim.c)

target_include_directories(maskim PRIVATE ${CFITSIO_INCLUDE_DIRS} src)
target_link_libraries(maskim PRIVATE ${CFITSIO_LIBRARIES} m)

# ---------------------------------------------------------------------------
# Install
# ---------------------------------------------------------------------------

install(TARGETS hotpants hotpants_cli maskim
    LIBRARY DESTINATION lib
    RUNTIME DESTINATION bin
)

install(FILES src/hotpants_api.h
    DESTINATION include
)
