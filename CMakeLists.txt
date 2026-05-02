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

find_package(LAPACK REQUIRED)
find_library(LAPACKE_LIBRARY NAMES lapacke REQUIRED
    DOC "LAPACKE C interface library")
find_path(LAPACKE_INCLUDE_DIR NAMES lapacke.h REQUIRED
    DOC "Directory containing lapacke.h")

# OpenMP is optional — wired up for future parallelism work.
find_package(OpenMP)

# ---------------------------------------------------------------------------
# hotpants (main executable)
# ---------------------------------------------------------------------------

add_executable(hotpants
    src/main.c
    src/alard.c
    src/functions.c
    src/vargs.c
)

target_include_directories(hotpants PRIVATE
    ${CFITSIO_INCLUDE_DIRS}
    ${LAPACKE_INCLUDE_DIR}
    src
)

target_link_libraries(hotpants PRIVATE
    ${CFITSIO_LIBRARIES}
    LAPACK::LAPACK
    ${LAPACKE_LIBRARY}
    m
)

if(OpenMP_C_FOUND)
    target_link_libraries(hotpants PRIVATE OpenMP::OpenMP_C)
endif()

# ---------------------------------------------------------------------------
# maskim (standalone mask application utility)
# ---------------------------------------------------------------------------

add_executable(maskim src/maskim.c)

target_include_directories(maskim PRIVATE ${CFITSIO_INCLUDE_DIRS} src)
target_link_libraries(maskim PRIVATE ${CFITSIO_LIBRARIES} m)

# ---------------------------------------------------------------------------
# Install
# ---------------------------------------------------------------------------

install(TARGETS hotpants maskim DESTINATION bin)
