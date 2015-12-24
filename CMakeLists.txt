cmake_minimum_required (VERSION 2.8.7) # >= 2.8.8 necessary for $<TARGET_OBJECTS:objlib>
include (CheckSymbolExists)
include (CheckCXXSymbolExists)
include (CheckIncludeFiles)
include (CheckIncludeFileCXX)
include (CheckCXXSourceCompiles)
include (GenerateExportHeader)
include (UsePkgConfig)
include (ProcessorCount)

project (oofem CXX C)
set (oofem_VERSION_MAJOR 2)
set (oofem_VERSION_MINOR 4)
set (oofem_VERSION_PATCH 0) # Unused
set (oofem_HOMEPAGE "http://www.oofem.org")
set (bp_HOMEPAGE "http://ksm.fsv.cvut.cz/%7Ebp/bp.html")
set (oofem_COPYRIGHT "Copyright (C) 1994-2014 Borek Patzak")
set (oofem_EMAIL "info@oofem.org")

# Default is debug build
if (NOT CMAKE_BUILD_TYPE)
    set (CMAKE_BUILD_TYPE Debug)
endif ()

if (MSVC)
    add_definitions (-D_USE_MATH_DEFINES)
    set (CMAKE_REQUIRED_DEFINITIONS "-D_USE_MATH_DEFINES")
    # I'm not sure we should block these warnings or not:
    # C4800 Conversion from int to bool performance warning (we should probably be using bool in most of these places though).
    # C4127 Conditional expression is constant: Most of these uses are intentional (debugging code etc.)
    # C4267 size_t to int conversion. Most of these are indices, but many places can be replaced by a for-each loop as well.
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -W4 -wd4786 -wd4100 -wd4996 -wd4706 -wd4512 -wd4251 -wd4800")
else ()
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -fPIC")
endif ()
# Check with more compilers as time go on
if (CMAKE_COMPILER_IS_GNUCXX AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 4.7)
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "Clang" OR CMAKE_COMPILER_IS_GNUCXX)
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
elseif (MSVC)
    # Nothing to do for MSVC.
endif ()

set (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DDEBUG")
set (CMAKE_CXX_FLAGS_PROFILING "${CMAKE_CXX_FLAGS_PROFILING} -O2 -pg")

#######################################################################
######################## Custom paths #################################
#######################################################################
set (SLEPC_DIR CACHE PATH "Manual SLEPc directory")
set (X11_DIR CACHE PATH "Manual X11 directory (OOFEG)")
set (TINYXML2_DIR CACHE PATH "Manual TinyXML-2 directory (CEMHYD only)")
set (MY_PETSC_DIR CACHE PATH "Manual PETSc directory")
set (MY_PETSC_ARCH CACHE STRING "Manual PETSc arch name (e.g. linux-gnu-cxx-opt)")
set (MPI_DIR CACHE PATH "Manual MPI directory")
set (VTK_DIR CACHE PATH "Manual VTK directory (For binary VTU export)")
set (LAPACK_DIR CACHE PATH "Manual LAPACK directory")
set (SPOOLES_DIR CACHE PATH "Manual SPOOLES directory")
set (ELIXIR_DIR CACHE PATH "Manual Elixir directory (OOFEG only)")
set (CKIT_DIR CACHE PATH "Manual CKit directory (OOFEG only)")
set (METIS_DIR CACHE PATH "Manual Metis directory (Load balancing)")
set (PARMETIS_DIR CACHE PATH "Manual ParMetis directory (Load balancing)")
#set (BOOST_DIR CACHE PATH "Manual Boost directory")

#######################################################################
######################## Options ######################################
#######################################################################
option (USE_OOFEG "Enable OOFEG" OFF)
if (${CMAKE_VERSION} VERSION_GREATER "2.8.7") # Compability with older version of cmake
    option (USE_SHARED_LIB "Build shared library" ON)
else ()
    set (USE_SHARED_LIB OFF)
endif ()

# Parallel stuff
option (USE_PARALLEL "Enable parallel support" OFF)
option (USE_METIS "Enable metis support" OFF)
option (USE_PARMETIS "Enable Parmetis support" OFF)
option (USE_OPENMP "Compile with OpenMP support (for parallel assembly)" OFF)
# Solvers and such
option (USE_DSS "Enable DSS module" OFF) # No reason to use this
option (USE_IML "Enable iml++ solvers" OFF) # or this
option (USE_SPOOLES "Enable SPOOLES-solver" OFF) # or this
option (USE_MKLPARDISO "Enable MKLPARDISO solver support" OFF)
option (USE_PETSC "Enable PETSc support (required for parallel)" OFF)
option (USE_SLEPC "Enable SLEPc support" OFF)
option (USE_LAPACK "Use LAPACK" OFF)
# Other external libraries
option (USE_TRIANGLE "Compile with Triangle bindings" OFF)
option (USE_VTK "Enable VTK (for exporting binary VTU-files)" OFF)
#option (USE_CGAL "CGAL" OFF)
# Internal modules
option (USE_SM "Enable structural mechanics module" ON)
option (USE_FM "Enable fluid mechanics module" ON)
option (USE_PFEM "Enable PFEM module" OFF)
option (USE_TM "Enable heat/mass transfer module" ON)
option (USE_CEMHYD "Enable CemHyd support" OFF)
#option (USE_BOOST "Enable Boost" OFF)
option (USE_PYTHON_BINDINGS "Enable Python bidings (OOFEM usable from Python). Implies USE_PYTHON." OFF)
option (USE_PYTHON "Enable Python to be used inside OOFEM (e.g. userdef Dirichlet BCs)" OFF)

if (USE_PARALLEL)
    add_definitions (-D__PARALLEL_MODE)
    set (USE_MPI ON)
endif ()

#######################################################################
######################## Internal libraries ###########################
#######################################################################
# I believe the subdirs should be used when including instead.
include_directories ("${oofem_SOURCE_DIR}/src/oofemlib")
include_directories ("${oofem_BINARY_DIR}")

# Is there really any reason not to compile all of them? Some has cross dependencies.
if (${CMAKE_VERSION} VERSION_LESS "2.8.8")
    list (APPEND LIBS core)
else ()
    list (APPEND LIBS $<TARGET_OBJECTS:core>)
endif ()

if (USE_FM)
    include_directories ("${oofem_SOURCE_DIR}/src/fm")
    add_definitions (-D__FM_MODULE)
    if (${CMAKE_VERSION} VERSION_LESS "2.8.8")
        list (APPEND LIBS fm)
    else ()
        list (APPEND LIBS $<TARGET_OBJECTS:fm>)
    endif ()
    list (APPEND MODULE_LIST "fm")
    if (USE_PFEM)
        add_definitions (-D__PFEM_MODULE)
    endif ()
endif ()

if (USE_TM)
    include_directories ("${oofem_SOURCE_DIR}/src/tm")
    add_definitions (-D__TM_MODULE)
    if (${CMAKE_VERSION} VERSION_LESS "2.8.8")
        list (APPEND LIBS tm)
    else ()
        list (APPEND LIBS $<TARGET_OBJECTS:tm>)
    endif ()
    list (APPEND MODULE_LIST "tm")
    if (USE_CEMHYD)
        add_definitions (-D__CEMHYD_MODULE)
        set (USE_TINYXML ON)
    endif ()
endif ()

if (USE_SM)
    include_directories ("${oofem_SOURCE_DIR}/src/sm")
    add_definitions (-D__SM_MODULE)
    if (${CMAKE_VERSION} VERSION_LESS "2.8.8")
        list (APPEND LIBS sm)
    else ()
        list (APPEND LIBS $<TARGET_OBJECTS:sm>)
    endif ()
    list (APPEND MODULE_LIST "sm")

    if (NOT MSVC AND NOT MINGW)
        find_library (DL_LIB dl) 
        list (APPEND EXT_LIBS ${DL_LIB})
    endif ()
    
#     if (MINGW)
#         find_library (DL_LIB ltdl HINTS "C:/MinGW/lib/" ) 
#         list (APPEND EXT_LIBS ${DL_LIB})
#     endif ()
    
endif ()

if (USE_DSS)
    include_directories ("${oofem_SOURCE_DIR}/src/dss")
    if (${CMAKE_VERSION} VERSION_LESS "2.8.8")
        list (APPEND LIBS dss)
    else ()
        list (APPEND LIBS $<TARGET_OBJECTS:dss>)
    endif ()
    list (APPEND MODULE_LIST "dss")
endif ()

if (USE_OPENMP)
    include (FindOpenMP)
    if (OPENMP_FOUND)
        set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
        set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
        set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
        list (APPEND MODULE_LIST "openmp")
    endif ()
endif ()

if (USE_OOFEG)
    add_definitions (-D__OOFEG)

    # Necessary X11 stuff;
    if (X11_DIR)
        set (CMAKE_LIBRARY_DIR ${X11_DIR})
        find_package (X11 REQUIRED PATHS ${X11_DIR})
        find_library (XAW_LIB Xaw PATHS ${X11_DIR})
    else ()
        find_package (X11 REQUIRED)
        find_library (XAW_LIB Xaw)
    endif ()

    # Krysls libraries:
    find_library (ELIXIR_LIB elixir PATH "${ELIXIR_DIR}/src")
    find_library (ESI_LIB esi PATH "${ELIXIR_DIR}/src")
    if (${ELIXIR_LIB} STREQUAL "ELIXIR_LIB-NOTFOUND")
        message (FATAL_ERROR "Elixir library not found (required for OOFEG)")
    endif ()
    find_library (CKIT_LIB ckit PATH "${CKIT_DIR}/src")
    if (${CKIT_LIB} STREQUAL "CKIT_LIB-NOTFOUND")
        message (FATAL_ERROR "Ckit library not found (required for OOFEG)")
    endif ()
    include_directories ("${ELIXIR_DIR}/include")
    include_directories ("${CKIT_DIR}/include")

    list (APPEND EXT_LIBS ${ESI_LIB} ${ELIXIR_LIB} ${CKIT_LIB} ${X11_LIBRARIES} ${X11_Xmu_LIB} ${X11_Xt_LIB} ${XAW_LIB})
    list (APPEND MODULE_LIST "oofeg")
endif ()

#######################################################################
######################## External libraries ###########################
#######################################################################

#if (USE_BOOST) 
#    set (CMAKE_LIBRARY_DIR ${BOOST_DIR})
#    find_package (BOOST REQUIRED)
#    list (APPEND EXT_LIBS ${BOOST_LIBRARIES})
#    if (BOOST_DIR)
#    	include_directories ("${BOOST_DIR}")
#    	list (APPEND MODULE_LIST "BOOST")
#     endif ()
# 
#    add_definitions (-D__BOOST_MODULE)
#endif ()

if (USE_LAPACK)
    set (CMAKE_LIBRARY_DIR ${LAPACK_DIR})
    find_package (LAPACK REQUIRED)
    list (APPEND EXT_LIBS ${LAPACK_LIBRARIES})
    list (APPEND MODULE_LIST "LAPACK")
    add_definitions (-D__LAPACK_MODULE)
endif ()

if (USE_VTK)
    set (CMAKE_LIBRARY_DIR ${VTK_DIR})
    find_package (VTK REQUIRED NO_MODULE)
    include (${VTK_USE_FILE}) # Sets necessary includes_directories
    add_definitions (-D__VTK_MODULE)
    list (APPEND EXT_LIBS ${VTK_LIBRARIES})
    list (APPEND MODULE_LIST "VTK")
endif ()

if (USE_PARMETIS)
    if (PARMETIS_DIR)
        find_library (PARMETIS_LIB parmetis PATH "${PARMETIS_DIR}/lib")
        include_directories ("${PARMETIS_DIR}/include")
    else ()
        find_library (PARMETIS_LIB parmetis)
    endif ()
    if (${PARMETIS_LIB} STREQUAL "PARMETIS_LIB-NOTFOUND")
        message (FATAL_ERROR "Parmetis library not found")
    endif ()
    list (APPEND EXT_LIBS ${PARMETIS_LIB})
    list (APPEND MODULE_LIST "parmetis")
    set (USE_MPI ON)
    set (USE_METIS ON)
endif ()

if (USE_METIS)
    if (METIS_DIR)
        find_library (METIS_LIB metis PATH "${METIS_DIR}/lib")
        include_directories ("${METIS_DIR}/include")
    else ()
        find_library (METIS_LIB metis)
    endif ()
    if (${METIS_LIB} STREQUAL "METIS_LIB-NOTFOUND")
        message (FATAL_ERROR "Metis library not found")
    endif ()
    list (APPEND EXT_LIBS ${METIS_LIB})
    list (APPEND MODULE_LIST "metis")
endif ()

if (USE_MKLPARDISO)
    # Find the actual libraries:
    list(APPEND MKL_LIBS "mkl_core" "mkl_gf_lp64" "mkl_gnu_thread" "pthread")
    foreach (f ${MKL_LIBS})
        if (MKL_DIR)
            find_library (${f}_LIB ${f} PATH "${MKL_DIR}/mkl/lib/intel64/" DOC "MKLPARDISO (library)")
            include_directories ("${MKL_DIR}/mkl/include")
        else ()
            find_library (${f}_LIB ${f} HINTS "/opt/intel/composerxe/mkl/lib/intel64/" DOC "MKLPARDISO (library)")
            include_directories ("/opt/intel/composerxe/mkl/include")
        endif ()
        if (${${f}_LIB} STREQUAL "${f}_LIB-NOTFOUND")
            message (FATAL_ERROR "MKLPARDISO library not found")
        endif ()
        list (APPEND EXT_LIBS ${${f}_LIB})
    endforeach ()

    #list (APPEND EXT_LIBS ${MKLPARDISO_LIB})
    list (APPEND MODULE_LIST "MKLPARDISO")
endif ()

if (USE_SLEPC)
    if (SLEPC_DIR)
        find_library (SLEPC_LIBS slepc PATH "${SLEPC_DIR}/lib")
        include_directories ("${SLEPC_DIR}/include")
    else ()
        find_library (SLEPC_LIBS slepc)
    endif ()
    if (${SLEPC_LIBS} STREQUAL "SLEPC_LIBS-NOTFOUND")
        message (FATAL_ERROR "SLEPc library not found")
    endif ()
    list (APPEND EXT_LIBS ${SLEPC_LIBS})
    list (APPEND MODULE_LIST "slepc")
    set (USE_PETSC ON)
    add_definitions (-D__SLEPC_MODULE)
endif ()

if (USE_PETSC)
    # A real PETScConfig.cmake file would have been nice to use:
    #find_package (PETSc 3.3 HINTS "${MY_PETSC_DIR}/${MY_PETSC_ARCH}")

    # Have to use pkgconfig information
    find_package (PkgConfig REQUIRED)
    
    set (ENV{PKG_CONFIG_PATH} "${MY_PETSC_DIR}/${MY_PETSC_ARCH}/lib64/pkgconfig:${MY_PETSC_DIR}/${MY_PETSC_ARCH}/lib/pkgconfig")
    pkg_check_modules (PC_PETSC REQUIRED PETSc)
    include_directories (${PC_PETSC_INCLUDE_DIRS})

    # Find the actual libraries:
    foreach (f ${PC_PETSC_LIBRARIES})
        find_library (${f}_LIBRARY NAMES ${f} lib${f}
            HINTS
            ${PC_PETSC_LIBRARY_DIRS}
        )
        list (APPEND EXT_LIBS ${${f}_LIBRARY})
    endforeach ()
    # Currently unused flags:
    #   PC_PETSC_LDFLAGS        ... all required linker flags
    #   PC_PETSC_LDFLAGS_OTHER  ... all other linker flags
    #   PC_PETSC_CFLAGS         ... all required cflags
    #   PC_PETSC_CFLAGS_OTHER   ... the other compiler flags

    #find_package (PETSc REQUIRED HINTS "${MY_PETSC_DIR}/${MY_PETSC_ARCH}/conf")
    #include_directories ("${MY_PETSC_DIR}/include" "${MY_PETSC_DIR}/${MY_PETSC_ARCH}/include" ${PETSC_INCLUDE_PATH})
    #find_library (PETSC_LIB petsc PATH "${MY_PETSC_DIR}/${MY_PETSC_ARCH}/lib/") 
    #if (${PETSC_LIB} STREQUAL "PETSC_LIB-NOTFOUND")
    #    message (FATAL_ERROR "PETSc library not found")
    #endif ()
    #list (APPEND EXT_LIBS ${PETSC_LIB} ${PETSC_PACKAGE_LIBS})
    list (APPEND MODULE_LIST "PETSc")
    set (USE_MPI ON)
    add_definitions (-D__PETSC_MODULE)
endif ()

if (USE_MPI)
    if (MPI_DIR)
        include_directories ("${MPI_DIR}/include")
        find_library (MPI_LIBRARY NAMES mpi PATHS ${MPI_DIR}/lib64 ${MPI_DIR}/lib)
        find_library (MPI_CXX_LIBRARY NAMES mpi_cxx PATHS ${MPI_DIR}/lib64 ${MPI_DIR}/lib)
        list (APPEND EXT_LIBS ${MPI_LIBRARY})
        if (MPI_CXX_LIBRARY)
            list (APPEND EXT_LIBS ${MPI_CXX_LIBRARY})
        endif ()
    else ()
        include (FindMPI)
        find_package (MPI REQUIRED)
        set (CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILES_FLAGS} ${MPI_COMPILE_FLAGS})
        set (CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${MPI_LINK_FLAGS})
        add_definitions (-DMPICH_IGNORE_CXX_SEEK)
        include_directories (${MPI_INCLUDE_PATH})
        list (APPEND EXT_LIBS ${MPI_LIBRARIES})
    endif (MPI_DIR)
    list (APPEND MODULE_LIST "MPI")
    # Compile and link flags?
endif ()

if (USE_TRIANGLE)
    find_library (TRIANGLE_LIB triangle DOC "Triangle (library)")
    if (${TRIANGLE_LIB} STREQUAL "TRIANGLE_LIB-NOTFOUND")
        message (FATAL_ERROR "Triangle library not found")
    endif ()
    add_definitions (-D__TRIANGLE_MODULE)
    list (APPEND EXT_LIBS ${TRIANGLE_LIB})
    list (APPEND MODULE_LIST "Triangle")
endif ()

if (USE_TINYXML)
    if (TINYXML2_DIR)
        include_directories (${TINYXML2_DIR})
        find_library (TINYXML2_LIB NAMES tinyxml2 HINTS ${TINYXML2_DIR} NO_DEFAULT_PATH)
    else ()
        find_library (TINYXML2_LIB NAMES tinyxml2 )
    endif ()
    if (${TINYXML2_LIB} STREQUAL "TINYXML2_LIB-NOTFOUND")
        message (FATAL_ERROR "TinyXML2 library not found")
    endif ()
    list (APPEND EXT_LIBS ${TINYXML2_LIB})
    list (APPEND MODULE_LIST "TinyXML-2")
endif ()

if (USE_IML)
    include_directories ("${oofem_SOURCE_DIR}")
    set (CMAKE_REQUIRED_INCLUDES "${oofem_SOURCE_DIR}")
    check_include_file_cxx ("iml/cg.h" HAVE_IML_CG_H)
    check_include_file_cxx ("iml/gmres.h" HAVE_IML_GMRES_H)
    if ((NOT HAVE_IML_CG_H) OR (NOT HAVE_IML_GMRES_H))
        message (FATAL_ERROR "Necessary IML++ headers not found")
    endif ()
    list (APPEND MODULE_LIST "IML++")
endif ()

if (USE_SPOOLES)
    if (SPOOLES_DIR)
        find_library (SPOOLES_LIB spooles DOC "SPOOLES library" PATH "${SPOOLES_DIR}/lib")
        include_directories ("${SPOOLES_DIR}/include")
    else ()
        find_library (SPOOLES_LIB spooles DOC "SPOOLES library")
    endif ()
    if (${SPOOLES_LIB} STREQUAL "SPOOLES_LIB-NOTFOUND")
        message (FATAL_ERROR "SPOOLES not found")
    endif ()
    list (APPEND EXT_LIBS ${SPOOLES_LIB})
    list (APPEND MODULE_LIST "SPOOLES")
endif ()

# We need the interpreter for the test suite.
find_package(PythonInterp "2.7")

if (USE_PYTHON_BINDINGS OR USE_PYTHON)
    # PythonLibs in cmake is bugged, have to try and rely on PythonInterp to get the suitable major and minor versions.
    set (CMAKE_LIBRARY_DIR ${PYTHON_DIR})
    find_package(PythonInterp "2.7" REQUIRED)
    find_package(PythonLibs "${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR}" REQUIRED)
    include_directories (${PYTHON_INCLUDE_DIRS})
    list (APPEND EXT_LIBS ${PYTHON_LIBRARIES})
    list (APPEND MODULE_LIST "Python")
    add_definitions(-D__PYTHON_MODULE)

    if (USE_PYTHON_BINDINGS)
        find_package(Boost COMPONENTS "python-py${PYTHON_VERSION_MAJOR}${PYTHON_VERSION_MINOR}" REQUIRED)
        include_directories (${Boost_INCLUDE_DIRS})
        list (APPEND EXT_LIBS ${Boost_LIBRARIES})
        add_definitions(-DBOOST_PYTHON)
    endif ()
endif ()


#######################################################################
######################## Configuration ################################
#######################################################################
site_name (HOST_NAME)
set (HOST_TYPE "${CMAKE_SYSTEM_PROCESSOR}-${CMAKE_SYSTEM_NAME}") # This should be renamed TARGET_TYPE.
set (BUILDNAME "OOFEM(${MODULE_LIST})-${HOST_TYPE}-${HOST_NAME}")

check_include_file_cxx ("unistd.h" HAVE_UNISTD_H)
check_include_file_cxx ("execinfo.h" HAVE_EXECINFO_H)
check_cxx_symbol_exists ("access" "unistd.h" HAVE_ACCESS)
check_cxx_symbol_exists ("cbrt" "cmath" HAVE_CBRT)
check_cxx_symbol_exists ("isnan" "cmath" HAVE_ISNAN)
if (USE_SHARED_LIB)
    set (OOFEM_EXPORT_HEADER "#include \"oofem_export.h\"")
    # For consistant library on linux and windows platforms.
    add_compiler_export_flags (LIBOOFEM_EXPORT_CFLAGS)
else ()
    set (OOFEM_EXPORT_HEADER "#define OOFEM_EXPORT\n#define OOFEM_NO_EXPORT")
endif ()
configure_file (
  "${oofem_SOURCE_DIR}/src/oofemcfg.h.in"
  "${oofem_BINARY_DIR}/oofemcfg.h"
  )

#######################################################################
######################## Targets ######################################
#######################################################################
# Subtargets
if (${CMAKE_VERSION} VERSION_LESS "2.8.8")
    # This won't work on windows, but windows users should have newer version than 2.8.7)
    set (LIB_TYPE SHARED)
else ()
    set (LIB_TYPE OBJECT)
endif ()
add_subdirectory ("${oofem_SOURCE_DIR}/src/oofemlib")
if (USE_DSS)
    add_subdirectory ("${oofem_SOURCE_DIR}/src/dss")
endif ()
if (USE_SM)
    add_subdirectory ("${oofem_SOURCE_DIR}/src/sm")
endif ()
if (USE_TM)
    add_subdirectory ("${oofem_SOURCE_DIR}/src/tm")
endif ()
if (USE_FM)
    add_subdirectory ("${oofem_SOURCE_DIR}/src/fm")
endif ()
set (oofempy)
if (USE_PYTHON_BINDINGS)
    set (oofempy ${oofem_SOURCE_DIR}/bindings/python/oofemlib.cpp)
endif ()

# Main targets
if (USE_SHARED_LIB)
    add_library (liboofem
        SHARED
        ${oofempy}
        ${LIBS} 
        )
    target_link_libraries (liboofem ${EXT_LIBS})
    set_target_properties (liboofem PROPERTIES OUTPUT_NAME oofem)
    if (NOT "${LIBOOFEM_EXPORT_CFLAGS}" STREQUAL "")
        set_target_properties (liboofem PROPERTIES COMPILE_FLAGS ${LIBOOFEM_EXPORT_CFLAGS})
    endif ()
    install (TARGETS liboofem DESTINATION lib)

    # This should be set automatically, but I think TARGET:objs confuses cmake:
    set_target_properties (core PROPERTIES COMPILE_DEFINITIONS "liboofem_EXPORTS")
    if (USE_FM)
        set_target_properties (fm PROPERTIES COMPILE_DEFINITIONS "liboofem_EXPORTS")
        if (NOT "${LIBOOFEM_EXPORT_CFLAGS}" STREQUAL "")
            set_target_properties (fm PROPERTIES COMPILE_FLAGS ${LIBOOFEM_EXPORT_CFLAGS})
        endif ()
    endif()
    if (USE_SM)
        set_target_properties (sm PROPERTIES COMPILE_DEFINITIONS "liboofem_EXPORTS")
        if (NOT "${LIBOOFEM_EXPORT_CFLAGS}" STREQUAL "")
            set_target_properties (sm PROPERTIES COMPILE_FLAGS ${LIBOOFEM_EXPORT_CFLAGS})
        endif ()
    endif()
    if (USE_TM)
        set_target_properties (tm PROPERTIES COMPILE_DEFINITIONS "liboofem_EXPORTS")
        if (NOT "${LIBOOFEM_EXPORT_CFLAGS}" STREQUAL "")
            set_target_properties (tm PROPERTIES COMPILE_FLAGS ${LIBOOFEM_EXPORT_CFLAGS})
        endif ()
    endif()

    generate_export_header (liboofem
        BASE_NAME OOFEM
        )
endif ()

# OOFEM target
if (USE_SHARED_LIB) # Prefer dynamic lib if available.
    add_executable (oofem ${oofem_SOURCE_DIR}/src/main/main.C)
    target_link_libraries (oofem liboofem)
elseif (${CMAKE_VERSION} VERSION_LESS "2.8.8") # Hackish backwards compatibility 
    add_executable (oofem ${oofem_SOURCE_DIR}/src/main/main.C)
    target_link_libraries (oofem ${LIBS} ${EXT_LIBS})
else ()
    add_executable (oofem ${oofem_SOURCE_DIR}/src/main/main.C ${LIBS})
    target_link_libraries (oofem ${EXT_LIBS})
endif ()
install (TARGETS oofem DESTINATION bin)

if (USE_OOFEG)
    #  OOFEG target:
    if (USE_SHARED_LIB)
        add_executable (oofeg ${oofem_SOURCE_DIR}/src/main/oofeg.C)
        target_link_libraries (oofeg liboofem)
    elseif (${CMAKE_VERSION} VERSION_LESS "2.8.8") # Hackish backwards compatibility
        add_executable (oofeg ${oofem_SOURCE_DIR}/src/main/oofeg.C)
        target_link_libraries (oofeg ${LIBS} ${EXT_LIBS})
    else ()
        add_executable (oofeg ${oofem_SOURCE_DIR}/src/main/oofeg.C ${LIBS})
        target_link_libraries (oofeg ${EXT_LIBS})
    endif ()
    install(TARGETS oofeg DESTINATION bin)
endif ()

# Example of using liboofem with dynamic input record:
add_executable(beam01 ${oofem_SOURCE_DIR}/bindings/oofemlib/beam01.C)
set_target_properties(beam01 PROPERTIES EXCLUDE_FROM_ALL TRUE)
target_link_libraries (beam01 liboofem)



# CppCheck target (not built by default)
add_custom_target(cppcheck)
set_target_properties(cppcheck PROPERTIES EXCLUDE_FROM_ALL TRUE)
add_custom_command(TARGET cppcheck 
    COMMAND cppcheck "--xml" "--template=gcc" "--enable=all" 
    "-I${oofem_BINARY_DIR}/"
    "-I${oofem_SOURCE_DIR}/src/oofemlib/"
    "-I${oofem_SOURCE_DIR}/src/sm/"
    "-I${oofem_SOURCE_DIR}/src/fm/"
    "-I${oofem_SOURCE_DIR}/src/tm/"
    "-DDEBUG"
    ${oofem_SOURCE_DIR}/src/oofemlib/
    ${oofem_SOURCE_DIR}/src/sm/
    ${oofem_SOURCE_DIR}/src/fm/
    ${oofem_SOURCE_DIR}/src/tm/
    ${oofem_SOURCE_DIR}/src/main/
    "2>"
    "cppcheck.xml"
    COMMENT "Running cppcheck on entire source"
    )


#######################################################################
######################## Documentation ################################
#######################################################################
find_package (Doxygen)
if (DOXYGEN_FOUND)
    set (DOXYGEN_INPUT ${oofem_BINARY_DIR}/doxygen.cfg)
    set (DOXYGEN_OUTPUT ${oofem_BINARY_DIR}/doc/api/)
    file (MAKE_DIRECTORY ${DOXYGEN_OUTPUT})

    configure_file (
        "${oofem_SOURCE_DIR}/doc/refman/footer.html.in"
        "${oofem_BINARY_DIR}/footer.html")
    configure_file (
        "${oofem_SOURCE_DIR}/doc/refman/doxygen.cfg.in"
        "${oofem_BINARY_DIR}/doxygen.cfg")

    add_custom_target (
        apidoc
        COMMAND ${CMAKE_COMMAND} -E echo_append "Building API Documentation..."
        COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_INPUT}
        COMMAND ${CMAKE_COMMAND} -E echo "Done."
        WORKING_DIRECTORY ${oofem_BINARY_DIR}
        DEPENDS)
endif (DOXYGEN_FOUND)

configure_file (
    "${oofem_SOURCE_DIR}/doc/include.tex.in"
    "${oofem_BINARY_DIR}")

# TODO: Manuals not done yet
add_custom_target (
    elementlibmanual
    COMMAND pdflatex elementlibmanual --output=${oofem_BINARY_DIR}/doc/
    COMMAND pdflatex elementlibmanual --output=${oofem_BINARY_DIR}/doc/
    #COMMAND tex4ht elementlibmanual --output=${oofem_BINARY_DIR}/doc/
    WORKING_DIRECTORY ${oofem_SOURCE_DIR}/doc/elementlibmanual/
    )

#######################################################################
######################## Tests ########################################
#######################################################################
set (oofem_TEST_DIR "${oofem_SOURCE_DIR}/tests/")
set (oofem_cmd  $<TARGET_FILE:oofem>)
file(TO_CMAKE_PATH "${oofem_TEST_DIR}/valgrind.supp" MEMORYCHECK_SUPPRESSIONS_FILE)
include (CTest)

#set (CTEST_MEMORYCHECK_SUPPRESSIONS_FILE ${oofem_TEST_DIR}/tests/valgrind.supp)

if (USE_PARALLEL)
    set (par_dir ${oofem_TEST_DIR}/partests)
    file (GLOB parallel_tests RELATIVE "${par_dir}" "${par_dir}/*/")
    foreach (case ${parallel_tests})
        file (GLOB files "${par_dir}/${case}/${case}.oofem.in.*")
        list (LENGTH files num_files)
        add_test (NAME "partest_${case}" WORKING_DIRECTORY ${par_dir}/${case} COMMAND "mpirun" "-np" ${num_files} ${oofem_cmd} "-p" "-f" ${case}.oofem.in "-ksp_type" "cg")
    endforeach (case)
    # partest_brazil_2d_nl7 requires more time on testing server than default 1500
    set_tests_properties(partest_brazil_2d_nl7 PROPERTIES TIMEOUT 2500)
endif ()

# Sequential test cases
if (USE_SM)
    file (GLOB sm_tests RELATIVE "${oofem_TEST_DIR}/sm" "${oofem_TEST_DIR}/sm/*.in")
    foreach (case ${sm_tests})
        add_test (NAME "test_${case}" WORKING_DIRECTORY ${oofem_TEST_DIR}/sm COMMAND ${oofem_cmd} "-f" ${case})
    endforeach (case)
endif ()

if (USE_FM)
    file (GLOB fm_tests RELATIVE "${oofem_TEST_DIR}/fm" "${oofem_TEST_DIR}/fm/*.in")
    foreach (case ${fm_tests})
        add_test (NAME "test_${case}" WORKING_DIRECTORY ${oofem_TEST_DIR}/fm COMMAND ${oofem_cmd} "-f" ${case})
    endforeach (case)
endif()

if (USE_TM)
    file (GLOB tm_tests RELATIVE "${oofem_TEST_DIR}/tm" "${oofem_TEST_DIR}/tm/*.in")
    foreach (case ${tm_tests})
        add_test (NAME "test_${case}" WORKING_DIRECTORY ${oofem_TEST_DIR}/tm COMMAND ${oofem_cmd} "-f" ${case})
    endforeach (case)
endif()

if (USE_TM AND USE_SM)
    file (GLOB tmsm_tests RELATIVE "${oofem_TEST_DIR}/tmsm" "${oofem_TEST_DIR}/tmsm/*.in")
    foreach (case ${tmsm_tests})
        add_test (NAME "test_${case}" WORKING_DIRECTORY ${oofem_TEST_DIR}/tmsm COMMAND ${oofem_cmd} "-f" ${case})
    endforeach (case)
endif()

if (USE_TM AND USE_FM)
    file (GLOB tmfm_tests RELATIVE "${oofem_TEST_DIR}/tmfm" "${oofem_TEST_DIR}/tmfm/*.in")
    foreach (case ${tmfm_tests})
        add_test (NAME "test_${case}" WORKING_DIRECTORY ${oofem_TEST_DIR}/tmfm COMMAND ${oofem_cmd} "-f" ${case})
    endforeach (case)
endif()

if (USE_TM AND USE_CEMHYD)
    file (GLOB tmcemhyd_tests RELATIVE "${oofem_TEST_DIR}/tmcemhyd" "${oofem_TEST_DIR}/tmcemhyd/*.in")
    foreach (case ${tmcemhyd_tests})
        add_test (NAME "test_${case}" WORKING_DIRECTORY ${oofem_TEST_DIR}/tmcemhyd COMMAND ${oofem_cmd} "-f" ${case})
    endforeach (case)
endif()

if (USE_FM AND USE_PFEM)
    file (GLOB fmpfem_tests RELATIVE "${oofem_TEST_DIR}/fmpfem" "${oofem_TEST_DIR}/fmpfem/*.in")
    foreach (case ${fmpfem_tests})
        add_test (NAME "test_${case}" WORKING_DIRECTORY ${oofem_TEST_DIR}/fmpfem COMMAND ${oofem_cmd} "-f" ${case})
    endforeach (case)
endif()


######################## Benchmarks ########################################

set (oofem_BENCHMARK_DIR "${oofem_SOURCE_DIR}/tests/benchmark")

if (USE_TM AND USE_SM AND USE_DSS AND USE_IML)
    file (GLOB tmsmdssiml_benchmark RELATIVE "${oofem_BENCHMARK_DIR}/tmsmdssiml" "${oofem_BENCHMARK_DIR}/tmsmdssiml/*.in")
    foreach (case ${tmsmdssiml_benchmark})
        add_test (NAME "benchmark_${case}" WORKING_DIRECTORY ${oofem_BENCHMARK_DIR}/tmsmdssiml COMMAND ${oofem_cmd} "-f" ${case})
    endforeach (case)
endif()

if (USE_FM)
    file (GLOB fm_benchmark RELATIVE "${oofem_BENCHMARK_DIR}/fm" "${oofem_BENCHMARK_DIR}/fm/*.in")
    foreach (case ${fm_benchmark})
        add_test (NAME "benchmark_${case}" WORKING_DIRECTORY ${oofem_BENCHMARK_DIR}/fm COMMAND ${oofem_cmd} "-f" ${case})
    endforeach (case)
endif()

if (USE_SM)
    file (GLOB sm_benchmark RELATIVE "${oofem_BENCHMARK_DIR}/sm" "${oofem_BENCHMARK_DIR}/sm/*.in")
    foreach (case ${sm_benchmark})
        add_test (NAME "benchmark_${case}" WORKING_DIRECTORY ${oofem_BENCHMARK_DIR}/sm COMMAND ${oofem_cmd} "-f" ${case})
    endforeach (case)
endif()

if (USE_TM)
    file (GLOB tm_benchmark RELATIVE "${oofem_BENCHMARK_DIR}/tm" "${oofem_BENCHMARK_DIR}/tm/*.in")
    foreach (case ${tm_benchmark})
        add_test (NAME "benchmark_${case}" WORKING_DIRECTORY ${oofem_BENCHMARK_DIR}/tm COMMAND ${oofem_cmd} "-f" ${case})
    endforeach (case)
endif()

if (USE_FM AND USE_SM AND USE_PFEM)
    file (GLOB fmpfem_benchmark RELATIVE "${oofem_BENCHMARK_DIR}/fmpfem" "${oofem_BENCHMARK_DIR}/fmpfem/*.in")
    foreach (case ${fmpfem_benchmark})
        add_test (NAME "benchmark_${case}" WORKING_DIRECTORY ${oofem_BENCHMARK_DIR}/fmpfem COMMAND ${oofem_cmd} "-f" ${case})
    endforeach (case)
endif()   
  
# Custom targets for test and benchmarks (which can conveniently be accessed from IDEs)

ProcessorCount(N)
if (N EQUAL 0)
  set (N 1)
endif ()

if (MSVC)
    add_custom_target (tests      COMMAND ${CMAKE_CTEST_COMMAND} "-j${N}" "-R" "^test" "-C" "Debug")
    add_custom_target (partests   COMMAND ${CMAKE_CTEST_COMMAND} "-R" "^partest" "-C" "Debug")
    add_custom_target (benchmarks COMMAND ${CMAKE_CTEST_COMMAND} "-j${N}" "-R" "^benchmark" "-C" "Debug")
else ()
    add_custom_target (tests      COMMAND ${CMAKE_CTEST_COMMAND} "-j${N}" "-R" "^test")
    add_custom_target (partests   COMMAND ${CMAKE_CTEST_COMMAND} "-R" "^partest")
    add_custom_target (benchmarks COMMAND ${CMAKE_CTEST_COMMAND} "-j${N}" "-R" "^benchmark")
endif ()


#######################################################################
######################## Packing ######################################
#######################################################################
### versions
set (CPACK_PACKAGE_VERSION_MAJOR ${oofem_VERSION_MAJOR})
set (CPACK_PACKAGE_VERSION_MINOR ${oofem_VERSION_MINOR})
set (CPACK_PACKAGE_VERSION_PATCH ${oofem_VERSION_PATCH})
set (CPACK_PACKAGE_VERSION "${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}")

### general settings
set (CPACK_PACKAGE_NAME "oofem")
set (CPACK_PACKAGE_VENDOR "OOFEM development team")
set (CPACK_PACKAGE_DESCRIPTION_SUMMARY "Object Oriented Finite Element")
#set (CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_SOURCE_DIR}/README")
set (CPACK_RESOURCE_FILE_README "${CMAKE_SOURCE_DIR}/README")
set (CPACK_RESOURCE_FILE_LICENSE "${CMAKE_SOURCE_DIR}/COPYING.LGPLv2.1")
set (CPACK_PACKAGE_INSTALL_DIRECTORY ${CPACK_PACKAGE_NAME})
set (CPACK_PACKAGE_EXECUTABLES "oofem")
set (CPACK_PACKAGE_FILE_NAME ${CPACK_PACKAGE_NAME}_${CPACK_PACKAGE_VERSION}_${CMAKE_SYSTEM_PROCESSOR})

### default generators
if (WIN32)
    set(CPACK_GENERATOR "ZIP")
else ()
    set (CPACK_GENERATOR "TGZ;DEB")
endif ()

### source package (type "make source_package" to produce a zip)
set (CPACK_SOURCE_GENERATOR "ZIP")
set (CPACK_SOURCE_PACKAGE_FILE_NAME "oofem-${CPACK_PACKAGE_VERSION}")
set (CPACK_SOURCE_IGNORE_FILES "~$;/build/;tags;cscope.*;.*\\\\.out$;\\\\.out\\\\.;/\\\\..*;\\\\.kdev4$;do_release;release_filter\\\\.pl") # Note: All hidden files are ignored

### deb package options
set (CPACK_DEBIAN_PACKAGE_MAINTAINER "Mikael Öhman <micketeer@gmail.com>")
set (CPACK_DEBIAN_PACKAGE_SECTION "Mathematics")
# autogenerate dependency information
set (CPACK_DEBIAN_PACKAGE_SHLIBDEPS ON)

set(CPACK_DEBIAN_PACKAGE_VERSION "${CPACK_PACKAGE_VERSION}.${CPACK_PACKAGE_VERSION_PATCH}+sid1")

include (CPack)
