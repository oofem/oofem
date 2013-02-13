cmake_minimum_required (VERSION 2.8) # >= 2.8.8 necessary for $<TARGET_OBJECTS:objlib>
include (CheckSymbolExists)
include (CheckCXXSymbolExists)
include (CheckIncludeFiles)
include (CheckIncludeFileCXX)
include (CheckCXXSourceCompiles)

project (oofem)
set (oofem_VERSION_MAJOR 2)
set (oofem_VERSION_MINOR 2)
set (oofem_VERSION_PATCH 0) # Unused
set (oofem_HOMEPAGE "http://www.oofem.org")
set (bp_HOMEPAGE "http://ksm.fsv.cvut.cz/%7Ebp/bp.html")
set (oofem_COPYRIGHT "2012")
set (oofem_EMAIL "info@oofem.org")

# Default is debug build
if (NOT CMAKE_BUILD_TYPE)
	set (CMAKE_BUILD_TYPE Debug)
endif ()

if (MSVC)
    add_definitions (-D_USE_MATH_DEFINES)
    set (CMAKE_REQUIRED_DEFINITIONS "-D_USE_MATH_DEFINES")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -W4 -wd4100 -wd4996 -wd4706 -wd4512")
else ()
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -fPIC")
endif ()

set (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DDEBUG")
set (CMAKE_CXX_FLAGS_PROFILING "${CMAKE_CXX_FLAGS_PROFILING} -O2 -pg")


#######################################################################
######################## Custom paths #################################
#######################################################################
set (SLEPC_DIR CACHE PATH "Manual SLEPc directory")
set (X11_DIR CACHE PATH "Manual X11 directory")
set (TINYXML2_DIR CACHE PATH "Manual TinyXML-2 directory")
set (MY_PETSC_DIR CACHE PATH "Manual PETSc directory")
set (MY_PETSC_ARCH CACHE STRING "Manual PETSc arch name (e.g. linux-gnu-cxx-opt)")
set (MPI_DIR CACHE PATH "Manual MPI directory")
set (VTK_DIR CACHE PATH "Manual VTK directory")
set (LAPACK_DIR CACHE PATH "Manual LAPACK directory")
set (IML_DIR CACHE PATH "Manual IML++ directory")
set (SPOOLES_DIR CACHE PATH "Manual SPOOLES directory")
set (ELIXIR_DIR CACHE PATH "Manual Elixir directory (OOFEG only)")
set (CKIT_DIR CACHE PATH "Manual CKit directory (OOFEG only)")

#######################################################################
######################## Options ######################################
#######################################################################
option (USE_OOFEG "Enable OOFEG" OFF)
option (USE_PYTHON_BINDINGS "Enable Python bidings" OFF)
if (${CMAKE_VERSION} VERSION_GREATER "2.8.7") # Compability with older version of cmake
    option (USE_SHARED_LIB "Build shared library" OFF)
    option (USE_STATIC_LIB "Build static library" OFF)
endif ()
# Parallel stuff
option (USE_PARALLEL "Enable parallel support" OFF)
option (USE_PARMETIS "Enable Parmetis support" OFF)
option (USE_OPENMP "Compile with OpenMP support (for parallel assembly)" OFF)
# Solvers and such
option (USE_DSS "Enable DSS module" OFF) # No reason to use this
option (USE_IML "Enable iml++ solvers" OFF) # or this
option (USE_SPOOLES "Enable SPOOLES-solver" OFF) # or this
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
option (USE_TM "Enable heat/mass transfer module" ON)
option (USE_CEMHYD "Enable CemHyd support" OFF)

if (USE_PARALLEL)
    add_definitions (-D__PARALLEL_MODE)
    set (USE_MPI ON)
endif ()


#######################################################################
######################## Internal libraries ###########################
#######################################################################
# I believe the subdirs should be used when including instead.
include_directories ("${oofem_SOURCE_DIR}/src/oofemlib")
include_directories ("${oofem_SOURCE_DIR}/src/oofemlib/xfem")
include_directories ("${oofem_SOURCE_DIR}/src/oofemlib/iga")
include_directories ("${oofem_BINARY_DIR}")

# Is there really any reason not to compile all of them? Some has cross dependencies.
if (${CMAKE_VERSION} VERSION_LESS "2.8.8")
    set (LIBS ${LIBS} core)
else ()
    set (LIBS ${LIBS} $<TARGET_OBJECTS:core>)
endif ()

if (USE_FM)
    include_directories ("${oofem_SOURCE_DIR}/src/fm")
    add_definitions (-D__FM_MODULE)
    if (${CMAKE_VERSION} VERSION_LESS "2.8.8")
        set (LIBS ${LIBS} fm)
    else ()
        set (LIBS ${LIBS} $<TARGET_OBJECTS:fm>)
    endif ()
    set (MODULE_LIST ${MODULE_LIST} "fm")
endif ()

if (USE_TM)
    include_directories ("${oofem_SOURCE_DIR}/src/tm")
    include_directories ("${oofem_SOURCE_DIR}/src/tm/cemhyd")
    add_definitions (-D__TM_MODULE)
    if (${CMAKE_VERSION} VERSION_LESS "2.8.8")
        set (LIBS ${LIBS} tm)
    else ()
        set (LIBS ${LIBS} $<TARGET_OBJECTS:tm>)
    endif ()
    set (MODULE_LIST ${MODULE_LIST} "tm")
    if (USE_CEMHYD)
        add_definitions (-D__CEMHYD_MODULE)
        set (USE_TINYXML ON)
    endif ()
endif ()

if (USE_SM)
    include_directories ("${oofem_SOURCE_DIR}/src/sm")
    add_definitions (-D__SM_MODULE)
    if (${CMAKE_VERSION} VERSION_LESS "2.8.8")
        set (LIBS ${LIBS} sm)
    else ()
        set (LIBS ${LIBS} $<TARGET_OBJECTS:sm>)
    endif ()
    set (MODULE_LIST ${MODULE_LIST} "sm")

    if (NOT MSVC AND NOT MINGW)
        find_library (DL_LIB dl) 
        set (EXT_LIBS ${EXT_LIBS} ${DL_LIB})
    endif ()
    
#     if (MINGW)
#         find_library (DL_LIB ltdl HINTS "C:/MinGW/lib/" ) 
#         set (EXT_LIBS ${EXT_LIBS} ${DL_LIB})
#     endif ()
    
endif ()

if (USE_DSS)
    include_directories ("${oofem_SOURCE_DIR}/src/dss")
    add_definitions (-D__DSS_MODULE)
    if (${CMAKE_VERSION} VERSION_LESS "2.8.8")
        set (LIBS ${LIBS} dss)
    else ()
        set (LIBS ${LIBS} $<TARGET_OBJECTS:dss>)
    endif ()
    set (MODULE_LIST ${MODULE_LIST} "dss")
endif ()

if (USE_OPENMP)
    include (FindOpenMP)
    if (OPENMP_FOUND)
        set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
        set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
        set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
        set (MODULE_LIST ${MODULE_LIST} "openmp")
    endif ()
endif ()

if (USE_OOFEG)
    add_definitions (-D__OOFEG)

    # Necessary X11 stuff;
    if (X11_DIR)
        set (CMAKE_LIBRARY_DIR ${CMAKE_LIBRARY_DIR} ${X11_DIR})
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

    set (EXT_LIBS ${EXT_LIBS} ${ESI_LIB} ${ELIXIR_LIB} ${CKIT_LIB} ${X11_LIBRARIES} ${X11_Xmu_LIB} ${X11_Xt_LIB} ${XAW_LIB})
    set (MODULE_LIST ${MODULE_LIST} "oofeg")
endif ()

if (BUILD_TESTING)
    find_package(PythonInterp)
endif ()

#######################################################################
######################## External libaries ############################
#######################################################################
if (USE_LAPACK)
    set (CMAKE_LIBRARY_DIR ${CMAKE_LIBRARY_DIR} ${LAPACK_DIR})
    find_package (LAPACK REQUIRED)
    set (EXT_LIBS ${EXT_LIBS} ${LAPACK_LIBRARIES})
    set (MODULE_LIST ${MODULE_LIST} "LAPACK")
    add_definitions (-D__LAPACK_MODULE)
endif ()

if (USE_VTK)
    set (CMAKE_LIBRARY_DIR ${CMAKE_LIBRARY_DIR} ${VTK_DIR})
    find_package (VTK REQUIRED)
    include (${VTK_USE_FILE}) # Sets necessary includes_directories
    add_definitions (-D__VTK_MODULE)
    set (EXT_LIBS ${EXT_LIBS} ${VTK_LIBRARIES})
    set (MODULE_LIST ${MODULE_LIST} "VTK")
endif ()

if (USE_PARMETIS)
    if (PARMETIS_DIR)
        find_library (PARMETIS_LIB parmetis PATH "${PARMETIS_DIR}/lib")
        include_directories ("${PARMETIS_DIR}/include")
    endif ()
    if (${PARMETIS_LIB} STREQUAL "PARMETIS_LIB-NOTFOUND")
        message (FATAL_ERROR "Parmetis library not found")
    endif ()
    set (EXT_LIBS ${EXT_LIBS} parmetis)
    set (MODULE_LIST ${MODULE_LIST} "parmetis")
    set (USE_MPI ON)
    add_definitions (-D__PARMETIS_MODULE)
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
    set (EXT_LIBS ${EXT_LIBS} ${SLEPC_LIBS})
    set (MODULE_LIST ${MODULE_LIST} "slepc")
    set (USE_PETSC ON)
    add_definitions (-D__SLEPC_MODULE)
endif ()

if (USE_PETSC)
    find_package (PETSc REQUIRED HINTS "${MY_PETSC_DIR}/${MY_PETSC_ARCH}/conf")
    include_directories ("${MY_PETSC_DIR}/include" "${MY_PETSC_DIR}/${MY_PETSC_ARCH}/include" ${PETSC_INCLUDE_PATH})
    find_library (PETSC_LIB petsc PATH "${MY_PETSC_DIR}/${MY_PETSC_ARCH}/lib/") 
    if (${PETSC_LIB} STREQUAL "PETSC_LIB-NOTFOUND")
        message (FATAL_ERROR "PETSc library not found")
    endif ()
    set (EXT_LIBS ${EXT_LIBS} ${PETSC_LIB} ${PETSC_PACKAGE_LIBS})
    set (MODULE_LIST ${MODULE_LIST} "PETSc")
    set (USE_MPI ON)
    add_definitions (-D__PETSC_MODULE)
endif ()

if (USE_MPI)
    if (MPI_DIR)
        include_directories ("${MPI_DIR}/include")
        link_directories ("${MPI_DIR}/lib")
        find_library (MPI_LIBRARY NAMES mpi PATHS ${MPI_DIR}/lib)
        set (EXT_LIBS ${EXT_LIBS} ${MPI_LIBRARY})
        find_library (MPI_CXX_LIBRARY NAMES mpi_cxx PATHS ${MPI_DIR}/lib)
        if (MPI_CXX_LIBRARY)
            set (EXT_LIBS ${EXT_LIBS} ${MPI_CXX_LIBRARY})
        endif ()
    else ()
        include (FindMPI)
        find_package (MPI REQUIRED)
        set (CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILES_FLAGS} ${MPI_COMPILE_FLAGS})
        set (CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${MPI_LINK_FLAGS})
        add_definitions (-DMPICH_IGNORE_CXX_SEEK)
        include_directories (${MPI_INCLUDE_PATH})
        set (EXT_LIBS ${EXT_LIBS} ${MPI_LIBRARIES})
    endif (MPI_DIR)
    set (MODULE_LIST ${MODULE_LIST} "MPI")
    # Compile and link flags?
endif ()

if (USE_TRIANGLE)
    find_library (TRIANGLE_LIB triangle DOC "Triangle (library)")
    if (${TRIANGLE_LIB} STREQUAL "TRIANGLE_LIB-NOTFOUND")
        message (FATAL_ERROR "Triangle library not found")
    endif ()
    add_definitions (-D__TRIANGLE_MODULE)
    set (EXT_LIBS ${EXT_LIBS} ${TRIANGLE_LIB})
    set (MODULE_LIST ${MODULE_LIST} "Triangle")
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
    set (EXT_LIBS ${EXT_LIBS} ${TINYXML2_LIB})
    set (MODULE_LIST ${MODULE_LIST} "TinyXML-2")
endif ()

if (USE_IML)
    if (IML_DIR)
        include_directories ("${IML_DIR}")
        set (CMAKE_REQUIRED_INCLUDES "${IML_DIR}")
    else ()
        include_directories ("${oofem_SOURCE_DIR}")
        set (CMAKE_REQUIRED_INCLUDES "${oofem_SOURCE_DIR}")
    endif ()
    check_include_file_cxx ("iml/cg.h" HAVE_IML_CG_H)
    check_include_file_cxx ("iml/gmres.h" HAVE_IML_GMRES_H)
    if ((NOT HAVE_IML_CG_H) OR (NOT HAVE_IML_GMRES_H))
        message (FATAL_ERROR "Necessary IML++ headers not found")
    endif ()
    add_definitions (-D__IML_MODULE)
    set (MODULE_LIST ${MODULE_LIST} "IML++")
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
    add_definitions (-D__SPOOLES_MODULE)
    set (EXT_LIBS ${EXT_LIBS} ${SPOOLES_LIB})
    set (MODULE_LIST ${MODULE_LIST} "SPOOLES")
endif ()

if (USE_PYTHON_BINDINGS)
    # PythonLibs in cmake is bugged, have to try and rely on PythonInterp to get the suitable major and minor versions.
    find_package(PythonInterp REQUIRED)
    find_package(PythonLibs ${PYTHON_VERSION_STRING} EXACT REQUIRED)
    find_package(Boost COMPONENTS "python-py${PYTHON_VERSION_MAJOR}${PYTHON_VERSION_MINOR}" REQUIRED)
    include_directories (${PYTHON_INCLUDE_DIRS})
    include_directories (${Boost_INCLUDE_DIRS})
    include_directories (${oofem_SOURCE_DIR}/src/main) # classfactory should probably be in oofemlib instead.
    set (EXT_LIBS ${EXT_LIBS} ${PYTHON_LIBRARIES} ${Boost_LIBRARIES})
    set (MODULE_LIST ${MODULE_LIST} "Python")
    add_definitions(-DBOOST_PYTHON)
endif ()


#######################################################################
######################## Configuration ################################
#######################################################################
site_name (HOST_NAME)
set (HOST_TYPE "${CMAKE_SYSTEM_PROCESSOR}-${CMAKE_SYSTEM_NAME}") # This should be renamed TARGET_TYPE.
set (BUILDNAME "OOFEM(${MODULE_LIST})-${HOST_TYPE}-${HOST_NAME}")

check_include_files ("sys/time.h" HAVE_SYS_TIME_H)
check_include_files ("time.h" HAVE_TIME_H)
if (HAVE_TIME_H AND HAVE_SYS_TIME_H)
    set (TIME_WITH_SYS_TIME true)
endif ()
check_include_files ("strings.h" HAVE_STRINGS_H)
check_include_files ("unistd.h" HAVE_UNISTD_H)
check_include_files ("execinfo.h" HAVE_EXECINFO_H)
check_symbol_exists ("strncasecmp" "strings.h" HAVE_STRNCASECMP)
check_symbol_exists ("access" "unistd.h" HAVE_ACCESS)
check_cxx_symbol_exists ("nearest" "cmath" HAVE_NEAREST) # In what header? I don't have it.
check_cxx_symbol_exists ("cbrt" "cmath" HAVE_CBRT)
check_cxx_symbol_exists ("isnan" "cmath" HAVE_ISNAN)
check_cxx_symbol_exists ("M_PI" "cmath" HAVE_M_PI)
check_cxx_symbol_exists ("M_LN2" "cmath" HAVE_M_LN2)
check_cxx_source_compiles ("
#include <stdio.h>
#define LOG_FATAL(file,...) fprintf (file, __VA_ARGS__)
int main(void) { return 0; }" HAVE_MACRO_VA_ARGS)

configure_file (
  "${oofem_SOURCE_DIR}/src/oofemcfg.h.in"
  "${oofem_BINARY_DIR}/oofemcfg.h"
  )


#######################################################################
######################## Targets ######################################
#######################################################################
# Subtargets
if (${CMAKE_VERSION} VERSION_LESS "2.8.8")
    set (LIB_TYPE STATIC)
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
if (${CMAKE_VERSION} VERSION_LESS "2.8.8")
    set (USE_SHARED_LIB OFF)
endif ()
if (USE_SHARED_LIB)
    add_library (liboofem
        SHARED
        # These doesn't really belong in main;
        ${oofem_SOURCE_DIR}/src/main/usrdefsub.C
        ${oofem_SOURCE_DIR}/src/main/classfactory.C
        ${oofempy}
        ${LIBS} 
        )
    target_link_libraries (liboofem ${EXT_LIBS})
    set_target_properties (liboofem PROPERTIES OUTPUT_NAME oofem)
    install (TARGETS liboofem DESTINATION lib)
else ()
    set (USE_STATIC_LIB ON)
endif ()

# Optional static lib (only compiled and installed if explicitly called):
if (USE_STATIC_LIB)
    if (${CMAKE_VERSION} VERSION_LESS "2.8.8")
        add_library (liboofem_static
            STATIC
            # These doesn't really belong in main;
            ${oofem_SOURCE_DIR}/src/main/usrdefsub.C
            ${oofem_SOURCE_DIR}/src/main/classfactory.C
            ${oofempy}
            )
        target_link_libraries (liboofem_static ${LIBS} ${LIBS} ${EXT_LIBS})
    else ()
        add_library (liboofem_static
            STATIC
            # These doesn't really belong in main;
            ${oofem_SOURCE_DIR}/src/main/usrdefsub.C
            ${oofem_SOURCE_DIR}/src/main/classfactory.C
            ${oofempy}
            ${LIBS} 
            )
        target_link_libraries (liboofem_static ${EXT_LIBS})
    endif ()
    set_target_properties (liboofem_static PROPERTIES OUTPUT_NAME oofem)
    install (TARGETS liboofem_static DESTINATION lib)
endif ()

# OOFEM target;
add_executable (oofem ${oofem_SOURCE_DIR}/src/main/main.C)
if (USE_SHARED_LIB) # Prefer dynamic lib if available.
    target_link_libraries (oofem liboofem)
else ()
    if (${CMAKE_VERSION} VERSION_LESS "2.8.8")
        target_link_libraries (oofem ${LIBS} liboofem_static ${LIBS} ${LIBS}) # Cyclic dependencies.
    else ()
        target_link_libraries (oofem liboofem_static)
    endif ()
endif ()
install (TARGETS oofem DESTINATION bin)

if (USE_OOFEG)
    #  OOFEG target:
    add_executable (oofeg ${oofem_SOURCE_DIR}/src/main/oofeg.C)
    if (USE_SHARED_LIB)
        target_link_libraries (oofeg liboofem)
    else ()
        target_link_libraries (oofeg liboofem_static liboofem_static)
    endif ()
    install(TARGETS oofeg DESTINATION bin)
endif ()

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
# This isn't perfect/elegant, partly due to the how tests work.
# TODO: Would like to run oofem from the binary dir (avoids creating *.out files in the source dir)
# unfortunately, that doesn't work for some tests that require sub-files using relative paths.
enable_testing ()
include (CTest)

set (python_path "python") # Check this
set (oofem_TEST_DIR "${oofem_SOURCE_DIR}/tests/")
set (oofem_cmd  $<TARGET_FILE:oofem>)
set (test_cmd  ${PYTHON_EXECUTABLE} "${CMAKE_CURRENT_BINARY_DIR}/testhelper.py")

# Helper script (perhaps one could do without this one)
# Some input files requier other files, so I didn't feel that there were any option to just running them in the main tree.
file (WRITE ${CMAKE_CURRENT_BINARY_DIR}/testhelper.py
"import subprocess,os,sys,glob
os.chdir(sys.argv[2])
if sys.argv[4] == 'p':
    files = glob.glob(sys.argv[1]+'.[0-9]*')
    mpicommand = 'mpirun -np '+str(len(files)) + ' '
    parallel = ' -p'
else:
    files = [sys.argv[1]]
    mpicommand = ''
    parallel = ''
if len(sys.argv) == 6:
    flags = ' ' + sys.argv[5]
else:
    flags = ''
command = mpicommand + sys.argv[3] + parallel + ' -f ' + sys.argv[1] + ' -l 2' + flags 
print('Executing: ' + command)
if os.system(command):
	exit(1)
else:
        for f in files:
	    python=\"${PYTHON_EXECUTABLE}\"
            extractor = \"${oofem_SOURCE_DIR}/tools/extractor.py\"
            print('Executing: ', [python, extractor, '-f', f, '-c'])
	    if subprocess.call([python, extractor, '-f', f, '-c']):
		exit(2)
exit(0)
")


if (USE_PARALLEL)
    set (par_dir ${oofem_TEST_DIR}/partests/)
    add_test (NAME partest_4test_np4 COMMAND ${test_cmd} 4test.oofem.in ${par_dir}/4test "${oofem_cmd}" "p" "-ksp_rtol 1.e-15")
    add_test (NAME partest_bar_np4 COMMAND ${test_cmd} bar.oofem.in ${par_dir}/bar "${oofem_cmd}" "p")
    add_test (NAME partest_barnl_np4 COMMAND ${test_cmd} barnl.oofem.in ${par_dir}/barnl "${oofem_cmd}" "p")
    add_test (NAME partest_barnl_np2 COMMAND ${test_cmd} barnl2.oofem.in ${par_dir}/barnl "${oofem_cmd}" "p")
    add_test (NAME partest_brazil_np2 COMMAND ${test_cmd} brazil_2d_nl2.oofem.in ${par_dir}/brazil "${oofem_cmd}" "p")
    add_test (NAME partest_brazil_np4 COMMAND ${test_cmd} brazil_2d_nl4.oofem.in ${par_dir}/brazil "${oofem_cmd}" "p")
    add_test (NAME partest_brazil_np7 COMMAND ${test_cmd} brazil_2d_nl7.oofem.in ${par_dir}/brazil "${oofem_cmd}" "p")
    add_test (NAME partest_brazil3d_np7 COMMAND ${test_cmd} brazil_3d_nl.oofem.in ${par_dir}/brazil3d "${oofem_cmd}" "p" "-ksp_rtol 1.e-15 -ksp_type cg")
    add_test (NAME partest_dyn_bar01_np2 COMMAND ${test_cmd} bar01.oofem.in ${par_dir}/dyn_bar01 "${oofem_cmd}" "p")
    add_test (NAME partest_dyn_bar02_np2 COMMAND ${test_cmd} bar02.oofem.in ${par_dir}/dyn_bar02 "${oofem_cmd}" "p")
    add_test (NAME partest_dyn_bar03_np2 COMMAND ${test_cmd} bar03.oofem.in ${par_dir}/dyn_bar03 "${oofem_cmd}" "p")
    # TODO: The tests requires different commands, and the files and named inconsistently. Some slight changes to naming could change this.
    # The flags are a bit trickier. It needs to be built into the input file.
    #file (GLOB parallel_tests RELATIVE "${par_dir}" "${par_dir}/*/")
    #foreach (case ${parallel_tests})
    #    add_test ("partest_${case}" ${test_cmd} ${case}.oofem.in "${par_dir}/${case}" ${oofem_cmd} "p")
    #endforeach (case)
endif ()
# Sequential test cases
file (GLOB sm_tests RELATIVE "${oofem_TEST_DIR}/sm" "${oofem_TEST_DIR}/sm/*.in")
file (GLOB fm_tests RELATIVE "${oofem_TEST_DIR}/fm" "${oofem_TEST_DIR}/fm/*.in")
file (GLOB tm_tests RELATIVE "${oofem_TEST_DIR}/tm" "${oofem_TEST_DIR}/tm/*.in")
file (GLOB tmsm_tests RELATIVE "${oofem_TEST_DIR}/tmsm" "${oofem_TEST_DIR}/tmsm/*.in")
foreach (case ${sm_tests})
    add_test (NAME "test_${case}" COMMAND ${test_cmd} ${case} "${oofem_TEST_DIR}/sm" ${oofem_cmd} "s")
endforeach (case)
foreach (case ${fm_tests})
    add_test (NAME "test_${case}" COMMAND ${test_cmd} ${case} "${oofem_TEST_DIR}/fm" ${oofem_cmd} "s")
endforeach (case)
foreach (case ${tm_tests})
    add_test (NAME "test_${case}" COMMAND ${test_cmd} ${case} "${oofem_TEST_DIR}/tm" ${oofem_cmd} "s")
endforeach (case)
foreach (case ${tmsm_tests})
    add_test (NAME "test_${case}" COMMAND ${test_cmd} ${case} "${oofem_TEST_DIR}/tmsm" ${oofem_cmd} "s")
endforeach (case)

# Benchmarks
file (GLOB benchmarks RELATIVE "${oofem_TEST_DIR}/benchmark" "${oofem_TEST_DIR}/benchmark/*.in")
foreach (case ${benchmarks})
    add_test (NAME "benchmark_${case}" COMMAND ${test_cmd} ${case} "${oofem_TEST_DIR}/benchmark" ${oofem_cmd} "s")
endforeach (case)


#file (GLOB par_tests RELATIVE "${oofem_TEST_DIR}/tests/partests" "${oofem_SOURCE_DIR}/tests/tm/*.in.*")
#runandcheck ("${oofem_TEST_DIR}/sm" ${case})
#add_test (test-${case} ${CMAKE_COMMAND} -E runcheckclean ${oofem_TEST_DIR} "sm" ${case} )
#add_test (test-${case} "${EXECUTABLE_OUTPUT_PATH}/oofem" "-f" ${case} "-l" "2")
#add_test (test-${case}-check "python" "${oofem_SOURCE_DIR}/tools/extractor.py" "-f" ${case} "-c")

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
set (CPACK_RESOURCE_FILE_LICENSE "${CMAKE_SOURCE_DIR}/gpl.txt")
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
set (CPACK_SOURCE_IGNORE_FILES "~$;/build/;tags;cscope.*;.*\\\\.out$;/\\\\..*;\\\\.kdev4$;do_release;release_filter\\\\.pl") # Note: All hidden files are ignored

### deb package options
set (CPACK_DEBIAN_PACKAGE_MAINTAINER "Mikael Öhman <micketeer@gmail.com>")
set (CPACK_DEBIAN_PACKAGE_SECTION "Mathematics")
if (CMAKE_SYSTEM_PROCESSOR STREQUAL "x86_64")
    set(CPACK_DEBIAN_PACKAGE_ARCHITECTURE "amd64")
else ()
    set(CPACK_DEBIAN_PACKAGE_ARCHITECTURE "i386") # Crude assumption. (Why can't cmake figure this out itself?)
endif ()
# autogenerate dependency information
set (CPACK_DEBIAN_PACKAGE_SHLIBDEPS ON)
# Comma separated list here is really bothersome.
#set (CPACK_DEBIAN_PACKAGE_DEPENDS "libpetsc3.3 (>=3.3), libslepc3.3 (>=3.3), libparmetis3.1 (>=3.1), libvtk5.8 (>=5.8), liblapack3 (>=3.2), libtriangle-1.6 (>=1.6)")

set(CPACK_DEBIAN_PACKAGE_VERSION "${CPACK_PACKAGE_VERSION}.${CPACK_PACKAGE_VERSION_PATCH}+sid1")

include (CPack)
