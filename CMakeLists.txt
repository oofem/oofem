cmake_minimum_required (VERSION 2.6)
include (CheckSymbolExists)
include (CheckCXXSymbolExists)
include (CheckIncludeFiles)
include (CheckIncludeFileCXX)
include (CheckCXXSourceCompiles)
include (FindPythonInterp)

project (oofem)
set (oofem_VERSION_MAJOR 2)
set (oofem_VERSION_MINOR 2)
set (oofem_VERSION_PATCH 0) # Unused
set (oofem_HOMEPAGE "http://www.oofem.org")
set (bp_HOMEPAGE "http://ksm.fsv.cvut.cz/%7Ebp/bp.html")
set (oofem_COPYRIGHT "2012")
set (oofem_EMAIL "info@oofem.org")

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
set (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DDEBUG")
#set (CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O2 -fPIE")
set (CMAKE_CXX_FLAGS_PROFILING "${CMAKE_CXX_FLAGS_PROFILING} -O2 -pg")

# Default is debug build
if(NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE Debug)
endif(NOT CMAKE_BUILD_TYPE)

if (MSVC)
    add_definitions (-D_USE_MATH_DEFINES)
    set (CMAKE_REQUIRED_DEFINITIONS "-D_USE_MATH_DEFINES")
endif ()

#######################################################################
######################## Custom paths #################################
#######################################################################
set (SLEPC_DIR CACHE PATH "Manual SLEPc directory")
set (MY_PETSC_DIR CACHE PATH "Manual PETSc directory")
set (MY_PETSC_ARCH CACHE STRING "Manual PETSc arch name (e.g. linux-gnu-cxx-opt)")
set (MPI_DIR CACHE PATH "Manual MPI directory")
set (IML_DIR CACHE PATH "Manual IML++ directory")
set (SPOOLES_DIR CACHE PATH "Manual SPOOLES directory")

#######################################################################
######################## Options ######################################
#######################################################################
option (USE_OOFEG "Enable OOFEG" OFF)
option (__ENABLE_COMPONENT_LABELS "Enable component labels" ON)
# Parallell stuff
option (USE_PARALLEL "Enable parallel support" OFF)
option (USE_PARMETIS "Enable Parmetis support" OFF)
option (USE_OPENMP "Compile with OpenMP support (for parallell assembly)" OFF)
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

if (USE_OOFEG)
    add_definitions (-D__OOFEG)
endif (USE_OOFEG)

if (USE_PARALLEL)
    add_definitions (-D__PARALLEL_MODE)
    set (USE_MPI ON)
endif (USE_PARALLEL)

#######################################################################
######################## Internal libraries ###########################
#######################################################################
# I believe the subdirs should be used when including instead.
include_directories ("${oofem_SOURCE_DIR}/src/oofemlib")
include_directories ("${oofem_SOURCE_DIR}/src/oofemlib/xfem")
include_directories ("${oofem_SOURCE_DIR}/src/oofemlib/iga")
include_directories ("${oofem_SOURCE_DIR}/src")
include_directories ("${oofem_BINARY_DIR}")

# Is there really any reason not to compile all of them? Some has cross dependencies.
if (USE_FM)
    include_directories ("${oofem_SOURCE_DIR}/src/fm")
    add_definitions (-D__FM_MODULE)
    set (LIBS ${LIBS} fm)
    set (MODULE_LIST ${MODULE_LIST} fm)
endif (USE_FM)

if (USE_TM)
    include_directories ("${oofem_SOURCE_DIR}/src/tm")
    include_directories ("${oofem_SOURCE_DIR}/src/tm/cemhyd")
    add_definitions (-D__TM_MODULE)
    set (LIBS ${LIBS} tm)
    set (MODULE_LIST ${MODULE_LIST} tm)
    if (USE_CEMHYD)
        add_definitions (-D__CEMHYD_MODULE)
        set (USE_TINYXML ON)
    endif (USE_CEMHYD)
endif (USE_TM)

if (USE_SM)
    include_directories ("${oofem_SOURCE_DIR}/src/sm")
    add_definitions (-D__SM_MODULE)
    set (LIBS ${LIBS} sm)
    set (MODULE_LIST ${MODULE_LIST} sm)
endif (USE_SM)

set (LIBS ${LIBS} oofemlib)

if (USE_DSS)
    include_directories ("${oofem_SOURCE_DIR}/src/dss")
    add_definitions (-D__DSS_MODULE)
    set (LIBS ${LIBS} dss)
    set (MODULE_LIST ${MODULE_LIST} dss)
endif (USE_DSS)

if (USE_OPENMP)
    include (FindOpenMP)
    if(OPENMP_FOUND)
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
    endif()
endif (USE_OPENMP)

#######################################################################
######################## External libaries ############################
#######################################################################
if (USE_LAPACK)
    find_package (LAPACK REQUIRED)
    set (LIBS ${LIBS} ${LAPACK_LIBRARIES})
    set (MODULE_LIST ${MODULE_LIST} lapack)
    add_definitions (-D__LAPACK_MODULE)
endif (USE_LAPACK)

if (USE_VTK)
    find_package(VTK REQUIRED)
    include (${VTK_USE_FILE}) # Sets necessary includes_directories
    add_definitions (-D__VTK_MODULE)
    #include_directories (${VTK_INCLUDE_DIRS})
    set (LIBS ${LIBS} ${VTK_LIBRARIES})
    set (MODULE_LIST ${MODULE_LIST} parmetis)
endif (USE_VTK)

if (USE_PARMETIS)
    if (PARMETIS_DIR)
        find_library (PARMETIS_LIB  NAMES parmetis PATH "${PARMETIS_DIR}/lib")
        include_directories ("${PARMETIS_DIR}/include")
    endif ()
    if (${PARMETIS_LIB} STREQUAL "parmetis-NOTFOUND")
        message (FATAL_ERROR "Parmetis library not found")
    endif ()
    set (LIBS ${LIBS} parmetis)
    set (MODULE_LIST ${MODULE_LIST} parmetis)
    set (USE_MPI ON)
    add_definitions (-D__PARMETIS_MODULE)
endif (USE_PARMETIS)

if (USE_SLEPC)
    if (SLEPC_DIR)
        find_library (SLEPC_LIBS NAMES slepc PATH "${SLEPC_DIR}/lib")
        include_directories ("${SLEPC_DIR}/include")
    else ()
        find_library (SLEPC_LIBS NAMES slepc)
    endif ()
    if (${SLEPC_LIB} STREQUAL "slepc-NOTFOUND")
        message (FATAL_ERROR "SLEPc library not found")
    endif ()
    set (LIBS ${LIBS} ${SLEPC_LIBS})
    set (MODULE_LIST ${MODULE_LIST} slepc)
    set (USE_PETSC ON)
    add_definitions (-D__SLEPC_MODULE)
endif (USE_SLEPC)

if (USE_PETSC)
    find_package (PETSc REQUIRED HINTS "${MY_PETSC_DIR}/${MY_PETSC_ARCH}/conf")
    include_directories ("${MY_PETSC_DIR}/include" "${MY_PETSC_DIR}/${MY_PETSC_ARCH}/include" ${PETSC_INCLUDE_PATH})
    set (LIBS ${LIBS} ${PETSC_PACKAGE_LIBS} "${MY_PETSC_DIR}/${MY_PETSC_ARCH}/lib/libpetsc.so")
    set (MODULE_LIST ${MODULE_LIST} PETSc)
    set (USE_MPI ON)
    add_definitions (-D__PETSC_MODULE)
endif (USE_PETSC)

if (USE_MPI)
    if (MPI_DIR)
        include_directories ("${MPI_DIR}/include")
        link_directories ("${MPI_DIR}/lib")
        set (LIBS ${LIBS} mpi_cxx)
    else ()
        include (FindMPI)
        find_package (MPI REQUIRED)
        set (CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILES_FLAGS} ${MPI_COMPILE_FLAGS})
        set (CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${MPI_LINK_FLAGS})
        add_definitions(-DMPICH_IGNORE_CXX_SEEK)
        include_directories (${MPI_INCLUDE_PATH})
        set (LIBS ${LIBS} ${MPI_LIBRARIES})
    endif (MPI_DIR)
    set (MODULE_LIST ${MODULE_LIST} MPI)
    # Compile and link flags?
endif (USE_MPI)

if (USE_TRIANGLE)
    find_library (TRIANGLE_LIB triangle DOC "Triangle (library)")
    if (${TRIANGLE_LIB} STREQUAL "triangle-NOTFOUND")
        message (FATAL_ERROR "Triangle not found")
    endif ()
    add_definitions (-D__TRIANGLE_MODULE)
    set (LIBS ${LIBS} ${TRIANGLE_LIB})
    set (MODULE_LIST ${MODULE_LIST} Triangle)
endif (USE_TRIANGLE)

if (USE_TINYXML)
    find_library (TINYXML_LIB tinyxml DOC "TinyXML library")
    if (${TINYXML_LIB} STREQUAL "tinyxml-NOTFOUND")
        message (FATAL_ERROR "TinyXML not found")
    endif ()
    add_definitions (-D__TINYXML_MODULE -DTIXML_USE_STL)
    set (LIBS ${LIBS} ${TINYXML_LIB})
    set (MODULE_LIST ${MODULE_LIST} TinyXML)
endif (USE_TINYXML)

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
    set (MODULE_LIST ${MODULE_LIST} iml++)
endif (USE_IML)

if (USE_SPOOLES)
    if (SPOOLES_DIR)
        find_library (SPOOLES_LIB spooles DOC "SPOOLES library" PATH "${SPOOLES_DIR/lib")
        include_directories ("${SPOOLES_DIR}/include")
    else ()
        find_library (SPOOLES_LIB spooles DOC "SPOOLES library")
    endif ()
    if (${SPOOLES_LIB} STREQUAL "spooles-NOTFOUND")
        message (FATAL_ERROR "SPOOLES not found")
    endif ()
    add_definitions (-D__SPOOLES_MODULE)
    set (LIBS ${LIBS} ${SPOOLES_LIB})
    set (MODULE_LIST ${MODULE_LIST} SPOOLES)
endif (USE_SPOOLES)

#######################################################################
######################## Configuration ################################
#######################################################################
site_name (HOST_NAME)
set (HOST_TYPE "${CMAKE_SYSTEM_PROCESSOR}-${CMAKE_SYSTEM_NAME}") # This should be renamed TARGET_TYPE.
check_include_files ("sys/time.h" HAVE_SYS_TIME_H)
check_include_files ("time.h" HAVE_TIME_H)
if (HAVE_TIME_H AND HAVE_SYS_TIME_H)
    set (TIME_WITH_SYS_TIME true)
endif ()
check_include_files ("strings.h" HAVE_STRINGS_H)
check_include_files ("unistd.h" HAVE_UNISTD_H)
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
add_subdirectory ("${oofem_SOURCE_DIR}/src/oofemlib")
if (USE_DSS)
    add_subdirectory ("${oofem_SOURCE_DIR}/src/dss")
endif (USE_DSS)
if (USE_SM)
    add_subdirectory ("${oofem_SOURCE_DIR}/src/sm")
endif (USE_SM)
if (USE_TM)
    add_subdirectory ("${oofem_SOURCE_DIR}/src/tm")
endif (USE_TM)
if (USE_FM)
    add_subdirectory ("${oofem_SOURCE_DIR}/src/fm")
endif (USE_FM)

# Main targets
# There should probably be a library + main program.
add_executable (oofem
    ${oofem_SOURCE_DIR}/src/main/main.C
    ${oofem_SOURCE_DIR}/src/main/usrdefsub.C
    ${oofem_SOURCE_DIR}/src/main/classfactory.C
    )
target_link_libraries (oofem ${LIBS})
install(TARGETS oofem RUNTIME DESTINATION bin)

if (USE_OOFEG)
    add_executable (oofeg ${oofem_SOURCE_DIR}/src/main/oofeg.C)
    target_link_libraries (oofeg ${LIBS})
    install(TARGETS oofeg RUNTIME DESTINATION bin)
endif ()

# With install paths

#######################################################################
######################## Documentation ################################
#######################################################################
# Not done yet, and i'm not sure if it's suitable.
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

add_custom_target (
    elementlibmanual
    COMMAND pdflatex elgmentlibmanual --output=${oofem_BINARY_DIR}/doc/
    COMMAND pdflatex elementlibmanual --output=${oofem_BINARY_DIR}/doc/
    #COMMAND tex4ht elementlibmanual --output=${oofem_BINARY_DIR}/doc/
    WORKING_DIRECTORY ${oofem_SOURCE_DIR}/doc/elementlibmanual/
    )

#add_custom_command(documentation command "cd doc; make pdf")
#add_custom_command(documentation command "cd doc/elementlibmanual; make pdf")

#######################################################################
######################## Tests ########################################
#######################################################################
# This isn't perfect/elegant, partly due to the how tests work.
# TODO: Would like to run oofem from the binary dir
# TODO: Inconsistent naming of parallel tests are cumbersome. Trivial to fix.
enable_testing ()
include (CTest)

set (python_path "python") # Check this
set (oofem_TEST_DIR "${oofem_SOURCE_DIR}/tests/")
set (oofem_cmd  $<TARGET_FILE:oofem>)
set (test_cmd  ${PYTHON_EXECUTABLE} "${CMAKE_CURRENT_BINARY_DIR}/testhelper.py")

# Helper script (perhaps one could do without this one)
# Some input files requier other files, so I didn't feel that there were any option to just running them in the main tree.
file (WRITE ${CMAKE_CURRENT_BINARY_DIR}/testhelper.py
"import os,sys,glob
os.chdir(sys.argv[2])
if sys.argv[4] == 'p':
    files = glob.glob(sys.argv[1]+'.[0-9]*')
    mpicommand = 'mpirun -np '+str(len(files)) + ' '
else:
    files = [sys.argv[1]]
    mpicommand = ''
if os.system(mpicommand + sys.argv[3] +' -f '+sys.argv[1]+' -l 2'):
    exit(1)
else:
    for f in files:
        if os.system('${PYTHON_EXECUTABLE} ${oofem_SOURCE_DIR}/tools/extractor.py -f '+f+' -c'):
            exit(2)
exit(0)
")

if (USE_PARALLEL)
    set (par_dir ${oofem_TEST_DIR}/partests/)
    set (parallel_tests 4test bar barnl brazil brazil3d dyn_bar01 dyn_bar02 dyn_bar03)
    #file (GLOB parallel_tests RELATIVE "${par_dir}" "${par_dir}/*/")
    add_test (NAME 4test_test COMMAND ${test_cmd} 4test.oofem.in ${par_dir}/4test "${oofem_cmd} -ksp_rtol 1.e-15" "p")
    add_test (NAME bar_test COMMAND ${test_cmd} bar.oofem.in ${par_dir}/bar "${oofem_cmd}" "p")
    add_test (NAME barnl_test COMMAND ${test_cmd} barnl.oofem.in ${par_dir}/barnl "${oofem_cmd}" "p")
    add_test (NAME brazil_test COMMAND ${test_cmd} brazil_2d_nl2.oofem.in ${par_dir}/brazil "${oofem_cmd}" "p")
    add_test (NAME brazil3d_test COMMAND ${test_cmd} brazil_3d_nl.oofem.in ${par_dir}/brazil3d "${oofem_cmd} -ksp_rtol 1.e-15 -ksp_type cg" "p")
    add_test (NAME dyn_bar01_test COMMAND ${test_cmd} bar01.oofem.in ${par_dir}/dyn_bar01 "${oofem_cmd}" "p")
    add_test (NAME dyn_bar02_test COMMAND ${test_cmd} bar02.oofem.in ${par_dir}/dyn_bar02 "${oofem_cmd}" "p")
    add_test (NAME dyn_bar03_test COMMAND ${test_cmd} bar03.oofem.in ${par_dir}/dyn_bar03 "${oofem_cmd}" "p")
    # The tests requires different commands, and the files and named inconsistently. Some slight changes to naming could change this.
    #foreach (case ${parallel_tests})
    #    add_test (${case}_test ${test_cmd} ${case}.oofem.in ${par_DIR}/${case} ${oofem_cmd} "p")
    #endforeach (case)

else () # Non-parallel test cases
    file (GLOB sm_tests RELATIVE "${oofem_TEST_DIR}/sm" "${oofem_TEST_DIR}/sm/*.in")
    file (GLOB fm_tests RELATIVE "${oofem_TEST_DIR}/fm" "${oofem_TEST_DIR}/fm/*.in")
    file (GLOB tm_tests RELATIVE "${oofem_TEST_DIR}/tm" "${oofem_TEST_DIR}/tm/*.in")
    file (GLOB tmsm_tests RELATIVE "${oofem_TEST_DIR}/tmsm" "${oofem_TEST_DIR}/tmsm/*.in")
    #file (GLOB tm_tests RELATIVE "${oofem_TEST_DIR}/partests/" "${oofem_TEST_DIR}/tm/*.in") # should run the parallel tests as well.
    foreach (case ${sm_tests})
        add_test (NAME ${case}_test COMMAND ${test_cmd} ${case} ${oofem_TEST_DIR}/sm ${oofem_cmd} "s")
    endforeach (case)
    foreach (case ${fm_tests})
        add_test (NAME ${case}_test COMMAND ${test_cmd} ${case} ${oofem_TEST_DIR}/fm ${oofem_cmd} "s")
    endforeach (case)
    foreach (case ${tm_tests})
        add_test (NAME ${case}_test COMMAND ${test_cmd} ${case} ${oofem_TEST_DIR}/tm ${oofem_cmd} "s")
    endforeach (case)
    foreach (case ${tmsm_tests})
        add_test (NAME ${case}_test COMMAND ${test_cmd} ${case} ${oofem_TEST_DIR}/tmsm ${oofem_cmd} "s")
    endforeach (case)
endif ()

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
set (CPACK_PACKAGE_VENDOR "OOFEM developement team")
set (CPACK_PACKAGE_DESCRIPTION_SUMMARY "Object Oriented Finite Element")
#set (CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_SOURCE_DIR}/README")
set (CPACK_RESOURCE_FILE_README "${CMAKE_SOURCE_DIR}/README")
set (CPACK_RESOURCE_FILE_LICENSE "${CMAKE_SOURCE_DIR}/gpl.txt")
set (CPACK_PACKAGE_INSTALL_DIRECTORY ${CPACK_PACKAGE_NAME})
set (CPACK_PACKAGE_EXECUTABLES "oofem")
set (CPACK_PACKAGE_FILE_NAME ${CPACK_PACKAGE_NAME}_${CPACK_PACKAGE_VERSION}_${CMAKE_SYSTEM_PROCESSOR})

### generators (TODO: Not sure if we need to specify this)
if (WIN32)
    set(CPACK_GENERATOR "ZIP")
else ()
    set (CPACK_GENERATOR "TGZ;DEB")
endif (WIN32)

### source package (type "make source_package" to produce a zip)
set(CPACK_SOURCE_GENERATOR "ZIP")
set(CPACK_SOURCE_PACKAGE_FILE_NAME "oofem-${CPACK_PACKAGE_VERSION}")
set(CPACK_SOURCE_IGNORE_FILES "*~;~$;[.]sw$;/[.]svn/;/[.]git/;.gitignore;/build/;tags;cscope.*;*.out")

### deb package options
set(CPACK_DEBIAN_PACKAGE_MAINTAINER "Mikael Öhman <micketeer@gmail.com>")
set(CPACK_DEBIAN_PACKAGE_SECTION "Mathematics")
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

include(CPack)
