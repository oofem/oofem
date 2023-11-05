/*
 *
 *                 #####    #####   ######  ######  ###   ###
 *               ##   ##  ##   ##  ##      ##      ## ### ##
 *              ##   ##  ##   ##  ####    ####    ##  #  ##
 *             ##   ##  ##   ##  ##      ##      ##     ##
 *            ##   ##  ##   ##  ##      ##      ##     ##
 *            #####    #####   ##      ######  ##     ##
 *
 *
 *             OOFEM : Object Oriented Finite Element Code
 *
 *               Copyright (C) 1993 - 2013   Borek Patzak
 *
 *
 *
 *       Czech Technical University, Faculty of Civil Engineering,
 *   Department of Structural Mechanics, 166 29 Prague, Czech Republic
 *
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with this library; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */

//  MAIN
//  Solves finite element problems.
//
#ifdef _PYTHON_EXTENSION
 #include <Python.h>
#endif

#include "engngm.h"
#include "oofemcfg.h"

#include "oofemtxtdatareader.h"
#include "datastream.h"
#include "util.h"
#include "error.h"
#include "logger.h"
#include "contextioerr.h"
#include "oofem_terminate.h"

#ifdef __PARALLEL_MODE
 #include "dyncombuff.h"
#endif

#ifdef __PETSC_MODULE
 #include <petsc.h>
#endif

#ifdef __SLEPC_MODULE
 #include <slepceps.h>
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef __OOFEG
 #include "oofeggraphiccontext.h"
#endif

#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <new>
#include <sstream>
// For passing PETSc/SLEPc arguments.
#include <fstream>
#include <iterator>

#include "classfactory.h"

using namespace oofem;

// debug
void oofem_debug(EngngModel &emodel);

void oofem_print_help();
void oofem_print_version();
void oofem_print_epilog();

// Finalize PETSc, SLEPc and MPI
void oofem_finalize_modules();

#define LOG_ERR_HEADER "_______________________________________________________"
#define LOG_ERR_TAIL   "_______________________________________________________\a\n"


// Handler for uncaught exceptions
void exception_handler() {
    try {
        auto eptr = std::current_exception();
        if (eptr) {
            std::rethrow_exception(eptr);
        }
    } catch (const RuntimeException& e) {

        fprintf(stderr, "%s\nOOFEM Error exception: %s\n%s", LOG_ERR_HEADER, e.what(), LOG_ERR_TAIL);
    #ifdef __GNUC__
        print_stacktrace();
    #endif
        oofem_logger.incrementErrorCounter();
        oofem_logger.printStatistics();

    } catch(const std::exception& e) {
        fprintf(stderr, "Caught exception: %s\n", e.what());
#ifdef __GNUC__
        print_stacktrace();
#endif
        exit(1);
    }
}


int main(int argc, char *argv[])
{
    // Stack trace on uncaught exceptions;
    std::set_terminate( exception_handler );

    int adaptiveRestartFlag = 0, restartStep = 0;
    bool parallelFlag = false, renumberFlag = false, debugFlag = false, contextFlag = false, restartFlag = false,
         inputFileFlag = false, outputFileFlag = false, errOutputFileFlag = false;
    std :: stringstream inputFileName, outputFileName, errOutputFileName;
    std :: vector< const char * >modulesArgs;

    int rank = 0;

#ifdef __PARALLEL_MODE
 #ifdef __USE_MPI
    MPI_Init(& argc, & argv);
    MPI_Comm_rank(MPI_COMM_WORLD, & rank);
    oofem_logger.setComm(MPI_COMM_WORLD);
 #endif
#endif

    //
    // check for options
    //
    if ( argc != 1 ) {
        // argv[0] is not read by PETSc and SLEPc.
        modulesArgs.push_back(argv [ 0 ]);
        for ( int i = 1; i < argc; i++ ) {
            if ( ( strcmp(argv [ i ], "-context") == 0 ) || ( strcmp(argv [ i ], "-c") == 0 ) ) {
                contextFlag = true;
            } else if ( strcmp(argv [ i ], "-v") == 0 ) {
                if ( rank == 0 ) {
                    oofem_print_version();
                }

                if ( argc == 2 ) {
#ifdef __USE_MPI
                    MPI_Finalize();
#endif
                    exit(EXIT_SUCCESS);     // exit if only "-v" option
                }
            } else if ( strcmp(argv [ i ], "-f") == 0 ) {
                if ( i + 1 < argc ) {
                    i++;
                    inputFileName << argv [ i ];
                    inputFileFlag = true;
                }
            } else if ( strcmp(argv [ i ], "-r") == 0 ) {
                if ( i + 1 < argc ) {
                    i++;
                    restartFlag = true;
                    restartStep = strtol(argv [ i ], NULL, 10);
                }
            } else if ( strcmp(argv [ i ], "-rn") == 0 ) {
                renumberFlag = true;
            } else if ( strcmp(argv [ i ], "-ar") == 0 ) {
                if ( i + 1 < argc ) {
                    i++;
                    adaptiveRestartFlag = strtol(argv [ i ], NULL, 10);
                }
            } else if ( strcmp(argv [ i ], "-l") == 0 ) {
                if ( i + 1 < argc ) {
                    i++;
                    int level = strtol(argv [ i ], NULL, 10);
                    oofem_logger.setLogLevel(level);
                }
            } else if ( strcmp(argv [ i ], "-qe") == 0 ) {
                if ( i + 1 < argc ) {
                    i++;
                    errOutputFileFlag = true;
                    errOutputFileName << argv [ i ];
                    i++;
                }
            } else if ( strcmp(argv [ i ], "-qo") == 0 ) {
                if ( i + 1 < argc ) {
                    i++;
                    outputFileFlag = true;
                    outputFileName << argv [ i ];
                }
            } else if ( strcmp(argv [ i ], "-d") == 0 ) {
                debugFlag = true;
            } else if ( strcmp(argv [ i ], "-p") == 0 ) {
#ifdef __PARALLEL_MODE
                parallelFlag = true;
#else
                fprintf(stderr, "\nCan't use -p, not compiled with parallel support\a\n\n");
                exit(EXIT_FAILURE);
#endif
            } else if ( strcmp(argv [i], "-t") == 0) {
#ifdef _OPENMP
                if ( i + 1 < argc ) {
                    i++;
                    int numberOfThreads = strtol(argv [ i ], NULL, 10);
                    omp_set_num_threads (numberOfThreads);
                }
#else
                fprintf(stderr, "\nCan't use -t, not compiled with OpenMP support\a\n\n");
                exit(EXIT_FAILURE);
#endif
            } else { // Arguments not handled by OOFEM is to be passed to PETSc
                modulesArgs.push_back(argv [ i ]);
            }
        }
    } else {
        if ( rank == 0 ) {
            oofem_print_help();
        }

#ifdef __USE_MPI
        MPI_Finalize();
#endif
        exit(EXIT_SUCCESS);
    }

    // check if input file given
    if ( !inputFileFlag ) {
        if ( rank == 0 ) {
            fprintf(stderr, "\nInput file not specified\a\n\n");
        }

#ifdef __USE_MPI
        MPI_Finalize();
#endif
        exit(EXIT_FAILURE);
    }

#if defined ( __PETSC_MODULE ) || defined ( __SLEPC_MODULE )
    int modulesArgc = modulesArgs.size();
    char **modulesArgv = const_cast< char ** >(& modulesArgs [ 0 ]);
#endif

#ifdef __PETSC_MODULE
    PetscInitialize(& modulesArgc, & modulesArgv, PETSC_NULL, PETSC_NULL);
#endif

#ifdef __SLEPC_MODULE
    SlepcInitialize(& modulesArgc, & modulesArgv, PETSC_NULL, PETSC_NULL);
#endif

#ifdef _PYTHON_EXTENSION
    Py_Initialize();
    // Adding . to the system path allows us to run Python functions stored in the working directory.
    PyRun_SimpleString("import sys");
    PyRun_SimpleString("sys.path.append(\".\")");
#endif

#ifdef __PARALLEL_MODE
    if ( parallelFlag ) {
        inputFileName << "." << rank;
        outputFileName << "." << rank;
        errOutputFileName << "." << rank;
    }
#endif
    if ( outputFileFlag ) {
        oofem_logger.appendLogTo( outputFileName.str() );
    }
    if ( errOutputFileFlag ) {
        oofem_logger.appendErrorTo( errOutputFileName.str() );
    }

    // print header to redirected output
    OOFEM_LOG_FORCED(PRG_HEADER_SM);

    OOFEMTXTDataReader dr( inputFileName.str() );
    auto problem = :: InstanciateProblem(dr, _processor, contextFlag, NULL, parallelFlag);
    dr.finish();
    if ( !problem ) {
        OOFEM_LOG_ERROR("Couldn't instanciate problem, exiting");
        exit(EXIT_FAILURE);
    }

    problem->checkProblemConsistency();
    problem->init();

    if ( renumberFlag ) {
        problem->setRenumberFlag();
    }

    if ( restartFlag ) {
        try {
            FileDataStream stream(problem->giveContextFileName(restartStep, 0), false);
            problem->restoreContext(stream, CM_State | CM_Definition);
        } catch ( const FileDataStream::CantOpen & e ) {
            printf("%s", e.what());
            exit(1);
        } catch ( ContextIOERR & c ) {
            c.print();
            exit(1);
        }
        problem->initStepIncrements();
    } else if ( adaptiveRestartFlag ) {
        problem->initializeAdaptive(adaptiveRestartFlag);
        problem->saveStepContext(problem->giveCurrentStep(),CM_State);
        // exit (1);
    }

    if ( debugFlag ) {
        oofem_debug(*problem);
    }


    try {
        problem->solveYourself();
    } catch(OOFEM_Terminate & c) {
        problem = nullptr;

        oofem_finalize_modules();

        return 1;
    }

    problem->terminateAnalysis();
#ifdef __PARALLEL_MODE
    if ( parallelFlag ) {
        DynamicCommunicationBuffer :: printInfo();
    }
#endif
    oofem_logger.printStatistics();
    problem = nullptr;

    oofem_finalize_modules();

    return 0;
}

void oofem_print_help()
{
    printf("\nOptions:\n\n");
    printf("  -v  prints oofem version\n");
    printf("  -f  (string) input file name\n");
    printf("  -r  (int) restarts analysis from given step\n");
    printf("  -ar (int) restarts adaptive analysis from given step\n");
    printf("  -l  (int) sets treshold for log messages (Errors=0, Warnings=1,\n");
    printf("            Relevant=2, Info=3, Debug=4)\n");
    printf("  -rn turns on renumbering\n");
    printf("  -qo (string) redirects the standard output stream to given file\n");
    printf("  -qe (string) redirects the standard error stream to given file\n");
    printf("  -c  creates context file for each solution step\n");
    printf("\n");
    oofem_print_epilog();
}

#ifndef HOST_TYPE
 #define HOST_TYPE "unknown"
#endif

void oofem_print_version()
{
    printf("\n%s (%s, %s)\nGit RepoURL: %s\n    Branch: %s, Hash: %s\n", PRG_VERSION, HOST_TYPE, MODULE_LIST, OOFEM_GIT_REPOURL, OOFEM_GIT_BRANCH, OOFEM_GIT_HASH);
    oofem_print_epilog();
}

void oofem_print_epilog()
{
    printf("\n%s\n", OOFEM_COPYRIGHT);
    printf("This is free software; see the source for copying conditions.  There is NO\n");
    printf("warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.\n\n");
}

void oofem_finalize_modules()
{
#ifdef __PETSC_MODULE
    PetscFinalize();
#endif

#ifdef __SLEPC_MODULE
    SlepcFinalize();
#endif

#ifdef __USE_MPI
    MPI_Finalize();
#endif

#ifdef _PYTHON_EXTENSION
    Py_Finalize();
#endif
}

//#include "loadbalancer.h"
//#include "xfem/iga.h"
#include "floatmatrix.h"
#include "domain.h"
#include "element.h"
void oofem_debug(EngngModel &emodel)
{
    FloatMatrix k;
    emodel.giveDomain(1)->giveElement(1)->giveCharacteristicMatrix(k, ConductivityMatrix, NULL);
    emodel.giveDomain(1)->giveElement(1)->giveCharacteristicMatrix(k, CapacityMatrix, NULL);

    //FloatMatrix k;
    //((BsplinePlaneStressElement*)emodel.giveDomain(1)->giveElement(1))->giveCharacteristicMatrix(k, StiffnessMatrix, NULL);

#ifdef __PARALLEL_MODE
    //LoadBalancer* lb = emodel.giveDomain(1)->giveLoadBalancer();
    //lb->calculateLoadTransfer();
    //lb->migrateLoad();
    //exit(1);
#endif
}

// Empty functions just so that we can link to the library even with oofeg compilation.
#ifdef __OOFEG
void ESICustomize(Widget parent_pane) { }
oofegGraphicContext gc [ OOFEG_LAST_LAYER ];
EView *myview;
void deleteLayerGraphics(int iLayer) { }
#endif
