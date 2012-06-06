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
 *               Copyright (C) 1993 - 2012   Borek Patzak
 *
 *
 *
 *       Czech Technical University, Faculty of Civil Engineering,
 *   Department of Structural Mechanics, 166 29 Prague, Czech Republic
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

//  MAIN
//  Solves finite element problems.
//

#ifndef __OOFEG

#include "engngm.h"
#include "freestor.h"
#include "compiler.h"

#include "oofemtxtdatareader.h"
#include "util.h"
#include "oofemdef.h"
#include "usrdefsub.h"
#include "error.h"
#include "logger.h"
#include "contextioerr.h"
#include "oofem_terminate.h"

#ifdef __PARALLEL_MODE
 #include "dyncombuff.h"
#endif

#ifdef __PETSC_MODULE
 #ifndef __MAKEDEPEND
  #include "petsc.h"
 #endif
#endif

#ifdef __SLEPC_MODULE
 #ifndef __MAKEDEPEND
  #include "slepceps.h"
 #endif
#endif

#ifndef __MAKEDEPEND
 #include <stdio.h>
 #include <string.h>
 #include <new>
 #include <sstream>
// For passing PETSc/SLEPc arguments.
 #include <fstream>
 #include <iterator>
#endif

#include "classfactory.h"

using namespace oofem;

// debug
void oofem_debug(EngngModel *emodel);

void oofem_print_help();
void oofem_print_version();
void oofem_print_epilog();

// Finalize PETSc, SLEPc and MPI
void oofem_finalize_modules();

/* Default oofem loggers */
Logger oofem :: oofem_logger(Logger :: LOG_LEVEL_INFO, stdout);
Logger oofem :: oofem_errLogger(Logger :: LOG_LEVEL_WARNING, stderr);

/* Global class factory */
ClassFactory oofem :: classFactory;

int main(int argc, char *argv[])
{
#ifndef _MSC_VER
    std :: set_new_handler(freeStoreError);   // prevents memory overflow
#endif

    int i;
    int adaptiveRestartFlag = 0, restartStepInfo [ 2 ];
    bool parallelFlag = false, renumberFlag = false, debugFlag = false, contextFlag = false, restartFlag = false,
         inputFileFlag = false, outputFileFlag = false, errOutputFileFlag = false;
    std :: stringstream inputFileName, outputFileName, errOutputFileName;
    std :: vector< const char * >modulesArgs;
    EngngModel *problem = 0;

    int rank = 0;

#ifdef __PARALLEL_MODE
 #ifdef __USE_MPI
    MPI_Init(& argc, & argv);
    MPI_Comm_rank(MPI_COMM_WORLD, & rank);
 #endif
    parallelFlag = true; ///@todo The default should be false even for parallel builds, and turned on with -p
#endif

    //
    // check for options
    //
    if ( argc != 1 ) {
        // argv[0] is not read by PETSc and SLEPc.
        modulesArgs.push_back(argv [ 0 ]);
        for ( i = 1; i < argc; i++ ) {
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
                    restartStepInfo [ 0 ] = strtol(argv [ i ] , NULL, 10);
                    restartStepInfo [ 1 ] = 0;
                }
            } else if ( strcmp(argv [ i ], "-rn") == 0 ) {
                renumberFlag = true;
            } else if ( strcmp(argv [ i ], "-ar") == 0 ) {
                if ( i + 1 < argc ) {
                    i++;
                    adaptiveRestartFlag = strtol( argv [ i ], NULL, 10);
                }
            } else if ( strcmp(argv [ i ], "-l") == 0 ) {
                if ( i + 1 < argc) {
                    i++;
                    int level = strtol( argv [ i ] , NULL, 10);
                    oofem_logger.setLogLevel(level);
                    oofem_errLogger.setLogLevel(level);
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
            } else { // Arguments not handled by OOFEM is to be passed to PETSc
                modulesArgs.push_back( argv [ i ] );
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

#ifdef __PARALLEL_MODE
    if ( parallelFlag ) {
        inputFileName << "." << rank;
        outputFileName << "." << rank;
        errOutputFileName << "." << rank;
    }
#endif
    if ( outputFileFlag ) {
        oofem_logger.appendlogTo( const_cast < char * > ( outputFileName.str().c_str() ) );
    }
    if ( errOutputFileFlag ) {
        oofem_errLogger.appendlogTo( const_cast < char * > ( errOutputFileName.str().c_str() ) );
    }

    // print header to redirected output
    LOG_FORCED_MSG(oofem_logger, PRG_HEADER_SM);

    OOFEMTXTDataReader dr( inputFileName.str().c_str() );
    problem = :: InstanciateProblem(& dr, _processor, contextFlag, NULL, parallelFlag);
    dr.finish();

    problem->checkProblemConsistency();
    problem->init();

    if ( renumberFlag ) {
        problem->setRenumberFlag();
    }

    if ( restartFlag ) {
        try {
            problem->restoreContext(NULL, CM_State, ( void * ) restartStepInfo);
        } catch(ContextIOERR & c) {
            c.print();
            exit(1);
        }
        problem->initStepIncrements();
    } else if ( adaptiveRestartFlag ) {
        problem->initializeAdaptive(adaptiveRestartFlag);
        problem->saveContext(NULL, CM_State);
        // exit (1);
    }

    if ( debugFlag ) {
        oofem_debug(problem);
    }


    try {
        problem->solveYourself();
    } catch(OOFEM_Terminate & c) {
        delete problem;

        oofem_finalize_modules();

        return 1;
    }

    problem->terminateAnalysis();
#ifdef __PARALLEL_MODE
    DynamicCommunicationBuffer :: printInfo();
#endif
    oofem_errLogger.printStatistics();
    delete problem;

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
    printf("\n%s (%s, %s)\nof %s on %s\n", PRG_VERSION, HOST_TYPE, MODULE_LIST, __DATE__, HOST_NAME);
    oofem_print_epilog();
}

void oofem_print_epilog()
{
    printf("\n%s\n", OOFEM_COPYRIGHT);
    printf("This is free software; see the source for copying conditions.  There is NO\n");
    printf("warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.\n\n");
}

void oofem_finalize_modules() {
#ifdef __PETSC_MODULE
    PetscFinalize();
#endif

#ifdef __SLEPC_MODULE
    SlepcFinalize();
#endif

#ifdef __USE_MPI
    MPI_Finalize();
#endif
}

#ifndef __MAKEDEPEND
 #include "loadbalancer.h"
 #include "iga.h"
#endif
void oofem_debug(EngngModel *emodel)
{
    //FloatMatrix k;
    //((BsplinePlaneStressElement*)emodel->giveDomain(1)->giveElement(1))->giveCharacteristicMatrix(k, StiffnessMatrix, NULL);

#ifdef __PARALLEL_MODE
    //LoadBalancer* lb = emodel->giveDomain(1)->giveLoadBalancer();
    //lb->calculateLoadTransfer();
    //lb->migrateLoad();
    //exit(1);
#endif
}

#endif
