/* $Header: /home/cvs/bp/oofem/main/src/main.C,v 1.4.4.1 2004/04/05 15:19:41 bp Exp $ */
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
 *               Copyright (C) 1993 - 2008   Borek Patzak
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
#endif

using namespace oofem;

// debug
void oofem_debug(EngngModel *emodel);


int gPF = 0;
void oofem_print_help();
void oofem_print_version();
void oofem_print_epilog();

/* Default oofem loggers */
Logger oofem::oofem_logger(Logger :: LOG_LEVEL_INFO, stdout);
Logger oofem::oofem_errLogger(Logger :: LOG_LEVEL_WARNING, stderr);


int main(int argc, char *argv[])
{
#ifndef _MSC_VER
    std :: set_new_handler(freeStoreError);   // prevents memory overflow
#endif

    int i;
    int contextFlag = 0, inputFileFlag = 0, restartFlag = 0, adaptiveRestartFlag = 0, restartStepInfo [ 2 ];
    int renumberFlag = 0;
    bool debugFlag = false;
    char inputFileName [ MAX_FILENAME_LENGTH + 10 ], buff [ MAX_FILENAME_LENGTH ];
    EngngModel *problem = 0;

    // print prg header on stdout
    printf("%s", PRG_HEADER_SM);

#ifdef __PARALLEL_MODE
    char fileName [ MAX_FILENAME_LENGTH ];
    int rank = 0;
#ifdef __USE_MPI
    MPI_Init(& argc, & argv);
    MPI_Comm_rank(MPI_COMM_WORLD, & rank);
#endif
#endif

#ifdef __PETSC_MODULE
    PetscInitialize(& argc, & argv, PETSC_NULL, PETSC_NULL);
#endif

#ifdef __SLEPC_MODULE
    SlepcInitialize(&argc, &argv, PETSC_NULL, PETSC_NULL);
#endif

    //
    // check for options
    //
    if ( argc != 1 ) {
        for ( i = 1; i < argc; i++ ) {
            if ( ( strcmp(argv [ i ], "-context") == 0 ) || ( strcmp(argv [ i ], "-c") == 0 ) ) {
                contextFlag = 1;
            } else if ( strcmp(argv [ i ], "-v") == 0 ) {
                oofem_print_version();
                if ( argc == 2 ) {
                    exit(1);     // exit if only "-v" option
                }
            } else if ( strcmp(argv [ i ], "-f") == 0 ) {
                if ( i + 1 < argc ) {
                    strcpy(inputFileName, argv [ i + 1 ]);
                    inputFileFlag = 1;
                }
            } else if ( strcmp(argv [ i ], "-r") == 0 ) {
                if ( i + 1 < argc ) {
                    strcpy(buff, argv [ i + 1 ]);
                    restartFlag = 1;
                    restartStepInfo [ 0 ] = strtol(buff, ( char ** ) NULL, 10);
                    restartStepInfo [ 1 ] = 0;
                }
            } else if ( strcmp(argv [ i ], "-rn") == 0 ) {
                renumberFlag = 1;
            } else if ( strcmp(argv [ i ], "-ar") == 0 ) {
                if ( i + 1 < argc ) {
                    strcpy(buff, argv [ i + 1 ]);
                    adaptiveRestartFlag = strtol(buff, ( char ** ) NULL, 10);
                }
            } else if ( strcmp(argv [ i ], "-l") == 0 ) {
                if ( i + 1 < argc ) {
                    strcpy(buff, argv [ i + 1 ]);
                    int level = strtol(buff, ( char ** ) NULL, 10);
                    oofem_logger.setLogLevel(level);
                    oofem_errLogger.setLogLevel(level);
                }
            } else if ( strcmp(argv [ i ], "-qe") == 0 ) {
                if ( i + 1 < argc ) {
                    strcpy(buff, argv [ i + 1 ]);
#ifdef __PARALLEL_MODE
                    sprintf(fileName, "%s.%d", buff, rank);
                    oofem_errLogger.appendlogTo(fileName);
#endif
#ifndef __PARALLEL_MODE
                    oofem_errLogger.appendlogTo(buff);
#endif
                }
            } else if ( strcmp(argv [ i ], "-qo") == 0 ) {
                if ( i + 1 < argc ) {
                    strcpy(buff, argv [ i + 1 ]);
#ifdef __PARALLEL_MODE
                    sprintf(fileName, "%s.%d", buff, rank);
                    oofem_logger.appendlogTo(fileName);
#endif
#ifndef __PARALLEL_MODE
                    oofem_logger.appendlogTo(buff);
#endif
                    // print header to redirected output
                    LOG_FORCED_MSG(oofem_logger, PRG_HEADER_SM);
                }
            } else if ( strcmp(argv [ i ], "-d") == 0 ) {
                debugFlag = true;
            }
        }
    } else {
        oofem_print_help();
        exit(1);
    }

    // check if input file given
    if ( !inputFileFlag ) {
        /*
         * ::giveInputDataFileName(inputFileName, MAX_FILENAME_LENGTH);
         */
        fprintf(stderr, "\nInput file not specified\a\n");
        exit(1);
    }

#ifdef __PARALLEL_MODE
    strcpy(fileName, inputFileName);
    sprintf(inputFileName, "%s.%d", fileName, rank);
#endif

    OOFEMTXTDataReader dr(inputFileName);
    problem = :: InstanciateProblem(& dr, _processor, contextFlag);
    dr.finish();

    problem->checkProblemConsistency();

    if ( renumberFlag ) {
        problem->setRenumberFlag();
    }

    if ( restartFlag ) {
        try {
            problem->restoreContext(NULL, CM_State, ( void * ) restartStepInfo);
        } catch ( ContextIOERR &c ) {
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
    } catch ( OOFEM_Terminate &c ) {
        delete problem;
#ifdef __PETSC_MODULE
        PetscFinalize();
#endif

#ifdef __SLEPC_MODULE
        SlepcFinalize();
#endif

#ifdef __PARALLEL_MODE
        MPI_Finalize();
#endif

        return 1;
    }

    problem->terminateAnalysis();
#ifdef __PARALLEL_MODE
    DynamicCommunicationBuffer :: printInfo();
#endif
    oofem_errLogger.printStatistics();
    delete problem;

#ifdef __PETSC_MODULE
    PetscFinalize();
#endif

#ifdef __SLEPC_MODULE
        SlepcFinalize();
#endif

#ifdef __PARALLEL_MODE
    MPI_Finalize();
#endif

    return 0;
}

void oofem_print_help() {
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

void oofem_print_version() {
    printf("\n%s (%s, %s)\nof %s on %s\n", PRG_VERSION, HOST_TYPE, MODULE_LIST, __DATE__, HOST_NAME);
    oofem_print_epilog();
}

void
oofem_print_epilog() {
    printf("\n%s\n", OOFEM_COPYRIGHT);
    printf("This is free software; see the source for copying conditions.  There is NO\n");
    printf("warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.\n\n");
}

#ifndef __MAKEDEPEND
#include "loadbalancer.h"
#endif
void oofem_debug(EngngModel *emodel)
{
#ifdef __PARALLEL_MODE
    //LoadBalancer* lb = emodel->giveDomain(1)->giveLoadBalancer();
    //lb->calculateLoadTransfer();
    //lb->migrateLoad();
    //exit(1);
#endif
}

#endif
