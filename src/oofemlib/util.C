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

#include "util.h"
#include "engngm.h"
#include "classfactory.h"
#include "inputrecord.h"
#include "datareader.h"
#include "error.h"

#include <cstring>
#if defined ( __GNUC__ ) && defined ( HAVE_EXECINFO_H )
#include <cxxabi.h>
#include <execinfo.h>
#include <cstdio>
#include <cstdlib>
#endif

namespace oofem {

void print_stacktrace(FILE *out, int skip, unsigned int max_frames)
{
#if defined ( __GNUC__ ) && defined ( HAVE_EXECINFO_H )
#if 0
    // Much simpler example from SO
    void *trace_elems[max_frames];
    int trace_elem_count(backtrace( trace_elems, max_frames ));
    char **stack_syms(backtrace_symbols( trace_elems, trace_elem_count ));
    for ( int i = 0 ; i < trace_elem_count ; ++i ) {
        std::cout << stack_syms[i] << "\n";
    }
    free( stack_syms );
#endif
    // Taken from https://idlebox.net/2008/0901-stacktrace-demangled/ which indicated free usage.
    int addrlen = 0;
    fprintf(out, "stack trace:\n");

    // storage array for stack trace address data
    void *addrlist [ max_frames + 1 ];

    // retrieve current stack addresses
    addrlen = backtrace( addrlist, sizeof( addrlist ) / sizeof( void * ) );
    if ( addrlen == 0 ) {
        fprintf(out, "  <empty, possibly corrupt>\n");
        return;
    }

    // resolve addresses into strings containing "filename(function+address)",
    // this array must be free()-ed
    char **symbollist;
    symbollist = backtrace_symbols(addrlist, addrlen);
    // allocate string which will be filled with the demangled function name
    size_t funcnamesize = 256;
    char *funcname = ( char * ) malloc(funcnamesize);

    // iterate over the returned symbol lines. skip the first, it is the
    // address of this function.
    for ( int i = skip+1; i < addrlen; i++ ) {
        char *begin_name = 0, *begin_offset = 0, *end_offset = 0;

        // find parentheses and +address offset surrounding the mangled name:
        // ./module(function+0x15c) [0x8048a6d]
        for ( char *p = symbollist [ i ]; * p; ++p ) {
            if ( * p == '(' ) {
                begin_name = p;
            } else if ( * p == '+' ) {
                begin_offset = p;
            } else if ( * p == ')' && begin_offset ) {
                end_offset = p;
                break;
            }
        }

        if ( begin_name && begin_offset && end_offset &&
            begin_name < begin_offset ) {
            * begin_name++ = '\0';
        * begin_offset++ = '\0';
        * end_offset = '\0';

        // mangled name is now in [begin_name, begin_offset) and caller
        // offset in [begin_offset, end_offset). now apply
        // __cxa_demangle():

        int status;
        char *ret = abi :: __cxa_demangle(begin_name,
                                            funcname, & funcnamesize, & status);
        if ( status == 0 ) {
            funcname = ret; // use possibly realloc()-ed string
            fprintf(out, "  %s : %s+%s\n",
                    symbollist [ i ], funcname, begin_offset);
        } else {
            // demangling failed. Output function name as a C function with
            // no arguments.
            fprintf(out, "  %s : %s()+%s\n",
                    symbollist [ i ], begin_name, begin_offset);
        }
            } else {
                // couldn't parse the line? print the whole line.
                fprintf(out, "  %s\n", symbollist [ i ]);
            }
    }

    free(funcname);
    free(symbollist);
#else
    fprintf(out, "No backtrace available\n");
#endif
}


std::unique_ptr<EngngModel> InstanciateProblem(DataReader &dr, problemMode mode, int contextFlag, EngngModel *_master, bool parallelFlag)
{
    std :: string problemName, dataOutputFileName, desc;

    dataOutputFileName = dr.giveOutputFileName();
    desc = dr.giveDescription();

    /* here we need copy of input record. The pointer returned by dr.giveInputRecord can (and will)
     * be updated as reading e-model components (nodes, etc). But we need this record being available
     * through the whole e-model instanciation
     */
    auto emodelir = dr.giveInputRecord(DataReader :: IR_emodelRec, 1).clone();
    emodelir->giveRecordKeywordField(problemName); ///@todo Make this function robust, it can't be allowed to fail (the record keyword is not a normal field-id)

    auto problem = classFactory.createEngngModel(problemName.c_str(), 1, _master);
    if ( !problem ) {
        OOFEM_WARNING( "Failed to construct engineering model of type \"%s\".\n", problemName.c_str() );
        return NULL;
    }

    problem->setProblemMode(mode);
    problem->setParallelMode(parallelFlag);

    if ( contextFlag ) {
        problem->setContextOutputMode(COM_Always);
    }

    problem->instanciateYourself( dr, *emodelir, dataOutputFileName.c_str(), desc.c_str() );
    emodelir->finish();

    return problem;
}
} // end namespace oofem
