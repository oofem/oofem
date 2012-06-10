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

#ifndef parallel_h
#define parallel_h

// turn on parallel verbose mode if uncommented
//#define __VERBOSE_PARALLEL

// turn on using MPI if uncommented
#define __USE_MPI

//
// ===========>  Do not modify bellow this line <==============
//

// macro for verbose print in parallel mode
//#define VERBOSEPARALLEL_PRINT(service, str,rank) fprintf(stderr, "\n[process rank %3d] %s - %s",rank, service, str);
#define VERBOSEPARALLEL_PRINT(service, str, rank) fprintf(stderr, "\n[%d] %s - %s", rank, service, str);


#ifdef __USE_MPI
// if MPI used, include headers
 #ifndef __MAKEDEPEND
  #include <mpi.h>
 #endif
#endif

namespace oofem {
#define PROCESSOR_NAME_LENGTH 60
} // end namespace oofem
#endif // parallel_h
