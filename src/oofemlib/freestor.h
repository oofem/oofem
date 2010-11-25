/* $Header: /home/cvs/bp/oofem/oofemlib/src/freestor.h,v 1.5 2003/04/06 14:08:24 bp Exp $ */
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


/*
 * //   ********************************************
 * //   *** DYNAMIC MEMORY ALLOCATION PROCEDURES ***
 * //   ********************************************
 */

#ifndef freestor_h
#define freestor_h

namespace oofem {
/*
 * This file does not define a class. Rather, it provides a few procedures
 * related with dynamic memory allocation.
 */

/**
 * Allocates space for given number of doubles and returns pointer to newly allocated space in memory.
 * Uses standard malloc. If not enough memory is available, freeStoreError procedure is called.
 * @return pointer to newly allocated space.
 */
double *allocDouble(int);
/**
 * Allocates space for given number of ints and returns pointer to newly allocated space in memory.
 * Uses standard malloc. If not enough memory is available, freeStoreError procedure is called.
 * @return pointer to newly allocated space.
 */
int *allocInt(int);
/**
 * This function is called whenever system is unable to allocate
 * required memory.
 */
void     freeStoreError();
/**
 * Deallocates the array of decimals 'a'.
 */
void     freeInt(int *);
/**
 * Deallocates the array of doubles 'a'.
 */
void     freeDouble(double *);
} // end namespace oofem
#endif // freestor_h







