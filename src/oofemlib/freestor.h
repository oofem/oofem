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
 *               Copyright (C) 1993 - 2011   Borek Patzak
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

#ifndef freestor_h
#define freestor_h

namespace oofem {
/*
 * This file does not define a class. Rather, it provides a few procedures
 * related with dynamic memory allocation.
 */

/**
 * Allocates space for given number of doubles and returns pointer to newly allocated space in memory.
 * Uses standard malloc. If not enough memory is available, error is printed then freeStoreError is called.
 * @param n Number of doubles to allocate.
 * @return Pointer to newly allocated space.
 */
double *allocDouble(int n);
/**
 * Allocates space for given number of ints and returns pointer to newly allocated space in memory.
 * Uses standard malloc. If not enough memory is available, error is printed then freeStoreError is called.
 * @param n Number of integers to allocate.
 * @return Pointer to newly allocated space.
 */
int *allocInt(int n);
/**
 * This function is called whenever system is unable to allocate
 * required memory.
 */
void freeStoreError();
/**
 * Deallocates the array of integers.
 * @param a Integers to deallocate.
 */
void freeInt(int *a);
/**
 * Deallocates the array of doubles.
 * @param a Doubles to deallocate.
 */
void freeDouble(double *a);
} // end namespace oofem
#endif // freestor_h







