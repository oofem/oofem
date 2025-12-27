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
 *               Copyright (C) 1993 - 2025   Borek Patzak
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

#include "floatmatrix.h"
#include "floatarray.h"
#include "mathfem.h"
#include "index.h"
#pragma once

using namespace FTensor;


/**
 * Functions providing basic operations with tensors
 * based on FTensor library
 * It provides transformation from the second and fourth order tensors in 3d to FloatArrayF<9> and FloatArray<9,9>
 * It provides calculation of determinant, cofactor, inverse, and tensor cross product
 * More to be added gradually
 * @author Martin Horak
 **/

namespace oofem {
class Tensor4_3d : public Tensor4< double, 3, 3, 3, 3 >
{
public:
    using Tensor4< double, 3, 3, 3, 3 >::Tensor4;

    /**
     * Creates fourth-order tensor in 3d
     * Initialized with zeros
     */
    Tensor4_3d() {
        // the first column
        this->data [ 0 ] [ 0 ] [ 0 ] [ 0 ] = 0;
        this->data [ 0 ] [ 1 ] [ 0 ] [ 0 ] = 0.;
        this->data [ 0 ] [ 2 ] [ 0 ] [ 0 ] = 0.;
        this->data [ 1 ] [ 0 ] [ 0 ] [ 0 ] = 0.;
        this->data [ 1 ] [ 1 ] [ 0 ] [ 0 ] = 0.;
        this->data [ 1 ] [ 2 ] [ 0 ] [ 0 ] = 0.;
        this->data [ 2 ] [ 0 ] [ 0 ] [ 0 ] = 0.;
        this->data [ 2 ] [ 1 ] [ 0 ] [ 0 ] = 0.;
        this->data [ 2 ] [ 2 ] [ 0 ] [ 0 ] = 0.;
        // the second column
        this->data [ 0 ] [ 0 ] [ 1 ] [ 1 ] = 0.;
        this->data [ 0 ] [ 1 ] [ 1 ] [ 1 ] = 0.;
        this->data [ 0 ] [ 2 ] [ 1 ] [ 1 ] = 0.;
        this->data [ 1 ] [ 0 ] [ 1 ] [ 1 ] = 0.;
        this->data [ 1 ] [ 1 ] [ 1 ] [ 1 ] = 0.;
        this->data [ 1 ] [ 2 ] [ 1 ] [ 1 ] = 0.;
        this->data [ 2 ] [ 0 ] [ 1 ] [ 1 ] = 0.;
        this->data [ 2 ] [ 1 ] [ 1 ] [ 1 ] = 0.;
        this->data [ 2 ] [ 2 ] [ 1 ] [ 1 ] = 0.;
        // the third column
        this->data [ 0 ] [ 0 ] [ 2 ] [ 2 ] = 0.;
        this->data [ 0 ] [ 1 ] [ 2 ] [ 2 ] = 0.;
        this->data [ 0 ] [ 2 ] [ 2 ] [ 2 ] = 0.;
        this->data [ 1 ] [ 0 ] [ 2 ] [ 2 ] = 0.;
        this->data [ 1 ] [ 1 ] [ 2 ] [ 2 ] = 0.;
        this->data [ 1 ] [ 2 ] [ 2 ] [ 2 ] = 0.;
        this->data [ 2 ] [ 0 ] [ 2 ] [ 2 ] = 0.;
        this->data [ 2 ] [ 1 ] [ 2 ] [ 2 ] = 0.;
        this->data [ 2 ] [ 2 ] [ 2 ] [ 2 ] = 0.;
        // the fourth column
        this->data [ 0 ] [ 0 ] [ 1 ] [ 2 ] = 0.;
        this->data [ 0 ] [ 1 ] [ 1 ] [ 2 ] = 0.;
        this->data [ 0 ] [ 2 ] [ 1 ] [ 2 ] = 0.;
        this->data [ 1 ] [ 0 ] [ 1 ] [ 2 ] = 0.;
        this->data [ 1 ] [ 1 ] [ 1 ] [ 2 ] = 0.;
        this->data [ 1 ] [ 2 ] [ 1 ] [ 2 ] = 0.;
        this->data [ 2 ] [ 0 ] [ 1 ] [ 2 ] = 0.;
        this->data [ 2 ] [ 1 ] [ 1 ] [ 2 ] = 0.;
        this->data [ 2 ] [ 2 ] [ 1 ] [ 2 ] = 0.;
        //the fifth column
        this->data [ 0 ] [ 0 ] [ 0 ] [ 2 ] = 0.;
        this->data [ 0 ] [ 1 ] [ 0 ] [ 2 ] = 0.;
        this->data [ 0 ] [ 2 ] [ 0 ] [ 2 ] = 0.;
        this->data [ 1 ] [ 0 ] [ 0 ] [ 2 ] = 0.;
        this->data [ 1 ] [ 1 ] [ 0 ] [ 2 ] = 0.;
        this->data [ 1 ] [ 2 ] [ 0 ] [ 2 ] = 0.;
        this->data [ 2 ] [ 0 ] [ 0 ] [ 2 ] = 0.;
        this->data [ 2 ] [ 1 ] [ 0 ] [ 2 ] = 0.;
        this->data [ 2 ] [ 2 ] [ 0 ] [ 2 ] = 0.;
        //the sixth column
        this->data [ 0 ] [ 0 ] [ 0 ] [ 1 ] = 0.;
        this->data [ 0 ] [ 1 ] [ 0 ] [ 1 ] = 0.;
        this->data [ 0 ] [ 2 ] [ 0 ] [ 1 ] = 0.;
        this->data [ 1 ] [ 0 ] [ 0 ] [ 1 ] = 0.;
        this->data [ 1 ] [ 1 ] [ 0 ] [ 1 ] = 0.;
        this->data [ 1 ] [ 2 ] [ 0 ] [ 1 ] = 0.;
        this->data [ 2 ] [ 0 ] [ 0 ] [ 1 ] = 0.;
        this->data [ 2 ] [ 1 ] [ 0 ] [ 1 ] = 0.;
        this->data [ 2 ] [ 2 ] [ 0 ] [ 1 ] = 0.;
        //the seventh column
        this->data [ 0 ] [ 0 ] [ 2 ] [ 1 ] = 0.;
        this->data [ 0 ] [ 1 ] [ 2 ] [ 1 ] = 0.;
        this->data [ 0 ] [ 2 ] [ 2 ] [ 1 ] = 0.;
        this->data [ 1 ] [ 0 ] [ 2 ] [ 1 ] = 0.;
        this->data [ 1 ] [ 1 ] [ 2 ] [ 1 ] = 0.;
        this->data [ 1 ] [ 2 ] [ 2 ] [ 1 ] = 0.;
        this->data [ 2 ] [ 0 ] [ 2 ] [ 1 ] = 0.;
        this->data [ 2 ] [ 1 ] [ 2 ] [ 1 ] = 0.;
        this->data [ 2 ] [ 2 ] [ 2 ] [ 1 ] = 0.;
        //the eight column
        this->data [ 0 ] [ 0 ] [ 2 ] [ 0 ] = 0.;
        this->data [ 0 ] [ 1 ] [ 2 ] [ 0 ] = 0.;
        this->data [ 0 ] [ 2 ] [ 2 ] [ 0 ] = 0.;
        this->data [ 1 ] [ 0 ] [ 2 ] [ 0 ] = 0.;
        this->data [ 1 ] [ 1 ] [ 2 ] [ 0 ] = 0.;
        this->data [ 1 ] [ 2 ] [ 2 ] [ 0 ] = 0.;
        this->data [ 2 ] [ 0 ] [ 2 ] [ 0 ] = 0.;
        this->data [ 2 ] [ 1 ] [ 2 ] [ 0 ] = 0.;
        this->data [ 2 ] [ 2 ] [ 2 ] [ 0 ] = 0.;
        //the nineght column
        this->data [ 0 ] [ 0 ] [ 1 ] [ 0 ] = 0.;
        this->data [ 0 ] [ 1 ] [ 1 ] [ 0 ] = 0.;
        this->data [ 0 ] [ 2 ] [ 1 ] [ 0 ] = 0.;
        this->data [ 1 ] [ 0 ] [ 1 ] [ 0 ] = 0.;
        this->data [ 1 ] [ 1 ] [ 1 ] [ 0 ] = 0.;
        this->data [ 1 ] [ 2 ] [ 1 ] [ 0 ] = 0.;
        this->data [ 2 ] [ 0 ] [ 1 ] [ 0 ] = 0.;
        this->data [ 2 ] [ 1 ] [ 1 ] [ 0 ] = 0.;
        this->data [ 2 ] [ 2 ] [ 1 ] [ 0 ] = 0.;
    }

    /**
     * Creates fourth-order tensor in 3d from floatmatrixf<9,9>
     */
    Tensor4_3d(const oofem::FloatMatrixF< 9, 9 > &mat) {
        // the first column
        this->data [ 0 ] [ 0 ] [ 0 ] [ 0 ] = mat.at(1, 1);
        this->data [ 0 ] [ 1 ] [ 0 ] [ 0 ] = mat.at(6, 1);
        this->data [ 0 ] [ 2 ] [ 0 ] [ 0 ] = mat.at(5, 1);
        this->data [ 1 ] [ 0 ] [ 0 ] [ 0 ] = mat.at(9, 1);
        this->data [ 1 ] [ 1 ] [ 0 ] [ 0 ] = mat.at(2, 1);
        this->data [ 1 ] [ 2 ] [ 0 ] [ 0 ] = mat.at(4, 1);
        this->data [ 2 ] [ 0 ] [ 0 ] [ 0 ] = mat.at(8, 1);
        this->data [ 2 ] [ 1 ] [ 0 ] [ 0 ] = mat.at(7, 1);
        this->data [ 2 ] [ 2 ] [ 0 ] [ 0 ] = mat.at(3, 1);
        // the second column
        this->data [ 0 ] [ 0 ] [ 1 ] [ 1 ] = mat.at(1, 2);
        this->data [ 0 ] [ 1 ] [ 1 ] [ 1 ] = mat.at(6, 2);
        this->data [ 0 ] [ 2 ] [ 1 ] [ 1 ] = mat.at(5, 2);
        this->data [ 1 ] [ 0 ] [ 1 ] [ 1 ] = mat.at(9, 2);
        this->data [ 1 ] [ 1 ] [ 1 ] [ 1 ] = mat.at(2, 2);
        this->data [ 1 ] [ 2 ] [ 1 ] [ 1 ] = mat.at(4, 2);
        this->data [ 2 ] [ 0 ] [ 1 ] [ 1 ] = mat.at(8, 2);
        this->data [ 2 ] [ 1 ] [ 1 ] [ 1 ] = mat.at(7, 2);
        this->data [ 2 ] [ 2 ] [ 1 ] [ 1 ] = mat.at(3, 2);
        // the third column
        this->data [ 0 ] [ 0 ] [ 2 ] [ 2 ] = mat.at(1, 3);
        this->data [ 0 ] [ 1 ] [ 2 ] [ 2 ] = mat.at(6, 3);
        this->data [ 0 ] [ 2 ] [ 2 ] [ 2 ] = mat.at(5, 3);
        this->data [ 1 ] [ 0 ] [ 2 ] [ 2 ] = mat.at(9, 3);
        this->data [ 1 ] [ 1 ] [ 2 ] [ 2 ] = mat.at(2, 3);
        this->data [ 1 ] [ 2 ] [ 2 ] [ 2 ] = mat.at(4, 3);
        this->data [ 2 ] [ 0 ] [ 2 ] [ 2 ] = mat.at(8, 3);
        this->data [ 2 ] [ 1 ] [ 2 ] [ 2 ] = mat.at(7, 3);
        this->data [ 2 ] [ 2 ] [ 2 ] [ 2 ] = mat.at(3, 3);
        // the fourth column
        this->data [ 0 ] [ 0 ] [ 1 ] [ 2 ] = mat.at(1, 4);
        this->data [ 0 ] [ 1 ] [ 1 ] [ 2 ] = mat.at(6, 4);
        this->data [ 0 ] [ 2 ] [ 1 ] [ 2 ] = mat.at(5, 4);
        this->data [ 1 ] [ 0 ] [ 1 ] [ 2 ] = mat.at(9, 4);
        this->data [ 1 ] [ 1 ] [ 1 ] [ 2 ] = mat.at(2, 4);
        this->data [ 1 ] [ 2 ] [ 1 ] [ 2 ] = mat.at(4, 4);
        this->data [ 2 ] [ 0 ] [ 1 ] [ 2 ] = mat.at(8, 4);
        this->data [ 2 ] [ 1 ] [ 1 ] [ 2 ] = mat.at(7, 4);
        this->data [ 2 ] [ 2 ] [ 1 ] [ 2 ] = mat.at(3, 4);
        //the fifth column
        this->data [ 0 ] [ 0 ] [ 0 ] [ 2 ] = mat.at(1, 5);
        this->data [ 0 ] [ 1 ] [ 0 ] [ 2 ] = mat.at(6, 5);
        this->data [ 0 ] [ 2 ] [ 0 ] [ 2 ] = mat.at(5, 5);
        this->data [ 1 ] [ 0 ] [ 0 ] [ 2 ] = mat.at(9, 5);
        this->data [ 1 ] [ 1 ] [ 0 ] [ 2 ] = mat.at(2, 5);
        this->data [ 1 ] [ 2 ] [ 0 ] [ 2 ] = mat.at(4, 5);
        this->data [ 2 ] [ 0 ] [ 0 ] [ 2 ] = mat.at(8, 5);
        this->data [ 2 ] [ 1 ] [ 0 ] [ 2 ] = mat.at(7, 5);
        this->data [ 2 ] [ 2 ] [ 0 ] [ 2 ] = mat.at(3, 5);
        //the sixth column
        this->data [ 0 ] [ 0 ] [ 0 ] [ 1 ] = mat.at(1, 6);
        this->data [ 0 ] [ 1 ] [ 0 ] [ 1 ] = mat.at(6, 6);
        this->data [ 0 ] [ 2 ] [ 0 ] [ 1 ] = mat.at(5, 6);
        this->data [ 1 ] [ 0 ] [ 0 ] [ 1 ] = mat.at(9, 6);
        this->data [ 1 ] [ 1 ] [ 0 ] [ 1 ] = mat.at(2, 6);
        this->data [ 1 ] [ 2 ] [ 0 ] [ 1 ] = mat.at(4, 6);
        this->data [ 2 ] [ 0 ] [ 0 ] [ 1 ] = mat.at(8, 6);
        this->data [ 2 ] [ 1 ] [ 0 ] [ 1 ] = mat.at(7, 6);
        this->data [ 2 ] [ 2 ] [ 0 ] [ 1 ] = mat.at(3, 6);
        //the seventh column
        this->data [ 0 ] [ 0 ] [ 2 ] [ 1 ] = mat.at(1, 7);
        this->data [ 0 ] [ 1 ] [ 2 ] [ 1 ] = mat.at(6, 7);
        this->data [ 0 ] [ 2 ] [ 2 ] [ 1 ] = mat.at(5, 7);
        this->data [ 1 ] [ 0 ] [ 2 ] [ 1 ] = mat.at(9, 7);
        this->data [ 1 ] [ 1 ] [ 2 ] [ 1 ] = mat.at(2, 7);
        this->data [ 1 ] [ 2 ] [ 2 ] [ 1 ] = mat.at(4, 7);
        this->data [ 2 ] [ 0 ] [ 2 ] [ 1 ] = mat.at(8, 7);
        this->data [ 2 ] [ 1 ] [ 2 ] [ 1 ] = mat.at(7, 7);
        this->data [ 2 ] [ 2 ] [ 2 ] [ 1 ] = mat.at(3, 7);
        //the eight column
        this->data [ 0 ] [ 0 ] [ 2 ] [ 0 ] = mat.at(1, 8);
        this->data [ 0 ] [ 1 ] [ 2 ] [ 0 ] = mat.at(6, 8);
        this->data [ 0 ] [ 2 ] [ 2 ] [ 0 ] = mat.at(5, 8);
        this->data [ 1 ] [ 0 ] [ 2 ] [ 0 ] = mat.at(9, 8);
        this->data [ 1 ] [ 1 ] [ 2 ] [ 0 ] = mat.at(2, 8);
        this->data [ 1 ] [ 2 ] [ 2 ] [ 0 ] = mat.at(4, 8);
        this->data [ 2 ] [ 0 ] [ 2 ] [ 0 ] = mat.at(8, 8);
        this->data [ 2 ] [ 1 ] [ 2 ] [ 0 ] = mat.at(7, 8);
        this->data [ 2 ] [ 2 ] [ 2 ] [ 0 ] = mat.at(3, 8);
        //the nineght column
        this->data [ 0 ] [ 0 ] [ 1 ] [ 0 ] = mat.at(1, 9);
        this->data [ 0 ] [ 1 ] [ 1 ] [ 0 ] = mat.at(6, 9);
        this->data [ 0 ] [ 2 ] [ 1 ] [ 0 ] = mat.at(5, 9);
        this->data [ 1 ] [ 0 ] [ 1 ] [ 0 ] = mat.at(9, 9);
        this->data [ 1 ] [ 1 ] [ 1 ] [ 0 ] = mat.at(2, 9);
        this->data [ 1 ] [ 2 ] [ 1 ] [ 0 ] = mat.at(4, 9);
        this->data [ 2 ] [ 0 ] [ 1 ] [ 0 ] = mat.at(8, 9);
        this->data [ 2 ] [ 1 ] [ 1 ] [ 0 ] = mat.at(7, 9);
        this->data [ 2 ] [ 2 ] [ 1 ] [ 0 ] = mat.at(3, 9);
    }

    /**
     * Transforms a fourth-order tensor into a floatmatrixf<9,9>, using the Voigt notation
     */

    inline FloatMatrixF< 9, 9 >to_voigt_form()
    {
        return { this->operator()(0, 0, 0, 0), this->operator()(1, 1, 0, 0), this->operator()(2, 2, 0, 0), this->operator()(1, 2, 0, 0), this->operator()(0, 2, 0, 0), this->operator()(0, 1, 0, 0), this->operator()(2, 1, 0, 0), this->operator()(2, 0, 0, 0), this->operator()(1, 0, 0, 0), this->operator()(0, 0, 1, 1), this->operator()(1, 1, 1, 1), this->operator()(2, 2, 1, 1), this->operator()(1, 2, 1, 1), this->operator()(0, 2, 1, 1), this->operator()(0, 1, 1, 1), this->operator()(2, 1, 1, 1), this->operator()(2, 0, 1, 1), this->operator()(1, 0, 1, 1), this->operator()(0, 0, 2, 2), this->operator()(1, 1, 2, 2), this->operator()(2, 2, 2, 2), this->operator()(1, 2, 2, 2), this->operator()(0, 2, 2, 2), this->operator()(0, 1, 2, 2), this->operator()(2, 1, 2, 2), this->operator()(2, 0, 2, 2), this->operator()(1, 0, 2, 2), this->operator()(0, 0, 1, 2), this->operator()(1, 1, 1, 2), this->operator()(2, 2, 1, 2), this->operator()(1, 2, 1, 2), this->operator()(0, 2, 1, 2), this->operator()(0, 1, 1, 2), this->operator()(2, 1, 1, 2), this->operator()(2, 0, 1, 2), this->operator()(1, 0, 1, 2), this->operator()(0, 0, 0, 2), this->operator()(1, 1, 0, 2), this->operator()(2, 2, 0, 2), this->operator()(1, 2, 0, 2), this->operator()(0, 2, 0, 2), this->operator()(0, 1, 0, 2), this->operator()(2, 1, 0, 2), this->operator()(2, 0, 0, 2), this->operator()(1, 0, 0, 2), this->operator()(0, 0, 0, 1), this->operator()(1, 1, 0, 1), this->operator()(2, 2, 0, 1), this->operator()(1, 2, 0, 1), this->operator()(0, 2, 0, 1), this->operator()(0, 1, 0, 1), this->operator()(2, 1, 0, 1), this->operator()(2, 0, 0, 1), this->operator()(1, 0, 0, 1), this->operator()(0, 0, 2, 1), this->operator()(1, 1, 2, 1), this->operator()(2, 2, 2, 1), this->operator()(1, 2, 2, 1), this->operator()(0, 2, 2, 1), this->operator()(0, 1, 2, 1), this->operator()(2, 1, 2, 1), this->operator()(2, 0, 2, 1), this->operator()(1, 0, 2, 1), this->operator()(0, 0, 2, 0), this->operator()(1, 1, 2, 0), this->operator()(2, 2, 2, 0), this->operator()(1, 2, 2, 0), this->operator()(0, 2, 2, 0), this->operator()(0, 1, 2, 0), this->operator()(2, 1, 2, 0), this->operator()(2, 0, 2, 0), this->operator()(1, 0, 2, 0), this->operator()(0, 0, 1, 0), this->operator()(1, 1, 1, 0), this->operator()(2, 2, 1, 0), this->operator()(1, 2, 1, 0), this->operator()(0, 2, 1, 0), this->operator()(0, 1, 1, 0), this->operator()(2, 1, 1, 0), this->operator()(2, 0, 1, 0), this->operator()(1, 0, 1, 0) };
    }

    //@todo, make this working
    /* attempt to implement noinline print, but it doesn't work, why?
     * void __attribute__((noinline)) printYourself() const
     * // Prints the receiver on screen.
     * {
     * printf("4th order tensor");
     * for ( int i = 1; i <= 3; i++ ) {
     *  for ( int j = 1; j <= 3; j++ ) {
     *    for ( int k = 1; k <= 3; k++ ) {
     *      for ( int l = 1; l <= 3; l++ ) {
     *        printf( "Component %d %d %d %d is: %10.3e \n", i,j,k,l,this->data[i][j][k][l] );
     *      }
     *    }
     *  }
     * }
     * }
     *
     * void __attribute__((noinline)) printComponent(int i, int j, int k, int l) const
     * // Prints the receiver on screen.
     * {
     * printf( "Component %d %d %d %d is: %10.3e\n  ", i,j,k,l,this->data[i][j][k][l] );
     * }
     */
};
} // end namespace oofem
