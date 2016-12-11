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

#include "unstructuredgridfield.h"
#include "feinterpol.h"
namespace oofem {

  FEI2dLineLin  UnstructuredGridField::Cell::i1 (1,2);
  FEI2dLineQuad UnstructuredGridField::Cell::i2 (1,2);
  FEI2dTrLin    UnstructuredGridField::Cell::i3 (1,2);
  FEI2dTrQuad   UnstructuredGridField::Cell::i4 (1,2);
  FEI2dQuadLin  UnstructuredGridField::Cell::i5 (1,2);
  FEI2dQuadQuad UnstructuredGridField::Cell::i6 (1,2);
  FEI3dTetLin   UnstructuredGridField::Cell::i7;

  FEInterpolation* UnstructuredGridField::Cell::interpTable[] = {
    &UnstructuredGridField::Cell::i1,
    &UnstructuredGridField::Cell::i2,
    &UnstructuredGridField::Cell::i3,
    &UnstructuredGridField::Cell::i4,
    &UnstructuredGridField::Cell::i5,
    &UnstructuredGridField::Cell::i6,
    &UnstructuredGridField::Cell::i7
  };
} // end namespace oofem
