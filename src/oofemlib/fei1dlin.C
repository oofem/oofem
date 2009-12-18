/* $Header: /home/cvs/bp/oofem/oofemlib/src/fei1dlin.C,v 1.1 2003/04/06 14:08:24 bp Exp $ */
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

#include "fei1dlin.h"
#include "mathfem.h"
#include "flotmtrx.h"
#include "dofmanager.h"
#include "node.h"

namespace oofem {

void
FEI1dLin :: evalN(FloatArray &answer, const FloatArray &lcoords, double time)
{
    double ksi = lcoords.at(1);
    answer.resize(2);

    answer.at(1) = ( 1. - ksi ) * 0.5;
    answer.at(2) = ( 1. + ksi ) * 0.5;

    return;
}

void
FEI1dLin :: evaldNdx(FloatMatrix &answer, const FloatArray **nc, const FloatArray &lcoords, double time)
{
    double l = this->computeLength(nc);
    answer.resize(2, 1);

    answer.at(1, 1) = -1.0 / l;
    answer.at(2, 1) =  1.0 / l;
}

void
FEI1dLin :: local2global(FloatArray &answer, const FloatArray **nc, const FloatArray &lcoords, double time)
{
    FloatArray n(2);
    answer.resize(1);

    this->evalN(n, lcoords, time);
    answer.at(1) = ( n.at(1) * nc [ 0 ]->at(cindx) +
                    n.at(2) * nc [ 1 ]->at(cindx) );
}

int
FEI1dLin :: global2local(FloatArray &answer, const FloatArray **nc, const FloatArray &coords, double time)
{
    double ksi, x1, x2;
    answer.resize(1);


    x1 = nc [ 0 ]->at(cindx);
    x2 = nc [ 1 ]->at(cindx);

    answer.at(1) = ksi = ( 2.0 * coords.at(1) - ( x1 + x2 ) ) / ( x2 - x1 );
    return ( fabs(ksi) <= 1.0 ) ? 1 : 0;
}

double
FEI1dLin :: giveTransformationJacobian(const FloatArray **nc, const FloatArray &lcoords, double time)
{
    double l = computeLength(nc);
    return 0.5 * l;
}

double
FEI1dLin :: computeLength(const FloatArray **nc)
{
    return ( nc [ 1 ]->at(cindx) - nc [ 0 ]->at(cindx) );
}

} // end namespace oofem
