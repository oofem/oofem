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

#include "../sm/Elements/PlaneStress/qplanstrssphf.h"
#include "classfactory.h"

namespace oofem {
REGISTER_Element(QPlaneStressPhF2d);

QPlaneStressPhF2d::QPlaneStressPhF2d( int n, Domain *aDomain ) : PhaseFieldElement( n, aDomain ), 
QPlaneStress2d(n, aDomain ) { }

void
QPlaneStressPhF2d :: giveDofManDofIDMask(int inode, IntArray &answer) const
{
    answer = { D_u, D_v, T_f }; ///@todo add damage dofID later
}

void
QPlaneStressPhF2d::giveDofManDofIDMask_u( IntArray &answer )
{
    answer = { D_u, D_v };
}

void
QPlaneStressPhF2d::giveDofManDofIDMask_d( IntArray &answer )
{
    answer = IntArray{ T_f };
}



} // end namespace oofem
