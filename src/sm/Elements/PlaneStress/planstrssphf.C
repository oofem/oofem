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

#include "../sm/Elements/PlaneStress/planstrssphf.h"
#include "classfactory.h"

namespace oofem {
REGISTER_Element(PlaneStressPhF2d);

PlaneStressPhF2d::PlaneStressPhF2d( int n, Domain *aDomain ) : PhaseFieldElement( n, aDomain ), 
PlaneStress2d(n, aDomain ) { }

void
PlaneStressPhF2d :: giveDofManDofIDMask(int inode, EquationID, IntArray &answer) const
{
    answer = {D_u, D_v, T_f}; ///@todo add damage dofID later
}

void
PlaneStressPhF2d :: giveDofManDofIDMask_u( IntArray &answer )
{
    answer = {D_u, D_v};
}

void
PlaneStressPhF2d :: giveDofManDofIDMask_d( IntArray &answer )
{
    answer = {T_f};
}



} // end namespace oofem
