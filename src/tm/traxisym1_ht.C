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

#include "traxisym1_ht.h"
#include "node.h"
#include "material.h"
#include "crosssection.h"
#include "gausspnt.h"
#include "gaussintegrationrule.h"
#include "flotmtrx.h"
#include "flotarry.h"
#include "intarray.h"
#include "domain.h"
#include "mathfem.h"
#include "engngm.h"
#include "structuralms.h"
#include "load.h"
#ifndef __MAKEDEPEND
 #include <stdio.h>
#endif

#ifdef __OOFEG
 #include "oofeggraphiccontext.h"
 #include "conTable.h"
#endif

namespace oofem {
TrAxisym1_ht :: TrAxisym1_ht(int n, Domain *aDomain) : Tr1_ht(n, aDomain)
// Constructor.
{ }

TrAxisym1_ht :: ~TrAxisym1_ht()
// Destructor
{ }

double
TrAxisym1_ht :: computeVolumeAround(GaussPoint *gp)
{
    double determinant, weight;
    determinant = fabs( this->interp.giveTransformationJacobian(* gp->giveCoordinates(), FEIElementGeometryWrapper(this)) );
    weight = gp->giveWeight();
    return determinant * weight * this->computeRadiusAt(gp); ///@todo What about 2*pi ?
}

double
TrAxisym1_ht :: computeEdgeVolumeAround(GaussPoint *gp, int iEdge)
{
    FloatArray gcoords;
    double determinant, radius;

    determinant = fabs( this->interp.edgeGiveTransformationJacobian(iEdge, *gp->giveCoordinates(), FEIElementGeometryWrapper(this)) );
    this->interp.edgeLocal2global(gcoords, iEdge, *gp->giveCoordinates(), FEIElementGeometryWrapper(this));
    radius = gcoords.at(1);
    return determinant * radius * gp->giveWeight();
}

double
TrAxisym1_ht :: computeRadiusAt(GaussPoint *gp)
{
    FloatArray gcoords;
    this->computeGlobalCoordinates(gcoords, *gp->giveCoordinates());
    return gcoords.at(1);
}
} // end namespace oofem
