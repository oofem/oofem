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

#include "xfem/inclusion.h"

#include "xfemmanager.h"
#include "element.h"
#include "feinterpol.h"
#include "gausspoint.h"

#include "classfactory.h"

#include <string>

namespace oofem {
REGISTER_EnrichmentItem(Inclusion)

Inclusion :: Inclusion(int n, XfemManager *xm, Domain *aDomain) :
    HybridEI(n, xm, aDomain),
    mpCrossSection(NULL)
{
    mpEnrichesDofsWithIdArray = {
        D_u, D_v, D_w
    };
}

Inclusion :: ~Inclusion()
{
    if ( mpCrossSection != NULL ) {
        mpCrossSection = NULL;
    }
}

bool Inclusion :: isMaterialModified(GaussPoint &iGP, Element &iEl, CrossSection * &opCS) const
{
    // Check if the point is located inside the inclusion

    FloatArray N;
    FEInterpolation *interp = iEl.giveInterpolation();
    interp->evalN( N, iGP.giveNaturalCoordinates(), FEIElementGeometryWrapper(& iEl) );

    const IntArray &elNodes = iEl.giveDofManArray();

    double levelSetGP = 0.0;
    interpLevelSet(levelSetGP, N, elNodes);

    if ( levelSetGP < 0.0 ) {
        opCS = mpCrossSection;
        return true;
    }

    return false;
}

IRResultType Inclusion :: initializeFrom(InputRecord *ir)
{
    IRResultType result;

    int crossSectionIndex = 0;
    IR_GIVE_FIELD(ir, crossSectionIndex, _IFT_Inclusion_CrossSection);
    mpCrossSection = this->giveDomain()->giveCrossSection(crossSectionIndex);

    return EnrichmentItem :: initializeFrom(ir);
}
} /* namespace oofem */
