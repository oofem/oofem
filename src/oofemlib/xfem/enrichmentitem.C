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

#include "xfemmanager.h"
#include "flotmtrx.h"
#include "enrichmentitem.h"
#include "geometry.h"
#include "element.h"
#include "enrichmentfunction.h"
#include "cltypes.h"
#include "conTable.h"
#include "oofem_limits.h"
#include "usrdefsub.h"
#include "mathfem.h"

namespace oofem {
EnrichmentItem :: EnrichmentItem(int n, XfemManager *xm, Domain *aDomain) : FEMComponent(n, aDomain)
{
    xmanager = xm;
    geometry = 0;
    enrichmentFunction = 0;
}

BasicGeometry *EnrichmentItem :: giveGeometry()
{
    return xmanager->giveGeometry(this->geometry);
}

EnrichmentFunction *EnrichmentItem :: giveEnrichmentFunction()
{
    return xmanager->giveEnrichmentFunction(this->enrichmentFunction);
}

IRResultType EnrichmentItem :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result; // Required by IR_GIVE_FIELD macro

    this->geometry = 0;
    this->enrichmentFunction = 0;

    IR_GIVE_FIELD(ir, geometry, IFT_EnrichmentItem_geometryItemNr, "geometryitem"); // Macro
    IR_GIVE_FIELD(ir, enrichmentFunction, IFT_EnrichmentItem_enrichmentFunctionNr, "enrichmentfunction"); // Macro
    // this->setEnrichmentFunction(enrItemFunction);
    // this should go into enrichmentfunction probably
    // enrItemFunction->insertEnrichmentItem(this);
    // enrItemFunction->setActive(this);
    return IRRT_OK;
}

bool EnrichmentItem :: interacts(Element *element)
{
    return this->giveGeometry()->intersects(element);
}

bool EnrichmentItem :: isOutside(BasicGeometry *bg)
{
    return this->giveGeometry()->isOutside(bg);
}

void EnrichmentItem :: computeIntersectionPoints(AList< FloatArray > *intersectionPoints, Element *element)
{
    this->giveGeometry()->computeIntersectionPoints(element, intersectionPoints);
}

int EnrichmentItem :: computeNumberOfIntersectionPoints(Element *element)
{
    return this->giveGeometry()->computeNumberOfIntersectionPoints(element);
}

bool EnrichmentItem :: isDofManEnriched(int nodeNumber)
{
    bool ret = false;
    // gets neighbouring elements of a node
    const IntArray *neighbours = domain->giveConnectivityTable()->giveDofManConnectivityArray(nodeNumber);
    for ( int i = 1; i <= neighbours->giveSize(); i++ ) {
        // for each of the neighbouring elements finds out whether it interacts with this EnrichmentItem
        if ( this->interacts( domain->giveElement( neighbours->at(i) ) ) ) {
            ret = true;
            break;
        }
    }

    return ret;
}

IRResultType Inclusion :: initializeFrom(InputRecord *ir)
{
    EnrichmentItem :: initializeFrom(ir);
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result; // Required by IR_GIVE_FIELD macro
    int material = 0;
    IR_GIVE_FIELD(ir, material, IFT_EnrichmentItem_materialNr, "material"); // Macro
    this->mat = this->giveDomain()->giveMaterial(material);
    return IRRT_OK;
}
} // end namespace oofem
