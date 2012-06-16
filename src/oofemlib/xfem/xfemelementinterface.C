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

#include "xfemelementinterface.h"
#include "enrichmentitem.h"
#include "engngm.h"
#include "gausspnt.h"
#include "materialmode.h"
#include "fei2dquadlin.h"
#include "patch.h"
#include "patchintegrationrule.h"
#include "delaunay.h"
#include "xfemmanager.h"
#include "flotarry.h"
#include "flotmtrx.h"

namespace oofem {
void XfemElementInterface :: XfemElementInterface_partitionElement(AList< Triangle > *answer, AList< FloatArray > *together)
{
    Delaunay dl;
    dl.triangulate(together, answer);
}

void XfemElementInterface :: XfemElementInterface_updateIntegrationRule()
{
    XfemManager *xf = this->element->giveDomain()->giveEngngModel()->giveXfemManager(1);
    if ( xf->isInteracted(element) ) {
        IntArray interactedEI;
        xf->getInteractedEI(interactedEI, element);
        AList< Triangle >triangles;
        AList< Triangle >triangles2;
        // all the points coming into triangulation
        AList< FloatArray >together1;
        AList< FloatArray >together2;
        this->XfemElementInterface_prepareNodesForDelaunay(& together1, & together2);
        this->XfemElementInterface_partitionElement(& triangles, & together1);
        this->XfemElementInterface_partitionElement(& triangles2, & together2);

        for ( int i = 1; i <= triangles2.giveSize(); i++ ) {
            int sz = triangles.giveSize();
            triangles.put( sz + 1, triangles2.at(i) );
            triangles2.unlink(i);
        }

        AList< IntegrationRule >irlist;
        for ( int i = 1; i <= triangles.giveSize(); i++ ) {
            int mat = 0;
            if ( xf->giveEnrichmentItem( interactedEI.at(1) )->giveGeometry()->isOutside( triangles.at(i) ) ) {
                mat = 1;
            } else   {
                mat = 2;
            }

            Patch *patch = new TrianglePatch(element, mat);
            for ( int j = 1; j <= triangles.at(i)->giveVertices()->giveSize(); j++ ) {
                FloatArray *nCopy = new FloatArray( *triangles.at( i )->giveVertex(j) );
                patch->setVertex(nCopy);
            }

            PatchIntegrationRule *pir = new PatchIntegrationRule(i, element, patch);
            int pointNr = 3;
            MaterialMode matMode = element->giveDefaultIntegrationRulePtr()->getIntegrationPoint(0)->giveMaterialMode();
            pir->setUpIntegrationPoints(_Triangle, pointNr, matMode);
            irlist.put(i, pir);
        }

        element->setIntegrationRules(& irlist);
    }
}

void XfemElementInterface :: XfemElementInterface_prepareNodesForDelaunay(AList< FloatArray > *answer1, AList< FloatArray > *answer2)
{
    XfemManager *xf = this->element->giveDomain()->giveEngngModel()->giveXfemManager(1);
    IntArray interactedEI;
    xf->getInteractedEI(interactedEI, element);
    // in intersecPoints the points of Element with interaction to EnrichmentItem will be stored
    AList< FloatArray >intersecPoints;
    for ( int i = 1; i <= interactedEI.giveSize(); i++ ) {
        xf->giveEnrichmentItem( interactedEI.at(i) )->computeIntersectionPoints(& intersecPoints, element);
    }

    // here the intersection points are copied in order to be put into two groups
    for ( int i = 1; i <= intersecPoints.giveSize(); i++ ) {
        int sz = answer1->giveSize();
        answer1->put( sz + 1, intersecPoints.at(i) );
        FloatArray *ip = intersecPoints.at(i);
        FloatArray *ipCopy = new FloatArray(*ip);
        int sz2 = answer2->giveSize();
        answer2->put(sz2 + 1, ipCopy);
    }

    if ( intersecPoints.giveSize() == 2 ) {
        // here the group is determined
        double x1 = intersecPoints.at(1)->at(1);
        double x2 = intersecPoints.at(2)->at(1);
        double y1 = intersecPoints.at(1)->at(2);
        double y2 = intersecPoints.at(2)->at(2);
        for ( int i = 1; i <= this->element->giveNumberOfDofManagers(); i++ ) {
            double x = element->giveDofManager(i)->giveCoordinates()->at(1);
            double y = element->giveDofManager(i)->giveCoordinates()->at(2);
            double det = ( x1 - x ) * ( y2 - y ) - ( x2 - x ) * ( y1 - y );
            FloatArray *node = element->giveDofManager(i)->giveCoordinates();
            FloatArray *nodesCopy = new FloatArray(*node);
            if ( det > 0.00001 ) {
                int sz = answer1->giveSize();
                answer1->put(sz + 1, nodesCopy);
            } else if ( det < ( -1 ) * 0.00001 ) {
                int sz = answer2->giveSize();
                answer2->put(sz + 1, nodesCopy);
            }
        }
    }

    // nodes of an element are copied to a different memory location
    // so that the whole container of points for triangulation can be dealt with
    // more easily (e.g. deleted)

    for ( int i = 1; i <= intersecPoints.giveSize(); i++ ) {
        intersecPoints.unlink(i);
    }
}
} // end namespace oofem
