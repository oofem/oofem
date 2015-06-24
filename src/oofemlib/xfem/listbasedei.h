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

#ifndef LISTBASEDEI_H_
#define LISTBASEDEI_H_

#define _IFT_ListBasedEI_Name "listbasedei"
#define _IFT_ListBasedEI_list "list"

#include "enrichmentitem.h"

#include <vector>

namespace oofem {
class XfemManager;
class Domain;

/**
 * EnrichmentItem with geometry defined by a set of nodes to be enriched.
 * @author Erik Svenning
 * @date Sep 9, 2014
 */
class OOFEM_EXPORT ListBasedEI : public EnrichmentItem
{
public:
    ListBasedEI(int n, XfemManager *xm, Domain *aDomain);
    virtual ~ListBasedEI();

    virtual const char *giveClassName() const { return "ListBasedEI"; }
    virtual const char *giveInputRecordName() const { return _IFT_ListBasedEI_Name; }

    virtual void updateGeometry();
    virtual void propagateFronts(bool &oFrontsHavePropagated);
    virtual void updateNodeEnrMarker(XfemManager &ixFemMan);

    virtual bool giveElementTipCoord(FloatArray &oCoord, double &oArcPos,  Element &iEl, const FloatArray &iElCenter) const;

    virtual void giveBoundingSphere(FloatArray &oCenter, double &oRadius) { OOFEM_ERROR("Not implemented.") }

protected:
    std :: vector< int >dofManList;
    double xi;
    int setNumber;
};
} /* namespace oofem */

#endif /* LISTBASEDEI_H_ */
