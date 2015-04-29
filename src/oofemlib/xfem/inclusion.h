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

#ifndef INCLUSIONEI_H_
#define INCLUSIONEI_H_

#define _IFT_Inclusion_Name "inclusion"
#define _IFT_Inclusion_CrossSection "crosssection"

#include "xfem/enrichmentitem.h"
#include "xfem/hybridei.h"

namespace oofem {
/**
 * EnrichmentItem describing an inclusion
 * @author Erik Svenning (among others)
 * @date Sep 9, 2014
 */
class OOFEM_EXPORT Inclusion : public HybridEI
{
protected:
    CrossSection *mpCrossSection;
public:
    Inclusion(int n, XfemManager *xm, Domain *aDomain);
    virtual ~Inclusion();

    // Returns true if the enrichment item can assign
    // a different material to any Gauss point.
    inline virtual bool canModifyMaterial() const { return true; }

    // Returns true if the enrichment item assigns a different material to the Gauss point
    virtual bool isMaterialModified(GaussPoint &iGP, Element &iEl, CrossSection * &opCS) const;


    virtual const char *giveClassName() const { return "Inclusion"; }
    virtual const char *giveInputRecordName() const { return _IFT_Inclusion_Name; }
    virtual IRResultType initializeFrom(InputRecord *ir);
    CrossSection *giveCrossSection() { return mpCrossSection; }
};
} /* namespace oofem */

#endif /* INCLUSIONEI_H_ */
