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

#ifndef patchintegrationrule_h
#define patchintegrationrule_h

#include "gaussintegrationrule.h"
#include "patch.h"

namespace oofem {

/**
 * Represents an IntegrationRule for an interacted element
 * the standard integration is replaced by an integration over a
 * patchset.
 */
class PatchIntegrationRule : public GaussIntegrationRule
{
protected:
    /// Patch.
    Patch *patch;

public:
    /// Constructor.
    PatchIntegrationRule(int n, Element *e, Patch *p);
    /// Destructor.
    virtual ~PatchIntegrationRule();
    virtual int SetUpPointsOnTriangle(int, MaterialMode);
    int giveMaterial() { return this->patch->giveMaterial(); }
    Patch *givePatch() { return this->patch; }

    virtual contextIOResultType saveContext(DataStream *stream, ContextMode mode, void *obj);
    virtual contextIOResultType restoreContext(DataStream *stream, ContextMode mode, void *obj);
    virtual classType giveClassID() const { return PatchIntegrationRuleClass; }
};
} // end namespace oofem
#endif // patchintegrationrule_h
