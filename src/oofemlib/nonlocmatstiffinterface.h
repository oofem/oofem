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

#ifndef nonlocmatstiffinterface_h
#define nonlocmatstiffinterface_h

#include "interface.h"
#include "nonlocalmaterialext.h"

namespace oofem {
class SparseMtrx;
class GaussPoint;
class TimeStep;
class oofegGraphicContext;
class UnknownNumberingScheme;

/**
 * Class Nonlocal Material Stiffness Interface. This is only abstract class.
 * This interface allows material model to add those services required to
 * compute and assemble nonlocal contribution to stiffness matrix.
 */
class OOFEM_EXPORT NonlocalMaterialStiffnessInterface : public Interface
{
public:
    /// Constructor.
    NonlocalMaterialStiffnessInterface() : Interface() { }

    /// Computes and adds IP contributions to destination matrix.
    virtual void NonlocalMaterialStiffnessInterface_addIPContribution(SparseMtrx &dest, const UnknownNumberingScheme &s,
                                                                      GaussPoint *gp, TimeStep *tStep) = 0;
    /**
     * Returns integration list of receiver. Contains localIntegrationRecord structures, containing
     * references to integration points and their weights that influence to nonlocal average in
     * receiver's associated integration point.
     */
    virtual std :: list< localIntegrationRecord > *NonlocalMaterialStiffnessInterface_giveIntegrationDomainList(GaussPoint *gp) = 0;

#ifdef __OOFEG
    /**
     * Plots the sparse structure of stiffness contribution.
     */
    virtual void NonlocalMaterialStiffnessInterface_showSparseMtrxStructure(GaussPoint *gp, oofegGraphicContext &gc, TimeStep *) { }
#endif
};
} // end namespace oofem
#endif // nonlocmatstiffinterface_h
