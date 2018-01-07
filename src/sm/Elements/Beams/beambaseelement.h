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

#ifndef beambaseelement_h
#define beambaseelement_h

#include "sm/Elements/structuralelement.h"
#include "sm/CrossSections/fiberedcs.h"
#include "sm/Materials/winklermodel.h"
#include "dofmanager.h"

namespace oofem {

/**
 * This class implements a base beam intented to be a base class for
 * beams based on lagrangian interpolation, where exact end forces can
 * be recovered.
 */
 class BeamBaseElement : public StructuralElement
{
protected:

 public:
    BeamBaseElement (int n, Domain *d);
    virtual ~BeamBaseElement();

protected:
    /** Computes element end force vector from applied loading in local coordinate system
     * @param answer computed end force vector due to non-nodal loading
     * @param tStep solution step
     * @param mode determines response mode
     */
    virtual void computeLocalForceLoadVector(FloatArray &answer, TimeStep *tStep, ValueModeType mode);

};
} // end namespace oofem
#endif // beambaseelement_h
