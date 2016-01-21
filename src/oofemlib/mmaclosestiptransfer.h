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

#ifndef mmaclosestiptransfer_h
#define mmaclosestiptransfer_h

#include "materialmappingalgorithm.h"

namespace oofem {
class Domain;
class Element;
class TimeStep;
class MaterialStatus;

/**
 * The class implements the closest integration point transfer of state variables.
 */
class OOFEM_EXPORT MMAClosestIPTransfer : public MaterialMappingAlgorithm
{
protected:
    GaussPoint *source;
    MaterialStatus *mpMaterialStatus;

public:
    /// Constructor
    MMAClosestIPTransfer();

    virtual void __init(Domain *dold, IntArray &type, FloatArray &coords, Set &sourceElemSet, TimeStep *tStep, bool iCohesiveZoneGP = false);

    virtual void finish(TimeStep *tStep) { }

    virtual int __mapVariable(FloatArray &answer, FloatArray &coords, InternalStateType type, TimeStep *tStep);

    virtual int mapStatus(MaterialStatus &oStatus) const;

    virtual const char *giveClassName() const { return "MMAClosestIPTransfer"; }
};
} // end namespace oofem
#endif // mmaclosestiptransfer_h
