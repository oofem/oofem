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

#ifndef PLMATERIALFORCE_H_
#define PLMATERIALFORCE_H_

#include "xfem/propagationlaw.h"

#include <memory>

#define _IFT_PLMaterialForce_Name "propagationlawmaterialforce"
#define _IFT_PLMaterialForce_Radius "radius" ///< Radius of region for domain integral
#define _IFT_PLMaterialForce_IncLength "incrementlength" ///< Increment length per time step
#define _IFT_PLMaterialForce_CrackPropThreshold "gc" ///<  Threshold for crack propagation

namespace oofem {

class Domain;
class EnrichmentDomain;
class DynamicInputRecord;
class MaterialForceEvaluator;


/**
 * Propagation law that propagates the crack in
 * the direction of the material force.
 *
 * @author Erik Svenning
 * @date Nov 14, 2014
 */
class OOFEM_EXPORT PLMaterialForce : public PropagationLaw
{
public:
    PLMaterialForce();
    virtual ~PLMaterialForce();

    const char *giveClassName() const override { return "PLMaterialForce"; }
    const char *giveInputRecordName() const override { return _IFT_PLMaterialForce_Name; }

    void initializeFrom(InputRecord &ir) override;
    void giveInputRecord(DynamicInputRecord &input) override;

    bool hasPropagation() const override { return mIncrementLength > 0.; } ///@todo Could this be done smarter? / Mikael
    bool propagateInterface(Domain &iDomain, EnrichmentFront &iEnrFront, TipPropagation &oTipProp) override;

    void setRadius(const double &iRadius) {mRadius = iRadius;}
    void setIncrementLength(const double &iIncrementLength) {mIncrementLength = iIncrementLength;}
    void setCrackPropThreshold(const double &iCrackPropThreshold) {mCrackPropThreshold = iCrackPropThreshold;}

protected:
    double mRadius, mIncrementLength, mCrackPropThreshold;

    std :: unique_ptr< MaterialForceEvaluator > mpMaterialForceEvaluator;
};

} /* namespace oofem */

#endif /* PLMATERIALFORCE_H_ */
