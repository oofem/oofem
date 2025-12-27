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
 *               Copyright (C) 1993 - 2025   Borek Patzak
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

#ifndef SRC_SM_XFEM_PROPAGATIONLAWS_PLPRINCIPALSTRAIN_H_
#define SRC_SM_XFEM_PROPAGATIONLAWS_PLPRINCIPALSTRAIN_H_

#include "xfem/propagationlaw.h"

#define _IFT_PLPrincipalStrain_Name "propagationlawprincipalstrain"
#define _IFT_PLPrincipalStrain_Radius "radius" ///< Radius away from tip used when picking sampling point
//#define _IFT_PLHoopStressCirc_AngleInc "angleinc" ///< Angle between sampling points on the circle
#define _IFT_PLPrincipalStrain_IncLength "incrementlength" ///< Increment length per time step
#define _IFT_PLPrincipalStrain_StrainThreshold "strainthreshold" ///< Threshold for crack propagation
#define _IFT_PLPrincipalStrain_RadialBasisFunc "useradialbasisfunc" ///< If radial basis functions should be used for strain interpolation

namespace oofem {
class Domain;
class EnrichmentDomain;
class DynamicInputRecord;

class OOFEM_EXPORT PLPrincipalStrain : public PropagationLaw {
public:
    PLPrincipalStrain();
    virtual ~PLPrincipalStrain();

    const char *giveClassName() const override { return "PLPrincipalStrain"; }
    const char *giveInputRecordName() const override { return _IFT_PLPrincipalStrain_Name; }

    void initializeFrom(InputRecord &ir) override;
    void giveInputRecord(DynamicInputRecord &input) override;

    bool hasPropagation() const override { return mIncrementLength > 0.; }
    bool propagateInterface(Domain &iDomain, EnrichmentFront &iEnrFront, TipPropagation &oTipProp) override;

    void setRadius(double iRadius) {mRadius = std::move(iRadius);}
    void setIncrementLength(double iIncrementLength) {mIncrementLength = std::move(iIncrementLength);}
    void setStrainThreshold(double iStrainThreshold) {mStrainThreshold = std::move(iStrainThreshold);}
    void setUseRadialBasisFunc(bool iUseRadialBasisFunc) {mUseRadialBasisFunc = std::move(iUseRadialBasisFunc);}

protected:
    double mRadius, mIncrementLength, mStrainThreshold;
    bool mUseRadialBasisFunc;

};
} // end namespace oofem

#endif /* SRC_SM_XFEM_PROPAGATIONLAWS_PLPRINCIPALSTRAIN_H_ */
