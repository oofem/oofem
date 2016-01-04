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

#ifndef PLHOOPSTRESSCIRC_H_
#define PLHOOPSTRESSCIRC_H_

#include "xfem/propagationlaw.h"

#define _IFT_PLHoopStressCirc_Name "propagationlawhoopstresscirc"
#define _IFT_PLHoopStressCirc_Radius "radius" ///< Radius of circle used for stress sampling points
#define _IFT_PLHoopStressCirc_AngleInc "angleinc" ///< Angle between sampling points on the circle
#define _IFT_PLHoopStressCirc_IncLength "incrementlength" ///< Increment length per time step
#define _IFT_PLHoopStressCirc_HoopStressThreshold "hoopstressthreshold" ///< Threshold for crack propagation
#define _IFT_PLHoopStressCirc_RadialBasisFunc "useradialbasisfunc" ///< If radial basis functions should be used for stress interpolation

namespace oofem {
class Domain;
class EnrichmentDomain;
class DynamicInputRecord;


/**
 * Propagation law that propagates the crack in the direction
 * that gives @f$ \sigma_{r\theta} = 0 @f$.
 * Based on
 * T.P. Fries and M. Baydoun:
 * "Crack propagation with the extended finite element method
 * and a hybrid explicit-implicit crack description",
 * Internat. J. Numer. Methods Engrg 89,
 * pp. 1527--1558 (2012)
 *
 * The stress is evaluated in several points on a circle
 * surrounding the crack tip.
 *
 * Compared to the paper above, the implementation has been extended
 * with a criterion for crack propagation instead of always
 * propagating a predefined increment length. Two options are
 * currently available for stress interpolation:
 * 1) Take stress of closest Gauss point
 * 2) Interpolate with radial basis functions
 *
 * @author Erik Svenning
 */
class OOFEM_EXPORT PLHoopStressCirc : public PropagationLaw
{
public:
    PLHoopStressCirc() : mRadius(0.0), mAngleInc(0.0), mIncrementLength(0.0), mHoopStressThreshold(0.0), mUseRadialBasisFunc(false) { }
    virtual ~PLHoopStressCirc() { }

    virtual const char *giveClassName() const { return "PLHoopStressCirc"; }
    virtual const char *giveInputRecordName() const { return _IFT_PLHoopStressCirc_Name; }

    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual void giveInputRecord(DynamicInputRecord &input);

    virtual bool hasPropagation() const { return mIncrementLength > 0.; } ///@todo Could this be done smarter? / Mikael
    virtual bool propagateInterface(Domain &iDomain, EnrichmentFront &iEnrFront, TipPropagation &oTipProp);

protected:
    double mRadius, mAngleInc, mIncrementLength, mHoopStressThreshold;
    bool mUseRadialBasisFunc;
};
} // end namespace oofem


#endif /* PLHOOPSTRESSCIRC_H_ */
