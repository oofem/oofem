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

#ifndef PLCZDAMAGERADIUS_H_
#define PLCZDAMAGERADIUS_H_

#include "xfem/propagationlaw.h"
#include "intarray.h"

#define _IFT_PLCZdamageRadius_Name "propagationlawczdamageradius"
#define _IFT_PLCZdamageRadius_IncRadius "incrementradius"           ///< Increment radius (from element nodes) per time step
#define _IFT_PLCZdamageRadius_DamageThreshold "damagethreshold"     ///< Damage threshold [0,1] for propagation
#define _IFT_PLCZdamageRadius_PropagationCS "propagationcs"         ///< Cross sections (must be part of csnum) viable for propagation

namespace oofem {
class Domain;
class EnrichmentDomain;
class DynamicInputRecord;


/**
 * Propagation law that propagates the (delamination) crack in a radius distance from element nodes
 * when the damage level in the associated cohesive zone reaces a defined value
 * Cracks w/o interface material as treated as fully damaged, thus will lead to propagation. 
 * cf. Främby, Fagerström & Bouzoulis, 'Adaptive modelling of delamination initiation and propagation using an equivalent single-layer shell approach', IJNME, 2016  
 *
 * @author Johannes Främby
 */
class OOFEM_EXPORT PLCZdamageRadius : public PropagationLaw
{
public:
    PLCZdamageRadius() : mIncrementRadius(0.0), mDamageThreshold(0.0), mPropCS(0) { }
    virtual ~PLCZdamageRadius() { }

    virtual const char *giveClassName() const { return "PLCZdamageRadius"; }
    virtual const char *giveInputRecordName() const { return _IFT_PLCZdamageRadius_Name; }

    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual void giveInputRecord(DynamicInputRecord &input);

    virtual bool hasPropagation() const { return mIncrementRadius > 0.; } ///@todo Could this be done smarter? / Mikael
    virtual bool propagateInterface(Domain &iDomain, EnrichmentFront &iEnrFront, TipPropagation &oTipProp);

    void setIncrementRadius(double iIncrementRadius) {mIncrementRadius = std::move(iIncrementRadius);}
    void setdamageThreshold(double idamageThreshold) {mDamageThreshold = std::move(idamageThreshold);}
    
    IntArray givePropagationCrossSections() {return this->mPropCS;}

protected:
    double mIncrementRadius, mDamageThreshold;
    IntArray mPropCS;
};
} // end namespace oofem


#endif /* PLHOOPSTRESSCIRC_H_ */
