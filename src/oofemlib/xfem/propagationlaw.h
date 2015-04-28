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


#ifndef PROPAGATIONLAW_H_
#define PROPAGATIONLAW_H_

#include "oofemcfg.h"
#include "inputrecord.h"

#define _IFT_PLDoNothing_Name "propagationlawdonothing"

#define _IFT_PLCrackPrescribedDir_Name "propagationlawcrackprescribeddir"
#define _IFT_PLCrackPrescribedDir_Dir "angle" // Angle in degrees
#define _IFT_PLCrackPrescribedDir_IncLength "incrementlength" // Increment per time step


namespace oofem {
class Domain;
class EnrichmentDomain;
class DynamicInputRecord;
class EnrichmentFront;
struct TipPropagation;

/**
 * Updates the geometry of evolving XFEM interfaces.
 * @author Erik Svenning
 */

class OOFEM_EXPORT PropagationLaw
{
public:
    PropagationLaw();
    virtual ~PropagationLaw();

    virtual const char *giveClassName() const = 0;
    virtual const char *giveInputRecordName() const = 0;

    virtual IRResultType initializeFrom(InputRecord *ir) = 0;
    virtual void giveInputRecord(DynamicInputRecord &input) = 0;

    virtual bool hasPropagation() const = 0;
    virtual bool propagateInterface(Domain &iDomain, EnrichmentFront &iEnrFront, TipPropagation &oTipProp) = 0;
};

/**
 * Dummy propagation law that does nothing.
 * @author Erik Svenning
 */
class OOFEM_EXPORT PLDoNothing : public PropagationLaw
{
public:
    PLDoNothing() { }
    virtual ~PLDoNothing() { }

    virtual const char *giveClassName() const { return "PLDoNothing"; }
    virtual const char *giveInputRecordName() const { return _IFT_PLDoNothing_Name; }

    virtual IRResultType initializeFrom(InputRecord *ir) { return IRRT_OK; }
    virtual void giveInputRecord(DynamicInputRecord &input);

    virtual bool hasPropagation() const { return false; }
    virtual bool propagateInterface(Domain &iDomain, EnrichmentFront &iEnrFront, TipPropagation &oTipProp) { return false; }
};

/**
 * Propagation law that propagates the crack in a predefined direction.
 * @author Erik Svenning
 */
class OOFEM_EXPORT PLCrackPrescribedDir : public PropagationLaw
{
public:
    PLCrackPrescribedDir() : mAngle(0.0), mIncrementLength(0.0) { }
    virtual ~PLCrackPrescribedDir() { }

    virtual const char *giveClassName() const { return "PLCrackPrescribedDir"; }
    virtual const char *giveInputRecordName() const { return _IFT_PLCrackPrescribedDir_Name; }

    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual void giveInputRecord(DynamicInputRecord &input);

    virtual bool hasPropagation() const { return mIncrementLength > 0.; }
    virtual bool propagateInterface(Domain &iDomain, EnrichmentFront &iEnrFront, TipPropagation &oTipProp);

protected:
    double mAngle, mIncrementLength;
};
} // end namespace oofem

#endif /* PROPAGATIONLAW_H_ */
