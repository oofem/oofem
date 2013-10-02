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

#include "inputrecord.h"

#define _IFT_PLDoNothing_Name "propagationlawdonothing"

#define _IFT_PLCrackPrescribedDir_Name "propagationlawcrackprescribeddir"
#define _IFT_PLCrackPrescribedDir_Dir "angle" // Angle in degrees
#define _IFT_PLCrackPrescribedDir_IncLength "incrementlength" // Increment per time step



namespace oofem {

class EnrichmentDomain;

/**
 * Updates the geometry of evolving XFEM interfaces.
 * @author Erik Svenning
 */

class PropagationLaw {
public:
	PropagationLaw();
	virtual ~PropagationLaw();

    virtual const char *giveClassName() const = 0;
    virtual const char *giveInputRecordName() const = 0;

    virtual IRResultType initializeFrom(InputRecord *ir) = 0;

	virtual void propagateInterfaces(EnrichmentDomain &iEnrDom) = 0;
};

class PLDoNothing: public PropagationLaw {
public:
	PLDoNothing() {};
	virtual ~PLDoNothing() {};

    virtual const char *giveClassName() const { return "PLDoNothing"; }
    virtual const char *giveInputRecordName() const { return _IFT_PLDoNothing_Name; }

    virtual IRResultType initializeFrom(InputRecord *ir) {return IRRT_OK;}

	virtual void propagateInterfaces(EnrichmentDomain &ioEnrDom) {};
};

class PLCrackPrescribedDir: public PropagationLaw {
public:
	PLCrackPrescribedDir():mAngle(0.0), mIncrementLength(0.0) {};
	virtual ~PLCrackPrescribedDir() {};

    virtual const char *giveClassName() const { return "PLCrackPrescribedDir"; }
    virtual const char *giveInputRecordName() const { return _IFT_PLCrackPrescribedDir_Name; }

    virtual IRResultType initializeFrom(InputRecord *ir);

	virtual void propagateInterfaces(EnrichmentDomain &ioEnrDom);

protected:
	double mAngle, mIncrementLength;
};

} // end namespace oofem

#endif /* PROPAGATIONLAW_H_ */
