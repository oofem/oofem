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

#ifndef INTMATELASTIC_H_
#define INTMATELASTIC_H_

#include "structuralinterfacematerial.h"
#include "structuralinterfacematerialstatus.h"

///@name Input fields for IntMatElastic
//@{
#define _IFT_IntMatElastic_Name "intmatelastic"
#define _IFT_IntMatElastic_kn "k"
//@}


namespace oofem {

/**
 * Linear elastic cohesive zone.
 * @author Erik Svenning
 */
class IntMatElastic : public StructuralInterfaceMaterial
{
public:
	IntMatElastic(int n, Domain * d);
	virtual ~IntMatElastic();

    virtual int hasNonLinearBehaviour()   { return 0; }

    virtual const char *giveClassName() const { return "IntMatElastic"; }
    virtual const char *giveInputRecordName() const { return _IFT_IntMatElastic_Name; }

    virtual void giveFirstPKTraction_3d(FloatArray &answer, GaussPoint *gp, const FloatArray &jumpVector,
                                        const FloatMatrix &F, TimeStep *tStep);

    virtual void give3dStiffnessMatrix_dTdj(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep);

    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual void giveInputRecord(DynamicInputRecord &input);

    virtual MaterialStatus *CreateStatus(GaussPoint *gp) const { return new StructuralInterfaceMaterialStatus(1, domain, gp); }

    virtual bool hasAnalyticalTangentStiffness() const { return true; }

protected:
	double k;
};

} /* namespace oofem */
#endif /* INTMATELASTIC_H_ */
