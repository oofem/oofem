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

/*
 * stvenantkirchhoff.C
 *
 *  Created on: May 24, 2013
 *      Author: svennine
 */

#include "stvenantkirchhoff.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "classfactory.h"


namespace oofem {

REGISTER_Material( StVenantKirchhoff );

StVenantKirchhoff :: StVenantKirchhoff(int n, Domain *d) : StructuralMaterial(n, d)
{}

StVenantKirchhoff::~StVenantKirchhoff()
{

}

int
StVenantKirchhoff :: hasMaterialModeCapability(MaterialMode mode)
//
// returns whether receiver supports the given mode
//
{
    if ( mode == _3dMat ) {
        return 1;
    }

    return 0;
}


void
StVenantKirchhoff :: give3dMaterialStiffnessMatrix(FloatMatrix &answer, MatResponseMode, GaussPoint *gp, TimeStep *atTime)

// returns the 6x6 tangent stiffness matrix

{

    double c1 = Lambda + 2.0*Mu;

    answer.resize(6, 6);

    answer.at(1, 1) = c1;
    answer.at(2, 2) = c1;
    answer.at(3, 3) = c1;
    answer.at(4, 4) = 1.0*Mu;
    answer.at(5, 5) = 1.0*Mu;
    answer.at(6, 6) = 1.0*Mu;
    answer.at(1, 2) = answer.at(2, 1) = Lambda;
    answer.at(1, 3) = answer.at(3, 1) = Lambda;
    answer.at(1, 4) = answer.at(4, 1) = 0.0;
    answer.at(1, 5) = answer.at(5, 1) = 0.0;
    answer.at(1, 6) = answer.at(6, 1) = 0.0;
    answer.at(2, 3) = answer.at(3, 2) = Lambda;
    answer.at(2, 4) = answer.at(4, 2) = 0.0;
    answer.at(2, 5) = answer.at(5, 2) = 0.0;
    answer.at(2, 6) = answer.at(6, 2) = 0.0;
    answer.at(3, 4) = answer.at(4, 3) = 0.0;
    answer.at(3, 5) = answer.at(5, 3) = 0.0;
    answer.at(3, 6) = answer.at(6, 3) = 0.0;
    answer.at(4, 5) = answer.at(5, 4) = 0.0;
    answer.at(4, 6) = answer.at(6, 4) = 0.0;
    answer.at(5, 6) = answer.at(6, 5) = 0.0;


}


void
StVenantKirchhoff :: giveRealStressVector(FloatArray &answer, MatResponseForm form, GaussPoint *gp, const FloatArray &totalStrain, TimeStep *atTime)

// returns 6 components of the stress corresponding to the given total strain

{

    FloatArray strainVector;

    StVenantKirchhoffMaterialStatus *status = static_cast< StVenantKirchhoffMaterialStatus * >( this->giveStatus(gp) );
    this->giveStressDependentPartOfStrainVector(strainVector, gp,
                                                totalStrain,
                                                atTime, VM_Total);

    answer.resize(6);
    double c1 = Lambda*( strainVector.at(1) + strainVector.at(2) + strainVector.at(3) );

    answer.at(1) = c1 + 2.0*Mu*strainVector.at(1);
    answer.at(2) = c1 + 2.0*Mu*strainVector.at(2);
    answer.at(3) = c1 + 2.0*Mu*strainVector.at(3);

    answer.at(4) = 		1.0*Mu*strainVector.at(4);
    answer.at(5) = 		1.0*Mu*strainVector.at(5);
    answer.at(6) = 		1.0*Mu*strainVector.at(6);

    // update gp
    status->letTempStrainVectorBe(totalStrain);
    status->letTempStressVectorBe(answer);

}


MaterialStatus *
StVenantKirchhoff :: CreateStatus(GaussPoint *gp) const
{
    StructuralMaterialStatus *status;

    status = new StructuralMaterialStatus(1, this->giveDomain(), gp);
    return status;
}


IRResultType
StVenantKirchhoff :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                // Required by IR_GIVE_FIELD macro


    StructuralMaterial :: initializeFrom(ir);

    // Read material properties here

    IR_GIVE_FIELD(ir, Lambda, _IFT_StVenantKirchhoff_lambda);
    IR_GIVE_FIELD(ir, Mu, _IFT_StVenantKirchhoff_mu);

    return IRRT_OK;
}


StVenantKirchhoffMaterialStatus :: StVenantKirchhoffMaterialStatus(int n, Domain *d, GaussPoint *g) : StructuralMaterialStatus(n, d, g)
{
    // init state variables
}


StVenantKirchhoffMaterialStatus :: ~StVenantKirchhoffMaterialStatus()
{}


void
StVenantKirchhoffMaterialStatus :: printOutputAt(FILE *file, TimeStep *tStep)
{
    // print state to output stream

    StructuralMaterialStatus :: printOutputAt(file, tStep);
    fprintf(file, "status { ");
    fprintf(file, "}\n");
}

// initialize temporary state variables according to equilibrated state vars
void
StVenantKirchhoffMaterialStatus :: initTempStatus()
{
    StructuralMaterialStatus :: initTempStatus();
}


// Called when equilibrium reached, set equilibrated vars according to temporary (working) ones.
void
StVenantKirchhoffMaterialStatus :: updateYourself(TimeStep *atTime)
{
    StructuralMaterialStatus :: updateYourself(atTime);
}
} // end namespace oofem
