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

#include "sm/Materials/InterfaceMaterials/structuralinterfacematerialphf.h"
#include "intmatphasefield.h"

#include "gausspoint.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "mathfem.h"
#include "datastream.h"
#include "contextioerr.h"
#include "classfactory.h"
#include "dynamicinputrecord.h"

namespace oofem {
REGISTER_Material(IntMatPhaseField);

IntMatPhaseField::IntMatPhaseField(int n, Domain *d) : StructuralInterfaceMaterialPhF(n, d) {}


bool
IntMatPhaseField :: hasMaterialModeCapability(MaterialMode mode) const
{
    // returns whether receiver supports given mode
    return mode == _3dInterface;
}

FloatArrayF<3>
IntMatPhaseField :: giveEngTraction_3d(const FloatArrayF<3> &jump, double damage, GaussPoint *gp, TimeStep *tStep) const
{
    IntMatPhaseFieldStatus *status = static_cast< IntMatPhaseFieldStatus * >( this->giveStatus(gp) );

    //status->initTempStatus();

    double g = compute_g(damage);

    FloatArrayF<3> answer;
    answer.at(1) = g * k*jump.at(1);
    answer.at(2) = g * k*jump.at(2);

    if ( jump.at(3) > 0.0 ) { // only degradation in tension
        answer.at(3) = g * k*jump.at(3);
    } else {
        answer.at(3) = k*jump.at(3);
    }

    double drivingEnergy = 0.5*k*jump.at(1)*jump.at(1) + 0.5*k*jump.at(2)*jump.at(2);
    if ( jump.at(3) > 0.0 ) { 
        drivingEnergy += 0.5*k*jump.at(3)*jump.at(3);
    }
    if ( drivingEnergy > status->giveTempDrivingEnergy() ) { // max val
        status->letTempDrivingEnergyBe( drivingEnergy ); 
    }

    status->tempDamage = damage; 

    status->letTempJumpBe(jump);
    status->letTempTractionBe(answer);

    return answer;
}

FloatMatrixF<3,3>
IntMatPhaseField :: give3dStiffnessMatrix_Eng(MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) const
{
    IntMatPhaseFieldStatus *status = static_cast< IntMatPhaseFieldStatus * >( this->giveStatus(gp) );

    const auto &jump = status->giveTempJump();

    double damage = status->giveDamage();
    double g = compute_g(damage);

    FloatMatrixF<3,3> answer;
    answer.at(1, 1) = g*k;
    answer.at(2, 2) = g*k;

    if ( jump.at(3) > 0.0 ) { // only degradation in tension
        answer.at(3,3) = g * k;
    } else {
        answer.at(3,3) = k;
    }
    return answer;
}


void
IntMatPhaseField :: giveTangents(FloatMatrix &Djj, FloatMatrix &Djd, FloatMatrix &Ddj, FloatMatrix &Ddd, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const
{ 
    IntMatPhaseFieldStatus *status = static_cast< IntMatPhaseFieldStatus * >( this->giveStatus(gp) );

    const auto &jump = status->giveTempJump();

    double damage = status->giveDamage();
    double g = compute_g(damage);

    Djj.resize(3, 3);
    Djj.zero();
    Djj.at(1, 1) = g*k;
    Djj.at(2, 2) = g*k;

    if ( jump.at(3) > 0.0 ) { // only degradation in tension
        Djj.at(3,3) = g * k;
    } else {
        Djj.at(3,3) = k;
    }

    double gPrime = compute_gPrime(damage);
    Djd.at(1, 1) = gPrime*k * jump.at(1);
    Djd.at(2, 2) = gPrime*k * jump.at(1);
    if ( jump.at(3) > 0.0 ) { // only degradation in tension
        Djd.at(3,3) = gPrime * k * jump.at(1);
    } else {
        Djd.at(3,3) = 0;
    }

    // Df/Dd, f= max[ g'*(psi_s+psi_n+) ]
    Ddj.beTranspositionOf(Djd);
    //Ddd = { g''* }
}

double
IntMatPhaseField :: giveDrivingForce(GaussPoint *gp) const
{
    IntMatPhaseFieldStatus *status = static_cast< IntMatPhaseFieldStatus * >( this->giveStatus(gp) );
    double gPrime = compute_gPrime(status->giveDamage());

    return gPrime / this->Gc * status->giveTempDrivingEnergy(); 
}

double
IntMatPhaseField :: giveDrivingForcePrime(GaussPoint *gp) const
{
    IntMatPhaseFieldStatus *status = static_cast< IntMatPhaseFieldStatus * >( this->giveStatus(gp) );    
    double gBis = compute_gBis(status->giveDamage());

    return gBis / this->Gc * status->giveTempDrivingEnergy();     
}



// double
// IntMatPhaseField :: computeFreeEnergy(GaussPoint *gp, TimeStep *tStep)
// {
//     //StructuralMaterialStatus *matStat = static_cast< StructuralMaterialStatus * >( gp->giveMaterialStatus() );
//     FloatArray strain, stress;
//     //stress = matStat->giveTempStressVector();
//     //strain = matStat->giveTempStrainVector();
//     //stress = matStat->giveStressVector();
//     //strain = matStat->giveStrainVector();
//     return 0.5 * stress.dotProduct( strain );
// }



double
IntMatPhaseField :: compute_g(double d) const
{
    // computes g = (1-d)^2 + r0
    double r0 = 1.0e-10;
    return (1.0 - d) * (1.0 - d) + r0;
}


double
IntMatPhaseField :: compute_gPrime(double d) const
{
    // compute Dg/Dd = -2*(1-d)
    return -2.0 * (1.0 - d);
}

double
IntMatPhaseField :: compute_gBis(double d) const
{
    // compute DDg/DDd = 2
    return 2.0;
}


void
IntMatPhaseField :: initializeFrom(InputRecord &ir)
{
    IR_GIVE_FIELD(ir, this->k, _IFT_IntMatPhaseField_kn);
    IR_GIVE_FIELD(ir, this->Gc, _IFT_IntMatPhaseField_gc);

    StructuralInterfaceMaterial :: initializeFrom(ir);
}


void IntMatPhaseField :: giveInputRecord(DynamicInputRecord &input)
{
    StructuralInterfaceMaterial :: giveInputRecord(input);

    input.setField(k, _IFT_IntMatPhaseField_kn);
}




IntMatPhaseFieldStatus :: IntMatPhaseFieldStatus(GaussPoint *g) : StructuralInterfaceMaterialStatus(g)
{
}


void
IntMatPhaseFieldStatus :: initTempStatus()
{
    StructuralInterfaceMaterialStatus :: initTempStatus();
    this->tempDrivingEnergy = 0.0;
    this->drivingEnergy = 0.0;
    this->tempDamage = 0.0;
}


void
IntMatPhaseFieldStatus :: updateYourself(TimeStep *atTime)
{
    StructuralInterfaceMaterialStatus :: updateYourself(atTime);
    this->drivingEnergy = this->tempDrivingEnergy;
}


} /* namespace oofem */
