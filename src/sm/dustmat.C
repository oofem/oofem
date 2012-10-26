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
 *               Copyright (C) 1993 - 2012   Borek Patzak
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

#include "dustmat.h"

#include "flotarry.h"
#include "flotmtrx.h"
#include "structuralms.h"
#include "gausspnt.h"
#include "intarray.h"
#include "structuralmaterial.h"
#include "isolinearelasticmaterial.h"
#include "structuralcrosssection.h"
#include "datastream.h"
#include "contextioerr.h"
#include "mathfem.h"

namespace oofem {
DustMaterialStatus :: DustMaterialStatus(int n, Domain *d, GaussPoint *gp) :
    StructuralMaterialStatus(n, d, gp),
    plasticStrainDeviator( gp->giveMaterialMode() ),
    tempPlasticStrainDeviator( gp->giveMaterialMode() )
{
	q = 0.;
}

DustMaterialStatus :: ~DustMaterialStatus()
{ }

void
DustMaterialStatus :: initTempStatus()
{
    // Call the function of the parent class to initialize the variables defined there.
    StructuralMaterialStatus :: initTempStatus();
    // tempVal = val
	 tempVolumetricPlasticStrain = volumetricPlasticStrain;
	 tempPlasticStrainDeviator = plasticStrainDeviator;
	 tempQ = q;
}

void
DustMaterialStatus :: updateYourself(TimeStep *atTime)
{
    // Call corresponding function of the parent class to update variables defined there.
    StructuralMaterialStatus :: updateYourself(atTime);
    // val = tempVal
	 volumetricPlasticStrain = tempVolumetricPlasticStrain;
	 plasticStrainDeviator = tempPlasticStrainDeviator;
	 q = tempQ;
}

void
DustMaterialStatus :: printOutputAt(FILE *file, TimeStep *tStep)
{
    // Call corresponding function of the parent class to print variables defined there.
    StructuralMaterialStatus :: printOutputAt(file, tStep);

    fprintf(file, "\tstatus { ");
    fprintf(file, "}\n");
}

contextIOResultType
DustMaterialStatus :: saveContext(DataStream *stream, ContextMode mode, void *obj)
{
    contextIOResultType iores;

    // save parent class status
    if ( ( iores = StructuralMaterialStatus :: saveContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    return CIO_OK;
}


contextIOResultType
DustMaterialStatus :: restoreContext(DataStream *stream, ContextMode mode, void *obj)
{
    contextIOResultType iores;

    // read parent class status
    if ( ( iores = StructuralMaterialStatus :: restoreContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    return CIO_OK;
}

//   *************************************************************
//   *** CLASS DUST MATERIAL   ***
//   *************************************************************


DustMaterial :: DustMaterial(int n, Domain *d) : StructuralMaterial(n, d)
{
    LEMaterial = new IsotropicLinearElasticMaterial(n, d);
}

DustMaterial :: ~DustMaterial()
{
    delete LEMaterial;
}

IRResultType
DustMaterial :: initializeFrom(InputRecord *ir)
{
    // Required by IR_GIVE_FIELD macro
    const char *__proc = "initializeFrom";
    IRResultType result;
    // call the corresponding service of structural material
    StructuralMaterial :: initializeFrom(ir);

    // call the corresponding service for the linear elastic material
    this->LEMaterial->initializeFrom(ir);
	 this->bulkModulus = this->LEMaterial->giveBulkModulus();
	 this->shearModulus = this->LEMaterial->giveShearModulus();

    // instanciate the variables defined in DustMaterial
    IR_GIVE_FIELD(ir, alpha, IFT_DustMaterial_alpha, "alpha");
    IR_GIVE_FIELD(ir, beta, IFT_DustMaterial_beta, "beta");
    IR_GIVE_FIELD(ir, lambda, IFT_DustMaterial_lambda, "lambda");
    IR_GIVE_FIELD(ir, theta, IFT_DustMaterial_theta, "theta");
    IR_GIVE_FIELD(ir, ft, IFT_DustMaterial_ft, "ft");
    IR_GIVE_FIELD(ir, rEllipse, IFT_DustMaterial_rEllipse, "rellipse");
	 hardeningType = 0;
	 IR_GIVE_FIELD(ir, hardeningType, IFT_DustMaterial_hardeningType, "ht");
	 if (hardeningType == 0) {
		 IR_GIVE_FIELD(ir, mHard, IFT_DustMaterial_mHard, "mhard");
	 } else if (hardeningType == 1) {
		 // TODO
	 }

    return IRRT_OK;
}

int
DustMaterial :: hasMaterialModeCapability(MaterialMode mMode)
{
    if ( ( mMode == _3dMat ) ||
        ( mMode == _PlaneStrain ) ||
        ( mMode == _3dRotContinuum ) ) {
        return 1;
    } else {
        return 0;
    }
}

void
DustMaterial :: giveRealStressVector(FloatArray &answer,
                                                  MatResponseForm form,
                                                  GaussPoint *gp,
                                                  const FloatArray &totalStrain,
                                                  TimeStep *atTime)
{
    FloatArray strainVectorR;

    DustMaterialStatus *status = ( DustMaterialStatus * ) ( this->giveStatus(gp) );

    // Initialize temp variables for this gauss point
    this->initTempStatus(gp);

    // subtract stress-independent part of strain
    this->giveStressDependentPartOfStrainVector(strainVectorR, gp, totalStrain, atTime, VM_Total);

    // perform the local stress return and update the history variables
    StrainVector strain( strainVectorR, gp->giveMaterialMode() );

    // copy total strain vector to the temp status
    status->letTempStrainVectorBe(totalStrain);

    // give back correct form of stressVector to giveRealStressVector
    if ( form == ReducedForm ) {
        answer = status->giveTempStressVector();
    } else {
        ( ( StructuralCrossSection * ) ( gp->giveElement()->giveCrossSection() ) )
        ->giveFullCharacteristicVector( answer, gp, status->giveTempStressVector() );
    }
}

void
DustMaterial :: give3dMaterialStiffnessMatrix(FloatMatrix &answer,
                                                           MatResponseForm form,
                                                           MatResponseMode mode,
                                                           GaussPoint *gp,
                                                           TimeStep *atTime)
{
    if ( mode == ElasticStiffness) {
        LEMaterial->giveCharacteristicMatrix(answer, form, mode, gp, atTime);
	 } else if (mode == SecantStiffness || mode == TangentStiffness) {
        LEMaterial->giveCharacteristicMatrix(answer, form, mode, gp, atTime);
	 } else {
        _error("\n");
    }
}

int
DustMaterial :: giveIPValue(FloatArray &answer,
                                         GaussPoint *gp,
                                         InternalStateType type,
                                         TimeStep *atTime)
{
    const DustMaterialStatus *status =
        ( DustMaterialStatus * ) giveStatus(gp);
    StrainVector plasticStrainVector(_3dMat);

    switch ( type ) {

        default:
            return StructuralMaterial :: giveIPValue(answer, gp, type, atTime);
        }

    return 0;
}

int
DustMaterial :: giveIPValueSize(InternalStateType type,
                                             GaussPoint *gp)
{
    switch ( type ) {

        default:
            return StructuralMaterial :: giveIPValueSize(type, gp);
        }
}

int
DustMaterial :: giveIntVarCompFullIndx(IntArray &answer,
                                                    InternalStateType type,
                                                    MaterialMode mmode)
{
    switch ( type ) {

        default:
            return StructuralMaterial :: giveIntVarCompFullIndx(answer, type, mmode);
        }
}

InternalStateValueType
DustMaterial :: giveIPValueType(InternalStateType type)
{
    switch ( type ) {

        default:
            return StructuralMaterial :: giveIPValueType(type);
        }
}

MaterialStatus *
DustMaterial :: CreateStatus(GaussPoint *gp) const
{
    DustMaterialStatus *status =
        new  DustMaterialStatus(1, StructuralMaterial :: giveDomain(), gp);
    return status;
}

double
DustMaterial :: functionFe(double i1)
{
	return alpha - lambda*exp(beta*i1) - theta*i1;
}

double
DustMaterial :: functionFc(double sn, double i1, double q)
{
	return sqrt(sn*sn + 1/rEllipse/rEllipse*(q-i1)*(q-i1));
}

double
DustMaterial :: yieldFunction1(double sn, double i1) {
	return sn - functionFe(i1);
}

double
DustMaterial :: yieldFunction2(double sn, double i1, double q)
{
	return functionFc(sn,i1,q) - functionFe(q);
}

double
DustMaterial :: yieldFunction3(double i1)
{
	return i1 - ft;
}

double
DustMaterial :: functionX(double q)
{
	return q - rEllipse*(alpha - lambda*exp(beta*q) - theta*q);
}

} // end namespace oofem

