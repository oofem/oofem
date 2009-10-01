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
 *               Copyright (C) 1993 - 2008   Borek Patzak
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
// file : misesmat.C


#include "misesmat.h"
#include "isolinearelasticmaterial.h"
#include "gausspnt.h"
#include "flotmtrx.h"
#include "flotarry.h"
#include "intarray.h"

#include "structuralcrosssection.h"
#include "mathfem.h"
//#include "datastream.h"
#include "contextioerr.h"

// constructor
MisesMat :: MisesMat(int n, Domain *d) : StructuralMaterial(n, d)
{    
    H = 0.;
    sig0 = 0.;
    linearElasticMaterial = new IsotropicLinearElasticMaterial(n, d);
}

// destructor
MisesMat :: ~MisesMat()
{ }

// specifies whether a given material mode is supported by this model
int 
MisesMat :: hasMaterialModeCapability (MaterialMode mode)
{
  /*
    if ( ( mode == _3dMat ) ||
        ( mode == _1dMat ) ||
        //<RESTRICTED_SECTION>
        ( mode == _PlaneStress )  ||
        //</RESTRICTED_SECTION>
        ( mode == _PlaneStrain )  ||
        ( mode == _2dPlateLayer ) ||
        ( mode == _2dBeamLayer )  ||
        ( mode == _3dShellLayer ) ||
        ( mode == _1dFiber ) ) {
        return 1;
    }
  */

  if ( ( mode == _3dMat ) )
    return 1;
  return 0;
}

// reads the model parameters from the input file
IRResultType
MisesMat :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // required by IR_GIVE_FIELD macro
    IRResultType result;                 // required by IR_GIVE_FIELD macro

    StructuralMaterial :: initializeFrom(ir);
    linearElasticMaterial->initializeFrom(ir); // takes care of elastic constants

    IR_GIVE_FIELD(ir, sig0, IFT_MisesMat_sig0, "sig0"); // uniaxial yield stress

    H = 0.;
    IR_GIVE_OPTIONAL_FIELD(ir, H, IFT_MisesMat_h, "h"); // hardening modulus

    return IRRT_OK;
}

// creates a new material status  corresponding to this class
MaterialStatus *
MisesMat :: CreateStatus(GaussPoint *gp) const
{
    MisesMatStatus *status;
    status = new MisesMatStatus(1, this->giveDomain(), gp);
    return status;
}

// computes the value of the yield function
double
MisesMat :: computeYieldValueAt(GaussPoint *gp, const FloatArray &stressVector, double kap)
{
    double f;

    // compute the square root of 3 times J2
    if ( stressVector.isNotEmpty() ) 
      f = sqrt(3. * this->computeJ2InvariantAt(stressVector) );
    else
      f = 0.;

    // subtract the current yield stress (given by the linear hardening rule)
    f -= (sig0 + H*kap);
    return f;
}

// computes the J2 invariant
double
MisesMat :: computeJ2InvariantAt(const FloatArray &stressVector)
{
    double answer;
    double v1, v2, v3;

    if ( stressVector.isEmpty() ) {
        return 0.0;
    }

    v1 = ( ( stressVector.at(1) - stressVector.at(2) ) * ( stressVector.at(1) - stressVector.at(2) ) );
    v2 = ( ( stressVector.at(2) - stressVector.at(3) ) * ( stressVector.at(2) - stressVector.at(3) ) );
    v3 = ( ( stressVector.at(3) - stressVector.at(1) ) * ( stressVector.at(3) - stressVector.at(1) ) );

    answer = ( 1. / 6. ) * ( v1 + v2 + v3 ) + stressVector.at(4) * stressVector.at(4) +
    stressVector.at(5) * stressVector.at(5) + stressVector.at(6) * stressVector.at(6);

    return answer;
}


// returns the stress vector in 3d stress space
// computed from the previous stress and current strain
void
MisesMat :: giveRealStressVector(FloatArray &answer,
                                          MatResponseForm form,
                                          GaussPoint *gp,
                                          const FloatArray &totalStrain,
                                          TimeStep *atTime)
{
    double kappa, yieldValue;
    FloatArray fullStressVector, reducedStressVector;
    FloatArray strainVectorR, elasticStrainVectorR, plasticStrainVectorR;
    FloatMatrix elasticStiffness, elasticCompliance;
    StructuralCrossSection *crossSection = ( StructuralCrossSection * )
                                           ( gp->giveElement()->giveCrossSection() );
    MisesMatStatus *status = ( MisesMatStatus * ) this->giveStatus(gp);

    this->initTempStatus(gp);
    this->initGpForNewStep(gp);

    // substract stress-independent part of strain (e.g. thermal strain)
    this->giveStressDependentPartOfStrainVector(strainVectorR, gp, totalStrain,
                                                atTime, VM_Total);

    // prepare the elastic stiffness and its inverse

    this->giveLinearElasticMaterial()->giveCharacteristicMatrix(elasticStiffness, ReducedForm, ElasticStiffness, gp, atTime);
    elasticCompliance.beInverseOf(elasticStiffness);

    // === radial return algorithm ===

    // get the initial plastic strain and initial kappa from the status
    status->givePlasticStrainVector(plasticStrainVectorR);
    status->giveCumulativePlasticStrain(kappa);
    
    // elastic predictor
    elasticStrainVectorR = totalStrain;
    elasticStrainVectorR.substract(plasticStrainVectorR);

    reducedStressVector.beProductOf(elasticStiffness, elasticStrainVectorR);
    crossSection->giveFullCharacteristicVector(fullStressVector, gp, reducedStressVector);

    // check the yield condition at the trial state
    yieldValue = this->computeYieldValueAt(gp, fullStressVector, kappa);
    if ( yieldValue > 0. ){

      // plastic return to the yield surface

    }

    // store the total strain in status
    status->letTempStrainVectorBe(totalStrain);

    // reduce the stress vector and store it in status
    crossSection->giveReducedCharacteristicVector(reducedStressVector, gp, fullStressVector);
    status->letTempStressVectorBe(reducedStressVector);

    // store the plastic strain and cumulative plastic strain
    status->letTempPlasticStrainVectorBe(plasticStrainVectorR);
    status->setTempCumulativePlasticStrain(kappa);

    // reduce the stress vector (if required) and return it
    if ( form == FullForm ) {
        answer = fullStressVector;
        return;
    } else                                                             {
        crossSection->giveReducedCharacteristicVector(answer, gp, fullStressVector);
        return;
    }
}

// returns the consistent (algorithmic) tangent stiffness matrix
void
MisesMat :: give3dMaterialStiffnessMatrix(FloatMatrix &answer, MatResponseForm form,
                                                   MatResponseMode mode,
                                                   GaussPoint *gp,
                                                   TimeStep *atTime)
{
    if ( mode == TangentStiffness ) {
        this->giveLinearElasticMaterial()->giveCharacteristicMatrix(answer, form, mode, gp, atTime);
    } else {
        this->giveLinearElasticMaterial()->giveCharacteristicMatrix(answer, form, mode, gp, atTime);
    }

    return;
}


//=============================================================================

MisesMatStatus :: MisesMatStatus(int n, Domain *d, GaussPoint *g) :
    StructuralMaterialStatus(n, d, g), plasticStrainVector(), tempPlasticStrainVector()
{
    kappa = tempKappa = 0.;
}

MisesMatStatus :: ~MisesMatStatus()
{ }

void
MisesMatStatus :: printOutputAt(FILE *file, TimeStep *tStep)
{
  int i, n;

    StructuralMaterialStatus :: printOutputAt(file, tStep);
    fprintf(file, "status { ");
    if ( tempKappa > kappa ) {
      fprintf(file, " Yielding, ");
    } else {
      fprintf(file, " Unloading, ");
    }

        n = plasticStrainVector.giveSize();
        fprintf(file, " plastic strains ");
        for ( i = 1; i <= n; i++ ) {
            fprintf( file, " % .4e", plasticStrainVector.at(i) );
        }

	fprintf(file, ", kappa ");
	fprintf( file, " % .4e", kappa );

	fprintf(file, "}\n");
}


// initializes temporary variables based on their values at the previous equlibrium state
void MisesMatStatus :: initTempStatus()
{
    StructuralMaterialStatus :: initTempStatus();

    if ( plasticStrainVector.giveSize() == 0 ) {
        plasticStrainVector.resize( ( ( StructuralMaterial * ) gp->giveMaterial() )->
                                   giveSizeOfReducedStressStrainVector( gp->giveMaterialMode() ) );
        plasticStrainVector.zero();
    }

    tempPlasticStrainVector = plasticStrainVector;
    tempKappa = kappa;
}


// updates internal variables when equilibrium is reached
void
MisesMatStatus :: updateYourself(TimeStep *atTime)
{
    StructuralMaterialStatus :: updateYourself(atTime);

    plasticStrainVector = tempPlasticStrainVector;
    kappa = tempKappa;
}


// saves full information stored in this status
// temporary variables are NOT stored
contextIOResultType
MisesMatStatus :: saveContext(DataStream *stream, ContextMode mode, void *obj)
{
    contextIOResultType iores;

    // save parent class status
    if ( ( iores = StructuralMaterialStatus :: saveContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    // write raw data
    if ( ( iores = plasticStrainVector.storeYourself(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }


    // !!!!!!!!!!!!!!!! kappa needs to be written !!!!!!!!!!!!!!!!!

    return CIO_OK;
}



contextIOResultType
MisesMatStatus :: restoreContext(DataStream *stream, ContextMode mode, void *obj)
//
// restores full information stored in stream to this Status
//
{
    contextIOResultType iores;

    // read parent class status
    if ( ( iores = StructuralMaterialStatus :: restoreContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = plasticStrainVector.restoreYourself(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    // !!!!!!!!!!!!!!!! kappa needs to be restored !!!!!!!!!!!!!!!!!

    return CIO_OK; // return succes
}

