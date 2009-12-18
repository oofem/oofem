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
#include "stressvector.h"
#include "strainvector.h"

#include "structuralcrosssection.h"
#include "mathfem.h"
#include "contextioerr.h"
#include "datastream.h"

namespace oofem {

// constructor
MisesMat :: MisesMat(int n, Domain *d) : StructuralMaterial(n, d)
{    
    linearElasticMaterial = new IsotropicLinearElasticMaterial(n, d);
    H = 0.;
    sig0 = 0.;
    G = 0.;
    K = 0.;
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

  if ( ( mode == _3dMat ) || ( mode == _3dMat_F ) )
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
    G = ((IsotropicLinearElasticMaterial*) linearElasticMaterial) -> giveShearModulus();
    K = ((IsotropicLinearElasticMaterial*) linearElasticMaterial) -> giveBulkModulus();

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

// returns the stress vector in 3d stress space
void
MisesMat :: giveRealStressVector(FloatArray &answer,
                                          MatResponseForm form,
                                          GaussPoint *gp,
                                          const FloatArray &totalStrain,
                                          TimeStep *atTime)
{
  MaterialMode mode = gp->giveMaterialMode();
  if ( mode == _3dMat )
    giveRealStressVectorComputedFromStrain(answer,form,gp,totalStrain,atTime);
  else if ( mode == _3dMat_F )
    giveRealStressVectorComputedFromDefGrad(answer,form,gp,totalStrain,atTime);
  else 
    OOFEM_ERROR("MisesMat::giveRealStressVector : unknown material response mode");
  return;
}

// returns the stress vector in 3d stress space
// computed from the previous plastic strain and current deformation gradient
void
MisesMat :: giveRealStressVectorComputedFromDefGrad(FloatArray &answer,
                                          MatResponseForm form,
                                          GaussPoint *gp,
                                          const FloatArray &totalDefGrad,
                                          TimeStep *atTime)
{
  FloatArray totalStrain;
  this->convertDefGradToGLStrain(totalDefGrad,totalStrain);
  giveRealStressVectorComputedFromStrain(answer,form,gp,totalStrain,atTime);
  return;
}

// converts the deformation gradient stored by columns in FloatArray F
// into the Green-Lagrange strain with components E11,E22,E33,E23,E31,E12 stored in FloatArray E
void
MisesMat :: convertDefGradToGLStrain(const FloatArray& F, FloatArray& E)
{
  if (F.giveSize() != 9)
    OOFEM_ERROR("MisesMat::convertDefGradToGLStrain : incorrect size of F (must be 9)");
  E.resize(6);
  // note: E(i) is equivalent to E.at(i+1)
  E(0) = F(0)*F(0) + F(1)*F(1) + F(2)*F(2); 
  E(1) = F(3)*F(3) + F(4)*F(4) + F(5)*F(5); 
  E(2) = F(6)*F(6) + F(7)*F(8) + F(8)*F(8); 
  E(3) = F(3)*F(6) + F(4)*F(7) + F(5)*F(8);
  E(4) = F(6)*F(0) + F(7)*F(1) + F(8)*F(2);
  E(5) = F(0)*F(3) + F(1)*F(4) + F(2)*F(5);
  return;
}

// returns the stress vector in 3d stress space
// computed from the previous plastic strain and current total strain
void
MisesMat :: giveRealStressVectorComputedFromStrain(FloatArray &answer,
                                          MatResponseForm form,
                                          GaussPoint *gp,
                                          const FloatArray &totalStrain,
                                          TimeStep *atTime)
{
    double kappa, yieldValue, dKappa;
    FloatArray reducedStress;
    FloatArray strain, plStrain;
    StructuralCrossSection *crossSection = ( StructuralCrossSection * )
                                           ( gp->giveElement()->giveCrossSection() );
    MisesMatStatus *status = ( MisesMatStatus * ) this->giveStatus(gp);

    this->initTempStatus(gp);
    this->initGpForNewStep(gp);

    // subtract the stress-independent part of strain (e.g. thermal strain)
    this->giveStressDependentPartOfStrainVector(strain, gp, totalStrain,
                                                atTime, VM_Total);

    // get the initial plastic strain and initial kappa from the status
    status->givePlasticStrain(plStrain);
    kappa = status->giveCumulativePlasticStrain();

    // === radial return algorithm ===
    
    // elastic predictor
    StrainVector elStrain(totalStrain,_3dMat);
    elStrain.substract(plStrain);
    StrainVector elStrainDev(_3dMat);
    double elStrainVol;
    elStrain.computeDeviatoricVolumetricSplit(elStrainDev,elStrainVol);
    StressVector trialStressDev(_3dMat);
    elStrainDev.applyDeviatoricElasticStiffness(trialStressDev, G);
    // store the deviatoric trial stress (reused by algorithmic stiffness) 
    status->letTrialStressDevBe(trialStressDev);

    // check the yield condition at the trial state
    double trialS = trialStressDev.computeStressNorm(); 
    double sigmaY = sig0 + H*kappa;
    yieldValue = sqrt(3./2.)*trialS - sigmaY;
    if ( yieldValue > 0. ){
      // increment of cumulative plastic strain
      dKappa = yieldValue / (H+3.*G);
      kappa += dKappa;
      StrainVector dPlStrain(_3dMat);
      // the following line is equivalent to multiplication by scaling matrix P
      trialStressDev.applyDeviatoricElasticCompliance(dPlStrain,0.5);
      // increment of plastic strain
      dPlStrain.times(sqrt(3./2.)*dKappa/trialS);
      plStrain.add(dPlStrain);
      // scaling of deviatoric trial stress 
      trialStressDev.times(1.-sqrt(6.)*G*dKappa/trialS);
    }

    // assemble the stress from the elastically computed volumetric part
    // and scaled deviatoric part
    StressVector fullStress(_3dMat);
    double stressVol = 3.*K*elStrainVol;
    trialStressDev.computeDeviatoricVolumetricSum(fullStress, stressVol);

    // store the total strain in status
    status->letTempStrainVectorBe(totalStrain);

    // reduce the stress vector and store it in status
    crossSection->giveReducedCharacteristicVector(reducedStress, gp, fullStress);
    status->letTempStressVectorBe(reducedStress);

    // store the plastic strain and cumulative plastic strain
    status->letTempPlasticStrainBe(plStrain);
    status->setTempCumulativePlasticStrain(kappa);

    // reduce the stress vector (if required) and return it
    if ( form == FullForm ) 
      answer = fullStress;
    else                                                             
      answer = reducedStress;
    return;
}

// returns the consistent (algorithmic) tangent stiffness matrix
void
MisesMat :: give3dMaterialStiffnessMatrix(FloatMatrix &answer, MatResponseForm form,
                                                   MatResponseMode mode,
                                                   GaussPoint *gp,
                                                   TimeStep *atTime)
{
  // start from the elastic stiffness
  this->giveLinearElasticMaterial()->giveCharacteristicMatrix(answer, form, mode, gp, atTime);
  if ( mode != TangentStiffness ) 
    return;

  MisesMatStatus *status = ( MisesMatStatus * ) this->giveStatus(gp);
  //StructuralCrossSection *crossSection = ( StructuralCrossSection * ) ( gp->giveElement()->giveCrossSection() );
  double kappa = status->giveCumulativePlasticStrain();
  // increment of cumulative plastic strain as an indicator of plastic loading
  double dKappa = status->giveTempCumulativePlasticStrain() - kappa;

  if (dKappa <= 0.0) // elastic loading - elastic stiffness plays the role of tangent stiffness
    return;

  // === plastic loading ===

  // yield stress at the beginning of the step
  double sigmaY = sig0 + H*kappa;

  // trial deviatoric stress and its norm
  StressVector trialStressDev(_3dMat);
  status -> giveTrialStressDev(trialStressDev);
  double trialS = trialStressDev.computeStressNorm(); 

  // one correction term
  FloatMatrix stiffnessCorrection(6,6);
  stiffnessCorrection.beDyadicProductOf(trialStressDev,trialStressDev);
  double factor = -2.*sqrt(6.)*G*G/trialS;
  double factor1 = factor*sigmaY/((H+3.*G)*trialS*trialS);
  stiffnessCorrection.times(factor1);
  answer.plus(stiffnessCorrection);

  // another correction term
  stiffnessCorrection.bePinvID();
  double factor2 = factor*dKappa;
  stiffnessCorrection.times(factor2);
  answer.plus(stiffnessCorrection);

  return;
}


//=============================================================================

MisesMatStatus :: MisesMatStatus(int n, Domain *d, GaussPoint *g) :
  StructuralMaterialStatus(n, d, g), plasticStrain(), tempPlasticStrain(), trialStress()
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
    /*
    // this would not be correct, since the status is already updated and kappa is now the "final" value
    if ( tempKappa > kappa ) {
      fprintf(file, " Yielding, ");
    } else {
      fprintf(file, " Unloading, ");
    }
    */

    // print the plastic strain
        n = plasticStrain.giveSize();
        fprintf(file, " plastic strains ");
        for ( i = 1; i <= n; i++ ) {
            fprintf( file, " % .4e", plasticStrain.at(i) );
        }

	// print the cumulative plastic strain
	fprintf(file, ", kappa ");
	fprintf( file, " % .4e", kappa );

	fprintf(file, "}\n");
}


// initializes temporary variables based on their values at the previous equlibrium state
void MisesMatStatus :: initTempStatus()
{
    StructuralMaterialStatus :: initTempStatus();

    if ( plasticStrain.giveSize() == 0 ) {
        plasticStrain.resize( ( ( StructuralMaterial * ) gp->giveMaterial() )->
                                   giveSizeOfReducedStressStrainVector( gp->giveMaterialMode() ) );
        plasticStrain.zero();
    }

    tempPlasticStrain = plasticStrain;
    tempKappa = kappa;
    trialStress.resize(0); // to indicate that it is not defined yet
}


// updates internal variables when equilibrium is reached
void
MisesMatStatus :: updateYourself(TimeStep *atTime)
{
    StructuralMaterialStatus :: updateYourself(atTime);

    plasticStrain = tempPlasticStrain;
    kappa = tempKappa;
    trialStress.resize(0); // to indicate that it is not defined any more
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

    // write plastic strain (vector)
    if ( ( iores = plasticStrain.storeYourself(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    // write cumulative plastic strain (scalar)
    if (!stream->write(&kappa,1)) THROW_CIOERR(CIO_IOERR);

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

    // read plastic strain (vector)
    if ( ( iores = plasticStrain.restoreYourself(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    // read cumulative plastic strain (scalar)
    if (!stream->read (&kappa,1)) THROW_CIOERR(CIO_IOERR);

    return CIO_OK; // return succes
}

} // end namespace oofem
