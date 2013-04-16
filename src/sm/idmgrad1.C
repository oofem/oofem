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

#include "idmgrad1.h"
#include "gausspnt.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "structuralcrosssection.h"
#include "mathfem.h"
#include "sparsemtrx.h"
#include "isolinearelasticmaterial.h"
#include "dynalist.h"
#include "error.h"
#include "nonlocalmaterialext.h"
#include "datastream.h"
#include "contextioerr.h"
#include "stressvector.h"
#include "strainvector.h"

#ifdef __PARALLEL_MODE
 #include "combuff.h"
#endif

namespace oofem {

IDGMaterial :: IDGMaterial(int n, Domain *d):IsotropicDamageMaterial1(n,d) 
    //
    // constructor
    //
{
  internalLength = 0;
  averType = 0;
}


IDGMaterial :: ~IDGMaterial()
//
// destructor
//
{ }

IRResultType
IDGMaterial :: initializeFrom(InputRecord *ir)
{
    //const char *__proc = "initializeFrom";     // Required by IR_GIVE_FIELD macro
    //IRResultType result;                // Required by IR_GIVE_FIELD macro
  
    IsotropicDamageMaterial1 :: initializeFrom(ir);

    //internal Length parameter    
    //IR_GIVE_OPTIONAL_FIELD(ir, internalLength, _IFT_IDGMaterial_internalLength, "internallength");
    averType = 1;
    // averaging type approach 0-standard else stress based averaging
    //IR_GIVE_OPTIONAL_FIELD(ir, averType, _IFT_IDGMaterial_averType, "avertype");
   
    // parameter for avetType = 1
    if(averType == 1) {
        beta = 0.5;
        //IR_GIVE_FIELD(ir, beta, _IFT_IDGMaterial_beta);
        t = 1;
        //IR_GIVE_FIELD(ir, t, _IFT_IDGMaterial_t, "t");
    }

    this->mapper.initializeFrom(ir);

    return IRRT_OK;
}

/////////////////////////////////////////////////////////////////////////////



int
IDGMaterial :: hasMaterialModeCapability(MaterialMode mode)
{
    if ( mode == _1dMatGrad || mode == _PlaneStressGrad ) {
        return 1;
    }

    return 0;
}

void
IDGMaterial :: giveCharacteristicMatrix(FloatMatrix &answer,
                                         MatResponseForm form, MatResponseMode rMode, GaussPoint *gp, TimeStep *atTime)
//
// Returns characteristic material stiffness matrix of the receiver
//
{
    MaterialMode mMode = gp->giveMaterialMode();
    switch ( mMode ) {
    case _1dMatGrad:
        if ( form == PDGrad_uu ) {
            give1dStressStiffMtrx(answer, form, rMode, gp, atTime);
        } else if ( form == PDGrad_ku ) {
            give1dKappaMatrix(answer, form, rMode, gp, atTime);
        } else if ( form == PDGrad_uk ) {
            give1dGprime(answer, form, rMode, gp, atTime);
        } else if ( form == PDGrad_kk ) {
            giveInternalLength(answer, form, rMode, gp, atTime);
        }
    case _PlaneStressGrad:
        if ( form == PDGrad_uu ) {
            givePlaneStressStiffMtrx(answer, form, rMode, gp, atTime);
        } else if ( form == PDGrad_ku ) {
            givePlaneStressKappaMatrix(answer, form, rMode, gp, atTime);
        } else if ( form == PDGrad_uk ) {
            givePlaneStressGprime(answer, form, rMode, gp, atTime);
        } else if ( form == PDGrad_kk ) {
            giveInternalLength(answer, form, rMode, gp, atTime);
        } else if ( form == PDGrad_LD ) {
            giveInternalLengthDerivative(answer, form, rMode, gp, atTime);
        }
        break;

    default:
        _error2( "giveCharacteristicMatrix : unknown mode (%s)", __MaterialModeToString(mMode) );
    }
}





/////////////////////////////////////////////////////////////////
// BEGIN: EVALUATION OF LOCAL STIFFNESS MATRIX

//compute derivative of the equivalent strain wrt strain
void
IDGMaterial :: computeEta(FloatMatrix &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *atTime)
{
    LinearElasticMaterial *lmat = this->giveLinearElasticMaterial();
    answer.resize(1,3);
    if ( strain.isEmpty() ) {
        answer.zero();
        return;
    }
 
    if ( this->equivStrainType == EST_Mazars ) {
        double posNorm = 0.0;
        FloatArray principalStrains, fullstrain;
        
        // if plane stress mode -> compute strain in z-direction from condition of zero stress in corresponding direction
        if ( gp->giveMaterialMode() == _PlaneStress || gp->giveMaterialMode() == _PlaneStressGrad  ) {
            double nu = lmat->give(NYxz, gp);
            FloatMatrix N,m,Eta(2,2);
            Eta.zero();
            FloatArray n(2);
            StrainVector fullStrain(strain,_PlaneStress);
            fullStrain.computePrincipalValDir(principalStrains, N);
            principalStrains.resize(3);
            principalStrains.at(3) = -nu * ( principalStrains.at(1) + principalStrains.at(2) ) / ( 1. - nu );

            for ( int i = 1; i <= 3; i++ ) {
                if ( i < 3 ) {
                    if ( principalStrains.at(i) > 0.0 ) {
                        double e = principalStrains.at(i);
                        for(int j = 1; j<3;j++)
                            n.at(j) = N.at(i,j);
                        m.beDyadicProductOf(n,n);
                        m.times(e);
                        Eta.add(m);
                    }
                }
                if ( principalStrains.at(i) > 0.0 ) 
                    posNorm += principalStrains.at(i) * principalStrains.at(i);
            }
            double kappa = sqrt(posNorm);
            Eta.times(1./kappa);
            answer.at(1,1) = Eta.at(1,1);
            answer.at(1,2) = Eta.at(2,2);
            answer.at(1,3) = Eta.at(1,2);
        }
    } else {
        _error("computeEta: unknown EquivStrainType");
    }
  
}


void
IDGMaterial ::  give1dStressStiffMtrx(FloatMatrix & answer,  MatResponseForm form, MatResponseMode mode, GaussPoint * gp,  TimeStep * tStep)
{
    IsotropicDamageMaterialStatus *status = static_cast< IsotropicDamageMaterialStatus * >( this->giveStatus(gp) );
    LinearElasticMaterial *lmat = this->giveLinearElasticMaterial();
    double om;
    om = status->giveTempDamage();
    om = min(om, maxOmega);
    answer.resize(1,1);
    answer.at(1,1) = lmat->give('E', gp);    
    answer.times(1.0 - om);

}

void
IDGMaterial :: give1dKappaMatrix(FloatMatrix &answer, MatResponseForm form, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
  IDGMaterialStatus *status = static_cast< IDGMaterialStatus * >( this->giveStatus(gp) );
  LinearElasticMaterial *lmat = this->giveLinearElasticMaterial();
  answer.resize(1,1);
  if(status->giveTempStrainVector().at(1) > 0)
    answer.at(1,1) = 1;
  else if(status->giveTempStrainVector().at(1) < 0)
    answer.at(1,1) = 2 * lmat -> give('n',gp) * lmat -> give('n',gp);
  else
    answer.at(1,1) = 0;

  
}

void
IDGMaterial :: give1dGprime(FloatMatrix &answer, MatResponseForm form, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
   IDGMaterialStatus *status = static_cast< IDGMaterialStatus * >( this->giveStatus(gp) );
   LinearElasticMaterial *lmat = this->giveLinearElasticMaterial();

   double damage = status ->giveDamage();
   double tempDamage = status ->giveTempDamage();
   double E = lmat->give('E', gp);

   answer.resize(1, 1);
   if ( ( tempDamage - damage ) > 0 ) {
     double nlKappa =  status->giveTempStrainVector().at(2);
     answer.at(1, 1) = E * status->giveTempStrainVector().at(1);
     double gPrime = damageFunctionPrime(nlKappa,gp);
     answer.times(gPrime);
    } else {
        answer.zero();
    }


}



void
IDGMaterial :: givePlaneStressStiffMtrx(FloatMatrix &answer, MatResponseForm form, MatResponseMode mode, GaussPoint *gp, TimeStep *atTime)
{
    IDGMaterialStatus *status = static_cast< IDGMaterialStatus * >( this->giveStatus(gp) );
    double tempDamage;
    if ( mode == ElasticStiffness ) {
        tempDamage = 0.0;
    } else {
        tempDamage = status->giveTempDamage();
    if ( tempDamage > 0.0 )
        tempDamage = min(tempDamage, maxOmega);
    }
    this->giveLinearElasticMaterial()->giveCharacteristicMatrix(answer, form, mode, gp, atTime);
    answer.times(1.0 - tempDamage);
    
#if 0
     double damage = status->giveDamage();
    if(tempDamage>damage) {
        double kappa = status->giveKappa();
        FloatArray strain;
        FloatArray stress;
        FloatArray eta;
        FloatMatrix correctionTerm;
        stress = status->giveTempStressVector();
        strain = status->giveTempStrainVector();
        double nlKappa = strain.at(4);
        strain.resize(3);
        stress.times(1./(1-tempDamage));
        stress.resize(3);
        this->computeEta(eta,strain,gp,atTime);
        double dDamage = damageFunctionPrime(nlKappa,gp);
        correctionTerm.beDyadicProductOf(stress,eta);
        correctionTerm.times(-dDamage);
        answer.add(correctionTerm);
      }
#endif
}


void
IDGMaterial :: givePlaneStressKappaMatrix(FloatMatrix &answer, MatResponseForm form, MatResponseMode mode, GaussPoint *gp, TimeStep *atTime)
{
    // only for Mazars equivalent deformation ...answer = <eps>/eps_eq 
    // only plane-stress case
    IDGMaterialStatus *status = static_cast< IDGMaterialStatus * >( this->giveStatus(gp) );
    //    double kappa = status->giveKappa();
    //    double nlKappa = totalStrain.at(4);
    FloatArray  totalStrain =  status->giveTempStrainVector();  
    double kappa = status->giveKappa();
    double tempKappa = status->giveTempKappa();
    totalStrain.resize(3);
    answer.resize(1,3);

    if ( tempKappa > kappa) {
        this->computeEta(answer,totalStrain,gp,atTime);
    } else 
        answer.zero();
  
}
  
void
IDGMaterial :: givePlaneStressGprime(FloatMatrix &answer, MatResponseForm form, MatResponseMode mode, GaussPoint *gp, TimeStep *atTime)
{

    IDGMaterialStatus *status = static_cast< IDGMaterialStatus * >( this->giveStatus(gp) );
    double damage = status->giveDamage();
    double tempDamage = status->giveTempDamage();
    answer.resize(3,1);
    if ( tempDamage > damage ) {
        double nlEquivStrain =  status->giveTempStrainVector().at(4);
        double gPrime =  this -> damageFunctionPrime(nlEquivStrain, gp);
        FloatArray stress =  status->giveTempStressVector();
        answer.at(1, 1) = stress.at(1)/(1-tempDamage);
        answer.at(2, 1) = stress.at(2)/(1-tempDamage);
        answer.at(3, 1) = stress.at(3)/(1-tempDamage);
        answer.times(gPrime);
    }
}


void
IDGMaterial :: giveInternalLength(FloatMatrix &answer, MatResponseForm form, MatResponseMode rMode, GaussPoint *gp, TimeStep *atTime)
{
    if ( averType == 0 ) {
        answer.resize(1, 1);
        answer.at(1, 1) = internalLength*internalLength;
    } else if ( averType == 1 ) {
        double distance = 0;
        //double distance = gp ->getDistanceToBoundary();
        answer.resize(1, 1);
        if ( distance < t*internalLength )
            answer.at(1,1) = (internalLength*beta + (1-beta)*distance/t)*(internalLength*beta + (1-beta)*distance/t);
        else
            answer.at(1, 1) = internalLength*internalLength;
    } else if ( averType == 2 ) {
        MaterialMode mMode = gp->giveMaterialMode();
        if ( mMode == _PlaneStressGrad ) {
            answer.resize(2,2);
            IDGMaterialStatus *status = static_cast< IDGMaterialStatus * >( this->giveStatus(gp) );
            StressVector finalStress(_PlaneStressGrad);
            FloatArray stress = status->giveTempStressVector();
            (FloatArray)finalStress = status->giveTempStressVector();
            FloatArray sigPrinc;
            FloatMatrix nPrinc;
            // get principal trial stresses (ordered) and principal stress directions
            finalStress.computePrincipalValDir(sigPrinc, nPrinc);
            if ( nPrinc.giveNumberOfRows() == 0 ) {
                nPrinc.resize(2,2);
                nPrinc.at(1,1) = 1;
                nPrinc.at(2,2) = 1;
            }
            // principal internal lengths
            double l1 = internalLength;
            double l2 = internalLength;
            double denominator = beta*sigPrinc.at(1)*sigPrinc.at(1)+(1-beta)*sigPrinc.at(2)*sigPrinc.at(2);
            if(denominator != 0)
                l2 *= sigPrinc.at(1)*sigPrinc.at(1)/denominator;
            // compose the internal Length matrix in global coordinates
            //   the first subscript refers to coordinate
            //   the second subscript refers to eigenvalue
            double n11 = nPrinc.at(1, 1);
            double n12 = nPrinc.at(1, 2);
            double n21 = nPrinc.at(2, 1);
            double n22 = nPrinc.at(2, 2);
            answer.at(1,1) = l1*n11 * n11 + l2 * n12 * n12;
            answer.at(1,2) = l1 * n11 * n21 + l2 * n12 * n22;
            answer.at(2,1) = l1 * n11 * n21 + l2 * n12 * n22;
            answer.at(2,2) = l1 * n21 * n21 + l2 * n22 * n22;
        }
    }
}


void
IDGMaterial :: giveInternalLengthDerivative(FloatMatrix &answer, MatResponseForm form, MatResponseMode mode, GaussPoint *gp, TimeStep *atTime)
{
  
    if ( averType == 1 ) {
        MaterialMode mMode = gp->giveMaterialMode();
        if ( mMode == _PlaneStressGrad ) {
            answer.resize(4,4);
            IDGMaterialStatus *status = static_cast< IDGMaterialStatus * >( this->giveStatus(gp) );
            StressVector finalStress(_PlaneStressGrad);
            (FloatArray)finalStress = status->giveTempStressVector();
            FloatArray sigPrinc;
            FloatMatrix nPrinc;
            // get principal trial stresses (ordered) and principal stress directions
            finalStress.computePrincipalValDir(sigPrinc, nPrinc);
            if(nPrinc.giveNumberOfRows() == 0) {
                nPrinc.resize(2,2);
                nPrinc.zero();
            }
        
            // principal internal lengths
            double denominator = (beta*sigPrinc.at(1)*sigPrinc.at(1)+(1-beta)*sigPrinc.at(2)*sigPrinc.at(2));
            double derivativeSig1 = 0;
            double derivativeSig2 = 0;
            if(denominator != 0) {
                derivativeSig1 = 2*internalLength*(sigPrinc.at(1)*denominator - beta*sigPrinc.at(1)*sigPrinc.at(1)*sigPrinc.at(1))/denominator/denominator;
                derivativeSig2 = 2*internalLength* (beta-1)*sigPrinc.at(1)*sigPrinc.at(1)*sigPrinc.at(2)/denominator/denominator;
            }
            // compose the internal Length matrix in global coordinates
            //   the first subscript refers to coordinate
            //   the second subscript refers to eigenvalue
            answer.at(1,1) = nPrinc.at(1, 1);
            answer.at(1,2) = nPrinc.at(1, 2);
            answer.at(2,1) = nPrinc.at(2, 1);
            answer.at(2,2) = nPrinc.at(2, 2);
            answer.at(3,3) = derivativeSig1;
            answer.at(4,4) = derivativeSig2;
        }
    }
}



////////////////////////////////////////////////////////////////////////////
void
IDGMaterial :: computeEquivalentStrain(double &kappa, const FloatArray &strain, GaussPoint *gp, TimeStep *atTime)
{
    LinearElasticMaterial *lmat = this->giveLinearElasticMaterial();
    StructuralCrossSection *crossSection = static_cast< StructuralCrossSection * >( gp->giveElement()->giveCrossSection() );

    if ( strain.isEmpty() ) {
        kappa = 0.;
        return;
    }

    if ( this->equivStrainType == EST_Mazars ) {
        double posNorm = 0.0;
        FloatArray principalStrains, fullstrain;

        crossSection->giveFullCharacteristicVector(fullstrain, gp, strain);
        fullstrain.resize(6);
        // if plane stress mode -> compute strain in z-direction from condition of zero stress in corresponding direction
        if ( gp->giveMaterialMode() == _PlaneStressGrad ) {
            double nu = lmat->give(NYxz, gp);
            fullstrain.at(3) = -nu * ( fullstrain.at(1) + fullstrain.at(2) ) / ( 1. - nu );
        } else if ( gp->giveMaterialMode() == _1dMatGrad ) {
            double nu = lmat->give(NYxz, gp);
            fullstrain.at(2) = -nu *fullstrain.at(1);
            fullstrain.at(3) = -nu *fullstrain.at(1);
        }

        this->computePrincipalValues(principalStrains, fullstrain, principal_strain);

        for ( int i = 1; i <= 3; i++ ) {
            if ( principalStrains.at(i) > 0.0 ) {
                posNorm += principalStrains.at(i) * principalStrains.at(i);
            }
        }

        kappa = sqrt(posNorm);
    } else if ( ( this->equivStrainType == EST_Rankine_Smooth ) || ( this->equivStrainType == EST_Rankine_Standard ) ) {
        // EST_Rankine equiv strain measure
        FloatMatrix de;
        FloatArray stress, fullStress, principalStress;
        double sum = 0.;

        lmat->giveCharacteristicMatrix(de, ReducedForm, SecantStiffness, gp, atTime);
        stress.beProductOf(de, strain);
        crossSection->giveFullCharacteristicVector(fullStress, gp, stress);
        this->computePrincipalValues(principalStress, fullStress, principal_stress);
        for ( int i = 1; i <= 3; i++ ) {
            if ( principalStress.at(i) > 0.0 ) {
                if ( this->equivStrainType == EST_Rankine_Smooth ) {
                    sum += principalStress.at(i) * principalStress.at(i);
                } else if ( sum < principalStress.at(i) ) {
                    sum = principalStress.at(i);
                }
            } else if ( sum < principalStress.at(i) ) {
                sum = principalStress.at(i);
            }
        }

        if ( this->equivStrainType == EST_Rankine_Smooth ) {
            sum = sqrt(sum);
        }

        kappa = sum / lmat->give('E', gp);
    } else if ( ( this->equivStrainType == EST_ElasticEnergy ) || ( this->equivStrainType == EST_ElasticEnergyPositiveStress ) || ( this->equivStrainType == EST_ElasticEnergyPositiveStrain ) ) {
        // equivalent strain expressions based on elastic energy
        FloatMatrix de;
        FloatArray stress;
        double sum;

        lmat->giveCharacteristicMatrix(de, ReducedForm, SecantStiffness, gp, atTime);
        if ( this->equivStrainType == EST_ElasticEnergy ) {
            // standard elastic energy
            stress.beProductOf(de, strain);
            sum = strain.dotProduct(stress);
        } else if ( this->equivStrainType == EST_ElasticEnergyPositiveStress ) {
            // elastic energy corresponding to positive part of stress
            FloatArray fullStress, principalStress;
            crossSection->giveFullCharacteristicVector(fullStress, gp, stress);
            this->computePrincipalValues(principalStress, fullStress, principal_stress);
            // TO BE FINISHED
            sum = 0.;
            OOFEM_ERROR("Elastic energy corresponding to positive part of stress not finished\n");
        } else {
            // elastic energy corresponding to positive part of strain
            // TO BE DONE
            sum = 0.;
            OOFEM_ERROR("Elastic energy corresponding to positive part of strain not finished\n");
        }

        kappa = sqrt( sum / lmat->give('E', gp) );
    } else if ( this->equivStrainType == EST_Mises ) {
        double nu = lmat->give(NYxz, NULL);
        FloatArray principalStrains, fullstrain;
        crossSection->giveFullCharacteristicVector(fullstrain, gp, strain);
        if ( gp->giveMaterialMode() == _PlaneStress ) {
            fullstrain.at(3) = -nu * ( fullstrain.at(1) + fullstrain.at(2) ) / ( 1. - nu );
        } else if ( gp->giveMaterialMode() == _1dMat ) {
            fullstrain.at(2) = -nu *fullstrain.at(1);
            fullstrain.at(3) = -nu *fullstrain.at(1);
        }

        this->computePrincipalValues(principalStrains, fullstrain, principal_strain);
        double I1e, J2e;
        this->computeStrainInvariants(principalStrains, I1e, J2e);
        double a, b, c;
        a = ( k - 1 ) * I1e / ( 2 * k * ( 1 - 2 * nu ) );
        b = ( k - 1 ) * ( k - 1 ) * I1e * I1e / ( ( 1 - 2 * nu ) * ( 1 - 2 * nu ) );
        c = 12 * k * J2e / ( ( 1 + nu ) * ( 1 + nu ) );
        kappa = a + 1 / ( 2 * k ) * sqrt(b + c);
    } else {
        _error("computeEquivalentStrain: unknown EquivStrainType");
    }
}


void
IDGMaterial :: initDamaged(double kappa, FloatArray &strainVector, GaussPoint *gp)
{
    int indx = 1;
    double le=0.;
    double E = this->giveLinearElasticMaterial()->give('E', gp);
    FloatArray principalStrains, crackPlaneNormal(3), fullstrain, crackVect(3);
    FloatMatrix principalDir(3, 3);
    IsotropicDamageMaterial1Status *status = static_cast< IsotropicDamageMaterial1Status * >( this->giveStatus(gp) );
    StructuralCrossSection *crossSection = static_cast< StructuralCrossSection * >( gp->giveElement()->giveCrossSection() );

    const double e0 = this->give(e0_ID, gp);
    const double ef = this->give(ef_ID, gp);
    const double gf = this->give(gf_ID, gp);
    double wf = this->give(wf_ID, gp);


    if ( softType == ST_Disable_Damage ) {
        return;
    }

    if ( gf != 0. ) { //cohesive crack model
        if ( softType == ST_Exponential_Cohesive_Crack ) { // exponential softening
            wf = gf / E / e0; // wf is the crack opening
        } else if ( softType == ST_Linear_Cohesive_Crack || softType == ST_BiLinear_Cohesive_Crack ) { // (bi) linear softening law
            wf = 2. * gf / E / e0; // wf is the crack opening
        } else {
            OOFEM_ERROR2("Gf unsupported for softening type softType = %d", softType);
        }
    }

    crossSection->giveFullCharacteristicVector(fullstrain, gp, strainVector);
    fullstrain.resize(6);
    if ( ( kappa > e0 ) && ( status->giveDamage() == 0. ) ) {
        this->computePrincipalValDir(principalStrains, principalDir, fullstrain, principal_strain);
        // find index of max positive principal strain
        for ( int i = 2; i <= 3; i++ ) {
            if ( principalStrains.at(i) > principalStrains.at(indx) ) {
                indx = i;
            }
        }

        for ( int i = 1; i <= 3; i++ ) {
            crackPlaneNormal.at(i) = principalDir.at(i, indx);
        }

        // find index with minimal value but non-zero for plane-stress condition - this is the crack direction
        indx = 1;
        for ( int i = 2; i <= 3; i++ ) {
            if ( principalStrains.at(i) < principalStrains.at(indx) && fabs( principalStrains.at(i) ) > 1.e-10 ) {
                indx = i;
            }
        }

        for ( int i = 1; i <= 3; i++ ) {
            crackVect.at(i) = principalDir.at(i, indx);
        }

        status->setCrackVector(crackVect);

        if ( isCrackBandApproachUsed() ) { // le needed only if the crack band approach is used
            // old approach (default projection method)
            // le = gp->giveElement()->giveCharacteristicLenght(gp, crackPlaneNormal);
            // new approach, with choice of method
            le = gp->giveElement()->giveCharacteristicSize(gp, crackPlaneNormal, ecsMethod);
            // remember le in corresponding status
            status->setLe(le);
        }

        // compute and store the crack angle (just for postprocessing)
        double ca = 3.1415926 / 2.;
        if ( crackPlaneNormal.at(1) != 0.0 ) {
            ca = atan( crackPlaneNormal.at(2) / crackPlaneNormal.at(1) );
        }

        status->setCrackAngle(ca);



        if ( this->gf != 0. && e0 >= ( wf / le ) ) { // case for a given fracture energy
            OOFEM_WARNING6("Fracturing strain %e is lower than the elastic strain e0=%e, possible snap-back. Element number %d, wf %e, le %e", wf / le, e0, gp->giveElement()->giveLabel(), wf, le);
            if ( checkSnapBack ) {
                OOFEM_ERROR1("\n");
            }
        } else if ( wf == 0. && e0 >= ef ) {
            OOFEM_WARNING5( "Fracturing strain ef=%e is lower than the elastic strain e0=%f, possible snap-back. Increase fracturing strain to %f. Element number %d", ef, e0, e0, gp->giveElement()->giveLabel() );
            if ( checkSnapBack ) {
                OOFEM_ERROR1("\n");
            }
        } else if ( ef == 0. && e0 * le >= wf ) {
            OOFEM_WARNING5( "Crack opening at zero stress wf=%f is lower than the elastic displacement w0=%f, possible snap-back. Increase crack opening wf to %f. Element number %d", wf, e0 * le, e0 * le, gp->giveElement()->giveLabel() );
            if ( checkSnapBack ) {
                OOFEM_ERROR1("\n");
            }
        }
    }
}




void
IDGMaterial :: giveRealStressVector(FloatArray &answer, MatResponseForm form, GaussPoint *gp,const FloatArray &totalStrain, TimeStep *atTime)
//
// returns real stress vector in 3d stress space of receiver according to
// previous level of stress and current
// strain increment, the only way, how to correctly update gp records
//
{
    IDGMaterialStatus *status = static_cast< IDGMaterialStatus * >( this->giveStatus(gp) );
    LinearElasticMaterial *lmat = this->giveLinearElasticMaterial();
    FloatArray strainVector, reducedTotalStrainVector,totalStrainVector, strain;

    FloatMatrix de;
    double f, equivStrain, tempKappa = 0.0, omega = 0.0;
    double nlEquivStrain;
    //////////////////////////////////////////////////////////////
    int size = totalStrain.giveSize();
    strain = totalStrain;
    strain.resize(size-1);
    nlEquivStrain = totalStrain.at(size);
    ////////////////////////////////////////////////////////////////

    this->initGpForNewStep(gp);

    // subtract stress independent part
    // note: eigenStrains (temperature) is not contained in mechanical strain stored in gp
    // therefore it is necessary to subtract always the total eigen strain value
    this->giveStressDependentPartOfStrainVector(reducedTotalStrainVector, gp, totalStrain, atTime, VM_Total);

       
    // compute equivalent strain
    this->computeEquivalentStrain(equivStrain, reducedTotalStrainVector, gp, atTime);

    if ( llcriteria == idm_strainLevelCR ) {
        // compute value of loading function if strainLevel crit apply
        f = nlEquivStrain - status->giveKappa();

        if ( f <= 0.0 ) {
            // damage does not grow
            tempKappa = status->giveKappa();
            omega     = status->giveDamage();
        } else {
            // damage grow
            tempKappa = nlEquivStrain;
            this->initDamaged(tempKappa, reducedTotalStrainVector, gp);
            // evaluate damage parameter
            this->computeDamageParam(omega, nlEquivStrain, reducedTotalStrainVector, gp);
        }
    } else if ( llcriteria == idm_damageLevelCR ) {
        // evaluate damage parameter first
        tempKappa = nlEquivStrain;
        this->initDamaged(tempKappa, strain, gp);
        this->computeDamageParam(omega, tempKappa, reducedTotalStrainVector, gp);
        if ( omega < status->giveDamage() ) {
            // unloading takes place
            omega = status->giveDamage();
            //printf (".");
        }
    } else {
        _error("giveRealStressVector: unsupported loading/uloading criteria");
    }


    lmat->giveCharacteristicMatrix(de, ReducedForm, SecantStiffness, gp, atTime);
    de.times(1.0 - omega);
    answer.beProductOf(de, strain);

    answer.resize(size);
    answer.at(size) = equivStrain;
    // update gp
    status->letTempStrainVectorBe(totalStrain);
    status->letTempStressVectorBe(answer);
    status->setTempKappa(tempKappa);
    status->setTempDamage(omega);
}




MaterialStatus *
IDGMaterial :: CreateStatus(GaussPoint *gp) const
{
    return new IDGMaterialStatus(1, IDGMaterial :: domain, gp);
}


IDGMaterialStatus :: IDGMaterialStatus(int n, Domain *d, GaussPoint *g) : IsotropicDamageMaterial1Status(n, d, g)
{

}


IDGMaterialStatus :: ~IDGMaterialStatus()
{ }




void
IDGMaterialStatus :: initTempStatus()
//
// initializes temp variables according to variables form previous equlibrium state.
// builds new crackMap
//
{
    IsotropicDamageMaterial1Status :: initTempStatus();
}



void
IDGMaterialStatus :: updateYourself(TimeStep *atTime)
//
// updates variables (nonTemp variables describing situation at previous equilibrium state)
// after a new equilibrium state has been reached
// temporary variables are having values corresponding to newly reched equilibrium.
//
{
    IsotropicDamageMaterial1Status :: updateYourself(atTime);
}



contextIOResultType
IDGMaterialStatus :: saveContext(DataStream *stream, ContextMode mode, void *obj)
//
// saves full information stored in this Status
// no temp variables stored
//
{
    contextIOResultType iores;
    // save parent class status
    if ( ( iores = IsotropicDamageMaterial1Status :: saveContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    // write a raw data
    if ( !stream->write(& le, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    return CIO_OK;
}

contextIOResultType
IDGMaterialStatus :: restoreContext(DataStream *stream, ContextMode mode, void *obj)
//
// restores full information stored in stream to this Status
//
{
    contextIOResultType iores;
    // read parent class status
    if ( ( iores = IsotropicDamageMaterialStatus :: restoreContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    // read raw data
    if ( !stream->read(& le, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    return CIO_OK;
}

}     // end namespace oofem
