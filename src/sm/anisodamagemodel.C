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

// This code is based on the anisotropic damage model proposed by Desmorat, Gatuingt and Ragueneau in
// their paper "Nonlocal anisotropic damage model and related computational aspects for quasi-brittle material"
// published in Engineering Fracture Mechanics 74 (2007) 1539-1560.

#include "anisodamagemodel.h"
#include "structuralmaterial.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "mathfem.h"
#include "datastream.h"
#include "contextioerr.h"
#include "dynamicinputrecord.h"
#include "isolinearelasticmaterial.h"
#include "gausspoint.h"
#include "classfactory.h"

namespace oofem {

REGISTER_Material(AnisotropicDamageMaterial);

AnisotropicDamageMaterial :: AnisotropicDamageMaterial(int n, Domain *d) : StructuralMaterial(n, d)
    //
    // constructor
    //
{
    linearElasticMaterial = new IsotropicLinearElasticMaterial(n, d);
    equivStrainType = EST_Unknown;
	a = 0.;
    A = 0.;
    kappa0 = 0.;
}

AnisotropicDamageMaterial :: ~AnisotropicDamageMaterial()
//
// destructor
//
{
    delete linearElasticMaterial;
}

int
AnisotropicDamageMaterial :: hasMaterialModeCapability(MaterialMode mode)
//
// returns whether receiver supports given mode
//
{
	return mode == _3dMat || mode == _PlaneStress ;
//	return mode == _3dMat || mode == _PlaneStress || mode == _PlaneStrain || mode == _1dMat;
}

void
AnisotropicDamageMaterial :: giveRealStressVector(FloatArray &answer, GaussPoint *gp,
                                                const FloatArray &totalStrain,
                                                TimeStep *atTime)
//
// returns real stress vector in 3d stress space of receiver according to
// previous level of stress and current strain increment, the only way,
// how to correctly update gp records
//
{
    AnisotropicDamageMaterialStatus *status = static_cast< AnisotropicDamageMaterialStatus * >( this->giveStatus(gp) );
    LinearElasticMaterial *lmat = this->giveLinearElasticMaterial();
    FloatArray strainVector, reducedTotalStrainVector;
    FloatMatrix de;
    double f, equivStrain, tempKappa = 0.0;
    MaterialMode mode = gp->giveMaterialMode();
    FloatMatrix  strainTensor;
//    double traceTempD = (this->A)*(equivStrain-this->kappa0);
    double traceTempD;
    FloatArray eVals;
    FloatMatrix eVecs;
	FloatMatrix damageTensor;

    this->initGpForNewStep(gp);

    // subtract stress independent part
    // note: eigenStrains (temperature) is not contained in mechanical strain stored in gp
    // therefore it is necessary to subtract always the total eigen strain value
    this->giveStressDependentPartOfStrainVector(reducedTotalStrainVector, gp, totalStrain, atTime, VM_Total);

    // compute equivalent strain
    this->computeEquivalentStrain(equivStrain, reducedTotalStrainVector, gp, atTime);
    FloatMatrix tempDamageTensor = status->giveDamage();
//    tempKappa = this->computeK(gp);
    tempKappa=(tempDamageTensor.at(1,1)+tempDamageTensor.at(2,2)+tempDamageTensor.at(3,3))/this->A + this->kappa0;

    FloatArray effectiveStressVector, fullEffectiveStressVector, stressVector;
    FloatMatrix effectiveStressTensor, stressTensor;
    lmat->giveStiffnessMatrix(de, ElasticStiffness, gp, atTime);
    effectiveStressVector.beProductOf(de,reducedTotalStrainVector);
    f=equivStrain-tempKappa;

    if ((tempKappa < this->kappa0) && (equivStrain < tempKappa)){ // no damage has been developed yet and the material behaves elastically
       	answer=effectiveStressVector;
//       	FloatMatrix tempDamageTensor = status->giveDamage();// changed, before: status->giveTempDamage();(18.2.14)
    }else{
		if(f <= 0.0){
//			FloatMatrix tempDamageTensor = status->giveDamage();// changed, before: status->giveTempDamage();(18.2.14)
		}else{
			AnisotropicDamageMaterial :: computeDamageTensor(tempDamageTensor, gp, totalStrain, atTime);
		}

		StructuralMaterial::giveFullSymVectorForm(fullEffectiveStressVector,effectiveStressVector, mode );
		effectiveStressTensor.beMatrixForm(fullEffectiveStressVector);
		double nu = lmat->give(NYxz, gp);
		strainTensor.resize(3,3);
		strainTensor.zero();
		if (mode == _PlaneStress)
		{
			strainTensor.at(1,1)=reducedTotalStrainVector.at(1);
			strainTensor.at(2,2)=reducedTotalStrainVector.at(2);
			strainTensor.at(3,3)=-nu * ( reducedTotalStrainVector.at(1) + reducedTotalStrainVector.at(2) ) / ( 1. - nu );
			strainTensor.at(2,3)=0.0;
			strainTensor.at(3,2)=0.0;
			strainTensor.at(1,3)=0.0;
			strainTensor.at(3,1)=0.0;
			strainTensor.at(1,2)=reducedTotalStrainVector.at(3)/2.0;
			strainTensor.at(2,1)=reducedTotalStrainVector.at(3)/2.0;
		}else{
			strainTensor.at(1,1)=reducedTotalStrainVector.at(1);
			strainTensor.at(2,2)=reducedTotalStrainVector.at(2);
			strainTensor.at(3,3)=reducedTotalStrainVector.at(3);
			strainTensor.at(2,3)=reducedTotalStrainVector.at(4)/2.0;
			strainTensor.at(3,2)=reducedTotalStrainVector.at(4)/2.0;
			strainTensor.at(1,3)=reducedTotalStrainVector.at(5)/2.0;
			strainTensor.at(3,1)=reducedTotalStrainVector.at(5)/2.0;
			strainTensor.at(1,2)=reducedTotalStrainVector.at(6)/2.0;
			strainTensor.at(2,1)=reducedTotalStrainVector.at(6)/2.0;
		}
		// Correct the trace of the damage tensor if necessary (see section 8.1 of the reference paper)
		traceTempD=computeTraceD(tempDamageTensor,strainTensor, gp);
		double effectiveStressTrace = effectiveStressTensor.at(1,1) + effectiveStressTensor.at(2,2) + effectiveStressTensor.at(3,3);
		// First term of the equation 53 of the reference paper
		FloatMatrix AuxMatrix;
		// Compute (1-D) (called ImD here)
		FloatMatrix ImD, sqrtImD;
		ImD.resize(3,3);
		ImD.zero();
		ImD.at(1,1) = ImD.at(2,2) = ImD.at(3,3) = 1;
		ImD.subtract(tempDamageTensor);
		// Compute the square root of (1-D), needed in the equation 53 of the reference paper
//		int checker1 = this->checkSymmetry(ImD);
		ImD.jaco_(eVals, eVecs, 40);
		sqrtImD.resize(3,3);
		sqrtImD.zero();
		for ( int i = 1; i <= 3; i++ ) {
			for ( int j = 1; j <= 3; j++ ) {
				for (int k =1; k<=3; k++){
					if (eVals.at(k)<0.0)
					{
						eVals.at(k)=0.0;
					}
				}
				sqrtImD.at(i, j) = sqrt(eVals.at(1)) * eVecs.at(i, 1) * eVecs.at(j, 1) + sqrt(eVals.at(2)) *eVecs.at(i, 2) * eVecs.at(j, 2) + sqrt(eVals.at(3))*eVecs.at(i, 3) * eVecs.at(j, 3);
			}
		}

		AuxMatrix.beProductOf(effectiveStressTensor,sqrtImD);
		stressTensor.beProductOf(sqrtImD,AuxMatrix);

		// Second term of the equation 53 of the reference paper
		double scalar = 0;
		AuxMatrix.zero();
		for (int i = 1; i <= 3; i++)
		{
			for (int j = 1; j <= 3; j++)
			{
				scalar += ImD.at(i,j) * effectiveStressTensor.at(i,j);
			}
		}
		scalar = scalar/(3.-traceTempD);
/*		if (fabs(3.-traceTempD)<0.001)
		{
			scalar=0.0;
		}*/
		AuxMatrix = ImD;
		AuxMatrix.times(scalar);
/*		if (equivStrain>100.0)
		{
			AuxMatrix.zero();
		}*/

		stressTensor.subtract(AuxMatrix);
		// Third term of the equation 53 of the reference paper
		AuxMatrix.zero();
		AuxMatrix.at(1,1) = AuxMatrix.at(2,2) = AuxMatrix.at(3,3) = 1./3.;
		if(effectiveStressTrace > 0) {
			AuxMatrix.times((1-traceTempD)*effectiveStressTrace);
		} else {
			AuxMatrix.times(effectiveStressTrace);
		}
/*		if (equivStrain>100.0)
		{
			AuxMatrix.zero();
		}*/
		stressTensor.add(AuxMatrix);
		stressVector.beSymVectorForm(stressTensor);
		StructuralMaterial :: giveReducedSymVectorForm(answer,stressVector,mode);

    }

    // update gp
    this->correctBigValues(stressTensor);
//    int checker20 = this->checkSymmetry(tempDamageTensor);
//    int checker21 = this->checkSymmetry(stressTensor);
//    int checker22 = this->checkSymmetry(strainTensor);
    status->letTempStrainVectorBe(totalStrain);
    status->letTempStressVectorBe(answer);
    status->setTempDamage(tempDamageTensor);
    status->setTempKappa(tempKappa);
	#ifdef keep_track_of_dissipated_energy
		status->computeWork(gp);
	#endif
}


void
AnisotropicDamageMaterial :: computeEquivalentStrain(double &kappa, const FloatArray &strain, GaussPoint *gp, TimeStep *atTime)
{
    LinearElasticMaterial *lmat = this->giveLinearElasticMaterial();

    if ( strain.isEmpty() ) {
        kappa = 0.;
        return;
    }

    if ( this->equivStrainType == EST_Mazars ) {
        double posNorm = 0.0;
        FloatArray principalStrains, fullstrain;

        StructuralMaterial :: giveFullSymVectorForm( fullstrain, strain, gp->giveMaterialMode() );

        // if plane stress mode -> compute strain in z-direction from condition of zero stress in corresponding direction
        if ( gp->giveMaterialMode() == _PlaneStress ) {
            double nu = lmat->give(NYxz, gp);
            fullstrain.at(3) = -nu * ( fullstrain.at(1) + fullstrain.at(2) ) / ( 1. - nu );
        } else if ( gp->giveMaterialMode() == _1dMat ) {
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
        double sum = 0.;
        FloatArray stress, fullStress, principalStress;
        FloatMatrix de;

        lmat->giveStiffnessMatrix(de, SecantStiffness, gp, atTime);
        stress.beProductOf(de, strain);
        StructuralMaterial :: giveFullSymVectorForm( fullStress, stress, gp->giveMaterialMode() );
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

        lmat->giveStiffnessMatrix(de, SecantStiffness, gp, atTime);
        if ( this->equivStrainType == EST_ElasticEnergy ) {
            // standard elastic energy
            stress.beProductOf(de, strain);
            sum = strain.dotProduct(stress);
        } else if ( this->equivStrainType == EST_ElasticEnergyPositiveStress ) {
            // elastic energy corresponding to positive part of stress
            FloatArray fullStress, principalStress;
            StructuralMaterial :: giveFullSymVectorForm( fullStress, stress, gp->giveMaterialMode() );
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
    } else if ( this->equivStrainType == EST_Griffith ) {
        double sum = 0.;
        FloatArray stress, fullStress, principalStress;
        FloatMatrix de;

        lmat->giveStiffnessMatrix(de, SecantStiffness, gp, atTime);
        stress.beProductOf(de, strain);
        StructuralMaterial :: giveFullSymVectorForm( fullStress, stress, gp->giveMaterialMode() );
        this->computePrincipalValues(principalStress, fullStress, principal_stress);
        for ( int i = 1; i <= 3; i++ ) {
            if ( principalStress.at(i) > 0.0 && sum < principalStress.at(i) ) {
                sum = principalStress.at(i);
            }
        }

        //Use Griffith criterion if Rankine not applied
        if (sum == 0.){
            sum = -pow(principalStress.at(1)-principalStress.at(3),2.)/8./(principalStress.at(1)+principalStress.at(3));
        }
        sum = max(sum,0.);
        kappa = sum / lmat->give('E', gp);
    } else {
        _error("computeEquivalentStrain: unknown EquivStrainType");
    }
}

// Computes Kappa according to the first damage law proposed in reference paper.
double
AnisotropicDamageMaterial :: computeK(GaussPoint *gp)
{
	double K, trace;
	FloatMatrix DamageTensor;
	AnisotropicDamageMaterialStatus *status = static_cast< AnisotropicDamageMaterialStatus * >( this->giveStatus(gp) );
	DamageTensor=status->giveDamage(); ///changed, before: status->giveTempDamage();(18.2.14)
	trace=DamageTensor.at(1,1)+DamageTensor.at(2,2)+DamageTensor.at(3,3);
	K=(1/this->A)*trace+this->kappa0;
	return K;
}

// Computes the percentage of deltaD that needs to be added so the first eigenvalue of (tempDamageTensor + alpha*deltaD) reaches Dc
// To do this, the middle point algorithm is used.
// @TODO: this algorithm is not particularly efficient and another algorithm could be implemented.
double
AnisotropicDamageMaterial :: obtainAlpha1(FloatMatrix tempDamageTensor, double deltaLambda, FloatMatrix positiveStrainTensor, double damageThreshold)
{
	double alpha_a, alpha_b, newAlpha, eps, maxDamage, size;
	FloatMatrix deltaD, positiveStrainTensorSquared, resultingDamageTensor, eVecs;
	FloatArray eVals;
	int cont;
	cont=1;
	alpha_a=0;
	alpha_b=1;
	newAlpha=(alpha_a+alpha_b)/2;
	positiveStrainTensorSquared.beProductOf(positiveStrainTensor, positiveStrainTensor);
	deltaD=positiveStrainTensorSquared;
	deltaD.times(newAlpha*deltaLambda);
	this->correctBigValues(deltaD);
	resultingDamageTensor=tempDamageTensor;
	resultingDamageTensor.add(deltaD);
//	int checker2 = this->checkSymmetry(resultingDamageTensor);
	resultingDamageTensor.jaco_(eVals, eVecs, 20);
	size=eVals.giveSize();
	maxDamage=eVals.at(1);
	for (int i=2; i<=size ;i++)
	{
		if(eVals.at(i)>maxDamage)
		{
			maxDamage=eVals.at(i);
		}
	}
	eps=maxDamage-damageThreshold;
	do
	{
		if (eps>0.0)
		{
//			alpha_a=alpha_a;
			alpha_b=newAlpha;
			newAlpha=(alpha_a+alpha_b)/2;
			deltaD=positiveStrainTensorSquared;
			deltaD.times(newAlpha*deltaLambda);
			this->correctBigValues(deltaD);
			resultingDamageTensor=tempDamageTensor;
			resultingDamageTensor.add(deltaD);
//			int checker3 = this->checkSymmetry(resultingDamageTensor);
			resultingDamageTensor.jaco_(eVals, eVecs, 20);
			size=eVals.giveSize();
			maxDamage=eVals.at(1);
			for (int i=2; i<=size ;i++)
			{
				if(eVals.at(i)>maxDamage)
				{
					maxDamage=eVals.at(i);
				}
			}
			eps=maxDamage-damageThreshold;
			cont=cont+1;
			if (cont==100)
				return newAlpha;
		}
		else
		{
			alpha_a=newAlpha;
//			alpha_b=alpha_b;
			newAlpha=(alpha_a+alpha_b)/2;
			deltaD=positiveStrainTensorSquared;
			deltaD.times(newAlpha*deltaLambda);
			this->correctBigValues(deltaD);
			resultingDamageTensor=tempDamageTensor;
			resultingDamageTensor.add(deltaD);
//			int checker4 = this->checkSymmetry(resultingDamageTensor);
			resultingDamageTensor.jaco_(eVals, eVecs, 20);
			size=eVals.giveSize();
			maxDamage=eVals.at(1);
			for (int i=2; i<=size ;i++)
			{
				if(eVals.at(i)>maxDamage)
				{
					maxDamage=eVals.at(i);
				}
			}
			eps=maxDamage-damageThreshold;
			cont=cont+1;
			if (cont==100)
				return newAlpha;
		}
	}while(fabs(eps)>1.0e-9);
	return newAlpha;
}

// Computes the percentage of deltaD that needs to be added so the second eigenvalue of (tempDamageTensor + alpha*deltaD) reaches Dc
// To do this, the middle point algorithm is used.
// @TODO: this algorithm is not particularly efficient and another algorithm could be implemented.
double
AnisotropicDamageMaterial :: obtainAlpha2(FloatMatrix tempDamageTensor, double deltaLambda, FloatMatrix positiveStrainTensor, FloatMatrix projPosStrainTensor, double damageThreshold)
{
	double alpha_a, alpha_b, newAlpha, eps, maxDamage, size, minVal;
	FloatMatrix deltaD, positiveStrainTensorSquared, resultingDamageTensor, eVecs;
	FloatArray eVals;
	int cont;
	cont=1;
	alpha_a=0;
	alpha_b=1;
	newAlpha=(alpha_a+alpha_b)/2;
	deltaD=projPosStrainTensor;
	deltaD.times(newAlpha*deltaLambda);
	this->correctBigValues(deltaD);
	resultingDamageTensor=tempDamageTensor;
	resultingDamageTensor.add(deltaD);
//	int checker5 = this->checkSymmetry(resultingDamageTensor);
	resultingDamageTensor.jaco_(eVals, eVecs, 40);
	size=eVals.giveSize();
	minVal=eVals.at(1);
	for (int i=2; i<=size ;i++)
	{
		if(eVals.at(i)<minVal)
		{
			minVal=eVals.at(i);
			maxDamage=((eVals.at(1)+eVals.at(2)+eVals.at(3))-eVals.at(i))/2;
		}
	}
	eps=maxDamage-damageThreshold;
	do
	{
		if (eps>0.0)
		{
//			alpha_a=alpha_a;
			alpha_b=newAlpha;
			newAlpha=(alpha_a+alpha_b)/2;
			deltaD=projPosStrainTensor;
			deltaD.times(newAlpha*deltaLambda);
			this->correctBigValues(deltaD);
			resultingDamageTensor=tempDamageTensor;
			resultingDamageTensor.add(deltaD);
//			int checker6 = this->checkSymmetry(resultingDamageTensor);
			resultingDamageTensor.jaco_(eVals, eVecs, 40);
			size=eVals.giveSize();
			minVal=eVals.at(1);
			for (int i=2; i<=size ;i++)
			{
				if(eVals.at(i)<minVal)
				{
					minVal=eVals.at(i);
					maxDamage=((eVals.at(1)+eVals.at(2)+eVals.at(3))-eVals.at(i))/2;
				}
			}
			eps=maxDamage-damageThreshold;
			cont=cont+1;
			if (cont==100)
				return newAlpha;
		}
		else
		{
			alpha_a=newAlpha;
//			alpha_b=alpha_b;
			newAlpha=(alpha_a+alpha_b)/2;
			deltaD=projPosStrainTensor;
			deltaD.times(newAlpha*deltaLambda);
			this->correctBigValues(deltaD);
			resultingDamageTensor=tempDamageTensor;
			resultingDamageTensor.add(deltaD);
//			int checker7 = this->checkSymmetry(resultingDamageTensor);
			resultingDamageTensor.jaco_(eVals, eVecs, 40);
			size=eVals.giveSize();
			minVal=eVals.at(1);
			for (int i=2; i<=size ;i++)
			{
				if(eVals.at(i)<minVal)
				{
					minVal=eVals.at(i);
					maxDamage=((eVals.at(1)+eVals.at(2)+eVals.at(3))-eVals.at(i))/2;
				}
			}
			eps=maxDamage-damageThreshold;
			cont=cont+1;
			if (cont==100)
				return newAlpha;
		}
	}while(fabs(eps)>1.0e-12);
	return newAlpha;
}

// Computes the percentage of deltaD that needs to be added so the third eigenvalue of (tempDamageTensor + alpha*deltaD) reaches Dc
// To do this, the middle point algorithm is used.
// @TODO: this algorithm is not particularly efficient and another algorithm could be implemented.
double
AnisotropicDamageMaterial :: obtainAlpha3(FloatMatrix tempDamageTensor, double deltaLambda, FloatMatrix positiveStrainTensor, FloatArray vec3, double damageThreshold)
{
	double alpha_a, alpha_b, newAlpha, eps, aux=0;
	FloatMatrix deltaD, positiveStrainTensorSquared, resultingDamageTensor, eVecs;
	FloatArray eVals, auxVec;
	int cont;
	cont=1;
	alpha_a=0;
	alpha_b=1;
	newAlpha=(alpha_a+alpha_b)/2;
	positiveStrainTensorSquared.beProductOf(positiveStrainTensor, positiveStrainTensor);
	for (int i=1; i<=3 ;i++)
	{
		aux=aux + vec3.at(i) * ( positiveStrainTensorSquared.at(i,1) * vec3.at(1) + positiveStrainTensorSquared.at(i,2) * vec3.at(2) + positiveStrainTensorSquared.at(i,3) * vec3.at(3) );
	}
	deltaD.beDyadicProductOf(vec3,vec3);
	deltaD.times(newAlpha*deltaLambda*aux);
	this->correctBigValues(deltaD);
	resultingDamageTensor=tempDamageTensor;
	resultingDamageTensor.add(deltaD);
//	int checker8 = this->checkSymmetry(resultingDamageTensor);
	resultingDamageTensor.jaco_(eVals, eVecs, 20);
	eps=eVals.at(1)+eVals.at(2)+eVals.at(3)-3*damageThreshold;
	do
	{
		if (eps>0.0)
		{
//			alpha_a=alpha_a;
			alpha_b=newAlpha;
			newAlpha=(alpha_a+alpha_b)/2;
			deltaD.beDyadicProductOf(vec3,vec3);
			deltaD.times(newAlpha*deltaLambda*aux);
			this->correctBigValues(deltaD);
			resultingDamageTensor=tempDamageTensor;
			resultingDamageTensor.add(deltaD);
//			int checker9 = this->checkSymmetry(resultingDamageTensor);
			resultingDamageTensor.jaco_(eVals, eVecs, 20);
			eps=eVals.at(1)+eVals.at(2)+eVals.at(3)-3*damageThreshold;
			cont=cont+1;
			if (cont==100)
				return newAlpha;
		}
		else
		{
			alpha_a=newAlpha;
//			alpha_b=alpha_b;
			newAlpha=(alpha_a+alpha_b)/2;
			deltaD.beDyadicProductOf(vec3,vec3);
			deltaD.times(newAlpha*deltaLambda*aux);
			this->correctBigValues(deltaD);
			resultingDamageTensor=tempDamageTensor;
			resultingDamageTensor.add(deltaD);
//			int checker10 = this->checkSymmetry(resultingDamageTensor);
			resultingDamageTensor.jaco_(eVals, eVecs, 20);
			eps=eVals.at(1)+eVals.at(2)+eVals.at(3)-3*damageThreshold;
			cont=cont+1;
			if (cont==100)
				return newAlpha;
		}
	}while(fabs(eps)>1.0e-9);
	return newAlpha;
}
//To check symmetry: delete this function when everything works fine
double
AnisotropicDamageMaterial :: checkSymmetry(FloatMatrix matrix)
{
	int a=0;
	int nRows=matrix.giveNumberOfRows();
	for (int i=1; i<=nRows ;i++)
	{
		for (int j=1; j<=nRows ;j++)
		{
			if (fabs(matrix.at(i,j)-matrix.at(j,i))<1.e-6) {;}
			else {a=1;}
		}
	}
	if (a==1)
	{a=1;}
	return a;
}

void
AnisotropicDamageMaterial :: correctBigValues(FloatMatrix &matrix)
{
	int nRows=matrix.giveNumberOfRows();
for ( int i = 1; i <= nRows; i++ ) {
	for ( int j = 1; j <= nRows; j++ ) {
		if (matrix.at(i,j)!=matrix.at(j,i))
		{
			double Aux=(matrix.at(i,j)+matrix.at(j,i))/2.0;
			matrix.at(i,j)=Aux;
			matrix.at(j,i)=Aux;
		}
	}
}
}

double
AnisotropicDamageMaterial :: computeTraceD(FloatMatrix tempDamageTensor, FloatMatrix strainTensor, GaussPoint *gp)
{
	AnisotropicDamageMaterialStatus *status = static_cast< AnisotropicDamageMaterialStatus * >( this->giveStatus(gp) );
	int flag = status->giveFlag();
	int tempFlag=status->giveTempFlag();

	// If flag = 0, the trace of the damage tensor has never been greater than 1 in any previous step
	if (flag==0)
	{
		if (tempFlag==0)
		{flag=0;}
		else
		{flag=1;}
	}
	else
	{
		flag=1;
	}


	double Dc=1, trD=0;

	// If flag = 0, the trace of the damage tensor has never been greater than 1 before
	if (flag==0)
	{
		if ( (strainTensor.at(1,1)+strainTensor.at(2,2)+strainTensor.at(3,3)) < 0 ) // Compression
		{
			trD=tempDamageTensor.at(1,1)+tempDamageTensor.at(2,2)+tempDamageTensor.at(3,3);
			if ( trD >= 1){	// The trace of the damage tensor is greater than 1 for the first time, then, flag turns into 1
				status->setTempFlag(1);}
		}
		else																		// Tension
		{
			if ((tempDamageTensor.at(1,1)+tempDamageTensor.at(2,2)+tempDamageTensor.at(3,3)) >= 1){
				trD=Dc;}
			else{
				trD=tempDamageTensor.at(1,1)+tempDamageTensor.at(2,2)+tempDamageTensor.at(3,3);}
		}
	}

	// If flag = 1, the trace of the damage tensor has become greater than 1 before
	if (flag==1)
	{
		if ( (strainTensor.at(1,1)+strainTensor.at(2,2)+strainTensor.at(3,3)) < 0 ) // Compression
		{
			trD=tempDamageTensor.at(1,1)+tempDamageTensor.at(2,2)+tempDamageTensor.at(3,3);
		}
		else
		{																		    // Tension
			trD=Dc;
/*			if ((tempDamageTensor.at(1,1)+tempDamageTensor.at(2,2)+tempDamageTensor.at(3,3)) >= 1){
				trD=Dc;}
			else{
				trD=tempDamageTensor.at(1,1)+tempDamageTensor.at(2,2)+tempDamageTensor.at(3,3);}
*/
		}
	}
	return trD;
}


void
AnisotropicDamageMaterial :: give3dMaterialStiffnessMatrix(FloatMatrix &answer,
                                                         MatResponseMode mode,
                                                         GaussPoint *gp,
                                                         TimeStep *atTime)
//
// Implementation of the 3D stiffness matrix, according to the equations 56 and 57 of the reference paper.
{
	AnisotropicDamageMaterialStatus *status = static_cast< AnisotropicDamageMaterialStatus * >( this->giveStatus(gp) );
	if ( mode == ElasticStiffness ) {
		this->giveLinearElasticMaterial()->give3dMaterialStiffnessMatrix(answer, mode, gp, atTime);
//		this->giveLinearElasticMaterial()->giveStiffnessMatrix(answer, mode, gp, atTime);
	} else {
		FloatArray strain = status->giveTempStrainVector();
		FloatArray totalStrain;
		FloatMatrix damageTensor, strainTensor;
		// The strain vector is turned into a tensor; for that, the elements that are out of the diagonal
		// must be divided by 2
		strainTensor.resize(3,3);
		strainTensor.zero();
		strainTensor.at(1,1)=strain.at(1);
		strainTensor.at(2,2)=strain.at(2);
		strainTensor.at(3,3)=strain.at(3);
		strainTensor.at(2,3)=strain.at(4)/2.0;
		strainTensor.at(3,2)=strain.at(4)/2.0;
		strainTensor.at(1,3)=strain.at(5)/2.0;
		strainTensor.at(3,1)=strain.at(5)/2.0;
		strainTensor.at(1,2)=strain.at(6)/2.0;
		strainTensor.at(2,1)=strain.at(6)/2.0;
		// The damage tensor is read
		damageTensor = status->giveTempDamage();
		totalStrain=status->giveTempStrainVector();
//		AnisotropicDamageMaterial :: computeDamageTensor(damageTensor, gp, totalStrain, atTime);
		AnisotropicDamageMaterial :: computeSecantOperator( answer, strainTensor, damageTensor, gp);
		for ( int j = 4; j <= 6; j++ ) {
			for (int i=1; i<=6; i++){
				answer.at(i,j)=answer.at(i,j)/2;
			}
		}
	}
}


void AnisotropicDamageMaterial :: givePlaneStressStiffMtrx(FloatMatrix &answer, MatResponseMode mode,
                                                         GaussPoint *gp, TimeStep *atTime)
{
	LinearElasticMaterial *lmat = this->giveLinearElasticMaterial();
	double nu = lmat->give(NYxz, gp);
	double E = lmat->give('E', gp);
	FloatArray totalStrain;
	AnisotropicDamageMaterialStatus *status = static_cast< AnisotropicDamageMaterialStatus * >( this->giveStatus(gp) );
	if ( mode == ElasticStiffness ) {
//		this->giveLinearElasticMaterial()->give3dMaterialStiffnessMatrix(answer, mode, gp, atTime);
		this->giveLinearElasticMaterial()->giveStiffnessMatrix(answer, mode, gp, atTime);
	} else {
		FloatArray strain = status->giveTempStrainVector();
		FloatMatrix damageTensor, strainTensor;
		// The strain vector is turned into a tensor; for that, the elements that are out of the diagonal
		// must be divided by 2
		damageTensor.resize(3,3);
		damageTensor.zero();
		strainTensor.resize(3,3);
		strainTensor.zero();
		strainTensor.at(1,1)=strain.at(1);
		strainTensor.at(2,2)=strain.at(2);
		strainTensor.at(3,3)=-nu * ( strain.at(1) + strain.at(2) ) / ( 1. - nu );
		strainTensor.at(2,3)=0.0;
		strainTensor.at(3,2)=0.0;
		strainTensor.at(1,3)=0.0;
		strainTensor.at(3,1)=0.0;
		strainTensor.at(1,2)=strain.at(3)/2.0;
		strainTensor.at(2,1)=strain.at(3)/2.0;
		// The damage tensor is read
		damageTensor = status->giveTempDamage();
		totalStrain=status->giveTempStrainVector();
		AnisotropicDamageMaterial :: computeDamageTensor(damageTensor, gp, totalStrain, atTime);
		FloatMatrix secantOperator;
		secantOperator.resize(6,6);
		AnisotropicDamageMaterial :: computeSecantOperator( secantOperator, strainTensor, damageTensor, gp);
		double C11, C12, C13, C16, C21, C22, C23, C26, C61, C62, C63, C66, q, r, s;
		C11=secantOperator.at(1,1);	C12=secantOperator.at(1,2);	C13=secantOperator.at(1,3);	C16=secantOperator.at(1,6);
		C21=secantOperator.at(2,1);	C22=secantOperator.at(2,2);	C23=secantOperator.at(2,3);	C26=secantOperator.at(2,6);
		C61=secantOperator.at(6,1);	C62=secantOperator.at(6,2);	C63=secantOperator.at(6,3);	C66=secantOperator.at(6,6);
		q=-nu/E;
		r=1/(1-C13*q);
		s=1/(1-q*C23-C23*q*q*r*C13);
		answer.resize(3,3);
		answer.at(2,1)=s*( C21 + C11 * C23 * q * r);
		answer.at(2,2)=s*(C22+C12*C23*q*r);
		answer.at(2,3)=s*(C26+C16*C23*q*r)*1/2;
		answer.at(1,1)=r*(C11+C13*q*answer.at(2,1));
		answer.at(1,2)=r*(C12+C13*q*answer.at(2,2));
		answer.at(1,3)=r*(C16+C13*q*answer.at(2,3))*1/2;
		answer.at(3,1)=C61+C63*q*(answer.at(1,1)+answer.at(2,1));
		answer.at(3,2)=C62+C63*q*(answer.at(1,2)+answer.at(2,2));
		answer.at(3,3)=(C66+C63*q*(answer.at(1,3)+answer.at(2,3)))*1/2;



	}
}

void AnisotropicDamageMaterial :: givePlaneStrainStiffMtrx(FloatMatrix &answer, MatResponseMode mode,
                                                         GaussPoint *gp, TimeStep *atTime)
{
	/*
	AnisotropicDamageMaterialStatus *status = static_cast< AnisotropicDamageMaterialStatus * >( this->giveStatus(gp) );
	if ( mode == ElasticStiffness ) {
		LinearElasticMaterial *lmat = this->giveLinearElasticMaterial();
		this->giveLinearElasticMaterial()->give3dMaterialStiffnessMatrix(answer, mode, gp, atTime);
	} else {
		FloatArray strain = status->giveTempStrainVector();
		FloatMatrix damageTensor, strainTensor;
		// The strain vector is turned into a tensor; for that, the elements that are out of the diagonal
		// must be divided by 2
		strainTensor.resize(3,3);
		strainTensor.zero();
		strainTensor.at(1,1)=strain.at(1);
		strainTensor.at(2,2)=strain.at(2);
		strainTensor.at(3,3)=0.0;
		strainTensor.at(2,3)=0.0;
		strainTensor.at(3,2)=0.0;
		strainTensor.at(1,3)=0.0;
		strainTensor.at(3,1)=0.0;
		strainTensor.at(1,2)=strain.at(3)/2.0;
		strainTensor.at(2,1)=strain.at(3)/2.0;
		// The damage tensor is read
		damageTensor = status->giveTempDamage();
		FloatMatrix secantOperator;
		AnisotropicDamageMaterial :: computeSecantOperator( secantOperator, strainTensor, damageTensor, gp);
		answer.at(1,1)=secantOperator.at(1,1);
		answer.at(1,2)=secantOperator.at(1,2);
		answer.at(1,3)=secantOperator.at(1,4);
		answer.at(2,1)=secantOperator.at(2,1);
		answer.at(2,2)=secantOperator.at(2,2);
		answer.at(2,3)=secantOperator.at(2,4);
		answer.at(3,1)=secantOperator.at(4,1);
		answer.at(3,2)=secantOperator.at(4,2);
		answer.at(3,3)=secantOperator.at(4,4);
	}
	*/
}

void AnisotropicDamageMaterial :: give1dStressStiffMtrx(FloatMatrix &answer, MatResponseMode mode,
                                                      GaussPoint *gp, TimeStep *atTime)
{
    // Implementation of the 3D stiffness matrix
    AnisotropicDamageMaterialStatus *status = static_cast< AnisotropicDamageMaterialStatus * >( this->giveStatus(gp) );
    if ( mode == ElasticStiffness ) {
		this->giveLinearElasticMaterial()->give3dMaterialStiffnessMatrix(answer, mode, gp, atTime);
    } else {
		FloatArray strain = status->giveTempStrainVector();
		if( ( strain.at(1) +strain.at(2) + strain.at(3) )  > 0) {
			//@todo eq 56
		} else {
			//@todo eq 57
		}
	}
}


void
AnisotropicDamageMaterial :: computeDamageTensor(FloatMatrix &answer, GaussPoint *gp,
	                                                const FloatArray &totalStrain,
	                                                TimeStep *atTime)
//
{

	//
	// returns real stress vector in 3d stress space of receiver according to
	// previous level of stress and current strain increment, the only way,
	// how to correctly update gp records
	//
	    AnisotropicDamageMaterialStatus *status = static_cast< AnisotropicDamageMaterialStatus * >( this->giveStatus(gp) );
	    LinearElasticMaterial *lmat = this->giveLinearElasticMaterial();
	    FloatArray strainVector, reducedTotalStrainVector;
	    FloatMatrix de;
		double Dc=1.00;
	    double f, equivStrain=0.0, tempKappa;
//	    MaterialMode mode = gp->giveMaterialMode();
	    FloatMatrix  strainTensor;
//	    double traceTempD = (this->A)*(equivStrain-this->kappa0);
	    FloatArray eVals;
	    FloatMatrix eVecs;

		FloatMatrix damageTensor;

	    this->initGpForNewStep(gp);

	    // subtract stress independent part
	    // note: eigenStrains (temperature) is not contained in mechanical strain stored in gp
	    // therefore it is necessary to subtract always the total eigen strain value
	    this->giveStressDependentPartOfStrainVector(reducedTotalStrainVector, gp, totalStrain, atTime, VM_Total);

	    // compute equivalent strain
	    this->computeEquivalentStrain(equivStrain, reducedTotalStrainVector, gp, atTime);

		// test of the damage condition: first we obtain kappa with the function computeK, and then we check the damage condition
/* 		if (this->computeK(gp) > status->giveKappa()) {
	    	tempKappa = this->computeK(gp);
	    }else{
	    	tempKappa = status->giveKappa();
	    }
*/
	    tempKappa = this->computeK(gp);
		f = equivStrain - tempKappa;
		FloatMatrix tempDamageTensor = status->giveDamage();// changed, before: status->giveTempDamage();(18.2.14)

	//	damageTensor = status->giveDamage();
		if ( f <= 0.0 ) {
			// damage does not grow
	//		tempKappa = status->giveKappa();
			answer.resize(3,3);
			answer.zero();
			answer=tempDamageTensor;
		} else {

			// damage grows
			tempKappa = equivStrain;
			double deltaLambda;
			FloatArray eVals, fullStrainVector;
			FloatMatrix eVecs, strainTensor, positiveStrainTensor, positiveStrainTensorSquared, tempDamageTensor0;
			MaterialMode mode = gp->giveMaterialMode();
			double nu = lmat->give(NYxz, gp);

			// Compute square of positive part of strain tensor
			//1.- converts strain vector to full form;
			StructuralMaterial::giveFullSymVectorForm(fullStrainVector,reducedTotalStrainVector, mode );
			// The strain vector is turned into a tensor; for that, the elements that are out of the diagonal
			// must be divided by 2
			strainTensor.resize(3,3);
			strainTensor.zero();
			if (mode == _PlaneStress)
				{
				strainTensor.at(1,1)=reducedTotalStrainVector.at(1);
				strainTensor.at(2,2)=reducedTotalStrainVector.at(2);
				strainTensor.at(3,3)=-nu * ( reducedTotalStrainVector.at(1) + reducedTotalStrainVector.at(2) ) / ( 1. - nu );
				strainTensor.at(2,3)=0.0;
				strainTensor.at(3,2)=0.0;
				strainTensor.at(1,3)=0.0;
				strainTensor.at(3,1)=0.0;
				strainTensor.at(1,2)=reducedTotalStrainVector.at(3)/2.0;
				strainTensor.at(2,1)=reducedTotalStrainVector.at(3)/2.0;
				}else{
				strainTensor.at(1,1)=reducedTotalStrainVector.at(1);
				strainTensor.at(2,2)=reducedTotalStrainVector.at(2);
				strainTensor.at(3,3)=reducedTotalStrainVector.at(3);
				strainTensor.at(2,3)=reducedTotalStrainVector.at(4)/2.0;
				strainTensor.at(3,2)=reducedTotalStrainVector.at(4)/2.0;
				strainTensor.at(1,3)=reducedTotalStrainVector.at(5)/2.0;
				strainTensor.at(3,1)=reducedTotalStrainVector.at(5)/2.0;
				strainTensor.at(1,2)=reducedTotalStrainVector.at(6)/2.0;
				strainTensor.at(2,1)=reducedTotalStrainVector.at(6)/2.0;
				}
			// computes polar decomposition and negative eigenvalues set to zero
//			int checker14 = this->checkSymmetry(strainTensor);
			strainTensor.jaco_(eVals, eVecs, 40);
			for (int i = 1; i <= 3; i++)
			{
				if(eVals.at(i) < 0)
					eVals.at(i) = 0;
			}
			// computes the positive part of the strain tensor
			positiveStrainTensor.resize(3,3);
			for ( int i = 1; i <= 3; i++ ) {
				for ( int j = 1; j <= 3; j++ ) {
					positiveStrainTensor.at(i, j) = eVals.at(1) * eVecs.at(i, 1) * eVecs.at(j, 1) + eVals.at(2) *eVecs.at(i, 2) * eVecs.at(j, 2) +  eVals.at(3) *eVecs.at(i, 3) * eVecs.at(j, 3);
				}
			}
			// computes the square of positiveStrainTensor
			positiveStrainTensorSquared.beProductOf(positiveStrainTensor, positiveStrainTensor);
			//compute delta Lambda
	//		double traceD = damageTensor.at(1,1) + damageTensor.at(2,2) + damageTensor.at(3,3);
			double traceD = tempDamageTensor.at(1,1) + tempDamageTensor.at(2,2) + tempDamageTensor.at(3,3);
			double traceTempD = (this->A)*(equivStrain-this->kappa0);
			// equation 50 of the reference paper
			deltaLambda =  (traceTempD - traceD)/equivStrain/equivStrain;
			//compute delta D: equation 48 of the reference paper
			FloatMatrix deltaD;
			deltaD = positiveStrainTensorSquared;
			deltaD.times(deltaLambda);
			this->correctBigValues(deltaD);
			// compute new damage tensor
			tempDamageTensor0=tempDamageTensor;
			tempDamageTensor0.add(deltaD);		//damage tensor is equal to tempDamageTensor+deltaD
			// The following loop implements the algorithm for a case in which the maximum damage threshold is reached,
			// in such a case, the remaining damage is projected on the other two directions available and, finally,
			// if the damage threshold is reached in two of the three possible directions, the remaining damage is projected
			// on the remaining third direction. If the threshold is reached in the three possible directions, the damage
			// tensor remains unchanged in the future and with all their eigenvalues equal to the damage threshold Dc.
			// This part of the code is based on the section 8.2 of the reference paper.
//			int checker15 = this->checkSymmetry(tempDamageTensor0);
			tempDamageTensor0.jaco_(eVals, eVecs, 20);
			if ((eVals.at(1)>(Dc)) || (eVals.at(2)>(Dc)) || (eVals.at(3)>(Dc)))
			{
				double alpha=0, deltaLambda1=0, Aux1=0, Aux2=0, Aux3=0;
				FloatMatrix deltaD1(3,3), positiveStrainTensorSquared(3,3), ProjMatrix(3,3),tempDamageTensor1(3,3), deltaD2(3,3), N11(3,3), N12(3,3), N13(3,3), N12sym(3,3), N13sym(3,3), projPosStrainTensor(3,3);
				FloatMatrix part1(3,3), part2(3,3), part3(3,3);
				FloatArray auxVals(3), auxVec1(3), auxVec2(3), auxVec3(3);

				// the percentage alpha of deltaD that needs to be added so the first eigenvalue reaches Dc is obtained
				alpha=obtainAlpha1(tempDamageTensor, deltaLambda, positiveStrainTensor, Dc);

				// deltaD1 is obtained --> deltaD1=alpha*deltaD
				positiveStrainTensorSquared.beProductOf(positiveStrainTensor, positiveStrainTensor);
				deltaD1=positiveStrainTensorSquared;
				deltaD1.times(alpha*deltaLambda);
				this->correctBigValues(deltaD1);
				tempDamageTensor1=tempDamageTensor;
				tempDamageTensor1.add(deltaD1);
				// The following lines describe the process to apply the equation 64 of the reference paper. First, the
				// eigenvalues and eigenvectors of the damage tensor resulting at the moment when the threshold is reached
				// are obtained
				// (note: the equation 64 is not correctly written in the paper, it should be as implemented here:
				// D_dot = lambda_dot * [ <e>^2-(nI·<e>^2 nI)(nI x nI) - 2(nII·<e>^2 nI)(nI x nII)_sym - 2(nIII·<e>^2 nI)(nI x nIII)_sym ]
//				int checker16 = this->checkSymmetry(tempDamageTensor1);
				tempDamageTensor1.jaco_(eVals, eVecs, 40);

				// The eigenvalues and eigenvectors are ordered, with the maximum eigenvalue being I, as its corresponding
				// eigenvector, and the other two being II and III. This is necessary so the equation 64 of the reference
				// paper can be applied
				if (eVals.at(1)>=eVals.at(2) &&  eVals.at(1)>=eVals.at(3))
				{
					auxVals.at(1)=eVals.at(1);
					auxVec1.at(1)=eVecs.at(1,1); auxVec1.at(2)=eVecs.at(2,1); auxVec1.at(3)=eVecs.at(3,1);
					auxVals.at(2)=eVals.at(2);
					auxVec2.at(1)=eVecs.at(1,2); auxVec2.at(2)=eVecs.at(2,2); auxVec2.at(3)=eVecs.at(3,2);
					auxVals.at(3)=eVals.at(3);
					auxVec3.at(1)=eVecs.at(1,3); auxVec3.at(2)=eVecs.at(2,3); auxVec3.at(3)=eVecs.at(3,3);
				}
				else if (eVals.at(2)>=eVals.at(1) &&  eVals.at(2)>=eVals.at(3))
				{
					auxVals.at(1)=eVals.at(2);
					auxVec1.at(1)=eVecs.at(1,2); auxVec1.at(2)=eVecs.at(2,2); auxVec1.at(3)=eVecs.at(3,2);
					auxVals.at(2)=eVals.at(1);
					auxVec2.at(1)=eVecs.at(1,1); auxVec2.at(2)=eVecs.at(2,1); auxVec2.at(3)=eVecs.at(3,1);
					auxVals.at(3)=eVals.at(3);
					auxVec3.at(1)=eVecs.at(1,3); auxVec3.at(2)=eVecs.at(2,3); auxVec3.at(3)=eVecs.at(3,3);
				}
				else
				{
					auxVals.at(1)=eVals.at(3);
					auxVec1.at(1)=eVecs.at(1,3); auxVec1.at(2)=eVecs.at(2,3); auxVec1.at(3)=eVecs.at(3,3);
					auxVals.at(2)=eVals.at(2);
					auxVec2.at(1)=eVecs.at(1,2); auxVec2.at(2)=eVecs.at(2,2); auxVec2.at(3)=eVecs.at(3,2);
					auxVals.at(3)=eVals.at(1);
					auxVec3.at(1)=eVecs.at(1,1); auxVec3.at(2)=eVecs.at(2,1); auxVec3.at(3)=eVecs.at(3,1);
				}

				// The symmetric part of the dyadic product of eigenvectors n1 and n2 is obtained
				N11.beDyadicProductOf(auxVec1,auxVec1);
				N12.beDyadicProductOf(auxVec1,auxVec2);
				for (int i=1; i<=3 ;i++)
				{
					for (int j=1; j<=3 ;j++)
					{
						N12sym.at(i,j)=0.5*(N12.at(i,j)+N12.at(j,i));
					}
				}

				// The symmetric part of the dyadic product of eigenvectors n1 and n3 is obtained
				N13.beDyadicProductOf(auxVec1,auxVec3);
				for (int i=1; i<=3 ;i++)
				{
					for (int j=1; j<=3 ;j++)
					{
						N13sym.at(i,j)=0.5*(N13.at(i,j)+N13.at(j,i));
					}
				}

				//The projected positive strain tensor is obtained
				for (int i=1; i<=3 ;i++)
				{
					Aux1=Aux1 + auxVec1.at(i) * ( positiveStrainTensorSquared.at(i,1) * auxVec1.at(1) + positiveStrainTensorSquared.at(i,2) * auxVec1.at(2) + positiveStrainTensorSquared.at(i,3) * auxVec1.at(3) ); //==(n1*(<e>^2*n1))     (eq. 64 of the reference paper)
					Aux2=Aux2 + auxVec2.at(i) * ( positiveStrainTensorSquared.at(i,1) * auxVec1.at(1) + positiveStrainTensorSquared.at(i,2) * auxVec1.at(2) + positiveStrainTensorSquared.at(i,3) * auxVec1.at(3) ); //==(n2*(<e>^2*n1))_sym (eq. 64 of the reference paper)
					Aux3=Aux3 + auxVec3.at(i) * ( positiveStrainTensorSquared.at(i,1) * auxVec1.at(1) + positiveStrainTensorSquared.at(i,2) * auxVec1.at(2) + positiveStrainTensorSquared.at(i,3) * auxVec1.at(3) ); //==(n3*(<e>^2*n1))_sym (eq. 64 of the reference paper)
				}
				N11.times(Aux1);
				N12sym.times(2*Aux2);
				N13sym.times(2*Aux3);

				// Finally, the expression between brackets in equation 64 is built up and called projPosStrainTensor
				projPosStrainTensor=positiveStrainTensorSquared;
				projPosStrainTensor.subtract(N11);
				projPosStrainTensor.subtract(N12sym);
				projPosStrainTensor.subtract(N13sym);

				// The following loop avoids numerical problems in the case that the trace of projPosStrainTensor is very small
				if ((projPosStrainTensor.at(1,1)+projPosStrainTensor.at(2,2)+projPosStrainTensor.at(3,3))<traceTempD*1e-10)
				{
					deltaLambda1=0.;
				}
				else
				{
					deltaLambda1=(traceTempD-(tempDamageTensor1.at(1,1)+tempDamageTensor1.at(2,2)+tempDamageTensor1.at(3,3)))/(projPosStrainTensor.at(1,1)+projPosStrainTensor.at(2,2)+projPosStrainTensor.at(3,3));
				}
				//projPosStrainTensor.symmetrized();
				deltaD2=projPosStrainTensor;
				deltaD2.times(deltaLambda1);
				this->correctBigValues(deltaD2);
				tempDamageTensor1.add(deltaD2);		//damage tensor is equal to tempDamageTensor+deltaD1+deltaD2

				// The following loop checks if after the addition of D2, any other eigenvalue of the damage tensor
				// has reached the threshold. If it has, it repeats the process, but this time projecting the
				// remaining damage on the direction of the remaining eigenvector
//				int checker17 = this->checkSymmetry(tempDamageTensor1);
				tempDamageTensor1.jaco_(eVals, eVecs, 40);
				if ((eVals.at(1)>(Dc)) || (eVals.at(2)>(Dc)) || (eVals.at(3)>(Dc)))
				{
					FloatMatrix deltaD3(3,3), projPosStrainTensor_new(3,3), tempDamageTensor2(3,3), deltaD4(3,3);
					FloatArray vec3(3);
					double alpha2=0, deltaLambda2=0, Aux4=0;
					// double val3=0;
					// Restoring the value of tempDamageTensor1 = tempDamageTensor + deltaD1
					tempDamageTensor1=tempDamageTensor;
					tempDamageTensor1.add(deltaD1);
					// the percentage alpha2 of deltaD2 that needs to be added to (tempDamageTensor+deltaD1) so the second eigenvalue
					// reaches Dc is obtained
					alpha2=obtainAlpha2(tempDamageTensor1, deltaLambda1, positiveStrainTensor, projPosStrainTensor, Dc);
					deltaD3=deltaD2;
					deltaD3.times(alpha2);
					tempDamageTensor2=tempDamageTensor1;
					tempDamageTensor2.add(deltaD3);
					// The smallest eigenvalue is detected and its eigenvector is used to build the new projPosStrainTensor
//					int checker18 = this->checkSymmetry(tempDamageTensor2);
					tempDamageTensor2.jaco_(eVals, eVecs, 40);
					if ( eVals.at(1)<=eVals.at(2) && eVals.at(1)<=eVals.at(3) )
					{
					  //val3=eVals.at(1);
						vec3.at(1)=eVecs.at(1,1);vec3.at(2)=eVecs.at(2,1);vec3.at(3)=eVecs.at(3,1);
					}
					else if ( eVals.at(2)<=eVals.at(1) && eVals.at(2)<=eVals.at(3) )
					{
					  //val3=eVals.at(2);
						vec3.at(1)=eVecs.at(1,2);vec3.at(2)=eVecs.at(2,2);vec3.at(3)=eVecs.at(3,2);
					}
					else
					{
					  //val3=eVals.at(3);
						vec3.at(1)=eVecs.at(1,3);vec3.at(2)=eVecs.at(2,3);vec3.at(3)=eVecs.at(3,3);
					}

					// The following loop computes nIII·<e>^2 nIII
					for (int i=1; i<=3 ;i++)
					{
						Aux4=Aux4 + vec3.at(i) * ( positiveStrainTensorSquared.at(i,1) * vec3.at(1) + positiveStrainTensorSquared.at(i,2) * vec3.at(2) + positiveStrainTensorSquared.at(i,3) * vec3.at(3) );
					}
					projPosStrainTensor_new.beDyadicProductOf(vec3,vec3);
					//
					projPosStrainTensor_new.times(Aux4);

					// The following loop avoids numerical problems in the case that the trace of projPosStrainTensor is very small
					if ((projPosStrainTensor_new.at(1,1)+projPosStrainTensor_new.at(2,2)+projPosStrainTensor_new.at(3,3))<traceTempD*1e-10)
					{
						deltaLambda2=0;
					}
					else
					{
						deltaLambda2=(traceTempD-(tempDamageTensor2.at(1,1)+tempDamageTensor2.at(2,2)+tempDamageTensor2.at(3,3)))/(projPosStrainTensor_new.at(1,1)+projPosStrainTensor_new.at(2,2)+projPosStrainTensor_new.at(3,3));
					}
				    deltaD4=projPosStrainTensor_new;
					deltaD4.times(deltaLambda2);
					tempDamageTensor2.add(deltaD4);		//damage tensor is equal to tempDamageTensor+deltaD1+deltaD3+deltaD4

					// The following loop checks if after the addition of D4, the remaining eigenvalue of the damage tensor
					// has reached the threshold. If it has, it computes a damage tensor with all its eigenvalues equal
					// to the damage threshold Dc
//					int checker19 = this->checkSymmetry(tempDamageTensor2);
					tempDamageTensor2.jaco_(eVals, eVecs, 40);
					if ((eVals.at(1)>(Dc)) || (eVals.at(2)>(Dc)) || (eVals.at(3)>(Dc)))
					{
						double alpha3=0;
						FloatMatrix deltaD5, tempDamageTensor3;
						tempDamageTensor2=tempDamageTensor;
						tempDamageTensor2.add(deltaD1);
						tempDamageTensor2.add(deltaD3);
						alpha3=obtainAlpha3(tempDamageTensor2, deltaLambda2, positiveStrainTensor, vec3 , Dc);
						deltaD5=deltaD4;
						deltaD5.times(alpha3);
						tempDamageTensor3=tempDamageTensor2;
						tempDamageTensor3.add(deltaD5);
						tempDamageTensor=tempDamageTensor3;
					}
					else
					{
						tempDamageTensor=tempDamageTensor2;
					}
				}
				else
				{
					tempDamageTensor=tempDamageTensor1;
				}
			}
			else
			{
				tempDamageTensor=tempDamageTensor0;
			}
			answer=tempDamageTensor;

		}
/*	    if (equivStrain>1.0)
	    {
	    	answer.resize(3,3);
	    	answer.zero();
	    	for (int i=1; i<=3 ;i++)
	    	{
	    		answer.at(i,i)=1.00;
	    	}
	    }*/
}

void
AnisotropicDamageMaterial :: computeSecantOperator(FloatMatrix &answer, FloatMatrix strainTensor, FloatMatrix damageTensor, GaussPoint *gp)
//
// Implementation of the 3D stiffness matrix, according to the equations 56 and 57 of the reference paper.
{
//    AnisotropicDamageMaterialStatus *status = static_cast< AnisotropicDamageMaterialStatus * >( this->giveStatus(gp) );
	LinearElasticMaterial *lmat = this->giveLinearElasticMaterial();
	FloatArray reducedStrainVector;
	double nu = lmat->give(NYxz, gp);
	double E = lmat->give('E', gp);
	double G, K;
	double traceD, Aux;
	FloatMatrix ImD, sqrtImD, eVecs, Imatrix;
	FloatArray eVals;
//	MaterialMode mode = gp->giveMaterialMode();
	G = E/(2.0*(1.0+nu));
	K = E/(3.0*(1.0-2.0*nu));
//	tempDamageTensor = status->giveTempDamage();
	//Compute the trace of the damage tensor, correcting it if necessary (see section 8.1 of the reference paper)
	traceD=computeTraceD(damageTensor,strainTensor, gp);
	if (fabs(3.-traceD)<0.001)
	{
		Aux=0.0;
	}else{
		Aux=(1./(3.-traceD));
	}
	// compute square root of (1-D)
	ImD.resize(3,3);
	ImD.zero();
	ImD.at(1,1) = ImD.at(2,2) = ImD.at(3,3) = 1.0;
	ImD.subtract(damageTensor);

	// computes square of positive part of strain tensor
//	int checker11=this->checkSymmetry(ImD);
	ImD.jaco_(eVals, eVecs, 40);
	sqrtImD.resize(3,3);
	for ( int i = 1; i <= 3; i++ ) {
		for ( int j = 1; j <= 3; j++ ) {
			for (int k =1; k<=3; k++){
				if (eVals.at(k)<0.0)
				{
					eVals.at(k)=0.0;
				}
			}
			sqrtImD.at(i, j) = sqrt(eVals.at(1)) * eVecs.at(i, 1) * eVecs.at(j, 1) + sqrt(eVals.at(2)) *eVecs.at(i, 2) * eVecs.at(j, 2) + sqrt(eVals.at(3))*eVecs.at(i, 3) * eVecs.at(j, 3);
		}
	}
	// To compute the expresions 56 and 57 of the reference paper, we need to work with fourth order tensors. To do this,
	// a structured called fourthOrderTensor is defined. This structure is composed by nine 3x3 FloatMatrix objects
	struct fourthOrderTensor
	{
		FloatMatrix Matrix_11kl;
		FloatMatrix Matrix_12kl;
		FloatMatrix Matrix_13kl;
		FloatMatrix Matrix_21kl;
		FloatMatrix Matrix_22kl;
		FloatMatrix Matrix_23kl;
		FloatMatrix Matrix_31kl;
		FloatMatrix Matrix_32kl;
		FloatMatrix Matrix_33kl;
	};
	// Four fourthOrderTensor structures are defined
	fourthOrderTensor Block1, Block2, Block3, secantOperator;
	Imatrix.resize(3,3);
	Imatrix.zero();
	Imatrix.at(1,1) = Imatrix.at(2,2) = Imatrix.at(3,3) = 1.0;

	// The fourthOrderTensor structures are initialised
	Block1.Matrix_11kl.resize(3,3);Block1.Matrix_12kl.resize(3,3);Block1.Matrix_13kl.resize(3,3);Block1.Matrix_21kl.resize(3,3);Block1.Matrix_22kl.resize(3,3);Block1.Matrix_23kl.resize(3,3);Block1.Matrix_31kl.resize(3,3);Block1.Matrix_32kl.resize(3,3);Block1.Matrix_33kl.resize(3,3);
	Block2.Matrix_11kl.resize(3,3);Block2.Matrix_12kl.resize(3,3);Block2.Matrix_13kl.resize(3,3);Block2.Matrix_21kl.resize(3,3);Block2.Matrix_22kl.resize(3,3);Block2.Matrix_23kl.resize(3,3);Block2.Matrix_31kl.resize(3,3);Block2.Matrix_32kl.resize(3,3);Block2.Matrix_33kl.resize(3,3);
	Block3.Matrix_11kl.resize(3,3);Block3.Matrix_12kl.resize(3,3);Block3.Matrix_13kl.resize(3,3);Block3.Matrix_21kl.resize(3,3);Block3.Matrix_22kl.resize(3,3);Block3.Matrix_23kl.resize(3,3);Block3.Matrix_31kl.resize(3,3);Block3.Matrix_32kl.resize(3,3);Block3.Matrix_33kl.resize(3,3);
	secantOperator.Matrix_11kl.resize(3,3);secantOperator.Matrix_12kl.resize(3,3);secantOperator.Matrix_13kl.resize(3,3);secantOperator.Matrix_21kl.resize(3,3);secantOperator.Matrix_22kl.resize(3,3);secantOperator.Matrix_23kl.resize(3,3);secantOperator.Matrix_31kl.resize(3,3);secantOperator.Matrix_32kl.resize(3,3);secantOperator.Matrix_33kl.resize(3,3);

	for (int k=1 ; k<=3 ; k++){
		for (int l=1 ; l<=3 ; l++){
			//The first block inside the brackets is obtained --> (1-D)^(1/2) _x_ (1-D)^(1/2)
			Block1.Matrix_11kl.at(k,l) = sqrtImD.at(1,k) * sqrtImD.at(1,l);
			Block1.Matrix_12kl.at(k,l) = sqrtImD.at(1,k) * sqrtImD.at(2,l);
			Block1.Matrix_13kl.at(k,l) = sqrtImD.at(1,k) * sqrtImD.at(3,l);
			Block1.Matrix_21kl.at(k,l) = sqrtImD.at(2,k) * sqrtImD.at(1,l);
			Block1.Matrix_22kl.at(k,l) = sqrtImD.at(2,k) * sqrtImD.at(2,l);
			Block1.Matrix_23kl.at(k,l) = sqrtImD.at(2,k) * sqrtImD.at(3,l);
			Block1.Matrix_31kl.at(k,l) = sqrtImD.at(3,k) * sqrtImD.at(1,l);
			Block1.Matrix_32kl.at(k,l) = sqrtImD.at(3,k) * sqrtImD.at(2,l);
			Block1.Matrix_33kl.at(k,l) = sqrtImD.at(3,k) * sqrtImD.at(3,l);
			//The second block inside the brackets is obtained --> ((1-D) x (1-D))/(3-trD)
			Block2.Matrix_11kl.at(k,l) = ImD.at(1,1) * ImD.at(k,l) * Aux;
			Block2.Matrix_12kl.at(k,l) = ImD.at(1,2) * ImD.at(k,l) * Aux;
			Block2.Matrix_13kl.at(k,l) = ImD.at(1,3) * ImD.at(k,l) * Aux;
			Block2.Matrix_21kl.at(k,l) = ImD.at(2,1) * ImD.at(k,l) * Aux;
			Block2.Matrix_22kl.at(k,l) = ImD.at(2,2) * ImD.at(k,l) * Aux;
			Block2.Matrix_23kl.at(k,l) = ImD.at(2,3) * ImD.at(k,l) * Aux;
			Block2.Matrix_31kl.at(k,l) = ImD.at(3,1) * ImD.at(k,l) * Aux;
			Block2.Matrix_32kl.at(k,l) = ImD.at(3,2) * ImD.at(k,l) * Aux;
			Block2.Matrix_33kl.at(k,l) = ImD.at(3,3) * ImD.at(k,l) * Aux;
			//The crossed-product of two identity tensors is obtained --> 1 x 1
			Block3.Matrix_11kl.at(k,l) = Imatrix.at(1,1) * Imatrix.at(k,l);
			Block3.Matrix_12kl.at(k,l) = Imatrix.at(1,2) * Imatrix.at(k,l);
			Block3.Matrix_13kl.at(k,l) = Imatrix.at(1,3) * Imatrix.at(k,l);
			Block3.Matrix_21kl.at(k,l) = Imatrix.at(2,1) * Imatrix.at(k,l);
			Block3.Matrix_22kl.at(k,l) = Imatrix.at(2,2) * Imatrix.at(k,l);
			Block3.Matrix_23kl.at(k,l) = Imatrix.at(2,3) * Imatrix.at(k,l);
			Block3.Matrix_31kl.at(k,l) = Imatrix.at(3,1) * Imatrix.at(k,l);
			Block3.Matrix_32kl.at(k,l) = Imatrix.at(3,2) * Imatrix.at(k,l);
			Block3.Matrix_33kl.at(k,l) = Imatrix.at(3,3) * Imatrix.at(k,l);

		}
	}
	// equation 56 of the reference paper
	if( ( strainTensor.at(1,1) +strainTensor.at(2,2) + strainTensor.at(3,3) )  > 0.0) {
		for (int k=1 ; k<=3 ; k++){
			for (int l=1 ; l<=3 ; l++){
			secantOperator.Matrix_11kl.at(k,l) = 2 * G * ( Block1.Matrix_11kl.at(k,l) - Block2.Matrix_11kl.at(k,l) ) + K * (1-traceD) * Block3.Matrix_11kl.at(k,l);
			secantOperator.Matrix_12kl.at(k,l) = 2 * G * ( Block1.Matrix_12kl.at(k,l) - Block2.Matrix_12kl.at(k,l) ) + K * (1-traceD) * Block3.Matrix_12kl.at(k,l);
			secantOperator.Matrix_13kl.at(k,l) = 2 * G * ( Block1.Matrix_13kl.at(k,l) - Block2.Matrix_13kl.at(k,l) ) + K * (1-traceD) * Block3.Matrix_13kl.at(k,l);
			secantOperator.Matrix_21kl.at(k,l) = 2 * G * ( Block1.Matrix_21kl.at(k,l) - Block2.Matrix_21kl.at(k,l) ) + K * (1-traceD) * Block3.Matrix_21kl.at(k,l);
			secantOperator.Matrix_22kl.at(k,l) = 2 * G * ( Block1.Matrix_22kl.at(k,l) - Block2.Matrix_22kl.at(k,l) ) + K * (1-traceD) * Block3.Matrix_22kl.at(k,l);
			secantOperator.Matrix_23kl.at(k,l) = 2 * G * ( Block1.Matrix_23kl.at(k,l) - Block2.Matrix_23kl.at(k,l) ) + K * (1-traceD) * Block3.Matrix_23kl.at(k,l);
			secantOperator.Matrix_31kl.at(k,l) = 2 * G * ( Block1.Matrix_31kl.at(k,l) - Block2.Matrix_31kl.at(k,l) ) + K * (1-traceD) * Block3.Matrix_31kl.at(k,l);
			secantOperator.Matrix_32kl.at(k,l) = 2 * G * ( Block1.Matrix_32kl.at(k,l) - Block2.Matrix_32kl.at(k,l) ) + K * (1-traceD) * Block3.Matrix_32kl.at(k,l);
			secantOperator.Matrix_33kl.at(k,l) = 2 * G * ( Block1.Matrix_33kl.at(k,l) - Block2.Matrix_33kl.at(k,l) ) + K * (1-traceD) * Block3.Matrix_33kl.at(k,l);
			}
		}
	// equation 57 of the reference paper
	} else {
		for (int k=1 ; k<=3 ; k++){
			for (int l=1 ; l<=3 ; l++){
			secantOperator.Matrix_11kl.at(k,l) = 2 * G * ( Block1.Matrix_11kl.at(k,l) - Block2.Matrix_11kl.at(k,l) ) + K * Block3.Matrix_11kl.at(k,l);
			secantOperator.Matrix_12kl.at(k,l) = 2 * G * ( Block1.Matrix_12kl.at(k,l) - Block2.Matrix_12kl.at(k,l) ) + K * Block3.Matrix_12kl.at(k,l);
			secantOperator.Matrix_13kl.at(k,l) = 2 * G * ( Block1.Matrix_13kl.at(k,l) - Block2.Matrix_13kl.at(k,l) ) + K * Block3.Matrix_13kl.at(k,l);
			secantOperator.Matrix_21kl.at(k,l) = 2 * G * ( Block1.Matrix_21kl.at(k,l) - Block2.Matrix_21kl.at(k,l) ) + K * Block3.Matrix_21kl.at(k,l);
			secantOperator.Matrix_22kl.at(k,l) = 2 * G * ( Block1.Matrix_22kl.at(k,l) - Block2.Matrix_22kl.at(k,l) ) + K * Block3.Matrix_22kl.at(k,l);
			secantOperator.Matrix_23kl.at(k,l) = 2 * G * ( Block1.Matrix_23kl.at(k,l) - Block2.Matrix_23kl.at(k,l) ) + K * Block3.Matrix_23kl.at(k,l);
			secantOperator.Matrix_31kl.at(k,l) = 2 * G * ( Block1.Matrix_31kl.at(k,l) - Block2.Matrix_31kl.at(k,l) ) + K * Block3.Matrix_31kl.at(k,l);
			secantOperator.Matrix_32kl.at(k,l) = 2 * G * ( Block1.Matrix_32kl.at(k,l) - Block2.Matrix_32kl.at(k,l) ) + K * Block3.Matrix_32kl.at(k,l);
			secantOperator.Matrix_33kl.at(k,l) = 2 * G * ( Block1.Matrix_33kl.at(k,l) - Block2.Matrix_33kl.at(k,l) ) + K * Block3.Matrix_33kl.at(k,l);
			}
		}
	}
	// The resulting material stiffness matrix is built
	answer.resize(6,6);
	answer.at(1,1)=secantOperator.Matrix_11kl.at(1,1);
	answer.at(1,2)=secantOperator.Matrix_11kl.at(2,2);
	answer.at(1,3)=secantOperator.Matrix_11kl.at(3,3);
	answer.at(1,4)=secantOperator.Matrix_11kl.at(2,3);
	answer.at(1,5)=secantOperator.Matrix_11kl.at(3,1);
	answer.at(1,6)=secantOperator.Matrix_11kl.at(1,2);
	answer.at(2,1)=secantOperator.Matrix_22kl.at(1,1);
	answer.at(2,2)=secantOperator.Matrix_22kl.at(2,2);
	answer.at(2,3)=secantOperator.Matrix_22kl.at(3,3);
	answer.at(2,4)=secantOperator.Matrix_22kl.at(2,3);
	answer.at(2,5)=secantOperator.Matrix_22kl.at(3,1);
	answer.at(2,6)=secantOperator.Matrix_22kl.at(1,2);
	answer.at(3,1)=secantOperator.Matrix_33kl.at(1,1);
	answer.at(3,2)=secantOperator.Matrix_33kl.at(2,2);
	answer.at(3,3)=secantOperator.Matrix_33kl.at(3,3);
	answer.at(3,4)=secantOperator.Matrix_33kl.at(2,3);
	answer.at(3,5)=secantOperator.Matrix_33kl.at(3,1);
	answer.at(3,6)=secantOperator.Matrix_33kl.at(1,2);
	answer.at(4,1)=secantOperator.Matrix_23kl.at(1,1);
	answer.at(4,2)=secantOperator.Matrix_23kl.at(2,2);
	answer.at(4,3)=secantOperator.Matrix_23kl.at(3,3);
	answer.at(4,4)=secantOperator.Matrix_23kl.at(2,3);
	answer.at(4,5)=secantOperator.Matrix_23kl.at(3,1);
	answer.at(4,6)=secantOperator.Matrix_23kl.at(1,2);
	answer.at(5,1)=secantOperator.Matrix_31kl.at(1,1);
	answer.at(5,2)=secantOperator.Matrix_31kl.at(2,2);
	answer.at(5,3)=secantOperator.Matrix_31kl.at(3,3);
	answer.at(5,4)=secantOperator.Matrix_31kl.at(2,3);
	answer.at(5,5)=secantOperator.Matrix_31kl.at(3,1);
	answer.at(5,6)=secantOperator.Matrix_31kl.at(1,2);
	answer.at(6,1)=secantOperator.Matrix_12kl.at(1,1);
	answer.at(6,2)=secantOperator.Matrix_12kl.at(2,2);
	answer.at(6,3)=secantOperator.Matrix_12kl.at(3,3);
	answer.at(6,4)=secantOperator.Matrix_12kl.at(2,3);
	answer.at(6,5)=secantOperator.Matrix_12kl.at(3,1);
	answer.at(6,6)=secantOperator.Matrix_12kl.at(1,2);

}

int
AnisotropicDamageMaterial :: giveIPValue(FloatArray &answer, GaussPoint *aGaussPoint, InternalStateType type, TimeStep *atTime)
{

	AnisotropicDamageMaterialStatus *status = static_cast< AnisotropicDamageMaterialStatus * >( this->giveStatus(aGaussPoint) );
    if ( type == IST_DamageTensor ) {
        answer.resize(6);
        answer.zero();
		answer.at(1)=status->giveDamage().at(1,1);
		answer.at(2)=status->giveDamage().at(2,2);
		answer.at(3)=status->giveDamage().at(3,3);
		answer.at(4)=status->giveDamage().at(2,3);
		answer.at(5)=status->giveDamage().at(1,3);
		answer.at(6)=status->giveDamage().at(1,2);
        return 1;
    } else if ( type == IST_PrincipalDamageTensor ) {
        FloatArray eVals;
		FloatMatrix eVecs;
		answer.resize(3);
        answer.zero();
//        int checker12=this->checkSymmetry(status->giveDamage());
		status->giveDamage().jaco_(eVals, eVecs, 20);
        answer.at(1) = eVals(1);
		answer.at(2) = eVals(2);
		answer.at(3) = eVals(3);
        return 1;
    } else if ( type == IST_DamageTensorTemp ) {
        answer.resize(6);
        answer.zero();
		answer.at(1)=status->giveTempDamage().at(1,1);
		answer.at(2)=status->giveTempDamage().at(2,2);
		answer.at(3)=status->giveTempDamage().at(3,3);
		answer.at(4)=status->giveTempDamage().at(2,3);
		answer.at(5)=status->giveTempDamage().at(1,3);
		answer.at(6)=status->giveTempDamage().at(1,2);
        return 1;
    } else if ( type == IST_PrincipalDamageTempTensor ) {
        FloatArray eVals;
		FloatMatrix eVecs;
		answer.resize(3);
        answer.zero();
//        int checker13=this->checkSymmetry(status->giveTempDamage());
		status->giveTempDamage().jaco_(eVals, eVecs, 20);
        answer.at(1) = eVals(1);
		answer.at(2) = eVals(2);
		answer.at(3) = eVals(3);
        return 1;
    } else if ( type == IST_MaxEquivalentStrainLevel ) {
        answer.resize(1);
        answer.at(1) = status->giveKappa();
        return 1;
    
#ifdef keep_track_of_dissipated_energy
    } else if ( type == IST_StressWorkDensity ) {
        answer.resize(1);
        answer.at(1) = status->giveStressWork();
        return 1;
    } else if ( type == IST_DissWorkDensity ) {
        answer.resize(1);
        answer.at(1) = status->giveDissWork();
    } else if ( type == IST_FreeEnergyDensity ) {
        answer.resize(1);
        answer.at(1) = status->giveStressWork() - status->giveDissWork();
        return 1;

#endif
    } else {
        return StructuralMaterial :: giveIPValue(answer, aGaussPoint, type, atTime);
    }

    return 1; // to make the compiler happy

	}
/*
InternalStateValueType
AnisotropicDamageMaterial :: giveIPValueType(InternalStateType type)
{
    if ( type == IST_DamageTensor || type == IST_DamageTensorTemp ) {
        return ISVT_TENSOR_S3;
    } else if ( type == IST_MaxEquivalentStrainLevel ) {
        return ISVT_SCALAR;
    } else if (  type == IST_PrincipalDamageTensor || type == IST_PrincipalDamageTempTensor ) {
        return ISVT_VECTOR;

#ifdef keep_track_of_dissipated_energy
    } else if ( type == IST_DissWorkDensity ) {
        return ISVT_SCALAR;
    } else if ( type == IST_StressWorkDensity ) {
        return ISVT_SCALAR;
    } else if ( type == IST_FreeEnergyDensity ) {
        return ISVT_SCALAR;

#endif
    } else {
        return StructuralMaterial :: giveIPValueType(type);
    }

}
*/
IRResultType
AnisotropicDamageMaterial :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    linearElasticMaterial->initializeFrom(ir);

    int equivStrainType;
    // specify the type of formula for equivalent strain
	IR_GIVE_OPTIONAL_FIELD(ir, equivStrainType, _IFT_AnisotropicDamageMaterial_equivStrainType);
    if ( equivStrainType == 1 ) {
        this->equivStrainType = EST_Rankine_Smooth;
    } else if ( equivStrainType == 2 ) {
        this->equivStrainType = EST_ElasticEnergy;
    } else if ( equivStrainType == 3 ) {
        this->equivStrainType = EST_Mises;
       // IR_GIVE_FIELD(ir, k, _IFT_IsotropicDamageMaterial1_k);
    } else if ( equivStrainType == 4 ) {
        this->equivStrainType = EST_Rankine_Standard;
    } else if ( equivStrainType == 5 ) {
        this->equivStrainType = EST_ElasticEnergyPositiveStress;
    } else if ( equivStrainType == 6 ) {
        this->equivStrainType = EST_ElasticEnergyPositiveStrain;
    } else if ( equivStrainType == 7 ) {
        this->equivStrainType = EST_Griffith;
    } else {
        this->equivStrainType = EST_Mazars;
	}

    IR_GIVE_FIELD(ir, A, _IFT_AnisotropicDamageMaterial_A);
    IR_GIVE_FIELD(ir, kappa0, _IFT_AnisotropicDamageMaterial_kappa0);

    return StructuralMaterial :: initializeFrom(ir);
}

void
AnisotropicDamageMaterial :: giveInputRecord(DynamicInputRecord &input)
{
    StructuralMaterial :: giveInputRecord(input);
    input.setField(this->kappa0, _IFT_AnisotropicDamageMaterial_kappa0);
    input.setField(this->A, _IFT_AnisotropicDamageMaterial_A);
}


AnisotropicDamageMaterialStatus :: AnisotropicDamageMaterialStatus(int n, Domain *d, GaussPoint *g) : StructuralMaterialStatus(n, d, g)
{
    kappa = tempKappa = 0.0;
	damage.resize(3,3);
	damage.zero();
	tempDamage.resize(3,3);
	tempDamage.zero();
	flag = tempFlag = 0;

#ifdef keep_track_of_dissipated_energy
    stressWork = tempStressWork = 0.0;
    dissWork = tempDissWork = 0.0;
#endif
}

AnisotropicDamageMaterialStatus :: ~AnisotropicDamageMaterialStatus()
{ }

void
AnisotropicDamageMaterialStatus :: printOutputAt(FILE *file, TimeStep *tStep)
{
    //StructuralMaterialStatus :: printOutputAt(file, tStep);
/*    fprintf(file, "status { ");
    if ( this->kappa > 0 && this->damage <= 0 ) {
        fprintf(file, "kappa %f", this->kappa);
    } else if ( this->damage > 0.0 ) {
//        fprintf( file, "kappa %f, damage %f crackVector %f %f %f", this->kappa, this->damage, this->crackVector.at(1), this->crackVector.at(2), this->crackVector.at(3) );

#ifdef keep_track_of_dissipated_energy
        fprintf(file, ", dissW %f, freeE %f, stressW %f ", this->dissWork, ( this->stressWork ) - ( this->dissWork ), this->stressWork);
    } else {
        fprintf(file, "stressW %f ", this->stressWork);
#endif
    }
	*/
//	fprintf( file, "kappa %f, damage %f %f %f", this->kappa, this->damage);

    //fprintf(file, "}\n");
}


void
AnisotropicDamageMaterialStatus :: initTempStatus()
{
    StructuralMaterialStatus :: initTempStatus();
    this->tempKappa = this->kappa;
    this->tempDamage = this->damage;
#ifdef keep_track_of_dissipated_energy
    this->tempStressWork = this->stressWork;
    this->tempDissWork = this->dissWork;
#endif
}


void
AnisotropicDamageMaterialStatus :: updateYourself(TimeStep *atTime)
{
    StructuralMaterialStatus :: updateYourself(atTime);
    this->kappa = this->tempKappa;
    this->damage = this->tempDamage;
    this->flag = this->tempFlag;
#ifdef keep_track_of_dissipated_energy
    this->stressWork = this->tempStressWork;
    this->dissWork = this->tempDissWork;
#endif
}


contextIOResultType
AnisotropicDamageMaterialStatus :: saveContext(DataStream *stream, ContextMode mode, void *obj)
{
/*    contextIOResultType iores;

    // save parent class status
    if ( ( iores = StructuralMaterialStatus :: saveContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    // write raw data
    if ( !stream->write(& kappa, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream->write(& damage, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

#ifdef keep_track_of_dissipated_energy
    if ( !stream->write(& stressWork, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream->write(& dissWork, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

#endif
*/
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    contextIOResultType iores;

    if ( stream == NULL ) {
        _error("saveContext : can't write into NULL stream");
    }

    // save parent class status
    if ( ( iores = StructuralMaterialStatus :: saveContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    // write raw data
    if ( !stream->write(& kappa, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }


    // write damage (vector)
    if ( ( iores = damage.storeYourself(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

#ifdef keep_track_of_dissipated_energy
    if ( !stream->write(& stressWork, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream->write(& dissWork, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }


#endif

    return CIO_OK;
}


contextIOResultType
AnisotropicDamageMaterialStatus :: restoreContext(DataStream *stream, ContextMode mode, void *obj)
{
/*    contextIOResultType iores;

    // read parent class status
    if ( ( iores = StructuralMaterialStatus :: restoreContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    // read raw data
    if ( !stream->read(& kappa, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream->read(& damage, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

#ifdef keep_track_of_dissipated_energy
    if ( !stream->read(& stressWork, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream->read(& dissWork, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

#endif
	*/

    contextIOResultType iores;


    // read parent class status
    if ( ( iores = StructuralMaterialStatus :: restoreContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    // read raw data
    if ( !stream->read(& kappa, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }


    // read damage (vector)
    if ( ( iores = damage.restoreYourself(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

#ifdef keep_track_of_dissipated_energy
    if ( !stream->read(& stressWork, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream->read(& dissWork, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

#endif


    return CIO_OK;
}

#ifdef keep_track_of_dissipated_energy
void
AnisotropicDamageMaterialStatus :: computeWork(GaussPoint *gp)
{
	/*
    // strain increment
    FloatArray deps;
    deps.beDifferenceOf(tempStrainVector, strainVector);

    // increment of stress work density
    double dSW = ( tempStressVector.dotProduct(deps) + stressVector.dotProduct(deps) ) / 2.;
    tempStressWork = stressWork + dSW;

    // elastically stored energy density
    double We = tempStressVector.dotProduct(tempStrainVector) / 2.;

    // dissipative work density
    tempDissWork = tempStressWork - We;
	*/
}
#endif
} // end namespace oofem
