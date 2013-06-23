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

#include "lsmastermatgrad.h"
#include "isolinearelasticmaterial.h"
#include "gausspoint.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "intarray.h"
#include "stressvector.h"
#include "strainvector.h"
#include "structuralcrosssection.h"
#include "mathfem.h"
#include "contextioerr.h"
#include "datastream.h"

namespace oofem {
// constructor
LsMasterMatGrad :: LsMasterMatGrad(int n, Domain *d) : LsMasterMat(n, d)
{
  slaveMat = 0;
}

// destructor
LsMasterMatGrad :: ~LsMasterMatGrad()
{ }

// specifies whether a given material mode is supported by this model
int
LsMasterMatGrad :: hasMaterialModeCapability(MaterialMode mode)
{
    if ( mode == _3dMatGrad_F || _3dMatGrad) {
        return 1;
    }

    return 0;
}

// creates a new material status  corresponding to this class
MaterialStatus *
LsMasterMatGrad :: CreateStatus(GaussPoint *gp) const
{
    LsMasterMatGradStatus *status;
    status = new LsMasterMatGradStatus(1, this->giveDomain(), gp, slaveMat);
    return status;
}


void
LsMasterMatGrad :: giveCharacteristicMatrix(FloatMatrix &answer,
                                         MatResponseForm form, MatResponseMode rMode, GaussPoint *gp, TimeStep *atTime)
//
// Returns characteristic material stiffness matrix of the receiver
//
{
    _error( "giveCharacteristicMatrix : shouldn't be called");
}


void
LsMasterMatGrad :: givePDGradMatrix_uu(FloatMatrix &answer, MatResponseForm form, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) 
{
    MaterialMode mMode = gp->giveMaterialMode();
    switch ( mMode ) {
    case _3dMatGrad_F:
    case _3dMatGrad:
        give3dMaterialStiffnessMatrix(answer, mode, gp, tStep);
        break;
    default:
        _error2( "giveCharacteristicMatrix : unknown mode (%s)", __MaterialModeToString(mMode) );
    }
}

void
LsMasterMatGrad :: givePDGradMatrix_ku(FloatMatrix &answer, MatResponseForm form, MatResponseMode mode, GaussPoint* gp, TimeStep* tStep)
{
    MaterialMode mMode = gp->giveMaterialMode();
    switch ( mMode ) {
    case _3dMatGrad_F:
    case _3dMatGrad:
        give3dKappaMatrix(answer, form, mode, gp, tStep);
        break;
    default:
        _error2( "giveCharacteristicMatrix : unknown mode (%s)", __MaterialModeToString(mMode) );
    }
}

void
LsMasterMatGrad :: givePDGradMatrix_uk(FloatMatrix &answer, MatResponseForm form, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
    MaterialMode mMode = gp->giveMaterialMode();
    switch ( mMode ) {
    case _3dMatGrad_F:
    case _3dMatGrad:
        give3dGprime(answer, form, mode, gp, tStep);
        break;
    default:
        _error2( "giveCharacteristicMatrix : unknown mode (%s)", __MaterialModeToString(mMode) );
    }
}

void
LsMasterMatGrad :: givePDGradMatrix_kk(FloatMatrix &answer, MatResponseForm form, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
    MaterialMode mMode = gp->giveMaterialMode();
    switch ( mMode ) {
    case _3dMatGrad_F:
    case _3dMatGrad:
        giveInternalLength(answer, form, mode, gp, tStep);
        break;
    default:
        _error2( "giveCharacteristicMatrix : unknown mode (%s)", __MaterialModeToString(mMode) );
    }
}

void
LsMasterMatGrad :: givePDGradMatrix_LD(FloatMatrix &answer, MatResponseForm form, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
    MaterialMode mMode = gp->giveMaterialMode();
    switch ( mMode ) {
    default:
        _error2( "giveCharacteristicMatrix : unknown mode (%s)", __MaterialModeToString(mMode) );
    }
}


void
LsMasterMatGrad :: giveInternalLength(FloatMatrix &answer, MatResponseForm form, MatResponseMode mode, GaussPoint *gp, TimeStep *atTime)
{
    this->initTempStatus(gp);
    GradDpMaterialExtensionInterface *graddpmat = dynamic_cast< GradDpMaterialExtensionInterface * >(domain->giveMaterial(slaveMat)->giveInterface(GradDpMaterialExtensionInterfaceType));
    if ( graddpmat == NULL ) {
        _warning2("checkConsistency: material %d has no Structural support", slaveMat);
        return;
    }
    graddpmat->givePDGradMatrix_kk(answer, form, mode, gp, atTime);
}

void
LsMasterMatGrad :: give3dGprime(FloatMatrix &answer, MatResponseForm form, MatResponseMode mode, GaussPoint *gp, TimeStep *atTime)
{
    LsMasterMatGradStatus *status = static_cast< LsMasterMatGradStatus * >( this->giveStatus(gp) );
    this->initTempStatus(gp);
    FloatMatrix P;
    FloatMatrix gPrime;
    MaterialMode mMode = gp->giveMaterialMode();
    GradDpMaterialExtensionInterface *graddpmat = dynamic_cast< GradDpMaterialExtensionInterface * >(domain->giveMaterial(slaveMat)->giveInterface(GradDpMaterialExtensionInterfaceType));
    if ( graddpmat == NULL ) {
        _warning2("checkConsistency: material %d has no Structural support", slaveMat);
        return;
    }

    if ( mMode == _3dMatGrad ) {
        graddpmat->givePDGradMatrix_uk(answer, form, mode, gp, atTime);
    }
    else {
        graddpmat->givePDGradMatrix_uk(gPrime, form, mode, gp, atTime);
        gPrime.at(4,1) =  2.* gPrime.at(4,1); 
        gPrime.at(5,1) =  2.* gPrime.at(5,1); 
        gPrime.at(6,1) =  2.* gPrime.at(6,1); 
        status->givePmatrix(P);
        answer.beProductOf(P,gPrime);
    }
   
}

 
void
LsMasterMatGrad :: give3dKappaMatrix(FloatMatrix &answer, MatResponseForm form, MatResponseMode mode, GaussPoint *gp, TimeStep *atTime)
{
    LsMasterMatGradStatus *status = static_cast< LsMasterMatGradStatus * >( this->giveStatus(gp) );
    this->initTempStatus(gp);
    FloatMatrix kappaMatrix,kappaTMatrix,P;
    GradDpMaterialExtensionInterface *graddpmat = dynamic_cast< GradDpMaterialExtensionInterface * >(domain->giveMaterial(slaveMat)->giveInterface(GradDpMaterialExtensionInterfaceType));
    if ( graddpmat == NULL ) {
        _warning2("checkConsistency: material %d has no Structural support", slaveMat);
        return;
    }
    MaterialMode mMode = gp->giveMaterialMode();
    if ( mMode == _3dMatGrad ) {
        graddpmat->givePDGradMatrix_ku(answer, form, mode, gp, atTime);
    } else {
        graddpmat->givePDGradMatrix_ku(kappaMatrix, form, mode, gp, atTime);
        status->givePmatrix(P);
        kappaMatrix.at(1,4) = 2.*kappaMatrix.at(1,4);
        kappaMatrix.at(1,5) = 2.*kappaMatrix.at(1,5);
        kappaMatrix.at(1,6) = 2.*kappaMatrix.at(1,6);
        kappaTMatrix.beTranspositionOf(kappaMatrix);
        kappaMatrix.beProductOf(P,kappaTMatrix);
        answer.beTranspositionOf(kappaMatrix);
    }
 }



// returns the stress vector in 3d stress space
void
LsMasterMatGrad :: giveRealStressVector(FloatArray &answer,
                                 MatResponseForm form,
                                 GaussPoint *gp,
                                 const FloatArray &totalStrain,
                                 TimeStep *atTime)
{

    LsMasterMatGradStatus *status = static_cast< LsMasterMatGradStatus * >( this->giveStatus(gp) );
    this->initTempStatus(gp);
    Material *mat;
    StructuralMaterial *sMat;
    mat = domain->giveMaterial(slaveMat);
    sMat = dynamic_cast< StructuralMaterial * >(mat);
    if ( sMat == NULL ) {
        _warning2("checkConsistency: material %d has no Structural support", slaveMat);
        return;
    }
    MaterialMode mMode = gp->giveMaterialMode();
    if ( mMode == _3dMatGrad ) {
        sMat->giveRealStressVector(answer, form, gp, totalStrain, atTime);
    }
    else {
        //store of deformation gradient into 3x3 matrix
        FloatMatrix F(3,3);
        F.at(1,1)  = totalStrain.at(1);
        F.at(2,1)  = totalStrain.at(2);
        F.at(3,1)  = totalStrain.at(3);
        F.at(1,2)  = totalStrain.at(4);
        F.at(2,2)  = totalStrain.at(5);
        F.at(3,2)  = totalStrain.at(6);
        F.at(1,3)  = totalStrain.at(7);
        F.at(2,3)  = totalStrain.at(8);
        F.at(3,3)  = totalStrain.at(9);

        //compute right Cauchy-Green, its eigenvalues and eigenvectors
        FloatMatrix C,Ft,eVecs,strain;
        FloatArray eVals;
        double lambda1,lambda2,lambda3,E1,E2,E3;
        Ft.beTranspositionOf(F);
        C.beProductOf(Ft,F);  
        C.jaco_(eVals,eVecs, 20);
        FloatMatrix GreenLagrangeStrain;
        FloatArray GreenLagrangeStrainVector(6);
        ///////////////////////////////////////////////////////////////
        GreenLagrangeStrain = C;
        GreenLagrangeStrain.at(1,1) =   GreenLagrangeStrain.at(1,1) -1;
        GreenLagrangeStrain.at(2,2) =   GreenLagrangeStrain.at(2,2) -1;
        GreenLagrangeStrain.at(3,3) =   GreenLagrangeStrain.at(3,3) -1;
        GreenLagrangeStrain.times(1./2.);
        GreenLagrangeStrainVector.at(1) = GreenLagrangeStrain.at(1,1);
        GreenLagrangeStrainVector.at(2) = GreenLagrangeStrain.at(2,2);
        GreenLagrangeStrainVector.at(3) = GreenLagrangeStrain.at(3,3);
        GreenLagrangeStrainVector.at(4) = 2*GreenLagrangeStrain.at(2,3);
        GreenLagrangeStrainVector.at(5) = 2*GreenLagrangeStrain.at(1,3);
        GreenLagrangeStrainVector.at(6) = 2*GreenLagrangeStrain.at(1,2);
        //////////////////////////////////////////////////////////

        lambda1 = eVals.at(1);
        lambda2 = eVals.at(2);
        lambda3 = eVals.at(3);
        if ( m == 0 ) {
            E1 = 1./2.*log(lambda1);
            E2 = 1./2.*log(lambda2);
            E3 = 1./2.*log(lambda3);
        }
        else {
            E1 = 1./(2*m) * ( pow(lambda1,m)-1);
            E2 = 1./(2*m) * ( pow(lambda2,m)-1);
            E3 = 1./(2*m) * ( pow(lambda3,m)-1);
        }
        /////////////////////////////////////////////////////////////////
        //logarithmic strain
        strain.resize(3,3);
        for ( int i = 1; i < 4 ; i++ ) {
            for(int j = 1; j < 4; j++ ) {
                strain.at(i,j) = E1 * eVecs.at(i,1)*eVecs.at(j,1) + E2 * eVecs.at(i,2)*eVecs.at(j,2) + E3 * eVecs.at(i,3)*eVecs.at(j,3);
            }
        }

        FloatArray strainArray(7);
        strainArray.at(1) = strain.at(1,1);
        strainArray.at(2) = strain.at(2,2);
        strainArray.at(3) = strain.at(3,3);
        strainArray.at(4) = 2*strain.at(2,3);
        strainArray.at(5) = 2*strain.at(1,3);
        strainArray.at(6) = 2*strain.at(1,2);
        strainArray.at(7) = totalStrain.at(10);
        ////////////////////////////////////////////////////////////////////////////////////////////
        //ask slave to compute stress

        FloatMatrix stress(3,3);
        FloatArray stressA;
        sMat->giveRealStressVector(stressA, form, gp, strainArray, atTime);
            /////////////////////////////////////////
        ////????????????????????????????????????????
        stressA.at(4) = 2*stressA.at(4);
        stressA.at(5) = 2*stressA.at(5);
        stressA.at(6) = 2*stressA.at(6);
        ////????????????????????????????????????????  
        ////////////////////////////Transformation matrix/////////////////////////////////////////////////////////
        FloatMatrix T(6,6);
        FloatMatrix tT(6,6);
        T.zero();
        T.at(1,1) = eVecs.at(1,1)*eVecs.at(1,1);
        T.at(1,2) = eVecs.at(2,1)*eVecs.at(2,1);
        T.at(1,3) = eVecs.at(3,1)*eVecs.at(3,1);
        T.at(1,4) = eVecs.at(2,1)*eVecs.at(3,1);
        T.at(1,5) = eVecs.at(1,1)*eVecs.at(3,1);
        T.at(1,6) = eVecs.at(1,1)*eVecs.at(2,1);

        T.at(2,1) = eVecs.at(1,2)*eVecs.at(1,2);
        T.at(2,2) = eVecs.at(2,2)*eVecs.at(2,2);
        T.at(2,3) = eVecs.at(3,2)*eVecs.at(3,2);
        T.at(2,4) = eVecs.at(2,2)*eVecs.at(3,2);
        T.at(2,5) = eVecs.at(1,2)*eVecs.at(3,2);
        T.at(2,6) = eVecs.at(1,2)*eVecs.at(2,2);

        T.at(3,1) = eVecs.at(1,3)*eVecs.at(1,3);
        T.at(3,2) = eVecs.at(2,3)*eVecs.at(2,3);
        T.at(3,3) = eVecs.at(3,3)*eVecs.at(3,3);
        T.at(3,4) = eVecs.at(2,3)*eVecs.at(3,3);
        T.at(3,5) = eVecs.at(1,3)*eVecs.at(3,3);
        T.at(3,6) = eVecs.at(1,3)*eVecs.at(2,3);

        T.at(4,1) = 2*eVecs.at(1,2)*eVecs.at(1,3);
        T.at(4,2) = 2*eVecs.at(2,2)*eVecs.at(2,3);
        T.at(4,3) = 2*eVecs.at(3,2)*eVecs.at(3,3);
        T.at(4,4) = eVecs.at(2,2)*eVecs.at(3,3) + eVecs.at(3,2)*eVecs.at(2,3);
        T.at(4,5) = eVecs.at(1,2)*eVecs.at(3,3) + eVecs.at(3,2)*eVecs.at(1,3);
        T.at(4,6) = eVecs.at(1,2)*eVecs.at(2,3) + eVecs.at(2,2)*eVecs.at(1,3);

        T.at(5,1) = 2*eVecs.at(1,1)*eVecs.at(1,3);
        T.at(5,2) = 2*eVecs.at(2,1)*eVecs.at(2,3);
        T.at(5,3) = 2*eVecs.at(3,1)*eVecs.at(3,3);
        T.at(5,4) = eVecs.at(2,1)*eVecs.at(3,3) + eVecs.at(3,1)*eVecs.at(2,3);
        T.at(5,5) = eVecs.at(1,1)*eVecs.at(3,3) + eVecs.at(3,1)*eVecs.at(1,3);
        T.at(5,6) = eVecs.at(1,1)*eVecs.at(2,3) + eVecs.at(2,1)*eVecs.at(1,3);

        T.at(6,1) = 2*eVecs.at(1,1)*eVecs.at(1,2);
        T.at(6,2) = 2*eVecs.at(2,1)*eVecs.at(2,2);
        T.at(6,3) = 2*eVecs.at(3,1)*eVecs.at(3,2);
        T.at(6,4) = eVecs.at(2,1)*eVecs.at(3,2) + eVecs.at(3,1)*eVecs.at(2,2);
        T.at(6,5) = eVecs.at(1,1)*eVecs.at(3,2) + eVecs.at(3,1)*eVecs.at(1,2);
        T.at(6,6) = eVecs.at(1,1)*eVecs.at(2,2) + eVecs.at(2,1)*eVecs.at(1,2);

        // this->constructTransformationMatrix(F,gp);

        tT.beTranspositionOf(T);

        FloatMatrix L2(6,6);
        double gamma12, gamma13, gamma23,gamma112,gamma221,gamma113,gamma331,gamma223,gamma332,gamma;
        FloatArray stressM;
        /// 7th component is the local cum pl. strain
        double kappa = stressA.at(7);
        stressA.resize(6);
        stressM.beProductOf(T,stressA);
        ////////////////////////////////// 
        stressM.at(4) = 1./2.*  stressM.at(4);
        stressM.at(5) = 1./2.*  stressM.at(5);
        stressM.at(6) = 1./2.*  stressM.at(6);
        //////////////?????????????????????????
        double lambda1P =  pow(lambda1,m-1);
        double lambda2P =  pow(lambda2,m-1);
        double lambda3P =  pow(lambda3,m-1);

        // three equal eigenvalues

        if ( ( lambda1 == lambda2 ) && ( lambda2 == lambda3 ) ) {
            gamma12 = gamma13 = gamma23  = 1./2.*lambda1P ;   

            L2.at(1,1) = 2*stressM.at(1)*(m-1)* pow(lambda1,m-2);
            L2.at(2,2) = 2*stressM.at(2)*(m-1)* pow(lambda2,m-2);
            L2.at(3,3) = 2*stressM.at(3)*(m-1)* pow(lambda3,m-2);
            
            L2.at(4,4) = 1./2.*(stressM.at(2)+ stressM.at(3))*(m-1)* pow(lambda2,m-2);
            L2.at(5,5) = 1./2.*(stressM.at(1) + stressM.at(3))*(m-1)* pow(lambda3,m-2);
            L2.at(6,6) = 1./2.*(stressM.at(1)+stressM.at(2))*(m-1)* pow(lambda1,m-2);
            
            L2.at(1,5) = L2.at(5,1) = stressM.at(5)*(m-1)* pow(lambda3,m-2);
            L2.at(1,6) = L2.at(6,1) = stressM.at(6)*(m-1)* pow(lambda1,m-2);
            L2.at(2,4) = L2.at(4,2) = stressM.at(4)*(m-1)* pow(lambda2,m-2);
            L2.at(2,6) = L2.at(6,2) = stressM.at(6)*(m-1)* pow(lambda1,m-2);
            L2.at(3,4) = L2.at(4,3) = stressM.at(4)*(m-1)* pow(lambda2,m-2);
            L2.at(3,5) = L2.at(5,3) = stressM.at(5)*(m-1)* pow(lambda3,m-2); 
            
            L2.at(4,5) = L2.at(5,4) =  1./2.*stressM.at(6)*(m-1)* pow(lambda1,m-2); 
            L2.at(4,6) = L2.at(6,4) = 1./2.*stressM.at(5)*(m-1)* pow(lambda1,m-2); ;
            L2.at(5,6) = L2.at(6,5) = 1./2.*stressM.at(4)*(m-1)* pow(lambda1,m-2); ;
        }

        //two equal eigenvalues
        else if ( lambda1 == lambda2 ) {
            gamma12  = 1./2.*lambda1P;
            gamma13 = (E1-E3)/(lambda1-lambda3);
            gamma23 = (E2-E3)/(lambda2-lambda3);
            gamma113 = (pow(lambda1,m-1)*(lambda1-lambda3)-2*(E1-E3))/((lambda1-lambda3)*(lambda1-lambda3));
            gamma331 = (pow(lambda3,m-1)*(lambda3-lambda1)-2*(E3-E1))/((lambda3-lambda1)*(lambda3-lambda1));
            

            L2.at(1,1) = 2*stressM.at(1)*(m-1)* pow(lambda1,m-2);
            L2.at(2,2) = 2*stressM.at(2)*(m-1)* pow(lambda2,m-2);
            L2.at(3,3) = 2*stressM.at(3)*(m-1)* pow(lambda3,m-2);

            L2.at(4,4) = stressM.at(2)*gamma113 + stressM.at(3)*gamma331;
            L2.at(5,5) = stressM.at(1)*gamma113 + stressM.at(3)*gamma331;
            L2.at(6,6) = 1./2.*(stressM.at(1)+stressM.at(2))*(m-1)* pow(lambda1,m-2);

            L2.at(1,5) = L2.at(5,1) = 2. * stressM.at(5)*gamma113;
            L2.at(1,6) = L2.at(6,1) = stressM.at(6)*(m-1)* pow(lambda1,m-2);
            L2.at(2,4) = L2.at(4,2) = 2. * stressM.at(4)*gamma113;
            L2.at(2,6) = L2.at(6,2) = stressM.at(6)*(m-1)* pow(lambda1,m-2);
            L2.at(3,4) = L2.at(4,3) = 2. * stressM.at(4)*gamma331;     
            L2.at(3,5) = L2.at(5,3) = 2. * stressM.at(5)*gamma331;
            L2.at(4,5) = L2.at(5,4) = stressM.at(6)*gamma113;
            L2.at(4,6) = L2.at(6,4) = stressM.at(5)*gamma113;
            L2.at(5,6) = L2.at(6,5) = stressM.at(4)*gamma113;
        }
        else if ( lambda2 == lambda3 ) {
            gamma23  = 1./2.* lambda2P;
            gamma12 = (E1-E2)/(lambda1-lambda2);
            gamma13 = (E1-E3)/(lambda1-lambda3);
            gamma112 = (pow(lambda1,m-1)*(lambda1-lambda2)-2*(E1-E2))/((lambda1-lambda2)*(lambda1-lambda2));
            gamma221 = (pow(lambda2,m-1)*(lambda2-lambda1)-2*(E2-E1))/((lambda2-lambda1)*(lambda2-lambda1));

            L2.at(1,1) = 2*stressM.at(1)*(m-1)* pow(lambda1,m-2);
            L2.at(2,2) = 2*stressM.at(2)*(m-1)* pow(lambda2,m-2);
            L2.at(3,3) = 2*stressM.at(3)*(m-1)* pow(lambda3,m-2);
            
            L2.at(4,4) = 1./2.*(stressM.at(2)+ stressM.at(3))*(m-1)* pow(lambda2,m-2);
            L2.at(5,5) = stressM.at(1)*gamma112 + stressM.at(3)*gamma221;
            L2.at(6,6) = stressM.at(1)*gamma112 + stressM.at(2)*gamma221;
        
            L2.at(1,5) = L2.at(5,1) = 2. * stressM.at(5)*gamma112;
            L2.at(1,6) = L2.at(6,1) = 2. * stressM.at(6)*gamma112;
            L2.at(2,4) = L2.at(4,2) = stressM.at(4)*(m-1)* pow(lambda2,m-2);

            L2.at(2,6) = L2.at(6,2) = 2. * stressM.at(6)*gamma221;
            L2.at(3,4) = L2.at(4,3) = stressM.at(4)*(m-1)* pow(lambda2,m-2);

            L2.at(3,5) = L2.at(5,3) = 2. * stressM.at(5)*gamma221;
            L2.at(4,5) = L2.at(5,4) = stressM.at(6)*gamma221;
            L2.at(4,6) = L2.at(6,4) = stressM.at(5)*gamma221;
            L2.at(5,6) = L2.at(6,5) = stressM.at(4)*gamma221;

        }
        else if ( lambda1 == lambda3 ) {
            gamma13 = 1./2.*lambda1P;
            gamma12 = (E1-E2)/(lambda1-lambda2);
            gamma23 = (E2-E3)/(lambda2-lambda3);
            gamma223 = (pow(lambda2,m-1)*(lambda2-lambda3)-2.*(E2-E3))/((lambda2-lambda3)*(lambda2-lambda3));
            gamma332 = (pow(lambda3,m-1)*(lambda3-lambda2)-2.*(E3-E2))/((lambda3-lambda2)*(lambda3-lambda2));
            
            L2.at(1,1) = 2. * stressM.at(1)*(m-1)* pow(lambda1,m-2);
            L2.at(2,2) = 2. * stressM.at(2)*(m-1)* pow(lambda2,m-2);
            L2.at(3,3) = 2. * stressM.at(3)*(m-1)* pow(lambda3,m-2);
            
            L2.at(4,4) = stressM.at(2)*gamma223 + stressM.at(3)*gamma332;
            L2.at(5,5) = 1./2. * (stressM.at(1) + stressM.at(3))*(m-1)* pow(lambda3,m-2);
            L2.at(6,6) = stressM.at(1)*gamma332 + stressM.at(2)*gamma223;

            L2.at(1,5) = L2.at(5,1) = stressM.at(5)*(m-1)* pow(lambda3,m-2);
            L2.at(1,6) = L2.at(6,1) = 2. * stressM.at(6)*gamma332;
            L2.at(2,4) = L2.at(4,2) = 2. * stressM.at(4)*gamma223;
            L2.at(2,6) = L2.at(6,2) = 2. * stressM.at(4)*gamma223;
            L2.at(3,4) = L2.at(4,3) = 2. * stressM.at(4)*gamma332;    
            L2.at(3,5) = L2.at(5,3) = stressM.at(5)*(m-1)* pow(lambda3,m-2);
            L2.at(4,5) = L2.at(5,4) = stressM.at(6)*gamma332;
            L2.at(4,6) = L2.at(6,4) = stressM.at(5)*gamma332;
            L2.at(5,6) = L2.at(6,5) = stressM.at(4)*gamma332;
        }
        //three different eigenvalues  
        else {
            gamma12 = (E1-E2)/(lambda1-lambda2);
            gamma13 = (E1-E3)/(lambda1-lambda3);
            gamma23 = (E2-E3)/(lambda2-lambda3);

            gamma112 = ( pow(lambda1,m-1)*(lambda1-lambda2)-2*(E1-E2))/((lambda1-lambda2)*(lambda1-lambda2));
            gamma221 = ( pow(lambda2,m-1)*(lambda2-lambda1)-2*(E2-E1))/((lambda2-lambda1)*(lambda2-lambda1));
            gamma113 = (pow(lambda1,m-1)*(lambda1-lambda3)-2*(E1-E3))/((lambda1-lambda3)*(lambda1-lambda3));
            gamma331 = (pow(lambda3,m-1)*(lambda3-lambda1)-2*(E3-E1))/((lambda3-lambda1)*(lambda3-lambda1));
            gamma223 = (pow(lambda2,m-1)*(lambda2-lambda3)-2*(E2-E3))/((lambda2-lambda3)*(lambda2-lambda3));
            gamma332 = (pow(lambda3,m-1)*(lambda3-lambda2)-2*(E3-E2))/((lambda3-lambda2)*(lambda3-lambda2));

            gamma = (lambda1*(E2-E3)+lambda2*(E3-E1)+lambda3*(E1-E2))/((lambda1-lambda2)*(lambda2-lambda3)*(lambda3-lambda1));            
            
            L2.at(1,1) = 2*stressM.at(1)*(m-1)* pow(lambda1,m-2);
            L2.at(2,2) = 2*stressM.at(2)*(m-1)* pow(lambda2,m-2);
            L2.at(3,3) = 2*stressM.at(3)*(m-1)* pow(lambda3,m-2);;

            L2.at(4,4) = stressM.at(2)*gamma223 + stressM.at(3)*gamma332;     
            L2.at(5,5) = stressM.at(1)*gamma113 + stressM.at(3)*gamma331;
            L2.at(6,6) = stressM.at(1)*gamma112 + stressM.at(2)*gamma221;

            L2.at(1,5) = L2.at(5,1) = 2. * stressM.at(5)*gamma113;
            L2.at(1,6) = L2.at(6,1) = 2. * stressM.at(6)*gamma112;
            L2.at(2,4) = L2.at(4,2) = 2. * stressM.at(4)*gamma223;
            L2.at(2,6) = L2.at(6,2) = 2. * stressM.at(6)*gamma221;
            L2.at(3,4) = L2.at(4,3) = 2. * stressM.at(4)*gamma332;
            L2.at(3,5) = L2.at(5,3) = 2. * stressM.at(5)*gamma331;
            L2.at(4,5) = L2.at(5,4) = 2. * stressM.at(6)*gamma;
            L2.at(4,6) = L2.at(6,4) = 2. * stressM.at(5)*gamma;
            L2.at(5,6) = L2.at(6,5) = 2. * stressM.at(4)*gamma;
        }
        /////////////////////////////////////////////////////////////////////////////

        FloatMatrix L1(6,6);
        L1.at(1,1) = lambda1P;
        L1.at(2,2) = lambda2P;
        L1.at(3,3) = lambda3P;
        L1.at(4,4) = gamma23;
        L1.at(5,5) = gamma13;
        L1.at(6,6) = gamma12;
        FloatMatrix junk, P,TL;
        junk.beProductOf(L1,T);
        P.beProductOf(tT,junk);
        //transformation of the stress to the 2PK stress
        answer.beProductOf(P,stressA);
        /////////////////////////////////////////////
        junk.zero();
        junk.beProductOf(L2,T);
        TL.beProductOf(tT,junk);
        GreenLagrangeStrainVector.resize(7);
        GreenLagrangeStrainVector.at(7) = totalStrain.at(9); 
        answer.resize(7);
        answer.at(7) = kappa;
        status->setPmatrix(P);
        status->setTLmatrix(TL);
        status->letTempStressVectorBe(answer);
        status->letTempStrainVectorBe(GreenLagrangeStrainVector);
    }
}

 void 
 LsMasterMatGrad :: give3dMaterialStiffnessMatrix(FloatMatrix &answer, MatResponseMode mode,GaussPoint * gp,TimeStep * atTime)
{
    LsMasterMatGradStatus *status = static_cast< LsMasterMatGradStatus * >( this->giveStatus(gp) );
    Material *mat;
    StructuralMaterial *sMat;
    FloatMatrix stiffness;
    MaterialMode mMode = gp->giveMaterialMode();
    mat = domain->giveMaterial(slaveMat);
    sMat = dynamic_cast< StructuralMaterial * >(mat);
    if ( sMat == NULL ) {
        _warning2("checkConsistency: material %d has no Structural support", slaveMat);
        return;
    }
    if ( mMode ==  _3dMatGrad ) {
        sMat->give3dMaterialStiffnessMatrix(answer, mode, gp, atTime);
    } else {
        sMat->give3dMaterialStiffnessMatrix(stiffness, mode, gp, atTime);
 

        stiffness.at(1,4) = 2.*stiffness.at(1,4);
        stiffness.at(4,1) = 2.*stiffness.at(4,1);
        stiffness.at(1,5) = 2.*stiffness.at(1,5);
        stiffness.at(5,1) = 2.*stiffness.at(5,1);
        stiffness.at(1,6) = 2.*stiffness.at(1,6);
        stiffness.at(6,1) = 2.*stiffness.at(6,1);
        stiffness.at(2,4) = 2.*stiffness.at(2,4);
        stiffness.at(4,2) = 2.*stiffness.at(4,2);
        stiffness.at(2,5) = 2.*stiffness.at(2,5);
        stiffness.at(5,2) = 2.*stiffness.at(5,2);
        stiffness.at(2,6) = 2.*stiffness.at(2,6);
        stiffness.at(6,2) = 2.*stiffness.at(6,2);
        stiffness.at(3,4) = 2.*stiffness.at(3,4);
        stiffness.at(4,3) = 2.*stiffness.at(4,3);
        stiffness.at(3,5) = 2.*stiffness.at(3,5);
        stiffness.at(5,3) = 2.*stiffness.at(5,3);
        stiffness.at(3,6) = 2.*stiffness.at(3,6);
        stiffness.at(6,3) = 2.*stiffness.at(6,3);
        stiffness.at(4,4) = 4.*stiffness.at(4,4);
        stiffness.at(4,5) = 4.*stiffness.at(4,5);
        stiffness.at(5,4) = 4.*stiffness.at(5,4);
        stiffness.at(4,6) = 4.*stiffness.at(4,6);
        stiffness.at(6,4) = 4.*stiffness.at(6,4);
        stiffness.at(5,5) = 4.*stiffness.at(5,5);
        stiffness.at(5,6) = 4.*stiffness.at(5,6);
        stiffness.at(6,5) = 4.*stiffness.at(6,5);
        stiffness.at(6,6) = 4.*stiffness.at(6,6);

        FloatMatrix P,TL,F,junk;
        FloatArray stress;
        stress =   status ->giveTempStressVector();
        status -> giveTransformationMatrix(F);
        status->givePmatrix(P);
        status->giveTLmatrix(TL);
        junk.resize(6,6);
        junk.zero();
        junk.beProductOf(stiffness,P);
        answer.beProductOf(P,junk);
        answer.add(TL); 
    }
}



IRResultType
LsMasterMatGrad :: initializeFrom(InputRecord *ir)
{
    //const char *__proc = "initializeFrom"; // required by IR_GIVE_FIELD macro
    //IRResultType result;                 // required by IR_GIVE_FIELD macro

    LsMasterMat :: initializeFrom(ir);
    
    return IRRT_OK;
}

//=============================================================================

LsMasterMatGradStatus :: LsMasterMatGradStatus(int n, Domain *d, GaussPoint *g, int s) : LsMasterMatStatus(n, d, g,s)
{ }

LsMasterMatGradStatus :: ~LsMasterMatGradStatus()
{ }


void
LsMasterMatGradStatus :: printOutputAt(FILE *file, TimeStep *tStep)
{
    LsMasterMatStatus :: printOutputAt(file, tStep);
    //LsMasterMatStatus:: printOutputAt(file, tStep);
}

// initializes temporary variables based on their values at the previous equlibrium state
void LsMasterMatGradStatus :: initTempStatus()
{
    LsMasterMatStatus :: initTempStatus();   
}


// updates internal variables when equilibrium is reached
void
LsMasterMatGradStatus :: updateYourself(TimeStep *atTime)
{
    LsMasterMatStatus :: updateYourself(atTime);
}


// saves full information stored in this status
// temporary variables are NOT stored
contextIOResultType
LsMasterMatGradStatus :: saveContext(DataStream *stream, ContextMode mode, void *obj)
{
    contextIOResultType iores;

    // save parent class status
    if ( ( iores = LsMasterMatGradStatus :: saveContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    // write raw data
 
    return CIO_OK;
}


contextIOResultType
LsMasterMatGradStatus :: restoreContext(DataStream *stream, ContextMode mode, void *obj)
//
// restores full information stored in stream to this Status
//
{
    contextIOResultType iores;

    // read parent class status
    if ( ( iores = LsMasterMatStatus :: restoreContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }
    

    return CIO_OK; // return succes
}

} // end namespace oofem
