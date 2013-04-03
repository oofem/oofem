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

#include "lsmastermat.h"
#include "isolinearelasticmaterial.h"
#include "gausspnt.h"
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
LsMasterMat :: LsMasterMat(int n, Domain *d) : StructuralMaterial(n, d)
{
  slaveMat = 0;
}

// destructor
LsMasterMat :: ~LsMasterMat()
{ }

// specifies whether a given material mode is supported by this model
int
LsMasterMat :: hasMaterialModeCapability(MaterialMode mode)
{
    if ( mode == _3dMat_F || mode == _3dMat) {
        return 1;
    }

    return 0;
}

// reads the model parameters from the input file
IRResultType
LsMasterMat :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // required by IR_GIVE_FIELD macro
    IRResultType result;                 // required by IR_GIVE_FIELD macro


    
    IR_GIVE_OPTIONAL_FIELD(ir, slaveMat, _IFT_LsMasterMat_slaveMat); // number of slave material
    IR_GIVE_OPTIONAL_FIELD(ir, m, _IFT_LsMasterMat_m); // type of Set-Hill strain tensor

    return IRRT_OK;
}

// creates a new material status  corresponding to this class
MaterialStatus *
LsMasterMat :: CreateStatus(GaussPoint *gp) const
{
    LsMasterMatStatus *status;
    status = new LsMasterMatStatus(1, this->giveDomain(), gp, slaveMat);
    return status;
}


  /*void
LsMasterMat :: constructTransformationMatrix(FloatMatrix F, GaussPoint *gp)
{
  LsMasterMatStatus *status = static_cast< LsMasterMatStatus * >( this->giveStatus(gp) );
  FloatMatrix answer;
  answer.resize(6,6);
  FloatMatrix invF;
  invF.beInverseOf(F);
 //first row of pull back transformation matrix
    answer.at(1, 1) = invF.at(1, 1) * invF.at(1, 1);
    answer.at(1, 2) = invF.at(1, 2) * invF.at(1, 2);
    answer.at(1, 3) = invF.at(1, 3) * invF.at(1, 3);
    answer.at(1, 4) = invF.at(1, 2) * invF.at(1, 3);
    answer.at(1, 5) = invF.at(1, 1) * invF.at(1, 3);
    answer.at(1, 6) = invF.at(1, 1) * invF.at(1, 2);
    //second row of pull back transformation matrix
    answer.at(2, 1) = invF.at(2, 1) * invF.at(2, 1);
    answer.at(2, 2) = invF.at(2, 2) * invF.at(2, 2);
    answer.at(2, 3) = invF.at(2, 3) * invF.at(2, 3);
    answer.at(2, 4) = invF.at(2, 2) * invF.at(2, 3);
    answer.at(2, 5) = invF.at(2, 1) * invF.at(2, 3);
    answer.at(2, 6) = invF.at(2, 1) * invF.at(2, 2);
    //third row of pull back transformation matrix
    answer.at(3, 1) = invF.at(3, 1) * invF.at(3, 1);
    answer.at(3, 2) = invF.at(3, 2) * invF.at(3, 2);
    answer.at(3, 3) = invF.at(3, 3) * invF.at(3, 3);
    answer.at(3, 4) = invF.at(3, 2) * invF.at(3, 3);
    answer.at(3, 5) = invF.at(3, 1) * invF.at(3, 3);
    answer.at(3, 6) = invF.at(3, 1) * invF.at(3, 2);
    //fourth row of pull back transformation matrix
    answer.at(4, 1) = 2.*invF.at(2, 1) * invF.at(3, 1);
    answer.at(4, 2) = 2.*invF.at(2, 2) * invF.at(3, 2);
    answer.at(4, 3) = 2.*invF.at(2, 3) * invF.at(3, 3);
    answer.at(4, 4) = ( invF.at(2, 2) * invF.at(3, 3) + invF.at(2, 3) * invF.at(3, 2) );
    answer.at(4, 5) = ( invF.at(2, 1) * invF.at(3, 3) + invF.at(2, 3) * invF.at(3, 1) );
    answer.at(4, 6) = ( invF.at(2, 1) * invF.at(3, 2) + invF.at(2, 2) * invF.at(3, 1) );
    //fifth row of pull back transformation matrix
    answer.at(5, 1) = 2.* invF.at(1, 1) * invF.at(3, 1);
    answer.at(5, 2) = 2.*invF.at(1, 2) * invF.at(3, 2);
    answer.at(5, 3) = 2.*invF.at(1, 3) * invF.at(3, 3);
    answer.at(5, 4) = ( invF.at(1, 2) * invF.at(3, 3) + invF.at(1, 3) * invF.at(3, 2) );
    answer.at(5, 5) = ( invF.at(1, 1) * invF.at(3, 3) + invF.at(1, 3) * invF.at(3, 1) );
    answer.at(5, 6) = ( invF.at(1, 1) * invF.at(3, 2) + invF.at(1, 2) * invF.at(3, 1) );
    //sixth row of pull back transformation matrix
    answer.at(6, 1) = 2.*invF.at(1, 1) * invF.at(2, 1);
    answer.at(6, 2) = 2.*invF.at(1, 2) * invF.at(2, 2);
    answer.at(6, 3) = 2.*invF.at(1, 3) * invF.at(2, 3);
    answer.at(6, 4) = ( invF.at(1, 2) * invF.at(2, 3) + invF.at(1, 3) * invF.at(2, 2) );
    answer.at(6, 5) = ( invF.at(1, 1) * invF.at(2, 3) + invF.at(1, 3) * invF.at(2, 1) );
    answer.at(6, 6) = ( invF.at(1, 1) * invF.at(2, 2) + invF.at(1, 2) * invF.at(2, 1) );
    ///////////////////////////////////////////
    status->setTransformationMatrix(answer);

    }*/


// returns the stress vector in 3d stress space
void
LsMasterMat :: giveRealStressVector(FloatArray &answer,
                                 MatResponseForm form,
                                 GaussPoint *gp,
                                 const FloatArray &totalStrain,
                                 TimeStep *atTime)
{
    LsMasterMatStatus *status = static_cast< LsMasterMatStatus * >( this->giveStatus(gp) );
    this->initTempStatus(gp);
    
    MaterialMode mode = gp->giveMaterialMode();
    if  ( mode == _3dMat ) {
        Material *mat;
        StructuralMaterial *sMat;
        mat = domain->giveMaterial(slaveMat);
        sMat = dynamic_cast< StructuralMaterial * >(mat);
        if ( sMat == NULL ) {
            _warning2("checkConsistency: material %d has no Structural support", slaveMat);
            return;
        }
        sMat->giveRealStressVector(answer, form, gp, totalStrain, atTime);
        
        status->letTempStressVectorBe(answer);
        status->letTempStrainVectorBe(totalStrain);
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
        C.jaco_(eVals,eVecs, 40);
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
            E1 = 1./(2.*m) * ( pow(lambda1,m)-1.);
            E2 = 1./(2.*m) * ( pow(lambda2,m)-1.);
            E3 = 1./(2.*m) * ( pow(lambda3,m)-1.);
        }
        /////////////////////////////////////////////////////////////////
        //logarithmic strain
        strain.resize(3,3);
        for ( int i = 1; i < 4; i++ ) {
            for ( int j = 1; j < 4; j++ ) {
                strain.at(i,j) = E1 * eVecs.at(i,1)*eVecs.at(j,1) + E2 * eVecs.at(i,2)*eVecs.at(j,2) + E3 * eVecs.at(i,3)*eVecs.at(j,3);
            }
        }

        FloatArray strainArray(6);
        strainArray.at(1) = strain.at(1,1);
        strainArray.at(2) = strain.at(2,2);
        strainArray.at(3) = strain.at(3,3);
        strainArray.at(4) = 2*strain.at(2,3);
        strainArray.at(5) = 2*strain.at(1,3);
        strainArray.at(6) = 2*strain.at(1,2);
        ////////////////////////////////////////////////////////////////////////////////////////////
        //ask slave to compute stress
        Material *mat;
        StructuralMaterial *sMat;
        mat = domain->giveMaterial(slaveMat);
        sMat = dynamic_cast< StructuralMaterial * >(mat);
        if ( sMat == NULL ) {
            _warning2("checkConsistency: material %d has no Structural support", slaveMat);
            return;
        }
        FloatMatrix stress(3,3);
        FloatArray stressA;
        sMat->giveRealStressVector(stressA, form, gp, strainArray, atTime);
        //////////////////////////////////////////////////////////
        //  this->constructTransformationMatrix(F,gp);
        
        
        //////////////////////////////////////////////////////
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
        
        tT.beTranspositionOf(T);
        
        FloatMatrix L2(6,6);
        double gamma12, gamma13, gamma23,gamma112,gamma221,gamma113,gamma331,gamma223,gamma332,gamma;
        FloatArray stressM;
        
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
            gamma12 = gamma13 = gamma23  = 1./2.*lambda1P;
            
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
        else if( lambda1 == lambda2 ) {
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
            
            L2.at(1,5) = L2.at(5,1) = 2.* stressM.at(5)*gamma113;
            L2.at(1,6) = L2.at(6,1) = stressM.at(6)*(m-1)* pow(lambda1,m-2);
            L2.at(2,4) = L2.at(4,2) = 2.* stressM.at(4)*gamma113;
            L2.at(2,6) = L2.at(6,2) = stressM.at(6)*(m-1)* pow(lambda1,m-2);
            L2.at(3,4) = L2.at(4,3) = 2. * stressM.at(4)*gamma331;     
            L2.at(3,5) = L2.at(5,3) = 2. * stressM.at(5)*gamma331;
            L2.at(4,5) = L2.at(5,4) = stressM.at(6)*gamma113;
            L2.at(4,6) = L2.at(6,4) = stressM.at(5)*gamma113;
            L2.at(5,6) = L2.at(6,5) = stressM.at(4)*gamma113;
        }
        
        else if( lambda2 == lambda3 ) {
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
            
            L2.at(1,5) = L2.at(5,1) = 2. *  stressM.at(5)*gamma112;
            L2.at(1,6) = L2.at(6,1) = 2. *  stressM.at(6)*gamma112;
            L2.at(2,4) = L2.at(4,2) = stressM.at(4)*(m-1)* pow(lambda2,m-2);
            L2.at(2,6) = L2.at(6,2) = 2. * stressM.at(6)*gamma221;
            L2.at(3,4) = L2.at(4,3) = stressM.at(4)*(m-1)* pow(lambda2,m-2);
            L2.at(3,5) = L2.at(5,3) = 2. * stressM.at(5)*gamma221;
            L2.at(4,5) = L2.at(5,4) = stressM.at(6)*gamma221;
            L2.at(4,6) = L2.at(6,4) = stressM.at(5)*gamma221;
            L2.at(5,6) = L2.at(6,5) = stressM.at(4)*gamma221;

        }
        else if( lambda1 == lambda3 ) {
            gamma13 = 1./2.*lambda1P;
            gamma12 = (E1-E2)/(lambda1-lambda2);
            gamma23 = (E2-E3)/(lambda2-lambda3);
            gamma223 = (pow(lambda2,m-1)*(lambda2-lambda3)-2*(E2-E3))/((lambda2-lambda3)*(lambda2-lambda3));
            gamma332 = (pow(lambda3,m-1)*(lambda3-lambda2)-2*(E3-E2))/((lambda3-lambda2)*(lambda3-lambda2));
            
            L2.at(1,1) = 2*stressM.at(1)*(m-1)* pow(lambda1,m-2);
            L2.at(2,2) = 2*stressM.at(2)*(m-1)* pow(lambda2,m-2);
            L2.at(3,3) = 2*stressM.at(3)*(m-1)* pow(lambda3,m-2);
            
            L2.at(4,4) = stressM.at(2)*gamma223 + stressM.at(3)*gamma332;
            L2.at(5,5) = 1./2.*(stressM.at(1) + stressM.at(3))*(m-1)* pow(lambda3,m-2);
            L2.at(6,6) = stressM.at(1)*gamma332 + stressM.at(2)*gamma223;
            
            L2.at(1,5) = L2.at(5,1) = stressM.at(5)*(m-1)* pow(lambda3,m-2);
            L2.at(1,6) = L2.at(6,1) = 2. *stressM.at(6)*gamma332;
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
            gamma113 = ( pow(lambda1,m-1)*(lambda1-lambda3)-2*(E1-E3))/((lambda1-lambda3)*(lambda1-lambda3));
            gamma331 = ( pow(lambda3,m-1)*(lambda3-lambda1)-2*(E3-E1))/((lambda3-lambda1)*(lambda3-lambda1));
            gamma223 = ( pow(lambda2,m-1)*(lambda2-lambda3)-2*(E2-E3))/((lambda2-lambda3)*(lambda2-lambda3));
            gamma332 = ( pow(lambda3,m-1)*(lambda3-lambda2)-2*(E3-E2))/((lambda3-lambda2)*(lambda3-lambda2));
            
            gamma = (lambda1*(E2-E3)+lambda2*(E3-E1)+lambda3*(E1-E2))/((lambda1-lambda2)*(lambda2-lambda3)*(lambda3-lambda1));            
        
            L2.at(1,1) = 2*stressM.at(1)*(m-1)* pow(lambda1,m-2);
            L2.at(2,2) = 2*stressM.at(2)*(m-1)* pow(lambda2,m-2);
            L2.at(3,3) = 2*stressM.at(3)*(m-1)* pow(lambda3,m-2);;
            
            L2.at(4,4) = stressM.at(2)*gamma223 + stressM.at(3)*gamma332;     
            L2.at(5,5) = stressM.at(1)*gamma113 + stressM.at(3)*gamma331;
            L2.at(6,6) = stressM.at(1)*gamma112 + stressM.at(2)*gamma221;
            
            L2.at(1,5) = L2.at(5,1) = 2.*stressM.at(5)*gamma113;
            L2.at(1,6) = L2.at(6,1) = 2.*stressM.at(6)*gamma112;
            L2.at(2,4) = L2.at(4,2) = 2.*stressM.at(4)*gamma223;
            L2.at(2,6) = L2.at(6,2) = 2.*stressM.at(6)*gamma221;
            L2.at(3,4) = L2.at(4,3) = 2.*stressM.at(4)*gamma332;
            L2.at(3,5) = L2.at(5,3) = 2.*stressM.at(5)*gamma331;
            L2.at(4,5) = L2.at(5,4) = 2.*stressM.at(6)*gamma;
            L2.at(4,6) = L2.at(6,4) = 2.*stressM.at(5)*gamma;
            L2.at(5,6) = L2.at(6,5) = 2.*stressM.at(4)*gamma;
            
            
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
        /////////////////////
        junk.zero();
        junk.beProductOf(L2,T);
        TL.beProductOf(tT,junk);
        status ->setTransformationMatrix(F);
        status -> setPmatrix(P);
        status -> setTLmatrix(TL);
        ////////////////////////////

        status->letTempStressVectorBe(answer);
        status->letTempStrainVectorBe(GreenLagrangeStrainVector);
    }
}

 void 
 LsMasterMat :: give3dMaterialStiffnessMatrix(FloatMatrix &answer,MatResponseForm form, MatResponseMode mode,GaussPoint * gp,TimeStep * atTime)
{
    LsMasterMatStatus *status = static_cast< LsMasterMatStatus * >( this->giveStatus(gp) );
    Material *mat;
    StructuralMaterial *sMat;
    FloatMatrix stiffness;
    MaterialMode mMode = gp->giveMaterialMode();
    if ( mMode == _3dMat ) {
        Material *mat;
        StructuralMaterial *sMat;
        mat = domain->giveMaterial(slaveMat);
        sMat = dynamic_cast< StructuralMaterial * >(mat);
        if ( sMat == NULL ) {
            _warning2("checkConsistency: material %d has no Structural support", slaveMat);
            return;
        }
        sMat->give3dMaterialStiffnessMatrix(answer,form,mode,gp,atTime);
    } else {
        mat = domain->giveMaterial(slaveMat);
        sMat = dynamic_cast< StructuralMaterial * >(mat);
        if ( sMat == NULL ) {
            _warning2("checkConsistency: material %d has no Structural support", slaveMat);
            return;
        }
        sMat->give3dMaterialStiffnessMatrix(stiffness,form,mode,gp,atTime);
        FloatMatrix P,TL,F,junk;
        FloatArray stress;
        stress =   status ->giveTempStressVector();
        status -> giveTransformationMatrix(F);
        ///////////////////////////////////////////////////////////
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
        /////////////////////////////////////////////////////////////
        status->givePmatrix(P);
        status->giveTLmatrix(TL);
        junk.resize(6,6);
        junk.zero();
        junk.beProductOf(stiffness,P);
        answer.beProductOf(P,junk);
        answer.add(TL);
    }
}


int
LsMasterMat :: giveIPValue(FloatArray &answer, GaussPoint *aGaussPoint, InternalStateType type, TimeStep *atTime)
{
    LsMasterMatStatus *status = static_cast< LsMasterMatStatus * >( this->giveStatus(aGaussPoint) );

    if ( type == IST_StressTensor ) {
        answer = status->giveStressVector();
        return 1;
        } 
    else if ( type == IST_StrainTensor ) {
        answer = status->giveStrainVector();
        return 1;
        }
    else {
        Material *mat;
        StructuralMaterial *sMat;
        mat = domain->giveMaterial(slaveMat);
        sMat = dynamic_cast< StructuralMaterial * >(mat);
        if ( sMat == NULL ) {
            _warning2("checkConsistency: material %d has no Structural support", slaveMat);
            return 0;
        }
      
        int result = sMat -> giveIPValue(answer, aGaussPoint, type, atTime);
        return result;
    }
   
}

InternalStateValueType
LsMasterMat :: giveIPValueType(InternalStateType type)
{
    Material *mat;
    StructuralMaterial *sMat;
    FloatMatrix stiffness;
    mat = domain->giveMaterial(slaveMat);
    sMat = dynamic_cast< StructuralMaterial * >(mat);
    if ( sMat == NULL ) {
        _warning2("checkConsistency: material %d has no Structural support", slaveMat);
    }
    
    InternalStateValueType result = sMat -> giveIPValueType(type);
    return result;
}


int
LsMasterMat :: giveIntVarCompFullIndx(IntArray &answer, InternalStateType type, MaterialMode mmode)
{
    Material *mat;
    StructuralMaterial *sMat;
    FloatMatrix stiffness;
    mat = domain->giveMaterial(slaveMat);
    sMat = dynamic_cast< StructuralMaterial * >(mat);
    if ( sMat == NULL ) {
        _warning2("checkConsistency: material %d has no Structural support", slaveMat);
        return 1;
    }
    
    int result = sMat ->  giveIntVarCompFullIndx(answer,type,mmode);
    return result;
}

int
LsMasterMat :: giveIPValueSize(InternalStateType type, GaussPoint *aGaussPoint)
{
    Material *mat;
    StructuralMaterial *sMat;
    mat = domain->giveMaterial(slaveMat);
    sMat = dynamic_cast< StructuralMaterial * >(mat);
    if ( sMat == NULL ) {
        _warning2("checkConsistency: material %d has no Structural support", slaveMat);
        return 1;
    }
    
    int result = sMat -> giveIPValueSize(type,aGaussPoint);
    return result;
}





//=============================================================================

LsMasterMatStatus :: LsMasterMatStatus(int n, Domain *d, GaussPoint *g, int s) : StructuralMaterialStatus(n, d, g)
{
    slaveMat = s;
    Pmatrix.resize(6,6);
    Pmatrix.at(1,1) =   Pmatrix.at(2,2) =   Pmatrix.at(3,3) =   Pmatrix.at(4,4) =   Pmatrix.at(5,5) =   Pmatrix.at(6,6) = 1;
    TLmatrix.resize(6,6);
    transformationMatrix.resize(6,6);
}

LsMasterMatStatus :: ~LsMasterMatStatus()
{ }



void
LsMasterMatStatus :: printOutputAt(FILE *file, TimeStep *tStep)
{
    Material *mat;
    StructuralMaterial *sMat;
    mat = domain->giveMaterial(slaveMat);
    sMat = static_cast< StructuralMaterial * >(mat);
    MaterialStatus* mS = sMat->giveStatus(gp);

    mS->printOutputAt(file,tStep);
    //  StructuralMaterialStatus :: printOutputAt(file, tStep);

}

// initializes temporary variables based on their values at the previous equlibrium state
void LsMasterMatStatus :: initTempStatus()
{
    Material *mat;
    StructuralMaterial *sMat;
    mat = domain->giveMaterial(slaveMat);
    sMat = static_cast< StructuralMaterial * >(mat);
    MaterialStatus* mS = sMat->giveStatus(gp);
    mS->initTempStatus();
    //StructuralMaterialStatus :: initTempStatus();   
}


// updates internal variables when equilibrium is reached
void
LsMasterMatStatus :: updateYourself(TimeStep *atTime)
{
    Material *mat;
    StructuralMaterial *sMat;
    mat = domain->giveMaterial(slaveMat);
    sMat = static_cast< StructuralMaterial * >(mat);
    MaterialStatus* mS = sMat->giveStatus(gp);
    mS->updateYourself(atTime);
    //  StructuralMaterialStatus :: updateYourself(atTime);
}


// saves full information stored in this status
// temporary variables are NOT stored
contextIOResultType
LsMasterMatStatus :: saveContext(DataStream *stream, ContextMode mode, void *obj)
{
    contextIOResultType iores;
    Material *mat;
    StructuralMaterial *sMat;
    mat = domain->giveMaterial(slaveMat);
    sMat = dynamic_cast< StructuralMaterial * >(mat);
    MaterialStatus* mS = sMat->giveStatus(gp);
    // save parent class status
    if ( ( iores = mS -> saveContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    // write raw data
 
    return CIO_OK;
}



contextIOResultType
LsMasterMatStatus :: restoreContext(DataStream *stream, ContextMode mode, void *obj)
//
// restores full information stored in stream to this Status
//
{
    contextIOResultType iores;

    // read parent class status
    if ( ( iores = StructuralMaterialStatus :: restoreContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    return CIO_OK; // return succes
}

} // end namespace oofem
