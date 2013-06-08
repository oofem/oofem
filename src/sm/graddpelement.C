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


#include "graddpelement.h"
#include "node.h"
#include "material.h"
#include "gausspoint.h"
#include "gaussintegrationrule.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "intarray.h"
#include "domain.h"
#include "cltypes.h"
#include "structuralms.h"
#include "mathfem.h"
#include "structuralcrosssection.h"
#include "nlstructuralelement.h"
#include "nonlocalbarrier.h"
#include "nlstructuralelement.h"

#include <cstdio>

namespace oofem {

GradDpElement :: GradDpElement () 
  // Constructor.
{
  nlGeo = 0;
  averType = 0;
}
 

void
GradDpElement :: setDisplacementLocationArray(IntArray &answer, int nPrimNodes, int nPrimVars, int nSecNodes, int nSecVars)
{
 
  answer.resize(locSize);
 
  for ( int i = 1; i <= totalSize; i++ ) {
      if(i<nSecNodes*nPrimVars+1)
            answer.at(i) = i +(int)(((i-1)/nPrimVars))*nSecVars;
      else if ( i > nSecNodes*(nPrimVars+nSecVars) )
            answer.at(i-nSecVars*nSecNodes) = i;
    }
 }

void
GradDpElement :: setNonlocalLocationArray(IntArray &answer, int nPrimNodes, int nPrimVars, int nSecNodes, int nSecVars)
{
    answer.resize(nlSize);
    for( int i = 1; i <= nlSize; i++ ) {
        answer.at(i) = i*nPrimVars+i;
    }
}


void
GradDpElement :: computeDisplacementDegreesOfFreedom(FloatArray &answer, TimeStep *stepN)
{
    StructuralElement* elem = this->giveStructuralElement();
    FloatArray u;
    answer.resize(locSize);
    answer.zero();
    
    elem -> computeVectorOf(EID_MomentumBalance, VM_Total, stepN, u);
    u.resize(totalSize);
    for ( int i = 1; i <= locSize; i++) {
        answer.at(i) = u.at(locU.at(i));
    }
}

void GradDpElement :: computeNonlocalDegreesOfFreedom(FloatArray &answer, TimeStep *stepN)
{
    StructuralElement* elem = this->giveStructuralElement();
    FloatArray u;
    answer.resize(nlSize);
    answer.zero();
    
    elem->computeVectorOf(EID_MomentumBalance, VM_Total, stepN, u);
    u.resize(totalSize);
    for  (int i = 1; i <= nlSize; i++ ) {
        answer.at(i) = u.at(locK.at(i));
    }
}

void
GradDpElement :: computeStressVector(FloatArray &answer, GaussPoint *gp, TimeStep *stepN)
{
    StructuralElement* elem = this->giveStructuralElement();
    FloatArray Epsilon;
    StructuralCrossSection *cs = static_cast< StructuralCrossSection * >( elem->giveCrossSection() );

    this->computeStrainVector(Epsilon, gp, stepN);
    cs->giveRealStresses(answer, ReducedForm, gp, Epsilon, stepN);
    int size = answer.giveSize()-1;
    answer.resize(size);
}


void
GradDpElement :: computeStrainVector(FloatArray &answer, GaussPoint *gp, TimeStep *stepN)
{
    FloatArray strain;
    double nlKappa;
    this->computeLocalStrainVector(strain,gp,stepN);
    this->computeNonlocalCumPlasticStrain(nlKappa,gp,stepN);

    answer = strain;
    int size = answer.giveSize();
    answer.resize(size+1);
    answer.at(size+1) = nlKappa;
}


void
GradDpElement :: computeLocalStrainVector(FloatArray &answer, GaussPoint *gp, TimeStep *stepN)
{
    int n;
    FloatMatrix b,A;
    FloatArray u,help;


    NLStructuralElement* elem = this->giveNLStructuralElement();
    nlGeo = elem->giveGeometryMode();
    
    this-> computeDisplacementDegreesOfFreedom(u,stepN);
    if ( nlGeo < 2 ) {
        elem->computeBmatrixAt(gp, b); 
        answer.beProductOf(b, u);
        // Green-Lagrange strain tensor (in vector form)
        // loop over all components of strain vector
        if ( nlGeo == 1 ) {
            n = answer.giveSize();
            for ( int i = 1; i <= n; i++ ) {
                // nonlinear part of the strain-displacement relation
                elem->computeNLBMatrixAt(A, gp, i);
                if ( A.isNotEmpty() ) {
                    help.beProductOf(A, u);
                    answer.at(i) += 0.5 * u.dotProduct(help);
                }
            }
        }
    }
    else if ( nlGeo == 2 ) {
        // deformation gradient will be used instead of strain
        elem->computeBFmatrixAt(gp, b);
        answer.beProductOf(b, u); // this gives the displacement gradient
        // unit matrix needs to be added
        // (needs to be adjusted if the mode is not 3d)
        answer.at(1) += 1.;
        answer.at(5) += 1.;
        answer.at(9) += 1.;
    }
}


void
GradDpElement :: computeNonlocalCumPlasticStrain(double &answer, GaussPoint *gp, TimeStep *stepN)
{
    FloatMatrix Nk;
    FloatArray u;
    FloatArray aux;

    this->computeNkappaMatrixAt(gp, Nk);
    this-> computeNonlocalDegreesOfFreedom(u,stepN); 
    aux.beProductOf(Nk, u);
    answer = aux.at(1);
}


void
GradDpElement :: computeNonlocalGradient(FloatArray &answer, GaussPoint *gp, TimeStep *stepN)
{
    FloatMatrix Bk;
    FloatArray u;
    FloatArray aux;

    this->computeBkappaMatrixAt(gp, Bk);
    this-> computeNonlocalDegreesOfFreedom(u,stepN); 
    answer.beProductOf(Bk, u);
}
    

void
GradDpElement :: giveNonlocalInternalForcesVector(FloatArray &answer,
                                              TimeStep *tStep, int useUpdatedGpRecord)
{
    //set displacement and nonlocal location array
    this-> setDisplacementLocationArray(locU,nPrimNodes, nPrimVars, nSecNodes, nSecVars);
    this-> setNonlocalLocationArray(locK,nPrimNodes, nPrimVars, nSecNodes, nSecVars);
    
    
    NLStructuralElement* elem = this->giveNLStructuralElement();
    // ?????????????
    MatResponseMode rMode = TangentStiffness;
    // ?????????????
    double tempKappa,dV;
    FloatMatrix stiffKappa, Nk;
    int size = nSecVars*nSecNodes;
    FloatArray fKappa(nlSize),aux(nlSize),dKappa,stress;
    aux.zero();
    GaussPoint *gp;
    Material *mat = elem->giveMaterial();
    IntegrationRule *iRule =  elem->giveIntegrationRule(0); 
    answer.resize(size);
    answer.zero();
    for ( int i = 0; i < iRule->giveNumberOfIntegrationPoints(); i++ ) {
        gp = iRule->getIntegrationPoint(i);
        this->computeNkappaMatrixAt(gp, Nk);
        for(int j = 1; j<=nlSize;j++)
        fKappa.at(j) = Nk.at(1,j);
        stress = static_cast< StructuralMaterialStatus * >( mat->giveStatus(gp) )->giveTempStressVector();    
        int size = stress.giveSize();
        tempKappa = stress.at(size);
        dV  = elem->computeVolumeAround(gp);
        fKappa.times(tempKappa);
        fKappa.times(-dV);
        aux.add(fKappa);
   }

    this->computeStiffnessMatrix_kk(stiffKappa,rMode,tStep);
    this-> computeNonlocalDegreesOfFreedom(dKappa,tStep);
    answer.beProductOf(stiffKappa,dKappa);
    answer.add(aux);
}


void
GradDpElement :: giveLocalInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord)
{
    NLStructuralElement* elem = this->giveNLStructuralElement();
    GaussPoint *gp;  
    IntegrationRule *iRule = elem->giveIntegrationRule(0);
    nlGeo = elem->giveGeometryMode();
    FloatMatrix b, bt, R, GNT,A;
    FloatArray bs, TotalStressVector;
    double dV;

    FloatMatrix  *ut = NULL;
    FloatArray u;

    answer.resize(0);
    if ( nlGeo ) {
        computeDisplacementDegreesOfFreedom(u,tStep);
        //      this->computeVectorOf(EID_MomentumBalance, VM_Total, tStep, u);
        if ( u.giveSize() ) {
            ut = new FloatMatrix( &u, 1);
            //delete u;
        }
        else {
            ut = NULL;
        }
    }
  
    for ( int i = 0; i < iRule->giveNumberOfIntegrationPoints(); i++ ) {
        gp = iRule->getIntegrationPoint(i);
        elem->computeBmatrixAt(gp, b);
        /////////////////////////////////
        if ( nlGeo ) {
            for ( int j = 1; j <= b.giveNumberOfRows(); j++ ) {
                // loop over each component of strain vector
                elem->computeNLBMatrixAt(A, gp, j);
                if (( A.isNotEmpty() ) && ( ut != NULL ) ) {
                    FloatMatrix b2;
                    b2.beProductOf(*ut,A);
                    for ( int k = 1; k <= b.giveNumberOfColumns(); k++ ) {
                        // add nonlinear contribution to each component
                        b.at(j, k) += b2.at(1, k); //mj
                    }
                }
            }
        } // end nlGeometry
        //////////////////////////////////////////////

        bt.beTranspositionOf(b);
        this->computeStressVector(TotalStressVector, gp, tStep);
        if ( TotalStressVector.giveSize() == 0 ) {
            break;
        }
    
        //
        // now every Gauss point has real stress vector
        //
        // compute nodal representation of internal forces using f = B^T*Sigma dV
        //
        dV  = elem->computeVolumeAround(gp);
        bs.beProductOf(bt, TotalStressVector);
        bs.times(dV);
        
        answer.add(bs); 
    }
    if ( nlGeo ) {
        delete ut;
    }
}

void
GradDpElement :: giveInternalForcesVector(FloatArray &answer,TimeStep *tStep, int useUpdatedGpRecord)
{
    //set displacement and nonlocal location array
    //this-> setDisplacementLocationArray(locU,nPrimNodes, nPrimVars, nSecNodes, nSecVars);
    //this-> setNonlocalLocationArray(locK,nPrimNodes, nPrimVars, nSecNodes, nSecVars);
    
    //StructuralElement* elem = this->giveStructuralElement();
    
    answer.resize(totalSize);
    answer.zero();
    FloatArray answerU;
    answerU.resize(locSize);
    answer.zero();
    FloatArray answerK(nlSize);
    answerK.zero();
    
    this->giveLocalInternalForcesVector(answerU,tStep, useUpdatedGpRecord);
    this->giveNonlocalInternalForcesVector(answerK,tStep, useUpdatedGpRecord);
    answer.assemble(answerU,locU);
    answer.assemble(answerK,locK);
}

void
GradDpElement :: computeForceLoadVector(FloatArray &answer, TimeStep *stepN, ValueModeType mode)
{
    //StructuralElement* elem = this->giveStructuralElement();
    
    //set displacement and nonlocal location array
    this-> setDisplacementLocationArray(locU,nPrimNodes, nPrimVars, nSecNodes, nSecVars);
    this-> setNonlocalLocationArray(locK,nPrimNodes, nPrimVars, nSecNodes, nSecVars);
    

    FloatArray localForces(locSize);
    FloatArray nlForces(nlSize);
    answer.resize(totalSize);
    answer.zero();
    this->computeLocForceLoadVector(localForces,stepN,mode); 

    answer.assemble(localForces,locU);
    answer.assemble(nlForces,locK);
}

void
GradDpElement :: computeNonForceLoadVector(FloatArray &answer, TimeStep *stepN, ValueModeType mode)
{
    //set displacement and nonlocal location array
    this-> setDisplacementLocationArray(locU,nPrimNodes, nPrimVars, nSecNodes, nSecVars);
    this-> setNonlocalLocationArray(locK,nPrimNodes, nPrimVars, nSecNodes, nSecVars);

    FloatArray localForces(locSize);
    FloatArray nlForces(nlSize);
    answer.resize(totalSize);
    answer.zero();

    this->computeLocNonForceLoadVector(answer,stepN,mode); 
}


  /************************************************************************/
void
GradDpElement :: computeLocForceLoadVector(FloatArray &answer, TimeStep *stepN, ValueModeType mode)
// computes the part of load vector, which is imposed by force loads acting
// on element volume (surface).
// When reactions forces are computed, they are computed from element::GiveRealStressVector
// in this vector a real forces are stored (temperature part is subtracted).
// so we need further sobstract part corresponding to non-nodeal loading.
{
    FloatMatrix T;
    NLStructuralElement* elem = this->giveNLStructuralElement();
    elem->computeLocalForceLoadVector(answer, stepN, mode);
    
    // transform result from global cs to nodal cs. if necessary
    if ( answer.isNotEmpty() ) {
        if ( elem->computeGtoLRotationMatrix(T) ) {
            // first back to global cs from element local
            answer.rotatedWith(T, 't');
        }
    }
    else {
        answer.resize(locSize);
        answer.zero();
    }
}


void
GradDpElement :: computeLocNonForceLoadVector(FloatArray &answer, TimeStep *stepN, ValueModeType mode)
// Computes the load vector of the receiver, at stepN.
{
    FloatArray helpLoadVector;
    StructuralElement* elem = this->giveStructuralElement();
    answer.resize(0);

    elem->computePrescribedStrainLoadVectorAt(helpLoadVector, stepN, mode);
    if ( helpLoadVector.giveSize() ) {
        answer.add(helpLoadVector);
    }
}


//tohle dela polylinenonlocalbarrier a budu to schovavat na gradientnim materialu
void GradDpElement::computeDistanceToBoundary()
{
    int nelem;
    StructuralElement* elem = this->giveStructuralElement();
    IntegrationRule *iRule = elem->giveIntegrationRule(0); 
    Domain *d = elem->giveMaterial()->giveDomain();
    nelem = d->giveNumberOfElements();
    
    for ( int i = 1; i <= nelem; i++ ) {
        for ( int j = 0 ; j < iRule->giveNumberOfIntegrationPoints(); j++ ) {
            GaussPoint *gp = iRule->getIntegrationPoint(j);
            FloatArray lcoord(2);
            FloatArray coord(2);
            lcoord(0) = gp->giveCoordinate(1);
            lcoord(1) = gp->giveCoordinate(2);
            elem ->computeGlobalCoordinates(coord,lcoord);
            //double distance = d->giveBoundary(1)->computeDistanceTo(coord.at(1),coord.at(2));
            //gp->setDistanceToBoundary(distance);
        }
    }
  
}

/*******************************************************************/

void
GradDpElement :: computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
{
    //set displacement and nonlocal location array
    this-> setDisplacementLocationArray(locU,nPrimNodes, nPrimVars, nSecNodes, nSecVars);
    this-> setNonlocalLocationArray(locK,nPrimNodes, nPrimVars, nSecNodes, nSecVars);
    
    
    answer.resize(totalSize,totalSize);
    answer.zero();


    if(averType == 1)
        this->computeDistanceToBoundary();
    

    FloatMatrix answer1,answer2,answer3,answer4;
    this->computeStiffnessMatrix_uu(answer1, rMode,tStep);
    this->computeStiffnessMatrix_uk(answer2, rMode,tStep);
    this->computeStiffnessMatrix_ku(answer3, rMode,tStep);
    this->computeStiffnessMatrix_kk(answer4, rMode,tStep);
    answer.assemble(answer1,locU);
    answer.assemble(answer2,locU,locK);
    answer.assemble(answer3,locK,locU);
    answer.assemble(answer4,locK);
}
 


void
GradDpElement :: computeStiffnessMatrix_uu(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
{
    double dV;
    GaussPoint *gp;
    NLStructuralElement* elem = this->giveNLStructuralElement();
    nlGeo = elem->giveGeometryMode();
    IntegrationRule *iRule = elem->giveIntegrationRule(0); 
    FloatMatrix A, *ut = NULL;
    FloatArray u;
    //MatResponseForm form = PDGrad_uu;
    bool matStiffSymmFlag = elem->giveCrossSection()->isCharacteristicMtrxSymmetric(rMode, elem->material);
    //////////////////////////////////////////////////// 
    
    if ( nlGeo ) {
        //    this->computeVectorOf(EID_MomentumBalance, VM_Total, tStep, u);
        computeDisplacementDegreesOfFreedom(u, tStep);
    
        if ( u.giveSize() ) {
            ut = new FloatMatrix( &u, 1);
            // delete u;
        }
        else {
            ut = NULL;
        }
    }
    //////////////////////////////


    FloatMatrix B,DB,D;
    answer.resize(locSize,locSize);
    answer.zero();
    if (!elem->isActivated(tStep)) return;
    
    Material *mat = elem->giveMaterial();
    for (int j = 0; j < iRule->giveNumberOfIntegrationPoints(); j++ ) {
        gp = iRule->getIntegrationPoint(j);
        elem->computeBmatrixAt(gp, B);
        if ( nlGeo ) {
            for ( int l = 1; l <=  B.giveNumberOfRows(); l++ ) {
                // loop over each component of strain vector
                elem->computeNLBMatrixAt(A, gp, l);
                if ( ( A.isNotEmpty() ) && ( ut != NULL ) ) {
                    FloatMatrix b2;
                    b2.beProductOf(* ut, A);
                    for ( int k = 1; k <= B.giveNumberOfColumns(); k++ ) {
                    // add nonlinear contribution to each component
                        B.at(l, k) += b2.at(1, k); //mj
                    }
                }
            }
        } // end nlGeometry
      
        mat-> giveCharacteristicMatrix(D,PDGrad_uu,rMode, gp, tStep);
        dV = elem->computeVolumeAround(gp);
        DB.beProductOf(D, B);
        if ( matStiffSymmFlag ) {
            answer.plusProductSymmUpper(B, DB, dV);
        } else {
            answer.plusProductUnsym(B, DB, dV);
        }
    }

    if ( nlGeo ) {
        delete ut;
    }

    if ( nlGeo ) {
        // assemble initial stress matrix
        for ( int i = 0; i < iRule->giveNumberOfIntegrationPoints(); i++ ) {
            gp = iRule->getIntegrationPoint(i);
            dV = elem->computeVolumeAround(gp);
            FloatArray stress = static_cast< StructuralMaterialStatus * >( mat->giveStatus(gp) )->giveTempStressVector();
            int n = 6;
            if ( n ) {
                for ( int j = 1; j <= n; j++ ) {
                    // loop over each component of strain vector
                    elem->computeNLBMatrixAt(A, gp, j);
                    if ( A.isNotEmpty() ) {
                        A.times(stress.at(j) * dV);
                        answer.add(A);
                    }
                }
            }
        }
    } // end nlGeometry

    if ( matStiffSymmFlag ) {
        answer.symmetrized();
    }
  
} 

void
GradDpElement :: computeStiffnessMatrix_ku(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
{
    NLStructuralElement* elem = this->giveNLStructuralElement();
    double dV;
    GaussPoint *gp;
    IntegrationRule *iRule = elem->giveIntegrationRule(0); 
    //MatResponseForm form = PDGrad_uu;
    ///////////////////////////////////////////////////////////////////
    FloatMatrix B,DB,D,Nk,NkT, NkDB;
    answer.resize(nlSize,locSize);
    answer.zero();
    Material *mat =elem->giveMaterial();
    // mat->getGradientFormulation(averType);
    FloatArray u;
    FloatMatrix A, *ut = NULL;
    if ( nlGeo ) {
        computeDisplacementDegreesOfFreedom(u, tStep);
        //this->computeVectorOf(EID_MomentumBalance, VM_Total, tStep, u);
      
        if ( u.giveSize() ) {
            ut = new FloatMatrix( &u, 1);
            // delete u;
        } else {
            ut = NULL;
        }
    }
 

    for ( int j = 0; j < iRule->giveNumberOfIntegrationPoints(); j++ ) {
        gp = iRule->getIntegrationPoint(j);
        elem->computeBmatrixAt(gp, B);
        if ( nlGeo ) {
            for ( int l = 1; l <=  B.giveNumberOfRows(); l++ ) {
                // loop over each component of strain vector
                elem->computeNLBMatrixAt(A, gp, l);
                if ( ( A.isNotEmpty() ) && ( ut != NULL ) ) {
                    FloatMatrix b2;
                    b2.beProductOf(* ut, A);
                    for ( int k = 1; k <= B.giveNumberOfColumns(); k++ ) {
                        // add nonlinear contribution to each component
                        B.at(l, k) += b2.at(1, k); //mj
                    }
                }
            }
        } // end nlGeometry
        mat->giveCharacteristicMatrix(D,PDGrad_ku,rMode, gp, tStep);
        this->computeNkappaMatrixAt(gp,Nk);
        NkT.beTranspositionOf(Nk);
        dV = elem->computeVolumeAround(gp);
        DB.beProductOf(D, B);
        NkDB.beProductOf(NkT,DB);
        NkDB.times(-dV);  
        answer.add(NkDB);

        if ( averType == 2 ) {
            FloatMatrix BtL;
            FloatMatrix lStiffDerivative1, dL1, dL2, B,Bk,Nk, D;
            FloatMatrix N1(2,2);
            FloatMatrix N2(2,2);
            FloatArray Gk,bl;
            FloatArray n1(3);
            FloatArray n2(3);
            elem ->computeBmatrixAt(gp,B);
            this ->computeBkappaMatrixAt(gp,Bk);

            mat-> giveCharacteristicMatrix(lStiffDerivative1,PDGrad_LD,rMode, gp, tStep);
            mat-> giveCharacteristicMatrix(D,PDGrad_uu,rMode, gp, tStep);
            this -> computeNonlocalGradient(Gk, gp,tStep);
            N1.at(1,1) = lStiffDerivative1.at(1,1)*lStiffDerivative1.at(1,1);
            N1.at(1,2) = lStiffDerivative1.at(2,1)*lStiffDerivative1.at(1,1);
            N1.at(2,1) = N1.at(1,2);
            N1.at(2,2) = lStiffDerivative1.at(2,1)*lStiffDerivative1.at(2,1);
            
            N2.at(1,1) = lStiffDerivative1.at(1,2)*lStiffDerivative1.at(1,2);
            N2.at(1,2) = lStiffDerivative1.at(2,2)*lStiffDerivative1.at(1,2);
            N2.at(2,1) = N2.at(1,2);
            N2.at(2,2) = lStiffDerivative1.at(2,2)*lStiffDerivative1.at(2,2);
            
            n1.at(1) = N1.at(1,1);
            n1.at(2) = N1.at(2,2);
            n1.at(3) = N1.at(1,2);
            
            n2.at(1) = N2.at(1,1);
            n2.at(2) = N2.at(2,2);
            n2.at(3) = N2.at(1,2);
            
            dL1 = N2;
            dL1.times(lStiffDerivative1.at(3,3));
            
            dL2 = N2;
            dL2.times(lStiffDerivative1.at(4,4));
            int nC = Bk.giveNumberOfColumns();
            int nR = B.giveNumberOfColumns();
            FloatMatrix result(nC,nR);
            result.zero();

            for ( int i = 1; i <= 3; i++ ) {
                for ( int q = 1; j <= 12; j++ ) {
                    for ( int m = 1; m <= 2; m++ ) {
                        for(int n = 1; n <= 2; n++ ) {
                            for(int p = 1; p <= 3; p++ ) {
                                for(int o = 1; o <= 3; o++ ) {
                                    result.at(i,q) += Bk.at(m,i)*(dL1.at(m,n)*n1.at(p)+dL2.at(m,n)*n2.at(p))*D.at(p,o)*B.at(o,q)*Gk.at(n); 
                                }
                            }
                        }
                    }
                }
            }
            result.times(dV);
            answer.add(result);
        }
    }
    if ( nlGeo ) {
      delete ut;
    }
}


void 
GradDpElement :: computeStiffnessMatrix_kk(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
{
    StructuralElement* elem = this->giveStructuralElement();
    double dV;
    double l;
    GaussPoint *gp;
    IntegrationRule *iRule = elem->giveIntegrationRule(0); 
    FloatMatrix lStiff;
    ///////////////////////////////////////////////////////////////////
    FloatMatrix Bk,Bt,BtB,N,Nt,NtN;

    Material *mat = elem->giveMaterial();
    answer.resize(nlSize,nlSize);
    answer.zero();
    //  mat->getGradientFormulation(averType);

    for ( int j = 0; j < iRule->giveNumberOfIntegrationPoints(); j++ ) {
        gp = iRule->getIntegrationPoint(j);
        this->computeNkappaMatrixAt(gp,N);
        Nt.beTranspositionOf(N);
        this->computeBkappaMatrixAt(gp,Bk);
        Bt.beTranspositionOf(Bk);
        dV = elem->computeVolumeAround(gp);
        mat-> giveCharacteristicMatrix(lStiff,PDGrad_kk,rMode, gp, tStep);
        NtN.beProductOf(Nt,N);
        NtN.times(dV);
        answer.add(NtN);
        if ( averType == 0 || averType == 1 ) {
            l = lStiff.at(1,1);
            BtB.beProductOf(Bt,Bk);
            BtB.times(l*l*dV);
            answer.add(BtB);
        }
        else if ( averType == 2 ) {
            FloatMatrix BtL;
            FloatMatrix lStiffDerivative1, dL1, dL2, B, D;
            FloatMatrix N1(2,2);
            FloatMatrix N2(2,2);
            FloatArray Gk,b;
            FloatArray n1(3);
            FloatArray n2(3);
            FloatMatrix Nk;
            this -> computeNkappaMatrixAt(gp,Nk);
            BtL.beProductOf(Bt,lStiff);
            BtB.beProductOf(BtL,Bk);
            BtB.times(dV);
            elem->computeBmatrixAt(gp,B);
            //mat-> giveCharacteristicMatrix(lStiffDerivative1,PDGrad_LD,rMode, gp, tStep);
            mat-> giveCharacteristicMatrix(D,PDGrad_ku,rMode, gp, tStep);
            this -> computeNonlocalGradient(Gk, gp,tStep);
            N1.at(1,1) = lStiffDerivative1.at(1,1)*lStiffDerivative1.at(1,1);
            N1.at(1,2) = lStiffDerivative1.at(2,1)*lStiffDerivative1.at(1,1);
            N1.at(2,1) = N1.at(1,2);
            N1.at(2,2) = lStiffDerivative1.at(2,1)*lStiffDerivative1.at(2,1);

            N2.at(1,1) = lStiffDerivative1.at(1,2)*lStiffDerivative1.at(1,2);
            N2.at(1,2) = lStiffDerivative1.at(2,2)*lStiffDerivative1.at(1,2);
            N2.at(2,1) = N2.at(1,2);
            N2.at(2,2) = lStiffDerivative1.at(2,2)*lStiffDerivative1.at(2,2);

            n1.at(1) = N1.at(1,1);
            n1.at(2) = N1.at(2,2);
            n1.at(3) = N1.at(1,2);

            n2.at(1) = N2.at(1,1);
            n2.at(2) = N2.at(2,2);
            n2.at(3) = N2.at(1,2);

            dL1 = N2;
            dL1.times(lStiffDerivative1.at(3,3));

            dL2 = N2;
            dL2.times(lStiffDerivative1.at(4,4));

            int nC = Bk.giveNumberOfColumns();
            int nR = Nk.giveNumberOfColumns();
            FloatMatrix result(nC,nR);

            for ( int i = 1; i <= 3; i++ ) {
                for ( int q = 1; j <= 3; j++ ) {
                    for ( int m = 1; m <= 2; m++ ) {
                        for ( int n = 1; n <= 2; n++ ) {
                            for ( int p = 1; p <= 3; p++ ) {
                                for ( int o = 1; o <= 3; o++ ) {
                                    result.at(i,q) += Bk.at(m,i)*(dL1.at(m,n)*n1.at(p)+dL2.at(m,n)*n2.at(p))*D.at(1,p)*Nk.at(1,q)*Gk.at(n); 
                                }
                            }
                        }
                    }
                }
            }

            result.times(dV);
            answer.add(BtB);
            answer.add(result);
        }
  } 
} 


void 
GradDpElement :: computeStiffnessMatrix_uk(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
{
    NLStructuralElement* elem = this->giveNLStructuralElement();
    double dV;
    GaussPoint *gp;
    Material *mat = elem->giveMaterial();
    IntegrationRule *iRule = elem->giveIntegrationRule(0); 
    ///////////////////////////////////////////////////////////////////
    FloatMatrix B,Bt,Nk,BSN;
    FloatMatrix BS;
    FloatMatrix gPSigma;
    FloatMatrix A, *ut = NULL;
    FloatArray u;
    answer.resize(locSize,nlSize);
    answer.zero();  
    
    if ( nlGeo ) {
        computeDisplacementDegreesOfFreedom(u, tStep);
        //this->computeVectorOf(EID_MomentumBalance, VM_Total, tStep, u);
      
        if ( u.giveSize() ) {
            ut = new FloatMatrix( &u, 1);
            // delete u;
        }
        else {
            ut = NULL;
        }
    }
 
    for ( int j = 0; j < iRule->giveNumberOfIntegrationPoints(); j++ ) {
        gp = iRule->getIntegrationPoint(j);
        mat-> giveCharacteristicMatrix(gPSigma,PDGrad_uk,rMode, gp, tStep);
        this->computeNkappaMatrixAt(gp,Nk);
        elem->computeBmatrixAt(gp,B);
        if ( nlGeo ) {
            for ( int l = 1; l <=  B.giveNumberOfRows(); l++ ) {
                // loop over each component of strain vector
                elem->computeNLBMatrixAt(A, gp, l);
                if ( ( A.isNotEmpty() ) && ( ut != NULL ) ) {
                    FloatMatrix b2;
                    b2.beProductOf(* ut, A);
                    
                    for ( int k = 1; k <= B.giveNumberOfColumns(); k++ ) {
                        // add nonlinear contribution to each component
                        B.at(l, k) += b2.at(1, k); //mj
                    }
                }
            }
        } // end nlGeometry
        Bt.beTranspositionOf(B);
        dV = elem->computeVolumeAround(gp);
        BS.beProductOf(Bt,gPSigma);
        BSN.beProductOf(BS,Nk);
        BSN.times(-dV);
        answer.add(BSN); 
    }
    if ( nlGeo ) {
      delete ut;
    }
}

IRResultType
GradDpElement :: initializeFrom(InputRecord *ir)
{
    //const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    //IRResultType result;                // Required by IR_GIVE_FIELD macro
    //nlGeo = 0;

    return IRRT_OK;
}

} // end namespace oofem

