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

#include "graddpelement.h"
#include "node.h"
#include "material.h"
#include "gausspnt.h"
#include "gaussintegrationrule.h"
#include "flotmtrx.h"
#include "flotarry.h"
#include "intarray.h"
#include "domain.h"
#include "cltypes.h"
#include "structuralms.h"
#include "mathfem.h"
#include "structuralcrosssection.h"
#include "structuralelement.h"

#ifndef __MAKEDEPEND
 #include <stdio.h>
#endif

namespace oofem {

GradDpElement :: GradDpElement ()
// Constructor.
{

}

void
GradDpElement :: setDisplacementLocationArray(IntArray &answer, int nPrimNodes, int nPrimVars, int nSecNodes, int nSecVars)
{

    answer.resize(locSize);

    for(int i=1; i<=totalSize; i++) {
        if ( i<nSecNodes*nPrimVars+1 )
            answer.at(i) = i +(int)(((i-1)/nPrimVars))*nSecVars;
        else if ( i>nSecNodes*(nPrimVars+nSecVars) )
            answer.at(i-nSecVars*nSecNodes) = i;
    }

}

void
GradDpElement :: setNonlocalLocationArray(IntArray &answer, int nPrimNodes, int nPrimVars, int nSecNodes, int nSecVars)
{
    answer.resize(nlSize);
    for (int i =1; i<=nlSize; i++) {
        answer.at(i) = i*nPrimVars+i;
    }
}


void
GradDpElement :: computeDisplacementDegreesOfFreedom(FloatArray &answer, GaussPoint *gp, TimeStep *stepN)
{

    StructuralElement* elem = this->giveStructuralElement();
    FloatArray u;
    answer.resize(locSize);
    answer.zero();

    elem -> computeVectorOf(EID_MomentumBalance, VM_Total, stepN, u);
    u.resize(totalSize);
    for (int i=1; i<=locSize; i++) {
        answer.at(i) = u.at(locU.at(i));
    }
}


void GradDpElement :: computeNonlocalDegreesOfFreedom(FloatArray &answer, GaussPoint *gp, TimeStep *stepN)
{
    StructuralElement* elem = this->giveStructuralElement();
    FloatArray u;
    answer.resize(nlSize);
    answer.zero();

    elem -> computeVectorOf(EID_MomentumBalance, VM_Total, stepN, u);
    u.resize(totalSize);
    for (int i=1; i<=nlSize; i++) {
        answer.at(i) = u.at(locK.at(i));
    }
}


void
GradDpElement :: computeStressVector(FloatArray &answer, GaussPoint *gp, TimeStep *stepN)
{
    StructuralElement* elem = this->giveStructuralElement();
    FloatArray Epsilon;
    StructuralCrossSection *cs = ( StructuralCrossSection * ) elem->giveCrossSection();

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
   FloatMatrix b;
   FloatArray u;
   StructuralElement* elem = this->giveStructuralElement();

#if 0
    if ( !this->isActivated(stepN) ) {
        answer.resize( this->giveCrossSection()->giveIPValueSize(IST_StrainTensor, gp) );
        answer.zero();
        return;
    }
#endif

    elem->computeBmatrixAt(gp, b);
    this->computeDisplacementDegreesOfFreedom(u, gp,stepN);
    answer.beProductOf(b, u);
}


void
GradDpElement :: computeNonlocalCumPlasticStrain(double &answer, GaussPoint *gp, TimeStep *stepN)
{
    FloatMatrix Nk;
    FloatArray u;
    FloatArray aux;

    /*if ( !this->isActivated(stepN) ) {
        answer = 0;
            return;
    }*/
    this->computeNkappaMatrixAt(gp, Nk);
    this->computeNonlocalDegreesOfFreedom(u, gp,stepN);
    aux.beProductOf(Nk, u);
    answer = aux.at(1);
}


void
GradDpElement :: giveNonlocalInternalForcesVector(FloatArray &answer,
                                              TimeStep *tStep, int useUpdatedGpRecord)
{
    //set displacement and nonlocal location array
    this-> setDisplacementLocationArray(locU,nPrimNodes, nPrimVars, nSecNodes, nSecVars);
    this-> setNonlocalLocationArray(locK,nPrimNodes, nPrimVars, nSecNodes, nSecVars);

    StructuralElement* elem = this->giveStructuralElement();
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
    for ( int i = 0; i < iRule->getNumberOfIntegrationPoints(); i++ ) {
        gp = iRule->getIntegrationPoint(i);
        this->computeNkappaMatrixAt(gp, Nk);
        for(int j = 1; j<=nlSize;j++)
        fKappa.at(j) = Nk.at(1,j);
        stress = ( (StructuralMaterialStatus *)mat->giveStatus(gp) )->giveTempStressVector();
        int size = stress.giveSize();
        tempKappa = stress.at(size);
        dV  = elem->computeVolumeAround(gp);
        fKappa.times(tempKappa);
        fKappa.times(-dV);
        aux.add(fKappa);
    }

    this->computeStiffnessMatrix_kk(stiffKappa,rMode,tStep);
    this-> computeNonlocalDegreesOfFreedom(dKappa,gp,tStep);
    answer.beProductOf(stiffKappa,dKappa);
    answer.add(aux);
}


void
GradDpElement :: giveLocalInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord)
{
    StructuralElement* elem = this->giveStructuralElement();
    GaussPoint *gp;
    IntegrationRule *iRule = elem->giveIntegrationRule(0);

    FloatMatrix b, bt, R, GNT;
    FloatArray bs, TotalStressVector;
    double dV;
    answer.resize(0);

    for ( int i = 0; i < iRule->getNumberOfIntegrationPoints(); i++ ) {
        gp = iRule->getIntegrationPoint(i);
        elem->computeBmatrixAt(gp, b);
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
    this->setDisplacementLocationArray(locU,nPrimNodes, nPrimVars, nSecNodes, nSecVars);
    this->setNonlocalLocationArray(locK,nPrimNodes, nPrimVars, nSecNodes, nSecVars);


    FloatArray localForces(locSize);
    FloatArray nlForces(nlSize);
    answer.resize(totalSize);
    answer.zero();
    this->computeLocForceLoadVector(localForces,stepN,mode);

    answer.assemble(localForces,locU);
    answer.assemble(nlForces,locK);
}


/************************************************************************/
void
GradDpElement :: computeLocForceLoadVector(FloatArray &answer, TimeStep *stepN, ValueModeType mode)
// computes the part of load vector, which is imposed by force loads acting
// on element volume (surface).
// Why is this function taken separately ?
// When reactions forces are computed, they are computed from element::GiveRealStressVector
// in this vector a real forces are stored (temperature part is subtracted).
// so we need further sobstract part corresponding to non-nodeal loading.
{
    FloatMatrix T;
    StructuralElement* elem = this->giveStructuralElement();
    elem->computeLocalForceLoadVector(answer, stepN, mode);

    if ( answer.isEmpty() ) {
        answer.resize(locSize);
        answer.zero();
    }
}


void
GradDpElement :: computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
{

    //set displacement and nonlocal location array
    this-> setDisplacementLocationArray(locU,nPrimNodes, nPrimVars, nSecNodes, nSecVars);
    this-> setNonlocalLocationArray(locK,nPrimNodes, nPrimVars, nSecNodes, nSecVars);


    answer.resize(totalSize,totalSize);
    answer.zero();

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
    StructuralElement* elem = this->giveStructuralElement();
    IntegrationRule *iRule = elem->giveIntegrationRule(0);
    //MatResponseForm form = PDGrad_uu;
    bool matStiffSymmFlag = elem->giveCrossSection()->isCharacteristicMtrxSymmetric(rMode, elem->material);
    ///////////////////////////////////////////////////////////////////
    FloatMatrix B,DB,D;
    answer.resize(locSize,locSize);
    answer.zero();
    if (!elem->isActivated(tStep)) return;

    Material *mat = elem->giveMaterial();
    for (int j = 0; j < iRule->getNumberOfIntegrationPoints(); j++ ) {
        gp = iRule->getIntegrationPoint(j);
        elem->computeBmatrixAt(gp, B);
        mat-> giveCharacteristicMatrix(D,PDGrad_uu,rMode, gp, tStep);
        dV = elem->computeVolumeAround(gp);
        DB.beProductOf(D, B);
        if ( matStiffSymmFlag ) {
            answer.plusProductSymmUpper(B, DB, dV);
        } else {
            answer.plusProductUnsym(B, DB, dV);
        }
    }
    if ( matStiffSymmFlag ) {
        answer.symmetrized();
    }
}


void
GradDpElement :: computeStiffnessMatrix_ku(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
{
    StructuralElement* elem = this->giveStructuralElement();
    double dV;
    GaussPoint *gp;
    IntegrationRule *iRule = elem->giveIntegrationRule(0);
    //MatResponseForm form = PDGrad_uu;
    ///////////////////////////////////////////////////////////////////
    FloatMatrix B,DB,D,Nk,NkT, NkDB;
    answer.resize(nlSize,locSize);
    answer.zero();
    Material *mat =elem->giveMaterial();
    for (int j = 0; j < iRule->getNumberOfIntegrationPoints(); j++ ) {
        gp = iRule->getIntegrationPoint(j);
        elem->computeBmatrixAt(gp, B);
        mat->giveCharacteristicMatrix(D,PDGrad_ku,rMode, gp, tStep);
        this->computeNkappaMatrixAt(gp,Nk);
        NkT.beTranspositionOf(Nk);
        dV = elem->computeVolumeAround(gp);
        DB.beProductOf(D, B);
        NkDB.beProductOf(NkT,DB);
        NkDB.times(-dV);
        answer.add(NkDB);
    }
}


void
GradDpElement :: computeStiffnessMatrix_kk(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
{
    StructuralElement* elem = this->giveStructuralElement();
    double dV;
    double R;
    GaussPoint *gp;
    IntegrationRule *iRule = elem->giveIntegrationRule(0);
    FloatMatrix lStiff;
    ///////////////////////////////////////////////////////////////////
    FloatMatrix B,Bt,BtB,N,Nt,NtN;
    Material *mat = elem->giveMaterial();
    answer.resize(nlSize,nlSize);
    answer.zero();

    for (int j = 0; j < iRule->getNumberOfIntegrationPoints(); j++ ) {
        gp = iRule->getIntegrationPoint(j);
        this->computeNkappaMatrixAt(gp,N);
        Nt.beTranspositionOf(N);
        this->computeBkappaMatrixAt(gp,B);
        Bt.beTranspositionOf(B);
        dV = elem->computeVolumeAround(gp);
        mat-> giveCharacteristicMatrix(lStiff,PDGrad_kk,rMode, gp, tStep);
        R = lStiff.at(1,1);
        NtN.beProductOf(Nt,N);
        NtN.times(dV);
        BtB.beProductOf(Bt,B);
        BtB.times(R*dV);
        answer.add(NtN);
        answer.add(BtB);

    }
}


void
GradDpElement :: computeStiffnessMatrix_uk(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
{
    StructuralElement* elem = this->giveStructuralElement();
    double dV;
    GaussPoint *gp;
    Material *mat = elem->giveMaterial();
    IntegrationRule *iRule = elem->giveIntegrationRule(0);
    ///////////////////////////////////////////////////////////////////
    FloatMatrix B,Bt,Nk,BSN;
    FloatMatrix BS;
    FloatMatrix gPSigma;
    answer.resize(locSize,nlSize);
    answer.zero();
    for (int j = 0; j < iRule->getNumberOfIntegrationPoints(); j++ ) {
        gp = iRule->getIntegrationPoint(j);
        mat-> giveCharacteristicMatrix(gPSigma,PDGrad_uk,rMode, gp, tStep);
        this->computeNkappaMatrixAt(gp,Nk);
        elem->computeBmatrixAt(gp,B);
        Bt.beTranspositionOf(B);
        dV = elem->computeVolumeAround(gp);
        BS.beProductOf(Bt,gPSigma);
        BSN.beProductOf(BS,Nk);
        BSN.times(-dV);
        answer.add(BSN);
    }
}

} // end namespace oofem
