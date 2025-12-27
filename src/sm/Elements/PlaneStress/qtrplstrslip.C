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
 *               Copyright (C) 1993 - 2025   Borek Patzak
 *
 *
 *
 *       Czech Technical University, Faculty of Civil Engineering,
 *   Department of Structural Mechanics, 166 29 Prague, Czech Republic
 *
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with this library; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */

#include "sm/Elements/PlaneStress/qtrplstrslip.h"
#include "fei2dtrquad.h"
#include "node.h"
#include "gausspoint.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "intarray.h"
#include "crosssection.h"
#include "gaussintegrationrule.h"
#include "mathfem.h"
#include "classfactory.h"
#include "Materials/structuralmaterial.h"
#include "Materials/structuralms.h"
#include "Materials/structuralslipfe2material.h"
#include "CrossSections/structuralcrosssection.h"


namespace oofem {
REGISTER_Element(QTrPlaneStress2dSlip);

QTrPlaneStress2dSlip :: QTrPlaneStress2dSlip(int n, Domain *aDomain) :
    QTrPlaneStress2d(n, aDomain)
{
    numberOfGaussPoints = 4;
}


void QTrPlaneStress2dSlip::giveDofManDofIDMask( int inode, IntArray &answer ) const
{
    answer = {
        D_u, D_v, S_u, S_v
    };
}


void QTrPlaneStress2dSlip::computeStiffnessMatrix( FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep )
{
    int ndof=this->computeNumberOfDofs();

    FloatMatrix Kaa, Kab, Kba, Kbb, daa, dab1, dab2, dba1, dba2, dbb1, dbb2, dbb3, dbb4;

    for ( auto &gp : *this->giveDefaultIntegrationRulePtr() ) {
        FloatMatrix Ba, Bb, Nb;
        double dV;
        this->computeBmatrixAt(gp, Ba);
        this->computeBHmatrixAt(gp, Bb);
        this->computeNmatrixAt(gp->giveNaturalCoordinates(), Nb);
        dV = this->computeVolumeAround(gp);

        //Get the sensitivities
        FloatMatrix dStressdEps, dStressdS, dStressdG, dBStressdEps, dBStressdS, dBStressdG, dRStressdEps, dRStressdS, dRStressdG;
        this->giveSensitivities(dStressdEps, dStressdS, dStressdG, dBStressdEps, dBStressdS, dBStressdG, dRStressdEps, dRStressdS, dRStressdG, rMode, gp, tStep );

        daa.beProductOf(dStressdEps, Ba);
        dab1.beProductOf(dStressdS,Nb);
        dab2.beProductOf(dStressdG,Bb);
        dba1.beProductOf(dBStressdEps,Ba);
        dba2.beProductOf(dRStressdEps,Ba);
        dbb1.beProductOf(dBStressdS,Nb);
        dbb2.beProductOf(dBStressdG,Bb);
        dbb3.beProductOf(dRStressdS,Nb);
        dbb4.beProductOf(dRStressdG,Bb);

        Kaa.plusProductUnsym(Ba, daa, dV);
        Kab.plusProductUnsym(Ba, dab1, dV);
        Kab.plusProductUnsym(Ba, dab2, dV);
        Kba.plusProductUnsym(Nb,dba1,dV);
        Kba.plusProductUnsym(Bb,dba2,dV);
        Kbb.plusProductUnsym(Nb,dbb1,dV);
        Kbb.plusProductUnsym(Nb,dbb2,dV);
        Kbb.plusProductUnsym(Bb,dbb3,dV);
        Kbb.plusProductUnsym(Bb,dbb4,dV);
    }

    answer.resize(ndof,ndof);
    answer.assemble(Kaa,aMask,aMask);
    answer.assemble(Kbb,bMask,bMask);
    answer.assemble(Kab,aMask,bMask);
    answer.assemble(Kba,bMask,aMask);

}


void QTrPlaneStress2dSlip::giveInternalForcesVector( FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord )
{
    int ndof=this->computeNumberOfDofs();

    //Compute local fields at nodes
    FloatArray u, a, b;
    this->computeVectorOf(VM_Total, tStep, u);

    a.beSubArrayOf(u,aMask);
    b.beSubArrayOf(u,bMask);

    FloatArray vStrain, vSlip, vSlipGradient;
    FloatArray Stress, bStress, rStress;
    FloatArray finta, fintb;

    for ( auto &gp : *this->giveDefaultIntegrationRulePtr() ) {
        FloatMatrix Ba, Bb, Nb;
        FloatArray f;
        double dV;
        this->computeBmatrixAt(gp, Ba);
        this->computeBHmatrixAt(gp, Bb);
        this->computeNmatrixAt(gp->giveNaturalCoordinates(), Nb);
        dV = this->computeVolumeAround(gp);

        vStrain.beProductOf(Ba, a);
        vSlip.beProductOf(Nb, b);
        vSlipGradient.beProductOf(Bb, b);

        this->giveHomogenizedFields(Stress, bStress, rStress, vStrain, vSlip, vSlipGradient, gp , tStep);

        finta.plusProduct(Ba,Stress,dV);
        fintb.plusProduct(Nb,bStress,dV);
        fintb.plusProduct(Bb,rStress,dV);
    }

    answer.resize(ndof);
    for (int i=1; i <= ndof/2 ; i++) {
        if (i % 2 == 1)  {
            answer.at(2*i-1) = finta.at(i);
            answer.at(2*i+1) = fintb.at(i);
        } else if (i % 2 == 0) {
            answer.at(2*i-2) = finta.at(i);
            answer.at(2*i) = fintb.at(i);
        }
    }
}


int QTrPlaneStress2dSlip::giveIPValue( FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep )
{
    if ( type == IST_DisplacementVector ) {
        FloatArray u, a;
        FloatMatrix N;
        this->computeVectorOf(VM_Total, tStep, u);
        a.beSubArrayOf(u, aMask);
        this->computeNmatrixAt(gp->giveSubPatchCoordinates(), N);
        answer.beProductOf(N, a);
        return 1;
    } else if ( type == IST_MacroSlipVector ) {
        FloatArray u, b;
        FloatMatrix N;
        this->computeVectorOf(VM_Total, tStep, u);
        b.beSubArrayOf(u,bMask);
        this->computeNmatrixAt(gp->giveSubPatchCoordinates(), N);
        answer.beProductOf(N, b);
        return 1;
    } else if ( type == IST_TransferStress ) {
        StructuralSlipFE2MaterialStatus* status=static_cast< StructuralSlipFE2MaterialStatus * >( gp->giveMaterialStatus() );
        answer = status->giveTransferStressVector();
        return 1;
    } else if ( type == IST_MacroSlipGradient ) {
        FloatArray u, b;
        FloatMatrix B;
        this->computeVectorOf(VM_Total, tStep, u);
        b.beSubArrayOf(u,bMask);
        this->computeBHmatrixAt(gp, B);
        answer.beProductOf(B, b);
        return 1;
    } else if (type == IST_ReinfMembraneStress ) {
        StructuralSlipFE2MaterialStatus* status=static_cast< StructuralSlipFE2MaterialStatus *>( gp->giveMaterialStatus() );
        answer = status->giveReinfStressVector();
        return 1;
    } else {
        return Element :: giveIPValue(answer, gp, type, tStep);
    }
}


void QTrPlaneStress2dSlip::giveHomogenizedFields( FloatArray &stress, FloatArray &bStress, FloatArray &rStress, const FloatArray &strain, const FloatArray &slip, const FloatArray &slipGradient, GaussPoint *gp, TimeStep *tStep )
{
    StructuralSlipFE2Material *mat = dynamic_cast<StructuralSlipFE2Material *>( this->giveStructuralCrossSection()->giveMaterial(gp) );
    if ( mat ) {
        mat->giveHomogenizedFields(stress, bStress, rStress, strain, slip, slipGradient, gp, tStep);
    } else {
        OOFEM_ERROR("Can't homogenize the fields. StructuralSlipFE2Material needed.");
    }
}


void QTrPlaneStress2dSlip::giveSensitivities( FloatMatrix &dStressdEps, FloatMatrix &dStressdS, FloatMatrix &dStressdG, FloatMatrix &dBStressdEps, FloatMatrix &dBStressdS,
    FloatMatrix &dBStressdG, FloatMatrix &dRStressdEps, FloatMatrix &dRStressdS, FloatMatrix &dRStressdG, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep )
{
    StructuralSlipFE2Material *mat = dynamic_cast<StructuralSlipFE2Material *>( this->giveStructuralCrossSection()->giveMaterial(gp) );
    if ( mat ) {
        mat->giveSensitivities(dStressdEps, dStressdS, dStressdG, dBStressdEps, dBStressdS, dBStressdG, dRStressdEps, dRStressdS, dRStressdG, mode, gp, tStep);
        /// for debugging
//        printf("dSdE \n");
//        dStressdEps.printYourself();
//        printf("dBSdE \n");
//        dBStressdEps.printYourself();
//        printf("dRSdE \n");
//        dRStressdEps.printYourself();
//        printf("dSdS \n");
//        dStressdS.printYourself();
//        printf("dBSdS \n");
//        dBStressdS.printYourself();
//        printf("dRSdS \n");
//        dRStressdS.printYourself();
//        printf("dSdG \n");
//        dStressdG.printYourself();
//        printf("dBSdG \n");
//        dBStressdG.printYourself();
//        printf("dRSdG \n");
//        dRStressdG.printYourself();
    } else {
        OOFEM_ERROR("Can't compute sensitivities. StructuralSlipFE2Material needed.");
    }
}
} // end namespace oofem