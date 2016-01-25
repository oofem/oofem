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

#include "lsmastermatgrad.h"
#include "Materials/isolinearelasticmaterial.h"
#include "gausspoint.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "intarray.h"
#include "stressvector.h"
#include "strainvector.h"
#include "mathfem.h"
#include "contextioerr.h"
#include "datastream.h"
#include "classfactory.h"

namespace oofem {

REGISTER_Material(LargeStrainMasterMaterial);

// constructor
LargeStrainMasterMaterialGrad :: LargeStrainMasterMaterialGrad(int n, Domain *d) : LargeStrainMasterMaterial(n, d), GradDpMaterialExtensionInterface(d)
{
    slaveMat = 0;
}

// destructor
LargeStrainMasterMaterialGrad :: ~LargeStrainMasterMaterialGrad()
{ }

// specifies whether a given material mode is supported by this model
int
LargeStrainMasterMaterialGrad :: hasMaterialModeCapability(MaterialMode mode)
{
    return mode == _3dMat;
}

// creates a new material status  corresponding to this class
MaterialStatus *
LargeStrainMasterMaterialGrad :: CreateStatus(GaussPoint *gp) const
{
    return new LargeStrainMasterMaterialStatus(1, this->giveDomain(), gp, slaveMat);
}


void
LargeStrainMasterMaterialGrad :: giveStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
//
// Returns characteristic material stiffness matrix of the receiver
//
{
    OOFEM_ERROR("shouldn't be called");
}


void
LargeStrainMasterMaterialGrad :: givePDGradMatrix_uu(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
    MaterialMode mMode = gp->giveMaterialMode();
    switch ( mMode ) {
    case _3dMat:
        LargeStrainMasterMaterial :: give3dMaterialStiffnessMatrix(answer, mode, gp, tStep);
        break;
    default:
        OOFEM_ERROR("unknown mode (%s)", __MaterialModeToString(mMode) );
    }
}

void
LargeStrainMasterMaterialGrad :: givePDGradMatrix_ku(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
    MaterialMode mMode = gp->giveMaterialMode();
    switch ( mMode ) {
    case _3dMat:
        give3dKappaMatrix(answer, mode, gp, tStep);
        break;
    default:
        OOFEM_ERROR("unknown mode (%s)", __MaterialModeToString(mMode) );
    }
}

void
LargeStrainMasterMaterialGrad :: givePDGradMatrix_uk(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
    MaterialMode mMode = gp->giveMaterialMode();
    switch ( mMode ) {
    case _3dMat:
        give3dGprime(answer, mode, gp, tStep);
        break;
    default:
        OOFEM_ERROR("unknown mode (%s)", __MaterialModeToString(mMode) );
    }
}

void
LargeStrainMasterMaterialGrad :: givePDGradMatrix_kk(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
    MaterialMode mMode = gp->giveMaterialMode();
    switch ( mMode ) {
    case _3dMat:
        giveInternalLength(answer, mode, gp, tStep);
        break;
    default:
        OOFEM_ERROR("unknown mode (%s)", __MaterialModeToString(mMode) );
    }
}

void
LargeStrainMasterMaterialGrad :: givePDGradMatrix_LD(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
    MaterialMode mMode = gp->giveMaterialMode();
    OOFEM_ERROR("unknown mode (%s)", __MaterialModeToString(mMode) );
}


void
LargeStrainMasterMaterialGrad :: giveInternalLength(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
    this->initTempStatus(gp);
    GradDpMaterialExtensionInterface *graddpmat = dynamic_cast< GradDpMaterialExtensionInterface * >( domain->giveMaterial(slaveMat)->giveInterface(GradDpMaterialExtensionInterfaceType) );
    if ( graddpmat == NULL ) {
        OOFEM_WARNING("material %d has no Structural support", slaveMat);
        return;
    }

    graddpmat->givePDGradMatrix_kk(answer, mode, gp, tStep);
}

void
LargeStrainMasterMaterialGrad :: give3dGprime(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
    LargeStrainMasterMaterialStatus *status = static_cast< LargeStrainMasterMaterialStatus * >( this->giveStatus(gp) );
    this->initTempStatus(gp);
    FloatMatrix gPrime;
    GradDpMaterialExtensionInterface *graddpmat = dynamic_cast< GradDpMaterialExtensionInterface * >( domain->giveMaterial(slaveMat)->giveInterface(GradDpMaterialExtensionInterfaceType) );
    if ( graddpmat == NULL ) {
        OOFEM_WARNING("material %d has no Structural support", slaveMat);
        return;
    }

    graddpmat->givePDGradMatrix_uk(gPrime, mode, gp, tStep);
    gPrime.at(4, 1) =  2. * gPrime.at(4, 1);
    gPrime.at(5, 1) =  2. * gPrime.at(5, 1);
    gPrime.at(6, 1) =  2. * gPrime.at(6, 1);
    answer.beProductOf(status->givePmatrix(), gPrime);
}


void
LargeStrainMasterMaterialGrad :: give3dKappaMatrix(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
    LargeStrainMasterMaterialStatus *status = static_cast< LargeStrainMasterMaterialStatus * >( this->giveStatus(gp) );
    this->initTempStatus(gp);
    FloatMatrix kappaMatrix;
    GradDpMaterialExtensionInterface *graddpmat = dynamic_cast< GradDpMaterialExtensionInterface * >( domain->giveMaterial(slaveMat)->giveInterface(GradDpMaterialExtensionInterfaceType) );
    if ( graddpmat == NULL ) {
        OOFEM_WARNING("material %d has no Structural support", slaveMat);
        return;
    }

    graddpmat->givePDGradMatrix_ku(kappaMatrix, mode, gp, tStep);
    kappaMatrix.at(1, 4) = 2. * kappaMatrix.at(1, 4);
    kappaMatrix.at(1, 5) = 2. * kappaMatrix.at(1, 5);
    kappaMatrix.at(1, 6) = 2. * kappaMatrix.at(1, 6);
    answer.beProductTOf(kappaMatrix, status->givePmatrix());
}



void
LargeStrainMasterMaterialGrad :: giveFirstPKStressVectorGrad(FloatArray &answer1, double &answer2, GaussPoint *gp, const FloatArray &vF, double nonlocalCumulatedStrain, TimeStep *tStep)
{
    LargeStrainMasterMaterialStatus *status = static_cast< LargeStrainMasterMaterialStatus * >( this->giveStatus(gp) );
    this->initTempStatus(gp);
    MaterialMode mode = gp->giveMaterialMode();
    if  ( mode == _3dMat ) {
        Material *mat;
        StructuralMaterial *sMat;
        mat = domain->giveMaterial(slaveMat);
        sMat = dynamic_cast< StructuralMaterial * >(mat);
        if ( sMat == NULL ) {
            OOFEM_WARNING("material %d has no Structural support", slaveMat);
            return;
        }

        GradDpMaterialExtensionInterface *dpmat = static_cast< GradDpMaterialExtensionInterface * >( sMat->giveInterface(GradDpMaterialExtensionInterfaceType) );
        if ( !dpmat ) {
            OOFEM_ERROR("Material doesn't implement the required DpGrad interface!");
        }


        double lambda1, lambda2, lambda3, E1, E2, E3;
        FloatArray eVals, SethHillStrainVector, stressVector, stressM;
        FloatMatrix F, C, eVecs, SethHillStrain;
        FloatMatrix L1, L2, T;
        //store of deformation gradient into 3x3 matrix
        F.beMatrixForm(vF);
        //compute right Cauchy-Green tensor(C), its eigenvalues and eigenvectors
        C.beTProductOf(F, F);
        // compute eigen values and eigen vectors of C
        C.jaco_(eVals, eVecs, 40);
        // compute Seth - Hill's strain measure, it depends on mParameter
        lambda1 = eVals.at(1);
        lambda2 = eVals.at(2);
        lambda3 = eVals.at(3);
        if ( m == 0 ) {
            E1 = 1. / 2. * log(lambda1);
            E2 = 1. / 2. * log(lambda2);
            E3 = 1. / 2. * log(lambda3);
        } else {
            E1 = 1. / ( 2. * m ) * ( pow(lambda1, m) - 1. );
            E2 = 1. / ( 2. * m ) * ( pow(lambda2, m) - 1. );
            E3 = 1. / ( 2. * m ) * ( pow(lambda3, m) - 1. );
        }

        SethHillStrain.resize(3, 3);
        for ( int i = 1; i < 4; i++ ) {
            for ( int j = 1; j < 4; j++ ) {
                SethHillStrain.at(i, j) = E1 * eVecs.at(i, 1) * eVecs.at(j, 1) + E2 *eVecs.at(i, 2) * eVecs.at(j, 2) + E3 *eVecs.at(i, 3) * eVecs.at(j, 3);
            }
        }



        SethHillStrainVector.beSymVectorFormOfStrain(SethHillStrain);
        dpmat->giveRealStressVectorGrad(stressVector, answer2, gp, SethHillStrainVector, nonlocalCumulatedStrain, tStep);
        this->constructTransformationMatrix(T, eVecs);

        stressVector.at(4) = 2 * stressVector.at(4);
        stressVector.at(5) = 2 * stressVector.at(5);
        stressVector.at(6) = 2 * stressVector.at(6);


        stressM.beProductOf(T, stressVector);
        stressM.at(4) = 1. / 2. *  stressM.at(4);
        stressM.at(5) = 1. / 2. *  stressM.at(5);
        stressM.at(6) = 1. / 2. *  stressM.at(6);

        this->constructL1L2TransformationMatrices(L1, L2, eVecs, stressM, E1, E2, E3);

        FloatMatrix junk, P, TL;
        FloatArray secondPK;
        junk.beProductOf(L1, T);
        P.beTProductOf(T, junk);
        //transformation of the stress to the 2PK stress and then to 1PK
        stressVector.at(4) = 0.5 * stressVector.at(4);
        stressVector.at(5) = 0.5 * stressVector.at(5);
        stressVector.at(6) = 0.5 * stressVector.at(6);
        secondPK.beProductOf(P, stressVector);
        answer1.beProductOf(F, secondPK); // P = F*S
        junk.zero();
        junk.beProductOf(L2, T);
        TL.beTProductOf(T, junk);

        status->setPmatrix(P);
        status->setTLmatrix(TL);
        status->letTempStressVectorBe(answer1);
    } else {
        OOFEM_ERROR("Unknown material mode.");
    }
}


IRResultType
LargeStrainMasterMaterialGrad :: initializeFrom(InputRecord *ir)
{
    return LargeStrainMasterMaterial :: initializeFrom(ir);
}

} // end namespace oofem
