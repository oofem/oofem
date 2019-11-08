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


#include "../sm/Elements/GradientDamage/graddamageelement.h"
#include "../sm/Materials/structuralms.h"
#include "../sm/CrossSections/structuralcrosssection.h"
#include "../sm/Elements/nlstructuralelement.h"
#include "../sm/Materials/graddamagematerialextensioninterface.h"
#include "node.h"
#include "material.h"
#include "gausspoint.h"
#include "gaussintegrationrule.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "intarray.h"
#include "domain.h"
#include "cltypes.h"
#include "mathfem.h"
#include "nonlocalbarrier.h"
#include "engngm.h"
#include "unknownnumberingscheme.h"

#include <cstdio>

namespace oofem {
GradientDamageElement :: GradientDamageElement()
{}



void
GradientDamageElement :: giveLocationArrayOfDofIDs(IntArray &locationArray_u, IntArray &locationArray_d, const UnknownNumberingScheme &s, const IntArray &dofIdArray_u, const IntArray &dofIdArray_d)
{
    // Routine to extract the location array of an element for given dofid array.
    locationArray_u.clear();
    locationArray_d.clear();
    NLStructuralElement *el = this->giveNLStructuralElement();
    int k = 0;
    IntArray nodalArray;
    for ( int i = 1; i <= el->giveNumberOfDofManagers(); i++ ) {
        DofManager *dMan = el->giveDofManager(i);
        int itt = 1;
        for ( int j = 1; j <= dofIdArray_u.giveSize(); j++ ) {
            if ( dMan->hasDofID( ( DofIDItem ) dofIdArray_u.at(j) ) ) {
                //  Dof *d = dMan->giveDofWithID( dofIdArray_u.at( j ) );
                locationArray_u.followedBy(k + itt);
            }
            itt++;
        }
        for ( int j = 1; j <= dofIdArray_d.giveSize(); j++ ) {
            if ( dMan->hasDofID( ( DofIDItem ) dofIdArray_d.at(j) ) ) {
                //Dof *d = dMan->giveDofWithID( dofIdArray_m.at( j ) );
                locationArray_d.followedBy(k + itt);
            }
            itt++;
        }
        k += dMan->giveNumberOfDofs();
    }
}





/*
 * void
 * GradientDamageElement :: setDisplacementLocationArray()
 * {
 * locU.resize(locSize);
 *
 * for ( int i = 1; i <= totalSize; i++ ) {
 *    if ( i < nSecNodes * nPrimVars + 1 ) {
 *        locU.at(i) = i + ( int ) ( ( ( i - 1 ) / nPrimVars ) ) * nSecVars;
 *    } else if ( i > nSecNodes * ( nPrimVars + nSecVars ) ) {
 *        locU.at(i - nSecVars * nSecNodes) = i;
 *    }
 * }
 * }
 *
 * void
 * GradientDamageElement :: setNonlocalLocationArray()
 * {
 * locK.resize(nlSize);
 * for ( int i = 1; i <= nlSize; i++ ) {
 *    locK.at(i) = i * nPrimVars + i;
 * }
 * }
 */

void
GradientDamageElement :: computeDisplacementDegreesOfFreedom(FloatArray &answer, TimeStep *tStep)
{
    this->giveStructuralElement()->computeVectorOf({ D_u, D_v, D_w }, VM_Total, tStep, answer);
}

void
GradientDamageElement :: computeNonlocalDegreesOfFreedom(FloatArray &answer, TimeStep *tStep,  ValueModeType mode)
{
    this->giveStructuralElement()->computeVectorOf({ G_0 }, mode, tStep, answer);
}

void
GradientDamageElement :: computeStressVector_and_localDamageDrivingVariable(FloatArray &answer, double &localDamageDrivingVariable, GaussPoint *gp, TimeStep *tStep)
{
    NLStructuralElement *elem = this->giveNLStructuralElement();

    double nlDamageDrivingVariable;

    int nlGeo = elem->giveGeometryMode();
    StructuralCrossSection *cs = elem->giveStructuralCrossSection();
    GradientDamageMaterialExtensionInterface *gdmat = static_cast< GradientDamageMaterialExtensionInterface * >(
        cs->giveMaterialInterface(GradientDamageMaterialExtensionInterfaceType, gp) );

    if ( !gdmat ) {
        OOFEM_ERROR("Material doesn't implement the required Gradient damage interface!");
    }

    this->computeNonlocalDamageDrivingVariable(nlDamageDrivingVariable, gp, tStep);
    if ( nlGeo == 0 ) {
        FloatArray Epsilon;
        this->computeStrainVector(Epsilon, gp, tStep);
        gdmat->giveRealStressVectorGradientDamage(answer, localDamageDrivingVariable, gp, Epsilon, nlDamageDrivingVariable, tStep);
    } else if ( nlGeo == 1 ) {
        if ( elem->giveDomain()->giveEngngModel()->giveFormulation() == TL ) {
            FloatArray vF;
            this->computeDeformationGradientVector(vF, gp, tStep);
            gdmat->giveFirstPKStressVectorGradientDamage(answer, localDamageDrivingVariable, gp, vF, nlDamageDrivingVariable, tStep);
        } else {
            FloatArray vF;
            this->computeDeformationGradientVector(vF, gp, tStep);
            gdmat->giveCauchyStressVectorGradientDamage(answer, localDamageDrivingVariable, gp, vF, nlDamageDrivingVariable, tStep);
        }
    }
}


void
GradientDamageElement :: computeStrainVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep)
{
    FloatArray u;
    FloatMatrix b;
    NLStructuralElement *elem = this->giveNLStructuralElement();

    this->computeDisplacementDegreesOfFreedom(u, tStep);
    elem->computeBmatrixAt(gp, b);
    answer.beProductOf(b, u);
}

void
GradientDamageElement :: computeDeformationGradientVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep)
{
    // Computes the deformation gradient in the Voigt format at the Gauss point gp of
    // the receiver at time step tStep.
    // Order of components: 11, 22, 33, 23, 13, 12, 32, 31, 21 in the 3D.

    // Obtain the current displacement vector of the element and subtract initial displacements (if present)
    FloatArray u;
    FloatMatrix B;
    NLStructuralElement *elem = this->giveNLStructuralElement();

    this->computeDisplacementDegreesOfFreedom(u, tStep);
    // Displacement gradient H = du/dX
    elem->computeBHmatrixAt(gp, B);
    answer.beProductOf(B, u);

    // Deformation gradient F = H + I
    MaterialMode matMode = gp->giveMaterialMode();
    if ( matMode == _3dMat || matMode == _PlaneStrain ) {
        answer.at(1) += 1.0;
        answer.at(2) += 1.0;
        answer.at(3) += 1.0;
    } else if ( matMode == _PlaneStress ) {
        answer.at(1) += 1.0;
        answer.at(2) += 1.0;
    } else if ( matMode == _1dMat ) {
        answer.at(1) += 1.0;
    } else {
        OOFEM_ERROR( "MaterialMode is not supported yet (%s)", __MaterialModeToString(matMode) );
    }
}

void
GradientDamageElement :: computeNonlocalDamageDrivingVariable(double &answer, GaussPoint *gp, TimeStep *tStep)
{
    FloatArray Nd, d;

    this->computeNdMatrixAt(gp, Nd);
    this->computeNonlocalDegreesOfFreedom(d, tStep);
    answer = Nd.dotProduct(d);
}


void
GradientDamageElement :: computeNonlocalDamageDrivingVariableGradient(FloatArray &answer, GaussPoint *gp, TimeStep *tStep)
{
    FloatMatrix Bd;
    FloatArray d;

    this->computeBdMatrixAt(gp, Bd);
    this->computeNonlocalDegreesOfFreedom(d, tStep);
    answer.beProductOf(Bd, d);
}


void
GradientDamageElement :: giveInternalForcesVector_d(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord)
{
    NLStructuralElement *elem = this->giveNLStructuralElement();

    double localDamageDrivingVariable = 0., f_dN = 0.,  nonlocalDamageDrivingVariable;
    FloatArray Nd, stress, d_d, nonlocalDamageDrivingVariable_grad, f_dB;
    FloatMatrix Bd;

    StructuralCrossSection *cs = elem->giveStructuralCrossSection();
    int size = nSecVars * nSecNodes;

    answer.resize(size);
    for ( GaussPoint *gp : *elem->giveIntegrationRule(0) ) {
        this->computeNdMatrixAt(gp, Nd);
        this->computeBdMatrixAt(gp, Bd);
        if ( useUpdatedGpRecord == 1 ) {
            StructuralMaterialStatus *ms =    static_cast< StructuralMaterialStatus * >( gp->giveMaterialStatus() );
            stress = ms->giveTempStressVector();
            /*	  StructuralCrossSection *cs = elem->giveStructuralCrossSection();
             * GradientDamageMaterialStatusExtensionInterface *gdms = dynamic_cast< GradientDamageMaterialStatusExtensionInterface * >(cs->giveMaterialInterface(GradientDamageMaterialStatusExtensionInterfaceType, gp) );
             * localCumulatedStrain = gdms->giveLocalCumulatedStrain();
             */
        } else {
            this->computeStressVector_and_localDamageDrivingVariable(stress, localDamageDrivingVariable, gp, tStep);
        }

        double dV = elem->computeVolumeAround(gp);
        this->computeNonlocalDegreesOfFreedom(d_d, tStep);

        nonlocalDamageDrivingVariable = Nd.dotProduct(d_d);
        nonlocalDamageDrivingVariable_grad.beProductOf(Bd, d_d);

        this->computeInternalForces_dN(f_dN, localDamageDrivingVariable, nonlocalDamageDrivingVariable, gp, tStep);
        this->computeInternalForces_dB(f_dB, localDamageDrivingVariable, nonlocalDamageDrivingVariable_grad, gp, tStep);


        GradientDamageMaterialExtensionInterface *gdmat = static_cast< GradientDamageMaterialExtensionInterface * >( cs->giveMaterialInterface(GradientDamageMaterialExtensionInterfaceType, gp) );
        if ( !gdmat ) {
            OOFEM_ERROR("Material doesn't implement the required Gradient damage interface!");
        }

        FloatArray N_f, B_f;
        N_f = Nd;
        N_f.times(f_dN);
        B_f.beTProductOf(Bd, f_dB);
        answer.add(dV, N_f);
        answer.add(dV, B_f);
    }


    // add penalty stiffness
    if ( penalty > 0. ) {
        FloatArray d;
        this->computeNonlocalDegreesOfFreedom(d, tStep, VM_Incremental);
        for ( int i = 1; i <= d.giveSize(); i++ ) {
            if ( d.at(i) <= 0. ) {
                answer.at(i) += penalty * d.at(i);
            }
        }
    }
}

void
GradientDamageElement :: computeInternalForces_dN(double &answer, double localDamageDrivingVariable, double nonlocalDamageDrivingVariable, GaussPoint *gp, TimeStep *tStep)
{
    NLStructuralElement *elem = this->giveNLStructuralElement();
    StructuralCrossSection *cs = elem->giveStructuralCrossSection();
    GradientDamageMaterialExtensionInterface *gdmat = static_cast< GradientDamageMaterialExtensionInterface * >(
        cs->giveMaterialInterface(GradientDamageMaterialExtensionInterfaceType, gp) );

    if ( !gdmat ) {
        OOFEM_ERROR("Material doesn't implement the required Gradient damage interface!");
    }

    gdmat->giveNonlocalInternalForces_N_factor(answer, nonlocalDamageDrivingVariable, gp, tStep);
}


void
GradientDamageElement :: computeInternalForces_dB(FloatArray &answer, double localDamageDrivingVariable, const FloatArray nonlocalDamageDrivingVariable_grad, GaussPoint *gp, TimeStep *tStep)
{
    NLStructuralElement *elem = this->giveNLStructuralElement();
    StructuralCrossSection *cs = elem->giveStructuralCrossSection();
    GradientDamageMaterialExtensionInterface *gdmat = static_cast< GradientDamageMaterialExtensionInterface * >(
        cs->giveMaterialInterface(GradientDamageMaterialExtensionInterfaceType, gp) );

    if ( !gdmat ) {
        OOFEM_ERROR("Material doesn't implement the required Gradient damage interface!");
    }

    gdmat->giveNonlocalInternalForces_B_factor(answer, nonlocalDamageDrivingVariable_grad, gp, tStep);
}


void
GradientDamageElement :: giveInternalForcesVector_u(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord)
{
    NLStructuralElement *elem = this->giveNLStructuralElement();
    int nlGeo = elem->giveGeometryMode();
    double localDamageDrivingVariable;
    FloatArray BS, vStress, vRedStress;
    FloatMatrix B;

    for ( GaussPoint *gp : *elem->giveIntegrationRule(0) ) {
        if ( nlGeo == 0 || elem->domain->giveEngngModel()->giveFormulation() == AL ) {
            elem->computeBmatrixAt(gp, B);
        } else if ( nlGeo == 1 ) {
            elem->computeBHmatrixAt(gp, B);
        }

        if ( useUpdatedGpRecord == 1 ) {
            StructuralMaterialStatus *ms =    static_cast< StructuralMaterialStatus * >( gp->giveMaterialStatus() );
            vStress = ms->giveTempStressVector();
            /*	  StructuralCrossSection *cs = elem->giveStructuralCrossSection();
             * GradientDamageMaterialStatusExtensionInterface *gdms = dynamic_cast< GradientDamageMaterialStatusExtensionInterface * >(cs->giveMaterialInterface(GradientDamageMaterialStatusExtensionInterfaceType, gp) );
             * localCumulatedStrain = gdms->giveLocalCumulatedStrain();
             */
        } else {
            this->computeStressVector_and_localDamageDrivingVariable(vStress, localDamageDrivingVariable, gp, tStep);
        }

        StructuralMaterial :: giveReducedSymVectorForm( vRedStress, vStress, gp->giveMaterialMode() );




        if ( vStress.giveSize() == 0 ) {         /// @todo is this check really necessary?
            break;
        }

        // Compute nodal internal forces at nodes as f = B^T*Stress dV
        double dV  = elem->computeVolumeAround(gp);
        BS.beTProductOf(B, vRedStress);
        answer.add(dV, BS);
    }
}


void
GradientDamageElement :: giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord)
{
    answer.resize(totalSize);
    answer.zero();
    FloatArray answerU;
    answerU.resize(locSize);
    answer.zero();
    FloatArray answerD(nlSize);
    answerD.zero();

    this->giveInternalForcesVector_u(answerU, tStep, useUpdatedGpRecord);
    this->giveInternalForcesVector_d(answerD, tStep, useUpdatedGpRecord);


    answer.assemble(answerU, locationArray_u);
    answer.assemble(answerD, locationArray_d);
}




void
GradientDamageElement :: computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
{
    answer.resize(totalSize, totalSize);
    answer.zero();

    FloatMatrix answer1, answer2, answer3, answer4;
    this->computeStiffnessMatrix_uu(answer1, rMode, tStep);
    this->computeStiffnessMatrix_ud(answer2, rMode, tStep);
    this->computeStiffnessMatrix_du(answer3, rMode, tStep);
    this->computeStiffnessMatrix_dd(answer4, rMode, tStep);
    answer.assemble(answer1, locationArray_u);
    answer.assemble(answer2, locationArray_u, locationArray_d);
    answer.assemble(answer3, locationArray_d, locationArray_u);
    answer.assemble(answer4, locationArray_d);
}


void
GradientDamageElement :: computeStiffnessMatrix_uu(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
{
    NLStructuralElement *elem = this->giveNLStructuralElement();
    StructuralCrossSection *cs = elem->giveStructuralCrossSection();
    FloatMatrix B, D, DB;

    int nlGeo = elem->giveGeometryMode();
    bool matStiffSymmFlag = elem->giveCrossSection()->isCharacteristicMtrxSymmetric(rMode);
    answer.clear();
    for ( GaussPoint *gp : *elem->giveIntegrationRule(0) ) {
        GradientDamageMaterialExtensionInterface *gdmat = dynamic_cast< GradientDamageMaterialExtensionInterface * >(
            cs->giveMaterialInterface(GradientDamageMaterialExtensionInterfaceType, gp) );
        if ( !gdmat ) {
            OOFEM_ERROR("Material doesn't implement the required DpGrad interface!");
        }
        if ( nlGeo == 0 ) {
            elem->computeBmatrixAt(gp, B);
        } else if ( nlGeo == 1 ) {
            if ( elem->domain->giveEngngModel()->giveFormulation() == AL ) {
                elem->computeBmatrixAt(gp, B);
            } else {
                elem->computeBHmatrixAt(gp, B);
            }
        }
        gdmat->giveGradientDamageStiffnessMatrix_uu(D, rMode, gp, tStep);
        double dV = elem->computeVolumeAround(gp);
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
GradientDamageElement :: computeStiffnessMatrix_du(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
{
    double dV;
    NLStructuralElement *elem = this->giveNLStructuralElement();
    FloatArray Nd;
    FloatMatrix B, DduB, Ddu;
    StructuralCrossSection *cs = elem->giveStructuralCrossSection();

    answer.clear();

    int nlGeo = elem->giveGeometryMode();

    for ( auto &gp : *elem->giveIntegrationRule(0) ) {
        GradientDamageMaterialExtensionInterface *gdmat = dynamic_cast< GradientDamageMaterialExtensionInterface * >(
            cs->giveMaterialInterface(GradientDamageMaterialExtensionInterfaceType, gp) );
        if ( !gdmat ) {
            OOFEM_ERROR("Material doesn't implement the required DpGrad interface!");
        }

        elem->computeBmatrixAt(gp, B);
        if ( nlGeo == 1 ) {
            if ( elem->domain->giveEngngModel()->giveFormulation() == AL ) {
                elem->computeBmatrixAt(gp, B);
            } else {
                elem->computeBHmatrixAt(gp, B);
            }
        }

        gdmat->giveGradientDamageStiffnessMatrix_du(Ddu, rMode, gp, tStep);
        this->computeNdMatrixAt(gp, Nd);
        dV = elem->computeVolumeAround(gp);
        if ( Ddu.giveNumberOfRows() > 0 ) {
            DduB.beTProductOf(Ddu, B);
            FloatMatrix Ndm(Nd, true);
            answer.plusProductUnsym(Ndm, DduB, dV);
        } else {
            if ( answer.giveNumberOfColumns() == 0 ) {
                answer.resize( Nd.giveSize(), B.giveNumberOfColumns() );
            }
        }
    }
}


void
GradientDamageElement :: computeStiffnessMatrix_dd(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
{
    StructuralElement *elem = this->giveStructuralElement();
    double dV;
    FloatMatrix Ddd_NN, Ddd_BB, Ddd_BN;
    FloatArray Nd;
    FloatMatrix Bd;
    FloatMatrix Ddd;
    StructuralCrossSection *cs = elem->giveStructuralCrossSection();

    answer.clear();
    for ( auto &gp : *elem->giveIntegrationRule(0) ) {
        GradientDamageMaterialExtensionInterface *gdmat = dynamic_cast< GradientDamageMaterialExtensionInterface * >(
            cs->giveMaterialInterface(GradientDamageMaterialExtensionInterfaceType, gp) );
        if ( !gdmat ) {
            OOFEM_ERROR("Material doesn't implement the required Gradient damage interface!");
        }

        this->computeNdMatrixAt(gp, Nd);
        FloatMatrix Ndm(Nd, true);
        this->computeBdMatrixAt(gp, Bd);
        dV = elem->computeVolumeAround(gp);

        gdmat->giveGradientDamageStiffnessMatrix_dd_NN(Ddd_NN, rMode, gp, tStep);
        gdmat->giveGradientDamageStiffnessMatrix_dd_BB(Ddd_BB, rMode, gp, tStep);
        gdmat->giveGradientDamageStiffnessMatrix_dd_BN(Ddd_BN, rMode, gp, tStep);



        if ( Ddd_NN.giveNumberOfRows() > 0 ) {
	  
            Ddd.beProductOf(Ddd_NN, Ndm);
            answer.plusProductUnsym(Ndm, Ddd, dV);
        } else {
            answer.plusProductUnsym(Ndm, Ndm, dV);
        }

        Ddd.beProductOf(Ddd_BB, Bd);
        answer.plusProductUnsym(Bd, Ddd, dV);



        if ( Ddd_BN.giveNumberOfRows() > 0 ) {
            Ddd.beProductOf(Ddd_BN, Ndm);
            answer.plusProductUnsym(Bd, Ddd, dV);
        }
    }

    // add penalty stiffness
    if ( penalty > 0. ) {
        FloatArray d;
        this->computeNonlocalDegreesOfFreedom(d, tStep, VM_Incremental);
        for ( int i = 1; i <= d.giveSize(); i++ ) {
            if ( d.at(i) <= 0. ) {
                answer.at(i, i) += penalty;
            }
        }
    }
}


void
GradientDamageElement :: computeStiffnessMatrix_ud(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
{
    NLStructuralElement *elem = this->giveNLStructuralElement();
    double dV;
    StructuralCrossSection *cs = elem->giveStructuralCrossSection();
    int nlGeo = elem->giveGeometryMode();
    FloatArray Nd;
    FloatMatrix B;
    FloatMatrix Dud, DudN;

    answer.clear();

    for ( auto &gp : *elem->giveIntegrationRule(0) ) {
        GradientDamageMaterialExtensionInterface *gdmat = dynamic_cast< GradientDamageMaterialExtensionInterface * >(
            cs->giveMaterialInterface(GradientDamageMaterialExtensionInterfaceType, gp) );
        if ( !gdmat ) {
            OOFEM_ERROR("Material doesn't implement the required gradient interface!");
        }
        gdmat->giveGradientDamageStiffnessMatrix_ud(Dud, rMode, gp, tStep);
        this->computeNdMatrixAt(gp, Nd);
        elem->computeBmatrixAt(gp, B);
        if ( nlGeo == 1 ) {
            if ( elem->domain->giveEngngModel()->giveFormulation() == AL ) {
                elem->computeBmatrixAt(gp, B);
            } else {
                elem->computeBHmatrixAt(gp, B);
            }
        }
        dV = elem->computeVolumeAround(gp);

        DudN.beProductTOf(Dud, Nd);
        answer.plusProductUnsym(B, DudN, dV);
    }
}


void
GradientDamageElement :: initializeFrom(InputRecord &ir)
{
    //nlGeo = 0;
    penalty = 0.;
    IR_GIVE_OPTIONAL_FIELD(ir, penalty, _IFT_GradientDamageElement_penalty);
}

void
GradientDamageElement :: postInitialize()
{
    this->giveLocationArray_u(locationArray_u);
    this->giveLocationArray_d(locationArray_d);
}


/*
 * void
 * GradientDamageElement :: postInitialize()
 * {
 * IntArray IdMask_u, IdMask_d;
 * this->giveDofManDofIDMask_u( IdMask_u );
 * this->giveDofManDofIDMask_d( IdMask_d );
 * this->giveLocationArrayOfDofIDs(locationArray_u,locationArray_d, EModelDefaultEquationNumbering(), IdMask_u, IdMask_d);
 * }
 */
} // end namespace oofem
