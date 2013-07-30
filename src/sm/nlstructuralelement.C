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

#include "nlstructuralelement.h"
#include "feinterpol.h"
#include "structuralms.h"
#include "domain.h"
#include "material.h"
#include "crosssection.h"
#include "integrationrule.h"
#include "intarray.h"
#include "floatarray.h"
#include "floatmatrix.h"
#include "structuralcrosssection.h"

namespace oofem {
NLStructuralElement :: NLStructuralElement(int n, Domain *aDomain) :
    StructuralElement(n, aDomain)
    // Constructor. Creates an element with number n, belonging to aDomain.
{
    nlGeometry = 0; // Geometrical nonlinearities disabled as default
}


void
NLStructuralElement :: computeDeformationGradientVector(FloatArray &answer, GaussPoint *gp, TimeStep *stepN)
{
    // Computes the deformation gradient in the Voigt format at the Gauss point gp of
    // the receiver at time step tStep.
    // Order of components: 11, 22, 33, 23, 13, 12, 32, 31, 21 in the 3D.

    // Obtain the current displacement vector of the element and subtract initial displacements (if present)
    FloatArray u;
    this->computeVectorOf(EID_MomentumBalance, VM_Total, stepN, u); // solution vector
    if ( initialDisplacements ) {
        u.subtract(* initialDisplacements);
    }

    // Displacement gradient H = du/dX
    FloatMatrix B;
    this->computeBHmatrixAt(gp, B);
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
        OOFEM_ERROR2( "computeDeformationGradientVector : MaterialMode is not supported yet (%s)", __MaterialModeToString(matMode) );
    }
}


void
NLStructuralElement :: computeFirstPKStressVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep)
{
    // Computes the first Piola-Kirchhoff stress vector containing the stresses at the  Gauss point
    // gp of the receiver at time step tStep. The deformation gradient F is computed and sent as
    // input to the material model.
    StructuralCrossSection *cs = this->giveStructuralCrossSection();

    FloatArray vF;
    this->computeDeformationGradientVector(vF, gp, tStep);
    cs->giveFirstPKStresses(answer, gp, vF, tStep);
}

void
NLStructuralElement :: computeCauchyStressVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep)
{
    // Computes the first Piola-Kirchhoff stress vector containing the stresses at the  Gauss point
    // gp of the receiver at time step tStep. The deformation gradient F is computed and sent as
    // input to the material model.
    StructuralCrossSection *cs = this->giveStructuralCrossSection();

    FloatArray vF;
    this->computeDeformationGradientVector(vF, gp, tStep);
    cs->giveCauchyStresses(answer, gp, vF, tStep);
}

void
NLStructuralElement :: giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord)
{
    FloatMatrix B;
    FloatArray vStress;

    // do not resize answer to computeNumberOfDofs(EID_MomentumBalance)
    // as this is valid only if receiver has no nodes with slaves
    // zero answer will resize accordingly when adding first contribution
    answer.resize(0);

    IntegrationRule *iRule = integrationRulesArray [ giveDefaultIntegrationRule() ];
    for ( int i = 0; i < iRule->giveNumberOfIntegrationPoints(); i++ ) {
        GaussPoint *gp = iRule->getIntegrationPoint(i);
        Material *mat = gp->giveMaterial();

        // Engineering (small strain) stress
        if ( nlGeometry == 0 ) {
            this->computeBmatrixAt(gp, B);
            if ( useUpdatedGpRecord == 1 ) {
                vStress = static_cast< StructuralMaterialStatus * >( mat->giveStatus(gp) )->giveStressVector();
            } else {
                this->computeStressVector(vStress, gp, tStep);
                ///@todo This is actaully inefficient since it constructs B and twice and collects the nodal unknowns over and over.
            }
        } else if ( nlGeometry == 1 ) {  // First Piola-Kirchhoff stress
            if ( this->domain->giveEngngModel()->giveFormulation() == AL ) { // Cauchy stress
                if ( useUpdatedGpRecord == 1 ) {
                    vStress = static_cast< StructuralMaterialStatus * >( mat->giveStatus(gp) )->giveCVector();
                } else {
                    this->computeCauchyStressVector(vStress, gp, tStep);
                }

                this->computeBmatrixAt(gp, B);
            } else { // First Piola-Kirchhoff stress
                if ( useUpdatedGpRecord == 1 ) {
                    vStress = static_cast< StructuralMaterialStatus * >( mat->giveStatus(gp) )->givePVector();
                } else {
                    this->computeFirstPKStressVector(vStress, gp, tStep);
                    ///@todo This is actaully inefficient since it constructs B and twice and collects the nodal unknowns over and over.
                }

                this->computeBHmatrixAt(gp, B);
            }
        }

        if ( vStress.giveSize() == 0 ) { /// @todo is this check really necessary?
            break;
        }

        // Compute nodal internal forces at nodes as f = B^T*Stress dV
        double dV  = this->computeVolumeAround(gp);
        answer.plusProduct(B, vStress, dV);
    }

    // If inactive: update fields but do not give any contribution to the internal forces
    if ( !this->isActivated(tStep) ) {
        answer.zero();
        return;
    }
}


void
NLStructuralElement :: giveInternalForcesVector_withIRulesAsSubcells(FloatArray &answer,
                                                                     TimeStep *tStep, int useUpdatedGpRecord)
{
    /**
     * Returns nodal representation of real internal forces computed from first Piola-Kirchoff stress
     * if useGpRecord == 1 then stresses stored in the gp are used, otherwise stresses are computed
     * this must be done if you want internal forces after element->updateYourself() has been called
     * for the same time step.
     * The integration procedure uses an integrationRulesArray for numerical integration.
     * Each integration rule is considered to represent a separate sub-cell/element. Typically this would be used when
     * integration of the element domain needs special treatment, e.g. when using the XFEM.
     */
    GaussPoint *gp;
    Material *mat = this->giveMaterial();
    IntegrationRule *iRule;

    FloatMatrix B;
    FloatArray BS, vP, vStress, u, BFu;

    IntArray irlocnum;
    FloatArray *m = & answer, temp;
    if ( this->giveInterpolation() && this->giveInterpolation()->hasSubPatchFormulation() ) {
        m = & temp;
    }

    // do not resize answer to computeNumberOfDofs(EID_MomentumBalance)
    // as this is valid only if receiver has no nodes with slaves
    // zero answer will resize accordingly when adding first contribution
    answer.resize(0);

    // loop over individual integration rules
    for ( int ir = 0; ir < numberOfIntegrationRules; ir++ ) {
        iRule = integrationRulesArray [ ir ];

        for ( int i = 0; i < iRule->giveNumberOfIntegrationPoints(); i++ ) {
            gp = iRule->getIntegrationPoint(i);

            if ( nlGeometry == 0 ) {
                this->computeBmatrixAt(gp, B);
                if ( useUpdatedGpRecord == 1 ) {
                    vStress = static_cast< StructuralMaterialStatus * >( mat->giveStatus(gp) )->giveStressVector();
                } else {
                    this->computeStressVector(vStress, gp, tStep);
                }
            } else if ( nlGeometry == 1 ) {
                if ( this->domain->giveEngngModel()->giveFormulation() == AL ) { // Cauchy stress
                    if ( useUpdatedGpRecord == 1 ) {
                        vStress = static_cast< StructuralMaterialStatus * >( mat->giveStatus(gp) )->giveCVector();
                    } else {
                        this->computeCauchyStressVector(vStress, gp, tStep);
                    }

                    this->computeBmatrixAt(gp, B);
                } else { // First Piola-Kirchhoff stress
                    if ( useUpdatedGpRecord == 1 ) {
                        vStress = static_cast< StructuralMaterialStatus * >( mat->giveStatus(gp) )->givePVector();
                    } else {
                        this->computeFirstPKStressVector(vStress, gp, tStep);
                    }

                    this->computeBHmatrixAt(gp, B);
                }
            }

            if ( vStress.giveSize() == 0 ) { //@todo is this really necessary?
                break;
            }

            // compute nodal representation of internal forces at nodes as f = B^T*stress dV
            double dV = this->computeVolumeAround(gp);
            m->plusProduct(B, vStress, dV);

            // localize irule contribution into element matrix
            if ( this->giveIntegrationRuleLocalCodeNumbers(irlocnum, iRule, EID_MomentumBalance) ) {
                answer.assemble(* m, irlocnum);
                m->resize(0);
            }
        }
    }

    // if inactive: update fields but do not give any contribution to the structure
    if ( !this->isActivated(tStep) ) {
        answer.zero();
        return;
    }
}






void
NLStructuralElement :: computeStiffnessMatrix(FloatMatrix &answer,
                                              MatResponseMode rMode, TimeStep *tStep)
{
    StructuralCrossSection *cs = this->giveStructuralCrossSection();
    bool matStiffSymmFlag = true;

    answer.resize( computeNumberOfDofs(EID_MomentumBalance), computeNumberOfDofs(EID_MomentumBalance) );
    if ( !this->isActivated(tStep) ) {
        return;
    }

    IntegrationRule *iRule;
    GaussPoint *gp;
    // Compute matrix from material stiffness (total stiffness for small def.) - B^T * dS/dE * B
    FloatMatrix B, D, DB;
    if ( numberOfIntegrationRules == 1 ) {
        iRule = integrationRulesArray [ giveDefaultIntegrationRule() ];
        for ( int j = 0; j < iRule->giveNumberOfIntegrationPoints(); j++ ) {
            gp = iRule->getIntegrationPoint(j);

            // Engineering (small strain) stiffness dSig/dEps
            if ( nlGeometry == 0 ) {
                this->computeBmatrixAt(gp, B);
                cs->giveCharMaterialStiffnessMatrix(D, rMode, gp, tStep);
                //cs->give_dSigdEps_StiffnessMatrix(D, rMode, gp, tStep);
            } else if ( nlGeometry == 1 ) {
                if ( this->domain->giveEngngModel()->giveFormulation() == AL ) { // Material stiffness dC/de
                    this->computeBmatrixAt(gp, B);
                    cs->giveStiffnessMatrix_dCde(D, rMode, gp, tStep);
                } else { // Material stiffness dP/dF
                    this->computeBHmatrixAt(gp, B);
                    cs->giveStiffnessMatrix_dPdF(D, rMode, gp, tStep);
                }
            }

            double dV = this->computeVolumeAround(gp);
            DB.beProductOf(D, B);
            if ( matStiffSymmFlag ) {
                answer.plusProductSymmUpper(B, DB, dV);
            } else {
                answer.plusProductUnsym(B, DB, dV);
            }
        }

        if ( this->domain->giveEngngModel()->giveFormulation() == AL ) {
            FloatMatrix initialStressMatrix;
            this->computeInitialStressMatrix(initialStressMatrix, tStep);
            answer.add(initialStressMatrix);
        }
    } else { /// @todo Verify that it works with large deformations
        if ( this->domain->giveEngngModel()->giveFormulation() == AL ) {
            OOFEM_ERROR("NLStructuralElement :: computeStiffnessMatrix - Updated lagrangian not supported yet");
        }

        int iStartIndx, iEndIndx, jStartIndx, jEndIndx;
        FloatMatrix Bi, BHi, Bj, BHj, Dij, DBj;
        for ( int i = 0; i < numberOfIntegrationRules; i++ ) {
            iStartIndx = integrationRulesArray [ i ]->getStartIndexOfLocalStrainWhereApply();
            iEndIndx   = integrationRulesArray [ i ]->getEndIndexOfLocalStrainWhereApply();
            for ( int j = 0; j < numberOfIntegrationRules; j++ ) {
                jStartIndx = integrationRulesArray [ j ]->getStartIndexOfLocalStrainWhereApply();
                jEndIndx   = integrationRulesArray [ j ]->getEndIndexOfLocalStrainWhereApply();
                if ( i == j ) {
                    iRule = integrationRulesArray [ i ];
                } else if ( integrationRulesArray [ i ]->giveNumberOfIntegrationPoints() < integrationRulesArray [ j ]->giveNumberOfIntegrationPoints() ) {
                    iRule = integrationRulesArray [ i ];
                } else {
                    iRule = integrationRulesArray [ j ];
                }

                for ( int k = 0; k < iRule->giveNumberOfIntegrationPoints(); k++ ) {
                    gp = iRule->getIntegrationPoint(k);

                    // Engineering (small strain) stiffness dSig/dEps
                    if ( nlGeometry == 0 ) {
                        this->computeBmatrixAt(gp, Bi, iStartIndx, iEndIndx);
                        cs->giveCharMaterialStiffnessMatrix(D, rMode, gp, tStep);
                        //cs->give_dSigdEps_StiffnessMatrix(D, rMode, gp, tStep);
                    } else if ( nlGeometry == 1 ) {
                        if ( this->domain->giveEngngModel()->giveFormulation() == AL ) {                         // Material stiffness dC/de
                            this->computeBmatrixAt(gp, B);
                            cs->giveStiffnessMatrix_dCde(D, rMode, gp, tStep);
                        } else {                        // Material stiffness dP/dF
                            this->computeBHmatrixAt(gp, B);
                            cs->giveStiffnessMatrix_dPdF(D, rMode, gp, tStep);
                        }
                    }


                    if ( i != j ) {
                        if ( nlGeometry == 0 ) {
                            this->computeBmatrixAt(gp, Bj, jStartIndx, jEndIndx);
                        } else if ( nlGeometry == 1 ) {
                            if ( this->domain->giveEngngModel()->giveFormulation() == AL ) {
                                this->computeBmatrixAt(gp, B);
                            } else {
                                this->computeBHmatrixAt(gp, Bj);
                            }
                        }
                    } else {
                        Bj  = Bi;
                        BHj = BHi;
                    }

                    Dij.beSubMatrixOf(D, iStartIndx, iEndIndx, jStartIndx, jEndIndx);
                    double dV = this->computeVolumeAround(gp);
                    DBj.beProductOf(Dij, Bj);
                    if ( matStiffSymmFlag ) {
                        answer.plusProductSymmUpper(Bi, DBj, dV);
                    } else {
                        answer.plusProductUnsym(Bi, DBj, dV);
                    }
                }
            }
        }
    }

    if ( matStiffSymmFlag ) {
        answer.symmetrized();
    }
}


void
NLStructuralElement :: computeStiffnessMatrix_withIRulesAsSubcells(FloatMatrix &answer,
                                                                   MatResponseMode rMode, TimeStep *tStep)
{
    GaussPoint *gp;
    IntegrationRule *iRule;
    StructuralCrossSection *cs = this->giveStructuralCrossSection();
    bool matStiffSymmFlag = cs->isCharacteristicMtrxSymmetric(rMode, this->material);
    answer.resize(0, 0);
    if ( !this->isActivated(tStep) ) {
        return;
    }

    FloatMatrix temp;
    FloatMatrix *m = & answer;
    if ( this->giveInterpolation() && this->giveInterpolation()->hasSubPatchFormulation() ) {
        m = & temp;
    }

    // Compute matrix from material stiffness
    FloatMatrix B, D, DB;
    FloatArray vS;
    FloatMatrix Smat, BH, SB;
    IntArray irlocnum;
    for ( int ir = 0; ir < numberOfIntegrationRules; ir++ ) {
        iRule = integrationRulesArray [ ir ];
        for ( int j = 0; j < iRule->giveNumberOfIntegrationPoints(); j++ ) {
            gp = iRule->getIntegrationPoint(j);

            if ( nlGeometry == 0 ) {
                this->computeBmatrixAt(gp, B);
                //cs->giveCharMaterialStiffnessMatrix(D, rMode, gp, tStep);

                //@todo This is a special treatment - asks the element instead of the cross section
                // This method is only used by one XFEM element which deals with inclusions
                this->computeConstitutiveMatrixAt(D, rMode, gp, tStep);
            } else if ( nlGeometry == 1 ) {
                if ( this->domain->giveEngngModel()->giveFormulation() == AL ) { // Material stiffness dC/de
                    this->computeBmatrixAt(gp, B);
                    cs->giveStiffnessMatrix_dCde(D, rMode, gp, tStep);
                } else { // Material stiffness dP/dF
                    this->computeBHmatrixAt(gp, B);
                    cs->giveStiffnessMatrix_dPdF(D, rMode, gp, tStep);
                }
            }

            double dV = this->computeVolumeAround(gp);
            DB.beProductOf(D, B);
            if ( matStiffSymmFlag ) {
                m->plusProductSymmUpper(B, DB, dV);
            } else {
                m->plusProductUnsym(B, DB, dV);
            }
        }

        // localize irule contribution into element matrix
        if ( this->giveIntegrationRuleLocalCodeNumbers(irlocnum, iRule, EID_MomentumBalance) ) {
            answer.assemble(* m, irlocnum);
            m->resize(0, 0);
        }
    }

    if ( matStiffSymmFlag ) {
        answer.symmetrized();
    }
}


void
NLStructuralElement :: computeInitialStressMatrix(FloatMatrix &answer, TimeStep *tStep)
{
    FloatArray stress, stressFull(6);
    FloatMatrix B, stress_ident, stress_identFull;
    IntArray indx;
    Material *mat = this->giveMaterial();

    answer.resize( computeNumberOfDofs(EID_MomentumBalance), computeNumberOfDofs(EID_MomentumBalance) );
    answer.zero();

    IntegrationRule *iRule = integrationRulesArray [ giveDefaultIntegrationRule() ];
    // assemble initial stress matrix
    for ( int i = 0; i < iRule->giveNumberOfIntegrationPoints(); i++ ) {
        GaussPoint *gp = iRule->getIntegrationPoint(i);
        stress = static_cast< StructuralMaterialStatus * >( mat->giveStatus(gp) )->giveStressVector();
        if ( stress.giveSize() ) {
            // Construct the stress_ident matrix
            StructuralMaterial :: giveVoigtSymVectorMask( indx, gp->giveMaterialMode() );
            stressFull.zero();
            stressFull.assemble(stress, indx);
            // The complicated part, the not-so-pretty product here: s_il delta_jk
            {
                stress_ident.at(1, 1) = stress.at(1);
                stress_ident.at(2, 2) = stress.at(2);
                stress_ident.at(3, 3) = stress.at(3);
                stress_ident.at(4, 4) = stress.at(2) + stress.at(3);
                stress_ident.at(5, 5) = stress.at(1) + stress.at(3);
                stress_ident.at(6, 6) = stress.at(1) + stress.at(2);

                stress_ident.at(1, 5) = stress.at(5);
                stress_ident.at(1, 6) = stress.at(6);

                stress_ident.at(2, 4) = stress.at(4);
                stress_ident.at(2, 5) = stress.at(6);

                stress_ident.at(3, 4) = stress.at(4);
                stress_ident.at(3, 5) = stress.at(5);

                stress_ident.at(4, 5) = stress.at(6);
                stress_ident.at(4, 6) = stress.at(5);
                stress_ident.at(5, 6) = stress.at(4);
            }
            stress_ident.beSubMatrixOf(stress_identFull, indx, indx);
            stress_ident.symmetrized();
            OOFEM_WARNING("NLStructuralElement :: computeInitialStressMatrix - Implementation not tested yet!");

            this->computeBmatrixAt(gp, B);
            answer.plusProductSymmUpper( B, stress_ident, this->computeVolumeAround(gp) );
        }
    }

    answer.symmetrized();
}



IRResultType
NLStructuralElement :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                   // Required by IR_GIVE_FIELD macro
    this->StructuralElement :: initializeFrom(ir);

    nlGeometry = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, nlGeometry, _IFT_NLStructuralElement_nlgeoflag);

    return IRRT_OK;
}

int
NLStructuralElement :: checkConsistency()
{
    if ( this->nlGeometry == 2 ) {
        OOFEM_ERROR("NLStructuralElement :: checkConsistency - nlGeometry = 2 is not supported anymore. If access to F is needed, then the material \n should overload giveFirstPKStressVector which has F as input.");
        return 0;
    }

    if ( this->nlGeometry != 0  &&  this->nlGeometry != 1 ) {
        OOFEM_ERROR2("NLStructuralElement :: checkConsistency - nlGeometry must be either 0 or 1 (%d not supported)", this->nlGeometry);
        return 0;
    } else {
        return 1;
    }
}
} // end namespace oofem
