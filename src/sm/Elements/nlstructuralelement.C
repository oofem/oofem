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

#include "../sm/Elements/nlstructuralelement.h"
#include "../sm/Materials/structuralms.h"
#include "../sm/CrossSections/structuralcrosssection.h"
#include "feinterpol.h"
#include "domain.h"
#include "material.h"
#include "crosssection.h"
#include "integrationrule.h"
#include "intarray.h"
#include "floatarray.h"
#include "floatmatrix.h"
#include "dynamicinputrecord.h"
#include "gausspoint.h"
#include "engngm.h"

namespace oofem {
NLStructuralElement :: NLStructuralElement(int n, Domain *aDomain) :
    StructuralElement(n, aDomain)
    // Constructor. Creates an element with number n, belonging to aDomain.
{
    nlGeometry = 0; // Geometrical nonlinearities disabled as default
}


double
NLStructuralElement :: computeCurrentVolume(TimeStep *tStep)
{
    double vol=0.0;
    for ( auto &gp: *this->giveDefaultIntegrationRulePtr() ) {
        FloatArray F;
        FloatMatrix Fm;


        computeDeformationGradientVector(F, gp, tStep);
        Fm.beMatrixForm(F);
        double J = Fm.giveDeterminant();

        FEInterpolation *interpolation = this->giveInterpolation();
        double detJ = fabs( ( interpolation->giveTransformationJacobian( gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) ) ) );

        vol += gp->giveWeight() * detJ * J;
    }

    return vol;
}


void
NLStructuralElement :: computeDeformationGradientVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep)
{
    // Computes the deformation gradient in the Voigt format at the Gauss point gp of
    // the receiver at time step tStep.
    // Order of components: 11, 22, 33, 23, 13, 12, 32, 31, 21 in the 3D.

    // Obtain the current displacement vector of the element and subtract initial displacements (if present)
    FloatArray u;
    this->computeVectorOf({D_u, D_v, D_w}, VM_Total, tStep, u); // solution vector
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
        OOFEM_ERROR("MaterialMode is not supported yet (%s)", __MaterialModeToString(matMode) );
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
    FloatArray vStress, vStrain, u;

    // This function can be quite costly to do inside the loops when one has many slave dofs.
    this->computeVectorOf(VM_Total, tStep, u);
    // subtract initial displacements, if defined
    if ( initialDisplacements ) {
        u.subtract(* initialDisplacements);
    }

    // zero answer will resize accordingly when adding first contribution
    answer.clear();

    for ( auto &gp: *this->giveDefaultIntegrationRulePtr() ) {
        StructuralMaterialStatus *matStat = static_cast< StructuralMaterialStatus * >( gp->giveMaterialStatus() );

        // Engineering (small strain) stress
        if ( nlGeometry == 0 ) {
            this->computeBmatrixAt(gp, B);
            if ( useUpdatedGpRecord == 1 ) {
                vStress = matStat->giveStressVector();
            } else {
                ///@todo Is this really what we should do for inactive elements?
                if ( !this->isActivated(tStep) ) {
                    vStrain.resize( StructuralMaterial :: giveSizeOfVoigtSymVector( gp->giveMaterialMode() ) );
                    vStrain.zero();
                }
                vStrain.beProductOf(B, u);
                this->computeStressVector(vStress, vStrain, gp, tStep);
            }
        } else if ( nlGeometry == 1 ) {  // First Piola-Kirchhoff stress
            if ( this->domain->giveEngngModel()->giveFormulation() == AL ) { // Cauchy stress
                if ( useUpdatedGpRecord == 1 ) {
                    vStress = matStat->giveCVector();
                } else {
                    this->computeCauchyStressVector(vStress, gp, tStep);
                }

                this->computeBmatrixAt(gp, B);
            } else { // First Piola-Kirchhoff stress
                if ( useUpdatedGpRecord == 1 ) {
                    vStress = matStat->givePVector();
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

        if ( nlGeometry == 1 ) {  // First Piola-Kirchhoff stress
            if ( vStress.giveSize() == 9 ) {
                FloatArray stressTemp;
                StructuralMaterial :: giveReducedVectorForm( stressTemp, vStress, gp->giveMaterialMode() );
                answer.plusProduct(B, stressTemp, dV);
            } else   {
                answer.plusProduct(B, vStress, dV);
            }
        } else   {
            if ( vStress.giveSize() == 6 ) {
                // It may happen that e.g. plane strain is computed
                // using the default 3D implementation. If so,
                // the stress needs to be reduced.
                // (Note that no reduction will take place if
                //  the simulation is actually 3D.)
                FloatArray stressTemp;
                StructuralMaterial :: giveReducedSymVectorForm( stressTemp, vStress, gp->giveMaterialMode() );
                answer.plusProduct(B, stressTemp, dV);
            } else   {
                answer.plusProduct(B, vStress, dV);
            }
        }
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
    /*
     * Returns nodal representation of real internal forces computed from first Piola-Kirchoff stress
     * if useGpRecord == 1 then stresses stored in the gp are used, otherwise stresses are computed
     * this must be done if you want internal forces after element->updateYourself() has been called
     * for the same time step.
     * The integration procedure uses an integrationRulesArray for numerical integration.
     * Each integration rule is considered to represent a separate sub-cell/element. Typically this would be used when
     * integration of the element domain needs special treatment, e.g. when using the XFEM.
     */

    FloatMatrix B;
    FloatArray vStress, vStrain;

    IntArray irlocnum;
    FloatArray *m = & answer, temp;
    if ( this->giveInterpolation() && this->giveInterpolation()->hasSubPatchFormulation() ) {
        m = & temp;
    }

    // zero answer will resize accordingly when adding first contribution
    answer.clear();


    // loop over individual integration rules
    for ( auto &iRule: integrationRulesArray ) {
        for ( GaussPoint *gp: *iRule ) {
            StructuralMaterialStatus *matStat = static_cast< StructuralMaterialStatus * >( gp->giveMaterialStatus() );

            if ( nlGeometry == 0 ) {
                this->computeBmatrixAt(gp, B);
                if ( useUpdatedGpRecord == 1 ) {
                    vStress = matStat->giveStressVector();
                } else {
                    this->computeStrainVector(vStrain, gp, tStep);
                    this->computeStressVector(vStress, vStrain, gp, tStep);
                }
            } else if ( nlGeometry == 1 ) {
                if ( this->domain->giveEngngModel()->giveFormulation() == AL ) { // Cauchy stress
                    if ( useUpdatedGpRecord == 1 ) {
                        vStress = matStat->giveCVector();
                    } else {
                        this->computeCauchyStressVector(vStress, gp, tStep);
                    }

                    this->computeBmatrixAt(gp, B);
                } else { // First Piola-Kirchhoff stress
                    if ( useUpdatedGpRecord == 1 ) {
                        vStress = matStat->givePVector();
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
            if ( this->giveIntegrationRuleLocalCodeNumbers(irlocnum, *iRule) ) {
                answer.assemble(* m, irlocnum);
                m->clear();
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
    bool matStiffSymmFlag = cs->isCharacteristicMtrxSymmetric(rMode);

    answer.clear();

    if ( !this->isActivated(tStep) ) {
        return;
    }

    // Compute matrix from material stiffness (total stiffness for small def.) - B^T * dS/dE * B
    if ( integrationRulesArray.size() == 1 ) {
        FloatMatrix B, D, DB;
        for ( auto &gp : *this->giveDefaultIntegrationRulePtr() ) {

            // Engineering (small strain) stiffness
            if ( nlGeometry == 0 ) {
                this->computeBmatrixAt(gp, B);
                this->computeConstitutiveMatrixAt(D, rMode, gp, tStep);
            } else if ( nlGeometry == 1 ) {
                if ( this->domain->giveEngngModel()->giveFormulation() == AL ) { // Material stiffness dC/de
                    this->computeBmatrixAt(gp, B);
                    /// @todo We probably need overloaded function (like above) here as well.
                    cs->giveStiffnessMatrix_dCde(D, rMode, gp, tStep);
                } else { // Material stiffness dP/dF
                    this->computeBHmatrixAt(gp, B);
                    /// @todo We probably need overloaded function (like above) here as well.
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
            OOFEM_ERROR("Updated lagrangian not supported yet");
        }

        int iStartIndx, iEndIndx, jStartIndx, jEndIndx;
        FloatMatrix Bi, Bj, D, Dij, DBj;
        for ( int i = 0; i < (int)integrationRulesArray.size(); i++ ) {
            iStartIndx = integrationRulesArray [ i ]->getStartIndexOfLocalStrainWhereApply();
            iEndIndx   = integrationRulesArray [ i ]->getEndIndexOfLocalStrainWhereApply();
            for ( int j = 0; j < (int)integrationRulesArray.size(); j++ ) {
                IntegrationRule *iRule;
                jStartIndx = integrationRulesArray [ j ]->getStartIndexOfLocalStrainWhereApply();
                jEndIndx   = integrationRulesArray [ j ]->getEndIndexOfLocalStrainWhereApply();
                if ( i == j ) {
                    iRule = integrationRulesArray [ i ].get();
                } else if ( integrationRulesArray [ i ]->giveNumberOfIntegrationPoints() < integrationRulesArray [ j ]->giveNumberOfIntegrationPoints() ) {
                    iRule = integrationRulesArray [ i ].get();
                } else {
                    iRule = integrationRulesArray [ j ].get();
                }

                for ( GaussPoint *gp: *iRule ) {

                    // Engineering (small strain) stiffness dSig/dEps
                    if ( nlGeometry == 0 ) {
                        this->computeBmatrixAt(gp, Bi, iStartIndx, iEndIndx);
                        this->computeConstitutiveMatrixAt(D, rMode, gp, tStep);
                    } else if ( nlGeometry == 1 ) {
                        if ( this->domain->giveEngngModel()->giveFormulation() == AL ) { // Material stiffness dC/de
                            this->computeBmatrixAt(gp, Bi);
                            cs->giveStiffnessMatrix_dCde(D, rMode, gp, tStep);
                        } else { // Material stiffness dP/dF
                            this->computeBHmatrixAt(gp, Bi);
                            cs->giveStiffnessMatrix_dPdF(D, rMode, gp, tStep);
                        }
                    }


                    if ( i != j ) {
                        if ( nlGeometry == 0 ) {
                            this->computeBmatrixAt(gp, Bj, jStartIndx, jEndIndx);
                        } else if ( nlGeometry == 1 ) {
                            if ( this->domain->giveEngngModel()->giveFormulation() == AL ) {
                                this->computeBmatrixAt(gp, Bj);
                            } else {
                                this->computeBHmatrixAt(gp, Bj);
                            }
                        }
                    } else {
                        Bj  = Bi;
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
    StructuralCrossSection *cs = this->giveStructuralCrossSection();
    bool matStiffSymmFlag = cs->isCharacteristicMtrxSymmetric(rMode);

    answer.clear();
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
    IntArray irlocnum;
    for ( auto &iRule: integrationRulesArray ) {
        for ( auto &gp : *iRule ) {

            if ( nlGeometry == 0 ) {
                this->computeBmatrixAt(gp, B);
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
        if ( this->giveIntegrationRuleLocalCodeNumbers(irlocnum, *iRule) ) {
            answer.assemble(* m, irlocnum);
            m->clear();
        }
    }

    if ( matStiffSymmFlag ) {
        answer.symmetrized();
    }
}


void
NLStructuralElement :: computeInitialStressMatrix(FloatMatrix &answer, TimeStep *tStep)
{
    FloatArray stress;
    FloatMatrix B, stress_ident, stress_identFull;
    IntArray indx;

    answer.clear();

    // assemble initial stress matrix
    for ( auto &gp : *this->giveDefaultIntegrationRulePtr() ) {
        // This function fetches the full form of the tensor
        this->giveIPValue(stress, gp, IST_StressTensor, tStep);
        if ( stress.giveSize() ) {
            // Construct the stress_ident matrix
            // The complicated part, the not-so-pretty product here: s_il delta_jk
            {
                stress_identFull.at(1, 1) = stress.at(1);
                stress_identFull.at(2, 2) = stress.at(2);
                stress_identFull.at(3, 3) = stress.at(3);
                stress_identFull.at(4, 4) = stress.at(2) + stress.at(3);
                stress_identFull.at(5, 5) = stress.at(1) + stress.at(3);
                stress_identFull.at(6, 6) = stress.at(1) + stress.at(2);

                stress_identFull.at(1, 5) = stress.at(5);
                stress_identFull.at(1, 6) = stress.at(6);

                stress_identFull.at(2, 4) = stress.at(4);
                stress_identFull.at(2, 5) = stress.at(6);

                stress_identFull.at(3, 4) = stress.at(4);
                stress_identFull.at(3, 5) = stress.at(5);

                stress_identFull.at(4, 5) = stress.at(6);
                stress_identFull.at(4, 6) = stress.at(5);
                stress_identFull.at(5, 6) = stress.at(4);
            }
            stress_ident.beSubMatrixOf(stress_identFull, indx, indx);
            stress_ident.symmetrized();
            OOFEM_WARNING("Implementation not tested yet!");

            this->computeBmatrixAt(gp, B);
            answer.plusProductSymmUpper( B, stress_ident, this->computeVolumeAround(gp) );
        }
    }

    answer.symmetrized();
}



IRResultType
NLStructuralElement :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                   // Required by IR_GIVE_FIELD macro

    nlGeometry = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, nlGeometry, _IFT_NLStructuralElement_nlgeoflag);

    return StructuralElement :: initializeFrom(ir);
}

void NLStructuralElement :: giveInputRecord(DynamicInputRecord &input)
{
    StructuralElement :: giveInputRecord(input);

    input.setField(nlGeometry, _IFT_NLStructuralElement_nlgeoflag);
}

int
NLStructuralElement :: checkConsistency()
{
    if ( this->nlGeometry == 2 ) {
        OOFEM_ERROR("nlGeometry = 2 is not supported anymore. If access to F is needed, then the material \n should overload giveFirstPKStressVector which has F as input.");
        return 0;
    }

    if ( this->nlGeometry != 0  &&  this->nlGeometry != 1 ) {
        OOFEM_ERROR("nlGeometry must be either 0 or 1 (%d not supported)", this->nlGeometry);
        return 0;
    } else {
        return 1;
    }
}
} // end namespace oofem
