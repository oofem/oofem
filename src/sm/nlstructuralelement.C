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
    nlGeometry = 0;
}


void
NLStructuralElement :: computeStrainVector(FloatArray &answer, GaussPoint *gp, TimeStep *stepN)
//
// Computes the vector containing the strains at the Gauss point gp of
// the receiver, at time step stepN. The nature of these strains depends
// on the element's type.
//
// element geometrical NON LINEARITY IS SUPPORTED ONLY IF BMAtrix is non-linear!!!!
// this version assumes Total Lagrange approach.
{
    int n, i;
    FloatMatrix b, A;
    FloatArray u, help;

    if ( !this->isActivated(stepN) ) {
        answer.resize( this->giveCrossSection()->giveIPValueSize(IST_StrainTensor, gp) );
        answer.zero();
        return;
    }

    fMode mode = domain->giveEngngModel()->giveFormulation();

    // === total Lagrangian formulation ===
    if ( mode == TL ) {
        // get the displacements
        this->computeVectorOf(EID_MomentumBalance, VM_Total, stepN, u);
        // subtract the initial displacements, if defined
        if ( initialDisplacements ) {
            u.subtract(*initialDisplacements);
        }

        if ( nlGeometry < 2 ) {
            // strain will be used as input for stress evaluation
            this->computeBmatrixAt(gp, b);

            // small strain tensor (in vector form)
            answer.beProductOf(b, u);

            // Green-Lagrange strain tensor (in vector form)
            // loop over all components of strain vector
            if ( nlGeometry == 1 || nlGeometry == -1) {


                n = answer.giveSize();
                for ( i = 1; i <= n; i++ ) {
                    // nonlinear part of the strain-displacement relation
                    this->computeNLBMatrixAt(A, gp, i);
                    if ( A.isNotEmpty() ) {
                        help.beProductOf(A, u);
                        answer.at(i) += 0.5 * u.dotProduct(help);
                    }
                }
            }
        } // end of nlGeometry = 0 or 1
        else { // nlGeometry = 2
            // deformation gradient will be used instead of strain
            this->computeBFmatrixAt(gp, b);
            answer.beProductOf(b, u); // this gives the displacement gradient
            // unit matrix needs to be added
            // (needs to be adjusted if the mode is not 3d)
            answer.at(1) += 1.;
            answer.at(5) += 1.;
            answer.at(9) += 1.;
        }
    } // end of total Lagrangian formulation
    else if ( mode == AL ) { // updated Lagrangian formulation
        OOFEM_ERROR("computeStrainVector : AL mode not supported now");
    }
}


void
NLStructuralElement :: computeDeformationGradientVector(FloatArray &answer, GaussPoint *gp, TimeStep *stepN)
{
    // Computes the deformation gradient stored as a vector (11, 21, 31, 12, 22, 32, 13, 23, 33) 
    // @todo rearange BF matrix to the order (11, 22, 33, 23, 13, 12, 32, 31, 21)
    FloatMatrix b;
    FloatArray u;

    this->computeVectorOf(EID_MomentumBalance, VM_Total, stepN, u); // solution vector
    // subtract the initial displacements, if defined
    if ( initialDisplacements ) {
        u.subtract(*initialDisplacements);
    }

    this->computeBFmatrixAt(gp, b);
    answer.beProductOf(b, u);   // displacement gradient H 
       
    // F = H + I 
    MaterialMode matMode = gp->giveMaterialMode();
    if ( matMode == _3dMat ) {
        answer.at(1) += 1.;
        answer.at(2) += 1.;
        answer.at(3) += 1.;
        
   } else {
    // handle plane stress/strain, axisym etc.
    OOFEM_ERROR("computeDeformationGradientVector : bad MaterialMode, only 3dMat mode is currently supported");
   }    
        

}

// remove, not in use
void 
NLStructuralElement :: computeGreenLagrangeStrainVector(FloatArray &answer, FloatArray &vF, MaterialMode matMode)
{
    // Computes the Green-Lagrange strain tensor: E=0.5(C-I)
    FloatMatrix F, E;
    F.beMatrixForm(vF);
    // 

    E.beTProductOf(F, F);    // C-Right Caucy-Green deformation tensor
    E.at(1, 1) += -1;
    E.at(2, 2) += -1;
    E.at(3, 3) += -1;
    E.times(0.5);

    FloatArray temp;
    answer.beReducedVectorFormOfStrain(E);       // Convert to Voight form Todo: add specific methods for strain/stress
}


void 
NLStructuralElement :: dyadicProductBelow(FloatMatrix &answer, FloatArray &A, FloatArray &B)
{
    // computes the non-standard "below" open product between two second order tensors
    // defined as C_ijkl = A_il * B_jk
    //@todo this should be improved by utilizing symmetry etc. C(9,9), A(9), B(6)(sym)
    answer.resize(9,9);
    answer(0,0) = A(0) * B(0); answer(0,1) = A(5) * B(5); answer(0,2) = A(4) * B(4); answer(0,3) = A(5) * B(4); 
    answer(0,4) = A(0) * B(4); answer(0,5) = A(0) * B(5); answer(0,6) = A(4) * B(5); answer(0,7) = A(4) * B(0); 
    answer(0,8) = A(5) * B(0); answer(1,0) = A(8) * B(5); answer(1,1) = A(1) * B(1); answer(1,2) = A(3) * B(3); 
    answer(1,3) = A(1) * B(3); answer(1,4) = A(8) * B(3); answer(1,5) = A(8) * B(1); answer(1,6) = A(3) * B(1); 
    answer(1,7) = A(3) * B(5); answer(1,8) = A(1) * B(5); answer(2,0) = A(7) * B(4); answer(2,1) = A(6) * B(3); 
    answer(2,2) = A(2) * B(2); answer(2,3) = A(6) * B(2); answer(2,4) = A(7) * B(2); answer(2,5) = A(7) * B(3); 
    answer(2,6) = A(2) * B(3); answer(2,7) = A(2) * B(4); answer(2,8) = A(6) * B(4); answer(3,0) = A(8) * B(4); 
    answer(3,1) = A(1) * B(3); answer(3,2) = A(3) * B(2); answer(3,3) = A(1) * B(2); answer(3,4) = A(8) * B(2); 
    answer(3,5) = A(8) * B(3); answer(3,6) = A(3) * B(3); answer(3,7) = A(3) * B(4); answer(3,8) = A(1) * B(4); 
    answer(4,0) = A(0) * B(4); answer(4,1) = A(5) * B(3); answer(4,2) = A(4) * B(2); answer(4,3) = A(5) * B(2); 
    answer(4,4) = A(0) * B(2); answer(4,5) = A(0) * B(3); answer(4,6) = A(4) * B(3); answer(4,7) = A(4) * B(4); 
    answer(4,8) = A(5) * B(4); answer(5,0) = A(0) * B(5); answer(5,1) = A(5) * B(1); answer(5,2) = A(4) * B(3); 
    answer(5,3) = A(5) * B(3); answer(5,4) = A(0) * B(3); answer(5,5) = A(0) * B(1); answer(5,6) = A(4) * B(1); 
    answer(5,7) = A(4) * B(5); answer(5,8) = A(5) * B(5); answer(6,0) = A(7) * B(5); answer(6,1) = A(6) * B(1); 
    answer(6,2) = A(2) * B(3); answer(6,3) = A(6) * B(3); answer(6,4) = A(7) * B(3); answer(6,5) = A(7) * B(1); 
    answer(6,6) = A(2) * B(1); answer(6,7) = A(2) * B(5); answer(6,8) = A(6) * B(5); answer(7,0) = A(7) * B(0); 
    answer(7,1) = A(6) * B(5); answer(7,2) = A(2) * B(4); answer(7,3) = A(6) * B(4); answer(7,4) = A(7) * B(4); 
    answer(7,5) = A(7) * B(5); answer(7,6) = A(2) * B(5); answer(7,7) = A(2) * B(0); answer(7,8) = A(6) * B(0); 
    answer(8,0) = A(8) * B(0); answer(8,1) = A(1) * B(5); answer(8,2) = A(3) * B(4); answer(8,3) = A(1) * B(4); 
    answer(8,4) = A(8) * B(4); answer(8,5) = A(8) * B(5); answer(8,6) = A(3) * B(5); answer(8,7) = A(3) * B(0); 
    answer(8,8) = A(1) * B(0); 


}


void
NLStructuralElement :: computeSecondPKStressVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep)
{
    // Computes the second Piola-Kirchoff stress vector containing the stresses at the Gauss point gp of
    // the receiver, at time step tStep. 
    StructuralCrossSection *cs = static_cast< StructuralCrossSection * >( this->giveCrossSection() );

    FloatArray vF;
    this->computeDeformationGradientVector(vF, gp, tStep);
    cs->giveSecondPKStresses(answer, ReducedForm, gp, vF, tStep);
}

// Old method
void
NLStructuralElement :: OLDgiveInternalForcesVector(FloatArray &answer,
                                                TimeStep *tStep, int useUpdatedGpRecord)
//
// returns nodal representation of real internal forces - necessary only for
// non-linear analysis.
// if useGpRecord == 1 then data stored in gp->giveStressVector() are used
// instead computing stressVector through this->ComputeStressVector();
// this must be done after you want internal forces after element->updateYourself()
// has been called for the same time step.
//
{
    GaussPoint *gp;
    Material *mat = this->giveMaterial();
    IntegrationRule *iRule = integrationRulesArray [ giveDefaultIntegrationRule() ];

    FloatMatrix b, A, *ut = NULL, b2;
    FloatArray bs, TotalStressVector, u;
    int i, j, k;
    double dV;

    // do not resize answer to computeNumberOfDofs(EID_MomentumBalance)
    // as this is valid only if receiver has no nodes with slaves
    // zero answer will resize accordingly when adding first contribution
    answer.resize(0);

    if ( nlGeometry ) {
        this->computeVectorOf(EID_MomentumBalance, VM_Total, tStep, u);
        if ( u.giveSize() ) {
            ut = new FloatMatrix( &u, 1);
        } else {
            ut = NULL;
        }
    }

    for ( i = 0; i < iRule->giveNumberOfIntegrationPoints(); i++ ) {
        gp = iRule->getIntegrationPoint(i);
        this->computeBmatrixAt(gp, b);
        if ( nlGeometry ) {
            for ( j = 1; j <= b.giveNumberOfRows(); j++ ) {
                // loop over each component of strain vector
                this->computeNLBMatrixAt(A, gp, j);
                if ( ( A.isNotEmpty() ) && ( ut != NULL ) ) {
                    b2.beProductOf(*ut,A);
                    for ( k = 1; k <= b.giveNumberOfColumns(); k++ ) {
                        // add nonlinear contribution to each component
                        b.at(j, k) += b2.at(1, k); //mj
                    }
                }
            }
        } // end nlGeometry

        if ( useUpdatedGpRecord == 1 ) {
            TotalStressVector = static_cast< StructuralMaterialStatus * >( mat->giveStatus(gp) )->giveStressVector();
        } else {
            this->computeStressVector(TotalStressVector, gp, tStep);
        }

        //
        // updates gp stress and strain record  acording to current
        // increment of displacement
        //
        if ( TotalStressVector.giveSize() == 0 ) {
            break;
        }

        //
        // now every gauss point has real stress vector
        //
        // compute nodal representation of internal forces using f = B^T*Sigma dV
        //
        dV  = this->computeVolumeAround(gp);
        bs.beTProductOf(b, TotalStressVector);
        
        answer.add(dV, bs);
    }

    if ( nlGeometry ) {
        delete ut;
    }

    // if inactive update fields; but do not contribute to structure
    if ( !this->isActivated(tStep) ) {
        answer.zero();
        return;
    }

    //answer.printYourself();

}


void
NLStructuralElement :: giveInternalForcesVector(FloatArray &answer,
                                                TimeStep *tStep, int useUpdatedGpRecord)
{
    // Returns nodal representation of real internal forces computed from first Piola-Kirchoff stress
    // if useGpRecord == 1 then data stored in gp->giveStressVector() are used
    // instead computing stressVector through this->ComputeStressVector();
    // this must be done after you want internal forces after element->updateYourself()
    // has been called for the same time step.

    if ( nlGeometry == 1 ) {
        OLDgiveInternalForcesVector(answer,tStep, useUpdatedGpRecord);
        return;
    }

    GaussPoint *gp;
    Material *mat = this->giveMaterial();
    IntegrationRule *iRule = integrationRulesArray [ giveDefaultIntegrationRule() ];

    FloatMatrix B;
    FloatArray BS, vP, vS, u, BFu;
    
    // do not resize answer to computeNumberOfDofs(EID_MomentumBalance)
    // as this is valid only if receiver has no nodes with slaves
    // zero answer will resize accordingly when adding first contribution
    answer.resize(0);

    for ( int i = 0; i < iRule->giveNumberOfIntegrationPoints(); i++ ) {
        gp = iRule->getIntegrationPoint(i);

        if ( nlGeometry == 0 ) {
            this->computeBmatrixAt(gp, B);
            if ( useUpdatedGpRecord == 1 ) {
                vS = static_cast< StructuralMaterialStatus * >( mat->giveStatus(gp) )->giveStressVector();
            } else {
                this->computeStressVector(vS, gp, tStep);
            }

        } else if ( nlGeometry == -1 ) {
         

            // updates gp stress and strain record  acording to current increment of displacement
            // now every gauss point has real stress vector
            if ( useUpdatedGpRecord == 1 ) {
                vS = static_cast< StructuralMaterialStatus * >( mat->giveStatus(gp) )->giveStressVector();
            } else {
                this->computeSecondPKStressVector(vS, gp, tStep); // currently gives S
            }

            // Compute stresses first in order to store the deformation gradient in the MaterialStatus
            // which is needed to compute GLB
            this->computeGLBMatrixAt(B, gp, tStep); 

        }

        if ( vS.giveSize() == 0 ) {
            break;
        }
        
        // compute nodal representation of internal forces at nodes as f = B^T*P dV
        double dV  = this->computeVolumeAround(gp);
        BS.beTProductOf(B, vS);
        answer.add(dV, BS);
         
    }



    // if inactive: update fields but do not give any contribution to the structure
    if ( !this->isActivated(tStep) ) {
        answer.zero();
        return;
    }

}

//@todo needs to be updated/rewritten /JB
void
NLStructuralElement :: giveInternalForcesVector_withIRulesAsSubcells(FloatArray &answer,
                                                                     TimeStep *tStep, int useUpdatedGpRecord)
//
// returns nodal representation of real internal forces - necessary only for
// non-linear analysis.
// if useGpRecord == 1 then data stored in gp->giveStressVector() are used
// instead computing stressVector through this->ComputeStressVector();
// this must be done after you want internal forces after element->updateYourself()
// has been called for the same time step.
//
{
    GaussPoint *gp;
    Material *mat = this->giveMaterial();
    IntegrationRule *iRule;

    FloatMatrix b, A, *ut = NULL, b2;
    FloatArray temp, bs, TotalStressVector, u;
    IntArray irlocnum;
    int ir, i, j, k;
    double dV;

    // do not resize answer to computeNumberOfDofs(EID_MomentumBalance)
    // as this is valid only if receiver has no nodes with slaves
    // zero answer will resize accordingly when adding first contribution
    answer.resize(0);

    if ( nlGeometry ) {
        this->computeVectorOf(EID_MomentumBalance, VM_Total, tStep, u);
        if ( u.giveSize() ) {
            ut = new FloatMatrix( &u, 1);
        } else {
            ut = NULL;
        }
    }

    FloatArray *m = & answer;
    if ( this->giveInterpolation() && this->giveInterpolation()->hasSubPatchFormulation() ) {
        m = & temp;
    }

    // loop over individual integration rules
    for ( ir = 0; ir < numberOfIntegrationRules; ir++ ) {
        iRule = integrationRulesArray [ ir ];

        for ( i = 0; i < iRule->giveNumberOfIntegrationPoints(); i++ ) {
            gp = iRule->getIntegrationPoint(i);
            this->computeBmatrixAt(gp, b);
            if ( nlGeometry ) {
                for ( j = 1; j <= b.giveNumberOfRows(); j++ ) {
                    // loop over each component of strain vector
                    this->computeNLBMatrixAt(A, gp, j);
                    if ( ( A.isNotEmpty() ) && ( ut != NULL ) ) {
                        b2.beProductOf(*ut,A);
                        for ( k = 1; k <= b.giveNumberOfColumns(); k++ ) {
                            // add nonlinear contribution to each component
                            b.at(j, k) += b2.at(1, k); //mj
                        }
                    }
                }
            } // end nlGeometry

            // TotalStressVector = gp->giveStressVector() ;
            if ( useUpdatedGpRecord == 1 ) {
                TotalStressVector = static_cast< StructuralMaterialStatus * >( mat->giveStatus(gp) )->giveStressVector();
            } else {
                this->computeStressVector(TotalStressVector, gp, tStep);
            }

            //
            // updates gp stress and strain record  acording to current
            // increment of displacement
            //
            if ( TotalStressVector.giveSize() == 0 ) {
                break;
            }

            //
            // now every gauss point has real stress vector
            //
            // compute nodal representation of internal forces using f = B^T*Sigma dV
            //
            dV  = this->computeVolumeAround(gp);
            bs.beTProductOf(b, TotalStressVector);

            m->add(dV, bs);

            // localize irule contribution into element matrix
            if ( this->giveIntegrationRuleLocalCodeNumbers(irlocnum, iRule, EID_MomentumBalance) ) {
                answer.assemble(* m, irlocnum);
                m->resize(0);
            }
        }
    } // end loop over irules

    if ( nlGeometry ) {
        delete ut;
    }

    // if inactive update fields; but do not contribute to structure
    if ( !this->isActivated(tStep) ) {
        answer.zero();
        return;
    }


}



// old method
void
NLStructuralElement :: OLDcomputeStiffnessMatrix(FloatMatrix &answer,
                                              MatResponseMode rMode, TimeStep *tStep)
//
// Computes numerically the stiffness matrix of the receiver.
// taking into account possible effects of nonlinear geometry
//
{
    int i, j, k, l, m, n, iStartIndx, iEndIndx, jStartIndx, jEndIndx;
    double dV;
    FloatMatrix d, A, *ut = NULL, b2;
    FloatMatrix bi, bj, dbj, dij;
    FloatArray u, stress;
    GaussPoint *gp;
    IntegrationRule *iRule;
    bool matStiffSymmFlag = this->giveCrossSection()->isCharacteristicMtrxSymmetric(rMode, this->material);

    answer.resize( computeNumberOfDofs(EID_MomentumBalance), computeNumberOfDofs(EID_MomentumBalance) );
    answer.zero();
    if ( !this->isActivated(tStep) ) {
        return;
    }

    Material *mat = this->giveMaterial();

    if ( nlGeometry ) {
        this->computeVectorOf(EID_MomentumBalance, VM_Total, tStep, u);
        if ( u.giveSize() ) {
            ut = new FloatMatrix( &u, 1);
        } else {
            ut = NULL;
        }
    }

    if ( numberOfIntegrationRules > 1 ) {
        for ( i = 0; i < numberOfIntegrationRules; i++ ) {
            iStartIndx = integrationRulesArray [ i ]->getStartIndexOfLocalStrainWhereApply();
            iEndIndx   = integrationRulesArray [ i ]->getEndIndexOfLocalStrainWhereApply();
            for ( j = 0; j < numberOfIntegrationRules; j++ ) {
                jStartIndx = integrationRulesArray [ j ]->getStartIndexOfLocalStrainWhereApply();
                jEndIndx   = integrationRulesArray [ j ]->getEndIndexOfLocalStrainWhereApply();
                if ( i == j ) {
                    iRule = integrationRulesArray [ i ];
                } else if ( integrationRulesArray [ i ]->giveNumberOfIntegrationPoints() < integrationRulesArray [ j ]->giveNumberOfIntegrationPoints() ) {
                    iRule = integrationRulesArray [ i ];
                } else {
                    iRule = integrationRulesArray [ j ];
                }

                for ( k = 0; k < iRule->giveNumberOfIntegrationPoints(); k++ ) {
                    gp = iRule->getIntegrationPoint(k);
                    this->computeBmatrixAt(gp, bi, iStartIndx, iEndIndx);
                    if ( i != j ) {
                        this->computeBmatrixAt(gp, bj, jStartIndx, jEndIndx);
                    } else {
                        bj = bi;
                    }

                    if ( nlGeometry ) {
                        for ( l = 0; l <  bi.giveNumberOfRows(); l++ ) {
                            // loop over each component of strain vector
                            this->computeNLBMatrixAt(A, gp, l + iStartIndx);
                            if ( ( A.isNotEmpty() ) && ( ut != NULL ) ) {
                                b2.beProductOf(* ut, A);
                                for ( m = 1; m <= bi.giveNumberOfColumns(); m++ ) {
                                    // add nonlinear contribution to each component
                                    bi.at(l + 1, m) += b2.at(1, m); //mj
                                }
                            }
                        }
                    }

                    if ( nlGeometry && ( i != j ) ) {
                        for ( l = 0; l <  bj.giveNumberOfRows(); l++ ) {
                            // loop over each component of strain vector
                            this->computeNLBMatrixAt(A, gp, l + jStartIndx);
                            if ( ( A.isNotEmpty() ) && ( ut != NULL ) ) {
                                b2.beProductOf(* ut, A);
                                for ( m = 1; m <= bj.giveNumberOfColumns(); m++ ) {
                                    // add nonlinear contribution to each component
                                    bj.at(l + 1, m) += b2.at(1, m); //mj
                                }
                            }
                        }
                    } // end nlGeometry

                    this->computeConstitutiveMatrixAt(d, rMode, gp, tStep);
                    dij.beSubMatrixOf(d, iStartIndx, iEndIndx, jStartIndx, jEndIndx);
                    dV  = this->computeVolumeAround(gp);
                    dbj.beProductOf(dij, bj);
                    if ( matStiffSymmFlag ) {
                        answer.plusProductSymmUpper(bi, dbj, dV);
                    } else {
                        answer.plusProductUnsym(bi, dbj, dV);
                    }
                }
            }
        }
    } else { // numberOfIntegrationRules == 1
        iRule = integrationRulesArray [ giveDefaultIntegrationRule() ];
        for ( j = 0; j < iRule->giveNumberOfIntegrationPoints(); j++ ) {
            gp = iRule->getIntegrationPoint(j);
            this->computeBmatrixAt(gp, bj);
            if ( nlGeometry ) {
                for ( l = 1; l <=  bj.giveNumberOfRows(); l++ ) {
                    // loop over each component of strain vector
                    this->computeNLBMatrixAt(A, gp, l);
                    if ( ( A.isNotEmpty() ) && ( ut != NULL ) ) {
                        b2.beProductOf(* ut, A);
                        for ( k = 1; k <= bj.giveNumberOfColumns(); k++ ) {
                            // add nonlinear contribution to each component
                            bj.at(l, k) += b2.at(1, k); //mj
                        }
                    }
                }
            } // end nlGeometry

            this->computeConstitutiveMatrixAt(d, rMode, gp, tStep);
            dV = this->computeVolumeAround(gp);
            dbj.beProductOf(d, bj);
            if ( matStiffSymmFlag ) {
                answer.plusProductSymmUpper(bj, dbj, dV);
            } else {
                answer.plusProductUnsym(bj, dbj, dV);
            }
        }
    }

    if ( nlGeometry ) {
        delete ut;
    }


    if ( nlGeometry ) {
        iRule = integrationRulesArray [ giveDefaultIntegrationRule() ];
        // assemble initial stress matrix
        for ( i = 0; i < iRule->giveNumberOfIntegrationPoints(); i++ ) {
            gp = iRule->getIntegrationPoint(i);
            dV = this->computeVolumeAround(gp);
            stress = static_cast< StructuralMaterialStatus * >( mat->giveStatus(gp) )->giveTempStressVector();
            n = stress.giveSize();
            if ( n ) {
                for ( j = 1; j <= n; j++ ) {
                    // loop over each component of strain vector
                    this->computeNLBMatrixAt(A, gp, j);
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
NLStructuralElement :: computeStiffnessMatrix(FloatMatrix &answer,
                                              MatResponseMode rMode, TimeStep *tStep)
{
    /* Computes the stiffness matrix of the receiver.
     * nlGeometry = 0   ->  small strain theory
     * nlGeometry = -1  ->  finite deformation theory based on Green-Lagrange strain E
     *                      and Second Piola-Kirchoff stress S in the virtual work
     */

    if ( nlGeometry == 1 ) { // old and soon redundant code
        OLDcomputeStiffnessMatrix(answer, rMode, tStep);
        return;
    }


    GaussPoint *gp;
    IntegrationRule *iRule;
    bool matStiffSymmFlag = this->giveCrossSection()->isCharacteristicMtrxSymmetric(rMode, this->material);

    answer.resize( computeNumberOfDofs(EID_MomentumBalance), computeNumberOfDofs(EID_MomentumBalance) );
    if ( !this->isActivated(tStep) ) {
        return;
    }

    // Compute matrix from material stiffness (total stiffness for small def.) - B^T * dS/dE * B
    FloatMatrix B, D, DB;
    if ( numberOfIntegrationRules == 1 ) {
        iRule = integrationRulesArray [ giveDefaultIntegrationRule() ];
        for ( int j = 0; j < iRule->giveNumberOfIntegrationPoints(); j++ ) {
            gp = iRule->getIntegrationPoint(j);
            
            if ( nlGeometry == 0 ) {
                this->computeBmatrixAt(gp, B);            
            } else if ( nlGeometry == -1 ) {
                this->computeGLBMatrixAt(B, gp, tStep); 
            }

            this->computeConstitutiveMatrixAt(D, rMode, gp, tStep); // dSdE
            double dV = this->computeVolumeAround(gp);
            DB.beProductOf(D, B);
            if ( matStiffSymmFlag ) {
                answer.plusProductSymmUpper(B, DB, dV);
            } else {
                answer.plusProductUnsym(B, DB, dV);
            }
        }

    } else { // @todo totally unchecked
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

                    if ( nlGeometry == 0 ) {
                        this->computeBmatrixAt(gp, Bi, iStartIndx, iEndIndx);

                    } else if ( nlGeometry == -1 ) { // @todo add support for Indx
                        this->computeGLBMatrixAt(Bi, gp, tStep);      
                    }

                    if ( i != j ) {
                        if ( nlGeometry == 0 ) {
                            this->computeBmatrixAt(gp, Bj, jStartIndx, jEndIndx);

                        } else if ( nlGeometry == -1 ) { // @todo add support for Indx
                            this->computeGLBMatrixAt(Bj, gp, tStep); 
                        }

                    } else {
                        Bj  = Bi;
                        BHj = BHi;
                    }

                    this->computeConstitutiveMatrixAt(D, rMode, gp, tStep);
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

    // Add geometric stiffness ("initial stress" stiffness)  - BH^T * S * BH
    // @todo The old code does not have any selective integration for this term, maybe it is never needed
    if ( nlGeometry == -1 ) {
        
        FloatArray vS, vEye;
        FloatMatrix Smat, BH, SB;
        vEye.setValues( 9, 1., 1., 1., 0., 0., 0., 0., 0., 0.);
        iRule = integrationRulesArray [ giveDefaultIntegrationRule() ];
        
        for ( int j = 0; j < iRule->giveNumberOfIntegrationPoints(); j++ ) {
            gp = iRule->getIntegrationPoint(j);

            this->computeBFmatrixAt(gp, BH);
            Material *mat = this->giveMaterial(); // shouldn't ask cross section?
            
            vS = static_cast< StructuralMaterialStatus * >( mat->giveStatus(gp) )->giveTempStressVector();
            dyadicProductBelow(Smat, vEye, vS);
            SB.beProductOf(Smat, BH);
            double dV = this->computeVolumeAround(gp);
            if ( matStiffSymmFlag ) {
                answer.plusProductSymmUpper(BH, SB, dV);
            } else {
                answer.plusProductUnsym(BH, SB, dV);
            }
        }
    }
        

    if ( matStiffSymmFlag ) {
        answer.symmetrized();
    }

}


int
NLStructuralElement :: giveVoigtIndexSym(int ind1, int ind2)
{
    // Returns the Voigt index corresponding to two given tensor indices.
    if ( ind1 == 1 && ind2 == 1 ) {
        return 1;
    } else if ( ind1 == 2 && ind2 == 2 ) {
        return 2;
    } else if ( ind1 == 3 && ind2 == 3 ) {
        return 3;
    } else if ( ( ind1 == 2 && ind2 == 3 ) || ( ind1 == 3 && ind2 == 2 ) ) {
        return 4;
    } else if ( ( ind1 == 1 && ind2 == 3 ) || ( ind1 == 3 && ind2 == 1 ) ) {
        return 5;
    } else if ( ( ind1 == 1 && ind2 == 2 ) || ( ind1 == 2 && ind2 == 1 ) ) {
        return 6;
    } else {
        OOFEM_ERROR("Error in giveVoigtIndex - bad indices");
        return -1;
    }
};



void
NLStructuralElement :: computeGLBMatrixAt(FloatMatrix &answer, GaussPoint *gp, TimeStep *tStep) 
{
    // Computes the B-matrix associated with the variation of the Green-Lagrange Strain E
    // @todo generalize to handle special case like plane strain etc

    FloatArray vH;
    this->computeBmatrixAt(gp, answer); // symmetric (small strain) B-matrix

    // Compute displacement gradient H
    Material *mat = this->giveMaterial();
    vH = static_cast< StructuralMaterialStatus * >( mat->giveStatus(gp) )->giveTempFVector();
    if ( !vH.giveSize() ) {
        this->computeDeformationGradientVector(vH, gp, tStep); 
    }
    vH.at(1) -= 1.0;
    vH.at(2) -= 1.0;
    vH.at(3) -= 1.0;

    /**
     Compute A(H) - 0.5*A(H)*H gives the quadratic terms in the GL strain 
     such that E = sym(H) + 0.5 * A(H)*H = ( B + 0.5 * A(H) * BH ) * a 

     Full 3d: 
        H = [ du/dx  dv/dy  dw/dz  dv/dz  du/dz  du/dy  dw/dy  dw/dx  dv/dx ]
        du/dx     0      0      0      0      0      0    dw/dx  dv/dx    
          0     dv/dy    0      0      0    du/dy  dw/dy    0      0
          0       0    dw/dz  dv/dz  du/dz    0      0      0      0
          0     dv/dz  dw/dy  dv/dy  du/dy  du/dz  dw/dz    0      0
        du/dz     0    dw/dx  dv/dx  du/dx    0      0    dw/dz  dv/dz
        du/dy   dv/dx    0      0      0    du/dx  dw/dx  dw/dy  dv/dy

     Plane strain - can be condensed from the above, keep: row 1,2,6, col 1, 2, 6, 9
        H = [ du/dx  dv/dy  du/dy  dv/dx ]
        du/dx     0      0     dv/dx    
          0     dv/dy  du/dy     0
        du/dy   dv/dx  du/dx   dv/dy

    */
    FloatMatrix A;
    A.resize(6,9);

    A.at(1,1) = vH.at(1); A.at(1,8) = vH.at(8); A.at(1,9) = vH.at(9);
    A.at(2,2) = vH.at(2); A.at(2,6) = vH.at(6); A.at(2,7) = vH.at(7);
    A.at(3,3) = vH.at(3); A.at(3,4) = vH.at(4); A.at(3,5) = vH.at(5);

    A.at(4,2) = vH.at(4); A.at(4,3) = vH.at(7);     
    A.at(4,4) = vH.at(2); A.at(4,5) = vH.at(6);
    A.at(4,6) = vH.at(5); A.at(4,7) = vH.at(3);
    
    A.at(5,1) = vH.at(5); A.at(5,3) = vH.at(8);     
    A.at(5,4) = vH.at(9); A.at(5,5) = vH.at(1);
    A.at(5,8) = vH.at(3); A.at(5,9) = vH.at(4);

    A.at(6,1) = vH.at(6); A.at(6,2) = vH.at(9);     
    A.at(6,6) = vH.at(1); A.at(6,7) = vH.at(8);
    A.at(6,8) = vH.at(7); A.at(6,9) = vH.at(2);

    // add A*BH to Blin
    FloatMatrix BF;
    this->computeBFmatrixAt(gp, BF);

    answer.addProductOf(A,BF);


}


//@todo rewrite for dPdF
void
NLStructuralElement :: computeStiffnessMatrix_withIRulesAsSubcells(FloatMatrix &answer,
                                                                   MatResponseMode rMode, TimeStep *tStep)
//
// Computes numerically the stiffness matrix of the receiver.
// taking into account possible effects of nonlinear geometry
//
{
    int i, j, k, l, n, ir;
    double dV;
    FloatMatrix temp, d, A, *ut = NULL, b2;
    FloatMatrix bi, bj, dbj, dij;
    FloatArray u, stress;
    GaussPoint *gp;
    IntegrationRule *iRule;
    IntArray irlocnum;
    int ndofs = computeNumberOfDofs(EID_MomentumBalance);
    answer.resize( ndofs, ndofs );
    answer.zero();
    if ( !this->isActivated(tStep) ) {
        return;
    }

    Material *mat = this->giveMaterial();

    if ( nlGeometry ) {
        this->computeVectorOf(EID_MomentumBalance, VM_Total, tStep, u);

        if ( u.giveSize() ) {
            ut = new FloatMatrix( &u, 1);
        } else {
            ut = NULL;
        }
    }

    FloatMatrix *m = & answer;
    if ( this->giveInterpolation() && this->giveInterpolation()->hasSubPatchFormulation() ) {
        m = & temp;
    }

    // loop over individual integration rules
    for ( ir = 0; ir < numberOfIntegrationRules; ir++ ) {
        iRule = integrationRulesArray [ ir ];
        for ( j = 0; j < iRule->giveNumberOfIntegrationPoints(); j++ ) {
            gp = iRule->getIntegrationPoint(j);
            this->computeBmatrixAt(gp, bj);
            if ( nlGeometry ) {
                for ( l = 1; l <=  bj.giveNumberOfRows(); l++ ) {
                    // loop over each component of strain vector
                    this->computeNLBMatrixAt(A, gp, l);
                    if ( ( A.isNotEmpty() ) && ( ut != NULL ) ) {
                        b2.beProductOf(* ut, A);
                        for ( k = 1; k <= bj.giveNumberOfColumns(); k++ ) {
                            // add nonlinear contribution to each component
                            bj.at(l, k) += b2.at(1, k); //mj
                        }
                    }
                }
            } // end nlGeometry

            this->computeConstitutiveMatrixAt(d, rMode, gp, tStep);
            dV = this->computeVolumeAround(gp);
            dbj.beProductOf(d, bj);
            m->plusProductSymmUpper(bj, dbj, dV);
        }

        if ( nlGeometry ) {
            delete ut;
        }

        // localize irule contribution into element matrix
        if ( this->giveIntegrationRuleLocalCodeNumbers(irlocnum, iRule, EID_MomentumBalance) ) {
            answer.assemble(* m, irlocnum);
            m->resize(0, 0);
        }
    }


    if ( nlGeometry ) {
        for ( ir = 0; ir < numberOfIntegrationRules; ir++ ) {
            m->resize(0, 0);
            iRule = integrationRulesArray [ ir ];

            // assemble initial stress matrix
            for ( i = 0; i < iRule->giveNumberOfIntegrationPoints(); i++ ) {
                gp = iRule->getIntegrationPoint(i);
                dV = this->computeVolumeAround(gp);
                stress = static_cast< StructuralMaterialStatus * >( mat->giveStatus(gp) )->giveStressVector();
                n = stress.giveSize();
                if ( n ) {
                    for ( j = 1; j <= n; j++ ) {
                        // loop over each component of strain vector
                        this->computeNLBMatrixAt(A, gp, j);
                        if ( A.isNotEmpty() ) {
                            A.times(stress.at(j) * dV);
                            m->add(A);
                        }
                    }
                }
            }

            // localize irule contribution into element matrix
            if ( this->giveIntegrationRuleLocalCodeNumbers(irlocnum, iRule, EID_MomentumBalance) ) {
                answer.assemble(* m, irlocnum);
                m->resize(0, 0);
            }
        }
    } // ens nlGeometry

    answer.symmetrized();
}



IRResultType
NLStructuralElement :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                // Required by IR_GIVE_FIELD macro
    this->StructuralElement :: initializeFrom(ir);

    nlGeometry = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, nlGeometry, _IFT_NLStructuralElement_nlgeoflag);

    return IRRT_OK;
}
} // end namespace oofem
