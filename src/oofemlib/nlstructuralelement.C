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

#include "nlstructuralelement.h"
#include "feinterpol.h"
#include "structuralms.h"
#include "domain.h"
#include "material.h"
#include "crosssection.h"
#include "integrationrule.h"
#include "intarray.h"
#include "flotarry.h"
#include "flotmtrx.h"

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
            if ( nlGeometry == 1 ) {
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
NLStructuralElement :: giveInternalForcesVector(FloatArray &answer,
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

    FloatMatrix b, bt, A, *ut = NULL, b2;
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

    for ( i = 0; i < iRule->getNumberOfIntegrationPoints(); i++ ) {
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

        bt.beTranspositionOf(b);
        if ( useUpdatedGpRecord == 1 ) {
            TotalStressVector = ( ( StructuralMaterialStatus * ) mat->giveStatus(gp) )
                                ->giveStressVector();
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
        bs.beProductOf(bt, TotalStressVector);
        bs.times(dV);

        answer.add(bs);
    }

    if ( nlGeometry ) {
        delete ut;
    }

    // if inactive update fields; but do not contribute to structure
    if ( !this->isActivated(tStep) ) {
        answer.zero();
        return;
    }
}


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

    FloatMatrix b, bt, A, *ut = NULL, b2;
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

        for ( i = 0; i < iRule->getNumberOfIntegrationPoints(); i++ ) {
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

            bt.beTranspositionOf(b);
            // TotalStressVector = gp->giveStressVector() ;
            if ( useUpdatedGpRecord == 1 ) {
                TotalStressVector = ( ( StructuralMaterialStatus * ) mat->giveStatus(gp) )
                                    ->giveStressVector();
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
            bs.beProductOf(bt, TotalStressVector);
            bs.times(dV);

            m->add(bs);

            // localize irule contribution into element matrix
            if ( this->giveIntegrationRuleLocalCodeNumbers(irlocnum, iRule, EID_MomentumBalance) ) {
                answer.assemble(* m, irlocnum);
                m->resize(0, 0);
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




void
NLStructuralElement :: computeStiffnessMatrix(FloatMatrix &answer,
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
                } else if ( integrationRulesArray [ i ]->getNumberOfIntegrationPoints() < integrationRulesArray [ j ]->getNumberOfIntegrationPoints() ) {
                    iRule = integrationRulesArray [ i ];
                } else {
                    iRule = integrationRulesArray [ j ];
                }

                for ( k = 0; k < iRule->getNumberOfIntegrationPoints(); k++ ) {
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
        for ( j = 0; j < iRule->getNumberOfIntegrationPoints(); j++ ) {
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
        for ( i = 0; i < iRule->getNumberOfIntegrationPoints(); i++ ) {
            gp = iRule->getIntegrationPoint(i);
            dV = this->computeVolumeAround(gp);
            stress = ( ( StructuralMaterialStatus * ) mat->giveStatus(gp) )->giveTempStressVector();
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

    FloatMatrix *m = & answer;
    if ( this->giveInterpolation() && this->giveInterpolation()->hasSubPatchFormulation() ) {
        m = & temp;
    }

    // loop over individual integration rules
    for ( ir = 0; ir < numberOfIntegrationRules; ir++ ) {
        iRule = integrationRulesArray [ ir ];
        for ( j = 0; j < iRule->getNumberOfIntegrationPoints(); j++ ) {
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
            for ( i = 0; i < iRule->getNumberOfIntegrationPoints(); i++ ) {
                gp = iRule->getIntegrationPoint(i);
                dV = this->computeVolumeAround(gp);
                stress = ( ( StructuralMaterialStatus * ) mat->giveStatus(gp) )->giveStressVector();
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


/*
 * void   NLStructuralElement :: updateYourself (TimeStep* stepN)
 * // Updates the receiver at end of step.
 * {
 *
 * StructuralElement :: updateYouself(stepN) ;
 *
 * }
 */


IRResultType
NLStructuralElement :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                // Required by IR_GIVE_FIELD macro
    this->StructuralElement :: initializeFrom(ir);

    nlGeometry = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, nlGeometry, IFT_NLStructuralElement_nlgeoflag, "nlgeo"); // Macro

    return IRRT_OK;
}
} // end namespace oofem
