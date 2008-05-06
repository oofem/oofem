/* $Header: /home/cvs/bp/oofem/oofemlib/src/nlstructuralelement.C,v 1.15 2003/04/06 14:08:25 bp Exp $ */
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
 *               Copyright (C) 1993 - 2008   Borek Patzak
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


//   file NLSTRUCTURALELEMENT.CC


#include "nlstructuralelement.h"
#include "structuralms.h"
#include "domain.h"
#include "timestep.h"
#include "node.h"
#include "dof.h"
#include "material.h"
#include "structuralcrosssection.h"
#include "bodyload.h"
#include "gausspnt.h"
#include "integrationrule.h"
#include "intarray.h"
#include "flotarry.h"
#include "flotmtrx.h"
#include "skyline.h"
#include "debug.h"
#include "verbose.h"

#ifndef __MAKEDEPEND
#include <stdlib.h>
#include <stdio.h>
#endif


NLStructuralElement :: NLStructuralElement(int n, Domain *aDomain) :
    StructuralElement(n, aDomain)
    // Constructor. Creates an element with number n, belonging to aDomain.
{
    nlGeometry = 0;
}

/*
 * void
 * NLStructuralElement ::  giveCharacteristicMatrix (FloatMatrix& answer,
 *                       CharType mtrx, TimeStep *tStep)
 * //
 * // returns characteristics matrix of receiver accordind to mtrx
 * //
 * {
 * if (mtrx == StiffnessMatrix)
 * this -> computeStiffnessMatrix(answer, TangentStiffness, tStep);
 * else if (mtrx == TangentStiffnessMatrix)
 * this -> computeStiffnessMatrix(answer, TangentStiffness, tStep);
 * else if (mtrx == SecantStiffnessMatrix)
 * this -> computeStiffnessMatrix(answer, SecantStiffness, tStep);
 * else if (mtrx == MassMatrix)
 * this ->  computeMassMatrix(answer, tStep);
 * else _error("giveCharacteristicMatrix: Unknown Type of characteristic mtrx.");
 *
 * return ;
 * }
 *
 *
 * void
 * NLStructuralElement ::  giveCharacteristicVector (FloatArray& answer, CharType mtrx,
 *                       TimeStep *tStep)
 * //
 * // returns characteristics vector of receiver accordind to mtrx
 * //
 * {
 * if (mtrx == ElementLoadVector) this -> computeLoadVectorAt (answer, tStep);
 * else if (mtrx == NodalInternalForcesVector ) this ->giveInternalForcesVector (answer, tStep) ;
 * else _error("giveCharacteristicVector: Unknown Type of characteristic mtrx.");
 *
 * return ;
 * }
 */


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
    int n, i, rot;
    FloatMatrix b, A;
    FloatArray u, help;

    rot     = this->updateRotationMatrix();
    fMode mode = domain->giveEngngModel()->giveFormulation();
    if ( mode == TL ) { // Total Lagrange formulation
        //b       = this -> ComputeBmatrixAt(gp) ;
        this->computeBmatrixAt(gp, b);
        this->computeVectorOf(EID_MomentumBalance, VM_Total, stepN, u);

        if ( rot ) {
            u.rotatedWith(this->rotationMatrix, 'n');
        }

        // linear part of strain tensor (in vector form)
        answer.beProductOf(b, u);
        //
        // nonlin part of strain vector
        // loop over all component of strain vector
        if ( nlGeometry ) {
            n = answer.giveSize();
            for ( i = 1; i <= n; i++ ) {
                // nonlin part of strain vector
                this->computeNLBMatrixAt(A, gp, i);
                if ( A.isNotEmpty() ) {
                    help.beProductOf(A, u);
                    answer.at(i) += 0.5 * dotProduct( u, help, u.giveSize() );
                    // delete help;
                    // delete A;
                }
            }
        }
    } else if ( mode == AL ) { // actualized Lagrange formulation
        _error("computeStrainVector : AL mode not supported now");
    }

    //delete b ;
    //delete u ;
    return;
}


/*
 * void
 * NLStructuralElement :: computeStressVector (FloatArray& answer, GaussPoint* gp, TimeStep* stepN)
 * // Computes the vector containing the stresses at the Gauss point gp of
 * // the receiver, at time step stepN. The nature of these stresses depends
 * // on the element's type.
 * // this version assumes TOTAL LAGRANGE APPROACH
 * {
 * FloatArray PrevEpsilon,Epsilon, incrementOfStrains ;
 * StructuralCrossSection* cs = (StructuralCrossSection*) this->giveCrossSection();
 * Material *mat = this->giveMaterial();
 *
 * this->computeStrainVector (Epsilon, gp,stepN) ;
 * PrevEpsilon = ((StructuralMaterialStatus*) mat->giveStatus(gp)) -> giveStrainVector ();
 *
 * if (PrevEpsilon.giveSize()) {
 * incrementOfStrains = PrevEpsilon;
 * incrementOfStrains.negated() ;
 * incrementOfStrains.add (Epsilon) ;
 * } else {
 *  incrementOfStrains = Epsilon;
 * }
 *
 * cs -> giveRealStresses (answer, ReducedForm, gp, strain ,stepN);
 *
 * //delete Epsilon;
 *
 * return  ;
 * }
 */




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

    FloatMatrix b, bt, A, *ut = NULL, *b2;
    FloatArray bs, TotalStressVector, u;
    int i, j, k, rot;
    double dV;

    answer.resize(0);
    rot = this->updateRotationMatrix();

    if ( nlGeometry ) {
        this->computeVectorOf(EID_MomentumBalance, VM_Total, tStep, u);

        if ( rot ) {
            u.rotatedWith(this->rotationMatrix, 'n');
        }

        if ( u.giveSize() ) {
            ut = new FloatMatrix(& u, 1);
            //delete u;
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
                    b2 = ut->Times(& A);
                    //delete A;
                    for ( k = 1; k <= b.giveNumberOfColumns(); k++ ) {
                        // add nonlinear contribution to each component
                        b.at(j, k) += b2->at(k, 1);
                    }

                    delete b2;
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
        if ( rot ) {
            bs.rotatedWith(this->rotationMatrix, 't');
        }

        answer.add(bs);
        //delete bs;
        //delete b;
        //delete bt;
        //delete TotalStressVector;
    }

    if ( nlGeometry ) {
        delete ut;
    }

    return;
}



void
NLStructuralElement :: computeStiffnessMatrix(FloatMatrix &answer,
                                              MatResponseMode rMode, TimeStep *tStep)
//
// Computes numerically the stiffness matrix of the receiver.
// taking into account possible effects of nonlinear geometry
//
{
    int i, j, k, l, m, n, iStartIndx, iEndIndx, jStartIndx, jEndIndx, rot;
    double dV;
    FloatMatrix d, A, *ut = NULL, b2;
    FloatMatrix bi, bj, dbj, dij;
    FloatArray u, stress;
    GaussPoint *gp;
    IntegrationRule *iRule;

    Material *mat = this->giveMaterial();
    rot = this->updateRotationMatrix();

    answer.resize( computeNumberOfDofs(EID_MomentumBalance), computeNumberOfDofs(EID_MomentumBalance) );
    answer.zero();

    if ( nlGeometry ) {
        this->computeVectorOf(EID_MomentumBalance, VM_Total, tStep, u);

        if ( rot ) {
            u.rotatedWith(this->rotationMatrix, 'n');
        }

        if ( u.giveSize() ) {
            ut = new FloatMatrix(& u, 1);
            // delete u;
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
                } else if ( integrationRulesArray [ i ]->getNumberOfIntegrationPoints() < integrationRulesArray [ j ]->getNumberOfIntegrationPoints() )      {
                    iRule = integrationRulesArray [ i ];
                } else                                                                                                                                                                             {
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
                                //delete A;
                                for ( m = 1; m <= bi.giveNumberOfColumns(); m++ ) {
                                    // add nonlinear contribution to each component
                                    bi.at(l + 1, m) += b2.at(m, 1);
                                }

                                // delete b2;
                            }
                        }
                    }

                    if ( nlGeometry && ( i != j ) ) {
                        for ( l = 0; l <  bj.giveNumberOfRows(); l++ ) {
                            // loop over each component of strain vector
                            this->computeNLBMatrixAt(A, gp, l + jStartIndx);
                            if ( ( A.isNotEmpty() ) && ( ut != NULL ) ) {
                                b2.beProductOf(* ut, A);
                                //delete A;
                                for ( m = 1; m <= bj.giveNumberOfColumns(); m++ ) {
                                    // add nonlinear contribution to each component
                                    bj.at(l + 1, m) += b2.at(m, 1);
                                }

                                // delete b2;
                            }
                        }
                    } // end nlGeometry

                    this->computeConstitutiveMatrixAt(d, rMode, gp, tStep);
                    dij.beSubMatrixOf(d, iStartIndx, iEndIndx, jStartIndx, jEndIndx);
                    dV  = this->computeVolumeAround(gp);
                    dbj.beProductOf(dij, bj);
                    answer.plusProductSymmUpper(bi, dbj, dV);
                    // delete bi; delete d; delete dij; delete dbj;
                    // delete d;
                    // if (i!=j) delete bj;
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
                        //delete A;
                        for ( k = 1; k <= bj.giveNumberOfColumns(); k++ ) {
                            // add nonlinear contribution to each component
                            bj.at(l, k) += b2.at(k, 1);
                        }

                        //delete b2;
                    }
                }
            } // end nlGeometry

            //      d  = this -> giveConstitutiveMatrix() ;
            this->computeConstitutiveMatrixAt(d, rMode, gp, tStep);
            dV = this->computeVolumeAround(gp);
            dbj.beProductOf(d, bj);
            answer.plusProductSymmUpper(bj, dbj, dV);

            // delete bj ;
            // delete dbj ;
            // delete d ;
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
            stress = ( ( StructuralMaterialStatus * ) mat->giveStatus(gp) )->giveStressVector();
            n = stress.giveSize();
            if ( n ) {
                for ( j = 1; j <= n; j++ ) {
                    // loop over each component of strain vector
                    this->computeNLBMatrixAt(A, gp, j);
                    if ( A.isNotEmpty() ) {
                        A.times(stress.at(j) * dV);
                        answer.plus(A);
                    }

                    //delete A;
                }
            }
        }
    } // end nlGeometry

    answer.symmetrized();
    if ( rot ) {
        answer.rotatedWith(* this->rotationMatrix);
    }

    return;
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

