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
 *               Copyright (C) 1993 - 2011   Borek Patzak
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

#include "structuralelementevaluator.h"
#include "flotarry.h"
#include "flotmtrx.h"
#include "domain.h"
#include "node.h"
#include "element.h"
#include "gausspnt.h"
#include "gaussintegrationrule.h"
#include "matresponsemode.h"
#include "crosssection.h"
#include "structuralcrosssection.h"
#include "mathfem.h"
#include "iga.h"

namespace oofem {
StructuralElementEvaluator :: StructuralElementEvaluator()
{
    this->rotationMatrix.beEmptyMtrx();
}

/*
 * int StructuralElementEvaluator :: giveIntegrationElementCodeNumbers(IntArray &answer, Element *elem,
 *                                                                  IntegrationRule *ie, EquationID ut) {
 *  int i;
 *  IntArray mask, nodeDofIDMask, nodalArray;
 *
 *  // first evaluate nonzero basis function mask
 *  if ( elem->giveInterpolation()->hasSubPatchFormulation() ) {
 *      IGAIntegrationElement *ee = ( IGAIntegrationElement * ) ie;
 *      elem->giveInterpolation()->giveKnotSpanBasisFuncMask(* ee->giveKnotSpan(), mask);
 *      // loop over nonzero shape functions and assemble localization array
 *      answer.resize(0);
 *      for ( i = 1; i <= mask.giveSize(); i++ ) {
 *          elem->giveDofManDofIDMask(mask.at(i), ut, nodeDofIDMask);
 *          elem->giveDofManager( mask.at(i) )->giveLocationArray(nodeDofIDMask, nodalArray);
 *          answer.followedBy(nodalArray);
 *      }
 *
 *      return 1;
 *  } else {
 *      return 0;
 *  }
 * }
 */

int StructuralElementEvaluator :: giveIntegrationElementLocalCodeNumbers(IntArray &answer, Element *elem,
                                                                         IntegrationRule *ie, EquationID ut)
{
    int i, j, nsd;
    IntArray mask, nodeDofIDMask, nodalArray;
    int dofmandof;

    // get number of dofs in node
    elem->giveDofManDofIDMask(1, ut, nodeDofIDMask);
    dofmandof = nodeDofIDMask.giveSize();

    nsd = elem->giveInterpolation()->giveNsd();

    // first evaluate nonzero basis function mask
    if ( elem->giveInterpolation()->hasSubPatchFormulation() ) {
        IGAIntegrationElement *ee = ( IGAIntegrationElement * ) ie;
        elem->giveInterpolation()->giveKnotSpanBasisFuncMask(* ee->giveKnotSpan(), mask);
        // loop over nonzero shape functions and assemble localization array
        answer.resize(0);
        for ( i = 1; i <= mask.giveSize(); i++ ) {
            nodalArray.resize( nodeDofIDMask.giveSize() );
            for ( j = 1; j <= nsd; j++ ) {
                nodalArray.at(j) = dofmandof * ( mask.at(i) - 1 ) + j;
            }

            answer.followedBy(nodalArray);
        }

        return 1;
    } else {
        return 0;
    }
}



void StructuralElementEvaluator :: giveCharacteristicMatrix(FloatMatrix &answer,
                                                            CharType mtrx, TimeStep *tStep)
//
// returns characteristics matrix of receiver accordind to mtrx
//
{
    if ( mtrx == StiffnessMatrix ) {
        this->computeStiffnessMatrix(answer, TangentStiffness, tStep);
    } else {
        OOFEM_ERROR2( "giveCharacteristicMatrix: Unknown Type of characteristic mtrx (%s)", __CharTypeToString(mtrx) );
    }

    return;
}

void StructuralElementEvaluator :: computeBcLoadVectorAt(FloatArray &answer, TimeStep *stepN, ValueModeType mode)
// Computes the load vector due to the boundary conditions acting on the
// receiver's nodes, at stepN. Returns NULL if this array contains only
// zeroes.
{
    FloatArray d, dp;
    FloatMatrix s;
    Element *elem = this->giveElement();
    int numberOfDofMans = elem->giveNumberOfDofManagers();
    /*
     * this -> computeVectorOfPrescribed(DisplacementVector,TotalMode,stepN, d) ;
     * if ((stepN->giveLoadResponseMode()==IncrementOfLoad) && (!stepN->isTheFirstStep())) {
     * this -> computeVectorOfPrescribed(DisplacementVector,TotalMode,stepN->givePreviousStep(), dp);
     * d.subtract (dp);
     * //delete dp;
     * }
     */
    this->computeVectorOfPrescribed(EID_MomentumBalance, mode, stepN, d);
    //this -> computeVectorOfPrescribed(DisplacementVector,umode,stepN, d) ;

    if ( d.containsOnlyZeroes() ) {
        answer.resize(0);
    } else {
        this->computeStiffnessMatrix(s, TangentStiffness, stepN);
        answer.beProductOf(s, d);
        answer.negated();
    }

    // delete d ;

    // if engngmodel supports dynamic change of static system
    // we must test if element has not been removed in previous step
    // if not, we must also test if there was previous BC on some DOF and now it is released.
    // if was, it is necessary to load it by reaction force.
    if ( elem->giveDomain()->giveEngngModel()->requiresUnknownsDictionaryUpdate() ) {
        FloatArray prevInternalForces;
        IntArray elementNodeMask, dofMask;
        DofManager *nodeI;
        Dof *dofJ;
        int nDofs, i, j, k = 0;

        if ( ( mode == VM_Incremental ) && ( !stepN->isTheFirstStep() ) ) {
            for ( i = 1; i <= numberOfDofMans; i++ ) {
                nodeI = elem->giveDofManager(i);
                elem->giveDofManDofIDMask(i, EID_MomentumBalance, elementNodeMask);
                nodeI->giveDofArray(elementNodeMask, dofMask);
                nDofs = dofMask.giveSize();
                for ( j = 1; j <= nDofs; j++ ) {
                    dofJ = nodeI->giveDof( dofMask.at(j) );
                    k++;
                    if ( !dofJ->hasBc(stepN) && dofJ->hasBc( stepN->givePreviousStep() ) ) {
                        if ( prevInternalForces.giveSize() == 0 ) {
                            // allocate and compute only if needed
                            // use updated gp record
                            this->giveInternalForcesVector(prevInternalForces,
                                                           stepN->givePreviousStep(), 1);
                        }

                        // check for allocated answer
                        if ( answer.giveSize() == 0 ) {
                            answer.resize( elem->computeNumberOfDofs(EID_MomentumBalance) );
                            answer.zero();
                        }

                        // add element part of reaction  to load vector
                        answer.at(k) -= prevInternalForces.at(k);
                    }
                }

                //delete elementNodeMask;
                //delete dofMask;
            }
        }
    }

    return;
}




void StructuralElementEvaluator :: computeStrainVector(FloatArray &answer, GaussPoint *gp, TimeStep *stepN)
{
    FloatArray u;
    Element *elem = this->giveElement();


    elem->computeVectorOf(EID_MomentumBalance, VM_Total, stepN, u);

    /*
     * // subtract initial displacements, if defined
     * if (initialDisplacements) u.subtract(initialDisplacements);
     */

    this->computeStrainVector(answer, gp, stepN, u);
}

void StructuralElementEvaluator :: computeStrainVector(FloatArray &answer, GaussPoint *gp, TimeStep *stepN, FloatArray &u)
// Computes the vector containing the strains at the Gauss point gp of
// the receiver, at time step stepN. The nature of these strains depends
// on the element's type.
{
    int i;
    FloatMatrix b;
    FloatArray ur;
    Element *elem = this->giveElement();

    if ( !this->isActivated(stepN) ) {
        answer.resize( elem->giveCrossSection()->giveIPValueSize(IST_StrainTensor, gp) );
        answer.zero();
        return;
    }

    this->computeBMatrixAt(b, gp);

    // get local code numbers corresponding to ir
    IntArray lc;
    this->giveIntegrationElementLocalCodeNumbers(lc, elem, gp->giveIntegrationRule(), EID_MomentumBalance);
    ur.resize( b.giveNumberOfColumns() );
    for ( i = 1; i <= lc.giveSize(); i++ ) {
        ur.at(i) = u.at( lc.at(i) );
    }

    answer.beProductOf(b, ur);
}


void StructuralElementEvaluator :: computeStressVector(FloatArray &answer, GaussPoint *gp, TimeStep *stepN)
// Computes the vector containing the stresses at the Gauss point gp of
// the receiver, at time step stepN. The nature of these stresses depends
// on the element's type.
// this version assumes TOTAL LAGRANGE APPROACH
{
    FloatArray Epsilon;
    Element *elem = this->giveElement();
    StructuralCrossSection *cs = ( StructuralCrossSection * ) elem->giveCrossSection();

    this->computeStrainVector(Epsilon, gp, stepN);
    cs->giveRealStresses(answer, ReducedForm, gp, Epsilon, stepN);
}


void StructuralElementEvaluator :: computeStressVector(FloatArray &answer, GaussPoint *gp, TimeStep *stepN, FloatArray &u)
// Computes the vector containing the stresses at the Gauss point gp of
// the receiver, at time step stepN. The nature of these stresses depends
// on the element's type.
// this version assumes TOTAL LAGRANGE APPROACH
{
    FloatArray Epsilon;
    Element *elem = this->giveElement();
    StructuralCrossSection *cs = ( StructuralCrossSection * ) elem->giveCrossSection();

    this->computeStrainVector(Epsilon, gp, stepN, u);
    cs->giveRealStresses(answer, ReducedForm, gp, Epsilon, stepN);
}

void StructuralElementEvaluator :: updateInternalState(TimeStep *stepN)
// Updates the receiver at end of step.
{
    FloatArray u;
    Element *elem = this->giveElement();

    elem->computeVectorOf(EID_MomentumBalance, VM_Total, stepN, u);

    /*
     * // subtract initial displacements, if defined
     * if (initialDisplacements) u.subtract(initialDisplacements);
     */

    int i, j;
    IntegrationRule *iRule;
    FloatArray stress;

    // force updating strains & stresses
    for ( i = 0; i < elem->giveNumberOfIntegrationRules(); i++ ) {
#ifdef __PARALLEL_MODE
      if (this->giveElement()->giveKnotSpanParallelMode(i) == Element_remote) continue;
#endif
        iRule = elem->giveIntegrationRule(i);
        for ( j = 0; j < iRule->getNumberOfIntegrationPoints(); j++ ) {
            computeStressVector(stress, iRule->getIntegrationPoint(j), stepN, u);
        }
    }

    /*
     * // Original unoptimized version
     * int i, j;
     * IntegrationRule *iRule;
     * FloatArray stress;
     * Element *elem = this->giveElement();
     *
     * // force updating strains & stresses
     * for ( i = 0; i < elem->giveNumberOfIntegrationRules(); i++ ) {
     * #ifdef __PARALLEL_MODE
     *     if (this->giveElement()->giveKnotSpanParallelMode(i) == Element_remote) continue;
     * #endif
     *  iRule = elem->giveIntegrationRule(i);
     *    for ( j = 0; j < iRule->getNumberOfIntegrationPoints(); j++ ) {
     *      computeStressVector(stress, iRule->getIntegrationPoint(j), stepN);
     *    }
     * }
     */
}


void StructuralElementEvaluator :: computeNonForceLoadVector(FloatArray &answer, TimeStep *stepN, ValueModeType mode)
// Computes the load vector of the receiver, at stepN.
{
    FloatArray helpLoadVector;

    answer.resize(0);

    // test for deactivation of receiver
    if ( ( mode == VM_Incremental ) && ( !stepN->isTheFirstStep() ) ) {
        if ( isActivated( stepN->givePreviousStep() ) && !isActivated(stepN) ) {
            // use updated gp record
            this->giveInternalForcesVector(answer, stepN->givePreviousStep(), 1);
        }
    }

    if ( !this->isActivated(stepN) ) {
        return;
    }

    /*
     * this->computePrescribedStrainLoadVectorAt(helpLoadVector, stepN, mode);
     * if ( helpLoadVector.giveSize() ) {
     * answer.add(helpLoadVector);
     * }
     */

    this->computeBcLoadVectorAt(helpLoadVector, stepN, mode);
    if ( helpLoadVector.giveSize() ) {
        answer.add(helpLoadVector);
    }
}


void StructuralElementEvaluator :: computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
{
    int ir, j, numberOfIntegrationRules;
    FloatMatrix temp, bj, d, dbj;
    IntegrationRule *iRule;
    GaussPoint *gp;
    Element *elem = this->giveElement();
    int ndofs = elem->computeNumberOfDofs(EID_MomentumBalance);
    bool matStiffSymmFlag = elem->giveCrossSection()->isCharacteristicMtrxSymmetric( rMode, elem->giveMaterial()->giveNumber() );
    IntArray irlocnum;
    double dV;

    answer.resize(ndofs, ndofs);
    answer.zero();

    FloatMatrix *m = & answer;
    if ( elem->giveInterpolation()->hasSubPatchFormulation() ) {
        m = & temp;
    }

    numberOfIntegrationRules = elem->giveNumberOfIntegrationRules();
    // loop over individual integration rules
    for ( ir = 0; ir < numberOfIntegrationRules; ir++ ) {

#ifdef __PARALLEL_MODE
      if (this->giveElement()->giveKnotSpanParallelMode(ir) == Element_remote) continue;
      //fprintf (stderr, "[%d] Computing element.knotspan %d.%d\n", elem->giveDomain()->giveEngngModel()->giveRank(), elem->giveNumber(), ir);
#endif
        m->resize(0, 0);
        iRule = elem->giveIntegrationRule(ir);
        // loop over individual integration points
        for ( j = 0; j < iRule->getNumberOfIntegrationPoints(); j++ ) {
            gp = iRule->getIntegrationPoint(j);
            this->computeBMatrixAt(bj, gp);
            //elem->computeConstitutiveMatrixAt(d, rMode, gp, tStep);
            ( ( StructuralCrossSection * ) elem->giveCrossSection() )
            ->giveCharMaterialStiffnessMatrix(d, rMode, gp, tStep);

            dV = this->computeVolumeAround(gp);
            dbj.beProductOf(d, bj);
            if ( matStiffSymmFlag ) {
                m->plusProductSymmUpper(bj, dbj, dV);
            } else {
                m->plusProductUnsym(bj, dbj, dV);
            }
        }

        if ( matStiffSymmFlag ) {
            m->symmetrized();
        }

        // localize irule contribution into element matrix
        if ( this->giveIntegrationElementLocalCodeNumbers(irlocnum, elem, iRule, EID_MomentumBalance) ) {
            answer.assemble(* m, irlocnum);
        }
    } // end loop over irules
}

} // end namespace oofem
