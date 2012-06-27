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
#include "structuralmaterial.h"
#include "structuralms.h"
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


void StructuralElementEvaluator :: giveCharacteristicVector(FloatArray &answer, CharType type, ValueModeType mode, TimeStep *tStep)
{
    if ( type == InternalForcesVector ) {
        this->giveInternalForcesVector(answer, tStep, false); /// @todo Only for total value mode type (?)
    } else if ( type == LastEquilibratedInternalForcesVector ) {
        this->giveInternalForcesVector(answer, tStep, true); /// @todo Only for total value mode type (?)
    } else {
        answer.resize(0);
    }
}


void StructuralElementEvaluator :: giveCharacteristicMatrix(FloatMatrix &answer,
                                                            CharType mtrx, TimeStep *tStep)
//
// returns characteristics matrix of receiver according to mtrx
//
{
    if ( mtrx == StiffnessMatrix ) {
        this->computeStiffnessMatrix(answer, TangentStiffness, tStep);
    } else if ( mtrx == MassMatrix ) {
        double mass;
        this->computeConsistentMassMatrix(answer, tStep, mass);
    } else if ( mtrx == LumpedMassMatrix ) {
        this->computeLumpedMassMatrix(answer, tStep);
    } else {
        OOFEM_ERROR2( "giveCharacteristicMatrix: Unknown Type of characteristic mtrx (%s)", __CharTypeToString(mtrx) );
    }
}


void StructuralElementEvaluator :: computeLumpedMassMatrix(FloatMatrix &answer, TimeStep *tStep)
// Returns the lumped mass matrix of the receiver.
{
    double mass;
    Element *elem = this->giveElement();
    int numberOfDofMans = elem->giveNumberOfDofManagers();
    IntArray nodeDofIDMask, dimFlag(3);
    int i, j, indx = 0, k, ldofs, dim;
    double summ;

    if ( !this->isActivated(tStep) ) {
        int ndofs = elem->computeNumberOfDofs(EID_MomentumBalance);
        answer.resize(ndofs, ndofs);
        answer.zero();
        return;
    }

    this->computeConsistentMassMatrix(answer, tStep, mass);
    ldofs = answer.giveNumberOfRows();

    for ( i = 1; i <= numberOfDofMans; i++ ) {
        elem->giveDofManDofIDMask(i, EID_MomentumBalance, nodeDofIDMask);
        for ( j = 1; j <= nodeDofIDMask.giveSize(); j++ ) {
            indx++;
            // zero all off-diagonal terms
            for ( k = 1; k <= ldofs; k++ ) {
                if ( k != indx ) {
                    answer.at(indx, k) = 0.;
                    answer.at(k, indx) = 0.;
                }
            }

            if ( ( nodeDofIDMask.at(j) != D_u ) && ( nodeDofIDMask.at(j) != D_v ) && ( nodeDofIDMask.at(j) != D_w ) ) {
                // zero corresponding diagonal member too <= no displacement dof
                answer.at(indx, indx) = 0.;
            } else if ( nodeDofIDMask.at(j) == D_u ) {
                dimFlag.at(1) = 1;
            } else if ( nodeDofIDMask.at(j) == D_v ) {
                dimFlag.at(2) = 1;
            } else if ( nodeDofIDMask.at(j) == D_w ) {
                dimFlag.at(3) = 1;
            }
        }
    }

    if ( indx != ldofs ) {
        OOFEM_ERROR("computeMassMatrix : internal consistency check failed");
    }

    dim = dimFlag.at(1) + dimFlag.at(2) + dimFlag.at(3);
    for ( summ = 0., k = 1; k <= ldofs; k++ ) {
        summ += answer.at(k, k);
    }

    answer.times(dim * mass / summ);
}


void StructuralElementEvaluator :: computeConsistentMassMatrix(FloatMatrix &answer, TimeStep *tStep, double &mass)
// Computes numerically the consistent (full) mass matrix of the receiver.
{
    Element *elem = this->giveElement();
    int ndofs = elem->computeNumberOfDofs(EID_MomentumBalance);
    double density, dV;
    FloatMatrix n;
    GaussPoint *gp;
    IntegrationRule *iRule;
    IntArray mask;

    answer.resize(ndofs, ndofs);
    answer.zero();
    if ( !this->isActivated(tStep) ) {
        return;
    }

    if ( ( iRule = this->giveMassMtrxIntegrationRule() ) ) {
        OOFEM_ERROR("computeConsistentMassMatrix no integration rule available");
    }

    this->giveMassMtrxIntegrationMask(mask);

    mass = 0.;

    for ( int ip = 0; ip < iRule->getNumberOfIntegrationPoints(); ip++ ) {
        gp      = iRule->getIntegrationPoint(ip);
        density = elem->giveMaterial()->give('d', gp);
        dV      = this->computeVolumeAround(gp);
        mass   += density * dV;
        this->computeNMatrixAt(n, gp);

        if ( mask.isEmpty() ) {
            answer.plusProductSymmUpper(n, n, density * dV);
        } else {
            double summ;

            for ( int i = 1; i <= ndofs; i++ ) {
                for ( int j = i; j <= ndofs; j++ ) {
                    summ = 0.;
                    for ( int k = 1; k <= n.giveNumberOfRows(); k++ ) {
                        if ( mask.at(k) == 0 ) {
                            continue;
                        }

                        summ += n.at(k, i) * n.at(k, j);
                    }

                    answer.at(i, j) += summ * density * dV;
                }
            }
        }
    }

    answer.symmetrized();
}


void StructuralElementEvaluator :: giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, bool useUpdatedGpRecord)
{
    Element *elem = this->giveElement();
    StructuralCrossSection *cs = ( StructuralCrossSection * ) elem->giveCrossSection();
    GaussPoint *gp;
    Material *mat = elem->giveMaterial();
    IntegrationRule *iRule;
    int ndofs = elem->computeNumberOfDofs(EID_MomentumBalance);
    FloatMatrix b;
    FloatArray bs, strain, stress, u, temp;
    IntArray irlocnum;
    double dV;

    elem->computeVectorOf(EID_MomentumBalance, VM_Total, tStep, u);

    answer.resize(ndofs);
    FloatArray *m = & answer;
    if ( elem->giveInterpolation()->hasSubPatchFormulation() ) {
        m = & temp;
    }

    int numberOfIntegrationRules = elem->giveNumberOfIntegrationRules();
    // loop over individual integration rules
    for ( int ir = 0; ir < numberOfIntegrationRules; ir++ ) {
        m->resize(0);
        iRule = elem->giveIntegrationRule(ir);
        for ( int i = 0; i < iRule->getNumberOfIntegrationPoints(); i++ ) {
            gp = iRule->getIntegrationPoint(i);
            this->computeBMatrixAt(b, gp);
            if ( useUpdatedGpRecord ) {
                stress = ( ( StructuralMaterialStatus * ) mat->giveStatus(gp) )->giveStressVector();
            } else {
                this->computeStrainVector(strain, gp, tStep, u); ///@todo This part computes the B matrix again; Inefficient.
                cs->giveRealStresses(stress, ReducedForm, gp, strain, tStep);
            }

            if ( stress.giveSize() == 0 ) {
                break;
            }

            // compute nodal representation of internal forces using f = B^T*Sigma dV
            dV  = this->computeVolumeAround(gp);
            bs.beTProductOf(b, stress);
            m->add(dV, bs);
        }
        // localize irule contribution into element matrix
        if ( this->giveIntegrationElementLocalCodeNumbers(irlocnum, elem, iRule, EID_MomentumBalance) ) {
            answer.assemble(* m, irlocnum);
        }
    } // end loop over irules

    // if inactive update state, but no contribution to global system
    if ( !this->isActivated(tStep) ) {
        answer.zero();
    }
}


void StructuralElementEvaluator :: computeStrainVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep, FloatArray &u)
// Computes the vector containing the strains at the Gauss point gp of
// the receiver, at time step tStep. The nature of these strains depends
// on the element's type.
{
    int i;
    FloatMatrix b;
    FloatArray ur;
    Element *elem = this->giveElement();

    if ( !this->isActivated(tStep) ) {
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

void StructuralElementEvaluator :: updateInternalState(TimeStep *tStep)
// Updates the receiver at end of step.
{
    FloatArray u;
    Element *elem = this->giveElement();
    StructuralCrossSection *cs = ( StructuralCrossSection * ) elem->giveCrossSection();

    elem->computeVectorOf(EID_MomentumBalance, VM_Total, tStep, u);

#if 0
     // subtract initial displacements, if defined
     if (initialDisplacements) u.subtract(initialDisplacements);
#endif

    int i, j;
    IntegrationRule *iRule;
    FloatArray strain, stress;
    GaussPoint *gp;

    // force updating strains & stresses
    for ( i = 0; i < elem->giveNumberOfIntegrationRules(); i++ ) {
#ifdef __PARALLEL_MODE
        if (this->giveElement()->giveKnotSpanParallelMode(i) == Element_remote) continue;
#endif
        iRule = elem->giveIntegrationRule(i);
        for ( j = 0; j < iRule->getNumberOfIntegrationPoints(); j++ ) {
            gp = iRule->getIntegrationPoint(j);
            this->computeStrainVector(strain, gp, tStep, u);
            cs->giveRealStresses(stress, ReducedForm, gp, strain, tStep);
        }
    }

#if 0
     // Original unoptimized version
     int i, j;
     IntegrationRule *iRule;
     FloatArray stress;
     Element *elem = this->giveElement();
     // force updating strains & stresses
     for ( i = 0; i < elem->giveNumberOfIntegrationRules(); i++ ) {
     #ifdef __PARALLEL_MODE
         if (this->giveElement()->giveKnotSpanParallelMode(i) == Element_remote) continue;
     #endif
        iRule = elem->giveIntegrationRule(i);
        for ( j = 0; j < iRule->getNumberOfIntegrationPoints(); j++ ) {
            computeStressVector(stress, iRule->getIntegrationPoint(j), tStep);
        }
     }
#endif
}


void StructuralElementEvaluator :: computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
{
    int ir, j, numberOfIntegrationRules;
    FloatMatrix temp, bj, d, dbj;
    IntegrationRule *iRule;
    GaussPoint *gp;
    Element *elem = this->giveElement();
    StructuralCrossSection *cs = ( StructuralCrossSection * ) elem->giveCrossSection();
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
            cs->giveCharMaterialStiffnessMatrix(d, rMode, gp, tStep);

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
