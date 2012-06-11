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

#include "libeam3dnl2.h"
#include "node.h"
#include "material.h"
#include "crosssection.h"
#include "gausspnt.h"
#include "gaussintegrationrule.h"
#include "structuralms.h"
#include "flotmtrx.h"
#include "intarray.h"
#include "flotarry.h"
#include "mathfem.h"
#include "timestep.h"
#include "contextioerr.h"

#ifdef __OOFEG
 #include "oofeggraphiccontext.h"
#endif

namespace oofem {
LIBeam3dNL2 :: LIBeam3dNL2(int n, Domain *aDomain) : NLStructuralElement(n, aDomain), q(4), tempQ(4) //, kappa (3)
{
    numberOfDofMans    = 2;
    l0                 = 0.;
    tempQCounter       = 0;
    referenceNode      = 0;
    // init kappa vector at centre
    // kappa.zero();
}


void
LIBeam3dNL2 :: computeSMtrx(FloatMatrix &answer, FloatArray &vec)
{
    if ( vec.giveSize() != 3 ) {
        _error("computeSMtrx: vec param size mismatch");
    }

    answer.resize(3, 3);

    answer.at(1, 1) = answer.at(2, 2) = answer.at(3, 3) = 0.;
    answer.at(1, 2) = -vec.at(3);
    answer.at(1, 3) =  vec.at(2);
    answer.at(2, 1) =  vec.at(3);
    answer.at(2, 3) = -vec.at(1);
    answer.at(3, 1) = -vec.at(2);
    answer.at(3, 2) =  vec.at(1);
}


void
LIBeam3dNL2 :: computeRotMtrx(FloatMatrix &answer, FloatArray &psi)
{
    FloatMatrix S(3, 3), SS(3, 3);
    double psiSize;

    if ( psi.giveSize() != 3 ) {
        _error("computeSMtrx: psi param size mismatch");
    }

    answer.resize(3, 3);
    answer.zero();

    psiSize = psi.computeNorm();
    answer.at(1, 1) = answer.at(2, 2) = answer.at(3, 3) = 1.;

    if ( psiSize <= 1.e-40 ) {
        return;
    }

    this->computeSMtrx(S, psi);
    SS.beProductOf(S, S);
    S.times(sin(psiSize) / psiSize);
    SS.times( ( 1. - cos(psiSize) ) / ( psiSize * psiSize ) );

    answer.add(S);
    answer.add(SS);
}


void
LIBeam3dNL2 :: updateTempQuaternion(TimeStep *tStep)
{
    // test if not previously updated temporary quaternion
    if ( tStep->giveSolutionStateCounter() != tempQCounter ) {
        // update temporary quaternion
        FloatArray u, centreSpin(3), q2(4);
        FloatMatrix dR(3, 3);
        double centreSpinSize;

        // ask element's displacement increments
        this->computeVectorOf(EID_MomentumBalance, VM_Incremental, tStep, u);

        // interpolate spin at the centre
        centreSpin.at(1) = 0.5 * ( u.at(4) + u.at(10) );
        centreSpin.at(2) = 0.5 * ( u.at(5) + u.at(11) );
        centreSpin.at(3) = 0.5 * ( u.at(6) + u.at(12) );

        centreSpinSize = centreSpin.computeNorm();
        if ( centreSpinSize > 1.e-30 ) {
            centreSpin.normalize();
            q2.at(1) = sin(centreSpinSize / 2.) * centreSpin.at(1);
            q2.at(2) = sin(centreSpinSize / 2.) * centreSpin.at(2);
            q2.at(3) = sin(centreSpinSize / 2.) * centreSpin.at(3);
            q2.at(4) = cos(centreSpinSize / 2.);

            // update temporary quaternion at the center
            tempQ.at(1) = q.at(4) * q2.at(1) + q2.at(4) * q.at(1) - ( q.at(2) * q2.at(3) - q.at(3) * q2.at(2) );
            tempQ.at(2) = q.at(4) * q2.at(2) + q2.at(4) * q.at(2) - ( q.at(3) * q2.at(1) - q.at(1) * q2.at(3) );
            tempQ.at(3) = q.at(4) * q2.at(3) + q2.at(4) * q.at(3) - ( q.at(1) * q2.at(2) - q.at(2) * q2.at(1) );
            tempQ.at(4) = q.at(4) * q2.at(4) - ( q.at(1) * q2.at(1) + q.at(2) * q2.at(2) + q.at(3) * q2.at(3) );
        } else {
            tempQ = q;
        }

        // remember timestamp
        tempQCounter = tStep->giveSolutionStateCounter();
    }
}


void
LIBeam3dNL2 :: computeRotMtrxFromQuaternion(FloatMatrix &answer, FloatArray &q)
{
    answer.resize(3, 3);

    answer.at(1, 1) = q.at(4) * q.at(4) + q.at(1) * q.at(1) - 0.5;
    answer.at(1, 2) = q.at(1) * q.at(2) - q.at(3) * q.at(4);
    answer.at(1, 3) = q.at(1) * q.at(3) + q.at(2) * q.at(4);

    answer.at(2, 1) = q.at(2) * q.at(1) + q.at(3) * q.at(4);
    answer.at(2, 2) = q.at(4) * q.at(4) + q.at(2) * q.at(2) - 0.5;
    answer.at(2, 3) = q.at(2) * q.at(3) - q.at(1) * q.at(4);

    answer.at(3, 1) = q.at(3) * q.at(1) - q.at(2) * q.at(4);
    answer.at(3, 2) = q.at(3) * q.at(2) + q.at(1) * q.at(4);
    answer.at(3, 3) = q.at(4) * q.at(4) + q.at(3) * q.at(3) - 0.5;

    answer.times(2.);
}


void
LIBeam3dNL2 :: computeQuaternionFromRotMtrx(FloatArray &answer, FloatMatrix &R)
{
    // Spurrier's algorithm

    int i, ii;
    double a, trR;

    answer.resize(4);

    trR = R.at(1, 1) + R.at(2, 2) + R.at(3, 3);
    a = trR;
    ii = 0;
    for ( i = 1; i <= 3; i++ ) {
        if ( R.at(i, i) > a ) {
            a = R.at(i, i);
            ii = i;
        }
    }

    if ( a == trR ) {
        //printf (".");
        answer.at(4) = 0.5 * sqrt(1. + a);
        answer.at(1) = ( R.at(3, 2) - R.at(2, 3) ) / ( 4. * answer.at(4) );
        answer.at(2) = ( R.at(1, 3) - R.at(3, 1) ) / ( 4. * answer.at(4) );
        answer.at(3) = ( R.at(2, 1) - R.at(1, 2) ) / ( 4. * answer.at(4) );
    } else {
        //printf (":");
        int jj, kk;
        if ( ii == 1 ) {
            jj = 2;
            kk = 3;
        } else if ( ii == 2 ) {
            jj = 3;
            kk = 1;
        } else {
            jj = 1;
            kk = 2;
        }

        answer.at(ii) = sqrt( 0.5 * a + 0.25 * ( 1. - trR ) );
        answer.at(4) = 0.25 * ( R.at(kk, jj) - R.at(jj, kk) ) / answer.at(ii);
        answer.at(jj) = 0.25 * ( R.at(jj, ii) + R.at(ii, jj) ) / answer.at(ii);
        answer.at(kk) = 0.25 * ( R.at(kk, ii) + R.at(ii, kk) ) / answer.at(ii);
    }
}


void
LIBeam3dNL2 :: computeStrainVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep)
{
    FloatArray ui(3), xd(3), eps(3), curv(3), ac(3);
    FloatMatrix sc(3, 3), tmid(3, 3), tempTc;
    IntArray mask(1);

    // update temp triad
    this->updateTempQuaternion(tStep);
    this->computeRotMtrxFromQuaternion(tempTc, this->tempQ);

    this->computeXdVector(xd, tStep);

    // compute eps_xl, gamma_l2, gamma_l3
    eps.beTProductOf(tempTc, xd);
    eps.times(1. / this->l0);
    eps.at(1) -= 1.0;

    // update curvature at midpoint
    this->computeTempCurv(curv, tStep);

    answer.resize(6);
    answer.at(1) = eps.at(1); // eps_xl
    answer.at(2) = eps.at(2); // gamma_l2
    answer.at(3) = eps.at(3); // gamma_l3
    answer.at(4) = curv.at(1); // kappa_1
    answer.at(5) = curv.at(2); // kappa_2
    answer.at(6) = curv.at(3); // kappa_3
}


void
LIBeam3dNL2 :: computeXMtrx(FloatMatrix &answer, TimeStep *tStep)
{
    int i, j;
    FloatArray xd(3);
    FloatMatrix s(3, 3);

    this->computeXdVector(xd, tStep);
    this->computeSMtrx(s, xd);

    answer.resize(12, 6);
    answer.zero();

    for ( i = 1; i < 4; i++ ) {
        answer.at(i, i)      = -1.0;
        answer.at(i + 6, i)   =  1.0;
        answer.at(i + 3, i + 3) = -1.0;
        answer.at(i + 9, i + 3) =  1.0;

        for ( j = 1; j < 4; j++ ) {
            answer.at(i + 3, j) = answer.at(i + 9, j) = 0.5 * s.at(j, i);
        }
    }
}


void
LIBeam3dNL2 :: giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord)
{
    int i, j;
    Material *mat = this->giveMaterial();
    IntegrationRule *iRule = integrationRulesArray [ giveDefaultIntegrationRule() ];
    GaussPoint *gp = iRule->getIntegrationPoint(0);
    FloatArray nm(6), TotalStressVector(6);
    FloatMatrix x, tempTc;
    double s1, s2;

    // update temp triad
    this->updateTempQuaternion(tStep);
    this->computeRotMtrxFromQuaternion(tempTc, this->tempQ);

    if ( useUpdatedGpRecord == 1 ) {
        TotalStressVector = ( ( StructuralMaterialStatus * ) mat->giveStatus(gp) )
                            ->giveStressVector();
    } else {
        this->computeStressVector(TotalStressVector, gp, tStep);
    }

    for ( i = 1; i <= 3; i++ ) {
        s1 = s2 = 0.0;
        for ( j = 1; j <= 3; j++ ) {
            s1 += tempTc.at(i, j) * TotalStressVector.at(j);
            s2 += tempTc.at(i, j) * TotalStressVector.at(j + 3);
        }

        nm.at(i)   = s1;
        nm.at(i + 3) = s2;
    }

    this->computeXMtrx(x, tStep);
    answer.beProductOf(x, nm);
}


void
LIBeam3dNL2 :: computeXdVector(FloatArray &answer, TimeStep *tStep)
{
    FloatArray u(3);

    answer.resize(3);
    // ask element's displacements
    this->computeVectorOf(EID_MomentumBalance, VM_Total, tStep, u);

    answer.at(1) = ( this->giveNode(2)->giveCoordinate(1) + u.at(7) ) -
                   ( this->giveNode(1)->giveCoordinate(1) + u.at(1) );
    answer.at(2) = ( this->giveNode(2)->giveCoordinate(2) + u.at(8) ) -
                   ( this->giveNode(1)->giveCoordinate(2) + u.at(2) );
    answer.at(3) = ( this->giveNode(2)->giveCoordinate(3) + u.at(9) ) -
                   ( this->giveNode(1)->giveCoordinate(3) + u.at(3) );
}


void
LIBeam3dNL2 :: computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
{
    int i, j, k;
    double s1, s2;
    FloatMatrix d, x, xt(12, 6), dxt, sn, sm, sxd, y, tempTc;
    FloatArray n(3), m(3), xd(3), TotalStressVector;
    IntegrationRule *iRule = integrationRulesArray [ giveDefaultIntegrationRule() ];
    GaussPoint *gp = iRule->getIntegrationPoint(0);

    answer.resize( this->computeNumberOfDofs(EID_MomentumBalance), this->computeNumberOfDofs(EID_MomentumBalance) );
    answer.zero();

    // linear part

    this->updateTempQuaternion(tStep);
    this->computeRotMtrxFromQuaternion(tempTc, this->tempQ);
    this->computeXMtrx(x, tStep);
    xt.zero();
    for ( i = 1; i <= 12; i++ ) {
        for ( j = 1; j <= 3; j++ ) {
            for ( k = 1; k <= 3; k++ ) {
                // compute x*Tbar, taking into account sparsity of Tbar
                xt.at(i, j)   += x.at(i, k) * tempTc.at(k, j);
                xt.at(i, j + 3) += x.at(i, k + 3) * tempTc.at(k, j);
            }
        }
    }

    this->computeConstitutiveMatrixAt(d, rMode, gp, tStep);
    dxt.beProductTOf(d, xt);
    answer.beProductOf(xt, dxt);
    answer.times(1. / this->l0);

    // geometric stiffness ks = ks1+ks2
    // ks1
    this->computeStressVector(TotalStressVector, gp, tStep);

    for ( i = 1; i <= 3; i++ ) {
        s1 = s2 = 0.0;
        for ( j = 1; j <= 3; j++ ) {
            s1 += tempTc.at(i, j) * TotalStressVector.at(j);
            s2 += tempTc.at(i, j) * TotalStressVector.at(j + 3);
        }

        n.at(i)   = s1;
        m.at(i)   = s2;
    }

    this->computeSMtrx(sn, n);
    this->computeSMtrx(sm, m);

    for ( i = 1; i <= 3; i++ ) {
        for ( j = 1; j <= 3; j++ ) {
            answer.at(i, j + 3)   += sn.at(i, j);
            answer.at(i, j + 9)   += sn.at(i, j);
            answer.at(i + 3, j + 3) += sm.at(i, j);
            answer.at(i + 3, j + 9) += sm.at(i, j);

            answer.at(i + 6, j + 3) -= sn.at(i, j);
            answer.at(i + 6, j + 9) -= sn.at(i, j);
            answer.at(i + 9, j + 3) -= sm.at(i, j);
            answer.at(i + 9, j + 9) -= sm.at(i, j);
        }
    }

    // ks2
    this->computeXdVector(xd, tStep);
    this->computeSMtrx(sxd, xd);

    y.beProductOf(sxd, sn);
    y.times(0.5);

    for ( i = 1; i <= 3; i++ ) {
        for ( j = 1; j <= 3; j++ ) {
            answer.at(i + 3, j)     -= sn.at(i, j);
            answer.at(i + 3, j + 3)   += y.at(i, j);
            answer.at(i + 3, j + 6)   += sn.at(i, j);
            answer.at(i + 3, j + 9)   += y.at(i, j);

            answer.at(i + 9, j)     -= sn.at(i, j);
            answer.at(i + 9, j + 3)   += y.at(i, j);
            answer.at(i + 9, j + 6)   += sn.at(i, j);
            answer.at(i + 9, j + 9)   += y.at(i, j);
        }
    }
}


void
LIBeam3dNL2 :: computeGaussPoints()
// Sets up the array of Gauss Points of the receiver.
{
    if ( !integrationRulesArray ) {
        numberOfIntegrationRules = 1;
        integrationRulesArray = new IntegrationRule * [ 1 ];
        integrationRulesArray [ 0 ] = new GaussIntegrationRule(1, this, 1, 2);
        integrationRulesArray [ 0 ]->setUpIntegrationPoints(_Line, 1, _3dBeam);
    }
}

IRResultType
LIBeam3dNL2 :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    // first call parent
    NLStructuralElement :: initializeFrom(ir);

    IR_GIVE_FIELD(ir, referenceNode, IFT_LIBeam3dNL2_refnode, "refnode"); // Macro
    if ( referenceNode == 0 ) {
        _error("instanciateFrom: wrong reference node specified");
    }

    /*
     * if (this->hasString (initString, "dofstocondense")) {
     *  dofsToCondense = this->ReadIntArray (initString, "dofstocondense");
     *  if (dofsToCondense->giveSize() >= 12)
     *    _error ("instanciateFrom: wrong input data for condensed dofs");
     * } else {
     *  dofsToCondense = NULL;
     * }
     */

    // compute initial triad at centre - requires nodal coordinates
    FloatMatrix lcs, tc;
    this->giveLocalCoordinateSystem(lcs);
    tc.beTranspositionOf(lcs);

    this->computeQuaternionFromRotMtrx(q, tc);
    this->computeGaussPoints();
    return IRRT_OK;
}


double
LIBeam3dNL2 :: giveLength()
// Returns the original length (l0) of the receiver.
{
    double dx, dy, dz;
    Node *nodeA, *nodeB;

    if ( l0 == 0. ) {
        nodeA   = this->giveNode(1);
        nodeB   = this->giveNode(2);
        dx      = nodeB->giveCoordinate(1) - nodeA->giveCoordinate(1);
        dy      = nodeB->giveCoordinate(2) - nodeA->giveCoordinate(2);
        dz      = nodeB->giveCoordinate(3) - nodeA->giveCoordinate(3);
        l0      = sqrt(dx * dx + dy * dy + dz * dz);
    }

    return l0;
}


void
LIBeam3dNL2 :: computeLumpedMassMatrix(FloatMatrix &answer, TimeStep *tStep)
// Returns the lumped mass matrix of the receiver. This expression is
// valid in both local and global axes.
{
    Material *mat;
    double halfMass;
    GaussPoint *gp = integrationRulesArray [ 0 ]->getIntegrationPoint(0);

    mat        = this->giveMaterial();
    halfMass   = mat->give('d', gp) * this->giveCrossSection()->give(CS_Area) * this->giveLength() / 2.;
    answer.resize(12, 12);
    answer.zero();
    answer.at(1, 1) = answer.at(2, 2) = answer.at(3, 3) = halfMass;
    answer.at(7, 7) = answer.at(8, 8) = answer.at(9, 9) = halfMass;
}


void
LIBeam3dNL2 :: computeNmatrixAt(GaussPoint *aGaussPoint, FloatMatrix &answer)
// Returns the displacement interpolation matrix {N} of the receiver, eva-
// luated at aGaussPoint.
{
    double ksi, n1, n2;

    ksi = aGaussPoint->giveCoordinate(1);
    n1  = ( 1. - ksi ) * 0.5;
    n2  = ( 1. + ksi ) * 0.5;

    answer.resize(6, 12);
    answer.zero();

    // u
    answer.at(1, 1) = n1;
    answer.at(1, 7) = n2;
    // v
    answer.at(2, 2) = n1;
    answer.at(2, 8) = n2;
    // w
    answer.at(3, 3) = n1;
    answer.at(3, 9) = n2;
    // fi_x
    answer.at(4, 4)  = n1;
    answer.at(4, 10) = n2;
    // fi_y
    answer.at(5, 5)  = n1;
    answer.at(5, 11) = n2;
    // fi_z
    answer.at(6, 6)  = n1;
    answer.at(6, 12) = n2;
}


double
LIBeam3dNL2 :: computeVolumeAround(GaussPoint *aGaussPoint)
// Returns the length of the receiver. This method is valid only if 1
// Gauss point is used.
{
    double weight  = aGaussPoint->giveWeight();
    return weight * 0.5 * this->giveLength();
}


void
LIBeam3dNL2 ::   giveDofManDofIDMask(int inode, EquationID, IntArray &answer) const
{
    answer.setValues(6, D_u, D_v, D_w, R_u, R_v, R_w);
}

int
LIBeam3dNL2 :: computeGlobalCoordinates(FloatArray &answer, const FloatArray &lcoords)
{
    double ksi, n1, n2;

    ksi = lcoords.at(1);
    n1  = ( 1. - ksi ) * 0.5;
    n2  = ( 1. + ksi ) * 0.5;

    answer.resize(3);
    answer.at(1) = n1 * this->giveNode(1)->giveCoordinate(1) + n2 *this->giveNode(2)->giveCoordinate(1);
    answer.at(2) = n1 * this->giveNode(1)->giveCoordinate(2) + n2 *this->giveNode(2)->giveCoordinate(2);
    answer.at(3) = n1 * this->giveNode(1)->giveCoordinate(3) + n2 *this->giveNode(2)->giveCoordinate(3);

    return 1;
}


void
LIBeam3dNL2 :: computeEgdeNMatrixAt(FloatMatrix &answer, GaussPoint *aGaussPoint)
{
    /*
     *
     * computes interpolation matrix for element edge.
     * we assemble locally this matrix for only nonzero
     * shape functions.
     * (for example only two nonzero shape functions for 2 dofs are
     * necessary for linear plane stress tringle edge).
     * These nonzero shape functions are then mapped to
     * global element functions.
     *
     * Uso of mapping technique will allow to assemble shape functions
     * without regarding particular side
     */

    this->computeNmatrixAt(aGaussPoint, answer);
}


void
LIBeam3dNL2 :: giveEdgeDofMapping(IntArray &answer, int iEdge) const
{
    /*
     * provides dof mapping of local edge dofs (only nonzero are taken into account)
     * to global element dofs
     */
    int i;

    if ( iEdge != 1 ) {
        _error("giveEdgeDofMapping: wrong edge number");
    }

    answer.resize(12);
    for ( i = 1; i <= 12; i++ ) {
        answer.at(i) = i;
    }
}


double
LIBeam3dNL2 :: computeEdgeVolumeAround(GaussPoint *aGaussPoint, int iEdge)
{
    if ( iEdge != 1 ) { // edge between nodes 1 2
        _error("computeEdgeVolumeAround: wrong egde number");
    }

    double weight  = aGaussPoint->giveWeight();
    return 0.5 * this->giveLength() * weight;
}


int
LIBeam3dNL2 :: giveLocalCoordinateSystem(FloatMatrix &answer)
//
// returns a unit vectors of local coordinate system at element
// stored rowwise (mainly used by some materials with ortho and anisotrophy)
//
{
    FloatArray lx(3), ly(3), lz(3), help(3);
    double length = this->giveLength();
    Node *nodeA, *nodeB, *refNode;
    int i;

    answer.resize(3, 3);
    answer.zero();
    nodeA  = this->giveNode(1);
    nodeB  = this->giveNode(2);
    refNode = this->giveDomain()->giveNode(this->referenceNode);

    for ( i = 1; i <= 3; i++ ) {
        lx.at(i) = ( nodeB->giveCoordinate(i) - nodeA->giveCoordinate(i) ) / length;
        help.at(i) = ( refNode->giveCoordinate(i) - nodeA->giveCoordinate(i) );
    }

    lz.beVectorProductOf(lx, help);
    lz.normalize();
    ly.beVectorProductOf(lz, lx);
    ly.normalize();

    for ( i = 1; i <= 3; i++ ) {
        answer.at(1, i) = lx.at(i);
        answer.at(2, i) = ly.at(i);
        answer.at(3, i) = lz.at(i);
    }

    return 1;
}


/*
 * int
 * LIBeam3dNL2 :: computeGtoLRotationMatrix (FloatMatrix& answer) // giveRotationMatrix ()
 * // Returns the rotation matrix of the receiver.
 * // rotation matrix is computed for original configuration only
 * {
 * FloatMatrix lcs;
 * int i,j;
 *
 * answer.resize(12,12);
 * answer.zero();
 *
 * this->giveLocalCoordinateSystem (lcs);
 *
 * for (i=1; i <= 3; i++)
 * for (j=1; j <= 3; j++) {
 * answer.at(i,j) = lcs.at(i,j);
 * answer.at(i+3, j+3) = lcs.at(i,j);
 * answer.at(i+6, j+6) = lcs.at(i,j);
 * answer.at(i+9, j+9) = lcs.at(i,j);
 * }
 *
 * //delete lcs;
 * return 1 ;
 * }
 */

int
LIBeam3dNL2 :: computeLoadGToLRotationMtrx(FloatMatrix &answer)
{
    /*
     * Returns transformation matrix from global coordinate system to local
     * element coordinate system for element load vector components.
     * If no transformation is necessary, answer is empty matrix (default);
     *
     * Does not support follower load
     */

    FloatMatrix lcs;
    int i, j;

    answer.resize(6, 6);
    answer.zero();

    this->giveLocalCoordinateSystem(lcs);

    for ( i = 1; i <= 3; i++ ) {
        for ( j = 1; j <= 3; j++ ) {
            answer.at(i, j) = lcs.at(i, j);
            answer.at(3 + i, 3 + j) = lcs.at(i, j);
        }
    }

    return 1;
}


int
LIBeam3dNL2 :: computeLoadLEToLRotationMatrix(FloatMatrix &answer, int iEdge, GaussPoint *gp)
{
    // returns transformation matrix from
    // edge local coordinate system
    // to element local c.s
    // (same as global c.s in this case)
    //
    // i.e. f(element local) = T * f(edge local)
    //
    answer.beEmptyMtrx();
    return 0;
}


void
LIBeam3dNL2 :: computeBodyLoadVectorAt(FloatArray &answer, Load *load, TimeStep *tStep, ValueModeType mode)
{
    NLStructuralElement::computeBodyLoadVectorAt(answer, load, tStep, mode);
    answer.times(this->giveCrossSection()->give(CS_Area));
}


void
LIBeam3dNL2 :: updateYourself(TimeStep *tStep)
// Updates the receiver at end of step.
{
    NLStructuralElement :: updateYourself(tStep);

    // update quaternion
    this->updateTempQuaternion(tStep);
    q = tempQ;

    // update curvature
    //FloatArray curv;
    //this->computeTempCurv (curv, tStep);
    //kappa = curv;
}

void
LIBeam3dNL2 :: initForNewStep()
// initializes receiver to new time step or can be used
// if current time step must be restarted
//
// call material->initGpForNewStep() for all GPs.
//
{
    NLStructuralElement :: initForNewStep();
    tempQ = q;
}


void
LIBeam3dNL2 :: computeTempCurv(FloatArray &answer, TimeStep *tStep)
{
    Material *mat = this->giveMaterial();
    IntegrationRule *iRule = integrationRulesArray [ giveDefaultIntegrationRule() ];
    GaussPoint *gp = iRule->getIntegrationPoint(0);
    FloatArray ui(3), xd(3), curv(3), ac(3), PrevEpsilon;
    FloatMatrix sc(3, 3), tmid(3, 3), tc(3, 3);

    answer.resize(3);

    /*
     * // update curvature at midpoint
     * // first, compute Tmid
     * // ask increments
     * this -> computeVectorOf(EID_MomentumBalance,UnknownMode_Incremental,tStep, ui) ;
     * ac.at(1) = 0.5*(ui.at(10) - ui.at(4));
     * ac.at(2) = 0.5*(ui.at(11) - ui.at(5));
     * ac.at(3) = 0.5*(ui.at(12) - ui.at(6));
     * this->computeSMtrx (sc, ac);
     * sc.times (1./2.);
     * // compute I+sc
     * sc.at(1,1) += 1.0;
     * sc.at(2,2) += 1.0;
     * sc.at(3,3) += 1.0;
     *
     * this->computeRotMtrxFromQuaternion(tc, this->q);
     * tmid.beProductOf (sc, tc);
     *
     * // update curvature at centre
     * ac.at(1) = (ui.at(10) - ui.at(4));
     * ac.at(2) = (ui.at(11) - ui.at(5));
     * ac.at(3) = (ui.at(12) - ui.at(6));
     *
     * answer.beTProductOf (tmid, ac);
     * answer.times (1/this->l0);
     * // ask for previous kappa
     * PrevEpsilon = ((StructuralMaterialStatus*) mat->giveStatus(gp)) -> giveStrainVector ();
     * if (PrevEpsilon.giveSize()) {
     * answer.at(1) += PrevEpsilon.at(4);
     * answer.at(2) += PrevEpsilon.at(5);
     * answer.at(3) += PrevEpsilon.at(6);
     * }
     */
    // update curvature
    // exact procedure due to Simo & Vu Quoc
    FloatMatrix dR, Rn, Ro;
    FloatArray e, om, omp, acp(3), kapgn1(3);
    double acSize, coeff;
    int i, j;

    this->computeVectorOf(EID_MomentumBalance, VM_Incremental, tStep, ui);

    ac.at(1) = 0.5 * ( ui.at(10) - ui.at(4) );
    ac.at(2) = 0.5 * ( ui.at(11) - ui.at(5) );
    ac.at(3) = 0.5 * ( ui.at(12) - ui.at(6) );

    this->computeRotMtrx(dR, ac);
    this->computeRotMtrxFromQuaternion(Ro, this->q);
    Rn.beProductOf(dR, Ro);

    acSize = ac.computeSquaredNorm();

    if ( acSize > 1.e-30 ) {
        FloatMatrix h(3, 3);
        ac.normalize();
        om = ac;
        om.times( 2. * tan(acSize / 2.) );

        coeff = ( 1. - ( acSize / sin(acSize) ) );
        for ( i = 1; i <= 3; i++ ) {
            for ( j = 1; j <= 3; j++ ) {
                h.at(i, j) = coeff * ac.at(i) * ac.at(j);
            }
        }

        for ( i = 1; i <= 3; i++ ) {
            h.at(i, i) = 1. - h.at(i, i);
        }

        // compute dPsi/ds
        acp.at(1) = ( ui.at(10) - ui.at(4) ) / this->l0;
        acp.at(2) = ( ui.at(11) - ui.at(5) ) / this->l0;
        acp.at(3) = ( ui.at(12) - ui.at(6) ) / this->l0;

        omp.beProductOf(h, acp);
        omp.times(2. * tan(acSize / 2.) / acSize);

        kapgn1.beVectorProductOf(om, omp);
        kapgn1.times(1. / 2.);
        kapgn1.add(omp);

        coeff = 1. / ( 1. + 0.25 * om.computeSquaredNorm() );
        kapgn1.times(coeff);

        answer.beTProductOf(Rn, kapgn1);
    }

    // ask for previous kappa
    PrevEpsilon = ( ( StructuralMaterialStatus * ) mat->giveStatus(gp) )->giveStrainVector();
    if ( PrevEpsilon.giveSize() ) {
        answer.at(1) += PrevEpsilon.at(4);
        answer.at(2) += PrevEpsilon.at(5);
        answer.at(3) += PrevEpsilon.at(6);
    }
}


contextIOResultType
LIBeam3dNL2 :: saveContext(DataStream *stream, ContextMode mode, void *obj)
//
// saves full element context (saves state variables, that completely describe
// current state)
//
{
    contextIOResultType iores;
    if ( ( iores = NLStructuralElement :: saveContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = q.storeYourself(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    return CIO_OK;
}


contextIOResultType
LIBeam3dNL2 :: restoreContext(DataStream *stream, ContextMode mode, void *obj)
//
// restores full element context (saves state variables, that completely describe
// current state)
//
{
    contextIOResultType iores;
    if ( ( iores = NLStructuralElement :: restoreContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = q.restoreYourself(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    return CIO_OK;
}


#ifdef __OOFEG
void LIBeam3dNL2 :: drawRawGeometry(oofegGraphicContext &gc)
{
    GraphicObj *go;

    if ( !gc.testElementGraphicActivity(this) ) {
        return;
    }

    //  if (!go) { // create new one
    WCRec p [ 2 ];   /* poin */
    EASValsSetLineWidth(OOFEG_RAW_GEOMETRY_WIDTH);
    EASValsSetColor( gc.getElementColor() );
    EASValsSetLayer(OOFEG_RAW_GEOMETRY_LAYER);
    p [ 0 ].x = ( FPNum ) this->giveNode(1)->giveCoordinate(1);
    p [ 0 ].y = ( FPNum ) this->giveNode(1)->giveCoordinate(2);
    p [ 0 ].z = ( FPNum ) this->giveNode(1)->giveCoordinate(3);
    p [ 1 ].x = ( FPNum ) this->giveNode(2)->giveCoordinate(1);
    p [ 1 ].y = ( FPNum ) this->giveNode(2)->giveCoordinate(2);
    p [ 1 ].z = ( FPNum ) this->giveNode(2)->giveCoordinate(3);
    go = CreateLine3D(p);
    EGWithMaskChangeAttributes(WIDTH_MASK | COLOR_MASK | LAYER_MASK, go);
    EGAttachObject(go, ( EObjectP ) this);
    EMAddGraphicsToModel(ESIModel(), go);
}


void LIBeam3dNL2 :: drawDeformedGeometry(oofegGraphicContext &gc, UnknownType type)
{
    GraphicObj *go;

    if ( !gc.testElementGraphicActivity(this) ) {
        return;
    }

    TimeStep *tStep = domain->giveEngngModel()->giveCurrentStep();
    double defScale = gc.getDefScale();
    //  if (!go) { // create new one
    WCRec p [ 2 ]; /* poin */
    const char *colors[] = {
        "red", "green", "blue"
    };

    EASValsSetLineWidth(OOFEG_DEFORMED_GEOMETRY_WIDTH);
    EASValsSetColor( gc.getDeformedElementColor() );
    EASValsSetLayer(OOFEG_DEFORMED_GEOMETRY_LAYER);
    p [ 0 ].x = ( FPNum ) this->giveNode(1)->giveUpdatedCoordinate(1, tStep, EID_MomentumBalance, defScale);
    p [ 0 ].y = ( FPNum ) this->giveNode(1)->giveUpdatedCoordinate(2, tStep, EID_MomentumBalance, defScale);
    p [ 0 ].z = ( FPNum ) this->giveNode(1)->giveUpdatedCoordinate(3, tStep, EID_MomentumBalance, defScale);

    p [ 1 ].x = ( FPNum ) this->giveNode(2)->giveUpdatedCoordinate(1, tStep, EID_MomentumBalance, defScale);
    p [ 1 ].y = ( FPNum ) this->giveNode(2)->giveUpdatedCoordinate(2, tStep, EID_MomentumBalance, defScale);
    p [ 1 ].z = ( FPNum ) this->giveNode(2)->giveUpdatedCoordinate(3, tStep, EID_MomentumBalance, defScale);
    go = CreateLine3D(p);
    EGWithMaskChangeAttributes(WIDTH_MASK | COLOR_MASK | LAYER_MASK, go);
    EMAddGraphicsToModel(ESIModel(), go);

    // draw centre triad
    FloatMatrix tc;
    int i, succ;
    double coeff = this->l0 / 3.;
    this->computeRotMtrxFromQuaternion(tc, this->q);

    p [ 0 ].x = 0.5 * ( ( FPNum ) this->giveNode(1)->giveUpdatedCoordinate(1, tStep, EID_MomentumBalance, defScale) + ( FPNum ) this->giveNode(2)->giveUpdatedCoordinate(1, tStep, EID_MomentumBalance, defScale) );
    p [ 0 ].y = 0.5 * ( ( FPNum ) this->giveNode(1)->giveUpdatedCoordinate(2, tStep, EID_MomentumBalance, defScale) + ( FPNum ) this->giveNode(2)->giveUpdatedCoordinate(2, tStep, EID_MomentumBalance, defScale) );
    p [ 0 ].z = 0.5 * ( ( FPNum ) this->giveNode(1)->giveUpdatedCoordinate(3, tStep, EID_MomentumBalance, defScale) + ( FPNum ) this->giveNode(2)->giveUpdatedCoordinate(3, tStep, EID_MomentumBalance, defScale) );

    // draw t1
    for ( i = 1; i <= 3; i++ ) {
        p [ 1 ].x = p [ 0 ].x + coeff *tc.at(1, i);
        p [ 1 ].y = p [ 0 ].y + coeff *tc.at(2, i);
        p [ 1 ].z = p [ 0 ].z + coeff *tc.at(3, i);

        EASValsSetColor( ColorGetPixelFromString(oofem_tmpstr(colors [ i - 1 ]), & succ) );

        go = CreateLine3D(p);
        EGWithMaskChangeAttributes(WIDTH_MASK | COLOR_MASK | LAYER_MASK, go);
        EMAddGraphicsToModel(ESIModel(), go);
    }
}
#endif
} // end namespace oofem
