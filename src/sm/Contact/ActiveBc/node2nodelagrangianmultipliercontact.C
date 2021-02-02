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


#include "sm/Contact/ActiveBc/node2nodelagrangianmultipliercontact.h"
#include "set.h"
#include "domain.h"
#include "node.h"
#include "masterdof.h"
#include "floatmatrix.h"
#include "unknownnumberingscheme.h"
#include "sparsemtrx.h"
#include "classfactory.h"

#ifdef _OPENMP
#include <omp.h>
#endif

namespace oofem {
REGISTER_BoundaryCondition(Node2NodeLagrangianMultiplierContact);


Node2NodeLagrangianMultiplierContact :: Node2NodeLagrangianMultiplierContact(int n, Domain *d) : ActiveBoundaryCondition(n, d), lmdm()
{}


void
Node2NodeLagrangianMultiplierContact :: initializeFrom(InputRecord &ir)
{
    ActiveBoundaryCondition :: initializeFrom(ir);
    this->useTangent = ir.hasField(_IFT_Node2NodeLagrangianMultiplierContact_useTangent);
    IR_GIVE_FIELD(ir, this->masterSet, _IFT_Node2NodeLagrangianMultiplierContact_masterSet);
    IR_GIVE_FIELD(ir, this->slaveSet, _IFT_Node2NodeLagrangianMultiplierContact_slaveSet);

    // these are internal lagrange multipliers used to enforce the no penetration for the contacting bodies
    for ( int pos = 0; pos < masterSet.giveSize(); ++pos ) {
        lmdm.push_back(std :: unique_ptr< Node >(new Node(0, domain) ) );
        lmdm.at(pos)->appendDof(new MasterDof(this->lmdm.at(pos).get(), ( DofIDItem ) ( this->giveDomain()->giveNextFreeDofID() ) ) );
    }
}


void
Node2NodeLagrangianMultiplierContact :: assemble(SparseMtrx &answer, TimeStep *tStep,
                                                 CharType type, const UnknownNumberingScheme &r_s, const UnknownNumberingScheme &c_s, 
                                                 double scale,
                                                 void* lock)
{
    if ( !this->useTangent || type != TangentStiffnessMatrix ) {
        return;
    }

    FloatMatrix K;
    IntArray loc, c_loc;

    IntArray dofIdArray = {
        D_u, D_v
    };

    if ( masterSet.giveSize() != slaveSet.giveSize() ) {
        OOFEM_ERROR("Number of master nodes does not match number of slave nodes")
    }

    std :: vector< IntArray >lambdaeq;
    this->giveLagrangianMultiplierLocationArray(r_s, lambdaeq);
    for ( int pos = 1; pos <= masterSet.giveSize(); ++pos ) {
        Node *masterNode = this->giveDomain()->giveNode(masterSet.at(pos) );
        Node *slaveNode = this->giveDomain()->giveNode(slaveSet.at(pos) );

        masterNode->giveLocationArray(dofIdArray, loc, r_s);
        slaveNode->giveLocationArray(dofIdArray, c_loc, c_s);
        loc.followedBy(c_loc);

        double gap = this->computeTangentFromContact(K, masterNode, slaveNode, tStep);
        if ( gap >= 0 ) {
            // to make the equation system regular in the case there is no contact we initialize the allocated equation to the following form 1*labmda = 0, forcing lagrange multiplier of inactive condition to be zero.
            FloatArray one(1);
            one.at(1) = 1;
#ifdef _OPENMP
            if (lock) omp_set_lock(static_cast<omp_lock_t*>(lock));
#endif
            answer.assemble(lambdaeq.at(pos - 1), one);
#ifdef _OPENMP
            if (lock) omp_unset_lock(static_cast<omp_lock_t*>(lock));
#endif
        } else {
#ifdef _OPENMP
            if (lock) omp_set_lock(static_cast<omp_lock_t*>(lock));
#endif
            answer.assemble(loc, lambdaeq.at(pos - 1), K);
#ifdef _OPENMP
            if (lock) omp_unset_lock(static_cast<omp_lock_t*>(lock));
#endif
            FloatMatrix Kt;
            Kt.beTranspositionOf(K);
#ifdef _OPENMP
            if (lock) omp_set_lock(static_cast<omp_lock_t*>(lock));
#endif
            answer.assemble(lambdaeq.at(pos - 1), loc, Kt);
#ifdef _OPENMP
            if (lock) omp_unset_lock(static_cast<omp_lock_t*>(lock));
#endif
        }
    }
}

void
Node2NodeLagrangianMultiplierContact :: assembleVector(FloatArray &answer, TimeStep *tStep,
                                                       CharType type, ValueModeType mode,
                                                       const UnknownNumberingScheme &s, 
                                                       FloatArray *eNorms,
                                                       void* lock)
{
    if ( masterSet.giveSize() != slaveSet.giveSize() ) {
        OOFEM_ERROR("Number of master nodes does not match number of slave nodes");
    }


    //IntArray dofIdArray = {D_u, D_v, D_w};
    IntArray dofIdArray = {
        D_u, D_v
    };


    if ( type == InternalForcesVector ) {
        // assemble lagrangian multiplier contribution to residuum
        int size = masterSet.giveSize();
        // assemble location array
        for ( int pos = 1; pos <= size; pos++ ) {
            IntArray loc, s_loc;
            FloatArray n;
            Node *masterNode = this->giveDomain()->giveNode(masterSet.at(pos) );
            Node *slaveNode = this->giveDomain()->giveNode(slaveSet.at(pos) );
            this->computeNormalMatrixAt(n, masterNode, slaveNode, tStep);
            Dof *mdof = * ( lmdm.at(pos - 1)->begin() );
            n.times(mdof->giveUnknown(mode, tStep) );
            masterNode->giveLocationArray(dofIdArray, loc, s);
            slaveNode->giveLocationArray(dofIdArray, s_loc, s);
            loc.followedBy(s_loc);
#ifdef _OPENMP
            if (lock) omp_set_lock(static_cast<omp_lock_t*>(lock));
#endif
            answer.assemble(n, loc);
#ifdef _OPENMP
            if (lock) omp_unset_lock(static_cast<omp_lock_t*>(lock));
#endif
        }
    } else if ( type == ExternalForcesVector ) {
        IntArray loc, c_loc;
        FloatArray fext;


        std :: vector< IntArray >lambdaeq;
        this->giveLagrangianMultiplierLocationArray(s, lambdaeq);

        for ( int pos = 1; pos <= masterSet.giveSize(); ++pos ) {
            Node *masterNode = this->giveDomain()->giveNode(masterSet.at(pos) );
            Node *slaveNode = this->giveDomain()->giveNode(slaveSet.at(pos) );
            this->computeExternalForcesFromContact(fext, masterNode, slaveNode, tStep);
#ifdef _OPENMP
            if (lock) omp_set_lock(static_cast<omp_lock_t*>(lock));
#endif
            answer.assemble(fext, lambdaeq.at(pos - 1) );
#ifdef _OPENMP
            if (lock) omp_unset_lock(static_cast<omp_lock_t*>(lock));
#endif
        }
    }
}



double
Node2NodeLagrangianMultiplierContact :: computeTangentFromContact(FloatMatrix &answer, Node *masterNode, Node *slaveNode, TimeStep *tStep)
{
    double gap;
    FloatArray Nv;
    this->computeGap(gap, masterNode, slaveNode, tStep);
    this->computeNormalMatrixAt(Nv, masterNode, slaveNode, tStep);
    answer.initFromVector(Nv, false);

    return gap;
    //answer.times(this->penalty);
    /*    if ( gap > 0.0 ) {
     *  int size = this->masterSet.giveSize();
     *  FloatMatrix help(2 * size, 2 * size);
     *  answer.beUnitMatrix();
     *  }*/
}

void
Node2NodeLagrangianMultiplierContact :: computeGap(double &answer, Node *masterNode, Node *slaveNode, TimeStep *tStep)
{
    FloatArray uS, uM;
    auto xs = slaveNode->giveCoordinates();
    auto xm = masterNode->giveCoordinates();
    FloatArray normal = xs - xm;
    double norm = normal.computeNorm();
    if ( norm < 1.0e-8 ) {
        OOFEM_ERROR("Couldn't compute normal between master node (num %d) and slave node (num %d), nodes are too close to each other.", masterNode->giveGlobalNumber(), slaveNode->giveGlobalNumber() );
    } else {
        normal.times(1.0 / norm);
    }

    slaveNode->giveUnknownVector(uS, { D_u, D_v }, VM_Total, tStep, true);
    masterNode->giveUnknownVector(uM, { D_u, D_v }, VM_Total, tStep, true);
    xs.add(uS);
    xm.add(uM);
    FloatArray dx = xs - xm;
    answer = dx.dotProduct(normal);
}


void
Node2NodeLagrangianMultiplierContact :: computeNormalMatrixAt(FloatArray &answer, Node *masterNode, Node *slaveNode, TimeStep *TimeStep)
{
    const auto &xs = slaveNode->giveCoordinates();
    const auto &xm = masterNode->giveCoordinates();
    auto normal = xs - xm;
    double norm = normal.computeNorm();
    if ( norm < 1.0e-8 ) {
        OOFEM_ERROR("Couldn't compute normal between master node (num %d) and slave node (num %d), nodes are too close to each other.", masterNode->giveGlobalNumber(), slaveNode->giveGlobalNumber() );
    } else {
        normal.times(1.0 / norm);
    }
    // The normal is not updated for node2node which is for small deformations only
    // C = {n -n}
    answer = {
        normal.at(1), normal.at(2),
        -normal.at(1), -normal.at(2)
    };
}


void
Node2NodeLagrangianMultiplierContact :: computeExternalForcesFromContact(FloatArray &answer, Node *masterNode, Node *slaveNode, TimeStep *tStep)
{
    answer.resize(1);
    this->computeGap(answer.at(1), masterNode, slaveNode, tStep);
    if ( answer.at(1) > 0.0 ) {
        answer.at(1) = 0.0;
    }
}


void
Node2NodeLagrangianMultiplierContact :: giveLagrangianMultiplierLocationArray(const UnknownNumberingScheme &r_s, std :: vector< IntArray > &answer)
{
    int size = this->masterSet.giveSize();
    answer.resize(size);
    IntArray dofIdArray = {
        D_u, D_v
    };

    // assemble location array
    IntArray l(1);
    for ( int i = 0; i < size; i++ ) {
        l.at(1) = r_s.giveDofEquationNumber(* lmdm.at(i)->begin() );
        answer.at(i) = l;
    }
}


void
Node2NodeLagrangianMultiplierContact :: giveLocationArrays(std :: vector< IntArray > &rows, std :: vector< IntArray > &cols, CharType type, const UnknownNumberingScheme &r_s, const UnknownNumberingScheme &c_s)
{
    IntArray r_loc, c_loc;
    rows.resize(3 * masterSet.giveSize() );
    cols.resize(3 * masterSet.giveSize() );
    IntArray dofIdArray = {
        D_u, D_v
    };
    std :: vector< IntArray >lambdaeq;
    this->giveLagrangianMultiplierLocationArray(r_s, lambdaeq);

    for ( int pos = 1; pos <= masterSet.giveSize(); pos++ ) {
        Node *masterNode = this->giveDomain()->giveNode(masterSet.at(pos) );
        Node *slaveNode = this->giveDomain()->giveNode(slaveSet.at(pos) );

        masterNode->giveLocationArray(dofIdArray, r_loc, r_s);
        slaveNode->giveLocationArray(dofIdArray, c_loc, c_s);

        // column block
        rows [ 0 + 3 * ( pos - 1 ) ] = r_loc;
        cols [ 0 + 3 * ( pos - 1 ) ] = lambdaeq.at(pos - 1);
        // row block
        cols [ 1 + 3 * ( pos - 1 ) ] = c_loc;
        rows [ 1 + 3 * ( pos - 1 ) ] = lambdaeq.at(pos - 1);
        // diagonal enry (some sparse mtrx implementation requaire this)
        rows [ 2 + 3 * ( pos - 1 ) ] = lambdaeq.at(pos - 1);
        cols [ 2 + 3 * ( pos - 1 ) ] = lambdaeq.at(pos - 1);
    }
}
} // namespace oofem
