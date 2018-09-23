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

#include "sm/Elements/Interfaces/structuralinterfaceelementphf.h"
#include "sm/Materials/InterfaceMaterials/structuralinterfacematerialphf.h"
#include "sm/Materials/InterfaceMaterials/structuralinterfacematerialstatus.h"

#include "sm/CrossSections/structuralinterfacecrosssection.h"
#include "feinterpol.h"
#include "domain.h"
#include "material.h"
#include "gausspoint.h"
#include "gaussintegrationrule.h"
#include "intarray.h"
#include "floatarray.h"
#include "floatmatrix.h"
#include "mathfem.h"
#include "dof.h"


namespace oofem {
StructuralInterfaceElementPhF :: StructuralInterfaceElementPhF(int n, Domain *aDomain) : StructuralInterfaceElement(n, aDomain)
{
    this->internalLength = 1.0e-1;
}


// remove later
void
StructuralInterfaceElementPhF :: computeDisplacementUnknowns(FloatArray &answer, ValueModeType valueMode, TimeStep *tStep)
{
    IntArray dofIdArray;
    this->giveDofManDofIDMask_u(dofIdArray);
    this->computeVectorOf(dofIdArray, valueMode, tStep, answer);
}


void
StructuralInterfaceElementPhF :: computeDamageUnknowns(FloatArray &answer, ValueModeType valueMode, TimeStep *tStep)
{
    IntArray dofIdArray;
    this->giveDofManDofIDMask_d(dofIdArray);
    this->computeVectorOf(dofIdArray, valueMode, tStep, answer);
}


void
StructuralInterfaceElementPhF :: computeLocationArrayOfDofIDs( const IntArray &dofIdArray, IntArray &answer )
{
    // Routine to compute the local ordering array an element given a dofid array.
    answer.resize( 0 );
    int k = 0;
    for(int i = 1; i <= this->giveNumberOfDofManagers(); i++) {
        DofManager *dMan = this->giveDofManager( i );
        for( int j = 1; j <= dofIdArray.giveSize( ); j++ ) {

            if ( dMan->hasDofID( (DofIDItem) dofIdArray.at( j ) ) ) {
                // hack
                answer.followedBy( k + j );
            }
        }
        k += dMan->giveNumberOfDofs( );
    }
}


void
StructuralInterfaceElementPhF :: giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord)
{
    // Computes the internal forces corresponding to the two fields u & d
    IntArray IdMask_u, IdMask_d;
    this->giveDofManDofIDMask_u( IdMask_u );
    this->giveDofManDofIDMask_d( IdMask_d );
    this->getLocationArray_u(loc_u);
    this->getLocationArray_d(loc_d);
    int ndofs = this->computeNumberOfDofs();
    answer.resize( ndofs);
    answer.zero();
    FloatArray answer_u(0);
    FloatArray answer_d(0);
    this->giveInternalForcesVector_u(answer_u, tStep, useUpdatedGpRecord);
    this->giveInternalForcesVector_d(answer_d, tStep, useUpdatedGpRecord);
    answer.assemble(answer_u, loc_u);
    answer.assemble(answer_d, loc_d);
}


void
StructuralInterfaceElementPhF :: giveInternalForcesVector_u(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord)
{
    FloatMatrix N;
    FloatArray u, traction, jump;

    IntArray dofIdArray;
    this->giveDofManDofIDMask_u(dofIdArray);
    this->computeVectorOf(dofIdArray, VM_Total, tStep, u);

    // subtract initial displacements, if defined
    if ( initialDisplacements.giveSize() ) {
        u.subtract(initialDisplacements);
    }

    // zero answer will resize accordingly when adding first contribution
    answer.clear();

    for ( auto &ip: *this->giveDefaultIntegrationRulePtr() ) {
        this->computeNmatrixAt(ip, N);

        jump.beProductOf(N, u);
        this->computeTraction(traction, ip, jump, tStep);
        //traction.resize(2);

        // compute internal cohesive forces as f = N^T*traction dA
        double dA = this->computeAreaAround(ip);
        answer.plusProduct(N, traction, dA);
    }
}

void
StructuralInterfaceElementPhF :: giveInternalForcesVector_d(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord)
{
    // computes int_A ( N^t * ( d + l*f ) + B^t * l^2 * gradd(d)  )* dA
    FloatArray BStress, Nd, BS, a_d, grad_d;
    FloatMatrix B;
    answer.clear();
    this->computeDamageUnknowns( a_d, VM_Total, tStep );

    double l = this->giveInternalLength();

    for ( auto &ip: *this->giveDefaultIntegrationRulePtr() ) {       

        // Part 1
        StructuralInterfaceMaterialPhF *mat = static_cast< StructuralInterfaceMaterialPhF * >( this->giveInterfaceCrossSection()->giveInterfaceMaterial() );

        double d = computeDamageAt( ip, VM_Total, tStep );
        double facN = d + l * mat->giveDrivingForce(ip);

        this->giveInterpolation()->evalN(Nd, ip->giveNaturalCoordinates(), FEIElementGeometryWrapper(this));       

        double dA = this->computeAreaAround(ip);
        answer.add(facN*dA, Nd);

        // Part 2
        this->computeBd_matrixAt(ip, B );
        grad_d.beProductOf(B, { a_d.at(1), a_d.at(2) });    
        BStress = grad_d * l * l;
        BS.beTProductOf(B, BStress);
        answer.add(dA, BS);
    }

    auto temp = answer;
    answer.resize(4);
    answer.zero();
    answer.assemble(temp,{1,2});
    answer.assemble(temp,{3,4});
}


double
StructuralInterfaceElementPhF :: computeDamageAt(GaussPoint *gp, ValueModeType valueMode, TimeStep *stepN)
{
    // d = N_d * a_d
    FloatArray dVec;
    computeDamageUnknowns(dVec, valueMode, stepN);
    FloatArray Nvec;
    this->giveInterpolation()->evalN(Nvec, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this));
    //return Nvec.dotProduct(dVec);
    return Nvec.dotProduct( {dVec.at(1), dVec.at(2) });
}


void
StructuralInterfaceElementPhF :: computeNd_matrixAt(const FloatArray &lCoords, FloatMatrix &N)
{
    FloatArray Nvec;
    this->giveInterpolation( )->evalN( Nvec, lCoords, FEIElementGeometryWrapper( this ) );
    //N.resize(1, Nvec.giveSize());
    N.beNMatrixOf(Nvec,1);
}

void
StructuralInterfaceElementPhF :: computeBd_matrixAt(GaussPoint *gp, FloatMatrix &answer)
{
    FloatMatrix dNdxi;
    FloatMatrix G;
    this->computeCovarBaseVectorsAt(gp, G);

    this->giveInterpolation( )->evaldNdxi( dNdxi, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper( this ) );
    answer.beProductTOf(G,dNdxi);
    //answer.beTranspositionOf( dNdx );
}


int
StructuralInterfaceElementPhF :: computeNumberOfDofs()
{
    //NLStructuralElement *el = static_cast< NLStructuralElement* >(this->giveElement( ) );
    int nDofs = 0;
    for( int i = 1; i <= this->giveNumberOfDofManagers(); i++) {
        nDofs += this->giveDofManager( i )->giveNumberOfDofs();
    }
    return nDofs;
}


void
StructuralInterfaceElementPhF :: computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
{
    //set displacement and nonlocal location array
    ///@todo this part is enough to do once
    //IntArray IdMask_u, IdMask_d;
    //this->giveDofManDofIDMask_u( IdMask_u );
    //this->giveDofManDofIDMask_d( IdMask_d );

    this->getLocationArray_u(loc_u);
    this->getLocationArray_d(loc_d);

    int nDofs = this->computeNumberOfDofs();
    answer.resize( nDofs, nDofs );
    answer.zero();

    FloatMatrix answer1, answer2, answer3, answer4;

    this->computeStiffnessMatrix_uu(answer1, rMode, tStep);
    //StructuralInterfaceElement :: computeStiffnessMatrix(answer1, rMode, tStep);

    //this->computeStiffnessMatrix_ud(answer2, rMode, tStep);
    //this->computeStiffnessMatrix_du(answer3, rMode, tStep); //symmetric
    //answer3.beTranspositionOf( answer2 );

    this->computeStiffnessMatrix_dd(answer4, rMode, tStep);

    answer.assemble( answer1, loc_u, loc_u );
    //answer.assemble( answer2, loc_u, loc_d );
    //answer.assemble( answer3, loc_d, loc_u );
    answer.assemble( answer4, loc_d, loc_d );
}


void
StructuralInterfaceElementPhF :: computeStiffnessMatrix_uu(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
{
    // This is the regular stiffness matrix times G

    // Computes the stiffness matrix of the receiver K_cohesive = int_A (  N^t * D * N ) dA
    // with the stiffness D = dT/dj
    FloatMatrix N, D, DN;
    bool matStiffSymmFlag = this->giveCrossSection()->isCharacteristicMtrxSymmetric(rMode);
    answer.clear();

    FloatMatrix rotationMatGtoL;
    for ( auto &ip: *this->giveDefaultIntegrationRulePtr() ) {

        this->giveStiffnessMatrix_Eng(D, rMode, ip, tStep);

        this->computeTransformationMatrixAt(ip, rotationMatGtoL);
        D.rotatedWith(rotationMatGtoL, 't');                      // transform stiffness to global coord system

        this->computeNmatrixAt(ip, N);
        DN.beProductOf(D, N);
        double dA = this->computeAreaAround(ip);
        if ( matStiffSymmFlag ) {
            answer.plusProductSymmUpper(N, DN, dA);
        } else {
            answer.plusProductUnsym(N, DN, dA);
        }
    }

    if ( matStiffSymmFlag ) {
        answer.symmetrized();
    }
}


void
StructuralInterfaceElementPhF :: computeStiffnessMatrix_ud(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
{
    FloatMatrix B, DB, N_d, DN, S(3,1);
    FloatArray stress;
#if 1
    answer.clear();

    for ( auto &ip: *this->giveDefaultIntegrationRulePtr() ) {
        double dA = this->computeAreaAround(ip);

        // compute int_V ( B^t * D_B * B )dV
        this->computeNd_matrixAt(ip->giveNaturalCoordinates(), N_d);
        // stress

        FloatMatrix B_u, N;
        this->computeNmatrixAt(ip, N);
        //stress = matStat->giveTempStressVector();
        //stress.times( this->computeGPrim( ip, VM_Total, tStep ) );
        S.setColumn(stress,1);
        DN.beProductOf(S, N_d);
        answer.plusProductUnsym(N, DN, dA);
    }
#endif
}


void
StructuralInterfaceElementPhF :: computeStiffnessMatrix_dd(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
{
    double l = this->giveInternalLength();

    FloatMatrix B_d, N_d;
    answer.clear();

    FloatMatrix rotationMatGtoL;
    for ( auto &ip: *this->giveDefaultIntegrationRulePtr() ) {
        StructuralInterfaceMaterialPhF *mat = static_cast< StructuralInterfaceMaterialPhF * >( this->giveInterfaceCrossSection()->giveInterfaceMaterial() );

        this->computeNd_matrixAt(ip->giveNaturalCoordinates(), N_d);
        this->computeBd_matrixAt(ip, B_d);

        // K_dd1 = ( 1 + l*f' ) * N^t*N;
        // K_dd2 = ( l*l ) * B^t*B;
        double dA = this->computeAreaAround(ip);

        double fPrime = mat->giveDrivingForcePrime(ip);
        double factorN = 1.0 + l* fPrime;
        double factorB = l*l;

        answer.plusProductUnsym(N_d, N_d, factorN * dA);
        answer.plusProductUnsym(B_d, B_d, factorB * dA);
    }

    auto temp = answer;
    answer.resize(4,4);
    answer.zero();
    answer.assemble(temp,{1,2},{1,2});
    answer.assemble(temp,{3,4},{3,4});
}


void
StructuralInterfaceElementPhF :: computeTraction(FloatArray &traction, IntegrationPoint *ip, FloatArray &jump, TimeStep *tStep)
{
    // Returns the traction in global coordinate system
    FloatMatrix rotationMatGtoL, F;
    this->computeTransformationMatrixAt(ip, rotationMatGtoL);
    jump.rotatedWith(rotationMatGtoL, 'n');      // transform jump to local coord system

    double damage = this->computeDamageAt(ip, VM_Total, tStep);

    this->giveEngTraction(traction, ip, jump, damage, tStep);

    traction.rotatedWith(rotationMatGtoL, 't');     // transform traction to global coord system
}


void
StructuralInterfaceElementPhF :: updateYourself(TimeStep *tStep)
{
    Element :: updateYourself(tStep);

    // record initial displacement if element not active
    if ( activityTimeFunction && !isActivated(tStep) ) {
        this->computeVectorOf(VM_Total, tStep, initialDisplacements);
    }
}


void
StructuralInterfaceElementPhF :: updateInternalState(TimeStep *tStep)
{
    // Updates the receiver at end of step.
    FloatArray tractionG, jumpL;

    // force updating strains & stresses
    for ( auto &iRule: integrationRulesArray ) {
        for ( auto &gp: *iRule ) {
            this->computeSpatialJump(jumpL, gp, tStep);
            this->computeTraction(tractionG, gp, jumpL, tStep);
        }
    }
}


} // end namespace oofem

