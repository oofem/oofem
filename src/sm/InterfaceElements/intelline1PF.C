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

#include "sm/InterfaceElements/intelline1PF.h"
#include "sm/Elements/Interfaces/intelline1.h"
#include "node.h"
#include "sm/CrossSections/structuralinterfacecrosssection.h"
#include "sm/Materials/InterfaceMaterials/structuralinterfacematerialstatus.h"
#include "gausspoint.h"
#include "gaussintegrationrule.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "intarray.h"
#include "mathfem.h"
#include "fei2dlinelin.h"
#include "classfactory.h"


namespace oofem {
REGISTER_Element(IntElLine1PF);

FEI2dLineLin IntElLine1PF :: interp(1, 1);


IntElLine1PF :: IntElLine1PF(int n, Domain *aDomain) :
    StructuralInterfaceElement(n, aDomain), PhaseFieldElement(n, aDomain)
{
    numberOfDofMans = 4;
    alpha.resize(2);
    alpha.zero();
    this->oldAlpha.resize(2);
    this->oldAlpha.zero();
    this->deltaAlpha.resize(2);
    this->deltaAlpha.zero();
    this->prescribed_damage = -1.0;
}



void
IntElLine1PF :: computeNmatrixAt(GaussPoint *ip, FloatMatrix &answer)
{
    // Returns the modified N-matrix which multiplied with u give the spatial jump.

    FloatArray N;
    FEInterpolation *interp = this->giveInterpolation();
    interp->evalN( N, ip->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );

    answer.resize(2, 8);
    answer.zero();
    answer.at(1, 1) = answer.at(2, 2) = -N.at(1);
    answer.at(1, 3) = answer.at(2, 4) = -N.at(2);

    answer.at(1, 5) = answer.at(2, 6) = N.at(1);
    answer.at(1, 7) = answer.at(2, 8) = N.at(2);
}


void
IntElLine1PF :: computeGaussPoints()
{
    if ( integrationRulesArray.size() == 0 ) {
        integrationRulesArray.resize( 1 );
        //integrationRulesArray[ 0 ] = std::make_unique<LobattoIntegrationRule>(1,domain, 1, 2);
        integrationRulesArray [ 0 ] = std::make_unique<GaussIntegrationRule>(1, this, 1, 2);
        integrationRulesArray [ 0 ]->SetUpPointsOnLine(this->numberOfGaussPoints, _2dInterface);
    }
}


FloatArrayF<2>
IntElLine1PF :: computeCovarBaseVectorAt(IntegrationPoint *ip) const
{
    FEInterpolation *interp = this->giveInterpolation();
    FloatMatrix dNdxi;
    interp->evaldNdxi( dNdxi, ip->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );
    FloatArrayF<2> G;
    int numNodes = this->giveNumberOfNodes();
    for ( int i = 1; i <= dNdxi.giveNumberOfRows(); i++ ) {
        double X1_i = 0.5 * ( this->giveNode(i)->giveCoordinate(1) + this->giveNode(i + numNodes / 2)->giveCoordinate(1) ); // (mean) point on the fictious mid surface
        double X2_i = 0.5 * ( this->giveNode(i)->giveCoordinate(2) + this->giveNode(i + numNodes / 2)->giveCoordinate(2) );
        G.at(1) += dNdxi.at(i, 1) * X1_i;
        G.at(2) += dNdxi.at(i, 1) * X2_i;
    }
    return G;
}


double
IntElLine1PF :: computeAreaAround(IntegrationPoint *ip)
{
    auto G = this->computeCovarBaseVectorAt(ip);
    auto ds = norm(G) * ip->giveWeight();
    return ds * this->giveCrossSection()->give(CS_Thickness, ip);
}


void
IntElLine1PF :: initializeFrom(InputRecord &ir)
{
    StructuralInterfaceElement :: initializeFrom(ir);

    if ( ir.hasField(_IFT_IntElLine1PF_prescribedDamage) ) {
        IR_GIVE_FIELD(ir, this->prescribed_damage, _IFT_IntElLine1PF_prescribedDamage);
        //this->alpha.setValues(2, damage, damage);
    }
}


void
IntElLine1PF :: giveDofManDofIDMask(int inode, IntArray &answer) const
{
    answer = {D_u, D_v, T_f};
}


void
IntElLine1PF :: computeTransformationMatrixAt(GaussPoint *gp, FloatMatrix &answer)
{
    // Transformation matrix to the local coordinate system
    FloatArray G;
    this->computeCovarBaseVectorAt(gp, G);
    G.normalize();

    answer.resize(2, 2);
    answer.at(1, 1) =  G.at(1);
    answer.at(2, 1) =  G.at(2);
    answer.at(1, 2) = -G.at(2);
    answer.at(2, 2) =  G.at(1);

}


FEInterpolation *
IntElLine1PF :: giveInterpolation() const
{
    return & interp;
}


void
IntElLine1PF :: computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
{
    //set displacement and nonlocal location array
    ///@todo this part is enough to do once
    IntArray IdMask_u, IdMask_d;
    this->giveDofManDofIDMask_u( IdMask_u );
    this->giveDofManDofIDMask_d( IdMask_d );
    this->computeLocationArrayOfDofIDs( IdMask_u, loc_u );
    this->computeLocationArrayOfDofIDs( IdMask_d, loc_d );

    int nDofs = this->computeNumberOfDofs();
    answer.resize( nDofs, nDofs );
    answer.zero();

    // Set the solution vectors (a_u, a_d) for the element
    this->computeDisplacementUnknowns( this->unknownVectorU, VM_Total, tStep );
    this->computeDamageUnknowns( this->unknownVectorD, VM_Total, tStep );
    this->computeDamageUnknowns( this->deltaUnknownVectorD, VM_Incremental, tStep );

    FloatMatrix temp;
    solveForLocalDamage(temp, tStep);

    FloatMatrix answer1, answer2, answer3, answer4;
    this->computeStiffnessMatrix_uu(answer1, rMode, tStep);
    //this->computeStiffnessMatrix_ud(answer2, rMode, tStep);
    //this->computeStiffnessMatrix_du(answer3, rMode, tStep); //symmetric
    answer3.beTranspositionOf( answer2 );
    this->computeStiffnessMatrix_dd(answer4, rMode, tStep);

    answer.assemble( answer1, loc_u, loc_u );
    //answer.assemble( answer2, loc_u, loc_d );
    //answer.assemble( answer3, loc_d, loc_u );
    answer.assemble( answer4, loc_d, loc_d );
}


void
IntElLine1PF :: computeStiffnessMatrix_uu(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
{
    // Computes the stiffness matrix of the receiver K_cohesive = int_A ( N^t * dT/dj * N ) dA
    FloatMatrix N, D, DN, GD, Gmat;
    bool matStiffSymmFlag = this->giveCrossSection()->isCharacteristicMtrxSymmetric(rMode);
    answer.resize(0, 0);
    FloatArray Nd;
    //IntegrationRule *iRule = integrationRulesArray [ giveDefaultIntegrationRule() ];
    FloatMatrix rotationMatGtoL;

    //for ( int j = 0; j < iRule->giveNumberOfIntegrationPoints(); j++ ) {
        //IntegrationPoint *ip = iRule->getIntegrationPoint(j);
    for ( auto &ip: *this->giveDefaultIntegrationRulePtr() ) {    

        if ( this->nlGeometry == 0 ) {
            this->giveStiffnessMatrix_Eng(D, rMode, ip, tStep);
        } else if ( this->nlGeometry == 1 ) {
            this->giveStiffnessMatrix_dTdj(D, rMode, ip, tStep);
        } else {
            OOFEM_ERROR("nlgeometry must be 0 or 1!")
        }

        this->computeNmatrixAt(ip, N);

        this->computeNd_vectorAt(ip->giveNaturalCoordinates(), Nd);
        double d = Nd.dotProduct(this->alpha);
        this->computeGMatrix(Gmat, d, ip, VM_Total, tStep);
        GD.beProductOf(Gmat,D);

        this->computeTransformationMatrixAt(ip, rotationMatGtoL);
        GD.rotatedWith(rotationMatGtoL, 't');                      // transform stiffness to global coord system

        DN.beProductOf(GD, N);

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
IntElLine1PF :: computeStiffnessMatrix_ud(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
{
    //IntegrationRule *iRule = integrationRulesArray [ giveDefaultIntegrationRule() ];

    FloatMatrix N, rotationMatGtoL;
    FloatArray a_u, traction, tractionTemp, jump, fu, fd(2), fd4(4);

    fd.zero();
    FloatArray a_d_temp, a_d, Bd, Nd;
    a_u = this->unknownVectorU;
    a_d_temp = this->unknownVectorD;
    a_d = { a_d_temp.at(1), a_d_temp.at(2) }; 

    FloatMatrix temp, Kud(8,2);
    fu.zero();
    Kud.zero();
    //for ( int i = 0; i < iRule->giveNumberOfIntegrationPoints(); i++ ) {
    for ( auto &ip: *this->giveDefaultIntegrationRulePtr() ) {      
        //IntegrationPoint *ip = iRule->getIntegrationPoint(i);
        this->computeNmatrixAt(ip, N);
        this->computeNd_vectorAt(ip->giveNaturalCoordinates(), Nd);
        jump.beProductOf(N, a_u);
        this->computeTraction(traction, ip, jump, tStep);

        // compute internal cohesive forces as f = N^T*traction dA
        double Gprim = computeGPrim(ip, VM_Total, tStep);
        double dA = this->computeAreaAround(ip);

        fu.plusProduct(N, traction, dA*Gprim);
        temp.beDyadicProductOf(fu,Nd);
        Kud.add(temp);

    }

    IntArray indxu, indx1, indx2;
    indxu = {1, 2, 3, 4, 5, 6, 7, 8};
    indx1 = {1, 2};
    indx2 = {3, 4};
    answer.resize(8,4);
    answer.zero();
    answer.assemble(Kud, indxu, indx1);    
    answer.assemble(Kud, indxu, indx2);


    // Numerical tangent
    // Set the solution vectors (a_u, a_d) for the element
    FloatArray a_u_ref, a_d_pert, a_d_ref, fu_ref, fd_ref, fu_pert, K_col, fd_temp, delta_a_d_ref, delta_a_d_pert;
    FloatMatrix K_num(8,2), Kdd_num;
    double eps = 1.0e-6;

    a_u_ref = this->unknownVectorU;
    a_d_ref = this->unknownVectorD;
    delta_a_d_ref = this->deltaUnknownVectorD;
    this->giveInternalForcesVectorUD(fu_ref, fd_ref, tStep, false);

    for ( int i = 1; i <=2; i++ ) {
        a_d_pert = a_d_ref;
        a_d_pert.at(i) += eps;
        this->unknownVectorD = a_d_pert;

        delta_a_d_pert = delta_a_d_ref;
        delta_a_d_pert.at(i) += eps;
        this->deltaUnknownVectorD = delta_a_d_pert;

        this->giveInternalForcesVectorUD(fu_pert, fd_temp, tStep, false);
        K_col = ( 1.0/eps ) * ( fu_pert - fu_ref );
        K_num.addSubVectorCol(K_col, 1, i);
    }

    answer.zero();
    answer.assemble(K_num, indxu, indx1);
    answer.assemble(K_num, indxu, indx2);

    this->unknownVectorU = a_u_ref;
    this->unknownVectorD = a_d_ref;
    this->deltaUnknownVectorD = delta_a_d_ref;
}


void
IntElLine1PF :: computeStiffnessMatrix_dd(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep) 
{
    // Computation of tangent: K_dd = \int Nd^t * ( -kp*neg_Mac'(alpha_dot)/delta_t + g_c/l + G''*Psibar) * Nd + 
    //                                \int Bd^t * (  g_c * l * [G^1 G^2]^t * [G^1 G^2] ) * Bd 
    //                              = K_dd1 + K_dd2
    int ndofs   = 8 + 4;
    int ndofs_u = 8;
    int ndofs_d = 4;

    double l       = this->giveInternalLength();
    double g_c     = this->giveCriticalEnergy();
    double kp      = this->givePenaltyParameter();
    double Delta_t = tStep->giveTimeIncrement();

    answer.resize( ndofs_d, ndofs_d );
    answer.zero();

    //IntegrationRule *iRule = integrationRulesArray [ giveDefaultIntegrationRule() ];
    FloatMatrix tempN(2,2), tempB(2,2), temp(2,2);
    FloatArray Nd, Bd;
    tempN.zero();
    tempB.zero();
    temp.zero();
    //for ( int j = 0; j < iRule->giveNumberOfIntegrationPoints(); j++ ) {
    //    IntegrationPoint *ip = iRule->getIntegrationPoint(j);
    for ( auto &ip: *this->giveDefaultIntegrationRulePtr() ) {
        computeBd_vectorAt(ip, Bd);
        computeNd_vectorAt(ip->giveNaturalCoordinates(), Nd);
        double dA = this->computeAreaAround(ip);

        //double Gbis   = this->computeGbis()
        double Gbis = 2.0;
        double Psibar  = this->computeFreeEnergy( ip, tStep );

        // K_dd1 = ( -kp*neg_Mac'(d_dot) / Delta_t + g_c/ l + G'' * Psibar ) * N^t*N; 
        double Delta_d = computeDamageAt(ip, VM_Incremental, tStep);
        double factorN = -kp * neg_MaCauleyPrime(Delta_d/Delta_t)/Delta_t +  g_c / l + Gbis * Psibar; 
        //double factorN = g_c / l + Gbis * Psibar; 

        //double Psibar0 = this->givePsiBar0();
        //double factorN = g_c / l + Gbis * this->MaCauley(Psibar-Psibar0);

        tempN.beDyadicProductOf(Nd, Nd);
        temp.add(factorN * dA, tempN);

        // K_dd2 =  g_c * l * Bd^t * Bd;
        double factorB = g_c * l;
        tempB.beDyadicProductOf(Bd, Bd);
        temp.add(factorB * dA, tempB);
    }

    IntArray indx1 = {1, 2}, indx2 = {3, 4};

    answer.assemble(temp, indx1, indx1);
    answer.assemble(temp, indx2, indx2);

#if 0
    // Numerical tangent
    // Set the solution vectors (a_u, a_d) for the element
    FloatArray a_u_ref, a_d_ref, a_d_pert, fu_ref, fd_ref, fd_pert, K_col, fu_temp;
    FloatMatrix K(4,4), Kdd_num;
    double eps = 1.0e-6;

    a_u_ref = this->unknownVectorU;
    a_d_ref = this->unknownVectorD;
    this->giveInternalForcesVectorUD(fu_ref, fd_ref, tStep, false);

    for ( int i = 1; i <=2; i++ ) {
        a_d_pert = a_d_ref;
        a_d_pert.at(i) += eps;
        this->unknownVectorD = a_d_pert;

        this->giveInternalForcesVectorUD(fu_temp, fd_pert, tStep, false);
        K_col = ( 1.0/eps ) * ( fd_pert - fd_ref );
        K.addSubVectorCol(K_col, 1, i);
    }

    Kdd_num.beSubMatrixOf(K,1,2,1,2);

    answer.zero();
    answer.assemble(Kdd_num, indx1, indx1);
    answer.assemble(Kdd_num, indx2, indx2);

    this->unknownVectorU = a_u_ref;
    this->unknownVectorD = a_d_ref;
#endif
}


void 
IntElLine1PF :: solveForLocalDamage(FloatMatrix &answer, TimeStep *tStep)
{
    // approach with dame defined over the element only

    if ( tStep->giveNumber() == 1 ) {
        return;
    }
    if ( this->prescribed_damage >-1.0e-8 ){
        return;
    }

    //IntegrationRule *iRule = integrationRulesArray [ giveDefaultIntegrationRule() ];

    FloatArray a_u, a_d_temp, a_d(2), traction, jump, fd(2), fd_ref(2), Nd, Bd;
    FloatMatrix Kdd(2,2), tempN, tempB;
    FloatArray delta_a_d;

    double kp      = this->givePenaltyParameter();
    double Delta_t = tStep->giveTimeIncrement();
    double l       = this->giveInternalLength();
    double g_c     = this->giveCriticalEnergy();

    fd.zero();
    Kdd.zero();
    //a_d.zero();
    a_d = this->alpha; // checked - evolves
    this->oldAlpha = this->alpha;
    for ( int nIter = 1; nIter <= 10; nIter++  ) {
        fd.zero();
        Kdd.zero();
        fd_ref.zero();
        //for ( int i = 0; i < iRule->giveNumberOfIntegrationPoints(); i++ ) {
        //    IntegrationPoint *ip = iRule->getIntegrationPoint(i);
        for ( auto &ip: *this->giveDefaultIntegrationRulePtr() ) {        

            this->computeBd_vectorAt(ip, Bd);
            this->computeNd_vectorAt(ip->giveNaturalCoordinates(), Nd);
            double dA = this->computeAreaAround(ip);

            // Internal force
            //double d       = computeDamageAt( ip, VM_Total, tStep);
            double d = Nd.dotProduct(a_d);
            double gradd   = Bd.dotProduct(a_d); // Dalpha/Ds
            //double Delta_d = computeDamageAt( ip, VM_Incremental, tStep);
            double Delta_d = Nd.dotProduct(a_d-oldAlpha);
            //double Gprim   = computeGPrim(ip, VM_Total, tStep);
            double Gprim   = -2.0 * (1.0 - d);
            double Psibar  = this->computeFreeEnergy( ip, tStep );

            double sectionalForcesScal = -kp * neg_MaCauley(Delta_d/Delta_t) + g_c / l * d + Gprim * Psibar;
	        double sectionalForcesVec  = g_c * l * gradd;
            fd = fd + ( Nd*sectionalForcesScal + Bd*sectionalForcesVec ) * dA;
            fd_ref = fd_ref + ( Nd * (Gprim * Psibar + 1.0e-8 ) ) * dA;

            // Tangent
            //double Gbis   = this->computeGbis()
            double Gbis = 2.0;

            // K_dd1 = ( -kp*neg_Mac'(d_dot) / Delta_t + g_c/ l + G'' * Psibar ) * N^t*N; 
            double factorN = -kp * neg_MaCauleyPrime(Delta_d/Delta_t)/Delta_t +  g_c / l + Gbis * Psibar; 
            //double factorN = g_c / l + Gbis * Psibar; 

            tempN.beDyadicProductOf(Nd, Nd);
            Kdd.add(factorN * dA, tempN);

            // K_dd2 =  g_c * l * Bd^t * Bd;
            double factorB = g_c * l;
            tempB.beDyadicProductOf(Bd, Bd);
            Kdd.add(factorB * dA, tempB);

        }
        //printf("norm %e \n", fd.computeNorm() );
        //if( fd.computeNorm() < 1.0e-5 ) {
        if ( fd.computeNorm()/fd_ref.computeNorm() < 1.0e-3 ) {
            this->alpha = a_d;
            return;
        }
        Kdd.solveForRhs(fd, delta_a_d);
        a_d.subtract(delta_a_d);
        this->deltaAlpha += delta_a_d;
        //this->alpha = a_d;
        //a_d.printYourself();
    }
    printf("norm %e \n", fd.computeNorm()/fd_ref.computeNorm() );
    OOFEM_ERROR("No convergence in phase field iterations")
}


void 
IntElLine1PF :: computeGMatrix(FloatMatrix &answer, const double damage,  GaussPoint *gp, ValueModeType valueMode, TimeStep *tStep) 
{
    double d = damage;
    if ( this->prescribed_damage > -1.0e-8 ) {
        d = prescribed_damage;
    }

    StructuralInterfaceMaterialStatus *matStat = static_cast< StructuralInterfaceMaterialStatus * >( gp->giveMaterialStatus() );
    const auto &strain = matStat->giveTempJump();
    double g2 = -1.0;
    if ( strain.at(3) < 0.0 ) { // Damage does not affect compressive stresses
        g2 = 1.0;
        //printf("compression \n");
    } else {
        g2 = (1.0 - d) * (1.0 - d);
        //printf("g %e \n", g2);
        //g2 = 1;
    }

    double g1 = (1.0 - d) * (1.0 - d);

    answer.resize(2,2);
    answer.zero();
    answer.at(1,1) = g1;
    answer.at(2,2) = g2;
}


double 
IntElLine1PF :: computeDamageAt(GaussPoint *gp, ValueModeType valueMode, TimeStep *stepN)
{
    // d = N_d * a_d
    //NLStructuralElement *el = static_cast< NLStructuralElement* >(this->giveElement( ) );
    StructuralInterfaceElement *el = this->giveElement( );
    FloatArray dVec;
    if ( valueMode == VM_Total ) {
        //dVec = this->unknownVectorD;
        //// use old value
        ////dVec.subtract(this->deltaUnknownVectorD);
        dVec = this->alpha;
    } else if ( valueMode == VM_Incremental ) {
        //dVec = this->deltaUnknownVectorD;
        //dVec = this->deltaAlpha;
        dVec = this->alpha - this->oldAlpha;
    }
    //dVec.resizeWithValues(2);
    FloatArray Nvec;
    el->giveInterpolation()->evalN(Nvec, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(el));
    return Nvec.dotProduct(dVec);
}


void
IntElLine1PF :: giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord)
{
    // Set the solution vectors (a_u, a_d) for the element
    this->computeDisplacementUnknowns( this->unknownVectorU, VM_Total, tStep );
    this->computeDamageUnknowns( this->unknownVectorD, VM_Total, tStep );
    this->computeDamageUnknowns( this->deltaUnknownVectorD, VM_Incremental, tStep );

    FloatArray fu, fd;
    this->giveInternalForcesVectorUD(fu, fd, tStep, useUpdatedGpRecord);

    // total vec
    IntArray IdMask_u, IdMask_d;
    this->giveDofManDofIDMask_u( IdMask_u );
    this->giveDofManDofIDMask_d( IdMask_d );
    this->computeLocationArrayOfDofIDs( IdMask_u, loc_u );
    this->computeLocationArrayOfDofIDs( IdMask_d, loc_d );

    int nDofs = this->computeNumberOfDofs();
    answer.resize( nDofs );
    answer.zero();
    answer.assemble(fu, loc_u);
    //answer.assemble(fd, loc_d);
}


void
IntElLine1PF :: giveInternalForcesVectorUD(FloatArray &fu, FloatArray &fd4, TimeStep *tStep, int useUpdatedGpRecord)
{
    // Computes internal forces
    // if useGpRecnord == 1 then data stored in ip->giveStressVector() are used
    // instead computing stressVector through this->ComputeStressVector();
    // this must be done after you want internal forces after element->updateYourself()
    // has been called for the same time step.

    //IntegrationRule *iRule = integrationRulesArray [ giveDefaultIntegrationRule() ];

    FloatMatrix Nu, rotationMatGtoL, Gmat;
    FloatArray a_u, a_d_temp, a_d, traction, jump, fd(2), Nd, Bd, GTraction;

    a_u      = this->unknownVectorU; 
    a_d_temp = this->unknownVectorD;
    a_d = { a_d_temp.at(1), a_d_temp.at(2) }; 

    fu.zero();
    fd.zero();

    for ( auto &ip: *this->giveDefaultIntegrationRulePtr() ) {    
    //for ( int i = 0; i < iRule->giveNumberOfIntegrationPoints(); i++ ) {
    //    IntegrationPoint *ip = iRule->getIntegrationPoint(i);
        this->computeNmatrixAt(ip, Nu);
        this->computeBd_vectorAt(ip, Bd);
        this->computeNd_vectorAt(ip->giveNaturalCoordinates(), Nd);

        // compute internal cohesive forces as f = N^T * g * traction dA
        jump.beProductOf(Nu, a_u);
        this->computeTraction(traction, ip, jump, tStep);

        //FloatMatrix temp;
        //solveForLocalDamage(temp, tStep);
        //double g  = this->computeG(ip, VM_Total, tStep);

        this->computeNd_vectorAt(ip->giveNaturalCoordinates(), Nd);
        double d = Nd.dotProduct(this->alpha);
        //    g = (1.0 - d) * (1.0 - d);

 
        this->computeGMatrix(Gmat, d, ip, VM_Total, tStep);
        GTraction.beProductOf(Gmat,traction);
        //double g  = this->computeOldG(ip, VM_Total, tStep);
        double dA = this->computeAreaAround(ip);

        //fu.plusProduct(Nu, traction, dA*g);
        fu.plusProduct(Nu, GTraction, dA);


        // damage field

        double kp      = this->givePenaltyParameter();
        double Delta_t = tStep->giveTimeIncrement();
        //double d       = computeDamageAt( ip, VM_Total, tStep);
        double Delta_d = computeDamageAt(ip, VM_Incremental, tStep);
        double l       = this->giveInternalLength();
        double g_c     = this->giveCriticalEnergy();
        double Gprim   = computeGPrim(ip, VM_Total, tStep);

        double Psibar  = this->computeFreeEnergy( ip, tStep );
        double gradd   = Bd.dotProduct(a_d); // Dalpha/Ds

        //double Psibar0 = this->givePsiBar0();
        //double sectionalForcesScal =  g_c / l * d + Gprim * this->MaCauley(Psibar-Psibar0);
        double sectionalForcesScal = -kp * neg_MaCauley(Delta_d/Delta_t) + g_c / l * d + Gprim * Psibar;
        //double sectionalForcesScal = g_c / l * d + Gprim * Psibar;
	    double sectionalForcesVec  = g_c * l * gradd;
        fd = fd + ( Nd*sectionalForcesScal + Bd*sectionalForcesVec ) * dA;
    }
    fd4.resize(4);
    fd4.zero();
    IntArray indx1, indx2;
    indx1 = {1, 2};
    indx2 = {3, 4};
    fd4.assemble(fd, indx1);
    fd4.assemble(fd, indx2);
}


double
IntElLine1PF :: computeFreeEnergy(GaussPoint *gp, TimeStep *tStep)
{
    StructuralInterfaceMaterialStatus *matStat = static_cast< StructuralInterfaceMaterialStatus * >( gp->giveMaterialStatus() );
    FloatArray strain, stress;
    stress = matStat->giveTempFirstPKTraction();
    strain = matStat->giveTempJump();
    //stress = matStat->giveFirstPKTraction();
    //strain = matStat->giveJump();
    //stress.printYourself();
    //strain.printYourself();

    // Damage only for positive normal traction
    if ( strain.at(3) < 0.0 ) {
        return 0.5 * ( stress.at(1)*strain.at(1) + stress.at(2)*strain.at(2) );
    } else {
        return 0.5 * stress.dotProduct( strain );
    }
}

double
IntElLine1PF :: computeOldG(GaussPoint *gp, ValueModeType valueMode, TimeStep *stepN)
{
    // computes Dg/Dd = (1-d)^2 + r0
    //double d = this->computeDamageAt(gp, VM_Total, stepN) - this->computeDamageAt(gp, VM_Incremental, stepN);
    double d = this->computeDamageAt(gp, VM_Total, stepN);
    double r0 = 1.0e-10;
    return (1.0 - d) * (1.0 - d) + r0;
}


void
IntElLine1PF :: giveDofManDofIDMask_u(IntArray &answer)
{
	//StructuralInterfaceElement :: giveDofManDofIDMask(-1, EID_MomentumBalance, answer); 
    //IntElLine1 :: giveDofManDofIDMask(-1, EID_MomentumBalance, answer);
    answer = {D_u, D_v};
}

void
IntElLine1PF :: giveDofManDofIDMask_d(IntArray &answer)
{
    answer = {T_f};
}


void
IntElLine1PF :: computeLocationArrayOfDofIDs( const IntArray &dofIdArray, IntArray &answer )
{
    // Routine to extract compute the location array an element given an dofid array.
    answer.resize( 0 );
    //NLStructuralElement *el = this->giveElement();
    StructuralInterfaceElement *el = this->giveElement();
    int k = 0;
    for(int i = 1; i <= el->giveNumberOfDofManagers(); i++) {
        DofManager *dMan = el->giveDofManager( i );
        for(int j = 1; j <= dofIdArray.giveSize( ); j++) {

            if(dMan->hasDofID( (DofIDItem) dofIdArray.at( j ) )) {
                Dof *d = dMan->giveDofWithID( dofIdArray.at( j ) );
                //answer.followedBy( k + d->giveNumber( ) );
            }
        }
        k += dMan->giveNumberOfDofs( );
    }
}

void
IntElLine1PF :: computeNd_vectorAt(const FloatArray &lCoords, FloatArray &N)
{
    StructuralInterfaceElement *el = this->giveElement();
    FloatArray Nvec;
    el->giveInterpolation( )->evalN( N, lCoords, FEIElementGeometryWrapper( el ) );

}

void
IntElLine1PF :: computeBd_vectorAt(GaussPoint *aGaussPoint, FloatArray &answer)
{
    // Returns the [numSpaceDim x nDofs] gradient matrix {B_d} of the receiver,
    // evaluated at gp.

    StructuralInterfaceElement *el = dynamic_cast< StructuralInterfaceElement* > ( this->giveElement() );
    FloatMatrix dNdxi;
    this->giveInterpolation()->evaldNdxi( dNdxi, aGaussPoint->giveNaturalCoordinates( ), FEIElementGeometryWrapper( this ) );

    answer.resize(2);
    for (int i = 1; i <= dNdxi.giveNumberOfRows(); i++) {
        answer.at(i) = dNdxi.at(i,1);
    }
    FloatArray G;
    this->computeCovarBaseVectorAt(aGaussPoint, G);
    answer.times( 1.0 / sqrt(G.dotProduct(G)) );
}

} // end namespace oofem
