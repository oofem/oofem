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

#include "shell7base.h"
#include "node.h"
#include "load.h"
#include "structuralms.h"
#include "mathfem.h"
#include "domain.h"
#include "equationid.h"
#include "gaussintegrationrule.h"
#include "gausspnt.h"
#include "feinterpol3d.h"
#include "boundaryload.h"
#include "constantpressureload.h"
#include "constantsurfaceload.h"

namespace oofem {
Shell7Base :: Shell7Base(int n, Domain *aDomain) : NLStructuralElement(n, aDomain),  LayeredCrossSectionInterface()
{}

IRResultType Shell7Base :: initializeFrom(InputRecord *ir)
{
    this->NLStructuralElement :: initializeFrom(ir);
    this->setupInitialNodeDirectors();
    return IRRT_OK;
}

int 
Shell7Base :: checkConsistency()
{
    NLStructuralElement :: checkConsistency();

    this->layeredCS = static_cast< LayeredCrossSection * >( this->giveCrossSection() );

    return 1;
}




Interface *Shell7Base :: giveInterface(InterfaceType it)
{
    switch ( it ) {
    case NodalAveragingRecoveryModelInterfaceType:
        return static_cast< NodalAveragingRecoveryModelInterface * >( this );

    case LayeredCrossSectionInterfaceType:
        return static_cast< LayeredCrossSectionInterface * >( this );

    default:
        return StructuralElement :: giveInterface(it);
    }
}


void
Shell7Base :: giveDofManDofIDMask(int inode, EquationID ut, IntArray &answer) const
{
    answer.setValues(7, D_u, D_v, D_w, W_u, W_v, W_w, Gamma);
}


int
Shell7Base :: giveNumberOfDofs() 
{
    return giveNumberOfFieldDofs(All);
}

int
Shell7Base :: computeGlobalCoordinates(FloatArray &answer, const FloatArray &lcoords)
{
    // should it return coord in reference or updated config?
    return 1;
}


double
Shell7Base :: giveGlobalZcoord(GaussPoint *gp)
{
    return ( * gp->giveCoordinates() ).at(3) * this->giveCrossSection()->give(CS_Thickness) * 0.5;
}




// Base vectors

#if 1

void
Shell7Base :: evalInitialCovarBaseVectorsAt(GaussPoint *gp, FloatArray &G1, FloatArray &G2, FloatArray &G3)
{
    double x, y, z, Mx, My, Mz, zeta;
    FloatArray lcoords = * gp->giveCoordinates();
    zeta = giveGlobalZcoord(gp);
    FloatArray N, M;
    FloatMatrix dNdxi, Mmat;

    // In plane base vectors
    FEInterpolation3d *fei = static_cast< FEInterpolation3d * >( this->giveInterpolation() );
    fei->evaldNdxi( dNdxi, lcoords, FEIElementGeometryWrapper(this) );

    G1.resize(3);
    G2.resize(3);
    G1.zero();
    G2.zero();

    for ( int i = 1; i <= this->giveNumberOfDofManagers(); i++ ) {
        FloatArray *nodeI = this->giveNode(i)->giveCoordinates();
        x = nodeI->at(1);
        y = nodeI->at(2);
        z = nodeI->at(3);

        M = this->giveInitialNodeDirector(i);
        Mx = M.at(1);
        My = M.at(2);
        Mz = M.at(3);

        G1.at(1) += dNdxi.at(i, 1) * ( x + zeta * Mx );
        G1.at(2) += dNdxi.at(i, 1) * ( y + zeta * My );
        G1.at(3) += dNdxi.at(i, 1) * ( z + zeta * Mz );
        G2.at(1) += dNdxi.at(i, 2) * ( x + zeta * Mx );
        G2.at(2) += dNdxi.at(i, 2) * ( y + zeta * My );
        G2.at(3) += dNdxi.at(i, 2) * ( z + zeta * Mz );
    }

    // Out of plane base vector = director
    this->evalInitialDirectorAt(gp, G3);     // G3=M
}


void
Shell7Base :: edgeEvalInitialCovarBaseVectorsAt(GaussPoint *gp, const int iedge, FloatArray &G1, FloatArray &G3)
{
    double x, y, z, Mx, My, Mz, zeta;
    FloatArray lcoords = * gp->giveCoordinates();
    zeta = 0.0;     // no variation i z (yet)
    FloatArray N, M, dNdxi;
    FloatMatrix Mmat;
    IntArray edgeNodes;
    FEInterpolation3d *fei = static_cast< FEInterpolation3d * >( this->giveInterpolation() );

    fei->computeLocalEdgeMapping(edgeNodes, iedge);
    fei->edgeEvaldNdxi( dNdxi, iedge, lcoords, FEIElementGeometryWrapper(this) );

    // Base vector along edge
    G1.resize(3);
    G1.zero();
    for ( int i = 1; i <= edgeNodes.giveSize(); i++ ) {
        FloatArray *nodeI = this->giveNode( edgeNodes.at(i) )->giveCoordinates();
        x = nodeI->at(1);
        y = nodeI->at(2);
        z = nodeI->at(3);
        M = this->giveInitialNodeDirector( edgeNodes.at(i) );
        Mx = M.at(1);
        My = M.at(2);
        Mz = M.at(3);
        G1.at(1) += dNdxi.at(i) * ( x + zeta * Mx );
        G1.at(2) += dNdxi.at(i) * ( y + zeta * My );
        G1.at(3) += dNdxi.at(i) * ( z + zeta * Mz );
    }

    // Director will be the second base vector
    this->edgeEvalInitialDirectorAt(gp, G3, iedge);
}

void
Shell7Base :: evalInitialContravarBaseVectorsAt(GaussPoint *gp, FloatArray &Gcon1, FloatArray &Gcon2, FloatArray &Gcon3)
{
    FloatArray Gcov1, Gcov2, Gcov3;
    this->evalInitialCovarBaseVectorsAt(gp, Gcov1, Gcov2, Gcov3);
    this->giveDualBase(Gcov1, Gcov2, Gcov3, Gcon1, Gcon2, Gcon3);
}

void
Shell7Base :: evalContravarBaseVectorsAt(GaussPoint *gp, FloatArray &gcon1, FloatArray &gcon2, FloatArray &gcon3, FloatArray &solVec)
{
    FloatArray gcov1, gcov2, gcov3;
    this->evalCovarBaseVectorsAt(gp, gcov1, gcov2, gcov3, solVec);
    this->giveDualBase(gcov1, gcov2, gcov3, gcon1, gcon2, gcon3);
}

void
Shell7Base :: giveDualBase(const FloatArray &G1, const FloatArray &G2, const FloatArray &G3, FloatArray &g1, FloatArray &g2, FloatArray &g3)
{
    FloatMatrix gmat, ginv, test;

    // Metric tensor
    gmat.resize(3, 3);
    gmat.at(1, 1) = G1.dotProduct(G1);
    gmat.at(1, 2) = G1.dotProduct(G2);
    gmat.at(1, 3) = G1.dotProduct(G3);
    gmat.at(2, 2) = G2.dotProduct(G2);
    gmat.at(2, 3) = G2.dotProduct(G3);
    gmat.at(3, 3) = G3.dotProduct(G3);
    gmat.symmetrized();

    ginv.beInverseOf(gmat);
    g1.resize(3);
    g1.zero();
    g2.resize(3);
    g2.zero();
    g3.resize(3);
    g3.zero();

    g1.add(ginv.at(1, 1), G1);
    g1.add(ginv.at(1, 2), G2);
    g1.add(ginv.at(1, 3), G3);
    g2.add(ginv.at(2, 1), G1);
    g2.add(ginv.at(2, 2), G2);
    g2.add(ginv.at(2, 3), G3);
    g3.add(ginv.at(3, 1), G1);
    g3.add(ginv.at(3, 2), G2);
    g3.add(ginv.at(3, 3), G3);
}

void
Shell7Base :: evalInitialDirectorAt(GaussPoint *gp, FloatArray &answer)
{       // Interpolates between the node directors
    FloatArray &lcoords = * gp->giveCoordinates();
    FloatArray N;
    FEInterpolation3d *fei = static_cast< FEInterpolation3d * >( this->giveInterpolation() );
    fei->evalN( N, lcoords, FEIElementGeometryWrapper(this) );

    answer.resize(3);
    answer.zero();
    for ( int i = 1; i <= this->giveNumberOfDofManagers(); i++ ) {
        answer.add( N.at(i), this->giveInitialNodeDirector(i) );
    }
}


void
Shell7Base :: edgeEvalInitialDirectorAt(GaussPoint *gp, FloatArray &answer, const int iEdge)
{       // Interpolates between the node directors
    FloatArray &lcoords = * gp->giveCoordinates();
    FloatArray N;
    IntArray edgeNodes;

    FEInterpolation3d *fei = static_cast< FEInterpolation3d * >( this->giveInterpolation() );

    fei->computeLocalEdgeMapping(edgeNodes, iEdge);
    fei->edgeEvalN( N, lcoords, FEIElementGeometryWrapper(this) );

    answer.resize(3);
    answer.zero();
    for ( int i = 1; i <= edgeNodes.giveSize(); i++ ) {
        answer.add( N.at(i), this->giveInitialNodeDirector( edgeNodes.at(i) ) );
    }
}


void
Shell7Base :: setupInitialNodeDirectors()
{   /* If the directors are not present in the input file, then they should be approximated as
     * normal to the initial surface. (No support in the input file at the moment.)
     */
    // Compute directors as normals to the surface

    FloatArray M(3), G1(3), G2(3), lcoords(2), nodeLocalXiCoords, nodeLocalEtaCoords;

    // Give the local coordinates for the element nodes - all at once
    this->giveLocalNodeCoords(nodeLocalXiCoords, nodeLocalEtaCoords);
    int nDofMan = this->giveNumberOfDofManagers();
    this->initialNodeDirectors.resize(nDofMan);

    for ( int node = 1; node <= nDofMan; node++ ) {
        this->initialNodeDirectors [ node - 1 ].resize(3);
        this->initialNodeDirectors [ node - 1 ].zero();

        lcoords.at(1) = nodeLocalXiCoords.at(node);
        lcoords.at(2) = nodeLocalEtaCoords.at(node);
        FEInterpolation3d *fei = static_cast< FEInterpolation3d * >( this->giveInterpolation() );
        FloatMatrix dNdxi;
        fei->evaldNdxi( dNdxi, lcoords, FEIElementGeometryWrapper(this) );

        G1.zero();
        G2.zero();
        M.zero();
        for ( int i = 1; i <= nDofMan; i++ ) {         // base vectors of the initial surface
            FloatArray *nodeI = this->giveNode(i)->giveCoordinates();
            G1.add(dNdxi.at(i, 1), * nodeI);
            G2.add(dNdxi.at(i, 2), * nodeI);
        }

        M.beVectorProductOf(G1, G2);
        M.normalize();
        this->initialNodeDirectors [ node - 1 ].add(M);
    }
}

void
Shell7Base :: evalCovarBaseVectorsAt(GaussPoint *gp, FloatArray &g1, FloatArray &g2, FloatArray &g3, FloatArray &genEps)
{
    FloatArray lcoords = * gp->giveCoordinates();
    double zeta = giveGlobalZcoord(gp);

    FloatArray dxdxi, dxdxi1, dxdxi2, m, dmdxi, dmdxi1, dmdxi2, dgamdxi,  test;
    double dgamdxi1, dgamdxi2, gam;
    this->giveGeneralizedStrainComponents(genEps, dxdxi1, dxdxi2, dmdxi1, dmdxi2, m, dgamdxi1, dgamdxi2, gam);

    double fac1 = ( zeta + 0.5 * gam * zeta * zeta );
    double fac2 = ( 0.5 * zeta * zeta );
    double fac3 = ( 1.0 + zeta * gam );

    g1 = dxdxi1;
    g1.add(fac1, dmdxi1);
    g1.add(fac2 * dgamdxi1, m);
    g2 = dxdxi2;
    g2.add(fac1, dmdxi2);
    g2.add(fac2 * dgamdxi2, m);
    g3 = m;
    g3.times(fac3);             //g3.add(fac3,m);
}

void
Shell7Base :: edgeEvalCovarBaseVectorsAt(GaussPoint *gp, const int iedge, FloatArray &g1, FloatArray &g3, TimeStep *tStep)
{
    double zeta = 0.0, gam;
    FloatArray lcoords = * gp->giveCoordinates();
    FloatArray a;
    FloatMatrix B;
    IntArray edgeNodes;
    FEInterpolation3d *fei = static_cast< FEInterpolation3d * >( this->giveInterpolation() );
    fei->computeLocalEdgeMapping(edgeNodes, iedge);
    this->edgeComputeBmatrixAt(gp, B, 1, ALL_STRAINS);
    this->edgeGiveUpdatedSolutionVector(a, iedge, tStep);

    FloatArray eps;           // generalized strain
    eps.beProductOf(B, a);    // [dxdxi, dmdxi, m, dgamdxi, gam]^T


    FloatArray dxdxi, m, dmdxi, test;
    double dgamdxi;
    dxdxi.setValues( 3, eps.at(1), eps.at(2), eps.at(3) );
    dmdxi.setValues( 3, eps.at(4), eps.at(5), eps.at(6) );
    m.setValues( 3, eps.at(7), eps.at(8), eps.at(9) );
    dgamdxi = eps.at(10);
    gam = eps.at(11);

    g1.resize(3),  g3.resize(3);
    double fac1 = ( zeta + 0.5 * gam * zeta * zeta );
    double fac2 = ( 0.5 * zeta * zeta );
    double fac3 = ( 1.0 + zeta * gam );

    g1.at(1) = dxdxi.at(1) + fac1 * dmdxi.at(1) + fac2 * m.at(1) * dgamdxi;
    g1.at(2) = dxdxi.at(2) + fac1 * dmdxi.at(2) + fac2 * m.at(2) * dgamdxi;
    g1.at(3) = dxdxi.at(3) + fac1 * dmdxi.at(3) + fac2 * m.at(3) * dgamdxi;

    g3.at(1) = fac3 * m.at(1);
    g3.at(2) = fac3 * m.at(2);
    g3.at(3) = fac3 * m.at(3);
}

#endif




// Tangent matrices

#if 1

void
Shell7Base :: computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
{
    //FloatMatrix massNum, mass;
    //this->computeMassMatrixNum(massNum, tStep);
    //this->computeMassMatrix(mass, tStep);

    int ndofs = this->giveNumberOfDofs();
    answer.resize(ndofs, ndofs);
    answer.zero();

    FloatMatrix temp;
    FloatArray solVec, totalSolVec;
    //this->giveUpdatedSolutionVector(solVec, tStep);

    this->giveUpdatedSolutionVector(totalSolVec, tStep);  //  a
    this->giveUpdatedSolutionVector(solVec, tStep); // da

    //this->computeBulkTangentMatrix(answer, solVec, rMode, tStep);
    //this->new_computeBulkTangentMatrix(temp, solVec, solVec, solVec, rMode, tStep);
    this->new_computeBulkTangentMatrix(answer, totalSolVec, solVec, solVec, rMode, tStep);

    // Add contribution due to pressure load
    int nLoads = this->boundaryLoadArray.giveSize() / 2;

    for ( int i = 1; i <= nLoads; i++ ) {     // For each pressure load that is applied
        int load_number = this->boundaryLoadArray.at(2 * i - 1);
        int iSurf = this->boundaryLoadArray.at(2 * i);         // load_id
        Load *load;
        load = this->domain->giveLoad(load_number);

        if ( dynamic_cast< ConstantPressureLoad * >( load ) ) {
            FloatMatrix K_pressure;
            this->computePressureTangentMatrix(K_pressure, load, iSurf, tStep);
            //answer.add(K_pressure); // Should assemble with ordering
        }
    }

}

void
Shell7Base :: computeLambdaMatrices(FloatMatrix lambda [ 3 ], FloatArray &genEps, double zeta)
{

    FloatArray m(3), dm1(3), dm2(3), temp1, temp2;
    double dgam1, dgam2, gam;
    this->giveGeneralizedStrainComponents(genEps, temp1, temp1, dm1, dm2, m, dgam1, dgam2, gam);

    // thickness coefficients
    double a = zeta + 0.5 * gam * zeta * zeta;
    double b = 0.5 * zeta * zeta;
    double c = 1. + gam * zeta;

    // lambda1 =  ( I,   0,  a*I,   0 ,  b*dgam1*I,  b*m,   0 ,  b*dm1 )
    // lambda2 =  ( 0,   I,   0 ,  a*I,  b*dgam2*I,   0 ,  b*m,  b*dm2 )

    FloatMatrix eye(3,3), aEye(3,3);
    eye.beUnitMatrix();
    aEye=eye;  aEye.times(a);
    lambda[ 0 ].resize(3,18);   lambda[ 0 ].zero();  
    lambda[ 1 ].resize(3,18);   lambda[ 1 ].zero();

    lambda[ 0 ].setSubMatrix(eye,1,1);  lambda[ 0 ].setSubMatrix(aEye,1,7);
    lambda[ 1 ].setSubMatrix(eye,1,4);  lambda[ 1 ].setSubMatrix(aEye,1,10);

    FloatMatrix bdg1eye(3,3), bdg2eye(3,3);
    bdg1eye = eye;  bdg1eye.times(b*dgam1);
    bdg2eye = eye;  bdg2eye.times(b*dgam2);
    lambda[ 0 ].setSubMatrix(bdg1eye,1,13);
    lambda[ 1 ].setSubMatrix(bdg2eye,1,13);

    FloatArray bm(3), bdm1(3), bdm2(3);
    bm = m; bm.times(b);
    bdm1 = dm1; bdm1.times(b);
    bdm2 = dm2; bdm2.times(b);
    lambda[ 0 ].setColumn(bm,16);   lambda[ 0 ].setColumn(bdm1,18);
    lambda[ 1 ].setColumn(bm,17);   lambda[ 1 ].setColumn(bdm2,18);

    // lambda3 =  ( 0,   0,   0 ,   0 ,     c*I   ,   0 ,   0 ,   xi*m )
    lambda[ 2 ].resize(3,18);   lambda[ 2 ].zero();
    lambda[ 2 ].at(1,13) = lambda[ 2 ].at(2,14) = lambda[ 2 ].at(3,15) = c;
    FloatArray zm(3);
    zm = m; zm.times(zeta);
    lambda[ 2 ].setColumn(zm,18);
}



void
Shell7Base :: new_computeBulkTangentMatrix(FloatMatrix &answer, FloatArray &solVec, FloatArray &solVecI, FloatArray &solVecJ, MatResponseMode rMode, TimeStep *tStep)
{
    FloatMatrix A [ 3 ] [ 3 ], lambdaI [ 3 ], lambdaJ [ 3 ];
    FloatMatrix L(18,18);
    FloatMatrix B11, B22, B32, B43, B53, B;
    FloatArray S1g(3), S2g(3), S3g(3);
    FloatMatrix K(42,42), tempAnswer(42,42);
    K.zero();
    tempAnswer.zero();

    int ndofs = Shell7Base :: giveNumberOfDofs();
    answer.resize(ndofs, ndofs);
    answer.zero();

    int numberOfLayers = this->layeredCS->giveNumberOfLayers();     

    for ( int layer = 1; layer <= numberOfLayers; layer++ ) {
        IntegrationRule *iRule = layerIntegrationRulesArray [ layer - 1 ];
        Material *mat = domain->giveMaterial( this->layeredCS->giveLayerMaterial(layer) );

        for ( int i = 0; i < iRule->getNumberOfIntegrationPoints(); i++ ) {
            GaussPoint *gp = iRule->getIntegrationPoint(i);


            this->computeBmatricesAt(gp, B11, B22, B32, B43, B53);
            FloatArray genEpsI, genEpsJ, genEps;
            this->computeGeneralizedStrainVector(genEpsI, solVecI, B11, B22, B32, B43, B53);
            this->computeGeneralizedStrainVector(genEpsJ, solVecJ, B11, B22, B32, B43, B53);
            this->computeGeneralizedStrainVector(genEps , solVec , B11, B22, B32, B43, B53);

            // Material stiffness
            Shell7Base :: computeLinearizedStiffness(gp, mat, tStep, S1g, S2g, S3g, A, genEps);

            double zeta = this->giveGlobalZcoord(gp);
            this->computeLambdaMatrices(lambdaI, genEpsI, zeta);
            this->computeLambdaMatrices(lambdaJ, genEpsJ, zeta);

            this->computeBmatrixAt(gp, B, 0, 0);
            // L = sum_{i,j} (lambdaI_i)^T * A^ij * lambdaJ_j
            // Naive implementation - should be optimized 
            // note: L will only be symmetric if lambdaI = lambdaJ
            FloatMatrix temp;
            L.zero();
            for ( int j = 0; j < 3; j++ ) {
                for ( int k = 0; k < 3; k++ ) {
                    this->computeTripleProduct(temp, lambdaI [ j ], A [ j ][ k ], lambdaJ [ k ]);
                    L.add(temp);
                }
            }
     
            FloatMatrix Ktemp, K;
            Ktemp.beProductOf(L, B);
            double dV = this->computeVolumeAroundLayer(gp, layer);
            K.beTProductOf(B,Ktemp);
            tempAnswer.add(dV, K);
        }
    }

    
    FloatMatrix test;
    test.beSubMatrixOf(tempAnswer,29,36,29,36);
    //test.printYourself();

    IntArray ordering = giveOrdering(All);
    answer.assemble(tempAnswer, ordering, ordering);


}



void
Shell7Base :: computeBulkTangentMatrix(FloatMatrix &answer, FloatArray &solVec, MatResponseMode rMode, TimeStep *tStep)
{
    FloatMatrix A [ 3 ] [ 3 ];
    FloatMatrix F1 [ 3 ];
    FloatArray t1(3), t2(3), t3(3);
    FloatArray f3 [ 2 ], f4, f5;

    FloatMatrix L11(6, 6), L12(6, 6), L13(6, 3), L14(6, 2), L15(6, 1),
    L22(6, 6), L23(6, 3), L24(6, 2), L25(6, 1),
    L33(3, 3), L34(3, 2), L35(3, 1),
    L44(2, 2), L45(2, 1),
    L55(1, 1);
    FloatMatrix B11, B22, B32, B43, B53;
    FloatArray S1g(3), S2g(3), S3g(3), m(3), dm1(3), dm2(3), temp1, temp2;
    FloatMatrix K11(18, 18), K12(18, 18), K13(18, 6), K22(18, 18), K23(18, 6), K33(6, 6);
    K11.zero(), K12.zero(), K13.zero(), K22.zero(), K23.zero(), K33.zero();

    FloatMatrix K11temp1, K12temp1, K13temp1, K22temp1, K22temp2, K23temp1, K23temp2, K33temp1, K33temp2;

    //FloatArray solVec;
    //this->giveUpdatedSolutionVector(solVec, tStep);

    //int ndofs = this->giveNumberOfDofs();
    //answer.resize(ndofs, ndofs);
    //answer.zero();

    LayeredCrossSection *layeredCS = dynamic_cast< LayeredCrossSection * >( this->giveCrossSection() );
    int numberOfLayers = layeredCS->giveNumberOfLayers();     // conversion from double to int - fix!

    for ( int layer = 1; layer <= numberOfLayers; layer++ ) {
        IntegrationRule *iRule = layerIntegrationRulesArray [ layer - 1 ];
        Material *mat = domain->giveMaterial( layeredCS->giveLayerMaterial(layer) );

        for ( int i = 0; i < iRule->getNumberOfIntegrationPoints(); i++ ) {
            GaussPoint *gp = iRule->getIntegrationPoint(i);

            double gam, dg1, dg2;

            this->computeBmatricesAt(gp, B11, B22, B32, B43, B53);
            FloatArray genEps;
            this->computeGeneralizedStrainVector(genEps, solVec, B11, B22, B32, B43, B53);
            this->giveGeneralizedStrainComponents(genEps, temp1, temp1, dm1, dm2, m, dg1, dg2, gam);


            // Material stiffness
            Shell7Base :: computeLinearizedStiffness(gp, mat, tStep, S1g, S2g, S3g, A, genEps);

            // Tangent stiffness
 #if 1
            // thickness coefficients
            double zeta = giveGlobalZcoord(gp);
            double a = zeta + 0.5 * gam * zeta * zeta;
            double b = 0.5 * zeta * zeta;
            double c = 1. + gam * zeta;

            // f1(alpha) = b*A(alpha,beta)*dg(beta) + c*A(alpha,3);
            F1 [ 0 ].resize(0, 0);
            F1 [ 1 ].resize(0, 0);
            F1 [ 2 ].resize(0, 0);
            F1 [ 0 ].add(dg1 * b, A [ 0 ] [ 0 ]);
            F1 [ 0 ].add(dg2 * b, A [ 0 ] [ 1 ]);
            F1 [ 0 ].add(c, A [ 0 ] [ 2 ]);
            F1 [ 1 ].add(dg1 * b, A [ 1 ] [ 0 ]);
            F1 [ 1 ].add(dg2 * b, A [ 1 ] [ 1 ]);
            F1 [ 1 ].add(c, A [ 1 ] [ 2 ]);
            F1 [ 2 ].add(dg1 * b, A [ 2 ] [ 0 ]);
            F1 [ 2 ].add(dg2 * b, A [ 2 ] [ 1 ]);
            F1 [ 2 ].add(c, A [ 2 ] [ 2 ]);

            //f3(alpha) = b*A(alpha,beta)*dm(beta) + zeta*A(alpha,3)*m;
            f3 [ 0 ].resize(0);
            f3 [ 1 ].resize(0);
            t1.beProductOf(A [ 0 ] [ 0 ], dm1);
            t2.beProductOf(A [ 0 ] [ 1 ], dm2);
            t3.beProductOf(A [ 0 ] [ 2 ], m);
            f3 [ 0 ].add(b, t1);
            f3 [ 0 ].add(b, t2);
            f3 [ 0 ].add(zeta, t3);

            //f32 = b*( A.at(2,1)*dm1 + A.at(2,2)*dm2 )  + zeta*A.at(2,3)*m;
            t1.beProductOf(A [ 1 ] [ 0 ], dm1);
            t2.beProductOf(A [ 1 ] [ 1 ], dm2);
            t3.beProductOf(A [ 1 ] [ 2 ], m);
            f3 [ 1 ].add(b, t1);
            f3 [ 1 ].add(b, t2);
            f3 [ 1 ].add(zeta, t3);


            //f4 = b*dm(alpha)*A(alpha,3) + zeta*A(3,3)*m;
            f4.resize(0);
            t1.beTProductOf(A [ 0 ] [ 2 ], dm1);
            t2.beTProductOf(A [ 1 ] [ 2 ], dm2);
            t3.beTProductOf(A [ 2 ] [ 2 ], m);
            f4.add(b, t1);
            f4.add(b, t2);
            f4.add(zeta, t3);

            // f5 = b*F1(alpha)*dm(alpha) + zeta*F2*m + zeta*N3
            f5.resize(0);
            t1.beTProductOf(F1 [ 0 ], dm1);
            t2.beTProductOf(F1 [ 1 ], dm2);
            t3.beTProductOf(F1 [ 2 ], m);
            f5.add(b, t1);
            f5.add(b, t2);
            f5.add(zeta, t3);
            f5.add(zeta, S3g);


            /*
             *
             *     A     a*A            F1           b*A*m             f3
             *    a*A    a^2*A          a*F1         a*b*A*m         a*f3+b*N
             *    F1     a*F1     b*F1*dgam+c*F2     b*(F1*m+N)        f5
             *    b*A*m  a*b*A*m    b*(F1*m+N)       b^2*m*A*m     b*m*f3
             *    f3     a*f3+b*N       f5           b*m*f3    b*dm*f3+zeta*m*f4
             */

            // L(1,1) = A(alpha,beta)
            L11.setSubMatrix(A [ 0 ] [ 0 ], 1, 1);
            L11.setSubMatrix(A [ 0 ] [ 1 ], 1, 4);
            L11.setSubMatrix(A [ 1 ] [ 0 ], 4, 1);
            L11.setSubMatrix(A [ 1 ] [ 1 ], 4, 4);

            // L(1,2) = a*A(alpha,beta)
            L12 = L11;
            L12.times(a);

            // L(1,3) = F1(alpha)
            L13.setSubMatrix(F1 [ 0 ], 1, 1);
            L13.setSubMatrix(F1 [ 1 ], 4, 1);

            // L(1,4) = b*A*m
            L14.zero();
            t1.beProductOf(A [ 0 ] [ 0 ], m);
            t2.beProductOf(A [ 0 ] [ 1 ], m);
            L14.addSubVectorCol(t1, 1, 1);
            L14.addSubVectorCol(t2, 1, 2);
            t1.beProductOf(A [ 1 ] [ 0 ], m);
            t2.beProductOf(A [ 1 ] [ 1 ], m);
            L14.addSubVectorCol(t1, 4, 1);
            L14.addSubVectorCol(t2, 4, 2);
            L14.times(b);

            // L(1,5) = f3(alpha)
            L15.zero();
            L15.addSubVectorCol(f3 [ 0 ], 1, 1);
            L15.addSubVectorCol(f3 [ 1 ], 4, 1);


            // L(2,2) = a^2*A(alpha,beta)= a*L12
            L22 = L12;
            L22.times(a);

            // L(2,3) = a*F1(alpha) = a*L13
            L23 = L13;
            L23.times(a);

            // L(2,4) = a*F1(alpha)=a*L14
            L24 = L14;
            L24.times(a);

            // L(2,5) = a*f3(alpha) + b*N(alpha)
            L25.zero();
            L25.addSubVectorCol(S1g, 1, 1);
            L25.addSubVectorCol(S2g, 4, 1);
            L25.times(b);
            L25.add(a, L15);


            // L(3,3) = b*F1(beta)*dgam(beta) + c*F2
            L33.resize(0, 0);
            L33.add(dg1 * b, F1 [ 0 ]);
            L33.add(dg2 * b, F1 [ 1 ]);
            L33.add(c, F1 [ 2 ]);


            // L(3,4) = b*( F1(beta)*m + N(beta) )
            t1.beTProductOf(F1 [ 0 ], m);
            t1.add(S1g);
            t1.times(b);
            t2.beTProductOf(F1 [ 1 ], m);
            t2.add(S2g);
            t2.times(b);
            L34.setColumn(t1, 1);
            L34.setColumn(t2, 2);

            // L(3,5) = f5
            L35.setColumn(f5, 1);

            // L(4,4) = b^2*m*A*m (2x2)
            t1.beProductOf(A [ 0 ] [ 0 ], m);
            L44.at(1, 1) = t1.dotProduct(m);
            t1.beProductOf(A [ 0 ] [ 1 ], m);
            L44.at(1, 2) = t1.dotProduct(m);
            t1.beProductOf(A [ 1 ] [ 0 ], m);
            L44.at(2, 1) = t1.dotProduct(m);
            t1.beProductOf(A [ 1 ] [ 1 ], m);
            L44.at(2, 2) = t1.dotProduct(m);
            L44.times(b * b);

            // L(4,5) = b*m*f3(alpha)
            L45.at(1, 1) = b * m.dotProduct(f3 [ 0 ]);
            L45.at(2, 1) = b * m.dotProduct(f3 [ 1 ]);

            // L(5,5) = b*m*f3(alpha)
            L55.at(1, 1) = b * dm1.dotProduct(f3 [ 0 ]) + b * dm2.dotProduct(f3 [ 1 ]) + zeta * m.dotProduct(f4);
#endif
            // K11 = BT11*L11*B11
            // K11 = BT11*K11temp1
            K11temp1.beProductOf(L11, B11);

            // K12 = BT11*(L12*B22+L13*B32)
            // K12 = BT11*K12temp1
            K12temp1.beProductOf(L12, B22);
            K12temp1.addProductOf(L13, B32);

            // K13 = BT11*(L14*B43 + L15*B53)
            // K13 = BT11*K13temp1
            K13temp1.beProductOf(L14, B43);
            K13temp1.addProductOf(L15, B53);

            // K22 = BT22*(L22*B22 + L23*B32) + BT32*(L32*B22 + L33*B32)
            // K22 = BT22*K22temp1            + BT32*K22temp2
            K22temp1.beProductOf(L22, B22);
            K22temp1.addProductOf(L23, B32);
            K22temp2.beTProductOf(L23, B22);
            K22temp2.addProductOf(L33, B32);

            // K23 = BT22*(L24*B43 + L25*B53) + BT32*(L34*B43 + L35*B53)
            // K23 = BT22*K23temp1            + BT32*K23temp2
            K23temp1.beProductOf(L24, B43);
            K23temp1.addProductOf(L25, B53);
            K23temp2.beProductOf(L34, B43);
            K23temp2.addProductOf(L35, B53);
            // K23 = (BT22*L24 + BT32*L34)*B43 + (BT22*L25 + BT32*L35)*B53 // TODO: This order would be faster
            //K23temp1.TbeProductOf(B22,L24);
            //K23temp1.addTProductOf(B32,L34);
            //K23temp2.beTProductOf(L22,L25);
            //K23temp2.addTProductOf(B32,L35);


            // K33 = BT43*(L44*B43 + L45*B53) + BT53*(L54*B43 + L55*B53)
            // K33 = BT43*K33temp1            + BT53*K33temp2
            K33temp1.beProductOf(L44, B43);
            K33temp1.addProductOf(L45, B53);
            K33temp2.beTProductOf(L45, B43);
            K33temp2.add(L55.at(1, 1), B53);            // L55 is just a scalar (and K33temp2 is just an array)

            double dV = this->computeVolumeAroundLayer(gp, layer);

            K11.plusProductSymmUpper(B11, K11temp1, dV);
            K12.plusProductUnsym(B11, K12temp1, dV);
            K13.plusProductUnsym(B11, K13temp1, dV);

            K22.plusProductSymmUpper(B22, K22temp1, dV);
            K22.plusProductSymmUpper(B32, K22temp2, dV);

            K23.plusProductUnsym(B22, K23temp1, dV);
            K23.plusProductUnsym(B32, K23temp2, dV);
            // Potential minor optimization, actually use FloatArray for anything that only has a single column
            //tmpA.beProductOf(K23temp1, B43); K23.add(dV, tmpA);
            //K23.plusDyadUnsymm(K33temp2, B53, dV);

            K33.plusProductSymmUpper(B43, K33temp1, dV);
            K33.plusProductSymmUpper(B53, K33temp2, dV);
            //K33.plusDyadSymmUpper(B53, K33temp2, dV); // Potential minor optimization (probably not worth it)
        }
    }

    K11.symmetrized();
    K22.symmetrized();
    K33.symmetrized();

    IntArray ordering_phibar = giveOrdering(Midplane);
    IntArray ordering_m = giveOrdering(Director);
    IntArray ordering_gam = giveOrdering(InhomStrain);

    answer.assemble(K11, ordering_phibar, ordering_phibar);
    answer.assemble(K12, ordering_phibar, ordering_m);
    answer.assemble(K13, ordering_phibar, ordering_gam);
    answer.assemble(K22, ordering_m,      ordering_m);
    answer.assemble(K23, ordering_m,      ordering_gam);
    answer.assemble(K33, ordering_gam,    ordering_gam);

    FloatMatrix K21, K31, K32;
    K21.beTranspositionOf(K12);
    K31.beTranspositionOf(K13);
    K32.beTranspositionOf(K23);
    answer.assemble(K21, ordering_m,      ordering_phibar);
    answer.assemble(K31, ordering_gam,    ordering_phibar);
    answer.assemble(K32, ordering_gam,    ordering_m);
}


void
Shell7Base :: computeLinearizedStiffness(GaussPoint *gp, Material *mat, TimeStep *tStep,
                                         FloatArray &S1g, FloatArray &S2g, FloatArray &S3g, FloatMatrix A [ 3 ] [ 3 ], FloatArray &genEps) {
    // Fix material for layered cross section

    FloatArray lcoords, cartStressVector, contravarStressVector, sectionalForces, BF, a;
    FloatArray g1, g2, g3, m(3), dm1(3), dm2(3), temp1, temp2;

    FloatMatrix D, Dcart, Mmat, S;

    //A = L^iklj * (g_k x g_l) + S^ij*I
    mat->giveCharacteristicMatrix(Dcart, FullForm, TangentStiffness, gp, tStep);     // L_ijkl - cartesian system (Voigt)
    this->transInitialCartesianToInitialContravar(gp, Dcart, D);      // L^ijkl - curvilinear system (Voigt)



    this->computeStressVector(cartStressVector, genEps, gp, mat, tStep);
    this->transInitialCartesianToInitialContravar(gp, cartStressVector, contravarStressVector);
    S.beMatrixForm(contravarStressVector);

    this->computeStressResultantsAt(gp, contravarStressVector, S1g, S2g, S3g, genEps);
    this->evalCovarBaseVectorsAt(gp, g1, g2, g3, genEps);

    FloatMatrix gg11, gg12, gg13, gg21, gg22, gg23, gg31, gg32, gg33;


    gg11.beDyadicProductOf(g1, g1);
    gg12.beDyadicProductOf(g1, g2);
    gg13.beDyadicProductOf(g1, g3);
    gg21.beDyadicProductOf(g2, g1);
    gg22.beDyadicProductOf(g2, g2);
    gg23.beDyadicProductOf(g2, g3);
    gg31.beDyadicProductOf(g3, g1);
    gg32.beDyadicProductOf(g3, g2);
    gg33.beDyadicProductOf(g3, g3);



    // position 11

    A [ 0 ] [ 0 ].resize(3, 3);
    A [ 0 ] [ 0 ].beUnitMatrix();
    A [ 0 ] [ 0 ].times( S.at(1, 1) );
    A [ 0 ] [ 0 ].add(D.at(1, 1), gg11);
    A [ 0 ] [ 0 ].add(D.at(1, 6), gg12);
    A [ 0 ] [ 0 ].add(D.at(1, 5), gg13);
    A [ 0 ] [ 0 ].add(D.at(6, 1), gg21);
    A [ 0 ] [ 0 ].add(D.at(6, 6), gg22);
    A [ 0 ] [ 0 ].add(D.at(6, 5), gg23);
    A [ 0 ] [ 0 ].add(D.at(5, 1), gg31);
    A [ 0 ] [ 0 ].add(D.at(5, 6), gg32);
    A [ 0 ] [ 0 ].add(D.at(5, 5), gg33);


    // position 21
    A [ 1 ] [ 0 ].resize(3, 3);
    A [ 1 ] [ 0 ].beUnitMatrix();
    A [ 1 ] [ 0 ].times( S.at(2, 1) );
    A [ 1 ] [ 0 ].add(D.at(6, 1), gg11);
    A [ 1 ] [ 0 ].add(D.at(6, 6), gg12);
    A [ 1 ] [ 0 ].add(D.at(6, 5), gg13);
    A [ 1 ] [ 0 ].add(D.at(2, 1), gg21);
    A [ 1 ] [ 0 ].add(D.at(2, 6), gg22);
    A [ 1 ] [ 0 ].add(D.at(2, 5), gg23);
    A [ 1 ] [ 0 ].add(D.at(4, 1), gg31);
    A [ 1 ] [ 0 ].add(D.at(4, 6), gg32);
    A [ 1 ] [ 0 ].add(D.at(4, 5), gg33);

    // position 31
    A [ 2 ] [ 0 ].resize(3, 3);
    A [ 2 ] [ 0 ].beUnitMatrix();
    A [ 2 ] [ 0 ].times( S.at(3, 1) );
    A [ 2 ] [ 0 ].add(D.at(5, 1), gg11);
    A [ 2 ] [ 0 ].add(D.at(5, 6), gg12);
    A [ 2 ] [ 0 ].add(D.at(5, 5), gg13);
    A [ 2 ] [ 0 ].add(D.at(4, 1), gg21);
    A [ 2 ] [ 0 ].add(D.at(4, 6), gg22);
    A [ 2 ] [ 0 ].add(D.at(4, 5), gg23);
    A [ 2 ] [ 0 ].add(D.at(3, 1), gg31);
    A [ 2 ] [ 0 ].add(D.at(3, 6), gg32);
    A [ 2 ] [ 0 ].add(D.at(3, 5), gg33);

    // position 12
    A [ 0 ] [ 1 ].resize(3, 3);
    A [ 0 ] [ 1 ].beUnitMatrix();
    A [ 0 ] [ 1 ].times( S.at(1, 2) );
    A [ 0 ] [ 1 ].add(D.at(1, 6), gg11);
    A [ 0 ] [ 1 ].add(D.at(1, 2), gg12);
    A [ 0 ] [ 1 ].add(D.at(1, 4), gg13);
    A [ 0 ] [ 1 ].add(D.at(6, 6), gg21);
    A [ 0 ] [ 1 ].add(D.at(6, 2), gg22);
    A [ 0 ] [ 1 ].add(D.at(6, 4), gg23);
    A [ 0 ] [ 1 ].add(D.at(5, 6), gg31);
    A [ 0 ] [ 1 ].add(D.at(5, 2), gg32);
    A [ 0 ] [ 1 ].add(D.at(5, 4), gg33);

    // position 22
    A [ 1 ] [ 1 ].resize(3, 3);
    A [ 1 ] [ 1 ].beUnitMatrix();
    A [ 1 ] [ 1 ].times( S.at(2, 2) );
    A [ 1 ] [ 1 ].add(D.at(6, 6), gg11);
    A [ 1 ] [ 1 ].add(D.at(6, 2), gg12);
    A [ 1 ] [ 1 ].add(D.at(6, 4), gg13);
    A [ 1 ] [ 1 ].add(D.at(2, 6), gg21);
    A [ 1 ] [ 1 ].add(D.at(2, 2), gg22);
    A [ 1 ] [ 1 ].add(D.at(2, 4), gg23);
    A [ 1 ] [ 1 ].add(D.at(4, 6), gg31);
    A [ 1 ] [ 1 ].add(D.at(4, 2), gg32);
    A [ 1 ] [ 1 ].add(D.at(4, 4), gg33);

    // position 32
    A [ 2 ] [ 1 ].resize(3, 3);
    A [ 2 ] [ 1 ].beUnitMatrix();
    A [ 2 ] [ 1 ].times( S.at(3, 2) );
    A [ 2 ] [ 1 ].add(D.at(5, 6), gg11);
    A [ 2 ] [ 1 ].add(D.at(5, 2), gg12);
    A [ 2 ] [ 1 ].add(D.at(5, 4), gg13);
    A [ 2 ] [ 1 ].add(D.at(4, 6), gg21);
    A [ 2 ] [ 1 ].add(D.at(4, 2), gg22);
    A [ 2 ] [ 1 ].add(D.at(4, 4), gg23);
    A [ 2 ] [ 1 ].add(D.at(3, 6), gg31);
    A [ 2 ] [ 1 ].add(D.at(3, 2), gg32);
    A [ 2 ] [ 1 ].add(D.at(3, 4), gg33);

    // position 13
    A [ 0 ] [ 2 ].resize(3, 3);
    A [ 0 ] [ 2 ].beUnitMatrix();
    A [ 0 ] [ 2 ].times( S.at(1, 3) );
    A [ 0 ] [ 2 ].add(D.at(1, 5), gg11);
    A [ 0 ] [ 2 ].add(D.at(1, 4), gg12);
    A [ 0 ] [ 2 ].add(D.at(1, 3), gg13);
    A [ 0 ] [ 2 ].add(D.at(6, 5), gg21);
    A [ 0 ] [ 2 ].add(D.at(6, 4), gg22);
    A [ 0 ] [ 2 ].add(D.at(6, 3), gg23);
    A [ 0 ] [ 2 ].add(D.at(5, 5), gg31);
    A [ 0 ] [ 2 ].add(D.at(5, 4), gg32);
    A [ 0 ] [ 2 ].add(D.at(5, 3), gg33);

    // position 23
    A [ 1 ] [ 2 ].resize(3, 3);
    A [ 1 ] [ 2 ].beUnitMatrix();
    A [ 1 ] [ 2 ].times( S.at(2, 3) );
    A [ 1 ] [ 2 ].add(D.at(6, 5), gg11);
    A [ 1 ] [ 2 ].add(D.at(6, 4), gg12);
    A [ 1 ] [ 2 ].add(D.at(6, 3), gg13);
    A [ 1 ] [ 2 ].add(D.at(2, 5), gg21);
    A [ 1 ] [ 2 ].add(D.at(2, 4), gg22);
    A [ 1 ] [ 2 ].add(D.at(2, 3), gg23);
    A [ 1 ] [ 2 ].add(D.at(4, 5), gg31);
    A [ 1 ] [ 2 ].add(D.at(4, 4), gg32);
    A [ 1 ] [ 2 ].add(D.at(4, 3), gg33);

    // position 33
    A [ 2 ] [ 2 ].resize(3, 3);
    A [ 2 ] [ 2 ].beUnitMatrix();
    A [ 2 ] [ 2 ].times( S.at(3, 3) );
    A [ 2 ] [ 2 ].add(D.at(5, 5), gg11);
    A [ 2 ] [ 2 ].add(D.at(5, 4), gg12);
    A [ 2 ] [ 2 ].add(D.at(5, 3), gg13);
    A [ 2 ] [ 2 ].add(D.at(4, 5), gg21);
    A [ 2 ] [ 2 ].add(D.at(4, 4), gg22);
    A [ 2 ] [ 2 ].add(D.at(4, 3), gg23);
    A [ 2 ] [ 2 ].add(D.at(3, 5), gg31);
    A [ 2 ] [ 2 ].add(D.at(3, 4), gg32);
    A [ 2 ] [ 2 ].add(D.at(3, 3), gg33);
}

void
Shell7Base :: computePressureTangentMatrix(FloatMatrix &answer, Load *load, const int iSurf, TimeStep *tStep)
{
    // Computes tangent matrix associated with the linearization of pressure loading. Assumes constant pressure.
    GaussPoint *gp;
    IntegrationRule *iRule = integrationRulesArray [ 1 ];   // rule #2 for mid-plane integration only

    double dA, a, b;
    FloatMatrix N, B, Nt, NtL, NtLB, L(7, 18);
    FloatArray lcoords, BF;
    double zeta = 1.e30;
    if ( iSurf == 1 ) {
        zeta = this->giveCrossSection()->give(CS_Thickness) * ( -0.5 );     // bottom surface -> iSurf = 1
    } else if ( iSurf == 2 ) {
        zeta = 0.0;                                                 // midplane surface -> iSurf = 2
    } else if ( iSurf == 3 ) {
        zeta = this->giveCrossSection()->give(CS_Thickness) * 0.5;       // top surface -> iSurf = 3
    } else {
        _error("computePressureForce: incompatible load surface must be 1, 2 or 3");
    }

    FloatArray solVec, pressure;
    this->giveUpdatedSolutionVector(solVec, tStep);


    int ndof = this->giveNumberOfDofs();
    answer.resize(ndof, ndof);
    answer.zero();
    for ( int i = 0; i < iRule->getNumberOfIntegrationPoints(); i++ ) {
        gp = iRule->getIntegrationPoint(i);
        lcoords = * gp->giveCoordinates();

        FloatArray g1, g2, g3, m(3), dm1(3), dm2(3), t1, t2;
        double gam, dg1, dg2;

        this->computeNmatrixAt(gp, N);
        this->computeBmatrixAt(gp, B);

        FloatArray genEps, temp1;
        genEps.beProductOf(B, solVec);        // [dxdxi, dmdxi, m, dgamdxi, gam]^T
        this->giveGeneralizedStrainComponents(genEps, temp1, temp1, dm1, dm2, m, dg1, dg2, gam);
        this->evalCovarBaseVectorsAt(gp, g1, g2, g3, genEps);
        FloatArray n;
        n.beVectorProductOf(g1, g2);                      // Compute normal (should not be normalized)

        // thickness coefficients
        a = zeta + 0.5 * gam * zeta * zeta;
        b = 0.5 * zeta * zeta;

        FloatMatrix W1(3, 3), W2(3, 3), temp3(3, 3);

        // W = skew-symmetric matrix
        W2.resize(3, 3);
        W2.zero();
        W2.at(2, 3) = -g2.at(1);
        W2.at(3, 2) =  g2.at(1);
        W2.at(1, 3) =  g2.at(2);
        W2.at(3, 1) = -g2.at(2);
        W2.at(1, 2) = -g2.at(3);
        W2.at(2, 1) =  g2.at(3);
        W2.negated();

        W1.resize(3, 3);
        W1.zero();
        W1.at(2, 3) = -g1.at(1);
        W1.at(3, 2) =  g1.at(1);
        W1.at(1, 3) =  g1.at(2);
        W1.at(3, 1) = -g1.at(2);
        W1.at(1, 2) = -g1.at(3);
        W1.at(2, 1) =  g1.at(3);
        W1.negated();



        // Tangent stiffness

        /*
         * [    W2     - W1      a*W2      -a*W1       b(W2*dg1-W1*dg2)          b*W2*m1      -b*W1*m2,      b(W2*dm1 - W1*dm2)       (3x18)
         * a*W2    -a*W1    a^2*W2      -a*W1      ab(W2*dg1-W1*dg2)        a*b*W2*m1    -a*b*W1*m2,     ab(W2*dm1 - W1*dm2)+bn    (3x18)
         * m^t*W2  -m^t*W1  a*m^t*W2  -a*m^t*W1   m^t*b(W2*dg1-W1*dg2)+n^t  b*m^t*W2*m1  -b*m^t*W1*m2  b*m^t*(W2*dm1 - W1*dm2)    ]  (1x18)
         */
        // L(1,1) = [W2, -W1]
        FloatMatrix L11(3, 6);
        L11.zero();
        temp3.add(-1, W1);
        L11.setSubMatrix(W2, 1, 1);
        L11.setSubMatrix(temp3, 1, 4);

        // L(1,2) = a*[W2, -W1] = a*L11
        FloatMatrix L12(3, 6);
        L12.zero();
        L12.add(a, L11);

        // L(1,3) =  b(W2*dg1-W1*dg2)
        FloatMatrix L13(3, 3);
        L13.zero();
        L13.add(dg1, W2);
        L13.add(-dg2, W1);
        L13.times(b);

        // L(1,4) = [b*W2*m, -b*W1*m]
        FloatMatrix L14(3, 2);
        L14.zero();
        t1.beProductOf(W2, m);
        t2.beProductOf(W1, m);
        t2.negated();
        L14.addSubVectorCol(t1, 1, 1);
        L14.addSubVectorCol(t2, 1, 2);
        L14.times(b);

        // L(1,5) = b(W2*dm1 - W1*dm2)
        FloatArray L15(3);
        L15.zero();
        t1.beProductOf(W2, dm1);
        t2.beProductOf(W1, dm2);
        L15.add(b, t1);
        L15.add(-b, t2);                       //L15.times(b);

        FloatMatrix L1(3, 18);
        L1.zero();
        L1.setSubMatrix(L11, 1, 1);
        L1.setSubMatrix(L12, 1, 7);
        L1.setSubMatrix(L13, 1, 13);
        L1.setSubMatrix(L14, 1, 16);
        L1.addSubVectorCol(L15, 1, 18);


        FloatMatrix L2(3, 18);
        L2.zero();
        L2.add(a, L1);
        t1.zero();
        t1.add(b, n);
        L2.addSubVectorCol(t1, 1, 18);


        // L(3,:) = m^t*L(:,1) + n^t*deltam
        FloatArray L3;
        L3.beTProductOf(L1, m);
        L3.times(b);
        L3.at(13) += b * n.at(1);
        L3.at(14) += b * n.at(2);
        L3.at(15) += b * n.at(3);

        L.zero();
        L.setSubMatrix(L1, 1, 1);
        L.setSubMatrix(L2, 4, 1);
        L.addSubVectorRow(L3, 7, 1);

        load->computeValueAt(pressure, tStep, * ( gp->giveCoordinates() ), VM_Total);        // pressure component
        L.times( -pressure.at(1) );


        // Tangent matrix (K = N^T*L*B*dA)
        NtL.beTProductOf(N, L);
        NtLB.beProductOf(NtL, B);
        dA = this->computeAreaAround(gp);
        NtLB.times(dA);
        answer.add(NtLB);
    }
}

#endif



// Strain and stress

#if 1

void
Shell7Base :: computeFAt(GaussPoint *gp, FloatMatrix &answer, FloatArray &genEps)
{
    // Compute deformation gradient as open product(g_i, G_i)
    FloatArray gcov1, gcov2, gcov3, Gcon1, Gcon2, Gcon3;
    this->evalCovarBaseVectorsAt(gp, gcov1, gcov2, gcov3, genEps);
    this->evalInitialContravarBaseVectorsAt(gp, Gcon1, Gcon2, Gcon3);

    FloatMatrix F1, F2, F3;
    F1.beDyadicProductOf(gcov1, Gcon1);
    F2.beDyadicProductOf(gcov2, Gcon2);
    F3.beDyadicProductOf(gcov3, Gcon3);
    answer.resize(3, 3);
    answer.zero();
    answer.add(F1);
    answer.add(F2);
    answer.add(F3);
}

void
Shell7Base :: computeStrainVector(FloatArray &answer, GaussPoint *gp, TimeStep *stepN, FloatArray &genEps)
{
    // Computes the Green-Lagrange strain tensor: E=0.5(C-I)
    FloatMatrix F, E;
    this->computeFAt(gp, F, genEps);     // Deformation gradient

    E.beTProductOf(F, F);    // C-Right Caucy-Green deformation tensor
    E.at(1, 1) += -1;
    E.at(2, 2) += -1;
    E.at(3, 3) += -1;
    E.times(0.5);

    FloatArray temp(6);
    temp.beReducedVectorForm(E);     // Convert to Voight form Todo: add enum strain/stress
    answer = temp;
    answer.at(4) = temp.at(4) * 2.0;   // correction of shear strains
    answer.at(5) = temp.at(5) * 2.0;
    answer.at(6) = temp.at(6) * 2.0;
}

void
Shell7Base :: computeStressVector(FloatArray &answer, FloatArray &genEps, GaussPoint *gp, Material *mat, TimeStep *stepN)
{
    FloatArray E;
    this->computeStrainVector(E, gp, stepN, genEps);     // Green-Lagrange strain vector
    static_cast< StructuralMaterial * >( mat )->giveRealStressVector(answer, ReducedForm, gp, E, stepN);
}

void
Shell7Base :: computeStressResultantsAt(GaussPoint *gp, FloatArray &Svec, FloatArray &S1g, FloatArray &S2g, FloatArray &S3g, FloatArray &genEps)
{
    FloatArray g1, g2, g3;
    this->evalCovarBaseVectorsAt(gp, g1, g2, g3, genEps);
    FloatMatrix S;
    S.beMatrixForm(Svec);

    // Sig =S(i,j)*g_j, - stress vectors on the surfaces given by g_j?
    for ( int j = 1; j <= 3; j++ ) {
        S1g.at(j) = S.at(1, 1) * g1.at(j) + S.at(1, 2) * g2.at(j) + S.at(1, 3) * g3.at(j);
        S2g.at(j) = S.at(2, 1) * g1.at(j) + S.at(2, 2) * g2.at(j) + S.at(2, 3) * g3.at(j);
        S3g.at(j) = S.at(3, 1) * g1.at(j) + S.at(3, 2) * g2.at(j) + S.at(3, 3) * g3.at(j);
    }
}

#endif




// Internal forces

#if 1


void
Shell7Base :: giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord)
//
// Computes internal forces as a summation of: sectional forces + convective mass force
{
    FloatArray totalSolVec, solVec;
    //this->giveUpdatedSolutionVector(solVec, tStep);
    this->giveUpdatedSolutionVector(solVec, tStep); // da
    answer.resize( this->giveNumberOfDofs() );
    answer.zero();

    FloatArray temp;
    this->computeSectionalForces(answer, tStep, solVec, useUpdatedGpRecord);

    ///@todo How to treat the convective force? Only active during dynamic simulations
}


void
Shell7Base :: computeSectionalForces(FloatArray &answer, TimeStep *tStep, FloatArray &solVec, int useUpdatedGpRecord)
//
{
    FloatMatrix B;
    FloatArray BtF, f, genEps;
    answer.resize( Shell7Base :: giveNumberOfDofs() );
    answer.zero();

    LayeredCrossSection *layeredCS = dynamic_cast< LayeredCrossSection * >( this->giveCrossSection() );
    int numberOfLayers = layeredCS->giveNumberOfLayers();     // conversion of types
    FloatArray f1(18), f2(18), f3(6);
    f1.zero();
    f2.zero();
    f3.zero();
    for ( int layer = 1; layer <= numberOfLayers; layer++ ) {
        IntegrationRule *iRuleL = layerIntegrationRulesArray [ layer - 1 ];
        Material *mat = domain->giveMaterial( layeredCS->giveLayerMaterial(layer) );

        for ( int j = 1; j <= iRuleL->getNumberOfIntegrationPoints(); j++ ) {
            GaussPoint *gp = iRuleL->getIntegrationPoint(j - 1);
            FloatMatrix B11, B22, B32, B43, B53;
            this->computeBmatricesAt(gp, B11, B22, B32, B43, B53);
            //this->computeGeneralizedStrainVector(genEps, solVec, B11, B22, B32, B43, B53);


            FloatArray genEpsD, totalSolVec;
            this->computeGeneralizedStrainVector(genEpsD, solVec, B11, B22, B32, B43, B53);

            this->giveUpdatedSolutionVector(totalSolVec, tStep); // a
            this->computeGeneralizedStrainVector(genEps, totalSolVec, B11, B22, B32, B43, B53);

            double zeta = giveGlobalZcoord(gp);
            FloatArray N, M, T, Ms;
            double Ts = 0.;
            this->computeSectionalForcesAt(N, M, T, Ms, Ts, gp, mat, tStep, genEps, genEpsD, zeta);

            // Computation of sectional forces: f = B^t*[N M T Ms Ts]^t
            FloatArray f1temp(18), f2temp(18), f3temp(6), temp;
            // f1 = BT11*N
            f1temp.beTProductOf(B11, N);

            // f2 = BT22*M + BT32*T
            f2temp.beTProductOf(B22, M);
            temp.beTProductOf(B32, T);
            f2temp.add(temp);

            // f3 = BT43*Ms + BT53*Ts
            f3temp.beTProductOf(B43, Ms);
            for ( int i = 1; i <= 6; i++ ) {
                f3temp.at(i) += B53.at(1, i) * Ts;
            }

            double dV = this->computeVolumeAroundLayer(gp, layer);
            f1.add(dV, f1temp);
            f2.add(dV, f2temp);
            f3.add(dV, f3temp);
        }
    }
    //f1.printYourself();
    IntArray ordering_phibar = giveOrdering(Midplane);
    IntArray ordering_m = giveOrdering(Director);
    IntArray ordering_gam = giveOrdering(InhomStrain);
    answer.assemble(f1, ordering_phibar);
    answer.assemble(f2, ordering_m);
    answer.assemble(f3, ordering_gam);
    
}


void
Shell7Base :: computeSectionalForcesAt(FloatArray &N, FloatArray  &M, FloatArray &T, FloatArray  &Ms, double &Ts, 
GaussPoint *gp, Material *mat, TimeStep *tStep, FloatArray &genEps, FloatArray &genEpsD, double zeta)
{

#if 1
    FloatArray S1g(3), S2g(3), S3g(3);
    FloatArray cartStressVector, contravarStressVector;
    FloatMatrix lambda[3];

    this->computeLambdaMatrices(lambda, genEpsD, zeta);


    this->computeStressVector(cartStressVector, genEps, gp, mat, tStep);
    this->transInitialCartesianToInitialContravar(gp, cartStressVector, contravarStressVector);
    this->computeStressResultantsAt(gp, contravarStressVector, S1g, S2g, S3g, genEps);

    
    FloatArray sectionalForces(18), temp;
    sectionalForces.zero();
    temp.beTProductOf(lambda [0], S1g);
    sectionalForces.add(temp);
    temp.beTProductOf(lambda [1], S2g);
    sectionalForces.add(temp);
    temp.beTProductOf(lambda [2], S3g);
    sectionalForces.add(temp);


    /* The following expressions are sectional forces (in the classic sense) if integrated over the thickness. However,
     * now they are integrated over the volume to produce f_int.*/

    
    N.resize(6);  // Membrane forces - N1, N2
    M.resize(6);  // Moments - M1, M2
    for ( int i = 1; i <=6; i++ ) {
        N.at(i) = sectionalForces.at(i);
        M.at(i) = sectionalForces.at(6+i);
    }
    
    // Shear force - T
    T.resize(3);
    T.at(1) = sectionalForces.at(13);
    T.at(2) = sectionalForces.at(14);
    T.at(3) = sectionalForces.at(15);

    // Higher order force - Ms
    Ms.resize(2);
    Ms.at(1) = sectionalForces.at(16);
    Ms.at(2) = sectionalForces.at(17);
    
    // Ts
    Ts = sectionalForces.at(18);

#else

    FloatArray g1, g2, g3, S1g(3), S2g(3), S3g(3), m(3), dm1(3), dm2(3), temp1, temp2;
    FloatArray lcoords, cartStressVector, contravarStressVector, sectionalForces, BF, a;
    double fac1, fac2, fac3, gam, dg1, dg2;

    lcoords = * gp->giveCoordinates();

    this->giveGeneralizedStrainComponents(genEps, temp1, temp1, dm1, dm2, m, dg1, dg2, gam);

    this->computeStressVector(cartStressVector, genEps, gp, mat, tStep);

    this->transInitialCartesianToInitialContravar(gp, cartStressVector, contravarStressVector);
    this->computeStressResultantsAt(gp, contravarStressVector, S1g, S2g, S3g, genEps);

    fac1 = zeta + 0.5 * gam * zeta * zeta;
    fac2 = 1. + gam * zeta;
    fac3 = 0.5 * zeta * zeta;

    /* The following expressions are sectional forces (in the classic sense) if integrated over the thickness. However,
     * now they are integrated over the volume to produce f_int.*/

    // Membrane forces - N1, N2
    N.resize(6);
    N.at(1) = S1g.at(1);
    N.at(2) = S1g.at(2);
    N.at(3) = S1g.at(3);
    N.at(4) = S2g.at(1);
    N.at(5) = S2g.at(2);
    N.at(6) = S2g.at(3);
    N.printYourself();
    // Moments - M1, M2
    M.resize(6);
    M.at(1) = fac1 * S1g.at(1);
    M.at(2) = fac1 * S1g.at(2);
    M.at(3) = fac1 * S1g.at(3);
    M.at(4) = fac1 * S2g.at(1);
    M.at(5) = fac1 * S2g.at(2);
    M.at(6) = fac1 * S2g.at(3);
    M.printYourself();
    // Shear force - T
    T.resize(3);
    T.at(1) = fac2 * S3g.at(1) + fac3 * ( S1g.at(1) * dg1 + S2g.at(1) * dg2 );
    T.at(2) = fac2 * S3g.at(2) + fac3 * ( S1g.at(2) * dg1 + S2g.at(2) * dg2 );
    T.at(3) = fac2 * S3g.at(3) + fac3 * ( S1g.at(3) * dg1 + S2g.at(3) * dg2 );


    // Higher order force - Ms
    Ms.resize(2);
    Ms.at(1) = fac3 * m.dotProduct(S1g);
    Ms.at(2) = fac3 * m.dotProduct(S2g);

    // Ts
    Ts = fac3 * ( dm1.dotProduct(S1g) + dm2.dotProduct(S2g) ) + zeta * m.dotProduct(S3g);
#endif
}

#endif




// Mass matrices

#if 1

void
Shell7Base :: computeThicknessMappingCoeff(GaussPoint *gp, FloatArray &answer)
{
    //thickness jacobian = ratio between volume and area: j0 = a3 + a2*zeta^2 + a1 * zeta
    // Returns array with a1-a3, used in expression for analytical integration of mass matrix.
    FloatArray lcoords = * gp->giveCoordinates();

    FloatMatrix dNdxi;
    FEInterpolation3d *fei = static_cast< FEInterpolation3d * >( this->giveInterpolation() );
    fei->evaldNdxi( dNdxi, lcoords, FEIElementGeometryWrapper(this) );

    FloatArray M, dM1(3), dM2(3), dX1(3), dX2(3);
    double gam, dg1, dg2;
    FloatMatrix B11, B22, B32, B43, B53;
    this->computeBmatricesAt(gp, B11, B22, B32, B43, B53);

    FloatArray initSolVec, genEps;
    this->giveInitialSolutionVector(initSolVec);

    this->computeGeneralizedStrainVector(genEps, initSolVec, B11, B22, B32, B43, B53);
    this->giveGeneralizedStrainComponents(genEps, dX1, dX2, dM1, dM2, M, dg1, dg2, gam);


    FloatArray temp, temp2;
    temp.beVectorProductOf(dX1, dX2);
    double sc = temp.computeNorm();
    answer.resize(3);
    answer.at(3) = M.dotProduct(temp) / sc;

    temp.beVectorProductOf(dX1, dM2);
    temp2.beVectorProductOf(dX2, dM1);
    temp.add(temp2);
    answer.at(2) = M.dotProduct(temp) / sc;

    temp.beVectorProductOf(dM1, dM2);
    //answer.at(1) = M.dotProduct(temp)/sc;
    answer.at(1) = temp.computeNorm() / sc;
}

void
Shell7Base :: computeLumpedMassMatrix(FloatMatrix &answer, TimeStep *tStep)
{
    // TODO: add algorithm for this
    OOFEM_ERROR("Shell7Base :: computeLumpedMassMatrix - No lumping algorithm implemented");
}


void
Shell7Base :: computeMassMatrix(FloatMatrix &answer, TimeStep *tStep)
{
    // Analytically integrated over the thickness. Constant density assumed.
    // => integration rule #2 = midplane only
    IntegrationRule *iRule = integrationRulesArray [ 1 ];

    //------------------------------
    FloatMatrix N, Nt, Ntm, NtmN, temp;
    FloatArray solVec;
    this->giveUpdatedSolutionVector(solVec, tStep);

    LayeredCrossSection *layeredCS = dynamic_cast< LayeredCrossSection * >( this->giveCrossSection() );
    Material *mat = domain->giveMaterial( layeredCS->giveLayerMaterial(1) );     // for now, while I don't have an analytical exp.


    for ( int i = 0; i < iRule->getNumberOfIntegrationPoints(); i++ ) {
        GaussPoint *gp = iRule->getIntegrationPoint(i);

        this->computeNmatrixAt(gp, N);
        FloatArray unknowns, m(3);
        unknowns.beProductOf(N, solVec);        // [x, m, gam]^T
        m.setValues( 3, unknowns.at(4), unknowns.at(5), unknowns.at(6) );
        double gam = unknowns.at(7);

        // Analytically integrated through the tickness
        /*     3    3    1
         * 3 [a*I  b*I   c*m      [A  B  C
         * 3       d*I   e*m    =     D  E
         * 1  sym       f*m.m]     sym   F]
         */

        double rho = mat->give('d', gp);
        FloatArray factors;
        this->giveMassFactorsAt(gp, factors, gam);

        FloatMatrix mass;
        mass.resize(7, 7);
        mass.at(1, 1) = mass.at(2, 2) = mass.at(3, 3) = factors.at(1);      // A
        mass.at(1, 4) = mass.at(2, 5) = mass.at(3, 6) = factors.at(2);      // B
        mass.at(1, 7) = factors.at(3) * m.at(1);
        mass.at(2, 7) = factors.at(3) * m.at(2);
        mass.at(3, 7) = factors.at(3) * m.at(3);        // C
        mass.at(4, 4) = mass.at(5, 5) = mass.at(6, 6) = factors.at(4);      // D
        mass.at(4, 7) = factors.at(5) * m.at(1);
        mass.at(5, 7) = factors.at(5) * m.at(2);
        mass.at(6, 7) = factors.at(5) * m.at(3);        // E
        mass.at(7, 7) = factors.at(6) * m.dotProduct(m);        // F
        mass.symmetrized();

        double dA = this->computeAreaAround(gp);
        Ntm.beTProductOf(N, mass);
        NtmN.beProductOf(Ntm, N);         // switch to sym prod upper something
        NtmN.times(dA * rho);
        temp.add(NtmN);
    }

    int ndofs = this->computeNumberOfDofs(EID_MomentumBalance);
    answer.resize(ndofs, ndofs);
    answer.zero();
    IntArray ordering_all = giveOrdering(All);
    answer.assemble(temp, ordering_all);
    answer.symmetrized();

 #if 1
    FloatMatrix A(18, 18), B(18, 18), C(18, 6), D(18, 18), E(18, 6), F(6, 6), eigVec, temp2;

    temp2.beSubMatrixOf(temp, 1, 18, 19, 36);
    //   temp2.printYourself();
    A.beSubMatrixOf(temp, 1, 18, 1, 18);
    B.beSubMatrixOf(temp, 1, 18, 19, 36);
    C.beSubMatrixOf(temp, 1, 18, 37, 42);
    D.beSubMatrixOf(temp, 19, 36, 19, 36);
    E.beSubMatrixOf(temp, 19, 36, 37, 42);
    F.beSubMatrixOf(temp, 37, 42, 37, 42);
    //    B.printYourself();

    FloatArray eig, sums(6);
    sums.zero();
    for ( int i = 1; i <= 18; i++ ) {
        sums.at(1) += A.at(i, i);
        sums.at(2) += B.at(i, i);
        sums.at(3) += C.at(i, 1) + C.at(i, 2) + C.at(i, 3) + C.at(i, 4) + C.at(i, 5) + C.at(i, 6);
        sums.at(4) += D.at(i, i);
        sums.at(5) += E.at(i, 1) + E.at(i, 2) + E.at(i, 3) + E.at(i, 4) + E.at(i, 5) + E.at(i, 6);
    }

    sums.at(6) += F.at(1, 1) + F.at(2, 2) + F.at(3, 3) + F.at(4, 4) + F.at(5, 5) + F.at(6, 6);
    printf("\n analytical \n");
    sums.printYourself();
 #endif
}


void
Shell7Base :: giveMassFactorsAt(GaussPoint *gp, FloatArray &factors, double &gam)
{
    double a1, a2, a3;
    FloatArray coeff;
    this->computeThicknessMappingCoeff(gp, coeff);
    a1 = coeff.at(1);
    a2 = coeff.at(2);
    a3 = coeff.at(3);

    double h, h2, h3, h5, gam2;
    h   = this->giveCrossSection()->give(CS_Thickness);
    h2  = h * h;
    h3 = h2 * h;
    h5 = h2 * h3;
    gam2 = gam * gam;

    factors.resize(6);
    factors.at(1) = a3 * h + ( a1 * h3 ) / 12.;
    factors.at(2) = h3 * ( 40. * a2 + 20. * a3 * gam + 3. * a1 * h2 * gam ) / 480.;
    factors.at(3) = h3 * ( 20. * a3 + 3. * a1 * h2 ) / 480. * 1.0;
    factors.at(4) = ( 28. * a3 * h3 * ( 80. + 3. * h2 * gam2 ) + 3. * h5 * ( 112. * a2 * gam + a1 * ( 112. + 5. * h2 * gam2 ) ) ) / 26880.;
    factors.at(5) = h5 * ( 56. * a2 + 28. * a3 * gam + 5. * a1 * h2 * gam ) / 8960.;
    factors.at(6) = h5 * ( 28. * a3 + 5. * a1 * h2 ) / 8960.;
}

void
Shell7Base :: computeMassMatrixNum(FloatMatrix &answer, TimeStep *tStep)
{
    // Num refers in this case to  numerical integration in both in-plane and through the thickness.
    // For analytically integrated throught he thickness, see computeMassMatrix


    FloatMatrix N, Nt, Ntm, NtmN, mass, temp;
    FloatArray solVec, unknowns;
    this->giveUpdatedSolutionVector(solVec, tStep);
    int ndofs = this->giveNumberOfDofs();
    temp.resize(ndofs, ndofs);
    temp.zero();


    LayeredCrossSection *layeredCS = dynamic_cast< LayeredCrossSection * >( this->giveCrossSection() );
    int numberOfLayers = layeredCS->giveNumberOfLayers();     // conversion of data

    FloatMatrix M11(18, 18), M12(18, 18), M13(18, 6), M22(18, 18), M23(18, 6), M33(6, 6);
    M11.zero();
    M12.zero();
    M13.zero();
    M22.zero();
    M23.zero();
    M33.zero();

    for ( int layer = 1; layer <= numberOfLayers; layer++ ) {
        IntegrationRule *iRuleL = layerIntegrationRulesArray [ layer - 1 ];
        Material *mat = domain->giveMaterial( layeredCS->giveLayerMaterial(layer) );

        for ( int j = 1; j <= iRuleL->getNumberOfIntegrationPoints(); j++ ) {
            GaussPoint *gp = iRuleL->getIntegrationPoint(j - 1);

            FloatMatrix N11, N22, N33;
            this->computeNmatricesAt(gp, N11, N22, N33);
            FloatArray xbar, m;
            double gam = 0.;
            this->computeSolutionFields(xbar, m, gam, solVec, N11, N22, N33);
            //this->computeNmatrixAt(gp, N);
            //unknowns.beProductOf(N,a); // [xbar, m, gam]^T
            //m.setValues(3, unknowns.at(4), unknowns.at(5), unknowns.at(6) );
            //double gam = unknowns.at(7);


            /* Consistent mass matrix M = int{N^t*mass*N}
             *
             *         3    3    1
             *         3 [a*I  b*I   c*m      [A  B  C
             * mass =   3       d*I   e*m    =     D  E
             *         1  sym       f*m.m]     sym   F]
             */


            double zeta = giveGlobalZcoord(gp);
            double fac1 = 4;
            double fac2 = 2.0 * zeta * ( 2.0 + gam * zeta );
            double fac3 = 2.0 * zeta * zeta;
            double fac4 = zeta * zeta * ( 2.0 + gam * zeta ) * ( 2.0 + gam * zeta );
            double fac5 = zeta * zeta * zeta * ( 2.0 + gam * zeta );
            double fac6 = zeta * zeta * zeta * zeta;
            FloatMatrix mass11(3, 3), mass12(3, 3), mass13(3, 1), mass21(3, 3), mass22(3, 3), mass23(3, 1), mass31(1, 3), mass32(1, 3), mass33(1, 1);
            mass11.zero();
            mass12.zero();
            mass13.zero();
            mass21.zero();
            mass22.zero();
            mass23.zero();
            mass31.zero();
            mass32.zero();
            mass33.zero();
            mass.resize(7, 7);
            mass11.at(1, 1) = mass11.at(2, 2) = mass11.at(3, 3) = fac1;          // A
            mass12.at(1, 1) = mass12.at(2, 2) = mass12.at(3, 3) = fac2;          // B
            mass13.at(1, 1) = fac3 * m.at(1);
            mass13.at(2, 1) = fac3 * m.at(2);
            mass13.at(3, 1) = fac3 * m.at(3);            // C
            mass22.at(1, 1) = mass22.at(2, 2) = mass22.at(3, 3) = fac4;          // D
            mass23.at(1, 1) = fac5 * m.at(1);
            mass23.at(2, 1) = fac5 * m.at(2);
            mass23.at(3, 1) = fac5 * m.at(3);            // E
            mass33.at(1, 1) = fac6 * m.dotProduct(m);            // F
            mass21.beTranspositionOf(mass12);
            mass31.beTranspositionOf(mass13);
            mass32.beTranspositionOf(mass23);
            //mass.symmetrized();

            double dV = this->computeVolumeAroundLayer(gp, layer);
            double rho = mat->give('d', gp);

            FloatMatrix M11temp, M12temp, M13temp, M22temp, M23temp, M33temp;
            this->computeTripleProduct(M11temp, N11, mass11, N11);
            this->computeTripleProduct(M12temp, N11, mass12, N22);
            this->computeTripleProduct(M13temp, N11, mass13, N33);
            this->computeTripleProduct(M22temp, N22, mass22, N22);
            this->computeTripleProduct(M23temp, N22, mass23, N33);
            this->computeTripleProduct(M33temp, N33, mass33, N33);
            M11.add(0.25 * rho * dV, M11temp);
            M12.add(0.25 * rho * dV, M12temp);
            M13.add(0.25 * rho * dV, M13temp);
            M22.add(0.25 * rho * dV, M22temp);
            M23.add(0.25 * rho * dV, M23temp);
            M33.add(0.25 * rho * dV, M33temp);
        }

        //M33.printYourself();
        answer.resize(ndofs, ndofs);
        answer.zero();

        IntArray ordering_phibar = giveOrdering(Midplane);
        IntArray ordering_m = giveOrdering(Director);
        IntArray ordering_gam = giveOrdering(InhomStrain);
        answer.assemble(M11, ordering_phibar, ordering_phibar);
        answer.assemble(M12, ordering_phibar, ordering_m);
        answer.assemble(M13, ordering_phibar, ordering_gam);
        answer.assemble(M22, ordering_m,      ordering_m);
        answer.assemble(M23, ordering_m,      ordering_gam);
        answer.assemble(M33, ordering_gam,    ordering_gam);

        FloatMatrix M21, M31, M32;
        M21.beTranspositionOf(M12);
        M31.beTranspositionOf(M13);
        M32.beTranspositionOf(M23);
        answer.assemble(M21, ordering_m,      ordering_phibar);
        answer.assemble(M31, ordering_gam,    ordering_phibar);
        answer.assemble(M32, ordering_gam,    ordering_m);
        answer.symmetrized();

 #if 1

        FloatArray eig, sums(6);
        sums.zero();
        for ( int i = 1; i <= 18; i++ ) {
            sums.at(1) += M11.at(i, i);
            sums.at(2) += M12.at(i, i);
            sums.at(3) += M13.at(i, 1) + M13.at(i, 2) + M13.at(i, 3) + M13.at(i, 4) + M13.at(i, 5) + M13.at(i, 6);
            sums.at(4) += M22.at(i, i);
            sums.at(5) += M23.at(i, 1) + M23.at(i, 2) + M23.at(i, 3) + M23.at(i, 4) + M23.at(i, 5) + M23.at(i, 6);
        }

        sums.at(6) += M33.at(1, 1) + M33.at(2, 2) + M33.at(3, 3) + M33.at(4, 4) + M33.at(5, 5) + M33.at(6, 6);
        printf("\n numerical \n");
        sums.printYourself();



 #endif
    }
}


void
Shell7Base :: computeConvectiveMassForce(FloatArray &answer, TimeStep *tStep)
{
    // Analytically integrated over the thickness. Constant density assumed.

    IntegrationRule *iRule = integrationRulesArray [ 1 ];   // rule 2 for mid-plane integration only
    GaussPoint *gp;


    //Material *mat = this->giveMaterial();
    FloatMatrix N, Nt;
    FloatArray lcoords, cartStressVector, contravarStressVector, sectionalForces, BF, a, da, unknowns, m, dm, aVec, daVec, fm(7), fM;
    double gam, dgam, dA;

    answer.resize(42);
    answer.zero();
    for ( int i = 0; i < iRule->getNumberOfIntegrationPoints(); i++ ) {
        gp = iRule->getIntegrationPoint(i);
        this->computeNmatrixAt(gp, N);
        this->giveUpdatedSolutionVector(aVec, tStep);
        // Fix for the new numbering in B & N
        this->computeVectorOf(EID_MomentumBalance, VM_Velocity, tStep, daVec);

        a.beProductOf(N, aVec);        // [ x,  m,  gam]^T
        da.beProductOf(N, daVec);        // [dx, dm, dgam]^T
        m.setValues( 3,  a.at(4),  a.at(5),  a.at(6) );
        gam =  a.at(7);
        dm.setValues( 3, da.at(4), da.at(5), da.at(6) );
        dgam = da.at(7);

        double a1, a2, a3, h, h2, h3, h5, fac1, fac2, fac3, rho;
        FloatArray coeff;


        rho = this->giveMaterial()->give('d', gp);
        h   = this->giveCrossSection()->give(CS_Thickness);
        h2  = h * h;
        h3 = h2 * h;
        h5 = h2 * h3;
        this->computeThicknessMappingCoeff(gp, coeff);
        a1 = coeff.at(1);
        a2 = coeff.at(2);
        a3 = coeff.at(3);

        // Convective mass "matrix"
        /*     3
         * 3 [ a*m
         * 3   b*m
         * 1  c*m.dm]*dgam
         */

        fac1 = rho * h3 * ( 20. * a3 + 3. * a1 * h2 ) / 240.;
        fac2 = rho * h5 * ( 56. * a2 + 28. * a3 * gam + 5. * a1 * h2 * gam ) / 4480.;
        fac3 = rho * h5 * ( 28. * a3 + 5. * a1 * h2 ) / 4480.;
        fm.at(1) = fac1 * dm.at(1) * dgam;
        fm.at(2) = fac1 * dm.at(2) * dgam;
        fm.at(3) = fac1 * dm.at(3) * dgam;
        fm.at(4) = fac2 * dm.at(1) * dgam;
        fm.at(5) = fac2 * dm.at(2) * dgam;
        fm.at(6) = fac2 * dm.at(3) * dgam;
        fm.at(7) = fac3 * m.dotProduct(dm) * gam;

        dA = this->computeAreaAround(gp);
        fM.beTProductOf(N, fm);
        fM.times(dA);
        answer.add(fM);
    }
}


void
Shell7Base :: computeTripleProduct(FloatMatrix &answer, const FloatMatrix &a, const FloatMatrix &b, const FloatMatrix &c)
{
    // Computes the product a^T*b*c
    FloatMatrix temp;
    temp.beTProductOf(a, b);
    answer.beProductOf(temp, c);
}


#endif // End mass matrices




// External forces

#if 1
void
Shell7Base :: computeEdgeLoadVectorAt(FloatArray &answer, Load *load, int iEdge, TimeStep *tStep, ValueModeType mode)
{
    BoundaryLoad *edgeLoad = dynamic_cast< BoundaryLoad * >( load );


    if ( edgeLoad ) {
        FloatArray fT, components, traction(3), temp;
        this->computeTractionForce(fT, iEdge, edgeLoad, tStep);
        IntArray mask;
        this->giveEdgeDofMapping(mask, iEdge);
        answer.resize( this->computeNumberOfDofs(EID_MomentumBalance) );
        answer.zero();
        answer.assemble(fT, mask);

        return;
    } else {
        _error("computeEdgeLoadVectorAt: incompatible load");
        return;
    }
}


// Surface
void
Shell7Base :: computeSurfaceLoadVectorAt(FloatArray &answer, Load *load,
                                         int iSurf, TimeStep *tStep, ValueModeType mode)
{
    BoundaryLoad *surfLoad = dynamic_cast< BoundaryLoad * >( load );
    if ( surfLoad ) {
        FloatArray solVec, force;
        this->giveUpdatedSolutionVector(solVec, tStep);
        this->computePressureForce(force, solVec, iSurf, surfLoad, tStep);

        IntArray mask;
        this->giveSurfaceDofMapping(mask, 1);         // same dofs regardless of iSurf
        answer.resize( this->computeNumberOfDofs(EID_MomentumBalance) );
        answer.zero();
        IntArray ordering_all = giveOrdering(All);
        answer.assemble(force, ordering_all);

        return;
    } else {
        _error("computeSurfaceLoadVectorAt: incompatible load");
        return;
    }
}

void
Shell7Base :: computePressureForce(FloatArray &answer, FloatArray solVec, const int iSurf, BoundaryLoad *surfLoad, TimeStep *tStep)
{
    // Computes pressure loading. Acts normal to the current (deformed) surface.

    // Should be special integration rule for top and bottom surface!!
    IntegrationRule *iRule = integrationRulesArray [ 1 ];   // rule #2 for mid-plane integration only
    GaussPoint *gp;

    FloatMatrix N, B;
    FloatArray Fp, fp, temp, n(3), m, aVec, a, load, traction, genEps;
    ;

    answer.resize( this->giveNumberOfDofs() );
    answer.zero();
    for ( int i = 0; i < iRule->getNumberOfIntegrationPoints(); i++ ) {
        gp = iRule->getIntegrationPoint(i);
        this->computeBmatrixAt(gp, B);
        this->computeNmatrixAt(gp, N);
        genEps.beProductOf(B, solVec);        // [dxdxi, dmdxi, m, dgamdxi, gam]^T

        this->computePressureForceAt(gp, fp, iSurf, genEps, surfLoad, tStep);

        double dA = this->computeAreaAround(gp);
        Fp.beTProductOf(N, fp);
        Fp.times(dA);
        answer.add(Fp);
    }
}


void
Shell7Base :: computePressureForceAt(GaussPoint *gp, FloatArray &answer, const int iSurf, FloatArray genEps, BoundaryLoad *surfLoad, TimeStep *tStep)
{
    // Computes pressure loading. Acts normal to the current (deformed) surface.

    FloatArray g1, g2, g3, m, load, traction;
    double gam;
    answer.resize(7);
    answer.zero();
    double zeta = 1.e30;
    if ( iSurf == 1 ) {
        zeta = this->giveCrossSection()->give(CS_Thickness) * ( -0.5 );     // bottom surface -> iSurf = 1
    } else if ( iSurf == 2 ) {
        zeta = 0.0;                                                 // midplane surface -> iSurf = 2
    } else if ( iSurf == 3 ) {
        zeta = this->giveCrossSection()->give(CS_Thickness) * 0.5;       // top surface -> iSurf = 3
    } else {
        _error("computePressureForceAt: incompatible load surface must be 1, 2 or 3");
    }

    m.setValues( 3,  genEps.at(13),  genEps.at(14),  genEps.at(15) );
    gam =  genEps.at(18);

    if ( dynamic_cast< ConstantPressureLoad * >( surfLoad ) ) {
        this->evalCovarBaseVectorsAt(gp, g1, g2, g3, genEps);         // m=g3
        surfLoad->computeValueAt(load, tStep, * ( gp->giveCoordinates() ), VM_Total);        // pressure components
        traction.beVectorProductOf(g1, g2);        // normal vector (unnormalized)
        traction.times( -load.at(1) );
    } else if ( dynamic_cast< ConstantSurfaceLoad * >( surfLoad ) ) {
        surfLoad->computeValueAt(traction, tStep, * ( gp->giveCoordinates() ), VM_Total);        // traction vector
    } else {
        _error("computePressureForceAt: incompatible load type");
    }


    double fac1 = zeta + 0.5 * zeta * zeta * gam;
    double fac2 = 0.5 * zeta * zeta;
    /*  Force
     * fp = [      n
     *         f1*n
     *    f2*dotprod(n,m)]
     */
    answer.at(1) = traction.at(1);
    answer.at(2) = traction.at(2);
    answer.at(3) = traction.at(3);
    answer.at(4) = fac1 * traction.at(1);
    answer.at(5) = fac1 * traction.at(2);
    answer.at(6) = fac1 * traction.at(3);
    answer.at(7) = fac2 * m.dotProduct(traction);
}


void
Shell7Base :: computeTractionForce(FloatArray &answer, const int iedge, BoundaryLoad *edgeLoad, TimeStep *tStep)
{
    // fix such that one can specify if the load should follow the deformed coord sys
    IntegrationRule *iRule = integrationRulesArray [ 2 ];   // rule #3 for edge integration of distributed loads given in [*/m]
    GaussPoint *gp;

    FloatMatrix N, Q;
    FloatArray g1, g2, g3, FT, fT(7), m, aVec, a, T(3), M, G1, components, lcoords;
    double dA;
    answer.resize( this->giveNumberOfEdgeDofs() );
    answer.zero();

    for ( int i = 0; i < iRule->getNumberOfIntegrationPoints(); i++ ) {
        gp = iRule->getIntegrationPoint(i);
        lcoords = * gp->giveCoordinates();

        edgeLoad->computeValueAt(components, tStep, lcoords, VM_Total);
        this->edgeComputeNmatrixAt(gp, N);
        this->edgeEvalCovarBaseVectorsAt(gp, iedge, g2, g3, tStep);
        g2.normalize();
        g3.normalize();
        g1.beVectorProductOf(g2, g3);
        g1.normalize();
        this->giveCoordTransMatrix(Q, g1, g2, g3);

        FloatArray c1(3), c2(3), t1, t2;
        c1.at(1) = components.at(1);
        c1.at(2) = components.at(2);
        c1.at(3) = components.at(3);
        c2.at(1) = components.at(4);
        c2.at(2) = components.at(5);
        c2.at(3) = components.at(6);
        t1.beTProductOf(Q, c1);
        t2.beTProductOf(Q, c2);

        fT.at(1) = t1.at(1);
        fT.at(2) = t1.at(2);
        fT.at(3) = t1.at(3);
        fT.at(4) = t2.at(1);
        fT.at(5) = t2.at(2);
        fT.at(6) = t2.at(3);
        fT.at(7) = components.at(7);
        fT = components;
        dA = this->edgeComputeLengthAround(gp, iedge);
        FT.beTProductOf(N, fT);
        FT.times(dA);
        answer.add(FT);
    }
}


void
Shell7Base :: computeBodyLoadVectorAt(FloatArray &answer, Load *forLoad, TimeStep *stepN, ValueModeType mode)
{
    OOFEM_ERROR("Error in computeBodyLoadVectorAt - currently not implemented");
}

#endif



// Integration
#if 1

double
Shell7Base :: computeAreaAround(GaussPoint *gp)
{
    FloatArray G1, G2, G3, temp;
    double detJ;
    this->evalInitialCovarBaseVectorsAt(gp, G1, G2, G3);
    temp.beVectorProductOf(G1, G2);
    detJ = temp.computeNorm();
    return detJ * gp->giveWeight() * 0.5;
}

double
Shell7Base :: edgeComputeLengthAround(GaussPoint *gp, const int iedge)
{
    FloatArray G1, G3, temp;
    double detJ;
    this->edgeEvalInitialCovarBaseVectorsAt(gp, iedge, G1, G3);
    detJ = G1.computeNorm();
    return detJ * gp->giveWeight();
}


double
Shell7Base :: computeVolumeAroundLayer(GaussPoint *gp, int layer)
{
    FloatArray G1, G2, G3, temp;
    double detJ;
    this->evalInitialCovarBaseVectorsAt(gp, G1, G2, G3);
    temp.beVectorProductOf(G1, G2);
    LayeredCrossSection *layeredCS = dynamic_cast< LayeredCrossSection * >( this->giveCrossSection() );
    detJ = temp.dotProduct(G3) * 0.5 * layeredCS->giveLayerThickness(layer);
    return detJ * gp->giveWeight();
}
#endif




// Transformation of bases and transformation matrices
# if 1
void
Shell7Base :: transInitialCartesianToInitialContravar(GaussPoint *gp, const FloatArray &VoightMatrix, FloatArray &answer)
{
    // Transform stress in cart system to curvilinear sytem
    // Uses Bond transformation matrix. No need to go from matrix to Voigt form and back.
    FloatArray Gcon1, Gcon2, Gcon3;
    this->evalInitialContravarBaseVectorsAt(gp, Gcon1, Gcon2, Gcon3);
    FloatMatrix GE, M;
    this->giveCoordTransMatrix(GE, Gcon1, Gcon2, Gcon3);
    this->giveBondTransMatrix(M, GE);
    answer.beProductOf(M, VoightMatrix);
}

void
Shell7Base :: transInitialCartesianToInitialContravar(GaussPoint *gp, const FloatMatrix &stiffness, FloatMatrix &answer)
{
    // Transform stifness tensor in cart system to curvilinear sytem
    // Uses Bond transformation matrix. No need to go from matrix to Voigt form and back.
    FloatArray Gcon1, Gcon2, Gcon3;
    this->evalInitialContravarBaseVectorsAt(gp, Gcon1, Gcon2, Gcon3);
    FloatMatrix GE, M, temp;
    this->giveCoordTransMatrix(GE, Gcon1, Gcon2, Gcon3);
    this->giveBondTransMatrix(M, GE);
    temp.beProductTOf(stiffness, M);
    answer.beProductOf(M, temp);
}

void
Shell7Base :: giveCoordTransMatrix(FloatMatrix &answer, FloatArray &g1, FloatArray &g2, FloatArray &g3,
                                   FloatArray &G1, FloatArray &G2, FloatArray &G3)
{   // Q_ij = dotproduct(g_i, G_j)
    answer.resize(3, 3);
    answer.at(1, 1) = g1.dotProduct(G1);
    answer.at(1, 2) = g1.dotProduct(G2);
    answer.at(1, 3) = g1.dotProduct(G3);
    answer.at(2, 1) = g2.dotProduct(G1);
    answer.at(2, 2) = g2.dotProduct(G2);
    answer.at(2, 3) = g2.dotProduct(G3);
    answer.at(3, 1) = g3.dotProduct(G1);
    answer.at(3, 2) = g3.dotProduct(G2);
    answer.at(3, 3) = g3.dotProduct(G3);
}

void
Shell7Base :: giveCoordTransMatrix(FloatMatrix &answer, FloatArray &g1, FloatArray &g2, FloatArray &g3)
{
    // Assumes transformation from a cartesian system
    // Q_ij = dotproduct(g_i, E_j) with E1 = [1 0 0]^T etc
    answer.resize(3, 3);
    answer.at(1, 1) = g1.at(1);
    answer.at(1, 2) = g1.at(2);
    answer.at(1, 3) = g1.at(3);
    answer.at(2, 1) = g2.at(1);
    answer.at(2, 2) = g2.at(2);
    answer.at(2, 3) = g2.at(3);
    answer.at(3, 1) = g3.at(1);
    answer.at(3, 2) = g3.at(2);
    answer.at(3, 3) = g3.at(3);
}

void
Shell7Base :: giveBondTransMatrix(FloatMatrix &answer, FloatMatrix &Q)
{
    /* Returns the Bond transformation matrix (M) between coordinate system 1 and 2, defined by the three base
     * vectors in each coordinate system g1i and g2i respectively.
     * Q is a transformation matrix defined as the direction cosines between the transformed system and the initial
     * system.*/


    answer.resize(6, 6);
    // block matrix 11 = Q^2
    answer.at(1, 1) = Q.at(1, 1) * Q.at(1, 1);
    answer.at(1, 2) = Q.at(1, 2) * Q.at(1, 2);
    answer.at(1, 3) = Q.at(1, 3) * Q.at(1, 3);
    answer.at(2, 1) = Q.at(2, 1) * Q.at(2, 1);
    answer.at(2, 2) = Q.at(2, 2) * Q.at(2, 2);
    answer.at(2, 3) = Q.at(2, 3) * Q.at(2, 3);
    answer.at(3, 1) = Q.at(3, 1) * Q.at(3, 1);
    answer.at(3, 2) = Q.at(3, 2) * Q.at(3, 2);
    answer.at(3, 3) = Q.at(3, 3) * Q.at(3, 3);

    // block matrix 12
    answer.at(1, 4) = 2. * Q.at(1, 2) * Q.at(1, 3);
    answer.at(1, 5) = 2. * Q.at(1, 3) * Q.at(1, 1);
    answer.at(1, 6) = 2. * Q.at(1, 1) * Q.at(1, 2);
    answer.at(2, 4) = 2. * Q.at(2, 2) * Q.at(2, 3);
    answer.at(2, 5) = 2. * Q.at(2, 3) * Q.at(2, 1);
    answer.at(2, 6) = 2. * Q.at(2, 1) * Q.at(2, 2);
    answer.at(3, 4) = 2. * Q.at(3, 2) * Q.at(3, 3);
    answer.at(3, 5) = 2. * Q.at(3, 3) * Q.at(3, 1);
    answer.at(3, 6) = 2. * Q.at(3, 1) * Q.at(3, 2);

    // block matrix 21
    answer.at(4, 1) = Q.at(2, 1) * Q.at(3, 1);
    answer.at(4, 2) = Q.at(2, 2) * Q.at(3, 2);
    answer.at(4, 3) = Q.at(2, 3) * Q.at(3, 3);
    answer.at(5, 1) = Q.at(3, 1) * Q.at(1, 1);
    answer.at(5, 2) = Q.at(3, 2) * Q.at(1, 2);
    answer.at(5, 3) = Q.at(3, 3) * Q.at(1, 3);
    answer.at(6, 1) = Q.at(1, 1) * Q.at(2, 1);
    answer.at(6, 2) = Q.at(1, 2) * Q.at(2, 2);
    answer.at(6, 3) = Q.at(1, 3) * Q.at(2, 3);

    // block matrix 22
    answer.at(4, 4) = Q.at(2, 2) * Q.at(3, 3) + Q.at(2, 3) * Q.at(3, 2);
    answer.at(4, 5) = Q.at(2, 1) * Q.at(3, 3) + Q.at(2, 3) * Q.at(3, 1);
    answer.at(4, 6) = Q.at(2, 2) * Q.at(3, 1) + Q.at(2, 1) * Q.at(3, 2);
    answer.at(5, 4) = Q.at(1, 2) * Q.at(3, 3) + Q.at(1, 3) * Q.at(3, 2);
    answer.at(5, 5) = Q.at(1, 3) * Q.at(3, 1) + Q.at(1, 1) * Q.at(3, 3);
    answer.at(5, 6) = Q.at(1, 1) * Q.at(3, 2) + Q.at(1, 2) * Q.at(3, 1);
    answer.at(6, 4) = Q.at(1, 2) * Q.at(2, 3) + Q.at(1, 3) * Q.at(2, 2);
    answer.at(6, 5) = Q.at(1, 3) * Q.at(2, 1) + Q.at(1, 1) * Q.at(2, 3);
    answer.at(6, 6) = Q.at(1, 1) * Q.at(2, 2) + Q.at(1, 2) * Q.at(2, 1);
}

#endif





// Recovery of nodal values

#if 1

void Shell7Base :: NodalAveragingRecoveryMI_computeSideValue(FloatArray &answer, int side, InternalStateType type, TimeStep *tStep)
{
    answer.resize(0);
}

int Shell7Base :: NodalAveragingRecoveryMI_giveDofManRecordSize(InternalStateType type)
{
    if ( type == IST_DirectorField ) {
        return 3;
    }

    return 0;
}

void Shell7Base :: NodalAveragingRecoveryMI_computeNodalValue(FloatArray &answer, int node, InternalStateType type, TimeStep *tStep)
{
    if ( type == IST_DirectorField ) {
        answer.resize(3);
        answer = this->giveInitialNodeDirector(node);
        answer.at(1) += this->giveNode(node)->giveDofWithID(W_u)->giveUnknown(EID_MomentumBalance, VM_Total, tStep);
        answer.at(2) += this->giveNode(node)->giveDofWithID(W_v)->giveUnknown(EID_MomentumBalance, VM_Total, tStep);
        answer.at(3) += this->giveNode(node)->giveDofWithID(W_w)->giveUnknown(EID_MomentumBalance, VM_Total, tStep);
        answer.times( this->giveCrossSection()->give(CS_Thickness) );
    } else {
        answer.resize(0);
    }
}

#endif




// Computation of solution vectors

#if 1

void
Shell7Base :: temp_computeVectorOf(IntArray &dofIdArray, ValueModeType u, TimeStep *stepN, FloatArray &answer)
{
    // Routine to extract vector given an array of dofid items
    // If a certain dofId does not exist a zero is used as value
    answer.resize( dofIdArray.giveSize() * numberOfDofMans );
    answer.zero();
    int k = 0;
    for ( int i = 1; i <= numberOfDofMans; i++ ) {
        DofManager *dMan = this->giveDofManager(i);        
        for (int j = 1; j <= dofIdArray.giveSize(); j++ ) {
            Dof *d = dMan->giveDof(j);
            k++;
            if ( dMan->hasDofID( (DofIDItem) dofIdArray.at(j) ) ) {
                answer.at(k) = d->giveUnknown(EID_MomentumBalance, VM_Total, stepN); ///@todo EID_MomentumBalance is just a dummy argument in this case and feels redundant
            }
        }
    }
}

//(int boundary, EquationID type, ValueModeType u, TimeStep *stepN, FloatArray &answer)
void
Shell7Base :: temp_computeBoundaryVectorOf(IntArray &dofIdArray, int boundary, ValueModeType u, TimeStep *stepN, FloatArray &answer)
{
    // Routine to extract vector given an array of dofid items
    // If a certain dofId does not exist a zero is used as value


    IntArray bNodes;
    this->giveInterpolation()->boundaryGiveNodes(bNodes, boundary);

    answer.resize( dofIdArray.giveSize() * bNodes.giveSize() );
    answer.zero();
    int k = 0;
    for ( int i = 1; i <= bNodes.giveSize(); i++ ) {
        DofManager *dMan = this->giveDofManager(bNodes.at(i));        
        for (int j = 1; j <= dofIdArray.giveSize(); j++ ) {
            Dof *d = dMan->giveDof(j);
            k++;
            if ( dMan->hasDofID( (DofIDItem) dofIdArray.at(j) ) ) {
                answer.at(k) = d->giveUnknown(EID_MomentumBalance, VM_Total, stepN); ///@todo EID_MomentumBalance is just a dummy argument in this case and feels redundant
            }
        }
    }


}

void
Shell7Base :: giveUpdatedSolutionVector(FloatArray &answer, TimeStep *tStep)
{
    // Computes updated solution as: x = X + dx, m = M + dM, gam = 0 + dgam
    this->giveInitialSolutionVector(answer);
    FloatArray temp;
    //this->computeVectorOf(EID_MomentumBalance, VM_Total, tStep, temp);

    int dummy = 0;
    IntArray dofIdArray;
    Shell7Base :: giveDofManDofIDMask(dummy, EID_MomentumBalance, dofIdArray);
    this->temp_computeVectorOf(dofIdArray, VM_Total, tStep, temp);
    
    answer.assemble( temp, this->giveOrdering(AllInv) );
}

void
Shell7Base :: giveUpdatedSolutionVectorC(FloatArray &answer, TimeStep *tStep)
{
    answer.resize( Shell7Base :: giveNumberOfDofs() );
    answer.zero();
    FloatArray temp;
    int dummy = 0;
    IntArray dofIdArray;
    Shell7Base :: giveDofManDofIDMask(dummy, EID_MomentumBalance, dofIdArray);
    this->temp_computeVectorOf(dofIdArray, VM_Total, tStep, temp);
    
    answer.assemble( temp, this->giveOrdering(AllInv) );
}

void
Shell7Base :: giveInitialSolutionVector(FloatArray &answer) {
    answer.resize( Shell7Base :: giveNumberOfDofs() );
    answer.zero();
    int ndofs_xm  = this->giveNumberOfFieldDofs(Midplane);
    // Reference position and directors
    for ( int i = 1, j = 0; i <= this->giveNumberOfDofManagers(); i++, j += 3 ) {
        FloatArray *Xi = this->giveNode(i)->giveCoordinates();
        FloatArray Mi = this->giveInitialNodeDirector(i);
        answer.at(1 + j) = Xi->at(1);
        answer.at(2 + j) = Xi->at(2);
        answer.at(3 + j) = Xi->at(3);
        answer.at(ndofs_xm + 1 + j) = Mi.at(1);
        answer.at(ndofs_xm + 2 + j) = Mi.at(2);
        answer.at(ndofs_xm + 3 + j) = Mi.at(3);
        // Assumes gam=0 at t=0
    }
}


void
Shell7Base :: edgeGiveUpdatedSolutionVector(FloatArray &answer, const int iedge, TimeStep *tStep)
{
    this->edgeGiveInitialSolutionVector(answer, iedge);
    FloatArray temp;
    //this->computeBoundaryVectorOf(iedge, EID_MomentumBalance, VM_Total, tStep, temp);
    int dummy = 0;
    IntArray dofIdArray;
    Shell7Base :: giveDofManDofIDMask(dummy, EID_MomentumBalance, dofIdArray);
    this->temp_computeBoundaryVectorOf(dofIdArray, iedge, VM_Total, tStep, temp);
    answer.assemble( temp, this->giveOrdering(EdgeInv) );
}

void
Shell7Base :: edgeGiveInitialSolutionVector(FloatArray &answer, const int iedge)
{
    FEInterpolation3d *fei = static_cast< FEInterpolation3d * >( this->giveInterpolation() );
    answer.resize( this->giveNumberOfEdgeDofs() );
    answer.zero();
    IntArray edgeNodes;
    fei->computeLocalEdgeMapping(edgeNodes, iedge);
    int ndofs_x = this->giveFieldSize(Midplane) * edgeNodes.giveSize();
    for ( int i = 1, j = 0; i <= edgeNodes.giveSize(); i++, j += 3 ) {
        FloatArray *Xi = this->giveNode( edgeNodes.at(i) )->giveCoordinates();
        FloatArray Mi  = this->giveInitialNodeDirector( edgeNodes.at(i) );
        answer.at(1 + j)   = Xi->at(1);
        answer.at(2 + j)   = Xi->at(2);
        answer.at(3 + j)   = Xi->at(3);
        answer.at(ndofs_x + 1 + j) = Mi.at(1);
        answer.at(ndofs_x + 2 + j) = Mi.at(2);
        answer.at(ndofs_x + 3 + j) = Mi.at(3);
        // gam(t=0)=0 is assumed
    }
}


void
Shell7Base :: giveGeneralizedStrainComponents(FloatArray genEps, FloatArray &dphidxi1, FloatArray &dphidxi2, FloatArray &dmdxi1,
                                              FloatArray &dmdxi2, FloatArray &m, double &dgamdxi1, double &dgamdxi2, double &gam) {
    // generealized strain vector  [dxdxi, dmdxi, m, dgamdxi, gam]^T
    dphidxi1.setValues( 3, genEps.at(1), genEps.at(2), genEps.at(3) );
    dphidxi2.setValues( 3, genEps.at(4), genEps.at(5), genEps.at(6) );
    dmdxi1.setValues( 3, genEps.at(7), genEps.at(8), genEps.at(9) );
    dmdxi2.setValues( 3, genEps.at(10), genEps.at(11), genEps.at(12) );
    m.setValues( 3, genEps.at(13), genEps.at(14), genEps.at(15) );
    dgamdxi1 = genEps.at(16);
    dgamdxi2 = genEps.at(17);
    gam = genEps.at(18);
}


void
Shell7Base :: computeGeneralizedStrainVector(FloatArray &answer, const FloatArray &solVec, const FloatMatrix &B11,
                                             const FloatMatrix &B22, const FloatMatrix &B32, const FloatMatrix &B43, const FloatMatrix  &B53) {
    answer.resize(18);
    answer.zero();
    int ndofs_xm  = this->giveNumberOfFieldDofs(Midplane);
    int ndofs_gam = this->giveNumberOfFieldDofs(InhomStrain);
    for ( int i = 1; i <= ndofs_xm; i++ ) {
        for ( int j = 1; j <= 6; j++ ) {
            answer.at(j)  += B11.at(j, i) * solVec.at(i);            // dx/dxi
            answer.at(6 + j)  += B22.at(j, i) * solVec.at(i + ndofs_xm);      // dm/dxi
        }
    }

    for ( int i = 1; i <= ndofs_xm; i++ ) {
        for ( int j = 1; j <= 3; j++ ) {
            answer.at(12 + j) += B32.at(j, i) * solVec.at(i + ndofs_xm);      // m
        }
    }

    for ( int i = 1; i <= ndofs_gam; i++ ) {
        for ( int j = 1; j <= 2; j++ ) {
            answer.at(15 + j) += B43.at(j, i) * solVec.at(i + ndofs_xm * 2);    // dgamma/dxi
        }
    }

    for ( int i = 1; i <= ndofs_gam; i++ ) {
        answer.at(18) += B53.at(1, i) * solVec.at(i + ndofs_xm * 2);  // gamma
    }
}


void
Shell7Base :: computeSolutionFields(FloatArray &xbar, FloatArray &m, double &gam, const FloatArray &solVec, const FloatMatrix &N11, const FloatMatrix &N22, const FloatMatrix &N33)
{
    // Returns xbar, m and gamma at a point corresponding to the N matrices
    xbar.resize(3);
    xbar.zero();
    m.resize(3);
    m.zero();
    int ndofs_xm  = this->giveNumberOfFieldDofs(Midplane);
    int ndofs_gam = this->giveNumberOfFieldDofs(InhomStrain);
    for ( int i = 1; i <= ndofs_xm; i++ ) {
        for ( int j = 1; j <= 3; j++ ) {
            xbar.at(j) += N11.at(j, i) * solVec.at(i);
            m.at(j) += N22.at(j, i) * solVec.at(i + ndofs_xm);
        }
    }

    gam = 0.;
    for ( int i = 1; i <= ndofs_gam; i++ ) {
        gam += N33.at(1, i) * solVec.at(i + ndofs_xm * 2);
    }
}


#endif




// N and B matrices

#if 1



void
Shell7Base :: edgeComputeNmatrixAt(GaussPoint *gp, FloatMatrix &answer)
// Returns the displacement interpolation matrix {N} of the receiver, eva-
// luated at aGaussPoint along one edge
{
    answer.resize( 7, this->giveNumberOfEdgeDofs() );
    answer.zero();

    FloatArray lcoords = * gp->giveCoordinates();
    FloatArray N;

    FEInterpolation3d *fei = static_cast< FEInterpolation3d * >( this->giveInterpolation() );
    fei->edgeEvalN( N, lcoords, FEIElementGeometryWrapper(this) );


    /*    9   9    3
     * 3 [N_x   0    0
     * 3   0   N_m   0
     * 1   0    0  N_gmm ]
     */
    int ndofs_xm = this->giveNumberOfEdgeDofs() / 7 * 3;   // numEdgeNodes * 3 dofs
    for ( int i = 1, j = 0; i <= this->giveNumberOfEdgeDofManagers(); i++, j += 3 ) {
        answer.at(1, 1 + j)   = N.at(i);
        answer.at(2, 2 + j)   = N.at(i);
        answer.at(3, 3 + j)   = N.at(i);
        answer.at(4, ndofs_xm + 1 + j) = N.at(i);
        answer.at(5, ndofs_xm + 2 + j) = N.at(i);
        answer.at(6, ndofs_xm + 3 + j) = N.at(i);
        answer.at(7, ndofs_xm * 2 + i)  = N.at(i);
    }
}


void
Shell7Base :: edgeComputeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int li, int ui)
/* Returns the  matrix {B} of the receiver, evaluated at aGaussPoint. Such that
 * B*a = [dxbar_dxi, dwdxi, w, dgamdxi, gam]^T, where a is the vector of unknowns
 */
{
    answer.resize( 11, this->giveNumberOfEdgeDofs() );
    answer.zero();
    FloatArray lcoords = * gp->giveCoordinates();
    FloatArray N, dNdxi;

    FEInterpolation3d *fei = static_cast< FEInterpolation3d * >( this->giveInterpolation() );
    fei->edgeEvalN( N, lcoords, FEIElementGeometryWrapper(this) );
    int iedge = 0;
    fei->edgeEvaldNdxi( dNdxi, iedge, lcoords, FEIElementGeometryWrapper(this) );

    /*
     * 3 [B_u   0    0
     * 3   0   B_w   0
     * 3   0   N_w   0
     * 1   0    0  B_gam
     * 1   0    0  N_gam]
     */
    int ndofs_xm = this->giveNumberOfEdgeDofs() / 7 * 3;   // numEdgeNodes * 3 dofs
    int ndofman = this->giveNumberOfEdgeDofManagers();
    // First row
    for ( int i = 1, j = 0; i <= ndofman; i++, j += 3  ) {
        answer.at(1, 1 + j) = dNdxi.at(i);
        answer.at(2, 2 + j) = dNdxi.at(i);
        answer.at(3, 3 + j) = dNdxi.at(i);
    }

    // Second row
    for ( int i = 1, j = 0; i <= ndofman; i++, j += 3  ) {
        answer.at(4, ndofs_xm + 1 + j) = dNdxi.at(i);
        answer.at(5, ndofs_xm + 2 + j) = dNdxi.at(i);
        answer.at(6, ndofs_xm + 3 + j) = dNdxi.at(i);
        answer.at(7, ndofs_xm + 1 + j) = N.at(i);
        answer.at(8, ndofs_xm + 2 + j) = N.at(i);
        answer.at(9, ndofs_xm + 3 + j) = N.at(i);
    }

    // Third row
    for ( int i = 1, j = 0; i <= ndofman; i++, j += 1  ) {
        answer.at(10, ndofs_xm * 2 + 1 + j) = dNdxi.at(i);
        answer.at(11, ndofs_xm * 2 + 1 + j) = N.at(i);
    }
}


void
Shell7Base :: computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int li, int ui)
/* Returns the  matrix {B} of the receiver, evaluated at aGaussPoint. Such that
 * B*a = [dxbar_dxi, dwdxi, w, dgamdxi, gam]^T, where a is the vector of unknowns
 */
{
    //int ndofs = this->giveNumberOfDofs();
    int ndofs = Shell7Base :: giveNumberOfDofs();
    int ndofs_xm  = this->giveNumberOfFieldDofs(Midplane);
    answer.resize(18, ndofs);
    answer.zero();
    FloatArray lcoords = * gp->giveCoordinates();
    FloatArray N;
    FloatMatrix dNdxi;
    this->giveInterpolation()->evalN( N, lcoords, FEIElementGeometryWrapper(this) );
    this->giveInterpolation()->evaldNdxi( dNdxi, lcoords, FEIElementGeometryWrapper(this) );

    /*    18   18   6
     * 6 [B_u   0   0
     * 6   0   B_w  0
     * 3   0   N_w  0
     * 2   0    0  B_gam
     * 1   0    0  N_gam]
     */
    int ndofman = this->giveNumberOfDofManagers();

    // First column
    for ( int i = 1, j = 0; i <= ndofman; i++, j += 3 ) {
        answer.at(1, 1 + j) = dNdxi.at(i, 1);
        answer.at(2, 2 + j) = dNdxi.at(i, 1);
        answer.at(3, 3 + j) = dNdxi.at(i, 1);
        answer.at(4, 1 + j) = dNdxi.at(i, 2);
        answer.at(5, 2 + j) = dNdxi.at(i, 2);
        answer.at(6, 3 + j) = dNdxi.at(i, 2);
    }

    // Second column
    for ( int i = 1, j = 0; i <= ndofman; i++, j += 3 ) {
        answer.at(7, ndofs_xm + 1 + j) = dNdxi.at(i, 1);
        answer.at(8, ndofs_xm + 2 + j) = dNdxi.at(i, 1);
        answer.at(9, ndofs_xm + 3 + j) = dNdxi.at(i, 1);
        answer.at(10, ndofs_xm + 1 + j) = dNdxi.at(i, 2);
        answer.at(11, ndofs_xm + 2 + j) = dNdxi.at(i, 2);
        answer.at(12, ndofs_xm + 3 + j) = dNdxi.at(i, 2);
        answer.at(13, ndofs_xm + 1 + j) = N.at(i);
        answer.at(14, ndofs_xm + 2 + j) = N.at(i);
        answer.at(15, ndofs_xm + 3 + j) = N.at(i);
    }

    // Third column
    for ( int i = 1, j = 0; i <= ndofman; i++, j += 1 ) {
        answer.at(16, ndofs_xm * 2 + 1 + j) = dNdxi.at(i, 1);
        answer.at(17, ndofs_xm * 2 + 1 + j) = dNdxi.at(i, 2);
        answer.at(18, ndofs_xm * 2 + 1 + j) = N.at(i);
    }
}


void
Shell7Base :: computeNmatrixAt(GaussPoint *gp, FloatMatrix &answer)
// Returns the displacement interpolation matrix {N} of the receiver, eva-
// luated at aGaussPoint.
{
    int ndofs = this->giveNumberOfDofs();
    int ndofs_xm  = this->giveNumberOfFieldDofs(Midplane);
    answer.resize(7, ndofs);
    answer.zero();
    FloatArray &lcoords = * gp->giveCoordinates();
    FloatArray N;
    this->giveInterpolation()->evalN( N, lcoords, FEIElementGeometryWrapper(this) );

    /*   nno*3 nno*3 nno
     * 3 [N_x   0    0
     * 3   0   N_m   0
     * 1   0    0  N_gmm ]
     */
    for ( int i = 1, j = 0; i <= this->giveNumberOfDofManagers(); i++, j += 3 ) {
        answer.at(1, 1 + j) = N.at(i);
        answer.at(2, 2 + j) = N.at(i);
        answer.at(3, 3 + j) = N.at(i);
        answer.at(4, ndofs_xm + 1 + j) = N.at(i);
        answer.at(5, ndofs_xm + 2 + j) = N.at(i);
        answer.at(6, ndofs_xm + 3 + j) = N.at(i);
        answer.at(7, ndofs_xm * 2 + i) = N.at(i);
    }
}



int
Shell7Base :: giveFieldSize(SolutionField fieldType)
{
    if ( fieldType == Midplane ) {
        return 3;
    } else if ( fieldType == Director  ) {
        return 3;
    } else if ( fieldType == InhomStrain  ) {
        return 1;
    } else if ( fieldType == All  ) {
        return 7;
    } else {
        _error("giveFieldSize: unknown fieldType");
        return 0;
    }
}

int
Shell7Base :: giveNumberOfFieldDofs(SolutionField fieldType)
{
    return this->giveNumberOfDofManagers() * giveFieldSize(fieldType);
}



void
Shell7Base :: computeFieldBmatrix(FloatMatrix &answer, FloatMatrix &dNdxi, SolutionField fieldType)
{
    // Returns the B matrix corresponding to a given solution field
    int ndofman = this->giveNumberOfDofManagers();
    int fieldSize = this->giveFieldSize(fieldType);
    answer.resize(2 * fieldSize, fieldSize * ndofman);
    answer.zero();

    for ( int row  = 1; row <= fieldSize; row++  ) {
        for ( int i = 1, j = 0; i <= ndofman; i++, j += fieldSize ) {
            answer.at(row, row + j) = dNdxi.at(i, 1);
            answer.at(row + fieldSize, row + j) = dNdxi.at(i, 2);
        }
    }
}


void
Shell7Base :: computeFieldNmatrix(FloatMatrix &answer, FloatArray &N, SolutionField fieldType)
{
    // Returns the N matrix corresponding to a given solution field
    int ndofman = this->giveNumberOfDofManagers();
    int fieldSize = this->giveFieldSize(fieldType);
    answer.resize(fieldSize, fieldSize * ndofman);
    answer.zero();

    for ( int i = 1, j = 0; i <= ndofman; i++, j += fieldSize ) {
        for ( int row  = 1; row <= fieldSize; row++  ) {
            answer.at(row, row + j) = N.at(i);
        }
    }
}


void
Shell7Base :: computeBmatricesAt(GaussPoint *gp, FloatMatrix &B11, FloatMatrix &B22, FloatMatrix &B32, FloatMatrix &B43, FloatMatrix &B53)
{
    /* Returns the submatrices of the B matrix, evaluated at gp.
     *   nno*3 nno*3 nno
     * 6 [B_u   0    0       [B11   0    0
     * 6   0   B_w   0         0   B22   0
     * 3   0   N_w   0      =  0   B32   0
     * 2   0    0   B_gam      0    0   B43
     * 1   0    0   N_gam]     0    0   B53]
     */
    FloatArray lcoords = * gp->giveCoordinates();
    FloatArray N;
    FloatMatrix dNdxi;
    this->giveInterpolation()->evalN( N, lcoords, FEIElementGeometryWrapper(this) );
    this->giveInterpolation()->evaldNdxi( dNdxi, lcoords, FEIElementGeometryWrapper(this) );

    this->computeFieldBmatrix(B11, dNdxi, Midplane);
    B22 = B11;
    this->computeFieldNmatrix(B32, N, Director);
    this->computeFieldBmatrix(B43, dNdxi, InhomStrain);
    this->computeFieldNmatrix(B53, N, InhomStrain);
}

void
Shell7Base :: computeNmatricesAt(GaussPoint *gp, FloatMatrix &N11, FloatMatrix &N22, FloatMatrix &N33)
{
    /* Returns the submatrices of the N matrix, evaluated at gp.
     *   nno*3 nno*3 nno
     * 3 [N_x   0    0        [N11   0    0
     * 3   0   N_m   0      =   0   N22   0
     * 1   0    0  N_gmm ]      0    0   N33]
     */
    FloatArray lcoords = * gp->giveCoordinates();
    FloatArray N;
    this->giveInterpolation()->evalN( N, lcoords, FEIElementGeometryWrapper(this) );

    this->computeFieldNmatrix(N11, N, Midplane);
    N22 = N11;
    this->computeFieldNmatrix(N33, N, InhomStrain);
}


#endif




// Misc functions

int
Shell7Base :: giveVoigtIndex(int ind1, int ind2)
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
} // end namespace oofem
