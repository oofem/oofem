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

#include "../sm/Elements/Shells/shell7base.h"
#include "../sm/Materials/structuralms.h"
#include "../sm/Loads/constantpressureload.h"
#include "node.h"
#include "load.h"
#include "mathfem.h"
#include "domain.h"
#include "gaussintegrationrule.h"
#include "gausspoint.h"
#include "boundaryload.h"
#include "constantsurfaceload.h"
#include "vtkxmlexportmodule.h"
#include "fracturemanager.h"
#include "dof.h"
#include <fstream>

namespace oofem {

FEI3dTrQuad  Shell7Base   :: interpolationForCZExport;
FEI3dWedgeQuad Shell7Base :: interpolationForExport;


Shell7Base :: Shell7Base(int n, Domain *aDomain) : NLStructuralElement(n, aDomain),  LayeredCrossSectionInterface(), 
    VTKXMLExportModuleElementInterface(), ZZNodalRecoveryModelInterface(this), FailureModuleElementInterface(){}

IRResultType Shell7Base :: initializeFrom(InputRecord *ir)
{
    return NLStructuralElement :: initializeFrom(ir);

}

int 
Shell7Base :: checkConsistency()
{
    NLStructuralElement :: checkConsistency();

    if ( this->layeredCS == NULL ) {
        OOFEM_ERROR("Elements derived from Shell7Base only supports layered cross section");
    }
    return ( this->layeredCS != NULL  &&  this->fei != NULL);
}


void
Shell7Base :: postInitialize()
{
    
    this->layeredCS = dynamic_cast< LayeredCrossSection * >( this->giveCrossSection()  );
    if ( this->layeredCS == NULL ) {
        OOFEM_ERROR("Shell7Base derived elements only supports layered cross section");
    }
    this->fei = dynamic_cast< FEInterpolation3d   * >( this->giveInterpolation() );
    this->setupInitialNodeDirectors();
    this->setupInitialSolutionVector();
    this->setupInitialEdgeSolutionVector();
    Element :: postInitialize();
    
        
    this->voigtIndices = {
        { 1, 6, 5 },
        { 9, 2, 4 },
        { 8, 7, 3 }
    };
    
    this->nlGeometry = 1;
}


Interface *Shell7Base :: giveInterface(InterfaceType it)
{
    switch ( it ) {
    case NodalAveragingRecoveryModelInterfaceType:
        return static_cast< NodalAveragingRecoveryModelInterface * >( this );

    case LayeredCrossSectionInterfaceType:
        return static_cast< LayeredCrossSectionInterface * >( this );

    case VTKXMLExportModuleElementInterfaceType:
        return static_cast< VTKXMLExportModuleElementInterface * >( this );

    case ZZNodalRecoveryModelInterfaceType:
        return static_cast< ZZNodalRecoveryModelInterface * >( this );

    case FailureModuleElementInterfaceType:
        return static_cast< FailureModuleElementInterface * >( this );

    default:
        return StructuralElement :: giveInterface(it);
    }
}


void
Shell7Base :: giveDofManDofIDMask(int inode, IntArray &answer) const
{
    answer = { D_u, D_v, D_w, W_u, W_v, W_w, Gamma };
}


int
Shell7Base :: giveNumberOfDofs() 
{
    return 7 * this->giveNumberOfDofManagers();
}

int
Shell7Base :: computeGlobalCoordinates(FloatArray &answer, const FloatArray &lcoords)
{
    // should it return coord in reference or updated config? -- Let it be initial so small def code remains the same
    //@todo move VTK method into this
    int layer = 1;
    vtkEvalInitialGlobalCoordinateAt(lcoords, layer, answer);

    return 1;
}

int
Shell7Base :: computeGlobalCoordinatesOnEdge(FloatArray &answer, const FloatArray &lcoords, const int iEdge)
{
    // should it return coord in reference or updated config? -- Let it be initial so small def code remains the same
//     double zeta = giveGlobalZcoordInLayer(lcoords.at(3), layer);
//     FloatArray N;
//     this->fei->evalN( N, lcoords, FEIElementGeometryWrapper(this) );
// 
//     globalCoords.clear();
//     for ( int i = 1; i <= this->giveNumberOfDofManagers(); i++ ) {
//         FloatArray &xbar = *this->giveNode(i)->giveCoordinates();
//         const FloatArray &M = this->giveInitialNodeDirector(i);
//         globalCoords.add(N.at(i), ( xbar + zeta * M ));
//     }
    return 1;
}



double
Shell7Base::giveGlobalZcoord( const FloatArray &lCoords )
{
    return lCoords.at(3) * this->layeredCS->give( CS_Thickness, lCoords, this, false ) * 0.5;
}


double // @todo move to layered crosssection
Shell7Base :: giveGlobalZcoordInLayer(double xi, int layer)
{
    // Mid position + part of the layer thickness
    return this->layeredCS->giveLayerMidZ(layer) + xi*this->layeredCS->giveLayerThickness(layer)*0.5 ;
}



// Base vectors

#if 1

void
Shell7Base :: evalInitialCovarBaseVectorsAt(const FloatArray &lcoords, FloatMatrix &Gcov)
{
    double zeta = giveGlobalZcoord( lcoords );
    FloatArray M;
    FloatMatrix dNdxi;

    // In plane base vectors
    this->fei->evaldNdxi( dNdxi, lcoords, FEIElementGeometryWrapper(this) );

    FloatArray G1, G2, nodeCoords; 
    for ( int i = 1; i <= this->giveNumberOfDofManagers(); i++ ) {
        FloatArray &xbar = * this->giveNode(i)->giveCoordinates();
        M = this->giveInitialNodeDirector(i);
        nodeCoords = (xbar + zeta*M);
        G1.add(dNdxi.at(i, 1), nodeCoords);
        G2.add(dNdxi.at(i, 2), nodeCoords);
    }

    // Out of plane base vector = director
    FloatArray G3;
    this->evalInitialDirectorAt(lcoords, G3);     // G3=M

    Gcov.resize(3,3);
    Gcov.setColumn(G1,1); Gcov.setColumn(G2,2); Gcov.setColumn(G3,3);
}


void
Shell7Base :: edgeEvalInitialCovarBaseVectorsAt(const FloatArray &lcoords, const int iedge, FloatArray &G1, FloatArray &G3)
{
    double zeta = 0.0;     // no variation i z (yet)
    FloatArray M, dNdxi, nodeCoords;

    IntArray edgeNodes;
    this->fei->computeLocalEdgeMapping(edgeNodes, iedge);
    this->fei->edgeEvaldNdxi( dNdxi, iedge, lcoords, FEIElementGeometryWrapper(this) );

    // Base vector along edge
    G1.clear();
    for ( int i = 1; i <= edgeNodes.giveSize(); i++ ) {
        FloatArray &xbar = * this->giveNode(edgeNodes.at(i))->giveCoordinates();
        M = this->giveInitialNodeDirector(edgeNodes.at(i));
        nodeCoords = (xbar + zeta*M);
        G1.add(dNdxi.at(i), nodeCoords);
    }

    // Director will be the second base vector
    this->edgeEvalInitialDirectorAt(lcoords, G3, iedge);
}


void
Shell7Base :: evalInitialContravarBaseVectorsAt(const FloatArray &lCoords, FloatMatrix &Gcon)
{
    FloatMatrix Gcov;
    this->evalInitialCovarBaseVectorsAt(lCoords, Gcov);
    this->giveDualBase(Gcov, Gcon);
}


void
Shell7Base :: giveDualBase( FloatMatrix &base1, FloatMatrix &base2)
{
    // Computes the dual base through thte inversion of the metric tensor
    FloatMatrix gMetric, ginv;  
    gMetric.beTProductOf(base1,base1); // Metric tensor
    ginv.beInverseOf(gMetric);
    base2.beProductTOf(base1,ginv);
}


void
Shell7Base :: evalInitialDirectorAt(const FloatArray &lcoords, FloatArray &answer)
{   
    // Interpolates between the node directors
    FloatArray N;
    this->fei->evalN( N, lcoords, FEIElementGeometryWrapper(this) );
    answer.clear();
    for ( int i = 1; i <= this->giveNumberOfDofManagers(); i++ ) {
        answer.add( N.at(i), this->giveInitialNodeDirector(i) );
    }
}


void
Shell7Base :: edgeEvalInitialDirectorAt(const FloatArray &lcoords, FloatArray &answer, const int iEdge)
{   
    // Interpolates between the node directors along an edge

    FloatArray N;
    IntArray edgeNodes;
    this->fei->computeLocalEdgeMapping(edgeNodes, iEdge);
    this->fei->edgeEvalN( N, iEdge, lcoords, FEIElementGeometryWrapper(this) );

    answer.clear();
    for ( int i = 1; i <= edgeNodes.giveSize(); i++ ) {
        answer.add( N.at(i), this->giveInitialNodeDirector( edgeNodes.at(i) ) );
    }
}



void
Shell7Base :: setupInitialNodeDirectors()
{   
    // Compute directors as normals to the surface
    FloatArray M(3), G1(3), G2(3), lcoords(2);
    FloatMatrix localNodeCoords;
    this->giveInterpolation()->giveLocalNodeCoords(localNodeCoords);
    
    int nDofMan = this->giveNumberOfDofManagers();
    this->initialNodeDirectors.resize(nDofMan);
    FloatMatrix dNdxi;
    for ( int node = 1; node <= nDofMan; node++ ) {
        lcoords.beColumnOf(localNodeCoords,node);      
        this->fei->evaldNdxi( dNdxi, lcoords, FEIElementGeometryWrapper(this) );

        G1.zero();
        G2.zero();
        // base vectors of the initial surface
        for ( int i = 1; i <= nDofMan; i++ ) {        
            FloatArray *nodeI = this->giveNode(i)->giveCoordinates();
            G1.add(dNdxi.at(i, 1), * nodeI);
            G2.add(dNdxi.at(i, 2), * nodeI);
        }

        M.beVectorProductOf(G1, G2);
        M.normalize();
        this->initialNodeDirectors [ node - 1 ] = M;
    }
}



void
Shell7Base :: evalCovarBaseVectorsAt(const FloatArray &lcoords, FloatMatrix &gcov, FloatArray &genEps, TimeStep *tStep)
{
    // Evaluates the covariant base vectors in the current configuration
    FloatArray g1; FloatArray g2; FloatArray g3;
    double zeta = giveGlobalZcoord( lcoords );

    FloatArray dxdxi1, dxdxi2, m, dmdxi1, dmdxi2;
    double dgamdxi1, dgamdxi2, gam;
    this->giveGeneralizedStrainComponents(genEps, dxdxi1, dxdxi2, dmdxi1, dmdxi2, m, dgamdxi1, dgamdxi2, gam);

    double fac1 = ( zeta + 0.5 * gam * zeta * zeta );
    double fac2 = ( 0.5 * zeta * zeta );
    double fac3 = ( 1.0 + zeta * gam );

    g1 = dxdxi1 + fac1*dmdxi1 + fac2*dgamdxi1*m;
    g2 = dxdxi2 + fac1*dmdxi2 + fac2*dgamdxi2*m;
    g3 = fac3*m;

    gcov.resize(3,3);
    gcov.setColumn(g1,1); gcov.setColumn(g2,2); gcov.setColumn(g3,3);

}


void
Shell7Base :: edgeEvalCovarBaseVectorsAt(const FloatArray &lcoords, const int iedge, FloatMatrix &gcov, TimeStep *tStep)
{
    // Evaluates the covariant base vectors in the current configuration for an edge
    double zeta = lcoords.at(3);

    FloatArray solVecEdge;
    FloatMatrix B;
    IntArray edgeNodes;
    this->fei->computeLocalEdgeMapping(edgeNodes, iedge);
    this->edgeComputeBmatrixAt(lcoords, B, 1, ALL_STRAINS);
    this->edgeGiveUpdatedSolutionVector(solVecEdge, iedge, tStep);

    FloatArray genEpsEdge;                 // generalized strain
    genEpsEdge.beProductOf(B, solVecEdge); // [dxdxi, dmdxi, m, dgamdxi, gam]^T

    FloatArray dxdxi, m, dmdxi;
    dxdxi = { 3, genEpsEdge.at(1), genEpsEdge.at(2), genEpsEdge.at(3) };
    dmdxi = { genEpsEdge.at(4), genEpsEdge.at(5), genEpsEdge.at(6) };
    m = { genEpsEdge.at(7), genEpsEdge.at(8), genEpsEdge.at(9) };
    double dgamdxi = genEpsEdge.at(10);
    double gam     = genEpsEdge.at(11);

    double fac1 = ( zeta + 0.5 * gam * zeta * zeta );
    double fac2 = ( 0.5 * zeta * zeta );
    double fac3 = ( 1.0 + zeta * gam );
    
    FloatArray g1, g2, g3;
    g2 = dxdxi + fac1*dmdxi + fac2*dgamdxi*m; // base vector along the edge
    g3 = fac3*m;                              // director field

    g2.normalize();
    g3.normalize();
    g1.beVectorProductOf(g2, g3);
    g1.normalize();
    gcov.resize(3,3);
    gcov.setColumn(g1,1);
    gcov.setColumn(g2,2);
    gcov.setColumn(g3,3);
}


#endif




// Tangent matrices

#if 1

void
Shell7Base :: computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
{

    FloatArray solVec;
    this->giveUpdatedSolutionVector(solVec, tStep); // a
    
    this->computeBulkTangentMatrix(answer, solVec, tStep);


    // Add contribution due to pressure load ///@todo should later be compted by the load
    int nLoads = this->boundaryLoadArray.giveSize() / 2;

    for ( int i = 1; i <= nLoads; i++ ) {     // For each pressure load that is applied
        int load_number = this->boundaryLoadArray.at(2 * i - 1);
        int iSurf = this->boundaryLoadArray.at(2 * i);         // load_id
        Load *load = this->domain->giveLoad(load_number);

        if ( dynamic_cast< ConstantPressureLoad * >( load ) ) {
            FloatMatrix K_pressure;
            this->computePressureTangentMatrix(K_pressure, load, iSurf, tStep);
            answer.add(K_pressure); // Should assemble with ordering
        }
    }

}

void
Shell7Base :: computeLambdaGMatrices(FloatMatrix lambda [ 3 ], FloatArray &genEps, double zeta)
{
    // computes the lambda^g matrices associated with the variation and linearization of the base vectors g_i.
    // \delta g_i = lambda_i * \delta n and \Delta g_i = lambda_i * \Delta n 
    // \delta n = B * \delta a and \Delta n = B * \Delta a
    // @todo optimize method
    FloatArray m(3), dm1(3), dm2(3), temp1;
    double dgam1, dgam2, gam;
    this->giveGeneralizedStrainComponents(genEps, temp1, temp1, dm1, dm2, m, dgam1, dgam2, gam);

    // thickness coefficients
    double a = zeta + 0.5 * gam * zeta * zeta;
    double b = 0.5 * zeta * zeta;
    double c = 1.0 + gam * zeta;

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
    bm   = b*m; 
    bdm1 = b*dm1; 
    bdm2 = b*dm2; 
    lambda[ 0 ].setColumn(bm,16);   lambda[ 0 ].setColumn(bdm1,18);
    lambda[ 1 ].setColumn(bm,17);   lambda[ 1 ].setColumn(bdm2,18);

    // lambda3 =  ( 0,   0,   0 ,   0 ,     c*I   ,   0 ,   0 ,   xi*m )
    lambda[ 2 ].resize(3,18);   lambda[ 2 ].zero();
    lambda[ 2 ].at(1,13) = lambda[ 2 ].at(2,14) = lambda[ 2 ].at(3,15) = c;
    FloatArray zm(3);
    zm = zeta*m; 
    lambda[ 2 ].setColumn(zm,18);
}

void
Shell7Base :: computeLambdaNMatrix(FloatMatrix &lambda, FloatArray &genEps, double zeta)
{
    // computes the lambda^n matrix associated with the variation and linearization of the position vector x.
    // \delta x = lambda * \delta \hat{x} with \hat{x} = [\bar{x}, m, \gamma]

    FloatArray m(3);
    m = { genEps.at(13), genEps.at(14), genEps.at(15) };
    double gam = genEps.at(18);

    // thickness coefficients
    double a = zeta + 0.5 * gam * zeta * zeta;
    double b = 0.5 * zeta * zeta;
    
    // lambda =  ( I, a*I, b*m )
    lambda.resize(3,7);
    lambda.zero();
    lambda.at(1,1) = lambda.at(2,2) = lambda.at(3,3) = 1.0;
    lambda.at(1,4) = lambda.at(2,5) = lambda.at(3,6) = a;
    lambda.setColumn(b*m,7); 

}


void
Shell7Base :: computeBulkTangentMatrix(FloatMatrix &answer, FloatArray &solVec, TimeStep *tStep)
{
    FloatMatrix A [ 3 ] [ 3 ], lambda[ 3 ], A_lambda(3,18), LB;
    FloatMatrix L(18,18), B;
    FloatMatrix tempAnswer;

    int ndofs = Shell7Base :: giveNumberOfDofs();
    answer.resize(ndofs, ndofs); tempAnswer.resize(ndofs, ndofs);
    answer.zero(); tempAnswer.zero();

    int numberOfLayers = this->layeredCS->giveNumberOfLayers();     
    FloatArray genEps;

    for ( int layer = 1; layer <= numberOfLayers; layer++ ) {
        StructuralMaterial *mat = static_cast< StructuralMaterial* >( domain->giveMaterial( this->layeredCS->giveLayerMaterial(layer) ) );

        for ( GaussPoint *gp : *integrationRulesArray [ layer - 1 ] ) {
            const FloatArray &lCoords = gp->giveNaturalCoordinates();

            this->computeBmatrixAt(lCoords, B);
            genEps.beProductOf(B, solVec);
            // Material stiffness
            Shell7Base :: computeLinearizedStiffness(gp, mat, tStep, A);

            double zeta = giveGlobalZcoord(lCoords);
            this->computeLambdaGMatrices(lambda, genEps, zeta);

            // L = sum_{i,j} (lambdaI_i)^T * A^ij * lambdaJ_j
            // note: L will only be symmetric if lambdaI = lambdaJ (not the case for xfem)
            L.zero();
            for (int i = 0; i < 3; i++) {
                A_lambda.zero();
                for (int j = 0; j < 3; j++) {
                    A_lambda.addProductOf(A[i][j], lambda[j]);
                }
                L.plusProductSymmUpper(lambda[i], A_lambda, 1.0);
            }
            L.symmetrized();
            LB.beProductOf(L, B);
            double dV = this->computeVolumeAroundLayer(gp, layer);
            tempAnswer.plusProductSymmUpper(B, LB, dV);
            
        }
    }
    tempAnswer.symmetrized();
    
    const IntArray &ordering = this->giveOrderingDofTypes();
    answer.assemble(tempAnswer, ordering, ordering);

}

void
Shell7Base :: computeLinearizedStiffness(GaussPoint *gp, StructuralMaterial *mat, TimeStep *tStep, FloatMatrix A [ 3 ] [ 3 ]) 
{
    const FloatArray &lcoords = gp->giveNaturalCoordinates();
    FloatMatrix D;

    // Material stiffness when internal work is formulated in terms of P and F:
    // \Delta(P*G^I) = L^IJ * \Delta g_J
    // A[I][J] = L^IJ = L_klmn * [G^I]_l * [G^J]_n
    mat->give3dMaterialStiffnessMatrix_dPdF(D, TangentStiffness, gp, tStep);    // D_ijkl - cartesian system (Voigt)
    FloatMatrix G;
    this->evalInitialContravarBaseVectorsAt(lcoords, G);
    for (int I = 1; I <= 3; I++) {
        for (int J = I; J <= 3; J++) {
            A[I - 1][J - 1].resize(3, 3);
            A[I - 1][J - 1].zero();

            for (int k = 1; k <= 3; k++) {
                for (int l = 1; l <= 3; l++) {
                    for (int m = 1; m <= 3; m++) {
                        for (int n = 1; n <= 3; n++) {
                    
                            A[I-1][J-1].at(k, m) += D.at( giveVoigtIndex(k, l), giveVoigtIndex(m, n) ) * G.at(l, I) * G.at(n, J);

                        }
                    }
                }
            }
        }
    }

    // position 21
    A [ 1 ] [ 0 ].beTranspositionOf( A [ 0 ] [ 1 ] );

    // position 31
    A [ 2 ] [ 0 ].beTranspositionOf( A [ 0 ] [ 2 ] );
    
    // position 32
    A [ 2 ] [ 1 ].beTranspositionOf( A [ 1 ] [ 2 ] );



}

void
Shell7Base :: computePressureTangentMatrix(FloatMatrix &answer, Load *load, const int iSurf, TimeStep *tStep)
{
    // Computes tangent matrix associated with the linearization of pressure loading. Assumes constant pressure.
    ConstantPressureLoad* pLoad = dynamic_cast< ConstantPressureLoad * >( load );
    
    //std :: unique_ptr< IntegrationRule >iRule( giveInterpolation()->giveBoundaryIntegrationRule(load->giveApproxOrder(), iSurf) );
    int nPointsTri = 6; //todo generalize
    GaussIntegrationRule iRule(1, this);
    iRule.SetUpPointsOnWedge(nPointsTri, 1, _3dMat); ///@todo replce with triangle which has a xi3-coord
    
    FloatMatrix N, B, LB, NLB, L(7, 18), gcov, W1, W2;
    FloatArray lcoords(3), solVec, pressure;
    FloatArray g1, g2, genEps;
    FloatMatrix lambdaG [ 3 ], lambdaN;

    double xi   = pLoad->giveLoadOffset();
    this->giveUpdatedSolutionVector(solVec, tStep);

    int ndof = Shell7Base :: giveNumberOfDofs();
    answer.resize(ndof, ndof);
    answer.zero();
    for ( auto *ip: iRule ) { // rule #2 for surface integration
        lcoords.at(1) = ip->giveNaturalCoordinate(1);
        lcoords.at(2) = ip->giveNaturalCoordinate(2);
        lcoords.at(3) = xi;     // local coord where load is applied
        double zeta = giveGlobalZcoord( lcoords );

        this->computeNmatrixAt(lcoords, N);
        this->computeBmatrixAt(lcoords, B);
        genEps.beProductOf(B, solVec);     

        // Traction tangent, L =  lambdaN * ( W2*lambdaG_1 - W1*lambdaG_2  ) 
        load->computeValueAt(pressure, tStep, ip->giveNaturalCoordinates(), VM_Total);        // pressure component   
        this->evalCovarBaseVectorsAt(lcoords, gcov, genEps, tStep);
        g1.beColumnOf(gcov,1);
        g2.beColumnOf(gcov,2);
        W1 = this->giveAxialMatrix(g1);
        W2 = this->giveAxialMatrix(g2);
        
        this->computeLambdaGMatrices(lambdaG, genEps, zeta);
        this->computeLambdaNMatrix(lambdaN, genEps, zeta);
        FloatMatrix W2L, W1L;
        W2L.beProductOf(W2,lambdaG[0]);
        W1L.beProductOf(W1,lambdaG[1]);
        W2L.subtract(W1L);
        L.beTProductOf(lambdaN, W2L);
        L.times( -pressure.at(1) );

        // Tangent matrix (K = N^T*L*B*dA)
        LB.beProductOf(L,B);
        NLB.beTProductOf(N,LB);
        double dA = this->computeAreaAround(ip, xi);
        answer.add(dA, NLB);
    }
}


FloatMatrix 
Shell7Base :: giveAxialMatrix(const FloatArray &v)
{
    // creates the skew-symmetric matrix W defined such that 
    // crossProduct(u,v) = W(v)*u
    // W = [   0   v(3)  -v(2)
    //      -v(3)   0     v(1)
    //       v(2) -v(1)    0  ]
    //
    FloatMatrix answer(3,3);
    answer.zero();
    answer.at(2, 3) =  v.at(1);
    answer.at(3, 2) = -v.at(1);
    answer.at(1, 3) = -v.at(2);
    answer.at(3, 1) =  v.at(2);
    answer.at(1, 2) =  v.at(3);
    answer.at(2, 1) = -v.at(3);
    return answer;
}

#endif



// Strain and stress

#if 1

void
Shell7Base :: computeFAt(const FloatArray &lCoords, FloatMatrix &answer, FloatArray &genEps, TimeStep *tStep)
{
    // Computes the deformation gradient in matrix form as open product(g_i, G^i) = gcov*Gcon^T
    FloatMatrix gcov, Gcon;
    this->evalCovarBaseVectorsAt(lCoords, gcov, genEps, tStep);
    this->evalInitialContravarBaseVectorsAt(lCoords, Gcon);
    answer.beProductTOf(gcov, Gcon);
}

void
Shell7Base :: computeStressMatrix(FloatMatrix &answer, FloatArray &genEps, GaussPoint *gp, Material *mat, TimeStep *tStep)
{
    FloatMatrix F;
    FloatArray vF, vP;
    computeFAt(gp->giveNaturalCoordinates(), F, genEps, tStep);
    vF.beVectorForm(F);
    static_cast< StructuralMaterial * >( mat )->giveFirstPKStressVector_3d(vP, gp, vF, tStep);
    answer.beMatrixForm(vP);
  
}


void
Shell7Base :: computeCauchyStressVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep)
{
    // Compute Cauchy stress from 2nd Piola stress
    FloatArray solVec;
    this->giveUpdatedSolutionVector(solVec, tStep); 
    FloatMatrix  B;
    const FloatArray &lCoords = gp->giveNaturalCoordinates();
    this->computeBmatrixAt(lCoords, B);
    FloatArray genEps;

    genEps.beProductOf(B, solVec);
    FloatMatrix F;
    this->computeFAt(lCoords, F, genEps, tStep);   

    FloatArray vS;
    giveIPValue(vS, gp, IST_StressTensor, tStep); // expects Second PK stress

    FloatMatrix S, temp, sigma;
    S.beMatrixFormOfStress(vS);
    temp.beProductTOf(S,F); 
    sigma.beProductOf(F,temp);
    sigma.times( 1.0/F.giveDeterminant() );
    answer.beSymVectorForm(sigma);
}


int 
Shell7Base :: giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep)
{
    // Compute special IST quantities this element should support
    switch (type) {
    case IST_CauchyStressTensor:
        this->computeCauchyStressVector(answer, gp, tStep);
        return 1;
    default:
        return Element :: giveIPValue(answer, gp, type, tStep);
    }
}

#endif




// Internal forces

#if 1


void
Shell7Base :: giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord)
{
    // Computes internal forces as a summation of: sectional forces + convective mass force

    FloatArray solVec;
    this->giveUpdatedSolutionVector(solVec, tStep); // placement vector x
    this->computeSectionalForces(answer, tStep, solVec, useUpdatedGpRecord);

    ///@todo How to treat the convective force? Only active during dynamic simulations
}


void
Shell7Base :: computeSectionalForces(FloatArray &answer, TimeStep *tStep, FloatArray &solVec, int useUpdatedGpRecord)
{
    // Computation of sectional forces: f = \int_V B^t * \hat{N} dV = \int_V B^t*[N M T Ms Ts]^t dV 
    int ndofs = Shell7Base :: giveNumberOfDofs();

    int numberOfLayers = this->layeredCS->giveNumberOfLayers();  
    FloatArray f, N;
    FloatArray genEps;
    FloatMatrix B;

    for ( int layer = 1; layer <= numberOfLayers; layer++ ) {
        Material *mat = domain->giveMaterial( this->layeredCS->giveLayerMaterial(layer) );

        for (GaussPoint *gp : *integrationRulesArray [ layer - 1 ]) {
            const FloatArray &lCoords = gp->giveNaturalCoordinates();
            this->computeBmatrixAt(lCoords, B);
            genEps.beProductOf(B, solVec);

            double zeta = giveGlobalZcoord(lCoords);
            this->computeSectionalForcesAt(N, gp, mat, tStep, genEps, zeta); // these are per unit volume
            
            double dV = this->computeVolumeAroundLayer(gp, layer);
            f.plusProduct(B, N, dV);
        }
    }

    answer.resize( ndofs );
    answer.zero();
    answer.assemble(f, this->giveOrderingDofTypes());
}


void
Shell7Base :: computeSectionalForcesAt(FloatArray &sectionalForces, IntegrationPoint *ip, Material *mat, TimeStep *tStep, FloatArray &genEps, double zeta)
{
    // New, in terms of PK1 stress
    // \Lambda_i * P * G^I
    FloatArray PG1(3), PG2(3), PG3(3);
    FloatMatrix lambda[3], Gcon, P, PG;
    this->computeStressMatrix(P, genEps, ip, mat, tStep);

    this->evalInitialContravarBaseVectorsAt(ip->giveNaturalCoordinates(), Gcon);
    PG.beProductOf(P,Gcon);
    PG1.beColumnOf(PG, 1);
    PG2.beColumnOf(PG, 2);
    PG3.beColumnOf(PG, 3);
    this->computeLambdaGMatrices(lambda, genEps, zeta); // associated with the variation of the test functions   

    // f = lambda_1^T * P*G^1 + lambda_2^T * P*G^2 + lambda_3^T * P*G^3
    sectionalForces.clear();
    sectionalForces.plusProduct(lambda[0], PG1, 1.0);
    sectionalForces.plusProduct(lambda[1], PG2, 1.0);
    sectionalForces.plusProduct(lambda[2], PG3, 1.0);
}

#endif




// Mass matrices

#if 1

void
Shell7Base :: computeThicknessMappingCoeff(GaussPoint *gp, FloatArray &answer)
{
    //thickness jacobian = ratio between volume and area: j0 = a3 + a2*zeta^2 + a1 * zeta
    // Returns array with a1-a3, used in expression for analytical integration of mass matrix.
    const FloatArray &lcoords = gp->giveNaturalCoordinates();

    FloatMatrix dNdxi;
    this->fei->evaldNdxi( dNdxi, lcoords, FEIElementGeometryWrapper(this) );

    FloatArray M, dM1(3), dM2(3), dX1(3), dX2(3);
    double gam, dg1, dg2;
    FloatMatrix B;
    this->computeBmatrixAt(lcoords, B);

    FloatArray initSolVec, genEps;
    initSolVec = this->giveInitialSolutionVector();
    genEps.beProductOf(B, initSolVec);
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
    int nPointsTri = 6; //todo generalize
    GaussIntegrationRule iRule(1, this);
    iRule.SetUpPointsOnWedge(nPointsTri, 1, _3dMat); ///@todo replce with triangle which has a xi3-coord
    
    
    //------------------------------
    FloatMatrix N, Ntm, NtmN, temp;
    FloatArray solVec;
    this->giveUpdatedSolutionVector(solVec, tStep);

    ///@todo: Should be changed to integration over the layers and use the corresponding material
    Material *mat = domain->giveMaterial( this->layeredCS->giveLayerMaterial(1) );     // for now, while I don't have an analytical exp.

    for ( auto &gp: iRule ) {
        const FloatArray &lCoords = gp->giveNaturalCoordinates();
        this->computeNmatrixAt(lCoords, N);
        FloatArray unknowns, m(3);
        unknowns.beProductOf(N, solVec);        // [x, m, gam]^T
        m = { unknowns.at(4), unknowns.at(5), unknowns.at(6) };
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
        double xi = 0.0;
        double dA = this->computeAreaAround(gp, xi);
        Ntm.beTProductOf(N, mass);
        NtmN.beProductOf(Ntm, N);         // switch to sym prod upper something
        NtmN.times(dA * rho);
        temp.add(NtmN);
    }

    int ndofs = this->computeNumberOfDofs();
    answer.resize(ndofs, ndofs);
    answer.zero();
    const IntArray &ordering = this->giveOrderingDofTypes();    
    answer.assemble(temp, ordering);
    answer.symmetrized();
}


void
Shell7Base :: giveMassFactorsAt(GaussPoint *gp, FloatArray &factors, double &gam)
{
    FloatArray coeff;
    this->computeThicknessMappingCoeff(gp, coeff);
    double a1 = coeff.at(1);
    double a2 = coeff.at(2);
    double a3 = coeff.at(3);

    double h  = this->giveCrossSection()->give(CS_Thickness,NULL);
    double h2 = h * h;
    double h3 = h2 * h;
    double h5 = h2 * h3;
    double gam2 = gam * gam;

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
    // For analytically integrated throught the thickness, see computeMassMatrix

    FloatMatrix mass;
    FloatArray solVec;
    this->giveUpdatedSolutionVector(solVec, tStep);
    int numberOfLayers = this->layeredCS->giveNumberOfLayers();     // conversion of data

    FloatMatrix M(42,42);

    for ( int layer = 1; layer <= numberOfLayers; layer++ ) {
        Material *mat = domain->giveMaterial( this->layeredCS->giveLayerMaterial(layer) );

        /* Consistent mass matrix M = int{N^t*mass*N}
         *
         *         3    3    1
         *         3 [a*I  b*I   c*m      [A  B  C
         * mass =  3       d*I   e*m    =     D  E
         *         1  sym       f*m.m]     sym   F]
         */

        for ( auto &gp: *integrationRulesArray [ layer - 1 ] ) {
            const FloatArray &lCoords = gp->giveNaturalCoordinates();
            FloatMatrix lambda, B, N, temp;
            FloatArray genEps;
            this->computeBmatrixAt(lCoords, B);
            genEps.beProductOf(B, solVec);    
            double zeta = giveGlobalZcoord(gp->giveNaturalCoordinates());
            this->computeLambdaNMatrix(lambda, genEps, zeta);
            
            // could also create lambda*N and then plusProdSymm - probably faster
            mass.beTProductOf(lambda,lambda);
            this->computeNmatrixAt(lCoords, N);
            temp.beProductOf(mass,N);
        
            double dV = this->computeVolumeAroundLayer(gp, layer);
            double rho = mat->give('d', gp);
            M.plusProductSymmUpper(N, temp, rho*dV);
        }
        M.symmetrized();
        const IntArray &ordering = this->giveOrderingDofTypes();
        answer.zero();
        answer.assemble(M, ordering, ordering);
#endif

    }
}


void
Shell7Base :: computeConvectiveMassForce(FloatArray &answer, TimeStep *tStep)
{
    ///@todo: very old version should be checked
    // Analytically integrated over the thickness. Constant density assumed.
    int nPointsTri = 6; //todo generalize
    GaussIntegrationRule iRule(1, this);
    iRule.SetUpPointsOnWedge(nPointsTri, 1, _3dMat); //@todo replce with triangle which has a xi3-coord
        
    FloatMatrix N;
    FloatArray a, da, m, dm, aVec, daVec, fm(7);
    double gam, dgam, dA;

    answer.resize(42);
    answer.zero();
    for ( auto &gp: iRule ) { // rule 2 for mid-plane integration only
        const FloatArray &lCoords = gp->giveNaturalCoordinates();
        this->computeNmatrixAt(lCoords, N);
        this->giveUpdatedSolutionVector(aVec, tStep);
        // Fix for the new numbering in B & N
        
        a.beProductOf(N, aVec);        // [ x,  m,  gam]^T
        da.beProductOf(N, daVec);      // [dx, dm, dgam]^T
        m = { a.at(4), a.at(5), a.at(6) };
        gam =  a.at(7);
        dm = { da.at(4), da.at(5), da.at(6) };
        dgam = da.at(7);

        double a1, a2, a3, h, h2, h3, h5, fac1, fac2, fac3, rho;
        FloatArray coeff;

        Material *mat = domain->giveMaterial( this->layeredCS->giveLayerMaterial(1) ); ///@todo fix this method
        rho = mat->give('d', gp);
        h   = this->giveCrossSection()->give(CS_Thickness, gp);
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
        double xi = 0.0;
        dA = this->computeAreaAround(gp, xi);
        answer.plusProduct(N, fm, dA);
    }
}


// External forces

#if 1
void
Shell7Base :: computeEdgeLoadVectorAt(FloatArray &answer, Load *load, int iEdge, TimeStep *tStep, ValueModeType mode)
{
    BoundaryLoad *edgeLoad = dynamic_cast< BoundaryLoad * >( load );
    if ( edgeLoad ) {
        this->computeTractionForce(answer, iEdge, edgeLoad, tStep, mode);
        return;
    } else {
        OOFEM_ERROR("Load type not supported");
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
        this->computePressureForce(force, solVec, iSurf, surfLoad, tStep, mode);
        IntArray mask;
        this->giveSurfaceDofMapping(mask, 1);         // same dofs regardless of iSurf->1
        answer.resize( this->computeNumberOfDofs() );
        answer.zero();
        answer.assemble(force, mask);

        return;
    } else {
        OOFEM_ERROR("Load type not supported");
        return;
    }
}

void
Shell7Base :: computePressureForce(FloatArray &answer, FloatArray solVec, const int iSurf, BoundaryLoad *surfLoad, TimeStep *tStep, ValueModeType mode)
{
    // Computes pressure loading. Acts normal to the current (deformed) surface.
    ConstantPressureLoad* pLoad = dynamic_cast< ConstantPressureLoad * >( surfLoad );
    int nPointsTri = 6; //todo generalize
    IntegrationRule *iRule = new GaussIntegrationRule(1, this);
    iRule->SetUpPointsOnWedge(nPointsTri, 1, _3dMat); //@todo replce with triangle which has a xi3-coord
    
    
    FloatMatrix N, B, lambda;
    FloatArray Fp, fp, genEps, genEpsC, lCoords(3), traction, solVecC;

    for ( auto *ip: *iRule ) { // rule #2 for surface integration
        lCoords.at(1) = ip->giveNaturalCoordinate(1);
        lCoords.at(2) = ip->giveNaturalCoordinate(2);
        lCoords.at(1) = pLoad->giveLoadOffset( );
        double zeta = giveGlobalZcoord( lCoords );
        this->computeBmatrixAt(lCoords, B);
        this->computeNmatrixAt(lCoords, N);
        this->giveUpdatedSolutionVector(solVecC, tStep);
        genEpsC.beProductOf(B, solVecC);       
        this->computePressureForceAt(ip, traction, iSurf, genEpsC, surfLoad, tStep, mode);
      
        genEps.beProductOf(B, solVec);
        this->computeLambdaNMatrix(lambda, genEps, zeta);
        fp.beTProductOf(lambda,traction);

        double xi = pLoad->giveLoadOffset();
        double dA = this->computeAreaAround(ip, xi);

        Fp.plusProduct(N, fp, dA);
    }
    answer.resize( Shell7Base :: giveNumberOfDofs() );
    answer.zero();
    answer.assemble(Fp, this->giveOrderingDofTypes());
}


void
Shell7Base :: computePressureForceAt(GaussPoint *gp, FloatArray &traction, const int iSurf, FloatArray genEps, BoundaryLoad *surfLoad, TimeStep *tStep, ValueModeType mode)
{
    // Computes pressure loading. Acts normal to the current (deformed) surface.
    //@todo traction load should be moved outside this method
    if ( iSurf != 0 ) {
        OOFEM_ERROR("incompatible load surface must be 0 for this element");
    }

    FloatArray load;
    if ( ConstantPressureLoad* pLoad = dynamic_cast< ConstantPressureLoad * >( surfLoad ) ) {
        FloatMatrix gcov; 
        FloatArray g1, g2;
        FloatArray lcoords(3);
        lcoords.at(1) = gp->giveNaturalCoordinate(1);
        lcoords.at(2) = gp->giveNaturalCoordinate(2);
        lcoords.at(3) = pLoad->giveLoadOffset();

        this->evalCovarBaseVectorsAt(lcoords, gcov, genEps, tStep); 
        g1.beColumnOf(gcov,1);
        g2.beColumnOf(gcov,2);
        surfLoad->computeValueAt(load, tStep, lcoords, mode);        // pressure component
        traction.beVectorProductOf(g1, g2);        // normal vector (should not be normalized due to integraton in reference config.)
        traction.times( -load.at(1) );
    } else if ( dynamic_cast< ConstantSurfaceLoad * >( surfLoad ) ) {
        surfLoad->computeValueAt(traction, tStep, gp->giveNaturalCoordinates(), mode);        // traction vector
    } else {
        OOFEM_ERROR("incompatible load type");
    }

}

void
Shell7Base :: evalCovarNormalAt(FloatArray &nCov, const FloatArray &lCoords, FloatArray &genEpsC, TimeStep *tStep)
{
    FloatMatrix gcov;
    this->evalCovarBaseVectorsAt(lCoords, gcov, genEpsC, tStep);
    FloatArray g1, g2;
    g1.beColumnOf(gcov,1);
    g2.beColumnOf(gcov,2);
    nCov.beVectorProductOf(g1, g2);
    nCov.normalize();
}

void
Shell7Base :: evalInitialCovarNormalAt(FloatArray &nCov, const FloatArray &lCoords)
{
    FloatMatrix Gcov;
    this->evalInitialCovarBaseVectorsAt(lCoords, Gcov);
    FloatArray G1, G2;
    G1.beColumnOf(Gcov,1);
    G2.beColumnOf(Gcov,2);
    nCov.beVectorProductOf(G1, G2);
    nCov.normalize();
}

void
Shell7Base :: computeTractionForce(FloatArray &answer, const int iEdge, BoundaryLoad *edgeLoad, TimeStep *tStep, ValueModeType mode)
{

    int approxOrder = edgeLoad->giveApproxOrder() + this->giveInterpolation()->giveInterpolationOrder();
    int numberOfGaussPoints = ( int ) ceil( ( approxOrder + 1. ) / 2. );
    GaussIntegrationRule iRule(1, this, 1, 1);
    iRule.SetUpPointsOnLine(numberOfGaussPoints, _Unknown);    
    FloatMatrix N, Q;
    FloatArray fT(7), Nf, components;
    
    Load :: CoordSystType coordSystType = edgeLoad->giveCoordSystMode();

    for ( GaussPoint *gp : iRule ) {    
        const FloatArray &lCoords = gp->giveNaturalCoordinates();
        edgeLoad->computeValueAt(components, tStep, lCoords, mode);
        this->edgeComputeNmatrixAt(lCoords, N);

        //if ( coordSystType ==  BoundaryLoad :: BL_UpdatedGlobalMode ) {
        if ( coordSystType ==  Load :: CST_UpdatedGlobal ) {
            
            // Updated global coord system
            FloatMatrix gcov;
            this->edgeEvalCovarBaseVectorsAt(lCoords, iEdge, gcov, tStep); 
            Q.beTranspositionOf(gcov);

            FloatArray distrForces(3), distrMoments(3), t1, t2;
            distrForces = { components.at(1), components.at(2), components.at(3)};
            distrMoments = { components.at(4), components.at(5), components.at(6) };
            t1.beTProductOf(Q, distrForces);
            t2.beTProductOf(Q, distrMoments);
            fT.addSubVector(t1,1);
            fT.addSubVector(t2,4);
            fT.at(7) = components.at(7); // don't do anything with the 'gamma'-load

        } else if( coordSystType == Load :: CST_Global ) { 
            // Undeformed global coord system
            for ( int i = 1; i <= 7; i++) {
                fT.at(i) = components.at(i);
            }
        } else {
            OOFEM_ERROR("Shell7Base :: computeTractionForce - does not support local coordinate system");
        }

        double dL = this->edgeComputeLengthAround(gp, iEdge);        
        
        Nf.plusProduct(N, fT, dL);
    }

    IntArray mask;
    this->giveEdgeDofMapping(mask, iEdge);
    answer.resize( Shell7Base :: giveNumberOfDofs()  );
    answer.zero();
    answer.assemble(Nf, mask);

}


void
Shell7Base :: computeBodyLoadVectorAt(FloatArray &answer, Load *forLoad, TimeStep *tStep, ValueModeType mode)
{
    OOFEM_ERROR("Shell7Base :: computeBodyLoadVectorAt - currently not implemented");
}

#endif



// Integration
#if 1

//@todo should be moved to tr-elements
double
Shell7Base :: edgeComputeLengthAround(GaussPoint *gp, const int iedge)
{
    FloatArray G1, G3;
    double detJ;
    const FloatArray &lcoords = gp->giveNaturalCoordinates();
    this->edgeEvalInitialCovarBaseVectorsAt(lcoords, iedge, G1, G3);
    detJ = G1.computeNorm();
    return detJ * gp->giveWeight();
}
#endif



// Recovery of nodal values

#if 1

void Shell7Base :: NodalAveragingRecoveryMI_computeSideValue(FloatArray &answer, int side, InternalStateType type, TimeStep *tStep)
{
    // what is this?
    answer.resize(0);
}

void Shell7Base :: NodalAveragingRecoveryMI_computeNodalValue(FloatArray &answer, int node, InternalStateType type, TimeStep *tStep)
{
    if ( type == IST_DirectorField ) {
        answer.resize(3);
        answer = this->giveInitialNodeDirector(node);
        answer.at(1) += this->giveNode(node)->giveDofWithID(W_u)->giveUnknown(VM_Total, tStep);
        answer.at(2) += this->giveNode(node)->giveDofWithID(W_v)->giveUnknown(VM_Total, tStep);
        answer.at(3) += this->giveNode(node)->giveDofWithID(W_w)->giveUnknown(VM_Total, tStep);
        answer.times( this->giveCrossSection()->give(CS_Thickness, NULL) );
    } else {
        answer.resize(0);
    }
}





void
Shell7Base :: NodalRecoveryMI_computeNValProduct(FloatMatrix &answer, int layer, InternalStateType type,
                                                                      TimeStep *tStep)
{  // evaluates N^T sigma over element volume
   // N(nsigma, nsigma*nnodes)
   // Definition : sigmaVector = N * nodalSigmaVector
    FloatArray stressVector, n;
    Element *elem = this->ZZNodalRecoveryMI_giveElement();

    //int size = ZZNodalRecoveryMI_giveDofManRecordSize(type);
    int size = 6; ///@todo this is hard coded for stress recovery

    answer.zero();
    for ( auto *gp: *integrationRulesArray [ layer - 1 ] ) {
        double dV = this->computeVolumeAroundLayer(gp, layer);

        if ( !elem->giveIPValue(stressVector, gp, type, tStep) ) {
            stressVector.resize(size);
            stressVector.zero();
        }

        this->ZZNodalRecoveryMI_ComputeEstimatedInterpolationMtrx(n, gp, type);
        answer.plusDyadUnsym(n, stressVector, dV);

    }
}

void
Shell7Base :: NodalRecoveryMI_computeNNMatrix(FloatArray &answer, int layer, InternalStateType type)
{
    //
    // Returns NTN matrix (lumped) for Zienkiewicz-Zhu
    // The size of N mtrx is (nstresses, nnodes*nstreses)
    // Definition : sigmaVector = N * nodalSigmaVector
    //
    FloatMatrix fullAnswer;
    FloatArray n;

    for ( auto &gp: *integrationRulesArray [ layer - 1 ] ) {
        double dV = this->computeVolumeAroundLayer(gp, layer);
        this->ZZNodalRecoveryMI_ComputeEstimatedInterpolationMtrx(n, gp, type);
        fullAnswer.plusDyadSymmUpper(n, dV);
    }
    fullAnswer.symmetrized();

    ///@todo Make a "lumpMatrix" function in FloatArray ? (summing over columns would also be faster)
    answer.resize(fullAnswer.giveNumberOfRows());
    for ( int i = 1; i <= fullAnswer.giveNumberOfRows(); i++ ) {
        double sum = 0.0;
        for ( int j = 1; j <= fullAnswer.giveNumberOfColumns(); j++ ) {
            sum += fullAnswer.at(i, j);
        }
        answer.at(i) = sum;
    }
}

void
Shell7Base :: ZZNodalRecoveryMI_ComputeEstimatedInterpolationMtrx(FloatArray &answer, GaussPoint *gp, InternalStateType type)
{
    // evaluates N matrix (interpolation estimated stress matrix)
    // according to Zienkiewicz & Zhu paper
    // N(nsigma, nsigma*nnodes)
    // Definition : sigmaVector = N * nodalSigmaVector
    FEInterpolation *interpol = static_cast< FEInterpolation * >( &this->interpolationForExport );

    // test if underlying element provides interpolation
    if ( interpol ) {
        ///@todo fix this whole compostie recovery thing in a better way

        interpol->evalN( answer, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );
    } else {
        // ok default implementation can not work, as element is not providing valid interpolation
        // to resolve this, one can overload this method for element implementing ZZNodalRecoveryModelInterface
        // or element should provide interpolation.
        OOFEM_ERROR("Element %d not providing valid interpolation", this->giveNumber() );
    }
}


void 
Shell7Base :: NodalRecoveryMI_recoverValues(std::vector<FloatArray> &recoveredValues, int layer, InternalStateType type, TimeStep *tStep)
{
    // ZZ recovery
    FloatArray nnMatrix;
    FloatMatrix nValProd;
    NodalRecoveryMI_computeNValProduct(nValProd, layer, type, tStep);
    NodalRecoveryMI_computeNNMatrix(nnMatrix, layer, type);
    int recoveredSize = nValProd.giveNumberOfColumns();
    int numNodes = nValProd.giveNumberOfRows();
    recoveredValues.resize(numNodes);
    
    for ( int i = 1; i <= numNodes; i++ ) {
        FloatArray temp(6);
        recoveredValues[i-1].resize(9);
        
        for ( int j = 1; j <= recoveredSize; j++ ) {
            temp.at(j) = nValProd.at(i,j)/nnMatrix.at(i);
        }
        
        recoveredValues[i-1] = convV6ToV9Stress(temp);
    }



}

#endif




// Computation of solution vectors

#if 1

void
Shell7Base :: temp_computeBoundaryVectorOf(IntArray &dofIdArray, int boundary, ValueModeType u, TimeStep *tStep, FloatArray &answer)
{
    ///@todo: NOT CHECKED!!!
    // Routine to extract vector given an array of dofid items
    // If a certain dofId does not exist a zero is used as value
    
    IntArray bNodes;
    this->fei->computeLocalEdgeMapping(bNodes, boundary);
    this->computeBoundaryVectorOf(bNodes, dofIdArray, u, tStep, answer); ///@todo uses new standard method

    //answer.resize( dofIdArray.giveSize() * bNodes.giveSize() );
    //answer.zero();
    //int k = 0;
    //for ( int i = 1; i <= bNodes.giveSize(); i++ ) {
    //    DofManager *dMan = this->giveDofManager(bNodes.at(i));
    //    for (int j = 1; j <= dofIdArray.giveSize(); j++ ) {
    //        Dof *d = dMan->giveDof(j);
    //        k++;
    //        if ( dMan->hasDofID( (DofIDItem) dofIdArray.at(j) ) ) {
    //            answer.at(k) = d->giveUnknown(VM_Total, tStep); 
    //        }
    //    }
    //}


}

void
Shell7Base :: giveUpdatedSolutionVector(FloatArray &answer, TimeStep *tStep)
{
    // Computes updated solution as: x = X + dx, m = M + dM, gam = 0 + dgam
    // This is because the element formulation is in terms of placement and not displacement.
    answer = this->giveInitialSolutionVector();
    FloatArray temp;
    int dummy = 0;
    IntArray dofIdArray;
    Shell7Base::giveDofManDofIDMask(dummy, dofIdArray);
    Element :: computeVectorOf(dofIdArray, VM_Total, tStep, temp, true);
    answer.assemble( temp, this->giveOrderingNodes() );
}



void
Shell7Base :: setupInitialSolutionVector() 
{
    // Gives the initial values of X, M and gamma which definies the initial position
    this->initialSolutionVector.resize( Shell7Base :: giveNumberOfDofs() );
    this->initialSolutionVector.zero();
    int ndofs_xm  = 3 * this->giveNumberOfDofManagers();
    
    // Reference position and directors
    for ( int i = 1, j = 0; i <= this->giveNumberOfDofManagers(); i++, j += 3 ) {
        FloatArray *Xi = this->giveNode(i)->giveCoordinates();
        FloatArray Mi  = this->giveInitialNodeDirector(i);
        this->initialSolutionVector.at(1 + j) = Xi->at(1);
        this->initialSolutionVector.at(2 + j) = Xi->at(2);
        this->initialSolutionVector.at(3 + j) = Xi->at(3);
        this->initialSolutionVector.at(ndofs_xm + 1 + j) = Mi.at(1);
        this->initialSolutionVector.at(ndofs_xm + 2 + j) = Mi.at(2);
        this->initialSolutionVector.at(ndofs_xm + 3 + j) = Mi.at(3);
        // Assumes gam=0 at t=0
    }
}


void
Shell7Base :: edgeGiveUpdatedSolutionVector(FloatArray &answer, const int iEdge, TimeStep *tStep)
{
    //this->edgeGiveInitialSolutionVector(answer, iedge);
    this->giveInitialEdgeSolutionVector(iEdge);
    FloatArray temp;
    int dummy = 0;
    IntArray dofIdArray;
    Shell7Base :: giveDofManDofIDMask(dummy, dofIdArray);
    this->temp_computeBoundaryVectorOf(dofIdArray, iEdge, VM_Total, tStep, temp);
    //answer.assemble( temp, this->giveOrdering(EdgeInv) );
    answer.assemble( temp, this->giveOrderingEdgeNodes());
}

void
Shell7Base :: setupInitialEdgeSolutionVector()
{
    //int numEdges = this->giveNumberOfBoundarySides();///TODO fix
    int numEdges = 3;
    this->initialEdgeSolutionVectors.resize( numEdges );
    for ( int iEdge = 1; iEdge <= numEdges; iEdge++ ) {
        FloatArray &solVec = this->initialEdgeSolutionVectors[iEdge-1];
        solVec.resize( this->giveNumberOfEdgeDofs() );
        solVec.zero();
        IntArray edgeNodes;
        this->fei->computeLocalEdgeMapping(edgeNodes, iEdge);
        int ndofs_x = 3 * edgeNodes.giveSize();
        for ( int i = 1, j = 0; i <= edgeNodes.giveSize(); i++, j += 3 ) {
            FloatArray *Xi = this->giveNode( edgeNodes.at(i) )->giveCoordinates();
            FloatArray Mi  = this->giveInitialNodeDirector( edgeNodes.at(i) );
            solVec.at(1 + j) = Xi->at(1);
            solVec.at(2 + j) = Xi->at(2);
            solVec.at(3 + j) = Xi->at(3);
            solVec.at(ndofs_x + 1 + j) = Mi.at(1);
            solVec.at(ndofs_x + 2 + j) = Mi.at(2);
            solVec.at(ndofs_x + 3 + j) = Mi.at(3);
            // gam(t=0)=0 is assumed
        }
    }
}

void
Shell7Base :: giveGeneralizedStrainComponents(FloatArray genEps, FloatArray &dphidxi1, FloatArray &dphidxi2, FloatArray &dmdxi1,
                                              FloatArray &dmdxi2, FloatArray &m, double &dgamdxi1, double &dgamdxi2, double &gam) {
    // generealized strain vector  [dxdxi, dmdxi, m, dgamdxi, gam]^T
    dphidxi1 = { genEps.at(1), genEps.at(2), genEps.at(3) };
    dphidxi2 = { genEps.at(4), genEps.at(5), genEps.at(6) };
    dmdxi1 = { genEps.at(7), genEps.at(8), genEps.at(9) };
    dmdxi2 = { genEps.at(10), genEps.at(11), genEps.at(12) };
    m = { genEps.at(13), genEps.at(14), genEps.at(15) };
    dgamdxi1 = genEps.at(16);
    dgamdxi2 = genEps.at(17);
         gam = genEps.at(18);
}


void
Shell7Base :: giveUnknownsAt(const FloatArray &lCoords, FloatArray &solVec, FloatArray &x, FloatArray &m, double &gam, TimeStep *tStep)
{
    // returns the unknowns evaluated at a point (xi1, xi2, xi3)
    FloatMatrix N;
    this->computeNmatrixAt(lCoords, N);
    FloatArray temp;
    temp.beProductOf(N, solVec);
    x = { temp.at(1), temp.at(2), temp.at(3) };
    m = { temp.at(4), temp.at(5), temp.at(6) };
    gam = temp.at(7);
}


#endif


// N and B matrices

#if 1


void
Shell7Base :: edgeComputeNmatrixAt(const FloatArray &lcoords, FloatMatrix &answer)
{
// Returns the displacement interpolation matrix {N} of the receiver 
// evaluated at gaussPoint along one edge.

    answer.resize( 7, this->giveNumberOfEdgeDofs() );
    answer.zero();

    FloatArray N;
    this->fei->edgeEvalN( N, 1, lcoords, FEIElementGeometryWrapper(this) );

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
Shell7Base :: edgeComputeBmatrixAt(const FloatArray &lcoords, FloatMatrix &answer, int li, int ui)
{
/* Returns the  matrix {B} of the receiver, evaluated at gp. Such that
 * B*a = [dxbar_dxi, dwdxi, w, dgamdxi, gam]^T, where a is the vector of unknowns
 */

    answer.resize( 11, this->giveNumberOfEdgeDofs() );
    answer.zero();
    FloatArray N, dNdxi;

    this->fei->edgeEvalN( N, 1, lcoords, FEIElementGeometryWrapper(this) );
    int iedge = 0;
    this->fei->edgeEvaldNdxi( dNdxi, iedge, lcoords, FEIElementGeometryWrapper(this) );

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
Shell7Base :: computeBmatrixAt(const FloatArray &lcoords, FloatMatrix &answer, int li, int ui)
{
    // Returns the  matrix {B} of the receiver, evaluated at gp. Such that
    // B*a = [dxbar_dxi, dwdxi, w, dgamdxi, gam]^T, where a is the vector of unknowns
 
    int ndofs = Shell7Base :: giveNumberOfDofs();
    int ndofs_xm  = 3 * this->giveNumberOfDofManagers();
    answer.resize(18, ndofs);
    answer.zero();
    FloatArray N;
    FloatMatrix dNdxi;
    this->fei->evalN( N, lcoords, FEIElementGeometryWrapper(this) );
    this->fei->evaldNdxi( dNdxi, lcoords, FEIElementGeometryWrapper(this) );

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
Shell7Base :: computeNmatrixAt(const FloatArray &iLocCoord, FloatMatrix &answer)
{
    // Returns the displacement interpolation matrix {N} of the receiver,
    // evaluated at gp.

    int ndofs = Shell7Base :: giveNumberOfDofs();
    int ndofs_xm  = 3 * this->giveNumberOfDofManagers();
    answer.resize(7, ndofs);
    answer.zero();
    FloatArray N;
    this->fei->evalN( N, iLocCoord, FEIElementGeometryWrapper(this) );
    
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


#endif


// VTK export
#if 1
void
Shell7Base :: vtkEvalInitialGlobalCoordinateAt(const FloatArray &localCoords, int layer, FloatArray &globalCoords)
{
    double zeta = giveGlobalZcoordInLayer(localCoords.at(3), layer);
    FloatArray N;
    this->fei->evalN( N, localCoords, FEIElementGeometryWrapper(this) );

    globalCoords.clear();
    for ( int i = 1; i <= this->giveNumberOfDofManagers(); i++ ) {
        FloatArray &xbar = *this->giveNode(i)->giveCoordinates();
        const FloatArray &M = this->giveInitialNodeDirector(i);
        globalCoords.add(N.at(i), ( xbar + zeta * M ));
    }

}

void
Shell7Base :: vtkEvalInitialGlobalCZCoordinateAt(const FloatArray &localCoords, int interface, FloatArray &globalCoords)
{
    double zeta = giveGlobalZcoordInLayer(1.0, interface);
    FloatArray N;
    this->fei->evalN( N, localCoords, FEIElementGeometryWrapper(this) );

    globalCoords.clear();
    for ( int i = 1; i <= this->giveNumberOfDofManagers(); i++ ) {
        FloatArray &xbar = *this->giveNode(i)->giveCoordinates();
        const FloatArray &M = this->giveInitialNodeDirector(i);
        globalCoords.add(N.at(i), ( xbar + zeta * M ));
    }

}

void
Shell7Base :: vtkEvalUpdatedGlobalCoordinateAt(const FloatArray &localCoords, int layer, FloatArray &globalCoords, TimeStep *tStep)
{
    FloatArray solVec;
    this->giveUpdatedSolutionVector(solVec, tStep);
    FloatArray x, m; double gam=0;
    this->giveUnknownsAt(localCoords, solVec, x, m, gam, tStep); 
    double zeta = giveGlobalZcoordInLayer(localCoords.at(3), layer);
    double fac = ( zeta + 0.5 * gam * zeta * zeta );
    globalCoords = x + fac*m;
}

void
Shell7Base :: giveCompositeExportData(std::vector< VTKPiece > &vtkPieces, IntArray &primaryVarsToExport, IntArray &internalVarsToExport, IntArray cellVarsToExport, TimeStep *tStep )
{
    vtkPieces.resize(1);
    this->giveShellExportData(vtkPieces[0], primaryVarsToExport, internalVarsToExport, cellVarsToExport, tStep );
    //this->giveCZExportData(vtkPieces[1], primaryVarsToExport, internalVarsToExport, cellVarsToExport, tStep );
    
}

void 
Shell7Base :: giveShellExportData(VTKPiece &vtkPiece, IntArray &primaryVarsToExport, IntArray &internalVarsToExport, IntArray cellVarsToExport, TimeStep *tStep )            
{   

    int numCells = this->layeredCS->giveNumberOfLayers();
    const int numCellNodes  = 15; // quadratic wedge
    int numTotalNodes = numCellNodes*numCells;

    vtkPiece.setNumberOfCells(numCells);
    vtkPiece.setNumberOfNodes(numTotalNodes);

    std::vector <FloatArray> nodeCoords;
    int val    = 1;
    int offset = 0;
    IntArray nodes(numCellNodes);

    // Compute fictious node coords
    int nodeNum = 1;
    for ( int layer = 1; layer <= numCells; layer++ ) {
        

        // Node coordinates
        this->giveFictiousNodeCoordsForExport(nodeCoords, layer);       
        
        for ( int node = 1; node <= numCellNodes; node++ ) {    
            vtkPiece.setNodeCoords(nodeNum, nodeCoords[node-1] );
            nodeNum += 1;
        }

        // Connectivity       
        for ( int i = 1; i <= numCellNodes; i++ ) {            
            nodes.at(i) = val++;
        }
        vtkPiece.setConnectivity(layer, nodes);
        
        // Offset
        offset += numCellNodes;
        vtkPiece.setOffset(layer, offset);

        // Cell types
        vtkPiece.setCellType(layer, 26); // Quadratic wedge
    }


    // Export nodal variables from primary fields        
    vtkPiece.setNumberOfPrimaryVarsToExport(primaryVarsToExport.giveSize(), numTotalNodes);

    std::vector<FloatArray> updatedNodeCoords;
    FloatArray u(3);
    std::vector<FloatArray> values;
    for ( int fieldNum = 1; fieldNum <= primaryVarsToExport.giveSize(); fieldNum++ ) {
        UnknownType type = ( UnknownType ) primaryVarsToExport.at(fieldNum);
        nodeNum = 1;
        for ( int layer = 1; layer <= numCells; layer++ ) {            
            
            if ( type == DisplacementVector ) { // compute displacement as u = x - X
                this->giveFictiousNodeCoordsForExport(nodeCoords, layer);
                this->giveFictiousUpdatedNodeCoordsForExport(updatedNodeCoords, layer, tStep);
                for ( int j = 1; j <= numCellNodes; j++ ) {
                    u = updatedNodeCoords[j-1];
                    u.subtract(nodeCoords[j-1]);
                    vtkPiece.setPrimaryVarInNode(fieldNum, nodeNum, u);
                    nodeNum += 1;        
                }

            } else {
                NodalRecoveryMI_recoverValues(values, layer, ( InternalStateType ) 1, tStep); // does not work well - fix
                for ( int j = 1; j <= numCellNodes; j++ ) {
                    vtkPiece.setPrimaryVarInNode(fieldNum, nodeNum, values[j-1]);
                    nodeNum += 1;
                }
            }
        }
    }

    // Export nodal variables from internal fields
    
    vtkPiece.setNumberOfInternalVarsToExport( internalVarsToExport.giveSize(), numTotalNodes );
    for ( int fieldNum = 1; fieldNum <= internalVarsToExport.giveSize(); fieldNum++ ) {
        InternalStateType type = ( InternalStateType ) internalVarsToExport.at(fieldNum);
        nodeNum = 1;
        //this->recoverShearStress(tStep);
        for ( int layer = 1; layer <= numCells; layer++ ) {            
            recoverValuesFromIP(values, layer, type, tStep);        
            for ( int j = 1; j <= numCellNodes; j++ ) {
                vtkPiece.setInternalVarInNode( fieldNum, nodeNum, values[j-1] );
                //ZZNodalRecoveryMI_recoverValues(el.nodeVars[fieldNum], layer, type, tStep);          
                nodeNum += 1;        
            }                                
        }  
    }


    // Export cell variables
    FloatArray average;
    vtkPiece.setNumberOfCellVarsToExport(cellVarsToExport.giveSize(), numCells);
    for ( int i = 1; i <= cellVarsToExport.giveSize(); i++ ) {
        InternalStateType type = ( InternalStateType ) cellVarsToExport.at(i);;
      
        for ( int layer = 1; layer <= numCells; layer++ ) {     
            std :: unique_ptr< IntegrationRule > &iRuleL = integrationRulesArray [ layer - 1 ];
            VTKXMLExportModule::computeIPAverage(average, iRuleL.get(), this, type, tStep);
            
            if ( average.giveSize() == 6 ) {
                vtkPiece.setCellVar(i, layer, convV6ToV9Stress(average) );
            } else {
                vtkPiece.setCellVar(i, layer, average );
            }

        }

    }




}


void 
Shell7Base :: recoverValuesFromIP(std::vector<FloatArray> &recoveredValues, int layer, InternalStateType type, TimeStep *tStep)
{
    // recover nodal values by coosing the ip closest to the node

    //FEInterpolation *interpol = static_cast< FEInterpolation * >( &this->interpolationForExport );

    // composite element interpolator
    FloatMatrix localNodeCoords;
    this->interpolationForExport.giveLocalNodeCoords(localNodeCoords);

    int numNodes = localNodeCoords.giveNumberOfColumns();
    recoveredValues.resize(numNodes);

    // Find closest ip to the nodes
    IntArray closestIPArray(numNodes);
    FloatArray nodeCoords, ipValues;

    for ( int i = 1; i <= numNodes; i++ ) {
        nodeCoords.beColumnOf(localNodeCoords, i);
        double distOld = 3.0; // should not be larger
        std :: unique_ptr< IntegrationRule > &iRule = integrationRulesArray [ layer - 1 ];
        for ( int j = 1; j <= iRule->giveNumberOfIntegrationPoints(); ++j ) {
            IntegrationPoint *ip = iRule->getIntegrationPoint(j);
            const FloatArray &ipCoords = ip->giveNaturalCoordinates();
            double dist = nodeCoords.distance(ipCoords);
            if ( dist < distOld ) {
                closestIPArray.at(i) = j;
                distOld = dist;
            }
        }
    }

   InternalStateValueType valueType =  giveInternalStateValueType(type);

    // recover ip values
    for ( int i = 1; i <= numNodes; i++ ) {
        IntegrationPoint *ip = integrationRulesArray [ layer - 1 ]->getIntegrationPoint( closestIPArray.at(i) );
        this->giveIPValue(ipValues, ip, type, tStep);
        if ( valueType == ISVT_TENSOR_S3 ) {
            recoveredValues[i-1].resize(9);
            recoveredValues[i-1] = convV6ToV9Stress(ipValues);
        //} else if ( ipValues.giveSize() == 0 && type == IST_AbaqusStateVector) {
        //    recoveredValues[i-1].resize(23);
        //    recoveredValues[i-1].zero();
        } else {
            recoveredValues[i-1] = ipValues;
        }
    }

}


void 
Shell7Base :: recoverShearStress(TimeStep *tStep)
{
    // Recover shear stresses at ip by numerical integration of the momentum balance through the thickness
    std::vector<FloatArray> recoveredValues;
    int numberOfLayers = this->layeredCS->giveNumberOfLayers();     // conversion of types
    
    GaussIntegrationRule iRuleThickness(1, this, 1, 1);
    iRuleThickness.SetUpPointsOnLine(this->layeredCS->giveNumIntegrationPointsInLayer(), _Unknown);    
    
    FloatArray dS, Sold;
    FloatMatrix B, Smat(2,6); // 2 stress components * num of in plane ip ///@todo generalize
    Smat.zero();
    FloatArray Tcon(6), Trec(6);  Tcon.zero(); Trec.zero();
    
     for ( int layer = 1; layer <= numberOfLayers; layer++ ) {
        std :: unique_ptr< IntegrationRule > &iRuleL = integrationRulesArray [ layer - 1 ];
        this->recoverValuesFromIP(recoveredValues, layer, IST_StressTensor, tStep);
        //this->ZZNodalRecoveryMI_recoverValues(recoveredValues, layer, IST_StressTensor, tStep);
        double thickness = this->layeredCS->giveLayerThickness(layer);

        //set up vector of stresses in the ip's = [S1_xx, S1_yy, S1_xy, ..., Sn_xx, Sn_yy, Sn_xy]
        int numNodes = 15;
        FloatArray aS(numNodes*3); 
        for ( int j = 1, pos = 0; j <= numNodes; j++, pos+=3 ) {
            aS.at(pos + 1) = recoveredValues[j-1].at(1);   // S_xx
            aS.at(pos + 2) = recoveredValues[j-1].at(2);   // S_yy
            aS.at(pos + 3) = recoveredValues[j-1].at(6);   // S_xy
        }
        int numInPlaneIP = 6;

        for ( int i = 0; i < iRuleThickness.giveNumberOfIntegrationPoints(); i++ ) { 
            double  dz = thickness * iRuleThickness.getIntegrationPoint(i)->giveWeight();
            
            for ( int j = 0; j < numInPlaneIP; j++ ) { 

                int point = i*numInPlaneIP + j; // integration point number
                GaussPoint *gp = iRuleL->getIntegrationPoint(point);

                this->computeBmatrixForStressRecAt(gp->giveNaturalCoordinates(), B, layer);
                dS.beProductOf(B,aS*(-dz)); // stress increment

                StructuralMaterialStatus* status = dynamic_cast< StructuralMaterialStatus* > ( gp->giveMaterialStatus() );
                Sold = status->giveStressVector();
                
                Smat.at(1,j+1) += dS.at(1); // add increment from each level
                Smat.at(2,j+1) += dS.at(2);

                //Tcon.at(j+1) += Sold.at(5)*dz;

                // Replace old stresses with  - this should probably not be done as it may affect the convergence in a nonlinear case
                Sold.at(5) = Smat.at(1,j+1); // S_xz
                Sold.at(4) = Smat.at(2,j+1); // S_yz


                status->letStressVectorBe(Sold);
                //Trec.at(j+1) += Sold.at(5)*dz;
            }
        }


    }

}


void
Shell7Base :: computeBmatrixForStressRecAt(const FloatArray &lcoords, FloatMatrix &answer, int layer)
{
    // Returns the  special matrix {B} of the receiver, evaluated at gp. Such that
    // B*a = [dS_xx/dx + dS_xy/dy, dS_yx/dx + dS_yy/dy ]^T, where a is the vector of in plane 
    // stresses [S_xx, S_yy, S_xy]
 
    // set up virtual cell geometry for an qwedge
    const int numNodes = 15;
    std::vector<FloatArray> nodes;
    giveFictiousNodeCoordsForExport(nodes, layer);

    FEInterpolation *interpol = static_cast< FEInterpolation * >( &this->interpolationForExport );
    FloatMatrix dNdx;
    interpol->evaldNdx( dNdx, lcoords, FEIVertexListGeometryWrapper( nodes ) ); 
    
    /*    
     * 1 [d/dx  0   d/dy
     * 1   0   d/dy d/dx]
     */
    int ndofs = numNodes*3;
    answer.resize(2, ndofs);
    for ( int i = 1, j = 0; i <= numNodes; i++, j += 3 ) {
        answer.at(1, j + 1) = dNdx.at(i, 1);
        answer.at(1, j + 3) = dNdx.at(i, 2);
        answer.at(2, j + 2) = dNdx.at(i, 2);
        answer.at(2, j + 3) = dNdx.at(i, 1);
    }
    
}





void 
Shell7Base :: giveFictiousNodeCoordsForExport(std::vector<FloatArray> &nodes, int layer)
{
    // compute fictious node coords
    FloatMatrix localNodeCoords;
    this->interpolationForExport.giveLocalNodeCoords(localNodeCoords);
    
    nodes.resize(localNodeCoords.giveNumberOfColumns());
    for ( int i = 1; i <= localNodeCoords.giveNumberOfColumns(); i++ ){
        FloatArray coords, localCoords(3);
        localCoords.beColumnOf(localNodeCoords,i);

        this->vtkEvalInitialGlobalCoordinateAt(localCoords, layer, coords);
        nodes[i-1].resize(3); 
        nodes[i-1] = coords;
    }

}


void 
Shell7Base :: giveFictiousCZNodeCoordsForExport(std::vector<FloatArray> &nodes, int interface)
{
    // compute fictious node coords
    FloatMatrix localNodeCoords;
    this->interpolationForCZExport.giveLocalNodeCoords(localNodeCoords);
    
    nodes.resize(localNodeCoords.giveNumberOfColumns());
    for ( int i = 1; i <= localNodeCoords.giveNumberOfColumns(); i++ ){
        FloatArray coords, localCoords(3);
        localCoords.beColumnOf(localNodeCoords,i);

        localCoords.at(3) = 1.0;
        this->vtkEvalInitialGlobalCoordinateAt(localCoords, interface, coords);
        nodes[i-1].resize(3); 
        nodes[i-1] = coords;
    }

}

void 
Shell7Base :: giveFictiousUpdatedNodeCoordsForExport(std::vector<FloatArray> &nodes, int layer, TimeStep *tStep)
{
    // compute fictious node coords

    FloatMatrix localNodeCoords;
    this->interpolationForExport.giveLocalNodeCoords(localNodeCoords);
    nodes.resize(localNodeCoords.giveNumberOfColumns());
    for ( int i = 1; i <= localNodeCoords.giveNumberOfColumns(); i++ ){
        FloatArray coords, localCoords(3);
        localCoords.beColumnOf(localNodeCoords,i);

        this->vtkEvalUpdatedGlobalCoordinateAt(localCoords, layer, coords, tStep);
        nodes[i-1].resize(3); 
        nodes[i-1] = coords;
    }
}


#endif



// Misc functions

FloatArray 
Shell7Base :: convV6ToV9Stress(const FloatArray &V6)
{
    FloatArray answer(9);
    answer.at(1) = V6.at(1);
    answer.at(5) = V6.at(2);
    answer.at(9) = V6.at(3);
    answer.at(6) = answer.at(8) = V6.at(4);
    answer.at(3) = answer.at(7) = V6.at(5);
    answer.at(2) = answer.at(4) = V6.at(6);
    return answer;
};


void 
Shell7Base :: computeInterLaminarStressesAt(int interfaceNum, TimeStep *tStep, std::vector < FloatArray > &interLamStresses)
{
    
    // Get integration rules on the upper and lower sides of the interface


    LayeredIntegrationRule *irLower = static_cast< LayeredIntegrationRule * >(this->giveIntegrationRule(interfaceNum-1)); // index from 0
    LayeredIntegrationRule *irUpper = static_cast< LayeredIntegrationRule * >(this->giveIntegrationRule(interfaceNum));
    IntegrationPoint *ip;

    // compute stresses
    int numInterfaceIP = irLower->upperInterfacePoints.giveSize();
    interLamStresses.resize(numInterfaceIP);
    FloatArray vSLower, vSUpper, nCov, vSMean(6), stressVec;
    FloatMatrix SMean;
    // compute mean interface stress vector for a gp-pair
    for ( int i = 1; i <= numInterfaceIP; i++ ) {
        ip = irUpper->getIntegrationPoint( irUpper->lowerInterfacePoints.at(i) );
        this->giveIPValue(vSUpper, ip, IST_CauchyStressTensor, tStep);
        ip = irLower->getIntegrationPoint( irLower->upperInterfacePoints.at(i) );
        this->giveIPValue(vSLower, ip, IST_CauchyStressTensor, tStep);

        this->evalInitialCovarNormalAt(nCov, ip->giveNaturalCoordinates());
        vSMean = 0.5 * ( vSUpper + vSLower );
        SMean.beMatrixFormOfStress(vSMean);
        stressVec.beProductOf(SMean,nCov);
        interLamStresses.at(i-1).resize( stressVec.giveSize() );
        interLamStresses.at(i-1) = stressVec;
    }

}


void
Shell7Base :: evaluateFailureCriteriaQuantities(FailureCriteriaStatus *fc, TimeStep *tStep) 
{
    //switch ( fc->giveType() ) {

    //case FC_MaxShearStress:
    //    // Stress ordering (1, 5, 9, 6, 3, 2) = (xx, yy, zz, yz, xz, xy)
    //    int numInterfaces = this->layeredCS->giveNumberOfLayers() - 1;
    //    std::vector < FloatArray > interLamStresses;
    //    fc->quantities.resize(numInterfaces); // will overwrite this every time
    //    for (int i = 1; i <= numInterfaces; i++ ) {    
    //        this->computeInterLaminarStressesAt(i, tStep, interLamStresses); // all 6 components in each evaluation point (ip)
    //        int numEvalPoints = interLamStresses.size();
    //        fc->quantities[i-1].resize(numEvalPoints); // most often = numIP
    //        
    //        for ( int j = 1; j <= numEvalPoints; j++) {
    //            FloatArray &values = fc->quantities[i-1][j-1]; // one resulting shear stress
    //            FloatArray &vS = interLamStresses[j-1];        // Stress in eval point
    //            values.resize(1);                              // scalar measure in this case
    //             
    //            values.at(1) = sqrt( vS.at(2)*vS.at(2) + vS.at(3)*vS.at(3) ); // components can't be right here? shouldn't it be a traction vector?

    //        }
    //    }

    //};
}

int
Shell7Base::giveSymVoigtIndex(int ind1, int ind2)
{
    // Returns the Voigt index corresponding to two given tensor indices for a symmetric tensor.
    // [11 12 13      [1 6 5
    //  21 22 23  = >  6 2 4
    //  31 32 33]      5 4 3]
    //
    std :: vector< std::vector<int> > voigtIndices =
        {
            { 1, 6, 5 },
            { 6, 2, 4 },
            { 5, 4, 3 }
        };
    
    return voigtIndices[ind1][ind2];

};


int
Shell7Base::giveVoigtIndex(int ind1, int ind2)
{
    // Returns the Voigt index corresponding to two given tensor indices for a general tensor.
    // [11 12 13      [1 6 5
    //  21 22 23  = >  9 2 4
    //  31 32 33]      8 7 3]
    //
    //std::vector< std::vector<int> > voigtIndices;
    //voigtIndices.resize(3);
    //voigtIndices[0] = { 1, 6, 5 };
    //voigtIndices[1] = { 9, 2, 4 };
    //voigtIndices[2] = { 8, 7, 3 };

    return this->voigtIndices[ind1-1][ind2-1];

};




} // end namespace oofem
