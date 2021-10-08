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

#include "sm/Elements/Shells/shell7base.h"
#include "sm/Materials/structuralms.h"
#include "sm/Loads/constantpressureload.h"
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
#include "floatarrayf.h"
#include "floatmatrixf.h"
#include "connectivitytable.h"
#include <fstream>

namespace oofem {

FEI3dTrQuad  Shell7Base   :: interpolationForCZExport;
FEI3dWedgeQuad Shell7Base :: interpolationForExport;


Shell7Base :: Shell7Base(int n, Domain *aDomain) : NLStructuralElement(n, aDomain),  LayeredCrossSectionInterface(), 
    VTKXMLExportModuleElementInterface(), ZZNodalRecoveryModelInterface(this), FailureModuleElementInterface(),
    recoverStress(false) 
    {}

void Shell7Base :: initializeFrom(InputRecord &ir)
{
    NLStructuralElement :: initializeFrom(ir);
    //IR_GIVE_OPTIONAL_FIELD(ir, this->recoverStress, _IFT_Shell7base_recoverStress);
    if ( ir.hasField(_IFT_Shell7base_recoverStress) ) {
        this->recoverStress = true;
    }

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

void 
Shell7Base::printOutputAt(FILE* file, TimeStep* tStep)
{
    if ( recoverStress ) {
        // Recover shear stresses
        this->recoverShearStress(tStep);
    }
    oofem::Element::printOutputAt(file, tStep);
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
    //int layer = 1;
    //vtkEvalInitialGlobalCoordinateAt(lcoords, layer, answer);
    
    
    double zeta = giveGlobalZcoord(lcoords);
    FloatArray N;
    this->fei->evalN( N, lcoords, FEIElementGeometryWrapper(this) );

    answer.clear();
    for ( int i = 1; i <= this->giveNumberOfDofManagers(); i++ ) {
        const auto &xbar = this->giveNode(i)->giveCoordinates();
        const auto &M = this->giveInitialNodeDirector(i);
        answer.add(N.at(i), ( xbar + zeta * M ));
    }

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
Shell7Base::giveGlobalZcoord( const FloatArrayF<3> &lCoords )
{
    return lCoords.at(3) * this->layeredCS->give( CS_Thickness, lCoords, this, false ) * 0.5;
}


double /// @todo move to layered crosssection
Shell7Base :: giveGlobalZcoordInLayer(double xi, int layer)
{
    // Mid position + part of the layer thickness
    return this->layeredCS->giveLayerMidZ(layer) + xi*this->layeredCS->giveLayerThickness(layer)*0.5 ;
}


// Base vectors

#if 1

FloatMatrixF<3,3>
Shell7Base :: evalInitialCovarBaseVectorsAt(const FloatArrayF<3> &lcoords)
{
    double zeta = giveGlobalZcoord( lcoords );
    FloatMatrix dNdxi;

    // In plane base vectors
    this->fei->evaldNdxi( dNdxi, lcoords, FEIElementGeometryWrapper(this) );

    FloatArrayF<3> G1, G2;
    for ( int i = 1; i <= this->giveNumberOfDofManagers(); i++ ) {
        const auto &xbar = FloatArrayF<3>(this->giveNode(i)->giveCoordinates());
        const auto &M = this->giveInitialNodeDirector(i);
        const auto &nodeCoords = (xbar + zeta*M);
        G1 += dNdxi.at(i, 1) * nodeCoords;
        G2 += dNdxi.at(i, 2) * nodeCoords;
    }

    // Out of plane base vector = director
    auto G3 = this->evalInitialDirectorAt(lcoords);     // G3=M

    FloatMatrixF<3,3> gcov;
    gcov.setColumn(G1, 0);
    gcov.setColumn(G2, 1);
    gcov.setColumn(G3, 2);
    return gcov;
}


std::pair<FloatArrayF<3>, FloatArrayF<3>>
Shell7Base :: edgeEvalInitialCovarBaseVectorsAt(const FloatArrayF<1> &lcoords, const int iedge)
{
    double zeta = 0.0;     // no variation i z (yet)

    const auto &edgeNodes = this->fei->computeLocalEdgeMapping(iedge);
    FloatArray dNdxi;
    this->fei->edgeEvaldNdxi( dNdxi, iedge, lcoords, FEIElementGeometryWrapper(this) );

    // Base vector along edge
    FloatArrayF<3> G1;
    for ( int i = 1; i <= edgeNodes.giveSize(); i++ ) {
        const auto &xbar = FloatArrayF<3>(this->giveNode(edgeNodes.at(i))->giveCoordinates());
        auto M = this->giveInitialNodeDirector(edgeNodes.at(i));
        auto nodeCoords = (xbar + zeta*M);
        G1 += dNdxi.at(i) * nodeCoords;
    }

    // Director will be the second base vector
    auto G3 = this->edgeEvalInitialDirectorAt(lcoords, iedge);
    
    return {G1, G3};
}


FloatMatrixF<3,3>
Shell7Base :: evalInitialContravarBaseVectorsAt(const FloatArrayF<3> &lCoords)
{
    auto Gcov = this->evalInitialCovarBaseVectorsAt(lCoords);
    return this->giveDualBase(Gcov);
}


FloatMatrixF<3,3>
Shell7Base :: giveDualBase( FloatMatrixF<3,3> &base1)
{
    // Computes the dual base through thte inversion of the metric tensor
    auto gMetric = Tdot(base1,base1); // Metric tensor
    auto ginv = inv(gMetric);
    return dotT(base1, ginv);
}


FloatArrayF<3>
Shell7Base :: evalInitialDirectorAt(const FloatArrayF<3> &lcoords)
{
    // Interpolates between the node directors
    FloatArray N;
    this->fei->evalN( N, lcoords, FEIElementGeometryWrapper(this) );
    FloatArrayF<3> g3;
    for ( int i = 1; i <= this->giveNumberOfDofManagers(); i++ ) {
        g3 += N.at(i) * this->giveInitialNodeDirector(i);
    }
    return g3;
}


FloatArrayF<3>
Shell7Base :: edgeEvalInitialDirectorAt(const FloatArrayF<1> &lcoords, const int iEdge)
{
    // Interpolates between the node directors along an edge

    FloatArray N;
    const auto &edgeNodes = this->fei->computeLocalEdgeMapping(iEdge);
    this->fei->edgeEvalN( N, iEdge, lcoords, FEIElementGeometryWrapper(this) );

    FloatArrayF<3> answer;
    for ( int i = 1; i <= edgeNodes.giveSize(); i++ ) {
        answer += N.at(i) * this->giveInitialNodeDirector( edgeNodes.at(i) );
    }
    return answer;
}


void
Shell7Base :: setupInitialNodeDirectors()
{
    // Compute directors as normals to the surface
    FloatArray lcoords;
    FloatMatrix localNodeCoords;
    this->giveInterpolation()->giveLocalNodeCoords(localNodeCoords);
    
    int nDofMan = this->giveNumberOfDofManagers();
    this->initialNodeDirectors.resize(nDofMan);
    FloatMatrix dNdxi;
    for ( int node = 1; node <= nDofMan; node++ ) {
        lcoords.beColumnOf(localNodeCoords,node);      
        this->fei->evaldNdxi( dNdxi, lcoords, FEIElementGeometryWrapper(this) );

        FloatArrayF<3> G1, G2;
        // base vectors of the initial surface
        for ( int i = 1; i <= nDofMan; i++ ) {        
            const auto &nodeI = FloatArrayF<3>(this->giveNode(i)->giveCoordinates());
            G1 += dNdxi.at(i, 1) * nodeI;
            G2 += dNdxi.at(i, 2) * nodeI;
        }
        this->initialNodeDirectors [ node - 1 ] = normalize(cross(G1, G2));
    }
}


FloatMatrixF<3,3>
Shell7Base :: evalCovarBaseVectorsAt(const FloatArrayF<3> &lcoords, FloatArray &genEps, TimeStep *tStep)
{
    // Evaluates the covariant base vectors in the current configuration
    double zeta = giveGlobalZcoord( lcoords );

    FloatArrayF<3> dxdxi1, dxdxi2, m, dmdxi1, dmdxi2;
    double dgamdxi1, dgamdxi2, gam;
    this->giveGeneralizedStrainComponents(genEps, dxdxi1, dxdxi2, dmdxi1, dmdxi2, m, dgamdxi1, dgamdxi2, gam);

    double fac1 = ( zeta + 0.5 * gam * zeta * zeta );
    double fac2 = ( 0.5 * zeta * zeta );
    double fac3 = ( 1.0 + zeta * gam );

    auto g1 = dxdxi1 + fac1*dmdxi1 + fac2*dgamdxi1*m;
    auto g2 = dxdxi2 + fac1*dmdxi2 + fac2*dgamdxi2*m;
    auto g3 = fac3*m;

    FloatMatrixF<3,3> gcov;
    gcov.setColumn(g1,0); gcov.setColumn(g2,1); gcov.setColumn(g3,2);
    return gcov;

}


FloatMatrixF<3,3>
Shell7Base :: edgeEvalCovarBaseVectorsAt(const FloatArrayF<3> &lcoords, const int iedge, TimeStep *tStep)
{
    // Evaluates the covariant base vectors in the current configuration for an edge
    double zeta = lcoords.at(3);

    FloatArray solVecEdge;
    FloatMatrix B;
    //const auto &edgeNodes = this->fei->computeLocalEdgeMapping(iedge);
    this->edgeComputeBmatrixAt(lcoords, B, 1, ALL_STRAINS);
    this->edgeGiveUpdatedSolutionVector(solVecEdge, iedge, tStep);

    FloatArray genEpsEdge;                 // generalized strain
    genEpsEdge.beProductOf(B, solVecEdge); // [dxdxi, dmdxi, m, dgamdxi, gam]^T

    FloatArrayF<3> dxdxi = { genEpsEdge.at(1), genEpsEdge.at(2), genEpsEdge.at(3) };
    FloatArrayF<3> dmdxi = { genEpsEdge.at(4), genEpsEdge.at(5), genEpsEdge.at(6) };
    FloatArrayF<3> m = { genEpsEdge.at(7), genEpsEdge.at(8), genEpsEdge.at(9) };
    double dgamdxi = genEpsEdge.at(10);
    double gam     = genEpsEdge.at(11);

    double fac1 = ( zeta + 0.5 * gam * zeta * zeta );
    double fac2 = ( 0.5 * zeta * zeta );
    double fac3 = ( 1.0 + zeta * gam );
    
    auto g2 = normalize(dxdxi + fac1*dmdxi + fac2*dgamdxi*m); // base vector along the edge
    auto g3 = normalize(fac3*m);                              // director field
    auto g1 = normalize(cross(g2, g3));

    FloatMatrixF<3,3> gcov;
    gcov.setColumn(g1, 0);
    gcov.setColumn(g2, 1);
    gcov.setColumn(g3, 2);
    return gcov;
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

std::array<FloatMatrixF<3,18>, 3>
Shell7Base :: computeLambdaGMatrices(FloatArray &genEps, double zeta)
{
    // computes the lambda^g matrices associated with the variation and linearization of the base vectors g_i.
    // \delta g_i = lambda_i * \delta n and \Delta g_i = lambda_i * \Delta n 
    // \delta n = B * \delta a and \Delta n = B * \Delta a
    // @todo optimize method
    FloatArrayF<3> m, dm1, dm2, temp1;
    double dgam1, dgam2, gam;
    this->giveGeneralizedStrainComponents(genEps, temp1, temp1, dm1, dm2, m, dgam1, dgam2, gam);

    // thickness coefficients
    double a = zeta + 0.5 * gam * zeta * zeta;
    double b = 0.5 * zeta * zeta;
    double c = 1.0 + gam * zeta;

    auto bm = b * m;
    auto bdm1 = b * dm1;
    auto bdm2 = b * dm2;

    std::array<FloatMatrixF<3,18>, 3> lambda;

    // lambda1 =  ( I,   0,  a*I,   0 ,  b*dgam1*I,  b*m,   0 ,  b*dm1 )
    lambda[ 0 ](0,0) = lambda[ 0 ](1,1) = lambda[ 0 ](2,2) = 1.;
    lambda[ 0 ](0,6) = lambda[ 0 ](1,7) = lambda[ 0 ](2,8) = a;
    lambda[ 0 ](0,12) = lambda[ 0 ](1,13) = lambda[ 0 ](2,14) = b * dgam1;
    lambda[ 0 ].setColumn(bm, 15);
    lambda[ 0 ].setColumn(bdm1, 17);
    
    // lambda2 =  ( 0,   I,   0 ,  a*I,  b*dgam2*I,   0 ,  b*m,  b*dm2 )
    lambda[ 1 ](0,3) = lambda[ 1 ](1,4) = lambda[ 1 ](2,5) = 1.;
    lambda[ 1 ](0,9) = lambda[ 1 ](1,10) = lambda[ 1 ](2,11) = a;
    lambda[ 1 ](0,12) = lambda[ 1 ](1,13) = lambda[ 1 ](2,14) = b * dgam2;
    lambda[ 1 ].setColumn(bm, 16);
    lambda[ 1 ].setColumn(bdm2, 17);

    // lambda3 =  ( 0,   0,   0 ,   0 ,     c*I   ,   0 ,   0 ,   xi*m )
    lambda[ 2 ](0,12) = lambda[ 2 ](1,13) = lambda[ 2 ](2,14) = c;
    auto zm = zeta * m;
    lambda[ 2 ].setColumn(zm, 17);

    return lambda;
}


FloatMatrixF<3,7>
Shell7Base :: computeLambdaNMatrix(FloatArray &genEps, double zeta)
{
    // computes the lambda^n matrix associated with the variation and linearization of the position vector x.
    // \delta x = lambda * \delta \hat{x} with \hat{x} = [\bar{x}, m, \gamma]

    FloatArrayF<3> m = { genEps.at(13), genEps.at(14), genEps.at(15) };
    double gam = genEps.at(18);

    // thickness coefficients
    double a = zeta + 0.5 * gam * zeta * zeta;
    double b = 0.5 * zeta * zeta;
    
    // lambda =  ( I, a*I, b*m )
    FloatMatrixF<3,7> lambda;
    lambda.at(1,1) = lambda.at(2,2) = lambda.at(3,3) = 1.0;
    lambda.at(1,4) = lambda.at(2,5) = lambda.at(3,6) = a;
    lambda.setColumn(b*m, 6); 
    return lambda;
}


void
Shell7Base :: computeBulkTangentMatrix(FloatMatrix &answer, FloatArray &solVec, TimeStep *tStep)
{
    FloatMatrix A [ 3 ] [ 3 ], A_lambda(3,18), LB;
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
            const auto &lCoords = gp->giveNaturalCoordinates();

            this->computeBmatrixAt(lCoords, B);
            genEps.beProductOf(B, solVec);
            // Material stiffness
            Shell7Base :: computeLinearizedStiffness(gp, mat, tStep, A);

            double zeta = giveGlobalZcoord(lCoords);
            auto lambda = this->computeLambdaGMatrices(genEps, zeta);

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
    
    const auto &ordering = this->giveOrderingDofTypes();
    answer.assemble(tempAnswer, ordering, ordering);

}


void
Shell7Base :: computeLinearizedStiffness(GaussPoint *gp, StructuralMaterial *mat, TimeStep *tStep, FloatMatrix A [ 3 ] [ 3 ]) 
{
    const auto &lcoords = gp->giveNaturalCoordinates();

    // Material stiffness when internal work is formulated in terms of P and F:
    // \Delta(P*G^I) = L^IJ * \Delta g_J
    // A[I][J] = L^IJ = L_klmn * [G^I]_l * [G^J]_n
    auto D = mat->give3dMaterialStiffnessMatrix_dPdF(TangentStiffness, gp, tStep);    // D_ijkl - cartesian system (Voigt)
    auto G = this->evalInitialContravarBaseVectorsAt(lcoords);
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
    
    FloatMatrix N, B, LB, NLB, L(7, 18);
    FloatArray lcoords(3), solVec, pressure;
    FloatArray genEps;

    double xi = pLoad->giveLoadOffset();
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
        auto gcov = this->evalCovarBaseVectorsAt(lcoords, genEps, tStep);
        auto g1 = gcov.column(0);
        auto g2 = gcov.column(1);
        auto W1 = this->giveAxialMatrix(g1);
        auto W2 = this->giveAxialMatrix(g2);
        
        auto lambdaG = this->computeLambdaGMatrices(genEps, zeta);
        auto lambdaN = this->computeLambdaNMatrix(genEps, zeta);
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


FloatMatrixF<3,3>
Shell7Base :: giveAxialMatrix(const FloatArrayF<3> &v)
{
    // creates the skew-symmetric matrix W defined such that 
    // crossProduct(u,v) = W(v)*u
    // W = [   0   v(3)  -v(2)
    //      -v(3)   0     v(1)
    //       v(2) -v(1)    0  ]
    //
    FloatMatrixF<3,3> answer;
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

FloatMatrixF<3,3>
Shell7Base :: computeFAt(const FloatArrayF<3> &lCoords, FloatArray &genEps, TimeStep *tStep)
{
    // Computes the deformation gradient in matrix form as open product(g_i, G^i) = gcov*Gcon^T
    auto gcov = this->evalCovarBaseVectorsAt(lCoords, genEps, tStep);
    auto Gcon = this->evalInitialContravarBaseVectorsAt(lCoords);
    return dotT(gcov, Gcon);
}


FloatMatrixF<3,3>
Shell7Base :: computeStressMatrix(FloatArray &genEps, GaussPoint *gp, Material *mat, TimeStep *tStep)
{
    auto F = computeFAt(gp->giveNaturalCoordinates(), genEps, tStep);
    auto vF = to_voigt_form(F);
    auto vP = static_cast< StructuralMaterial * >( mat )->giveFirstPKStressVector_3d(vF, gp, tStep);
    return from_voigt_form(vP);
}


void
Shell7Base :: computeCauchyStressVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep)
{
    // Compute Cauchy stress from 2nd Piola stress
    FloatArray solVec;
    this->giveUpdatedSolutionVector(solVec, tStep); 
    const auto &lCoords = gp->giveNaturalCoordinates();
    FloatMatrix B;
    this->computeBmatrixAt(lCoords, B);
    FloatArray genEps;
    genEps.beProductOf(B, solVec);

    auto F = this->computeFAt(lCoords, genEps, tStep);   

    FloatArray vS;
    giveIPValue(vS, gp, IST_StressTensor, tStep); // expects Second PK stress

    FloatMatrix S, temp, sigma;
    S.beMatrixFormOfStress(vS);
    temp.beProductTOf(S,F); 
    sigma.beProductOf(F,temp);
    sigma.times( 1.0/ det(F) );
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
    auto P = this->computeStressMatrix(genEps, ip, mat, tStep);
    auto Gcon = this->evalInitialContravarBaseVectorsAt(ip->giveNaturalCoordinates());
    auto PG = dot(P, Gcon);
    auto PG1 = PG.column(0);
    auto PG2 = PG.column(1);
    auto PG3 = PG.column(2);
    auto lambda = this->computeLambdaGMatrices(genEps, zeta); // associated with the variation of the test functions   

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
    const auto &lcoords = gp->giveNaturalCoordinates();

    FloatMatrix dNdxi;
    this->fei->evaldNdxi( dNdxi, lcoords, FEIElementGeometryWrapper(this) );

    FloatArrayF<3> M, dM1, dM2, dX1, dX2;
    double gam, dg1, dg2;
    FloatMatrix B;
    this->computeBmatrixAt(lcoords, B);

    FloatArray initSolVec, genEps;
    initSolVec = this->giveInitialSolutionVector();
    genEps.beProductOf(B, initSolVec);
    this->giveGeneralizedStrainComponents(genEps, dX1, dX2, dM1, dM2, M, dg1, dg2, gam);

    auto temp = cross(dX1, dX2);
    double sc = norm(temp);
    answer.resize(3);
    answer.at(3) = dot(M, temp) / sc;

    temp = cross(dX1, dM2) + cross(dX2, dM1);
    answer.at(2) = dot(M, temp) / sc;

    temp = cross(dM1, dM2);
    //answer.at(1) = dot(M, temp)/sc;
    answer.at(1) = norm(temp) / sc;
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
            FloatMatrix B, N, temp;
            FloatArray genEps;
            this->computeBmatrixAt(lCoords, B);
            genEps.beProductOf(B, solVec);    
            double zeta = giveGlobalZcoord(gp->giveNaturalCoordinates());
            auto lambda = this->computeLambdaNMatrix(genEps, zeta);

            // could also create lambda*N and then plusProdSymm - probably faster
            mass = Tdot(lambda, lambda);
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
    iRule.SetUpPointsOnWedge(nPointsTri, 1, _3dMat); ///@todo replace with triangle which has a xi3-coord

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
void Shell7Base :: computeBoundaryEdgeLoadVector(FloatArray &answer, BoundaryLoad *load, int boundary, CharType type, ValueModeType mode, TimeStep *tStep, bool global)
{
    answer.clear();
    if ( type != ExternalForcesVector ) {
        return;
    }
    this->computeTractionForce(answer, boundary, load, tStep, mode);
}


/*
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
    } else {
        OOFEM_ERROR("Load type not supported");
    }
}
*/
  
void
Shell7Base :: computePressureForce(FloatArray &answer, FloatArray solVec, const int iSurf, BoundaryLoad *surfLoad, TimeStep *tStep, ValueModeType mode)
{
    // Computes pressure loading. Acts normal to the current (deformed) surface.
//     ConstantPressureLoad* pLoad = dynamic_cast< ConstantPressureLoad * >( surfLoad );
//     ConstantSurfaceLoad*  sLoad = dynamic_cast< ConstantSurfaceLoad *  >( surfLoad );
    double xi = 0.0;
    if ( ConstantPressureLoad* pLoad = dynamic_cast< ConstantPressureLoad * >( surfLoad ) ) {
        xi = pLoad->giveLoadOffset();
        //printf("pLoad xi = %e \n",xi);
    } 
    else if ( ConstantSurfaceLoad* sLoad = dynamic_cast< ConstantSurfaceLoad * >( surfLoad ) ) {
        xi = sLoad->giveLoadOffset();
        //printf("sLoad xi = %e \n",xi);
    } 
    
    int nPointsTri = 6; //todo generalize
    auto iRule = std::make_unique<GaussIntegrationRule>(1, this);
//     iRule->SetUpPointsOnWedge(nPointsTri, 1, _3dMat); //@todo replce with triangle which has a xi3-coord
    iRule->SetUpPointsOnTriangle(nPointsTri, _3dMat);
    
    
    FloatMatrix N, B, lambda;
    FloatArray Fp, fp, genEps, genEpsC, traction, solVecC;

    for ( auto *ip: *iRule ) { // rule #2 for surface integration
        FloatArrayF<3> lCoords = {
            ip->giveNaturalCoordinate(1),
            ip->giveNaturalCoordinate(2),
            xi,
        };
        double zeta = giveGlobalZcoord( lCoords );
        this->computeBmatrixAt(lCoords, B);
        this->computeNmatrixAt(lCoords, N);
        this->giveUpdatedSolutionVector(solVecC, tStep);
        genEpsC.beProductOf(B, solVecC);
        this->computePressureForceAt(ip, traction, iSurf, genEpsC, surfLoad, tStep, mode);

        genEps.beProductOf(B, solVec);
        lambda = this->computeLambdaNMatrix(genEps, zeta);
        fp.beTProductOf(lambda,traction);

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

        gcov = this->evalCovarBaseVectorsAt(lcoords, genEps, tStep); 
        g1.beColumnOf(gcov,1);
        g2.beColumnOf(gcov,2);
        surfLoad->computeValueAt(load, tStep, lcoords, mode);        // pressure component
        traction.beVectorProductOf(g1, g2);        // normal vector (should not be normalized due to integraton in reference config.)
        traction.times( -load.at(1) );
    } else if ( ConstantSurfaceLoad* sLoad = dynamic_cast< ConstantSurfaceLoad * >( surfLoad ) ) {
        FloatArray lcoords(3), tempTraction;
        lcoords.at(1) = gp->giveNaturalCoordinate(1);
        lcoords.at(2) = gp->giveNaturalCoordinate(2);
        lcoords.at(3) = sLoad->giveLoadOffset();
        surfLoad->computeValueAt(tempTraction, tStep, lcoords, mode);        // traction vector
        traction.resize(3);
        traction.zero();
        //IntArray sLoadDofIDs = sLoad->giveDofIDs();
        //sLoadDofIDs.printYourself("sLoadDofIDs");
        traction.assemble(tempTraction,sLoad->giveDofIDs());
    } else {
        OOFEM_ERROR("incompatible load type");
    }

}

FloatArrayF<3>
Shell7Base :: evalCovarNormalAt(const FloatArrayF<3> &lCoords, FloatArray &genEpsC, TimeStep *tStep)
{
    auto gcov = this->evalCovarBaseVectorsAt(lCoords, genEpsC, tStep);
    auto g1 = gcov.column(0);
    auto g2 = gcov.column(1);
    return normalize(cross(g1, g2));
}

FloatArrayF<3>
Shell7Base :: evalInitialCovarNormalAt(const FloatArrayF<3> &lCoords)
{
    auto Gcov = this->evalInitialCovarBaseVectorsAt(lCoords);
    auto G1 = Gcov.column(0);
    auto G2 = Gcov.column(1);
    return normalize(cross(G1, G2));
}

void
Shell7Base :: computeTractionForce(FloatArray &answer, const int iEdge, BoundaryLoad *edgeLoad, TimeStep *tStep, ValueModeType mode, bool map2elementDOFs)
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
        fT.zero();
        edgeLoad->computeValueAt(components, tStep, lCoords, mode);
        this->edgeComputeNmatrixAt(lCoords, N);

        //if ( coordSystType ==  BoundaryLoad :: BL_UpdatedGlobalMode ) {
        if ( coordSystType ==  Load :: CST_UpdatedGlobal ) {
            
            // Updated global coord system
            auto gcov = this->edgeEvalCovarBaseVectorsAt(lCoords, iEdge, tStep); 
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
    /* Note: now used from computeBoundaryEdgeLoadVector, so result should be in global cs for edge dofs */
    if (map2elementDOFs) {
      IntArray mask;
      this->giveEdgeDofMapping(mask, iEdge);
      answer.resize( Shell7Base :: giveNumberOfDofs()  );
      answer.zero();
      answer.assemble(Nf, mask);
    } else {
      answer = Nf;
    }

}


void
Shell7Base :: computeBodyLoadVectorAt(FloatArray &answer, Load *forLoad, TimeStep *tStep, ValueModeType mode)
{
    OOFEM_ERROR("Shell7Base :: computeBodyLoadVectorAt - currently not implemented");
}

#endif



// Integration
#if 1

///@todo should be moved to tr-elements
double
Shell7Base :: edgeComputeLengthAround(GaussPoint *gp, const int iedge)
{
    const auto &lcoords = gp->giveNaturalCoordinates();
    //auto [G1, G3] = this->edgeEvalInitialCovarBaseVectorsAt(lcoords, iedge);
    auto tmp = this->edgeEvalInitialCovarBaseVectorsAt(lcoords, iedge);
    auto G1 = tmp.first;
    double detJ = norm(G1);
    return detJ * gp->giveWeight();
}
#endif



// Recovery of nodal values

#if 1

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
{
    // evaluates N^T sigma over element volume
    // N(nsigma, nsigma*nnodes)
    // Definition : sigmaVector = N * nodalSigmaVector
    FloatArray stressVector, n;

    int size = 6; ///@todo this is hard coded for stress recovery

    answer.zero();
    for ( auto *gp: *integrationRulesArray [ layer - 1 ] ) {
        double dV = this->computeVolumeAroundLayer(gp, layer);

        if ( !this->giveIPValue(stressVector, gp, type, tStep) ) {
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
            temp.at(j) = nValProd.at(i,j) / nnMatrix.at(i);
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

    const auto &bNodes = this->fei->computeLocalEdgeMapping(boundary);
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
        const auto &Xi = this->giveNode(i)->giveCoordinates();
        const auto &Mi = this->giveInitialNodeDirector(i);
        this->initialSolutionVector.at(1 + j) = Xi.at(1);
        this->initialSolutionVector.at(2 + j) = Xi.at(2);
        this->initialSolutionVector.at(3 + j) = Xi.at(3);
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
        const auto &edgeNodes = this->fei->computeLocalEdgeMapping(iEdge);
        int ndofs_x = 3 * edgeNodes.giveSize();
        for ( int i = 1, j = 0; i <= edgeNodes.giveSize(); i++, j += 3 ) {
            const auto &Xi = this->giveNode( edgeNodes.at(i) )->giveCoordinates();
            const auto &Mi = this->giveInitialNodeDirector( edgeNodes.at(i) );
            solVec.at(1 + j) = Xi.at(1);
            solVec.at(2 + j) = Xi.at(2);
            solVec.at(3 + j) = Xi.at(3);
            solVec.at(ndofs_x + 1 + j) = Mi.at(1);
            solVec.at(ndofs_x + 2 + j) = Mi.at(2);
            solVec.at(ndofs_x + 3 + j) = Mi.at(3);
            // gam(t=0)=0 is assumed
        }
    }
}

void
Shell7Base :: giveGeneralizedStrainComponents(FloatArray genEps, FloatArrayF<3> &dphidxi1, FloatArrayF<3> &dphidxi2, FloatArrayF<3> &dmdxi1,
                                              FloatArrayF<3> &dmdxi2, FloatArrayF<3> &m, double &dgamdxi1, double &dgamdxi2, double &gam)
{
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
Shell7Base :: giveUnknownsAt(const FloatArrayF<3> &lCoords, const FloatArray &solVec, FloatArrayF<3> &x, FloatArrayF<3> &m, double &gam, TimeStep *tStep)
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
    // Returns the  matrix {B} of the receiver, evaluated at gp. Such that
    // B*a = [dxbar_dxi, dwdxi, w, dgamdxi, gam]^T, where a is the vector of unknowns

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
FloatArrayF<3>
Shell7Base :: vtkEvalInitialGlobalCoordinateAt(const FloatArrayF<3> &localCoords, int layer)
{
    double zeta = giveGlobalZcoordInLayer(localCoords.at(3), layer);
    FloatArray N;
    this->fei->evalN( N, localCoords, FEIElementGeometryWrapper(this) );

    FloatArrayF<3> globalCoords;
    for ( int i = 1; i <= this->giveNumberOfDofManagers(); i++ ) {
        const auto &xbar = FloatArrayF<3>(this->giveNode(i)->giveCoordinates());
        const auto &M = this->giveInitialNodeDirector(i);
        globalCoords += N.at(i) * ( xbar + zeta * M );
    }
    return globalCoords;
}

FloatArrayF<3>
Shell7Base :: vtkEvalInitialGlobalCZCoordinateAt(const FloatArrayF<3> &localCoords, int interface)
{
    double zeta = giveGlobalZcoordInLayer(1.0, interface);
    FloatArray N;
    this->fei->evalN( N, localCoords, FEIElementGeometryWrapper(this) );

    FloatArrayF<3> globalCoords;
    for ( int i = 1; i <= this->giveNumberOfDofManagers(); i++ ) {
        const auto &xbar = FloatArrayF<3>(this->giveNode(i)->giveCoordinates());
        const auto &M = this->giveInitialNodeDirector(i);
        globalCoords += N.at(i) * ( xbar + zeta * M );
    }
    return globalCoords;
}

FloatArrayF<3>
Shell7Base :: vtkEvalUpdatedGlobalCoordinateAt(const FloatArrayF<3> &localCoords, int layer, TimeStep *tStep)
{
    FloatArray solVec;
    this->giveUpdatedSolutionVector(solVec, tStep);
    FloatArrayF<3> x, m; double gam=0;
    this->giveUnknownsAt(localCoords, solVec, x, m, gam, tStep); 
    double zeta = giveGlobalZcoordInLayer(localCoords.at(3), layer);
    double fac = ( zeta + 0.5 * gam * zeta * zeta );
    return x + fac*m;
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

    int val    = 1;
    int offset = 0;
    IntArray nodes(numCellNodes);

    // Compute fictious node coords
    int nodeNum = 1;
    for ( int layer = 1; layer <= numCells; layer++ ) {

        // Node coordinates
        auto nodeCoords = this->giveFictiousNodeCoordsForExport(layer);
        
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
    vtkPiece.setNumberOfPrimaryVarsToExport(primaryVarsToExport, numTotalNodes);

    FloatArray u(3);
    std::vector<FloatArray> values;
    for ( int fieldNum = 1; fieldNum <= primaryVarsToExport.giveSize(); fieldNum++ ) {
        UnknownType type = ( UnknownType ) primaryVarsToExport.at(fieldNum);
        nodeNum = 1;
        for ( int layer = 1; layer <= numCells; layer++ ) {
            
            if ( type == DisplacementVector ) { // compute displacement as u = x - X
                auto nodeCoords = this->giveFictiousNodeCoordsForExport(layer);
                auto updatedNodeCoords = this->giveFictiousUpdatedNodeCoordsForExport(layer, tStep);
                for ( int j = 1; j <= numCellNodes; j++ ) {
                    u = updatedNodeCoords[j-1];
                    u.subtract(nodeCoords[j-1]);
                    vtkPiece.setPrimaryVarInNode(type, nodeNum, u);
                    nodeNum += 1;        
                }

            } else {
                NodalRecoveryMI_recoverValues(values, layer, ( InternalStateType ) 1, tStep); // does not work well - fix
                for ( int j = 1; j <= numCellNodes; j++ ) {
                    vtkPiece.setPrimaryVarInNode(type, nodeNum, values[j-1]);
                    nodeNum += 1;
                }
            }
        }
    }

    // Export nodal variables from internal fields
    
    vtkPiece.setNumberOfInternalVarsToExport( internalVarsToExport, numTotalNodes );
    for ( int fieldNum = 1; fieldNum <= internalVarsToExport.giveSize(); fieldNum++ ) {
        InternalStateType type = ( InternalStateType ) internalVarsToExport.at(fieldNum);
        nodeNum = 1;
        
//         if ( recoverStress ) {
//         // Recover shear stresses
//         this->recoverShearStress(tStep);
//         }
        
        for ( int layer = 1; layer <= numCells; layer++ ) {            
            recoverValuesFromIP(values, layer, type, tStep);        
            for ( int j = 1; j <= numCellNodes; j++ ) {
                vtkPiece.setInternalVarInNode( type, nodeNum, values[j-1] );
                //ZZNodalRecoveryMI_recoverValues(el.nodeVars[fieldNum], layer, type, tStep);
                nodeNum += 1;
            }
        }
    }


    // Export cell variables
    FloatArray average;
    vtkPiece.setNumberOfCellVarsToExport(cellVarsToExport, numCells);
    for ( int i = 1; i <= cellVarsToExport.giveSize(); i++ ) {
        InternalStateType type = ( InternalStateType ) cellVarsToExport.at(i);

        for ( int layer = 1; layer <= numCells; layer++ ) {
            std :: unique_ptr< IntegrationRule > &iRuleL = integrationRulesArray [ layer - 1 ];
            VTKBaseExportModule::computeIPAverage(average, iRuleL.get(), this, type, tStep);
            
            if ( average.giveSize() == 6 ) {
                vtkPiece.setCellVar(type, layer, convV6ToV9Stress(average) );
            } else {
                vtkPiece.setCellVar(type, layer, average );
            }
        }
    }
}


void 
Shell7Base :: recoverValuesFromIP(std::vector<FloatArray> &recoveredValues, int layer, InternalStateType type, TimeStep *tStep, stressRecoveryType SRtype)
{
    // recover nodal values from IP:s
    
    switch (SRtype) {
        case copyIPvalue:
            // Find closest ip to the nodes
            CopyIPvaluesToNodes(recoveredValues, layer, type, tStep);
            break;
        case LSfit:
            // Do least square fit
            nodalLeastSquareFitFromIP(recoveredValues, layer, type, tStep);
            break;
        case L2fit:
            // Do L2 fit
            OOFEM_ERROR("L2-fit not implemented.");
            break;
        default:
            OOFEM_ERROR("Incorrect stress recovery type.");
    }
#if 0
    if (this->giveGlobalNumber() == 20) {
        for (FloatArray inod: recoveredValues) {
            inod.printYourself("recovered-värden");
        }
    }
#endif
    
}


void 
Shell7Base :: CopyIPvaluesToNodes(std::vector<FloatArray> &recoveredValues, int layer, InternalStateType type, TimeStep *tStep)
{
    // Method for copying value in closest IP to nodes
    
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
        for ( int j = 0; j < iRule->giveNumberOfIntegrationPoints(); ++j ) {
            IntegrationPoint *ip = iRule->getIntegrationPoint(j);
            const FloatArray &ipCoords = ip->giveNaturalCoordinates();
            ///TODO z-coord in parent shell. 
            double dist = distance(nodeCoords, ipCoords);
            if ( dist < distOld ) {
                closestIPArray.at(i) = j;
                distOld = dist;
            }
        }
    }

    // recover ip values
    for ( int i = 1; i <= numNodes; i++ ) {
        IntegrationPoint *ip = integrationRulesArray [ layer - 1 ]->getIntegrationPoint( closestIPArray.at(i) );
        this->giveIPValue(ipValues, ip, type, tStep);
        if ( giveInternalStateValueType(type) == ISVT_TENSOR_S3 ) {
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
Shell7Base :: nodalLeastSquareFitFromIP(std::vector<FloatArray> &recoveredValues, int layer, InternalStateType type, TimeStep *tStep)
{   
    // Method for computing nodal values from IP:s by least square fit
    
    // composite element interpolator
    FloatMatrix localNodeCoords;
    this->interpolationForExport.giveLocalNodeCoords(localNodeCoords);

    int numNodes = localNodeCoords.giveNumberOfColumns();
    recoveredValues.resize(numNodes);
    
    std :: unique_ptr< IntegrationRule > &iRule = integrationRulesArray [ layer - 1 ];
    IntegrationPoint *ip;
    int numIP = iRule->giveNumberOfIntegrationPoints();
    FloatArray nodeCoords, ipCoords;
    InternalStateValueType valueType =  giveInternalStateValueType(type);
    
    if (numNodes > numIP) {
        OOFEM_ERROR("Least square fit not possible for more nodes than IP:s per layer.");
    }
    
    // Find IP values and set up matrix of base functions
    FloatMatrix Nbar, NbarTNbar, NbarTNbarInv, Nhat, ipValues, temprecovedValues, temprecovedValuesT; 
    Nbar.resize(numIP,numNodes);
    //int numSC; 
    //if ( valueType == ISVT_TENSOR_S3 ) { numSC = 6; } else { numSC = 9; };
    //ipValues.resize(numIP,numSC);
    
    for ( int i = 0; i < numIP; i++ ) {
        ip = iRule->getIntegrationPoint(i);
        FloatArray tempIPvalues;
        this->giveIPValue(tempIPvalues, ip, type, tStep);
#if 0
        // Test of analytical dummy values in IPs
        FloatArray ipGlobalCoords;
        ipCoords = *ip->giveNaturalCoordinates();
        this->vtkEvalInitialGlobalCoordinateAt(ipCoords, layer, ipGlobalCoords);
//         this->computeGlobalCoordinates( ipGlobalCoords, ipCoords );
        tempIPvalues.resize(9);
        tempIPvalues.zero();
        tempIPvalues.at(1) = 7.8608e+09 * ( 0.2 - ipGlobalCoords.at(1) ) * ipGlobalCoords.at(3);
#endif
        ipValues.addSubVectorRow(tempIPvalues,i+1,1);
        
        // set up virtual cell geometry for an qwedge
        auto nodes = giveFictiousNodeCoordsForExport(layer);
        FEInterpolation *interpol = static_cast< FEInterpolation * >( &this->interpolationForExport );
        FloatArray N;
        interpol->evalN( N, ip->giveNaturalCoordinates(), FEIVertexListGeometryWrapper( nodes ) ); 
        Nbar.addSubVectorRow(N,i+1,1);
    }
    // Nhat = inv(Nbar^T*Nbar)*Nbar^T
    NbarTNbar.beTProductOf(Nbar,Nbar);
    NbarTNbarInv.beInverseOf(NbarTNbar);
    Nhat.beProductTOf(NbarTNbarInv,Nbar);

    temprecovedValues.beProductOf(Nhat,ipValues);
    temprecovedValuesT.beTranspositionOf(temprecovedValues);
    ///TODO se över nodnumrering.
    IntArray strangeNodeNumbering = {2, 1, 3, 5, 4, 6, 7, 9, 8, 10, 12, 11, 13, 14, 15 };
    for (int i = 0; i < numNodes; i++ ) {
        if ( valueType == ISVT_TENSOR_S3 ) {
            FloatArray nodalStresses; 
            nodalStresses.beColumnOf(temprecovedValuesT,i+1);
            recoveredValues[strangeNodeNumbering[i]-1] = convV6ToV9Stress(nodalStresses);
        } else {
            recoveredValues[strangeNodeNumbering[i]-1].beColumnOf(temprecovedValuesT,i+1);
        }
    }
}

void 
Shell7Base :: giveL2contribution(FloatMatrix &ipValues, FloatMatrix &Nbar, int layer, InternalStateType type, TimeStep *tStep)
{ 
    // composite element interpolator
    FloatMatrix localNodeCoords;
    this->interpolationForExport.giveLocalNodeCoords(localNodeCoords);

    int numNodes = localNodeCoords.giveNumberOfColumns();
    
    std :: unique_ptr< IntegrationRule > &iRule = integrationRulesArray [ layer - 1 ];
    int numIP = iRule->giveNumberOfIntegrationPoints();
    FloatArray nodeCoords, ipCoords;
    InternalStateValueType valueType =  giveInternalStateValueType(type);
    
    // Find IP values and set up matrix of base functions
    // FloatMatrix Nbar, NbarTNbar, NbarTNbarInv, Nhat, ipValues, temprecovedValues, temprecovedValuesT; 
    Nbar.resize(numIP,numNodes);
    int numSC; 
    if ( valueType == ISVT_TENSOR_S3 ) { numSC = 6; } else { numSC = 9; };
    ipValues.resize(numIP,numSC);
    
    for ( int i = 0; i < numIP; i++ ) {
        auto *ip = iRule->getIntegrationPoint(i);
        FloatArray tempIPvalues;
        this->giveIPValue(tempIPvalues, ip, type, tStep);
        ipValues.addSubVectorRow(tempIPvalues,i+1,1);
        
        // set up virtual cell geometry for an qwedge
        auto nodes = giveFictiousNodeCoordsForExport(layer);
        FEInterpolation *interpol = static_cast< FEInterpolation * >( &this->interpolationForExport );
        FloatArray N;
        interpol->evalN( N, ip->giveNaturalCoordinates(), FEIVertexListGeometryWrapper( nodes ) ); 
        Nbar.addSubVectorRow(N,i+1,1);
    }
}

void
Shell7Base :: giveSPRcontribution(FloatMatrix &eltIPvalues, FloatMatrix &eltPolynomialValues, int layer, InternalStateType type, TimeStep *tStep)
{
    /* Returns the integration point values and polynomial values of the element interpolation wedge at layer, such that
       [GPvalue] = [P(x_GP,y_GP,z_GP)]*a, where 
       P(x,y,z) = [1 x y z yz xz xz x^2 y^2] (NB: z^2 term omitted)
       a = [a1 ... a9] is the vector of coefficients, which is to be found by the polynomial least square fit.
    */
    
    std :: unique_ptr< IntegrationRule > &iRule = this->integrationRulesArray[layer-1];                                     // Stämmer det här???
    IntegrationPoint *ip;
    int numEltIP = iRule->giveNumberOfIntegrationPoints();
    eltIPvalues.clear();
    eltPolynomialValues.clear();
    
    // Loop over IP:s in wedge interpolation
    for ( int iIP = 0; iIP < numEltIP; iIP++ ) {
        ip = iRule->getIntegrationPoint(iIP);
        
        // Collect IP-value
        FloatArray IPvalues;
        this->giveIPValue(IPvalues, ip, type, tStep);
        eltIPvalues.addSubVectorRow(IPvalues,iIP+1,1);
        
        // Collect global coordinates for IP and assemble to P.
        FloatArray IpCoords;
        IpCoords = ip->giveGlobalCoordinates();
        FloatArray iRowP = {1,IpCoords.at(1),IpCoords.at(2),IpCoords.at(3),
                            IpCoords.at(2)*IpCoords.at(3),IpCoords.at(1)*IpCoords.at(3),IpCoords.at(1)*IpCoords.at(2),
                            IpCoords.at(1)*IpCoords.at(1),IpCoords.at(2)*IpCoords.at(2)};
        eltPolynomialValues.addSubVectorRow(iRowP,iIP+1,1);
        
    }
}

void
Shell7Base :: giveTractionBC(FloatMatrix &tractionTop, FloatMatrix &tractionBtm, TimeStep *tStep)
{
    int numInPlaneIP = tractionTop.giveNumberOfColumns(); 
    tractionTop.zero(); tractionBtm.zero();
    
    int numLoads = this->boundaryLoadArray.giveSize() / 2;
    
    for ( int i = 1; i <= numLoads; i++ ) {     // For each pressure load that is applied
        
        int load_number = this->boundaryLoadArray.at(2 * i - 1);
        //int iSurf = this->boundaryLoadArray.at(2 * i);         // load_id
        Load *load = this->domain->giveLoad(load_number);
        
        if ( ConstantSurfaceLoad* sLoad = dynamic_cast< ConstantSurfaceLoad * >( load ) ) {
            
//             FloatArray tractionBC, tempTractionTop, tempTractionBtm; tempTractionTop.clear(); tempTractionBtm.clear();
            FloatArray tractionBC(3), tempBC;
            load->computeComponentArrayAt(tempBC,tStep,VM_Total);   ///TODO: update VM-mode?
            double xi = sLoad->giveLoadOffset();
            tractionBC.assemble(tempBC,sLoad->giveDofIDs());
            if ( tractionBC.giveSize() != tractionBtm.giveNumberOfRows() ) {
                OOFEM_ERROR("Number of stress components don't match");
            }
            
            // Assume linear variation of placement
            //tempTractionTop.add(((1.0+xi)/2.0),tractionBC);        // positive z-normal
            //tempTractionBtm.add((-(1.0-xi)/2.0),tractionBC);       // negative z-normal            
            for (int iIP = 1 ; iIP <= numInPlaneIP ; iIP++ ) {
                
                tractionTop.addSubVectorCol(((1.0+xi)/2.0)*tractionBC,1,iIP);
                tractionBtm.addSubVectorCol((-(1.0-xi)/2.0)*tractionBC,1,iIP);
            }
        } else {
            OOFEM_ERROR("Load type not supported");
        }
    }
    
//     if (this->giveGlobalNumber() == 48) {
//         tractionTop.printYourself("tractionTop");
//         tractionBtm.printYourself("tractionBtm");
//     }
}
  
void 
Shell7Base :: recoverShearStress(TimeStep *tStep)
{
    // Recover shear stresses at ip by integration of the momentum balance through the thickness    
    
    int numberOfLayers = this->layeredCS->giveNumberOfLayers();     // conversion of types
    
    int numThicknessIP = this->layeredCS->giveNumIntegrationPointsInLayer();        
    if (numThicknessIP < 2) {
        // Polynomial fit involves linear z-component at the moment
        OOFEM_ERROR("To few thickness IP per layer to do polynomial fit");
    }
    
    int numInPlaneIP = 6; ///TODO generalise this;
//     int numInPlaneIP = this->giveIntegrationRule()->giveNumberOfIntegrationPoints();
    int numWedgeIP = numInPlaneIP*numThicknessIP;
    double totalThickness = this->layeredCS->computeIntegralThick(); 
//     double integralThickness = 0.0;
    double zeroThicknessLevel = - 0.5 * totalThickness;     // assumes midplane is the geometric midplane of layered structure.
//     double zeroThicknessLevel = this->layeredCS->give( CS_BottomZCoord, &lCoords, this, false );  
//     FloatArray Sold;
    FloatMatrix dSmatLayer(3,numInPlaneIP), SmatOld(3,numInPlaneIP);  // 3 stress components (S_xz, S_yz, S_zz) * num of in plane ip 
    FloatMatrix dSmatLayerIP(3,numWedgeIP);
    
    // Fitting of Szz to traction BC
    std::vector <FloatMatrix> dSmatIP; dSmatIP.resize(numberOfLayers);              //recovered stress values in wedge IPs
    std::vector <FloatMatrix> dSmat; dSmat.resize(numberOfLayers);              //recovered stress values at layer top
    std::vector <FloatMatrix> dSmatIPupd; dSmatIPupd.resize(numberOfLayers);     //recovered stress values in wedge IPs fitted to  top and btm traction BC. 
    std::vector <FloatMatrix> dSmatupd; dSmatupd.resize(numberOfLayers);              //recovered stress values at delamination top, adjusted for top and btm traction BC
    
    // Find top and bottom BC NB: assumes constant over element surface.
    FloatMatrix tractionTop(3,numInPlaneIP), tractionBtm(3,numInPlaneIP); 
    giveTractionBC(tractionTop, tractionBtm, tStep);
    SmatOld = tractionBtm;        // integration performed from bottom
//     SmatOld.zero();    
//     for ( int i = 1; i <= numInPlaneIP ; i++) {
//         SmatOld.at(1,i) = tractionBtm.at(1);
//         SmatOld.at(2,i) = tractionBtm.at(2);
//         SmatOld.at(3,i) = tractionBtm.at(3);
//     }
    
    // Integration from the bottom
    for ( int layer = 1; layer <= numberOfLayers; layer++ ) {

        /* Recover values by a polynomial fit to the stress values in a patch of elements closest to the this element.
            * The vector of GPvalues is [GPvalue] = [P(x_GP,y_GP,z_GP)]*a, where 
            * P(x,y,z) = [1 x y z yz xz xz x^2 y^2] (NB: z^2 term omitted)
            * a = [a1 ... a9] is the vector of coefficients of P
            * a = inv(A)*b, where 
            * A = [P]^T*[P], b = [P]^T*[GPvalue], calculated over the appropriate patch.
            * the gradient of the GP value (used in the stress recovery) is then directly calculated by
            * d[GPvalue]/di = [dP/di|GP]*a, i = x,y,z.
            * dP/dx = [0 1 0 0 0 z y 2x 0]
            * dP/dy = [0 0 1 0 z 0 x 0 2y]
            * dP/dz = [0 0 0 1 y x 0 0 0 ]
            * which is then analytically integrated over z (see giveZintegratedPolynomialGradientForStressRecAt)
            * The recovery of the transverse normal stress in performed using the same analytical expression 
            * The values in the GP is overwritten by the recovery
            */
        
#if 0
        // Integration from bottom and adjusting Szz to top/btm traction BC. Shear stress not adjusted to top BC.
        giveLayerContributionToSR(dSmatLayer, dSmatLayerIP, layer, zeroThicknessLevel, tStep);
        
        
        // Save layer stresses (IP- and interface values) to be able to fit to traction BC
        // Only Szz fitted. ///TODO fit shear stresses to traction BC?
        dSmatIP[layer-1].beSubMatrixOf(dSmatLayerIP, 3, 3, 1, numWedgeIP);
        dSmat[layer-1].beSubMatrixOf(dSmatLayer, 3, 3, 1, numInPlaneIP);
        //if (this->giveGlobalNumber() == 48 ) { dSmatLayerIP.printYourself(); dSmatIP[layer -1].printYourself(); }
        
        updateLayerTransvShearStressesSR(dSmatLayerIP, SmatOld, layer);          
        
//         updateLayerTransvStressesSR(dSmatLayerIP, SmatOld, layer);

        SmatOld.add(dSmatLayer); // Integrated stress over laminate
#endif
        
        // Perform integration and distribute the integration error to top and btm 
        // Integration of Szz requires both top and btm traction BC
        // Integration of shear stresses only require one BC, however the error is distributed across thickness (ie the same as doing integration from top and btm and taking average)
        // fulfillment of top and btm traction BC optional for shear stresses.
        
        giveLayerContributionToSR(dSmat[layer-1], dSmatIP[layer-1], layer, zeroThicknessLevel, tStep);
        SmatOld.add(dSmat[layer-1]); // Integrated stress over laminate
        
        zeroThicknessLevel += this->layeredCS->giveLayerThickness(layer);      
        
//         integralThickness += thickness;
        
    } 
    
    // Fitting of stresses in IPs to traction BC
    zeroThicknessLevel = - 0.5 * totalThickness;
    
#if 0
    // Only Szz
    FloatMatrix integratedStress;
    integratedStress.beSubMatrixOf(SmatOld, 3, 3, 1, numInPlaneIP);
    FloatArray tempTractionBtm = {tractionBtm.at(3)};
    FloatArray tempTractionTop = {tractionTop.at(3)};
    fitRecoveredStress2BC(dSmatIPupd, dSmat, dSmatIP, integratedStress, tempTractionBtm, tempTractionTop, zeroThicknessLevel, {1});
//     fitRecoveredStress2BC(dSmatIPupd, dSmat, dSmatIP, integratedStress, {tractionBtm.at(3)}, {tractionBtm.at(3)}, zeroThicknessLevel, {1});
    
    FloatArray zeros(numInPlaneIP); zeros.zero();
    for ( int layer = 1; layer <= numberOfLayers; layer++ ) {
//         if (this->giveGlobalNumber() == 48 ) { dSzzmatIP[layer -1].printYourself(); }
        updateLayerTransvNormalStressSR( dSmatIPupd[layer -1], zeros, layer); // traction on bottom already taken into account in the integration. 
    }
#endif

    // All transverse stress components
    fitRecoveredStress2BC(dSmatIPupd, dSmatupd, dSmat, dSmatIP, SmatOld, tractionBtm, tractionTop, zeroThicknessLevel, {0.0,0.0,1.0}, 1, numberOfLayers);    // {0.0,0.0,1.0}: only Szz fulfills BC, shear stress integration error is only distributed to top and btm
    
    for ( int layer = 1; layer <= numberOfLayers; layer++ ) {
        //if (this->giveGlobalNumber() == 84 ) { dSmatIPupd[layer -1].printYourself(); }
        updateLayerTransvStressesSR( dSmatIPupd[layer -1], layer); 
    }
    
    
}

void 
Shell7Base :: giveRecoveredTransverseInterfaceStress(std::vector<FloatMatrix> &transverseStress, TimeStep *tStep)
{
    // Recover shear stresses at ip by integration of the momentum balance through the thickness  
    transverseStress.clear();
    
    int numberOfLayers = this->layeredCS->giveNumberOfLayers();     // conversion of types
    int numThicknessIP = this->layeredCS->giveNumIntegrationPointsInLayer();        
    if (numThicknessIP < 2) {
        // Polynomial fit involves linear z-component at the moment
        OOFEM_ERROR("To few thickness IP per layer to do polynomial fit");
    }
    
    int numInPlaneIP = 6; ///TODO generalise this;
    double totalThickness = this->layeredCS->computeIntegralThick(); 
    double zeroThicknessLevel = - 0.5 * totalThickness;     // assumes midplane is the geometric midplane of layered structure.
    
    FloatMatrix SmatOld(3,numInPlaneIP);  // 3 stress components (S_xz, S_yz, S_zz) * num of in plane ip 
    
    // Fitting of Szz to traction BC
    std::vector <FloatMatrix> dSmatIP; dSmatIP.resize(numberOfLayers);              //recovered stress values in wedge IPs
    std::vector <FloatMatrix> dSmat; dSmat.resize(numberOfLayers);              //recovered stress values at layer top
    std::vector <FloatMatrix> dSmatIPupd; dSmatIPupd.resize(numberOfLayers);     //recovered stress values in wedge IPs fitted to  top and btm traction BC. 
    std::vector <FloatMatrix> dSmatupd; dSmatupd.resize(numberOfLayers);              //recovered stress values at layer top, fitted to  top and btm traction BC
    
    // Find top and bottom BC NB: assumes constant over element surface.
    FloatMatrix tractionTop(3,numInPlaneIP), tractionBtm(3,numInPlaneIP); 
    giveTractionBC(tractionTop, tractionBtm, tStep);
    SmatOld = tractionBtm;        // integration performed from bottom
//     if (this->giveGlobalNumber() == 10) {
//         printf("Target time: %i \n",tStep->giveNumber());
//         printf("Intrinsic time: %f \n",tStep->giveIntrinsicTime());
//         printf("Target time: %f \n",tStep->giveTargetTime());
//         tractionTop.printYourself("tractionTop"); 
//         tractionBtm.printYourself("tractionBtm");
//     }  
    
    // Integration from the bottom
    for ( int layer = 1; layer <= numberOfLayers; layer++ ) {

        /* Recover values by a polynomial fit to the stress values in a patch of elements closest to the this element.
            * The vector of GPvalues is [GPvalue] = [P(x_GP,y_GP,z_GP)]*a, where 
            * P(x,y,z) = [1 x y z yz xz xz x^2 y^2] (NB: z^2 term omitted)
            * a = [a1 ... a9] is the vector of coefficients of P
            * a = inv(A)*b, where 
            * A = [P]^T*[P], b = [P]^T*[GPvalue], calculated over the appropriate patch.
            * the gradient of the GP value (used in the stress recovery) is then directly calculated by
            * d[GPvalue]/di = [dP/di|GP]*a, i = x,y,z.
            * dP/dx = [0 1 0 0 0 z y 2x 0]
            * dP/dy = [0 0 1 0 z 0 x 0 2y]
            * dP/dz = [0 0 0 1 y x 0 0 0 ]
            * which is then analytically integrated over z (see giveZintegratedPolynomialGradientForStressRecAt)
            * The recovery of the transverse normal stress in performed using the same analytical expression 
            * The values in the GP is overwritten by the recovery
            */
        
        // Perform integration and distribute the integration error to top and btm 
        // Integration of Szz requires both top and btm traction BC
        // Integration of shear stresses only require one BC, however the error is distributed across thickness (ie the same as doing integration from top and btm and taking average)
        // fulfillment of top and btm traction BC optional for shear stresses.
        
        giveLayerContributionToSR(dSmat[layer-1], dSmatIP[layer-1], layer, zeroThicknessLevel, tStep);
        SmatOld.add(dSmat[layer-1]); // Integrated stress over laminate
        
        zeroThicknessLevel += this->layeredCS->giveLayerThickness(layer);      
        
    } 
    
    // Fitting of stresses in IPs to traction BC
    zeroThicknessLevel = - 0.5 * totalThickness;

    // All transverse stress components
    transverseStress.resize(numberOfLayers-1);              //recovered stress values at layer interfaces 
    fitRecoveredStress2BC(dSmatIPupd, dSmatupd, dSmat, dSmatIP, SmatOld, tractionBtm, tractionTop, zeroThicknessLevel, {0.0,0.0,1.0}, 1, numberOfLayers);    // {0.0,0.0,1.0}: only Szz fulfills BC, shear stress integration error is only distributed to top and btm
    
    for ( int layer = 1 ; layer < numberOfLayers ; layer++ ) {
        transverseStress[layer-1] = dSmatupd[layer-1];
//         if ( this->giveGlobalNumber() == 61 ) {
//             transverseStress[layer-1].printYourself();
//         }
    }
}

void 
Shell7Base :: giveLayerContributionToSR(FloatMatrix &dSmatLayer, FloatMatrix &dSmatLayerIP, int layer, double zeroThicknessLevel, TimeStep *tStep)
{
    /* Recover values by a polynomial fit to the stress values in a patch of elements closest to the this element.
    * The vector of GPvalues is [GPvalue] = [P(x_GP,y_GP,z_GP)]*a, where 
    * P(x,y,z) = [1 x y z yz xz xz x^2 y^2] (NB: z^2 term omitted)
    * a = [a1 ... a9] is the vector of coefficients of P
    * a = inv(A)*b, where 
    * A = [P]^T*[P], b = [P]^T*[GPvalue], calculated over the appropriate patch.
    * the gradient of the GP value (used in the stress recovery) is then directly calculated by
    * d[GPvalue]/di = [dP/di|GP]*a, i = x,y,z.
    * dP/dx = [0 1 0 0 0 z y 2x 0]
    * dP/dy = [0 0 1 0 z 0 x 0 2y]
    * dP/dz = [0 0 0 1 y x 0 0 0 ]
    * which is then analytically integrated over z (see giveZintegratedPolynomialGradientForStressRecAt)
    * The recovery of the transverse normal stress in performed using the same analytical expression 
    */
    std :: unique_ptr< IntegrationRule > &iRuleL = this->integrationRulesArray [ layer - 1 ];   
       
    int numInPlaneIP = 6; ///TODO generalise this; 
    int numThicknessIP = this->layeredCS->giveNumIntegrationPointsInLayer();     
    dSmatLayer.clear(); dSmatLayer.resize(3,numInPlaneIP);                  // 3 stress components (S_xz, S_yz, S_zz) * num of in plane ip 
    dSmatLayerIP.clear(); dSmatLayerIP.resize(3,numInPlaneIP*numThicknessIP); // 3 stress components (S_xz, S_yz, S_zz) * num of wedge IP
    
    //recoveredValues.resize(numWedgeNodes);
    IntArray centreElNum = {this->giveGlobalNumber()};
    IntArray centreElNodes = this->giveDofManArray(); 
    IntArray patchEls; 
    int numSampledIP = 0, numCoefficents = 11;
//         int numPatchEls = patchEls.giveSize();
//         FloatMatrix patchIpValues(numPatchEls*numThicknessIP*numInPlaneIP,6), P(numPatchEls*numThicknessIP*numInPlaneIP,numCoefficents);
    FloatMatrix patchIpValues, P;
    patchIpValues.zero();
    P.zero();
    
    Domain *d = this->giveDomain();
    d->giveConnectivityTable()->giveElementNeighbourList(patchEls,centreElNum); 
    
    for (int patchElNum : patchEls) { 
        
        
        // Get (pointer to) current element in patch and calculate its addition to the patch
        Shell7Base *patchEl = static_cast<Shell7Base*>(d->giveElement(patchElNum));
        
        std :: unique_ptr< IntegrationRule > &iRule = patchEl->integrationRulesArray[layer-1];                                     // Stämmer det här???
        IntegrationPoint *ip;
        int numEltIP = iRule->giveNumberOfIntegrationPoints();
        
        // Loop over IP:s in wedge interpolation
        for ( int iIP = 0; iIP < numEltIP; iIP++ ) {
            ip = iRule->getIntegrationPoint(iIP);    
//                 printf("blajjja \n");
            
            // Collect IP-value
            FloatArray tempIPvalues;
            patchEl->giveIPValue(tempIPvalues, ip, IST_StressTensor, tStep);
            patchIpValues.addSubVectorRow(tempIPvalues,numSampledIP + iIP+1,1);
            
            // Collect global coordinates for IP and assemble to P.
            FloatArray IpCoords;
            IpCoords = ip->giveGlobalCoordinates();
//             IpCoords.at(3) -= zeroThicknessLevel;    // make bottom of layer ref point for polynimial fit, ie using layer z-coords (not shell)    
/*            FloatArray iRowP = {1,IpCoords.at(1),IpCoords.at(2),IpCoords.at(3),
                                IpCoords.at(2)*IpCoords.at(3),IpCoords.at(1)*IpCoords.at(3),IpCoords.at(1)*IpCoords.at(2),
                                IpCoords.at(1)*IpCoords.at(1),IpCoords.at(2)*IpCoords.at(2)};  */ 
            // P = [1 x y z yz xz xy x² y² x²y x²z]
            FloatArray iRowP = {1.0,IpCoords.at(1),IpCoords.at(2),IpCoords.at(3),
                                IpCoords.at(2)*IpCoords.at(3),IpCoords.at(1)*IpCoords.at(3),IpCoords.at(1)*IpCoords.at(2),
                                IpCoords.at(1)*IpCoords.at(1),IpCoords.at(2)*IpCoords.at(2),
                                IpCoords.at(1)*IpCoords.at(1)*IpCoords.at(3),IpCoords.at(2)*IpCoords.at(2)*IpCoords.at(3)}; 
            // P = [x y yz xz xy x² y² x²y x²z]
//             FloatArray iRowP = {IpCoords.at(1),IpCoords.at(2),
//                                 IpCoords.at(2)*IpCoords.at(3),IpCoords.at(1)*IpCoords.at(3),IpCoords.at(1)*IpCoords.at(2),
//                                 IpCoords.at(1)*IpCoords.at(1),IpCoords.at(2)*IpCoords.at(2),
//                                 IpCoords.at(1)*IpCoords.at(1)*IpCoords.at(3),IpCoords.at(2)*IpCoords.at(2)*IpCoords.at(3)};
            // P = [x y yz xz xy x² y² x³ y³ xyz x²z y²z]
//             FloatArray iRowP = {IpCoords.at(1),IpCoords.at(2),
//                                 IpCoords.at(2)*IpCoords.at(3),IpCoords.at(1)*IpCoords.at(3),IpCoords.at(1)*IpCoords.at(2),
//                                 IpCoords.at(1)*IpCoords.at(1),IpCoords.at(2)*IpCoords.at(2),
//                                 IpCoords.at(1)*IpCoords.at(1)*IpCoords.at(1),IpCoords.at(2)*IpCoords.at(2)*IpCoords.at(2),IpCoords.at(1)*IpCoords.at(2)*IpCoords.at(3),
//                                 IpCoords.at(1)*IpCoords.at(1)*IpCoords.at(3),IpCoords.at(2)*IpCoords.at(2)*IpCoords.at(3)};
            P.addSubVectorRow(iRowP,numSampledIP + iIP+1,1);
            
        }
        
#if 0       
        // Replace above with function
        // Get contribution from element in patch
        FloatMatrix eltIPvalues, eltPolynomialValues;
        patchEl->giveSPRcontribution(eltIPvalues,eltPolynomialValues,layer,IST_StressTensor,tStep);
        int numEltIP = eltIPvalues.giveNumberOfRows();
        
        // Add to patch
        P.setSubMatrix(eltPolynomialValues,numSampledIP+1,1);
        patchIpValues.setSubMatrix(eltIPvalues,numSampledIP+1,1);
#endif

        // increase number of sampled IPs 
        numSampledIP += numEltIP;
    }
    
    if (numSampledIP < numCoefficents*numThicknessIP) {
        // Polynomial fit involves quadratic x- and y-components
        OOFEM_ERROR("To few in-plane IP to do polynomial fit");
    }
    
    // Find polynomial coefficients for patch
    FloatMatrix A, invA, Abar, a;
    // A = P^T*P
    A.beTProductOf(P,P);
    //A.writeCSV("Amatrix");
    invA.beInverseOf(A);
    // Abar = inv(A)*P^T
    Abar.beProductTOf(invA,P);
    // a = Abar*ipValues, a.size = 11 x StressComponents (11 = number of coefficients)
    a.beProductOf(Abar,patchIpValues);
    
    
    // Assemble appropriate coefficients matrix for computing S_xz and S_yz and S_zz
    // aSiz = [a(:,1)^T a(:,6)^T; a(:,6)^T a(:,2)^T]
    // aSzz = [a(:,1)^T 2*a(:,6)^T a(:,2)^T]
    // aSzz = sum[a81 a76 a82]
    FloatMatrix aSiz, aSzz; 
    aSiz.resize(2,numCoefficents*2);
    aSzz.resize(1,numCoefficents*3);
    for (int iCol = 1; iCol <= numCoefficents; iCol++) {
        aSiz.at(1,iCol) = a.at(iCol,1);
        aSiz.at(1,iCol+numCoefficents) = a.at(iCol,6);
        aSiz.at(2,iCol) = a.at(iCol,6);
        aSiz.at(2,iCol+numCoefficents) = a.at(iCol,2);
        aSzz.at(1,iCol) = a.at(iCol,1);
        aSzz.at(1,iCol+numCoefficents) = 2.0*a.at(iCol,6);
        aSzz.at(1,iCol+2*numCoefficents) = a.at(iCol,2);
    }
//     double aSzz = a.at(8,1) + a.at(7,6) + a.at(8,2); 
    
    
    // Integrate gradient of S_xz, S_yz and S_zz in all IP-stacks
    dSmatLayer.zero(); 
    dSmatLayerIP.zero();
    double thickness = this->layeredCS->giveLayerThickness(layer);
    
    // Compute recovered stresses
    for ( int j = 0; j < numInPlaneIP; j++ ) { 
        FloatArray dS, GPcoords(3), dSzz;
        
        // Approximation at bottom of layer (to fulfill integration BC)
        FloatArray dSBtm, dSzzBtm;
        GaussPoint *gp = iRuleL->getIntegrationPoint(j);
        GPcoords.zero();
        GPcoords = gp->giveGlobalCoordinates();
        GPcoords.at(3) = zeroThicknessLevel;
        // S_xz and S_yz
        dSBtm.zero();
        FloatArray intGradPBtm;
        this->giveZintegratedPolynomialGradientForStressRecAt(intGradPBtm,GPcoords);
        dSBtm.beProductOf(aSiz,intGradPBtm);
        // S_zz and dSzz/dz
        dSzzBtm.zero();
        FloatArray int2Grad2PBtm;
        this->giveZ2integratedPolynomial2GradientForStressRecAt(int2Grad2PBtm,GPcoords);
        dSzzBtm.beProductOf(aSzz,int2Grad2PBtm*(-1.0));
        
        // thickness GPs
        for ( int i = 0; i < numThicknessIP; i++ ) {

            int point = i*numInPlaneIP + j; // wedge integration point number
            
            GaussPoint *gp = iRuleL->getIntegrationPoint(point);
            GPcoords.zero();
            GPcoords = gp->giveGlobalCoordinates();
            #if 0
                if (this->giveGlobalNumber() == 126) {
                    if (layer == 1) {
                        if (i == 0) {
                            GPcoords.printYourself("wedge GPs");
                        }
                    }
                }
            #endif
//             GPcoords.at(3) -= zeroThicknessLevel;                      // Using layer z-coords (not shell)
            
            // Calculate z-integration of fitted stress variation in position of GP such that 
            // Siz = - Integral( dS_ij/dx + dS_ij/dy )dz = - [a(ij) a(ij)]*Integral( dP/dx dP/dy )dz = - [a(ij) a(ij)]*IntGradP + Sij^k-1
            
            // S_xz and S_yz
            dS.zero();
            FloatArray intGradP;
            this->giveZintegratedPolynomialGradientForStressRecAt(intGradP,GPcoords);
            dS.beProductOf(aSiz,intGradP*(-1.0));
//             dSmatLayerIP.at(1,point+1) = dS.at(1);                             // S_xz
//             dSmatLayerIP.at(2,point+1) = dS.at(2);                             // S_yz
            dSmatLayerIP.at(1,point+1) = dS.at(1) + dSBtm.at(1);                             // S_xz Using polynomial fit to shell coordinates
            dSmatLayerIP.at(2,point+1) = dS.at(2) + dSBtm.at(2);                             // S_yz Using polynomial fit to shell coordinates
            
            // S_zz
            dSzz.zero();
            FloatArray int2Grad2P;
            this->giveZ2integratedPolynomial2GradientForStressRecAt(int2Grad2P,GPcoords);
            dSzz.beProductOf(aSzz,int2Grad2P);
            // Using only polynomial fit
    //         dSmatLayerIP.at(3,point+1) = dSzz.at(1);   
            // Using polynomial fit to shell coordinates
            dSmatLayerIP.at(3,point+1) = dSzz.at(1) + dSzzBtm.at(1);   
            // Old polynomial fit
    //         dSmatLayerIP.at(3,point+1) = aSzz*GPcoords.at(3)*GPcoords.at(3);  

        }
        
        // Calculate stresses at upper interface of layer (use the x-y-coords of the upper GP of the layer)
        
//         GPcoords.at(3) = thickness;                      // Using layer z-coords (not shell)
        GPcoords.at(3) = thickness + zeroThicknessLevel;
        
        // S_xz and S_yz
        dS.zero();
        FloatArray intGradP;
        this->giveZintegratedPolynomialGradientForStressRecAt(intGradP,GPcoords);
        dS.beProductOf(aSiz,intGradP*(-1.0));
//         dSmatLayer.at(1,j+1) = dS.at(1);      // S_xz
//         dSmatLayer.at(2,j+1) = dS.at(2);      // S_yz
        dSmatLayer.at(1,j+1) = dS.at(1) + dSBtm.at(1);      // S_xz
        dSmatLayer.at(2,j+1) = dS.at(2) + dSBtm.at(2);      // S_yz
        
        // S_zz
        dSzz.zero();
        FloatArray int2Grad2P;
        this->giveZ2integratedPolynomial2GradientForStressRecAt(int2Grad2P,GPcoords);
        dSzz.beProductOf(aSzz,int2Grad2P);
        // Using only polynomial fit
//         dSmatLayer.at(3,j+1) = dSzz.at(1);  
        // Using polynomial fit to shell coordinates
        dSmatLayer.at(3,j+1) = dSzz.at(1) + dSzzBtm.at(1);     
        // Old polynomial fit
//         dSmatLayer.at(3,j+1) = aSzz*GPcoords.at(3)*GPcoords.at(3);                
        
    }
    
# if 0
    // Fit P2(x,y) = [x y xy]*a to S_xz and S_xy in element layer and add to integration of Szz. 
    int numCoefficents2 = 2;
    FloatMatrix eltIpValues, P2;
    eltIpValues.zero();
    P2.zero();
    
    // Loop over IP:s in wedge interpolation
    IntegrationPoint *ip;
    int numEltIP = iRuleL->giveNumberOfIntegrationPoints();
    for ( int iIP = 0; iIP < numEltIP; iIP++ ) {
        ip = iRuleL->getIntegrationPoint(iIP);    
        
        // Collect IP-value
        FloatArray tempIPvalues;
        this->giveIPValue(tempIPvalues, ip, IST_StressTensor, tStep);
        FloatArray tempSaz(2); tempSaz.beSubArrayOf(tempIPvalues,{4,5});     // extract Syz and Sxz
        eltIpValues.addSubVectorRow(tempSaz,iIP+1,1);
        
        // Collect global coordinates for IP and assemble to P.
        FloatArray IpCoords;
        IpCoords = ip->giveGlobalCoordinates();
        IpCoords.at(3) -= zeroThicknessLevel;    // make bottom of layer ref point for polynimial fit. 
        // P = [x y]
        FloatArray iRowP = {IpCoords.at(1),IpCoords.at(2)}; 
        // P = [x y xy]
//         FloatArray iRowP = {IpCoords.at(1),IpCoords.at(2),IpCoords.at(1)*IpCoords.at(2)};
        P2.addSubVectorRow(iRowP,iIP+1,1);
        
    }
    
    // Find polynomial coefficients for patch
    FloatMatrix A2, invA2, A2bar, a2;
    // A2 = P2^T*P2
    A2.beTProductOf(P2,P2);
    invA2.beInverseOf(A2);
    // A2bar = inv(A2)*P2^T
    A2bar.beProductTOf(invA2,P2);
    // a2 = A2bar*ipValues, a2.size = number of coefficients x 2
    a2.beProductOf(A2bar,eltIpValues);
    
    // Assemble appropriate coefficients matrix for computing S_zz
    // a2Szz = [a2(:,1)^T a2(:,2)^T]
    FloatMatrix a2Szz; 
    a2Szz.resize(1,numCoefficents2*2);
    FloatArray temp; 
    temp.beColumnOf(a2,1);
    a2Szz.addSubVectorRow(temp,1,1);
    temp.beColumnOf(a2,2);
    a2Szz.addSubVectorRow(temp,1,numCoefficents2+1);
//     if ( this->giveGlobalNumber() == 50) {
//         a2.printYourself("a2");
//         a2Szz.printYourself("a2Szz");
//     }   
    
    for ( int j = 0; j < numInPlaneIP; j++ ) { 
        
        FloatArray GPcoords(3), dSzz2;
        
        // thickness GPs
        for ( int i = 0; i < numThicknessIP; i++ ) {

            int point = i*numInPlaneIP + j; // wedge integration point number
            
            GaussPoint *gp = iRuleL->getIntegrationPoint(point);
            GPcoords.zero();
            GPcoords = gp->giveGlobalCoordinates();
            GPcoords.at(3) -= zeroThicknessLevel;
            
            dSzz2.zero();
            // P = [x y]
            FloatArray intGradP2 = {GPcoords.at(3),0,
                                    0,GPcoords.at(3)};
            // P = [x y xy]
    //         FloatArray intGradP2 = {GPcoords.at(3),0,GPcoords.at(2)*GPcoords.at(3),
    //                                 0,GPcoords.at(3),GPcoords.at(1)*GPcoords.at(3)};
            dSzz2.beProductOf(a2Szz,intGradP2*(-1.0));
            dSmatLayerIP.at(3,point+1) += dSzz2.at(1);  
        }
        
        GPcoords.at(3) = thickness;
        dSzz2.zero();
        // P = [x y]
        FloatArray intGradP2 = {GPcoords.at(3),0,
                                0,GPcoords.at(3)};
        // P = [x y xy]
//         FloatArray intGradP2 = {GPcoords.at(3),0,GPcoords.at(2)*GPcoords.at(3),
//                                 0,GPcoords.at(3),GPcoords.at(1)*GPcoords.at(3)};
        dSzz2.beProductOf(a2Szz,intGradP2*(-1.0));
        if ( this->giveGlobalNumber() == 50) {
            a2Szz.printYourself("a2Szz");
            intGradP2.printYourself("intGradP2");
            dSzz2.printYourself("dSzz2");
        } 
        dSmatLayer.at(3,j+1) += dSzz2.at(1);   
    }
# endif
    
}

void 
Shell7Base :: fitRecoveredStress2BC(std::vector<FloatMatrix> &answer1, std::vector<FloatMatrix> &answer2, std::vector<FloatMatrix> &dSmat, std::vector<FloatMatrix> &dSmatIP, FloatMatrix &SmatOld, FloatMatrix &tractionBtm, FloatMatrix &tractionTop, double zeroThicknessLevel, FloatArray fulfillBC, int startLayer, int endLayer)
{
    
    // Adjust the integrated stress values to take top and bottom surface BC into account
    // Optional fulfillment of shear traction BC (default is distribution of integration error across thickness)
    // NB: recovery of normal stress (Szz) requires both top and bottom BC fulfillment.
    
    int numInPlaneIP = 6; ///TODO generalise this.
    int numberOfLayers = startLayer - endLayer + 1;
    int numThicknessIP = this->layeredCS->giveNumIntegrationPointsInLayer();    
    int numSC = SmatOld.giveNumberOfRows();    // assume same number of stress component to be fitted in each layer
    double totalThickness = 0.0;
    if ( numberOfLayers == this->layeredCS->giveNumberOfLayers() ) {
        totalThickness = this->layeredCS->computeIntegralThick(); 
    } else {
        for ( int layer = startLayer ; layer <= endLayer; layer++ ) {
            totalThickness += this->layeredCS->giveLayerThickness(layer);
        }
    }
    
    if (numSC != tractionBtm.giveNumberOfRows() || numSC != tractionTop.giveNumberOfRows() ) {
        OOFEM_ERROR("Number of stress components don't match");
    }
    
    // Find integration error and correction term
    // Assumes constant integration error across thickness
    FloatMatrix C1(numSC,numInPlaneIP);
    FloatMatrix intError(numSC,numInPlaneIP);
    FloatMatrix SOld(numSC,numInPlaneIP);
    for ( int iSC = 1; iSC <= numSC; iSC++) {
        
        for ( int j = 0 ; j < numInPlaneIP ; j++ ) { 
        
            intError.at(iSC,j+1) = tractionTop.at(iSC,j+1) - SmatOld.at(iSC,j+1);
            C1.at(iSC,j+1) = (1.0/totalThickness)*intError.at(iSC,j+1);
            SOld.at(iSC,j+1) = tractionBtm.at(iSC,j+1) - 0.5*(1.0 - fulfillBC.at(iSC))*intError.at(iSC,j+1);    // last term distributes integration error on top and btm if fulfillBC = 0;
            
        }
    }
            
    FloatArray GPcoords;
        
    for ( int layer = startLayer ; layer <= endLayer; layer++ ) {
        
        std :: unique_ptr< IntegrationRule > &iRuleL = this->integrationRulesArray [ layer - 1 ];  
        answer1[layer-startLayer].resize(numSC,numInPlaneIP*numThicknessIP);
        answer2[layer-startLayer].resize(numSC,numInPlaneIP);
        
        for ( int iSC = 1; iSC <= numSC; iSC++) {
        
            for ( int j = 0 ; j < numInPlaneIP ; j++ ) { 
                
//                 if ( layer == 1 ) {
//                     SOld.at(iSC,j+1) = tractionBtm.at(iSC,j+1);
//                 }
            
                // Assume constant integration error across thickness
                //double C1 = (1.0/totalThickness)*(tractionTop.at(iSC,j+1) - SmatOld.at(iSC,j+1));
                
                for ( int i = 0 ; i < numThicknessIP ; i++ ) {
                    
                    int point = i*numInPlaneIP + j; // wedge integration point number
                
                    GaussPoint *gp = iRuleL->getIntegrationPoint(point);
                    GPcoords.zero();
                    GPcoords = gp->giveGlobalCoordinates();
                    GPcoords.at(3) -= zeroThicknessLevel; 
                    
                    // Stress in wedge GP:s
                    answer1[layer-startLayer].at(iSC,point+1) = dSmatIP[layer-startLayer].at(iSC,point+1) + C1.at(iSC,j+1)*GPcoords.at(3) + SOld.at(iSC,j+1);
                }
                
                GPcoords.at(3) = this->layeredCS->giveLayerThickness(layer);
                SOld.at(iSC,j+1) += dSmat[layer-startLayer].at(iSC,j+1) + C1.at(iSC,j+1)*GPcoords.at(3);
                
            }  
            // Stress in interface GP:s
            answer2[layer-startLayer] = SOld;
        }
        zeroThicknessLevel += this->layeredCS->giveLayerThickness(layer);  
    }
}

void
Shell7Base :: updateLayerTransvStressesSR(FloatMatrix &dSmatLayerIP, int layer)
{
    // add stresses from lower interface AND replace stresses in wedge GP
    
    std :: unique_ptr< IntegrationRule > &iRuleL = this->integrationRulesArray [ layer - 1 ];  
    
    int numInPlaneIP = 6; ///TODO generalise this; 
    int numThicknessIP = this->layeredCS->giveNumIntegrationPointsInLayer();         
    
    for ( int j = 0; j < numInPlaneIP; j++ ) { 
        for ( int i = 0; i < numThicknessIP; i++ ) {
            int point = i*numInPlaneIP + j; // wedge integration point number
            //dSmatLayerIP.at(1,point+1) += SmatOld.at(1,j+1);  // S_xz
            //dSmatLayerIP.at(2,point+1) += SmatOld.at(2,j+1);  // S_yz
            //dSmatLayerIP.at(3,point+1) += SmatOld.at(3,j+1);  // S_zz
            
            // Replace stresses
            GaussPoint *gp = iRuleL->getIntegrationPoint(point);
            StructuralMaterialStatus* status = dynamic_cast< StructuralMaterialStatus* > ( gp->giveMaterialStatus() );
            FloatArray Sold = status->giveStressVector(); 
            
            Sold.at(5) = dSmatLayerIP.at(1,point+1); // S_xz
            Sold.at(4) = dSmatLayerIP.at(2,point+1); // S_yz
            Sold.at(3) = dSmatLayerIP.at(3,point+1); // S_zz        
            if ( Sold.giveSize() > 6 ) {
                Sold.at(8) = Sold.at(5); // S_xz
                Sold.at(7) = Sold.at(4); // S_yz
            }
            status->letStressVectorBe(Sold);
        }
    } 
}

# if 1
void
Shell7Base :: updateLayerTransvShearStressesSR(FloatMatrix &dSmatLayerIP, FloatMatrix &SmatOld, int layer)
{
    // add stresses from lower interface AND replace stresses in wedge GP
    
    std :: unique_ptr< IntegrationRule > &iRuleL = this->integrationRulesArray [ layer - 1 ];  
    
    int numInPlaneIP = 6; ///TODO generalise this; 
    int numThicknessIP = this->layeredCS->giveNumIntegrationPointsInLayer();         
    
    for ( int j = 0; j < numInPlaneIP; j++ ) { 
        for ( int i = 0; i < numThicknessIP; i++ ) {
            int point = i*numInPlaneIP + j; // wedge integration point number
            //dSmatLayerIP.at(1,point+1) += SmatOld.at(1,j+1);  // S_xz
            //dSmatLayerIP.at(2,point+1) += SmatOld.at(2,j+1);  // S_yz
            
            // Replace stresses
            GaussPoint *gp = iRuleL->getIntegrationPoint(point);
            StructuralMaterialStatus* status = dynamic_cast< StructuralMaterialStatus* > ( gp->giveMaterialStatus() );
            FloatArray Sold = status->giveStressVector(); 
            #if 0
                if (this->giveGlobalNumber() == 126) {
                    if (j == 3) {
                        if (i == 1) {
                            Sold.printYourself("Original");
                        }
                    }
                }
            #endif
            Sold.at(5) = dSmatLayerIP.at(1,point+1); // S_xz
            Sold.at(4) = dSmatLayerIP.at(2,point+1); // S_yz   
            if ( Sold.giveSize() > 6 ) {
                Sold.at(8) = Sold.at(5); // S_xz
                Sold.at(7) = Sold.at(4); // S_yz
            }
            status->letStressVectorBe(Sold);
        }
    } 
}

void
Shell7Base :: updateLayerTransvNormalStressSR(FloatMatrix &dSzzMatLayerIP, FloatArray &SzzMatOld, int layer)
{
    // add stresses from lower interface AND replace stresses in wedge GP
    
    std :: unique_ptr< IntegrationRule > &iRuleL = this->integrationRulesArray [ layer - 1 ];  
    
    int numInPlaneIP = 6; ///TODO generalise this;  
    int numThicknessIP = this->layeredCS->giveNumIntegrationPointsInLayer();         
    
    for ( int j = 0; j < numInPlaneIP; j++ ) { 
        for ( int i = 0; i < numThicknessIP; i++ ) {
            int point = i*numInPlaneIP + j; // wedge integration point number
            dSzzMatLayerIP.at(1,point+1) += SzzMatOld.at(j+1);  // S_zz
            
            // Replace stresses
            GaussPoint *gp = iRuleL->getIntegrationPoint(point);
            StructuralMaterialStatus* status = dynamic_cast< StructuralMaterialStatus* > ( gp->giveMaterialStatus() );
            FloatArray Sold = status->giveStressVector(); 
            
            Sold.at(3) = dSzzMatLayerIP.at(1,point+1); // S_zz        
            status->letStressVectorBe(Sold);
        }
    } 
}
# endif

#if 0
void
Shell7Base :: givePolynomialGradientForStressRecAt(FloatArray &answer, FloatArray &coords)
{
    /* Return the special matrix {gradP} of the reciever evaluated at a global coordinate (of a Gauss point)
     * Used for integration of S_xz and S_yz such that
     * a(ij)*gradP = dS_ij/dx + dS_ij/dy, where a(ij) is a vector of coefficients to the polynomial
     * P(x,y,z) = [1 x y z yz xz xy x² y² x²z y²z] (NB: z^2 term omitted)
     * associated with stress component ij.
     * The gradient of P(x,y,z) is given by
     * dP/dx = [0 1 0 0 0 z y 2x 0]
     * dP/dy = [0 0 1 0 z 0 x 0 2y]
     * dP/dz = [0 0 0 1 y x 0 0 0 ] (NB not returned atm)
     * answer = gradP = [dP/dx dP/dy]^T
     */
    
    answer.zero();
    answer = {0,1,0,0,           0,coords.at(3),coords.at(2),2.0*coords.at(1),             0,
              0,0,1,0,coords.at(3),           0,coords.at(1),             0,2.0*coords.at(2)};
}
#endif


void
Shell7Base :: giveZintegratedPolynomialGradientForStressRecAt(FloatArray &answer, FloatArray &coords)
{
    /* Return the special matrix {gradP} of the reciever evaluated at a global coordinate (of a Gauss point)
     * Used for integration of S_xz and S_yz such that
     * a(ij)*gradP = dS_ij/dx + dS_ij/dy, where a(ij) is a vector of coefficients to the polynomial
     * P(x,y,z) = [1 x y z yz xz xy x² y² x²z y²z] (NB: z^2 term omitted)
     * associated with stress component ij.
     * The gradient of P(x,y,z) is given by
     * dP/dx = [0 1 0 0 0 z y 2x 0 2xz 0]
     * dP/dy = [0 0 1 0 z 0 x 0 2y 0 2yz]
     * and the z-integration is then:
     * I[dP/dx]dz = [0 z 0 0 0    z²/2 yz 2xz 0   xz² 0  ]
     * I[dP/dy]dz = [0 0 z 0 z²/2 0    xz 0   2yz 0   yz²]
     * answer = IntgradP = [I[dP/dx]dz I[dP/dy]yz]^T
     */
    
    answer.zero();
    // P = [1 x y z yz xz xy x² y²]
//     answer = {0,coords.at(3),           0,0,                                  0,0.5*coords.at(3)*coords.at(3),coords.at(2)*coords.at(3),2.0*coords.at(1)*coords.at(3),                            0,
//               0,           0,coords.at(3),0,0.5*coords.at(3)*coords.at(3),                                  0,coords.at(1)*coords.at(3),                            0,2.0*coords.at(2)*coords.at(3)};
    
    // P = [1 x y z yz xz xy x² y² x²z y²z]
    answer = {0,coords.at(3),0,0,0,0.5*coords.at(3)*coords.at(3),coords.at(2)*coords.at(3),2.0*coords.at(1)*coords.at(3),0,coords.at(1)*coords.at(3)*coords.at(3),0,
              0,0,coords.at(3),0,0.5*coords.at(3)*coords.at(3),0,coords.at(1)*coords.at(3),0,2.0*coords.at(2)*coords.at(3),0,coords.at(2)*coords.at(3)*coords.at(3)};
    
    // P = [x y yz xz xy x² y² x²z y²z]
//     answer = {coords.at(3),0,0,0.5*coords.at(3)*coords.at(3),coords.at(2)*coords.at(3),2.0*coords.at(1)*coords.at(3),0,coords.at(1)*coords.at(3)*coords.at(3),0,
//               0,coords.at(3),0.5*coords.at(3)*coords.at(3),0,coords.at(1)*coords.at(3),0,2.0*coords.at(2)*coords.at(3),0,coords.at(2)*coords.at(3)*coords.at(3)};
    
    // P = [x y yz xz xy x² y² x²z y²z xyz]
//     answer = {coords.at(3),0,0,0.5*coords.at(3)*coords.at(3),coords.at(2)*coords.at(3),2.0*coords.at(1)*coords.at(3),0,coords.at(1)*coords.at(3)*coords.at(3),0,0.5*coords.at(2)*coords.at(3)*coords.at(3),
//               0,coords.at(3),0.5*coords.at(3)*coords.at(3),0,coords.at(1)*coords.at(3),0,2.0*coords.at(2)*coords.at(3),0,coords.at(2)*coords.at(3)*coords.at(3),0.5*coords.at(1)*coords.at(3)*coords.at(3)};
        
    // P = [x y yz xz xy x² y² x³ y³ xyz x²z y²z]
//     answer = {coords.at(3),0,0,0.5*coords.at(3)*coords.at(3),coords.at(2)*coords.at(3),2.0*coords.at(1)*coords.at(3),0,3.0*coords.at(1)*coords.at(1)*coords.at(3),0,0.5*coords.at(2)*coords.at(3)*coords.at(3),coords.at(1)*coords.at(3)*coords.at(3),0,
//               0,coords.at(3),0.5*coords.at(3)*coords.at(3),0,coords.at(1)*coords.at(3),0,2.0*coords.at(2)*coords.at(3),0,3.0*coords.at(2)*coords.at(2)*coords.at(3),0.5*coords.at(1)*coords.at(3)*coords.at(3),0,coords.at(2)*coords.at(3)*coords.at(3)};

    
}

void
Shell7Base :: giveZ2integratedPolynomial2GradientForStressRecAt(FloatArray &answer, FloatArray &coords)
{
    /* 
     * I[d(I[dP/dx]dz)/dx]dz = [0 0 0 0 0 0 0 z² 0  z³/3 0   ]
     * I[d(I[dP/dy]dz)/dy]dz = [0 0 0 0 0 0 0 0  z² 0    z³/3]
     * I[d(I[dP/dx]dz)/dy]dz = I[d(I[dP/dy]dz)/dx]dz = [0 0 0 0 0 0 z²/2 0 0 0 0]
     * answer = int2grad2P = [I[d(I[dP/dx]dz)/dx]dz I[d(I[dP/dx]dz)/dy]dz I[d(I[dP/dy]dz)/dy]dz]^T
     */
    
    answer.zero();
    
    // P = [1 x y z yz xz xy x² y² x²z y²z]
    answer = {0,0,0,0,0,0,0,coords.at(3)*coords.at(3),0,(1.0/3.0)*coords.at(3)*coords.at(3)*coords.at(3),0,
              0,0,0,0,0,0,0.5*coords.at(3)*coords.at(3),0,0,0,0,
              0,0,0,0,0,0,0,0,coords.at(3)*coords.at(3),0,(1.0/3.0)*coords.at(3)*coords.at(3)*coords.at(3)};
    
    // P = [x y yz xz xy x² y² x²z y²z]
//     answer = {0,0,0,0,0,coords.at(3)*coords.at(3),0,(1.0/3.0)*coords.at(3)*coords.at(3)*coords.at(3),0,
//               0,0,0,0,0.5*coords.at(3)*coords.at(3),0,0,0,0,
//               0,0,0,0,0,0,coords.at(3)*coords.at(3),0,(1.0/3.0)*coords.at(3)*coords.at(3)*coords.at(3)};
//               
    // P = [x y yz xz xy x² y² x²z y²z]
//     answer = {0,0,0,0,0,coords.at(3)*coords.at(3),0,3.0*coords.at(1)*coords.at(3)*coords.at(3),0,0,(1.0/3.0)*coords.at(3)*coords.at(3)*coords.at(3),0,
//               0,0,0,0,0.5*coords.at(3)*coords.at(3),0,0,0,0,(1.0/6.0)*coords.at(3)*coords.at(3)*coords.at(3),0,0,
//               0,0,0,0,0,0,coords.at(3)*coords.at(3),0,3.0*coords.at(2)*coords.at(3)*coords.at(3),0,0,(1.0/3.0)*coords.at(3)*coords.at(3)*coords.at(3)};
              
}

#if 0
void
Shell7Base :: computeBmatrixForStressRecAt(FloatArray &lcoords, FloatMatrix &answer, int layer, bool intSzz)
{
    /* Returns the  special matrix {B} of the receiver, evaluated at aGaussPoint. 
     * For integration of S_xz and S_yz such that
     * B*a = [dS_xx/dx + dS_xy/dy, dS_yx/dx + dS_yy/dy ]^T, where a is the vector of in plane 
     * stresses a = [S1_xx, S1_yy, S1_xy, ..., Sn_xx, Sn_yy, Sn_xy] for n layer nodes
     * For integration of S_zz such that
     * B*a = [dS_yz/dx + dS_yz/dy ]^T, where a is the vector of in plane 
     * stresses a = [S1_xz, S1_yz, ..., Sn_xz, Sn_yz] for n layer nodes
     */
 
    // set up virtual cell geometry for an qwedge
    const int numNodes = this->interpolationForExport.giveNumberOfNodes();  ///TODO lägg till interpolationForExport för SR
    std::vector<FloatArray> nodes;
    giveFictiousNodeCoordsForExport(nodes, layer);

    FEInterpolation *interpol = static_cast< FEInterpolation * >( &this->interpolationForExport );
    FloatMatrix dNdx;
    interpol->evaldNdx( dNdx, lcoords, FEIVertexListGeometryWrapper( nodes ) ); 
    
    /*    
    * 1 [d/dx  0   d/dy
    * 1   0   d/dy d/dx]
    */
    if (!intSzz) {
        int ndofs = numNodes*3;
        answer.resize(2, ndofs);
        for ( int i = 1, j = 0; i <= numNodes; i++, j += 3 ) {
            answer.at(1, j + 1) = dNdx.at(i, 1);
            answer.at(1, j + 3) = dNdx.at(i, 2);
            answer.at(2, j + 2) = dNdx.at(i, 2);
            answer.at(2, j + 3) = dNdx.at(i, 1);
        }
    } else {
        int ndofs = numNodes*2;
        answer.resize(1, ndofs);
        for ( int i = 1, j = 0; i <= numNodes; i++, j += 2 ) {
            answer.at(1, j + 1) = dNdx.at(i, 1);
            answer.at(1, j + 2) = dNdx.at(i, 2);
        }
    }
    
}
#endif




std::vector<FloatArray> 
Shell7Base :: giveFictiousNodeCoordsForExport(int layer)
{
    // compute fictious node coords
    FloatMatrix localNodeCoords;
    this->interpolationForExport.giveLocalNodeCoords(localNodeCoords);
    
    std::vector<FloatArray> nodes(localNodeCoords.giveNumberOfColumns());
    for ( int i = 1; i <= localNodeCoords.giveNumberOfColumns(); i++ ){
        FloatArray localCoords(3);
        localCoords.beColumnOf(localNodeCoords,i);

        nodes[i-1] = this->vtkEvalInitialGlobalCoordinateAt(localCoords, layer);
    }
    return nodes;
}


std::vector<FloatArray>
Shell7Base :: giveFictiousCZNodeCoordsForExport(int interface)
{
    // compute fictious node coords
    FloatMatrix localNodeCoords;
    this->interpolationForCZExport.giveLocalNodeCoords(localNodeCoords);
    
    std::vector<FloatArray> nodes(localNodeCoords.giveNumberOfColumns());
    for ( int i = 1; i <= localNodeCoords.giveNumberOfColumns(); i++ ){
        FloatArray localCoords(3);
        localCoords.beColumnOf(localNodeCoords,i);

        localCoords.at(3) = 1.0;
        nodes[i-1] = this->vtkEvalInitialGlobalCoordinateAt(localCoords, interface);
    }
    return nodes;
}

std::vector<FloatArray>
Shell7Base :: giveFictiousUpdatedNodeCoordsForExport(int layer, TimeStep *tStep)
{
    // compute fictious node coords

    FloatMatrix localNodeCoords;
    this->interpolationForExport.giveLocalNodeCoords(localNodeCoords);
    std::vector<FloatArray> nodes(localNodeCoords.giveNumberOfColumns());
    for ( int i = 1; i <= localNodeCoords.giveNumberOfColumns(); i++ ){
        FloatArray localCoords(3);
        localCoords.beColumnOf(localNodeCoords,i);

        nodes[i-1] = this->vtkEvalUpdatedGlobalCoordinateAt(localCoords, layer, tStep);
    }
    return nodes;
}


#endif



// Misc functions

FloatArrayF<9>
Shell7Base :: convV6ToV9Stress(const FloatArrayF<6> &V6)
{
    FloatArrayF<9> answer;
    answer.at(1) = V6.at(1);
    answer.at(2) = V6.at(2);
    answer.at(3) = V6.at(3);
    answer.at(4) = answer.at(7) = V6.at(4);
    answer.at(5) = answer.at(8) = V6.at(5);
    answer.at(6) = answer.at(9) = V6.at(6);
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

        nCov = this->evalInitialCovarNormalAt(ip->giveNaturalCoordinates());
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
