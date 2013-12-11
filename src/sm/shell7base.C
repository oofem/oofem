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

#include "shell7base.h"
#include "node.h"
#include "load.h"
#include "structuralms.h"
#include "mathfem.h"
#include "domain.h"
#include "equationid.h"
#include "gaussintegrationrule.h"
#include "gausspoint.h"
#include "feinterpol3d.h"
#include "boundaryload.h"
#include "constantpressureload.h"
#include "constantsurfaceload.h"
#include "vtkxmlexportmodule.h"
#include "fracturemanager.h"

namespace oofem {
FEI3dWedgeQuad Shell7Base :: interpolationForExport;

Shell7Base :: Shell7Base(int n, Domain *aDomain) : NLStructuralElement(n, aDomain),  LayeredCrossSectionInterface(),
    VTKXMLExportModuleElementInterface(), ZZNodalRecoveryModelInterface(), FailureModuleElementInterface() {}

IRResultType Shell7Base :: initializeFrom(InputRecord *ir)
{
    this->NLStructuralElement :: initializeFrom(ir);
    return IRRT_OK;
}

int
Shell7Base :: checkConsistency()
{
    NLStructuralElement :: checkConsistency();

    if ( this->layeredCS == NULL ) {
        OOFEM_ERROR("Elements derived from Shell7Base only supports layered cross section");
    }
    return ( this->layeredCS != NULL  &&  this->fei != NULL );
}


void
Shell7Base :: postInitialize()
{
    this->layeredCS = dynamic_cast< LayeredCrossSection * >( this->giveCrossSection() );
    this->fei       = dynamic_cast< FEInterpolation3d * >( this->giveInterpolation() );
    this->setupInitialNodeDirectors();

    Element :: postInitialize();
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
    //@todo move VTK method into this
    return 1;
}


double
Shell7Base :: giveGlobalZcoord(double xi, FloatArray &lc)
{
    return xi * this->layeredCS->give(CS_Thickness, & lc, NULL, this) * 0.5;
}

double // @todo move to layered crosssection
Shell7Base :: giveGlobalZcoordInLayer(double xi, int layer)
{
    // Mid position + part of the layer thickness
    return this->layeredCS->giveLayerMidZ(layer) + xi *this->layeredCS->giveLayerThickness(layer) * 0.5;
}



// Base vectors

#if 1

void
Shell7Base :: evalInitialCovarBaseVectorsAt(FloatArray &lcoords, FloatMatrix &Gcov)
{
    double zeta = giveGlobalZcoord(lcoords.at(3), lcoords);
    FloatArray M;
    FloatMatrix dNdxi;

    // In plane base vectors
    this->fei->evaldNdxi( dNdxi, lcoords, FEIElementGeometryWrapper(this) );

    FloatArray G1(3), G2(3);
    G1.zero();
    G2.zero();
    for ( int i = 1; i <= this->giveNumberOfDofManagers(); i++ ) {
        FloatArray &xbar = * this->giveNode(i)->giveCoordinates();
        M = this->giveInitialNodeDirector(i);
        FloatArray nodeCoords = ( xbar + zeta * M );
        G1 += dNdxi.at(i, 1) * nodeCoords;
        G2 += dNdxi.at(i, 2) * nodeCoords;
    }

    // Out of plane base vector = director
    FloatArray G3;
    this->evalInitialDirectorAt(lcoords, G3);     // G3=M

    Gcov.resize(3, 3);
    Gcov.setColumn(G1, 1);
    Gcov.setColumn(G2, 2);
    Gcov.setColumn(G3, 3);
}

void
Shell7Base :: edgeEvalInitialCovarBaseVectorsAt(FloatArray &lcoords, const int iedge, FloatArray &G1, FloatArray &G3)
{
    double zeta = 0.0;     // no variation i z (yet)
    //double zeta = lcoords.at(3);
    FloatArray M, dNdxi;

    IntArray edgeNodes;
    this->fei->computeLocalEdgeMapping(edgeNodes, iedge);
    this->fei->edgeEvaldNdxi( dNdxi, iedge, lcoords, FEIElementGeometryWrapper(this) );

    // Base vector along edge
    G1.resize(3);
    G1.zero();
    for ( int i = 1; i <= edgeNodes.giveSize(); i++ ) {
        FloatArray &xbar = * this->giveNode( edgeNodes.at(i) )->giveCoordinates();
        M = this->giveInitialNodeDirector( edgeNodes.at(i) );
        FloatArray nodeCoords = ( xbar + zeta * M );
        G1 += dNdxi.at(i) * nodeCoords;
    }

    // Director will be the second base vector
    this->edgeEvalInitialDirectorAt(lcoords, G3, iedge);
}


void
Shell7Base :: evalInitialContravarBaseVectorsAt(FloatArray &lCoords, FloatMatrix &Gcon)
{
    FloatMatrix Gcov;
    this->evalInitialCovarBaseVectorsAt(lCoords, Gcov);
    this->giveDualBase(Gcov, Gcon);
}


void
Shell7Base :: evalContravarBaseVectorsAt(FloatArray &lcoords, FloatMatrix &gcon, FloatArray &genEps)
{
    // Not in use at the moment!
    FloatMatrix gcov;
    this->evalCovarBaseVectorsAt(lcoords, gcov, genEps);
    this->giveDualBase(gcov, gcon);
}

void
Shell7Base :: giveDualBase(FloatMatrix &base1, FloatMatrix &base2)
{
    // Computes the dual base through thte inversion of the metric tensor
    FloatMatrix gMetric, ginv;
    gMetric.beTProductOf(base1, base1); // Metric tensor
    ginv.beInverseOf(gMetric);
    base2.beProductTOf(base1, ginv);
}

void
Shell7Base :: evalInitialDirectorAt(GaussPoint *gp, FloatArray &answer)
{
    // Interpolates between the node directors
    FloatArray &lcoords = * gp->giveCoordinates();
    FloatArray N;
    this->fei->evalN( N, lcoords, FEIElementGeometryWrapper(this) );

    answer.resize(3);
    answer.zero();
    for ( int i = 1; i <= this->giveNumberOfDofManagers(); i++ ) {
        answer.add( N.at(i), this->giveInitialNodeDirector(i) );
    }
}

void
Shell7Base :: evalInitialDirectorAt(FloatArray &lcoords, FloatArray &answer)
{
    // Interpolates between the node directors
    FloatArray N;
    this->fei->evalN( N, lcoords, FEIElementGeometryWrapper(this) );
    answer.resize(3);
    answer.zero();
    for ( int i = 1; i <= this->giveNumberOfDofManagers(); i++ ) {
        answer.add( N.at(i), this->giveInitialNodeDirector(i) );
    }
}


void
Shell7Base :: edgeEvalInitialDirectorAt(FloatArray &lcoords, FloatArray &answer, const int iEdge)
{
    // Interpolates between the node directors along an edge

    FloatArray N;
    IntArray edgeNodes;
    this->fei->computeLocalEdgeMapping(edgeNodes, iEdge);
    this->fei->edgeEvalN( N, iEdge, lcoords, FEIElementGeometryWrapper(this) );

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
    //@todo average field for smooth geometry
    FloatArray M(3), G1(3), G2(3), lcoords(2), nodeLocalXiCoords, nodeLocalEtaCoords;

    // Give the local coordinates for the element nodes - all at once
    this->giveLocalNodeCoords(nodeLocalXiCoords, nodeLocalEtaCoords);
    int nDofMan = this->giveNumberOfDofManagers();
    this->initialNodeDirectors.resize(nDofMan);
    FloatMatrix dNdxi;
    for ( int node = 1; node <= nDofMan; node++ ) {
        this->initialNodeDirectors [ node - 1 ].resize(3);
        this->initialNodeDirectors [ node - 1 ].zero();
        lcoords.at(1) = nodeLocalXiCoords.at(node);
        lcoords.at(2) = nodeLocalEtaCoords.at(node);

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
        this->initialNodeDirectors [ node - 1 ].add(M);
    }
}



void
Shell7Base :: evalCovarBaseVectorsAt(FloatArray &lcoords, FloatMatrix &gcov, FloatArray &genEps)
{
    // Evaluates the covariant base vectors in the current configuration
    FloatArray g1;
    FloatArray g2;
    FloatArray g3;
    double zeta = giveGlobalZcoord(lcoords.at(3), lcoords);

    FloatArray dxdxi1, dxdxi2, m, dmdxi1, dmdxi2;
    double dgamdxi1, dgamdxi2, gam;
    this->giveGeneralizedStrainComponents(genEps, dxdxi1, dxdxi2, dmdxi1, dmdxi2, m, dgamdxi1, dgamdxi2, gam);

    double fac1 = ( zeta + 0.5 * gam * zeta * zeta );
    double fac2 = ( 0.5 * zeta * zeta );
    double fac3 = ( 1.0 + zeta * gam );

    g1 = dxdxi1 + fac1 * dmdxi1 + fac2 * dgamdxi1 * m;
    g2 = dxdxi2 + fac1 * dmdxi2 + fac2 * dgamdxi2 * m;
    g3 = fac3 * m;

    gcov.resize(3, 3);
    gcov.setColumn(g1, 1);
    gcov.setColumn(g2, 2);
    gcov.setColumn(g3, 3);
}


void
Shell7Base :: edgeEvalCovarBaseVectorsAt(FloatArray &lcoords, const int iedge, FloatMatrix &gcov, TimeStep *tStep)
{
    // Evaluates the covariant base vectors in the current configuration for en edge

    //double zeta = 0.0; //@todo fix integration rule for arbitrary z-coord
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
    dxdxi.setValues( 3, genEpsEdge.at(1), genEpsEdge.at(2), genEpsEdge.at(3) );
    dmdxi.setValues( 3, genEpsEdge.at(4), genEpsEdge.at(5), genEpsEdge.at(6) );
    m.setValues( 3, genEpsEdge.at(7), genEpsEdge.at(8), genEpsEdge.at(9) );
    double dgamdxi = genEpsEdge.at(10);
    double gam     = genEpsEdge.at(11);

    double fac1 = ( zeta + 0.5 * gam * zeta * zeta );
    double fac2 = ( 0.5 * zeta * zeta );
    double fac3 = ( 1.0 + zeta * gam );

    FloatArray g1, g2, g3;
    g2 = dxdxi + fac1 * dmdxi + fac2 * dgamdxi * m; // base vector along the edge
    g3 = fac3 * m;                              // director field

    g2.normalize();
    g3.normalize();
    g1.beVectorProductOf(g2, g3);
    g1.normalize();
    gcov.resize(3, 3);
    gcov.setColumn(g1, 1);
    gcov.setColumn(g2, 2);
    gcov.setColumn(g3, 3);
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

    FloatArray solVec, totalSolVec;
    this->giveUpdatedSolutionVector(totalSolVec, tStep);
    this->giveUpdatedSolutionVector(solVec, tStep); // a

    // first 'solVec' corresponds to the point where the tangent i evaluated and solVecLeft, solVecRight
    // corresponds to the solution used for evaluation of the lambda matrices
    //@todo rewrite this method since the XFEM part does not use this anymore
    // new_computeBulkTangentMatrix(answer, solVecPoint, solVecLeft, solVecRight, rMode, tStep);
    this->new_computeBulkTangentMatrix(answer, solVec, solVec, solVec, rMode, tStep);

    // Add contribution due to pressure load
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

    FloatMatrix eye(3, 3), aEye(3, 3);
    eye.beUnitMatrix();
    aEye = eye;
    aEye.times(a);
    lambda [ 0 ].resize(3, 18);
    lambda [ 0 ].zero();
    lambda [ 1 ].resize(3, 18);
    lambda [ 1 ].zero();

    lambda [ 0 ].setSubMatrix(eye, 1, 1);
    lambda [ 0 ].setSubMatrix(aEye, 1, 7);
    lambda [ 1 ].setSubMatrix(eye, 1, 4);
    lambda [ 1 ].setSubMatrix(aEye, 1, 10);

    FloatMatrix bdg1eye(3, 3), bdg2eye(3, 3);
    bdg1eye = eye;
    bdg1eye.times(b * dgam1);
    bdg2eye = eye;
    bdg2eye.times(b * dgam2);
    lambda [ 0 ].setSubMatrix(bdg1eye, 1, 13);
    lambda [ 1 ].setSubMatrix(bdg2eye, 1, 13);

    FloatArray bm(3), bdm1(3), bdm2(3);
    bm   = b * m;
    bdm1 = b * dm1;
    bdm2 = b * dm2;
    lambda [ 0 ].setColumn(bm, 16);
    lambda [ 0 ].setColumn(bdm1, 18);
    lambda [ 1 ].setColumn(bm, 17);
    lambda [ 1 ].setColumn(bdm2, 18);

    // lambda3 =  ( 0,   0,   0 ,   0 ,     c*I   ,   0 ,   0 ,   xi*m )
    lambda [ 2 ].resize(3, 18);
    lambda [ 2 ].zero();
    lambda [ 2 ].at(1, 13) = lambda [ 2 ].at(2, 14) = lambda [ 2 ].at(3, 15) = c;
    FloatArray zm(3);
    zm = zeta * m;
    lambda [ 2 ].setColumn(zm, 18);
}

void
Shell7Base :: computeLambdaNMatrix(FloatMatrix &lambda, FloatArray &genEps, double zeta)
{
    // computes the lambda^n matrix associated with the variation and linearization of the position vector x.
    FloatArray m(3);
    m.setValues( 3, genEps.at(13), genEps.at(14), genEps.at(15) );
    double gam = genEps.at(18);

    // thickness coefficients
    double a = zeta + 0.5 * gam * zeta * zeta;
    double b = 0.5 * zeta * zeta;

    // lambda =  ( I, a*I, b*m )
    lambda.resize(3, 7);
    lambda.zero();
    lambda.at(1, 1) = lambda.at(2, 2) = lambda.at(3, 3) = 1.0;
    lambda.at(1, 4) = lambda.at(2, 5) = lambda.at(3, 6) = a;
    lambda.setColumn(b * m, 7);
}


void
Shell7Base :: new_computeBulkTangentMatrix(FloatMatrix &answer, FloatArray &solVec, FloatArray &solVecI, FloatArray &solVecJ, MatResponseMode rMode, TimeStep *tStep)
{
    //@todo optimize method - remove I,J-stuff since XFEM does not need this anymore
    FloatMatrix A [ 3 ] [ 3 ], lambdaI [ 3 ], lambdaJ [ 3 ];
    FloatMatrix L(18, 18), B;
    FloatMatrix K, tempAnswer;

    int ndofs = Shell7Base :: giveNumberOfDofs();
    answer.resize(ndofs, ndofs);
    tempAnswer.resize(ndofs, ndofs);
    answer.zero();
    tempAnswer.zero();

    int numberOfLayers = this->layeredCS->giveNumberOfLayers();
    FloatMatrix temp;
    FloatArray genEpsI, genEpsJ, genEps, lCoords;

    for ( int layer = 1; layer <= numberOfLayers; layer++ ) {
        IntegrationRule *iRule = integrationRulesArray [ layer - 1 ];
        StructuralMaterial *mat = static_cast< StructuralMaterial * >( domain->giveMaterial( this->layeredCS->giveLayerMaterial(layer) ) );

        for ( int i = 0; i < iRule->giveNumberOfIntegrationPoints(); i++ ) {
            GaussPoint *gp = iRule->getIntegrationPoint(i);
            lCoords = * gp->giveCoordinates();

            this->computeBmatrixAt(lCoords, B, 0, 0);
            this->computeGeneralizedStrainVectorNew(genEpsI, solVecI, B);
            this->computeGeneralizedStrainVectorNew(genEpsJ, solVecJ, B);
            this->computeGeneralizedStrainVectorNew(genEps, solVec, B);

            // Material stiffness
            Shell7Base :: computeLinearizedStiffness(gp, mat, tStep, A, genEps);

            double zeta = giveGlobalZcoord(gp->giveCoordinate(3), lCoords);
            this->computeLambdaGMatrices(lambdaI, genEpsI, zeta);
            this->computeLambdaGMatrices(lambdaJ, genEpsJ, zeta);

            // L = sum_{i,j} (lambdaI_i)^T * A^ij * lambdaJ_j
            // @todo Naive implementation - should be optimized
            // note: L will only be symmetric if lambdaI = lambdaJ
            L.zero();
            for ( int j = 0; j < 3; j++ ) {
                for ( int k = 0; k < 3; k++ ) {
                    this->computeTripleProduct(temp, lambdaI [ j ], A [ j ] [ k ], lambdaJ [ k ]);
                    L.add(temp);
                }
            }

            this->computeTripleProduct(K, B, L, B);
            double dV = this->computeVolumeAroundLayer(gp, layer);
            tempAnswer.add(dV, K);
        }
    }

    const IntArray &ordering = this->giveOrdering(All);
    answer.assemble(tempAnswer, ordering, ordering);
}

void
Shell7Base :: computeLinearizedStiffness(GaussPoint *gp, StructuralMaterial *mat, TimeStep *tStep,
                                         FloatMatrix A [ 3 ] [ 3 ], FloatArray &genEps)
{
    FloatArray cartStressVector, contravarStressVector;
    FloatMatrix D, Dcart, S;

    //A = L^iklj * (g_k x g_l) + S^ij*I
    //mat->giveStiffnessMatrix(Dcart, TangentStiffness, gp, tStep);     // L_ijkl - cartesian system (Voigt)
    mat->give3dMaterialStiffnessMatrix(Dcart, TangentStiffness, gp, tStep);     // L_ijkl - cartesian system (Voigt)
    this->transInitialCartesianToInitialContravar(gp, Dcart, D);      // L^ijkl - curvilinear system (Voigt)

    FloatArray lcoords = * gp->giveCoordinates();

    //this->computeStressMatrix(cartStressVector, genEps, gp, mat, tStep);
    this->computeStressVectorInMaterial(cartStressVector, genEps, gp, mat, tStep);
    this->transInitialCartesianToInitialContravar(gp, cartStressVector, contravarStressVector);
    S.beMatrixFormOfStress(contravarStressVector);

    FloatMatrix gcov;
    this->evalCovarBaseVectorsAt(lcoords, gcov, genEps);

    FloatMatrix gg11, gg12, gg13, gg21, gg22, gg23, gg31, gg32, gg33;
    FloatArray g1, g2, g3;
    g1.beColumnOf(gcov, 1);
    g2.beColumnOf(gcov, 2);
    g3.beColumnOf(gcov, 3);

    gg11.beDyadicProductOf(g1, g1);
    gg12.beDyadicProductOf(g1, g2);
    gg13.beDyadicProductOf(g1, g3);
    gg22.beDyadicProductOf(g2, g2);
    gg23.beDyadicProductOf(g2, g3);
    gg33.beDyadicProductOf(g3, g3);

    gg21.beTranspositionOf(gg12);
    gg31.beTranspositionOf(gg13);
    gg32.beTranspositionOf(gg23);

    // position 11
    A [ 0 ] [ 0 ].resize(3, 3);
    A [ 0 ] [ 0 ].beUnitMatrix();
    A [ 0 ] [ 0 ].times( S.at(1, 1) );
    A [ 0 ] [ 0 ].add(D.at(1, 1), gg11);
    A [ 0 ] [ 0 ].add(D.at(1, 6), gg12);
    A [ 0 ] [ 0 ].add(D.at(1, 5), gg13); // = (5,1)^T ?
    A [ 0 ] [ 0 ].add(D.at(6, 1), gg21); // = (1,6)^T ?
    A [ 0 ] [ 0 ].add(D.at(6, 6), gg22);
    A [ 0 ] [ 0 ].add(D.at(6, 5), gg23);
    A [ 0 ] [ 0 ].add(D.at(5, 1), gg31);
    A [ 0 ] [ 0 ].add(D.at(5, 6), gg32); // = (6,5)^T ?
    A [ 0 ] [ 0 ].add(D.at(5, 5), gg33);


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


    // position 22
    A [ 1 ] [ 1 ].resize(3, 3);
    A [ 1 ] [ 1 ].beUnitMatrix();
    A [ 1 ] [ 1 ].times( S.at(2, 2) );
    A [ 1 ] [ 1 ].add(D.at(6, 6), gg11);
    A [ 1 ] [ 1 ].add(D.at(6, 2), gg12);
    A [ 1 ] [ 1 ].add(D.at(6, 4), gg13);
    A [ 1 ] [ 1 ].add(D.at(2, 6), gg21); // = (6,2)^T ?
    A [ 1 ] [ 1 ].add(D.at(2, 2), gg22);
    A [ 1 ] [ 1 ].add(D.at(2, 4), gg23);
    A [ 1 ] [ 1 ].add(D.at(4, 6), gg31); // = (6,4)^T ?
    A [ 1 ] [ 1 ].add(D.at(4, 2), gg32); // = (2,4)^T ?
    A [ 1 ] [ 1 ].add(D.at(4, 4), gg33);


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
    A [ 2 ] [ 2 ].add(D.at(4, 5), gg21); // = (5,4)^T ?
    A [ 2 ] [ 2 ].add(D.at(4, 4), gg22);
    A [ 2 ] [ 2 ].add(D.at(4, 3), gg23);
    A [ 2 ] [ 2 ].add(D.at(3, 5), gg31); // = (5,3)^T ?
    A [ 2 ] [ 2 ].add(D.at(3, 4), gg32); // = (4,3)^T ?
    A [ 2 ] [ 2 ].add(D.at(3, 3), gg33);

    // position 21
    A [ 1 ] [ 0 ].beTranspositionOf(A [ 0 ] [ 1 ]);

    // position 31
    A [ 2 ] [ 0 ].beTranspositionOf(A [ 0 ] [ 2 ]);

    // position 32
    A [ 2 ] [ 1 ].beTranspositionOf(A [ 1 ] [ 2 ]);
}

void
Shell7Base :: computePressureTangentMatrix(FloatMatrix &answer, Load *load, const int iSurf, TimeStep *tStep)
{
    // Computes tangent matrix associated with the linearization of pressure loading. Assumes constant pressure.
    IntegrationRule *iRule = specialIntegrationRulesArray [ 1 ];   // rule #2 for surface integration
    ConstantPressureLoad *pLoad = dynamic_cast< ConstantPressureLoad * >( load );

    FloatMatrix N, B, NLB, L(7, 18), gcov, W1, W2;
    FloatArray lcoords(3), solVec, pressure;
    FloatArray g1, g2, genEps;
    FloatMatrix lambdaG [ 3 ], lambdaN;

    double xi   = pLoad->giveLoadOffset();
    this->giveUpdatedSolutionVector(solVec, tStep);

    int ndof = Shell7Base :: giveNumberOfDofs();
    answer.resize(ndof, ndof);
    answer.zero();
    for ( int i = 0; i < iRule->giveNumberOfIntegrationPoints(); i++ ) {
        IntegrationPoint *ip = iRule->getIntegrationPoint(i);
        lcoords.at(1) = ip->giveCoordinate(1);
        lcoords.at(2) = ip->giveCoordinate(2);
        lcoords.at(3) = xi;     // local coord where load is applied

        double zeta = this->giveGlobalZcoord( xi, * ip->giveCoordinates() );
        this->computeNmatrixAt(lcoords, N);
        this->computeBmatrixAt(lcoords, B);
        genEps.beProductOf(B, solVec);

        // Traction tangent, L =  lambdaN * ( W2*lambdaG_1 - W1*lambdaG_2  )
        load->computeValueAt(pressure, tStep, * ( ip->giveCoordinates() ), VM_Total);        // pressure component
        this->evalCovarBaseVectorsAt(lcoords, gcov, genEps);
        g1.beColumnOf(gcov, 1);
        g2.beColumnOf(gcov, 2);
        W1 = this->giveAxialMatrix(g1);
        W2 = this->giveAxialMatrix(g2);

        this->computeLambdaGMatrices(lambdaG, genEps, zeta);
        this->computeLambdaNMatrix(lambdaN, genEps, zeta);
        FloatMatrix W2L, W1L;
        W2L.beProductOf(W2, lambdaG [ 0 ]);
        W1L.beProductOf(W1, lambdaG [ 1 ]);
        W2L.subtract(W1L);
        L.beTProductOf(lambdaN, W2L);
        L.times( -pressure.at(1) );


        // Tangent matrix (K = N^T*L*B*dA)
        this->computeTripleProduct(NLB, N, L, B);
        double dA = this->computeAreaAround(ip, xi);
        answer.add(dA, NLB);
    }
}


FloatMatrix
Shell7Base :: giveAxialMatrix(const FloatArray &v)
{
    // creates the skew-symmetric matrix W defined such that
    // crossProduct(u,v) = W(u)*v
    // W = [   0   v(3)  -v(2)
    //      -v(3)   0     v(1)
    //       v(2) -v(1)    0  ]
    //
    FloatMatrix answer(3, 3);
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
Shell7Base :: computeFAt(FloatArray &lCoords, FloatMatrix &answer, FloatArray &genEps)
{
    // Compute deformation gradient as open product(g_i, G^i) = gcov*Gcon^T
    // Output is a 3x3 matrix
    FloatMatrix gcov, Gcon;
    //FloatArray lcoords;
    //lcoords = *gp->giveCoordinates();
    this->evalCovarBaseVectorsAt(lCoords, gcov, genEps);
    this->evalInitialContravarBaseVectorsAt(lCoords, Gcon);
    answer.beProductTOf(gcov, Gcon);
}

void
Shell7Base :: computeE(FloatMatrix &answer, FloatMatrix &F)
{
    // Computes the Green-Lagrange strain tensor: E=0.5(C-I)
    // Output is a 3x3 matrix
    // utilize the symmetry to save time?
    answer.beTProductOf(F, F);    // C-Right Caucy-Green deformation tensor F^T*F
    answer.at(1, 1) += -1;
    answer.at(2, 2) += -1;
    answer.at(3, 3) += -1;
    answer.times(0.5);
}

void
Shell7Base :: computeStrainVectorF(FloatArray &answer, GaussPoint *gp, TimeStep *tStep, FloatArray &genEps)
//Shell7Base :: computeStrainVectorFrom(FloatArray &answer, GaussPoint *gp, TimeStep *tStep, FloatArray &genEps)
{
    // Computes the Green-Lagrange strain tensor: E=0.5(C-I)
    FloatMatrix F, E;
    FloatArray lcoords;
    lcoords = * gp->giveCoordinates();
    this->computeFAt(lcoords, F, genEps);       // Deformation gradient
    this->computeE(E, F);                  // Green-Lagrange strain tensor
    answer.beSymVectorFormOfStrain(E); // Convert to reduced Voight form (6 components)
}

void
Shell7Base :: computeStressMatrix(FloatArray &answer, FloatArray &genEps, GaussPoint *gp, Material *mat, TimeStep *tStep)
//computeStrainVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep)
{
    FloatArray vE;
    this->computeStrainVectorF(vE, gp, tStep, genEps);     // Green-Lagrange strain vector in Voigt form
    static_cast< StructuralMaterial * >( mat )->giveRealStressVector(answer, gp, vE, tStep);
}


void
Shell7Base :: computeCauchyStressVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep)
{
    // Compute Cauchy stress from 2nd Piola stress
    FloatArray solVec, lCoords;
    this->giveUpdatedSolutionVector(solVec, tStep);
    FloatMatrix B;
    lCoords = * gp->giveCoordinates();
    this->computeBmatrixAt(lCoords, B);
    FloatArray genEps;
    this->computeGeneralizedStrainVectorNew(genEps, solVec, B);
    FloatMatrix F;
    this->computeFAt(lCoords, F, genEps);

    FloatArray vS;
    giveIPValue(vS, gp, IST_StressTensor, tStep); // expects Second PK stress

    FloatMatrix S, temp, sigma;
    S.beMatrixFormOfStress(vS);
    temp.beProductTOf(S, F);
    sigma.beProductOf(F, temp);
    sigma.times( 1.0 / F.giveDeterminant() );
    answer.beSymVectorForm(sigma);
}


void
Shell7Base :: computeStressResultantsAt(GaussPoint *gp, FloatArray &Svec, FloatArray &S1g, FloatArray &S2g, FloatArray &S3g, FloatArray &genEps)
{
    // Computes the stress resultants in the covariant system Sig =S(i,j)*g_j
    FloatMatrix gcov;
    FloatArray lCoords = * gp->giveCoordinates();
    this->evalCovarBaseVectorsAt(lCoords, gcov, genEps);
    FloatMatrix S, Sig;
    S.beMatrixFormOfStress(Svec);
    Sig.beProductTOf(gcov, S);       // Stress resultants stored in each column
    S1g.beColumnOf(Sig, 1);
    S2g.beColumnOf(Sig, 2);
    S3g.beColumnOf(Sig, 3);
}

int
Shell7Base :: giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep)
{
    // Compute special IST quantities this element should support
    switch ( type ) {
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
    answer.resize( this->giveNumberOfDofs() );
    answer.zero();
    this->computeSectionalForces(answer, tStep, solVec, useUpdatedGpRecord);

    ///@todo How to treat the convective force? Only active during dynamic simulations
}


void
Shell7Base :: computeSectionalForces(FloatArray &answer, TimeStep *tStep, FloatArray &solVec, int useUpdatedGpRecord)
{
    //
    int ndofs = Shell7Base :: giveNumberOfDofs();

    int numberOfLayers = this->layeredCS->giveNumberOfLayers();     // conversion of types
    FloatArray f(ndofs), ftemp, sectionalForces;
    FloatArray genEps, genEpsD, totalSolVec, lCoords;
    FloatMatrix B;

    f.zero();
    for ( int layer = 1; layer <= numberOfLayers; layer++ ) {
        IntegrationRule *iRuleL = integrationRulesArray [ layer - 1 ];
        Material *mat = domain->giveMaterial( this->layeredCS->giveLayerMaterial(layer) );

        for ( int j = 1; j <= iRuleL->giveNumberOfIntegrationPoints(); j++ ) {
            GaussPoint *gp = iRuleL->getIntegrationPoint(j - 1);
            lCoords = * gp->giveCoordinates();
            this->computeBmatrixAt(lCoords, B);

            this->computeGeneralizedStrainVectorNew(genEpsD, solVec, B);
            this->giveUpdatedSolutionVector(totalSolVec, tStep);
            this->computeGeneralizedStrainVectorNew(genEps, totalSolVec, B);

            double zeta = giveGlobalZcoord(gp->giveCoordinate(3), lCoords);
            this->computeSectionalForcesAt(sectionalForces, gp, mat, tStep, genEps, genEpsD, zeta);

            // Computation of sectional forces: f = B^t*[N M T Ms Ts]^t
            ftemp.beTProductOf(B, sectionalForces);
            double dV = this->computeVolumeAroundLayer(gp, layer);
            f.add(dV, ftemp);
        }
    }

    answer.resize(ndofs);
    answer.zero();
    const IntArray &ordering_all = this->giveOrdering(All);
    answer.assemble(f, ordering_all);
}


void
Shell7Base :: computeSectionalForcesAt(FloatArray &sectionalForces, IntegrationPoint *ip, Material *mat, TimeStep *tStep, FloatArray &genEpsC, FloatArray &genEpsD, double zeta)
{
    FloatArray S1g(3), S2g(3), S3g(3);
    FloatArray cartStressVector, contravarStressVector;
    FloatMatrix lambda [ 3 ];

    this->computeStressVectorInMaterial(cartStressVector, genEpsC, ip, mat, tStep);
    //this->computeStressMatrix(cartStressVector, genEpsC, ip, mat, tStep);
    this->transInitialCartesianToInitialContravar(ip, cartStressVector, contravarStressVector);
    this->computeStressResultantsAt(ip, contravarStressVector, S1g, S2g, S3g, genEpsC);
    this->computeLambdaGMatrices(lambda, genEpsD, zeta);

    // f = lambda_1^T * S1g + lambda_2^T * S2g + lambda_3^T * S3g
    sectionalForces.resize(18);
    FloatArray temp;
    sectionalForces.zero();
    temp.beTProductOf(lambda [ 0 ], S1g);
    sectionalForces.add(temp);
    temp.beTProductOf(lambda [ 1 ], S2g);
    sectionalForces.add(temp);
    temp.beTProductOf(lambda [ 2 ], S3g);
    sectionalForces.add(temp);
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
    this->fei->evaldNdxi( dNdxi, lcoords, FEIElementGeometryWrapper(this) );

    FloatArray M, dM1(3), dM2(3), dX1(3), dX2(3);
    double gam, dg1, dg2;
    FloatMatrix B;
    this->computeBmatrixAt(gp, B);

    FloatArray initSolVec, genEps;
    this->giveInitialSolutionVector(initSolVec);

    this->computeGeneralizedStrainVectorNew(genEps, initSolVec, B);
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
    IntegrationRule *iRule = specialIntegrationRulesArray [ 1 ];

    //------------------------------
    FloatMatrix N, Ntm, NtmN, temp;
    FloatArray solVec, lCoords;
    this->giveUpdatedSolutionVector(solVec, tStep);

    ///@todo: Should be changed to integration over the layers and use the corresponding material
    Material *mat = domain->giveMaterial( this->layeredCS->giveLayerMaterial(1) );     // for now, while I don't have an analytical exp.

    for ( int i = 0; i < iRule->giveNumberOfIntegrationPoints(); i++ ) {
        GaussPoint *gp = iRule->getIntegrationPoint(i);
        lCoords = * gp->giveCoordinates();
        this->computeNmatrixAt(lCoords, N);
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
    const IntArray &ordering_all = this->giveOrdering(All);
    answer.assemble(temp, ordering_all);
    answer.symmetrized();

 #if 0
    FloatMatrix A(18, 18), B(18, 18), C(18, 6), D(18, 18), E(18, 6), F(6, 6), eigVec, temp2;

    temp2.beSubMatrixOf(temp, 1, 18, 19, 36);
    A.beSubMatrixOf(temp, 1, 18, 1, 18);
    B.beSubMatrixOf(temp, 1, 18, 19, 36);
    C.beSubMatrixOf(temp, 1, 18, 37, 42);
    D.beSubMatrixOf(temp, 19, 36, 19, 36);
    E.beSubMatrixOf(temp, 19, 36, 37, 42);
    F.beSubMatrixOf(temp, 37, 42, 37, 42);

    FloatArray sums(6);
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
 #endif
}


void
Shell7Base :: giveMassFactorsAt(GaussPoint *gp, FloatArray &factors, double &gam)
{
    FloatArray coeff;
    this->computeThicknessMappingCoeff(gp, coeff);
    double a1 = coeff.at(1);
    double a2 = coeff.at(2);
    double a3 = coeff.at(3);

    double h  = this->giveCrossSection()->give(CS_Thickness, gp->giveCoordinates(), NULL, this);
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
    // For analytically integrated throught he thickness, see computeMassMatrix

    FloatMatrix mass, temp;
    FloatArray solVec, lCoords;
    this->giveUpdatedSolutionVector(solVec, tStep);
    int ndofs = this->giveNumberOfDofs();
    temp.resize(ndofs, ndofs);
    temp.zero();

    int numberOfLayers = this->layeredCS->giveNumberOfLayers();     // conversion of data

    FloatMatrix M11(18, 18), M12(18, 18), M13(18, 6), M22(18, 18), M23(18, 6), M33(6, 6), M(42, 42);
    M.zero();
    M11.zero();
    M12.zero();
    M13.zero();
    M22.zero();
    M23.zero();
    M33.zero();

    for ( int layer = 1; layer <= numberOfLayers; layer++ ) {
        IntegrationRule *iRuleL = integrationRulesArray [ layer - 1 ];
        Material *mat = domain->giveMaterial( this->layeredCS->giveLayerMaterial(layer) );
 #if 0
        for ( int j = 1; j <= iRuleL->giveNumberOfIntegrationPoints(); j++ ) {
            GaussPoint *gp = iRuleL->getIntegrationPoint(j - 1);

            FloatMatrix N11, N22, N33;
            this->computeNmatricesAt(gp, N11, N22, N33);
            FloatArray xbar, m;
            double gam = 0.;
            //this->computeSolutionFields(xbar, m, gam, solVec, N11, N22, N33);
            FloatArray localCoords = * gp->giveCoordinates();
            this->giveUnknownsAt(localCoords, solVec, xbar, m, gam, tStep);
            /* Consistent mass matrix M = int{N^t*mass*N}
             *
             *         3    3    1
             *         3 [a*I  b*I   c*m      [A  B  C
             * mass =  3       d*I   e*m    =     D  E
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

        answer.resize(ndofs, ndofs);
        answer.zero();

        const IntArray &ordering_phibar = this->giveOrdering(Midplane);
        const IntArray &ordering_m = this->giveOrdering(Director);
        const IntArray &ordering_gam = this->giveOrdering(InhomStrain);
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
 #else
        for ( int j = 1; j <= iRuleL->giveNumberOfIntegrationPoints(); j++ ) {
            GaussPoint *gp = iRuleL->getIntegrationPoint(j - 1);
            lCoords = * gp->giveCoordinates();
            FloatMatrix lambda, B, N, temp;
            FloatArray genEps;
            this->computeBmatrixAt(lCoords, B);
            genEps.beProductOf(B, solVec);
            //double zeta = giveGlobalZcoord(gp);
            double zeta = giveGlobalZcoord(gp->giveCoordinate(3), lCoords);
            this->computeLambdaNMatrix(lambda, genEps, zeta);

            // could also create lambda*N and then plusProdSymm - probably faster
            mass.beTProductOf(lambda, lambda);
            this->computeNmatrixAt(lCoords, N);
            temp.beProductOf(mass, N);

            double dV = this->computeVolumeAroundLayer(gp, layer);
            double rho = mat->give('d', gp);
            M.plusProductSymmUpper(N, temp, rho * dV);
        }
        M.symmetrized();
        const IntArray &ordering_all = this->giveOrdering(All);
        answer.zero();
        answer.assemble(M, ordering_all, ordering_all);
 #endif
 #if 0

        FloatArray sums(6);
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



 #endif
    }
}


void
Shell7Base :: computeConvectiveMassForce(FloatArray &answer, TimeStep *tStep)
{
    //@todo: very old version should be checked
    // Analytically integrated over the thickness. Constant density assumed.

    IntegrationRule *iRule = specialIntegrationRulesArray [ 1 ];   // rule 2 for mid-plane integration only
    GaussPoint *gp;
    FloatMatrix N;
    FloatArray lCoords, a, da, m, dm, aVec, daVec, fm(7), fM;
    double gam, dgam, dA;

    answer.resize(42);
    answer.zero();
    for ( int i = 0; i < iRule->giveNumberOfIntegrationPoints(); i++ ) {
        gp = iRule->getIntegrationPoint(i);
        lCoords = * gp->giveCoordinates();
        this->computeNmatrixAt(lCoords, N);
        this->giveUpdatedSolutionVector(aVec, tStep);
        // Fix for the new numbering in B & N
        this->computeVectorOfDofIDs(EID_MomentumBalance, VM_Velocity, tStep, daVec);

        a.beProductOf(N, aVec);        // [ x,  m,  gam]^T
        da.beProductOf(N, daVec);        // [dx, dm, dgam]^T
        m.setValues( 3,  a.at(4),  a.at(5),  a.at(6) );
        gam =  a.at(7);
        dm.setValues( 3, da.at(4), da.at(5), da.at(6) );
        dgam = da.at(7);

        double a1, a2, a3, h, h2, h3, h5, fac1, fac2, fac3, rho;
        FloatArray coeff;

        Material *mat = domain->giveMaterial( this->layeredCS->giveLayerMaterial(1) ); ///@todo fix this method
        rho = mat->give('d', gp);
        h   = this->giveCrossSection()->give(CS_Thickness, & lCoords, NULL, this);
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
        fM.beTProductOf(N, fm);
        fM.times(dA);
        answer.add(fM);
    }
}


void
Shell7Base :: computeTripleProduct(FloatMatrix &answer, const FloatMatrix &a, const FloatMatrix &b, const FloatMatrix &c)
{
    // Computes the product a^T*b*c
    // @todo room for optimization + add scalar product for integration
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
        this->computeTractionForce(answer, iEdge, edgeLoad, tStep, mode);
        return;
    } else {
        _error("Shell7Base :: computeEdgeLoadVectorAt: load type not supported");
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
        _error("Shell7Base :: computeSurfaceLoadVectorAt: load type not supported");
        return;
    }
}

void
Shell7Base :: computePressureForce(FloatArray &answer, FloatArray solVec, const int iSurf, BoundaryLoad *surfLoad, TimeStep *tStep, ValueModeType mode)
{
    // Computes pressure loading. Acts normal to the current (deformed) surface.
    ConstantPressureLoad *pLoad = dynamic_cast< ConstantPressureLoad * >( surfLoad );
    IntegrationRule *iRule = specialIntegrationRulesArray [ 1 ];   // rule #2 for surface integration

    FloatMatrix N, B, lambda;
    FloatArray Fp, fp, genEps, genEpsC, lCoords, traction, solVecC;
    answer.resize( Shell7Base :: giveNumberOfDofs() );
    answer.zero();
    for ( int i = 0; i < iRule->giveNumberOfIntegrationPoints(); i++ ) {
        IntegrationPoint *ip = iRule->getIntegrationPoint(i);
        lCoords = * ip->giveCoordinates();
        double zeta = this->giveGlobalZcoord(pLoad->giveLoadOffset(), lCoords);
        this->computeBmatrixAt(lCoords, B);
        this->computeNmatrixAt(lCoords, N);
        this->giveUpdatedSolutionVector(solVecC, tStep);
        genEpsC.beProductOf(B, solVecC);
        this->computePressureForceAt(ip, traction, iSurf, genEpsC, surfLoad, tStep, mode);

        genEps.beProductOf(B, solVec);
        this->computeLambdaNMatrix(lambda, genEps, zeta);
        fp.beTProductOf(lambda, traction);

        double xi = pLoad->giveLoadOffset();
        double dA = this->computeAreaAround(ip, xi);

        Fp.beTProductOf(N, fp);
        Fp.times(dA);
        answer.assemble( Fp, this->giveOrdering(All) );
    }
}


void
Shell7Base :: computePressureForceAt(GaussPoint *gp, FloatArray &traction, const int iSurf, FloatArray genEps, BoundaryLoad *surfLoad, TimeStep *tStep, ValueModeType mode)
{
    // Computes pressure loading. Acts normal to the current (deformed) surface.
    //@todo traction load should be moved outside this method
    if ( iSurf != 0 ) {
        _error("computePressureForceAt: incompatible load surface must be 0 for this element");
    }

    FloatArray load;
    if ( ConstantPressureLoad * pLoad = dynamic_cast< ConstantPressureLoad * >( surfLoad ) ) {
        FloatMatrix gcov;
        FloatArray g1, g2;
        FloatArray lcoords(3);
        lcoords.at(1) = gp->giveCoordinate(1);
        lcoords.at(2) = gp->giveCoordinate(2);
        lcoords.at(3) = pLoad->giveLoadOffset();

        this->evalCovarBaseVectorsAt(lcoords, gcov, genEps);
        g1.beColumnOf(gcov, 1);
        g2.beColumnOf(gcov, 2);
        surfLoad->computeValueAt(load, tStep, lcoords, mode);        // pressure component
        traction.beVectorProductOf(g1, g2);        // normal vector (should not be normalized due to integraton in reference config.)
        traction.times( -load.at(1) );
    } else if ( dynamic_cast< ConstantSurfaceLoad * >( surfLoad ) ) {
        surfLoad->computeValueAt(traction, tStep, * ( gp->giveCoordinates() ), mode);        // traction vector
    } else {
        _error("computePressureForceAt: incompatible load type");
    }
}

void
Shell7Base :: evalCovarNormalAt(FloatArray &nCov, FloatArray &lCoords, FloatArray &genEpsC)
{
    FloatMatrix gcov;
    this->evalCovarBaseVectorsAt(lCoords, gcov, genEpsC);
    FloatArray g1, g2;
    g1.beColumnOf(gcov, 1);
    g2.beColumnOf(gcov, 2);
    nCov.beVectorProductOf(g1, g2);
    nCov.normalize();
}

void
Shell7Base :: evalInitialCovarNormalAt(FloatArray &nCov, FloatArray &lCoords)
{
    FloatMatrix Gcov;
    this->evalInitialCovarBaseVectorsAt(lCoords, Gcov);
    FloatArray G1, G2;
    G1.beColumnOf(Gcov, 1);
    G2.beColumnOf(Gcov, 2);
    nCov.beVectorProductOf(G1, G2);
    nCov.normalize();
}

void
Shell7Base :: computeTractionForce(FloatArray &answer, const int iEdge, BoundaryLoad *edgeLoad, TimeStep *tStep, ValueModeType mode)
{
    //
    IntegrationRule *iRule = specialIntegrationRulesArray [ 2 ];   // rule #3 for edge integration of distributed loads given in [*/m]
    GaussPoint *gp;

    FloatMatrix N, Q;
    FloatArray fT(7), components, lcoords;

    //BoundaryLoad :: CoordSystType coordSystType = edgeLoad->giveCoordSystMode();
    Load :: CoordSystType coordSystType = edgeLoad->giveCoordSystMode();

    FloatArray Nftemp(21), Nf(21);
    Nf.zero();
    for ( int i = 0; i < iRule->giveNumberOfIntegrationPoints(); i++ ) {
        gp = iRule->getIntegrationPoint(i);
        lcoords = * gp->giveCoordinates();

        edgeLoad->computeValueAt(components, tStep, lcoords, mode);
        FloatArray lCoords = * gp->giveCoordinates();
        this->edgeComputeNmatrixAt(lCoords, N);

        //if ( coordSystType ==  BoundaryLoad :: BL_UpdatedGlobalMode ) {
        if ( coordSystType ==  Load :: CST_UpdatedGlobal ) {
            // Updated global coord system
            FloatMatrix gcov;
            FloatArray lCoords = * gp->giveCoordinates();
            this->edgeEvalCovarBaseVectorsAt(lCoords, iEdge, gcov, tStep);
            Q.beTranspositionOf(gcov);

            FloatArray distrForces(3), distrMoments(3), t1, t2;
            distrForces.setValues( 3, components.at(1), components.at(2), components.at(3) );
            distrMoments.setValues( 3, components.at(4), components.at(5), components.at(6) );
            t1.beTProductOf(Q, distrForces);
            t2.beTProductOf(Q, distrMoments);
            fT.addSubVector(t1, 1);
            fT.addSubVector(t2, 4);
            fT.at(7) = components.at(7); // don't do anything with the 'gamma'-load

            //} else if( coordSystType == BoundaryLoad :: BL_GlobalMode ) {
        } else if ( coordSystType == Load :: CST_Global ) {
            // Undeformed global coord system
            for ( int i = 1; i <= 7; i++ ) {
                fT.at(i) = components.at(i);
            }
        } else {
            OOFEM_ERROR("Shell7Base :: computeTractionForce - does not support local coordinate system");
        }

        double dL = this->edgeComputeLengthAround(gp, iEdge);

        Nftemp.beTProductOf(N, fT * dL);
        Nf.add(Nftemp);
    }

    IntArray mask;
    this->giveEdgeDofMapping(mask, iEdge);
    answer.resize( Shell7Base :: giveNumberOfDofs() );
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
    FloatArray lcoords;
    lcoords = * gp->giveCoordinates();
    this->edgeEvalInitialCovarBaseVectorsAt(lcoords, iedge, G1, G3);
    detJ = G1.computeNorm();
    return detJ * gp->giveWeight();
}
#endif




// Transformation of bases and transformation matrices

# if 1
// general function
void
Shell7Base :: transInitialCartesianToInitialContravar(GaussPoint *gp, const FloatArray &VoightMatrix, FloatArray &answer)
{
    // Transform stress in cart system to curvilinear sytem
    // Uses Bond transformation matrix. No need to go from matrix to Voigt form and back.
    FloatMatrix Gcon;
    FloatArray lcoords;
    lcoords = * gp->giveCoordinates();
    this->evalInitialContravarBaseVectorsAt(lcoords, Gcon);
    FloatMatrix GE, M;
    GE.beTranspositionOf(Gcon); // Transformation matrix
    this->giveBondTransMatrix(M, GE);

    answer.beProductOf(M, VoightMatrix);
}

// general function
void
Shell7Base :: transInitialCartesianToInitialContravar(GaussPoint *gp, const FloatMatrix &stiffness, FloatMatrix &answer)
{
    // Transform stifness tensor in cart system to curvilinear sytem
    // Uses Bond transformation matrix. No need to go from matrix to Voigt form and back.
    FloatMatrix Gcon;
    FloatArray lcoords;
    lcoords = * gp->giveCoordinates();
    this->evalInitialContravarBaseVectorsAt(lcoords, Gcon);
    FloatMatrix GE, M, temp;
    GE.beTranspositionOf(Gcon); // Transformation matrix
    this->giveBondTransMatrix(M, GE);
    temp.beProductTOf(stiffness, M);
    answer.beProductOf(M, temp);
}

// general function
void
Shell7Base :: giveCoordTransMatrix(FloatMatrix &answer, FloatArray &g1, FloatArray &g2, FloatArray &g3,
                                   FloatArray &G1, FloatArray &G2, FloatArray &G3)
{   // Q_ij = dotproduct(g_i, G_j) = g_i^T*G_j
    // currently not in use
    //answer.beTProductOf(g,G);
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

// general function
void
Shell7Base :: giveCoordTransMatrix(FloatMatrix &answer, FloatArray &g1, FloatArray &g2, FloatArray &g3)
{
    // Assumes transformation from a cartesian system
    // Q_ij = dotproduct(g_i, E_j) with E1 = [1 0 0]^T etc
    // Q = g_i^T
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

// general function
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
        answer.times( this->giveCrossSection()->give(CS_Thickness, this->giveNode(node)->giveCoordinates(), this, false) );
    } else {
        answer.resize(0);
    }
}





void
Shell7Base :: ZZNodalRecoveryMI_computeNValProduct(FloatMatrix &answer, int layer, InternalStateType type,
                                                   TimeStep *tStep)
{  // evaluates N^T sigma over element volume
   // N(nsigma, nsigma*nnodes)
   // Definition : sigmaVector = N * nodalSigmaVector
    FloatArray stressVector, n;
    Element *elem  = this->ZZNodalRecoveryMI_giveElement();
    IntegrationRule *iRule = integrationRulesArray [ layer - 1 ];
    GaussPoint *gp;

    //int size = ZZNodalRecoveryMI_giveDofManRecordSize(type);
    int size = 6; ///@todo this is hard coded for stress recovery
    int numDofMans = 15;
    answer.resize(numDofMans, size);

    answer.zero();
    for ( int i = 0; i < iRule->giveNumberOfIntegrationPoints(); i++ ) {
        gp  = iRule->getIntegrationPoint(i);
        double dV = this->computeVolumeAroundLayer(gp, layer);

        if ( !elem->giveIPValue(stressVector, gp, type, tStep) ) {
            stressVector.resize(size);
            stressVector.zero();
        }

        this->ZZNodalRecoveryMI_ComputeEstimatedInterpolationMtrx(n, gp, type);
        for ( int j = 1; j <= numDofMans; j++ ) {
            for ( int k = 1; k <= size; k++ ) {
                answer.at(j, k) += n.at(j) * stressVector.at(k) * dV;
            }
        }
    }
}

void
Shell7Base :: ZZNodalRecoveryMI_computeNNMatrix(FloatArray &answer, int layer, InternalStateType type)
{
    //
    // Returns NTN matrix (lumped) for Zienkiewicz-Zhu
    // The size of N mtrx is (nstresses, nnodes*nstreses)
    // Definition : sigmaVector = N * nodalSigmaVector
    //
    //double volume = 0.0;
    FloatMatrix fullAnswer;
    FloatArray n;
    IntegrationRule *iRule = integrationRulesArray [ layer - 1 ];
    GaussPoint *gp;

    int size = 15;
    fullAnswer.resize(size, size);
    fullAnswer.zero();
    for ( int i = 0; i < iRule->giveNumberOfIntegrationPoints(); i++ ) {
        gp  = iRule->getIntegrationPoint(i);
        double dV = this->computeVolumeAroundLayer(gp, layer);
        this->ZZNodalRecoveryMI_ComputeEstimatedInterpolationMtrx(n, gp, type);
        //fullAnswer.plusDyadSymmUpper(n, n, dV);
        fullAnswer.plusDyadSymmUpper(n, dV);
        //volume += dV;
    }


    fullAnswer.symmetrized();
    answer.resize(size);
    for ( int i = 1; i <= size; i++ ) {
        double sum = 0.0;
        for ( int j = 1; j <= size; j++ ) {
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
    FEInterpolation *interpol = static_cast< FEInterpolation * >( & this->interpolationForExport );

    // test if underlying element provides interpolation
    if ( interpol ) {
        ///@todo fix this whole compostie recovery thing in a better way
        //if ( !this->giveIPValueSize(type, gp) ) {
        //    OOFEM_ERROR3("ZZNodalRecoveryMI_computeNNMatrix: Element %d not supporting type %d", this->giveNumber(), type);
        //    return;
        //}

        interpol->evalN( answer, * gp->giveCoordinates(), FEIElementGeometryWrapper(this) );
    } else {
        // ok default implementation can not work, as element is not providing valid interpolation
        // to resolve this, one can overload this method for element implementing ZZNodalRecoveryModelInterface
        // or element should provide interpolation.
        OOFEM_ERROR2( "ZZNodalRecoveryMI_computeNNMatrix: Element %d not providing valid interpolation", this->giveNumber() );
    }
}


void
Shell7Base :: ZZNodalRecoveryMI_recoverValues(std :: vector< FloatArray > &recoveredValues, int layer, InternalStateType type, TimeStep *tStep)
{
    // ZZ recovery
    FloatArray nnMatrix;
    FloatMatrix nValProd;
    ZZNodalRecoveryMI_computeNValProduct(nValProd, layer, type, tStep);
    ZZNodalRecoveryMI_computeNNMatrix(nnMatrix, layer, type);
    int recoveredSize = nValProd.giveNumberOfColumns();
    int numNodes = nValProd.giveNumberOfRows();
    recoveredValues.resize(numNodes);

    for ( int i = 1; i <= numNodes; i++ ) {
        //recoveredValues[i-1].resize(recoveredSize);
        FloatArray temp(6);
        recoveredValues [ i - 1 ].resize(9);

        for ( int j = 1; j <= recoveredSize; j++ ) {
            temp.at(j) = nValProd.at(i, j) / nnMatrix.at(i);
            //recoveredValues[i-1].at(j) = nValProd.at(i,j)/nnMatrix.at(i);
        }

        recoveredValues [ i - 1 ] = convV6ToV9Stress(temp);
    }
}

#endif




// Computation of solution vectors

#if 1

void
Shell7Base :: giveSolutionVector(FloatArray &answer, const IntArray &dofIdArray, TimeStep *tStep)
{
    // Returns the solution vector corresponding to all the dofs
    FloatArray temp;
    computeVectorOfDofIDs(dofIdArray, VM_Total, tStep, temp);
    answer.resize( Shell7Base :: giveNumberOfDofs() );
    answer.zero();
    answer.assemble( temp, this->giveOrdering(AllInv) );
}


void
Shell7Base :: computeVectorOfDofIDs(const IntArray &dofIdArray, ValueModeType u, TimeStep *tStep, FloatArray &answer)
{
    // Routine to extract the solution vector for an element given an dofid array.
    // Size will be numberOfDofs and if a certain dofId does not exist a zero is used as value.
    // This method may e.g. be used to obtain the enriched part of the solution vector
    ///@todo: generalize so it can be used by all XFEM elements
    answer.resize( Shell7Base :: giveNumberOfDofs() );
    answer.zero();
    int k = 0;
    for ( int i = 1; i <= numberOfDofMans; i++ ) {
        DofManager *dMan = this->giveDofManager(i);
        for ( int j = 1; j <= dofIdArray.giveSize(); j++ ) {
            //k++;
            if ( dMan->hasDofID( ( DofIDItem ) dofIdArray.at(j) ) ) {
                Dof *d = dMan->giveDofWithID( dofIdArray.at(j) );
                /// @todo: will fail if any other dof then gamma is excluded from enrichment
                /// since I only add "j". Instead I should skip certain dof numbers when incrementing
                answer.at(k + j) = d->giveUnknown(VM_Total, tStep);
            }
        }
        k += 7;
    }
}

void
Shell7Base :: temp_computeBoundaryVectorOf(IntArray &dofIdArray, int boundary, ValueModeType u, TimeStep *tStep, FloatArray &answer)
{
    ///@todo: NOT CHECKED!!!
    // Routine to extract vector given an array of dofid items
    // If a certain dofId does not exist a zero is used as value

    IntArray bNodes;
    this->fei->computeLocalEdgeMapping(bNodes, boundary);

    answer.resize( dofIdArray.giveSize() * bNodes.giveSize() );
    answer.zero();
    int k = 0;
    for ( int i = 1; i <= bNodes.giveSize(); i++ ) {
        DofManager *dMan = this->giveDofManager( bNodes.at(i) );
        for ( int j = 1; j <= dofIdArray.giveSize(); j++ ) {
            Dof *d = dMan->giveDof(j);
            k++;
            if ( dMan->hasDofID( ( DofIDItem ) dofIdArray.at(j) ) ) {
                answer.at(k) = d->giveUnknown(VM_Total, tStep);
            }
        }
    }
}

void
Shell7Base :: giveUpdatedSolutionVector(FloatArray &answer, TimeStep *tStep)
{
    // Computes updated solution as: x = X + dx, m = M + dM, gam = 0 + dgam
    // This is because the element formulation is in terms of placement and not displacement.
    this->giveInitialSolutionVector(answer); // X & M
    FloatArray temp;
    int dummy = 0;
    IntArray dofIdArray;
    Shell7Base :: giveDofManDofIDMask(dummy, EID_MomentumBalance, dofIdArray);
    this->computeVectorOfDofIDs(dofIdArray, VM_Total, tStep, temp);

    answer.assemble( temp, this->giveOrdering(AllInv) );
}


void
Shell7Base :: giveInitialSolutionVector(FloatArray &answer)
{
    // Gives the initial values of X, M and gamma which definies the initial position
    answer.resize( Shell7Base :: giveNumberOfDofs() );
    answer.zero();

    int ndofs_xm  = this->giveNumberOfFieldDofs(Midplane);
    // Reference position and directors
    for ( int i = 1, j = 0; i <= this->giveNumberOfDofManagers(); i++, j += 3 ) {
        FloatArray *Xi = this->giveNode(i)->giveCoordinates();
        FloatArray Mi  = this->giveInitialNodeDirector(i);
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
    int dummy = 0;
    IntArray dofIdArray;
    Shell7Base :: giveDofManDofIDMask(dummy, EID_MomentumBalance, dofIdArray);
    this->temp_computeBoundaryVectorOf(dofIdArray, iedge, VM_Total, tStep, temp);
    answer.assemble( temp, this->giveOrdering(EdgeInv) );
}

void
Shell7Base :: edgeGiveInitialSolutionVector(FloatArray &answer, const int iedge)
{
    answer.resize( this->giveNumberOfEdgeDofs() );
    answer.zero();
    IntArray edgeNodes;
    this->fei->computeLocalEdgeMapping(edgeNodes, iedge);
    int ndofs_x = this->giveFieldSize(Midplane) * edgeNodes.giveSize();
    for ( int i = 1, j = 0; i <= edgeNodes.giveSize(); i++, j += 3 ) {
        FloatArray *Xi = this->giveNode( edgeNodes.at(i) )->giveCoordinates();
        FloatArray Mi  = this->giveInitialNodeDirector( edgeNodes.at(i) );
        answer.at(1 + j) = Xi->at(1);
        answer.at(2 + j) = Xi->at(2);
        answer.at(3 + j) = Xi->at(3);
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
Shell7Base :: computeGeneralizedStrainVectorNew(FloatArray &answer, const FloatArray &solVec, const FloatMatrix &B) {
    // Returns an array genEps = [dxdxi, dmdxi, m, dgamdxi, gam]^T with size 18
    answer.resize(18);
    answer.zero();


    /*    18   18   6
     * 6 [B_u   0   0
     * 6   0   B_w  0
     * 3   0   N_w  0
     * 2   0    0  B_gam
     * 1   0    0  N_gam]
     */
    int ndofs_xm  = this->giveNumberOfFieldDofs(Midplane);
    int ndofs_gam = this->giveNumberOfFieldDofs(InhomStrain);
    for ( int i = 1; i <= ndofs_xm; i++ ) {
        for ( int j = 1; j <= 6; j++ ) { // 3dofs*2
            answer.at(j)      += B.at(j, i) * solVec.at(i);                                  // dx/dxi
            answer.at(6 + j)  += B.at(6 + j, i + ndofs_xm) * solVec.at(i + ndofs_xm);        // dm/dxi
        }
    }

    for ( int i = 1; i <= ndofs_xm; i++ ) {
        for ( int j = 1; j <= 3; j++ ) { // 3 dofs
            answer.at(12 + j) += B.at(12 + j, i + ndofs_xm) * solVec.at(i + ndofs_xm);        // m
        }
    }

    for ( int i = 1; i <= ndofs_gam; i++ ) {
        for ( int j = 1; j <= 2; j++ ) { // 2 dofs
            answer.at(15 + j) += B.at(15 + j, i + ndofs_xm * 2) * solVec.at(i + ndofs_xm * 2);    // dgamma/dxi
        }
    }

    for ( int i = 1; i <= ndofs_gam; i++ ) { //1 dof
        answer.at(18) += B.at(18, i + ndofs_xm * 2) * solVec.at(i + ndofs_xm * 2);            // gamma
    }
}

// two copies of the same function?
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


void
Shell7Base :: giveUnknownsAt(FloatArray &lcoords, FloatArray &solVec, FloatArray &x, FloatArray &m, double gam, TimeStep *tStep)
{
    // returns the unknowns evaluated at a point (xi1, xi2, xi3)
    FloatArray N;
    this->fei->evalN( N, lcoords, FEIElementGeometryWrapper(this) );

    x.resize(3);
    x.zero();
    m.resize(3);
    m.zero();
    int ndofs_xm  = this->giveNumberOfFieldDofs(Midplane);
    for ( int i = 1, j = 0; i <= this->giveNumberOfDofManagers(); i++, j += 3 ) {
        x.at(1) += N.at(i) * solVec.at(1 + j);
        x.at(2) += N.at(i) * solVec.at(2 + j);
        x.at(3) += N.at(i) * solVec.at(3 + j);
        m.at(1) += N.at(i) * solVec.at(ndofs_xm + 1 + j);
        m.at(2) += N.at(i) * solVec.at(ndofs_xm + 2 + j);
        m.at(3) += N.at(i) * solVec.at(ndofs_xm + 3 + j);
        gam     += N.at(i) * solVec.at(2 * ndofs_xm + i);
    }
}


#endif




// N and B matrices

#if 1


void
Shell7Base :: edgeComputeNmatrixAt(FloatArray &lcoords, FloatMatrix &answer)
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
Shell7Base :: edgeComputeBmatrixAt(FloatArray &lcoords, FloatMatrix &answer, int li, int ui)
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
Shell7Base :: computeBmatrixAt(FloatArray &lcoords, FloatMatrix &answer, int li, int ui)
{
    // Returns the  matrix {B} of the receiver, evaluated at gp. Such that
    // B*a = [dxbar_dxi, dwdxi, w, dgamdxi, gam]^T, where a is the vector of unknowns

    int ndofs = Shell7Base :: giveNumberOfDofs();
    int ndofs_xm  = this->giveNumberOfFieldDofs(Midplane);
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
    int ndofs_xm  = this->giveNumberOfFieldDofs(Midplane);
    answer.resize(7, ndofs);
    answer.zero();
    FloatArray N;
    this->fei->evalN( N, iLocCoord, FEIElementGeometryWrapper(this) );
    //N.printYourself();
    //lcoords.printYourself();
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



#endif


// VTK export
#if 1
void
Shell7Base :: vtkEvalInitialGlobalCoordinateAt(FloatArray &localCoords, int layer, FloatArray &globalCoords)
{
    double zeta = giveGlobalZcoordInLayer(localCoords.at(3), layer);
    FloatArray N;
    this->fei->evalN( N, localCoords, FEIElementGeometryWrapper(this) );

    globalCoords.resize(3);
    globalCoords.zero();
    for ( int i = 1; i <= this->giveNumberOfDofManagers(); i++ ) {
        FloatArray &xbar = * this->giveNode(i)->giveCoordinates();
        FloatArray M = this->giveInitialNodeDirector(i);
        globalCoords += N.at(i) * ( xbar + zeta * M );
    }
}


void
Shell7Base :: vtkEvalUpdatedGlobalCoordinateAt(FloatArray &localCoords, int layer, FloatArray &globalCoords, TimeStep *tStep)
{
    FloatArray solVec;
    this->giveUpdatedSolutionVector(solVec, tStep);
    FloatArray x, m;
    double gam = 0;
    this->giveUnknownsAt(localCoords, solVec, x, m, gam, tStep);
    double zeta = giveGlobalZcoordInLayer(localCoords.at(3), layer);
    double fac = ( zeta + 0.5 * gam * zeta * zeta );
    globalCoords = x + fac * m;
}


 #if 0
void
Shell7Base :: giveCompositeExportData(CompositeCell &compositeCell, IntArray &primaryVarsToExport, IntArray &cellVarsToExport, VTKPiece &vtkPiece, TimeStep *tStep)
{
    int numCells = this->layeredCS->giveNumberOfLayers();
    const int numCellNodes  = 15; // quadratic wedge
    int numTotalNodes = numCellNodes * numCells;
    vtkPiece.elCellTypes.resize(numCells);
    vtkPiece.elOffsets.resize(numCells);
    vtkPiece.nodeCoords.resize(numTotalNodes);


    compositeCell.elements.resize(numCells);
    compositeCell.numSubEl = numCells;
    compositeCell.numTotalNodes = numTotalNodes;

    int val    = 0;
    int offset = 0;
    // Compute fictious node coords
    int nodeNum = 0;
    for ( int layer = 1; layer <= numCells; layer++ ) {
        Cell &el = compositeCell.elements [ layer - 1 ];

        // Node coordinates
        this->giveFictiousNodeCoordsForExport(el.nodeCoords, layer);

        for ( int node = 1; node <= numCellNodes; node++ ) {
            vtkPiece.nodeCoords [ nodeNum ] = el.nodeCoords [ node - 1 ];
            nodeNum += 1;
        }

        // Connectivity
        el.connectivity.resize(numCellNodes);
        vtkPiece.connectivity.resize(numTotalNodes);
        for ( int i = 1; i <= numCellNodes; i++ ) {
            el.connectivity.at(i) = val++;
            vtkPiece.connectivity.at(i) = val++;
        }

        // Offset
        offset += numCellNodes;
        el.offset = offset;
        vtkPiece.elOffsets.at(layer) = offset;

        // Cell types
        el.cellType = 26; // Quadratic wedge
        vtkPiece.elCellTypes.at(layer) = 26;

        // Export nodal variables
        el.nodeVars.resize( primaryVarsToExport.giveSize() + cellVarsToExport.giveSize() );



        vtkPiece.nodeVars.resize( primaryVarsToExport.giveSize() );
        int nodeVarNum = 0;
        for ( int i = 1; i <= primaryVarsToExport.giveSize(); i++ ) {
            UnknownType type = ( UnknownType ) primaryVarsToExport.at(i);

            vtkPiece.nodeVars [ i - 1 ].resize(numTotalNodes);
            if ( type == DisplacementVector ) {
                std :: vector< FloatArray >updatedNodeCoords;
                giveFictiousUpdatedNodeCoordsForExport(updatedNodeCoords, layer, tStep);
                FloatArray u(3);
                el.nodeVars [ nodeVarNum ].resize(numCellNodes);
                for ( int j = 1; j <= numCellNodes; j++ ) {
                    u = updatedNodeCoords [ j - 1 ];
                    u.subtract(el.nodeCoords [ j - 1 ]);
                    el.nodeVars [ nodeVarNum ] [ j - 1 ].resize(3);
                    el.nodeVars [ nodeVarNum ] [ j - 1 ] = u;

                    //vtkPiece.nodeVars[i-1][inode-1] = valueArray;
                }
            } else {
                ZZNodalRecoveryMI_recoverValues(el.nodeVars [ i - 1 ], layer, ( InternalStateType ) 1, tStep);
            }
            nodeVarNum++;
        }

        vtkPiece.nodeVarsFromIS.resize( cellVarsToExport.giveSize() );
        for ( int i = 1; i <= cellVarsToExport.giveSize(); i++ ) {
            InternalStateType type = ( InternalStateType ) cellVarsToExport.at(i);
            ;

            this->recoverShearStress(tStep);
            recoverValuesFromIP(el.nodeVars [ nodeVarNum ], layer, type, tStep);
            //ZZNodalRecoveryMI_recoverValues(el.nodeVars[nodeVarNum], layer, type, tStep);

            nodeVarNum++;
        }


        // Export element (cell) variables
        IntegrationRule *iRuleL = integrationRulesArray [ layer - 1 ];
        int numElVars = cellVarsToExport.giveSize();
        el.elVars.resize(numElVars);
        for ( int i = 1; i <= numElVars; i++ ) {
            InternalStateType type = ( InternalStateType ) cellVarsToExport.at(i);
            GaussPoint *gp;

            // stress
            FloatArray average(6), temp;
            average.zero();
            double gptot = 0.0;
            // Avergae gp values
            for ( int j = 0; j < iRuleL->giveNumberOfIntegrationPoints(); j++ ) {
                gp = iRuleL->getIntegrationPoint(j);
                this->giveIPValue(temp, gp, type, tStep);
                gptot += gp->giveWeight();
                average.add(gp->giveWeight(), temp);
            }
            average.times(1. / gptot);

            el.elVars [ i - 1 ] = convV6ToV9Stress(average);
        }
    }
}

 #else

void
Shell7Base :: giveCompositeExportData(VTKPiece &vtkPiece, IntArray &primaryVarsToExport, IntArray &internalVarsToExport, IntArray cellVarsToExport, TimeStep *tStep)
{
    int numCells = this->layeredCS->giveNumberOfLayers();
    const int numCellNodes  = 15; // quadratic wedge
    int numTotalNodes = numCellNodes * numCells;

    //int numCells = this->layeredCS->giveNumberOfLayers();
    //const int numCellNodes  = 15; // quadratic wedge
    //int numTotalNodes = numCellNodes*numCells;
    vtkPiece.setNumberOfCells(numCells);
    vtkPiece.setNumberOfNodes(numTotalNodes);

    //vtkPiece.elCellTypes.resize(numCells);
    //vtkPiece.elOffsets.resize(numCells);
    //vtkPiece.nodeCoords.resize(numTotalNodes);
    //vtkPiece.connectivity.resize(numCells);

    std :: vector< FloatArray >nodeCoords;
    int val    = 1;
    int offset = 0;
    IntArray nodes(numCellNodes);

    // Compute fictious node coords
    int nodeNum = 1;
    for ( int layer = 1; layer <= numCells; layer++ ) {
        // Node coordinates
        this->giveFictiousNodeCoordsForExport(nodeCoords, layer);

        for ( int node = 1; node <= numCellNodes; node++ ) {
            vtkPiece.setNodeCoords(nodeNum, nodeCoords [ node - 1 ]);
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

    std :: vector< FloatArray >updatedNodeCoords;
    FloatArray u(3);
    std :: vector< FloatArray >values;
    for ( int fieldNum = 1; fieldNum <= primaryVarsToExport.giveSize(); fieldNum++ ) {
        UnknownType type = ( UnknownType ) primaryVarsToExport.at(fieldNum);
        nodeNum = 1;
        for ( int layer = 1; layer <= numCells; layer++ ) {
            if ( type == DisplacementVector ) { // compute displacement as u = x - X
                this->giveFictiousNodeCoordsForExport(nodeCoords, layer);
                this->giveFictiousUpdatedNodeCoordsForExport(updatedNodeCoords, layer, tStep);
                for ( int j = 1; j <= numCellNodes; j++ ) {
                    u = updatedNodeCoords [ j - 1 ];
                    u.subtract(nodeCoords [ j - 1 ]);
                    vtkPiece.setPrimaryVarInNode(fieldNum, nodeNum, u);
                    nodeNum += 1;
                }
            } else {
                ZZNodalRecoveryMI_recoverValues(values, layer, ( InternalStateType ) 1, tStep); // does not work well - fix
                for ( int j = 1; j <= numCellNodes; j++ ) {
                    vtkPiece.setPrimaryVarInNode(fieldNum, nodeNum, values [ j - 1 ]);
                    nodeNum += 1;
                }
            }
        }
    }

    // Export nodal variables from internal fields

    vtkPiece.setNumberOfInternalVarsToExport(internalVarsToExport.giveSize(), numTotalNodes);
    for ( int fieldNum = 1; fieldNum <= internalVarsToExport.giveSize(); fieldNum++ ) {
        InternalStateType type = ( InternalStateType ) internalVarsToExport.at(fieldNum);
        nodeNum = 1;
        this->recoverShearStress(tStep);
        for ( int layer = 1; layer <= numCells; layer++ ) {
            recoverValuesFromIP(values, layer, type, tStep);
            for ( int j = 1; j <= numCellNodes; j++ ) {
                vtkPiece.setInternalVarInNode(fieldNum, nodeNum, values [ j - 1 ]);
                //ZZNodalRecoveryMI_recoverValues(el.nodeVars[fieldNum], layer, type, tStep);
                nodeNum += 1;
            }
        }
    }


    // Export cell variables
    FloatArray average;
    vtkPiece.setNumberOfCellVarsToExport(cellVarsToExport.giveSize(), numCells);
    for ( int i = 1; i <= cellVarsToExport.giveSize(); i++ ) {
        InternalStateType type = ( InternalStateType ) cellVarsToExport.at(i);
        ;

        for ( int layer = 1; layer <= numCells; layer++ ) {
            IntegrationRule *iRuleL = integrationRulesArray [ layer - 1 ];
            VTKXMLExportModule :: computeIPAverage(average, iRuleL, this, type, tStep);
            vtkPiece.setCellVar( i, layer, convV6ToV9Stress(average) );
        }
    }
}

 #endif




void
Shell7Base :: recoverValuesFromIP(std :: vector< FloatArray > &recoveredValues, int layer, InternalStateType type, TimeStep *tStep)
{
    // recover nodal values by coosing the ip closest to the node

    //FEInterpolation *interpol = static_cast< FEInterpolation * >( &this->interpolationForExport );

    // composite element interpolator

    int numNodes = 15;
    recoveredValues.resize(numNodes);

    IntegrationRule *iRule = integrationRulesArray [ layer - 1 ];
    IntegrationPoint *ip;

    // Find closest ip to the nodes
    IntArray closestIPArray(numNodes);
    FloatArray nodeCoords, ipCoords, ipValues;
    FloatArray nodeLocalXi1Coords, nodeLocalXi2Coords, nodeLocalXi3Coords;
    giveLocalNodeCoordsForExport(nodeLocalXi1Coords, nodeLocalXi2Coords, nodeLocalXi3Coords);

    for ( int i = 1; i <= numNodes; i++ ) {
        nodeCoords.setValues( 3, nodeLocalXi1Coords.at(i), nodeLocalXi2Coords.at(i), nodeLocalXi3Coords.at(i) );
        double distOld = 3.0; // should not be larger
        for ( int j = 0; j < iRule->giveNumberOfIntegrationPoints(); j++ ) {
            ip = iRule->getIntegrationPoint(j);
            ipCoords = * ip->giveCoordinates();
            double dist = nodeCoords.distance(ipCoords);
            if ( dist < distOld ) {
                closestIPArray.at(i) = j;
                distOld = dist;
            }
        }
    }

    // recover ip values
    for ( int i = 1; i <= numNodes; i++ ) {
        ip = iRule->getIntegrationPoint( closestIPArray.at(i) );
        this->giveIPValue(ipValues, ip, type, tStep);
        recoveredValues [ i - 1 ].resize(9);
        recoveredValues [ i - 1 ] = convV6ToV9Stress(ipValues);
    }
}


void
Shell7Base :: recoverShearStress(TimeStep *tStep)
{
    // Recover shear stresses at ip by numerical integration of the momentum balance through the thickness
    std :: vector< FloatArray >recoveredValues;
    int numberOfLayers = this->layeredCS->giveNumberOfLayers();     // conversion of types
    IntegrationRule *iRuleThickness = specialIntegrationRulesArray [ 0 ];
    FloatArray dS, Sold;
    FloatMatrix B, Smat(2, 6); // 2 stress components * num of in plane ip ///@todo generalize
    Smat.zero();
    FloatArray Tcon(6), Trec(6);
    Tcon.zero();
    Trec.zero();

    for ( int layer = 1; layer <= numberOfLayers; layer++ ) {
        IntegrationRule *iRuleL = integrationRulesArray [ layer - 1 ];
        this->recoverValuesFromIP(recoveredValues, layer, IST_StressTensor, tStep);
        //this->ZZNodalRecoveryMI_recoverValues(recoveredValues, layer, IST_StressTensor, tStep);
        double thickness = this->layeredCS->giveLayerThickness(layer);

        //set up vector of stresses in the ip's = [S1_xx, S1_yy, S1_xy, ..., Sn_xx, Sn_yy, Sn_xy]
        int numNodes = 15;
        FloatArray aS(numNodes * 3);
        for ( int j = 1, pos = 0; j <= numNodes; j++, pos += 3 ) {
            aS.at(pos + 1) = recoveredValues [ j - 1 ].at(1);   // S_xx
            aS.at(pos + 2) = recoveredValues [ j - 1 ].at(2);   // S_yy
            aS.at(pos + 3) = recoveredValues [ j - 1 ].at(6);   // S_xy
        }
        int numInPlaneIP = 6;

        for ( int i = 0; i < iRuleThickness->giveNumberOfIntegrationPoints(); i++ ) {
            double dz = thickness * iRuleThickness->getIntegrationPoint(i)->giveWeight();

            for ( int j = 0; j < numInPlaneIP; j++ ) {
                int point = i * numInPlaneIP + j; // integration point number
                GaussPoint *gp = iRuleL->getIntegrationPoint(point);

                this->computeBmatrixForStressRecAt(* gp->giveCoordinates(), B, layer);
                dS.beProductOf( B, aS * ( -dz ) ); // stress increment

                //StructuralMaterialStatus* status = dynamic_cast< StructuralMaterialStatus* > (gp->giveMaterialStatus(1));
                StructuralMaterialStatus *status = dynamic_cast< StructuralMaterialStatus * >( gp->giveMaterialStatus() );
                Sold = status->giveStressVector();

                Smat.at(1, j + 1) += dS.at(1); // add increment from each level
                Smat.at(2, j + 1) += dS.at(2);

                //Tcon.at(j+1) += Sold.at(5)*dz;

                // Replace old stresses with  - this should probably not be done as it may affect the convergence in a nonlinear case
                Sold.at(5) = Smat.at(1, j + 1); // S_xz
                Sold.at(4) = Smat.at(2, j + 1); // S_yz


                status->letStressVectorBe(Sold);
                //Trec.at(j+1) += Sold.at(5)*dz;
            }
        }
    }

    //printf("Shear force from continuum stress [N/M] \n");
    //Tcon.printYourself();
    //printf("Shear force from recovered stress [N/M] \n");
    //Trec.printYourself();
}


void
Shell7Base :: computeBmatrixForStressRecAt(FloatArray &lcoords, FloatMatrix &answer, int layer)
{
    // Returns the  special matrix {B} of the receiver, evaluated at gp. Such that
    // B*a = [dS_xx/dx + dS_xy/dy, dS_yx/dx + dS_yy/dy ]^T, where a is the vector of in plane
    // stresses [S_xx, S_yy, S_xy]

    // set up virtual cell geometry for an qwedge
    const int numNodes = 15;
    std :: vector< FloatArray >nodes;
    giveFictiousNodeCoordsForExport(nodes, layer);

    //int VTKWedge2EL [] = { 3, 1, 2, 6, 4, 5, 9, 7, 8, 12, 10, 11, 15, 13, 14 };
    int VTKWedge2EL [] = {
        1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15
    };
    FloatArray *coords [ numNodes ];


    for ( int i = 1; i <= numNodes; i++ ) {
        //int pos = i;
        int pos = VTKWedge2EL [ i - 1 ];
        //coords[ i - 1 ] = new FloatArray();
        coords [ i - 1 ] = & nodes [ pos - 1 ];
    }

    FEInterpolation *interpol = static_cast< FEInterpolation * >( & this->interpolationForExport );
    FloatMatrix dNdx;
    interpol->evaldNdx( dNdx, lcoords, FEIVertexListGeometryWrapper(numNodes, ( const FloatArray ** ) coords) );





    /*
     * 1 [d/dx  0   d/dy
     * 1   0   d/dy d/dx]
     */
    int ndofs = numNodes * 3;
    answer.resize(2, ndofs);
    for ( int i = 1, j = 0; i <= numNodes; i++, j += 3 ) {
        answer.at(1, j + 1) = dNdx.at(i, 1);
        answer.at(1, j + 3) = dNdx.at(i, 2);
        answer.at(2, j + 2) = dNdx.at(i, 2);
        answer.at(2, j + 3) = dNdx.at(i, 1);
    }

    // how to destruct?
    for ( int i = 1; i <= numNodes; i++ ) {
        //coords [ i - 1 ]->~FloatArray();
    }
}





void
Shell7Base :: giveFictiousNodeCoordsForExport(std :: vector< FloatArray > &nodes, int layer)
{
    // compute fictious node coords
    FloatArray nodeLocalXi1Coords, nodeLocalXi2Coords, nodeLocalXi3Coords;
    giveLocalNodeCoordsForExport(nodeLocalXi1Coords, nodeLocalXi2Coords, nodeLocalXi3Coords);
    nodes.resize(15);
    for ( int i = 1; i <= 15; i++ ) {
        FloatArray coords, localCoords(3);
        localCoords.at(1) = nodeLocalXi1Coords.at(i);
        localCoords.at(2) = nodeLocalXi2Coords.at(i);
        localCoords.at(3) = nodeLocalXi3Coords.at(i);

        this->vtkEvalInitialGlobalCoordinateAt(localCoords, layer, coords);
        nodes [ i - 1 ].resize(3);
        nodes [ i - 1 ] = coords;
        //nodes[i-1] = localCoords;
    }
}

void
Shell7Base :: giveFictiousUpdatedNodeCoordsForExport(std :: vector< FloatArray > &nodes, int layer, TimeStep *tStep)
{
    // compute fictious node coords
    FloatArray nodeLocalXi1Coords, nodeLocalXi2Coords, nodeLocalXi3Coords;
    giveLocalNodeCoordsForExport(nodeLocalXi1Coords, nodeLocalXi2Coords, nodeLocalXi3Coords);
    nodes.resize(15);
    for ( int i = 1; i <= 15; i++ ) {
        FloatArray coords, localCoords(3);
        localCoords.at(1) = nodeLocalXi1Coords.at(i);
        localCoords.at(2) = nodeLocalXi2Coords.at(i);
        localCoords.at(3) = nodeLocalXi3Coords.at(i);

        this->vtkEvalUpdatedGlobalCoordinateAt(localCoords, layer, coords, tStep);
        nodes [ i - 1 ].resize(3);
        nodes [ i - 1 ] = coords;
    }
}


void
Shell7Base :: giveLocalNodeCoordsForExport(FloatArray &nodeLocalXi1Coords, FloatArray &nodeLocalXi2Coords, FloatArray &nodeLocalXi3Coords) {
    // Local coords for a quadratic wedge element (VTK cell type 26)
    double z = 0.999;
    nodeLocalXi1Coords.setValues(15, 1., 0., 0., 1., 0., 0., .5, 0., .5, .5, 0., .5, 1., 0., 0.);
    nodeLocalXi2Coords.setValues(15, 0., 1., 0., 0., 1., 0., .5, .5, 0., .5, .5, 0., 0., 1., 0.);
    nodeLocalXi3Coords.setValues(15, -z, -z, -z,  z,  z,  z, -z, -z, -z,  z,  z,  z, 0., 0., 0.);
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


int
Shell7Base :: giveVoigtIndex(int ind1, int ind2)
{
    // Returns the Voigt index corresponding to two given tensor indices.
    //@todo which component ordering does this assume?
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
}

void
Shell7Base :: computeStrainVectorInLayer(FloatArray &answer, const FloatArray &masterGpStrain, GaussPoint *masterGp, GaussPoint *slaveGp, TimeStep *tStep)
{
    OOFEM_ERROR("Shell7Base :: computeStrainVectorInLayer - Should not be called! Not meaningful for this element.");
}



void
Shell7Base :: computeInterLaminarStressesAt(int interfaceNum, TimeStep *tStep, std :: vector< FloatArray > &interLamStresses)
{
    // Get integration rules on the upper and lower sides of the interface


    LayeredIntegrationRule *irLower = static_cast< LayeredIntegrationRule * >( this->giveIntegrationRule(interfaceNum - 1) ); // index from 0
    LayeredIntegrationRule *irUpper = static_cast< LayeredIntegrationRule * >( this->giveIntegrationRule(interfaceNum) );
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

        this->evalInitialCovarNormalAt( nCov, * ip->giveCoordinates() );
        vSMean = 0.5 * ( vSUpper + vSLower );
        SMean.beMatrixFormOfStress(vSMean);
        stressVec.beProductOf(SMean, nCov);
        interLamStresses.at(i - 1).resize( stressVec.giveSize() );
        interLamStresses.at(i - 1) = stressVec;
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
} // end namespace oofem
