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

#include "../sm/Elements/Bars/truss1d.h"
#include "../sm/CrossSections/structuralcrosssection.h"
#include "fei1dlin.h"
#include "node.h"
#include "material.h"
#include "gausspoint.h"
#include "gaussintegrationrule.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "intarray.h"
#include "mathfem.h"
#include "classfactory.h"

#ifdef __OOFEG
 #include "oofeggraphiccontext.h"
#endif

namespace oofem {
REGISTER_Element(Truss1d);

FEI1dLin Truss1d :: interp(1); // Initiates the static interpolator


Truss1d :: Truss1d(int n, Domain *aDomain) :
    StructuralElement(n, aDomain),
    ZZNodalRecoveryModelInterface(this), NodalAveragingRecoveryModelInterface(),
    SpatialLocalizerInterface(this),
    ZZErrorEstimatorInterface(this),
    HuertaErrorEstimatorInterface()
{
    numberOfDofMans = 2;
}


void Truss1d :: computeGaussPoints()
// Sets up the array of Gauss Points of the receiver.
{
    if ( integrationRulesArray.size() == 0 ) {
        integrationRulesArray.resize( 1 );
        integrationRulesArray [ 0 ].reset( new GaussIntegrationRule(1, this, 1, 2) );
        this->giveCrossSection()->setupIntegrationPoints(* integrationRulesArray [ 0 ], 1, this);
    }
}

FEInterpolation *
Truss1d :: giveInterpolation() const { return & interp; }

void
Truss1d :: computeLumpedMassMatrix(FloatMatrix &answer, TimeStep *tStep)
// Returns the lumped mass matrix of the receiver. This expression is
// valid in both local and global axes.
{
    GaussPoint *gp = integrationRulesArray [ 0 ]->getIntegrationPoint(0);
    double density = this->giveStructuralCrossSection()->give('d', gp);
    double halfMass = density * this->giveCrossSection()->give(CS_Area, gp) * this->computeLength() * 0.5;
    answer.resize(2, 2);
    answer.zero();
    answer.at(1, 1) = halfMass;
    answer.at(2, 2) = halfMass;
}


void
Truss1d :: computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int li, int ui)
//
// Returns linear part of geometrical equations of the receiver at gp.
// Returns the linear part of the B matrix
//
{
    FloatMatrix dN;
    this->interp.evaldNdx( dN, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );
    answer.beTranspositionOf(dN); ///@todo It would be more suitable to follow the column-major version as done in FEI-classes
}


void
Truss1d :: computeBHmatrixAt(GaussPoint *gp, FloatMatrix &answer)
//
// Returns the [1x2] displacement gradient matrix {BH} of the receiver,
// evaluated at gp.
// @todo not checked if correct
{
    this->computeBmatrixAt(gp, answer);
}



void
Truss1d :: computeNmatrixAt(const FloatArray &iLocCoord, FloatMatrix &answer)
// Returns the displacement interpolation matrix {N} of the receiver,
// evaluated at gp.
{
    FloatArray n;
    this->interp.evalN( n, iLocCoord, FEIElementGeometryWrapper(this) );
    // Reshape
    answer.resize(1, 2); ///@todo It would be more suitable to follow the column-major order and just do answer.setColumn(...)
    answer.at(1, 1) = n.at(1);
    answer.at(1, 2) = n.at(2);
}


double
Truss1d :: computeVolumeAround(GaussPoint *gp)
// Returns the length of the receiver. This method is valid only if 1
// Gauss point is used.
{
    double detJ = fabs( this->interp.giveTransformationJacobian( gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) ) );
    return detJ *gp->giveWeight() * this->giveCrossSection()->give(CS_Area, gp);
}


void
Truss1d :: computeStressVector(FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep)
{
    this->giveStructuralCrossSection()->giveRealStress_1d(answer, gp, strain, tStep);
}

void
Truss1d :: computeConstitutiveMatrixAt(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
{
    this->giveStructuralCrossSection()->giveStiffnessMatrix_1d(answer, rMode, gp, tStep);
}

void
Truss1d :: giveDofManDofIDMask(int inode, IntArray &answer) const
{
    answer = {D_u};
}


void
Truss1d :: HuertaErrorEstimatorI_setupRefinedElementProblem(RefinedElement *refinedElement, int level, int nodeId,
                                                            IntArray &localNodeIdArray, IntArray &globalNodeIdArray,
                                                            HuertaErrorEstimatorInterface :: SetupMode sMode, TimeStep *tStep,
                                                            int &localNodeId, int &localElemId, int &localBcId,
                                                            IntArray &controlNode, IntArray &controlDof,
                                                            HuertaErrorEstimator :: AnalysisMode aMode)
{
    int inode, nodes = 2;
    FloatArray *corner [ 2 ], midNode, cor [ 2 ];
    double x = 0.0;

    if ( sMode == HuertaErrorEstimatorInterface :: NodeMode ||
        ( sMode == HuertaErrorEstimatorInterface :: BCMode && aMode == HuertaErrorEstimator :: HEE_linear ) ) {
        for ( inode = 0; inode < nodes; inode++ ) {
            corner [ inode ] = this->giveNode(inode + 1)->giveCoordinates();
            if ( corner [ inode ]->giveSize() != 3 ) {
                cor [ inode ].resize(3);
                cor [ inode ].at(1) = corner [ inode ]->at(1);
                cor [ inode ].at(2) = 0.0;
                cor [ inode ].at(3) = 0.0;

                corner [ inode ] = & ( cor [ inode ] );
            }

            x += corner [ inode ]->at(1);
        }

        midNode.resize(3);

        midNode.at(1) = x / nodes;
        midNode.at(2) = 0.0;
        midNode.at(3) = 0.0;
    }

    this->setupRefinedElementProblem1D(this, refinedElement, level, nodeId, localNodeIdArray, globalNodeIdArray,
                                       sMode, tStep, nodes, corner, midNode,
                                       localNodeId, localElemId, localBcId,
                                       controlNode, controlDof, aMode, "Truss1d");
}


void Truss1d :: HuertaErrorEstimatorI_computeNmatrixAt(GaussPoint *gp, FloatMatrix &answer)
{
    computeNmatrixAt(gp->giveSubPatchCoordinates(), answer);
}


#ifdef __OOFEG
void Truss1d :: drawRawGeometry(oofegGraphicContext &gc, TimeStep *tStep)
{
    GraphicObj *go;
    //  if (!go) { // create new one
    WCRec p [ 2 ]; /* point */
    if ( !gc.testElementGraphicActivity(this) ) {
        return;
    }

    EASValsSetLineWidth(OOFEG_RAW_GEOMETRY_WIDTH);
    EASValsSetColor( gc.getElementColor() );
    EASValsSetLayer(OOFEG_RAW_GEOMETRY_LAYER);
    p [ 0 ].x = ( FPNum ) this->giveNode(1)->giveCoordinate(1);
    p [ 0 ].y = 0.;
    p [ 0 ].z = 0.0;
    p [ 1 ].x = ( FPNum ) this->giveNode(2)->giveCoordinate(1);
    p [ 1 ].y = 0.;
    p [ 1 ].z = 0.0;
    go = CreateLine3D(p);
    EGWithMaskChangeAttributes(WIDTH_MASK | COLOR_MASK | LAYER_MASK, go);
    EGAttachObject(go, ( EObjectP ) this);
    EMAddGraphicsToModel(ESIModel(), go);
}


void Truss1d :: drawDeformedGeometry(oofegGraphicContext &gc, TimeStep *tStep, UnknownType type)
{
    GraphicObj *go;
    double defScale = gc.getDefScale();
    //  if (!go) { // create new one
    WCRec p [ 2 ]; /* point */
    if ( !gc.testElementGraphicActivity(this) ) {
        return;
    }

    EASValsSetLineWidth(OOFEG_DEFORMED_GEOMETRY_WIDTH);
    EASValsSetColor( gc.getDeformedElementColor() );
    EASValsSetLayer(OOFEG_DEFORMED_GEOMETRY_LAYER);
    p [ 0 ].x = ( FPNum ) this->giveNode(1)->giveUpdatedCoordinate(1, tStep, defScale);
    p [ 0 ].y = 0.;
    p [ 0 ].z = 0.;

    p [ 1 ].x = ( FPNum ) this->giveNode(2)->giveUpdatedCoordinate(1, tStep, defScale);
    p [ 1 ].y = 0.;
    p [ 1 ].z = 0.;
    go = CreateLine3D(p);
    EGWithMaskChangeAttributes(WIDTH_MASK | COLOR_MASK | LAYER_MASK, go);
    EMAddGraphicsToModel(ESIModel(), go);
}


void Truss1d :: drawScalar(oofegGraphicContext &gc, TimeStep *tStep)
{
    int i, indx, result = 0;
    WCRec p [ 2 ];
    GraphicObj *tr;
    FloatArray v1, v2;
    double s [ 2 ], defScale;

    if ( !gc.testElementGraphicActivity(this) ) {
        return;
    }

    if ( gc.giveIntVarMode() == ISM_recovered ) {
        result += this->giveInternalStateAtNode(v1, gc.giveIntVarType(), gc.giveIntVarMode(), 1, tStep);
        result += this->giveInternalStateAtNode(v2, gc.giveIntVarType(), gc.giveIntVarMode(), 2, tStep);
    } else if ( gc.giveIntVarMode() == ISM_local ) {
        GaussPoint *gp = integrationRulesArray [ 0 ]->getIntegrationPoint(0);
        result += giveIPValue(v1, gp, gc.giveIntVarType(), tStep);
        v2 = v1;
        result *= 2;
    }

    if ( result != 2 ) {
        return;
    }

    indx = gc.giveIntVarIndx();

    s [ 0 ] = v1.at(indx);
    s [ 1 ] = v2.at(indx);

    EASValsSetLayer(OOFEG_VARPLOT_PATTERN_LAYER);

    if ( ( gc.getScalarAlgo() == SA_ISO_SURF ) || ( gc.getScalarAlgo() == SA_ISO_LINE ) ) {
        for ( i = 0; i < 2; i++ ) {
            if ( gc.getInternalVarsDefGeoFlag() ) {
                // use deformed geometry
                defScale = gc.getDefScale();
                p [ i ].x = ( FPNum ) this->giveNode(i + 1)->giveUpdatedCoordinate(1, tStep, defScale);
                p [ i ].y = 0.;
                p [ i ].z = 0.;
            } else {
                p [ i ].x = ( FPNum ) this->giveNode(i + 1)->giveCoordinate(1);
                p [ i ].y = 0.;
                p [ i ].z = 0.;
            }
        }

        //EASValsSetColor(gc.getYieldPlotColor(ratio));
        tr =  CreateLine3D(p);
        EGWithMaskChangeAttributes(LAYER_MASK, tr);
        EMAddGraphicsToModel(ESIModel(), tr);
    } else if ( ( gc.getScalarAlgo() == SA_ZPROFILE ) || ( gc.getScalarAlgo() == SA_COLORZPROFILE ) ) {
        double landScale = gc.getLandScale();

        for ( i = 0; i < 2; i++ ) {
            if ( gc.getInternalVarsDefGeoFlag() ) {
                // use deformed geometry
                defScale = gc.getDefScale();
                p [ i ].x = ( FPNum ) this->giveNode(i + 1)->giveUpdatedCoordinate(1, tStep, defScale);
                p [ i ].y = 0.0;
                p [ i ].z = s [ i ] * landScale;
            } else {
                p [ i ].x = ( FPNum ) this->giveNode(i + 1)->giveCoordinate(1);
                p [ i ].y = 0.0;
                p [ i ].z = s [ i ] * landScale;
            }
        }

        if ( gc.getScalarAlgo() == SA_ZPROFILE ) {
            /*
             * EASValsSetColor(gc.getDeformedElementColor());
             * EASValsSetLineWidth(OOFEG_DEFORMED_GEOMETRY_WIDTH);
             * tr =  CreateLine3D(p);
             * EGWithMaskChangeAttributes(WIDTH_MASK | COLOR_MASK | LAYER_MASK, tr);
             */
            WCRec pp [ 4 ];
            pp [ 0 ].x = p [ 0 ].x;
            pp [ 0 ].y = 0.0;
            pp [ 0 ].z = 0.0;
            pp [ 1 ].x = p [ 0 ].x;
            pp [ 1 ].y = 0.0;
            pp [ 1 ].z = p [ 0 ].z;
            pp [ 2 ].x = p [ 1 ].x;
            pp [ 2 ].y = 0.0;
            pp [ 2 ].z = p [ 1 ].z;
            pp [ 3 ].x = p [ 1 ].x;
            pp [ 3 ].y = 0.0;
            pp [ 3 ].z = 0.0;
            tr = CreateQuad3D(pp);
            EASValsSetLineWidth(OOFEG_DEFORMED_GEOMETRY_WIDTH);
            EASValsSetColor( gc.getDeformedElementColor() );
            //EASValsSetLayer(OOFEG_DEFORMED_GEOMETRY_LAYER);
            EASValsSetFillStyle(FILL_HOLLOW);
            EGWithMaskChangeAttributes(WIDTH_MASK | FILL_MASK | COLOR_MASK | LAYER_MASK, tr);
            EMAddGraphicsToModel(ESIModel(), tr);
        } else {
            //tr =  CreateTriangleWD3D(p, s[0], s[1], s[2]);
            EASValsSetColor( gc.getDeformedElementColor() );
            tr =  CreateLine3D(p);
            EGWithMaskChangeAttributes(WIDTH_MASK | COLOR_MASK | LAYER_MASK, tr);
            EMAddGraphicsToModel(ESIModel(), tr);
        }
    }
}

#endif


Interface *
Truss1d :: giveInterface(InterfaceType interface)
{
    if ( interface == ZZNodalRecoveryModelInterfaceType ) {
        return static_cast< ZZNodalRecoveryModelInterface * >(this);
    } else if ( interface == NodalAveragingRecoveryModelInterfaceType ) {
        return static_cast< NodalAveragingRecoveryModelInterface * >(this);
    } else if ( interface == SpatialLocalizerInterfaceType ) {
        return static_cast< SpatialLocalizerInterface * >(this);
    } else if ( interface == ZZErrorEstimatorInterfaceType ) {
        return static_cast< ZZErrorEstimatorInterface * >(this);
    } else if ( interface == HuertaErrorEstimatorInterfaceType ) {
        return static_cast< HuertaErrorEstimatorInterface * >(this);
    }

    return NULL;
}


void
Truss1d :: NodalAveragingRecoveryMI_computeNodalValue(FloatArray &answer, int node,
                                                      InternalStateType type, TimeStep *tStep)
{
    GaussPoint *gp = integrationRulesArray [ 0 ]->getIntegrationPoint(0);
    this->giveIPValue(answer, gp, type, tStep);
}

} // end namespace oofem
