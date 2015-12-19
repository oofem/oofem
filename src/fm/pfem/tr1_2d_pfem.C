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

#include "tr1_2d_pfem.h"
#include "node.h"
#include "material.h"
#include "crosssection.h"
#include "gausspoint.h"
#include "gaussintegrationrule.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "intarray.h"
#include "domain.h"
#include "mathfem.h"
#include "engngm.h"
#include "pfem.h"
#include "fluiddynamicmaterial.h"
#include "fluidcrosssection.h"
#include "load.h"
#include "timestep.h"
#include "boundaryload.h"
#include "mathfem.h"
#include "contextioerr.h"

#ifdef __OOFEG
 #include "oofeggraphiccontext.h"
 #include "connectivitytable.h"
#endif

namespace oofem {
FEI2dTrLin TR1_2D_PFEM :: velocityInterpolation(1, 2);
FEI2dTrLin TR1_2D_PFEM :: pressureInterpolation(1, 2);

IntArray TR1_2D_PFEM :: edge_ordering [ 3 ] = {
    IntArray(6), IntArray(6), IntArray(6)
};

IntArray TR1_2D_PFEM :: velocityDofMask = {
    1, 2, 4, 5, 7, 8
};
IntArray TR1_2D_PFEM :: pressureDofMask = {
    3, 6, 9
};


TR1_2D_PFEM :: TR1_2D_PFEM(int n, Domain *aDomain, int particle1, int particle2, int particle3, int mat, int cs) :
    PFEMElement2d(n, aDomain)
{
    // Constructor.
    numberOfDofMans  = 3;
    IntArray aBodyLoadArry(1);
    aBodyLoadArry.at(1) = 3;
    this->setDofManagers({ particle1, particle2, particle3 });
    // CHECK THIS - NOT NICE
    this->setBodyLoads(aBodyLoadArry);
    this->material = mat;
    this->setCrossSection(cs);
    this->postInitialize();
}

TR1_2D_PFEM :: ~TR1_2D_PFEM()
// Destructor
{ }

int
TR1_2D_PFEM :: computeNumberOfDofs()
{
    return 9;
}

void
TR1_2D_PFEM ::   giveDofManDofIDMask(int inode, IntArray &answer) const
{
    answer = {
        V_u, V_v, P_f
    };
}

// NOT IN USE
void
TR1_2D_PFEM ::   giveElementDofIDMask(IntArray &answer) const
{
    this->giveDofManDofIDMask(1, answer);
}


IRResultType
TR1_2D_PFEM :: initializeFrom(InputRecord *ir)
{
    IRResultType ret = this->PFEMElement :: initializeFrom(ir);

    this->computeGaussPoints();
    return ret;
}

void
TR1_2D_PFEM :: computeGaussPoints()
// Sets up the array containing the four Gauss points of the receiver.
{
    if ( integrationRulesArray.size() == 0 ) {
        integrationRulesArray.resize(1);
        integrationRulesArray [ 0 ].reset( new GaussIntegrationRule(1, this, 1, 3) );
        integrationRulesArray [ 0 ]->setUpIntegrationPoints(_Triangle, 3, _2dFlow);
    }
}

void
TR1_2D_PFEM :: computeForceVector(FloatArray &answer, TimeStep *atTime) //F
{
    FloatArray vec;

    int nLoads = this->boundaryLoadArray.giveSize() / 2;
    answer.resize(6);
    answer.zero();
    for ( int i = 1; i <= nLoads; i++ ) {  // For each Neumann boundary condition
        int load_number = this->boundaryLoadArray.at(2 * i - 1);
        int load_id = this->boundaryLoadArray.at(2 * i);
        Load *load = this->domain->giveLoad(load_number);
        bcGeomType ltype = load->giveBCGeoType();

        if ( ltype == EdgeLoadBGT ) {
            this->computeEdgeBCSubVectorAt(vec, load, load_id, atTime);
            answer.add(vec);
        }
    }

    nLoads = this->giveBodyLoadArray()->giveSize();
    for ( int i = 1; i <= nLoads; i++ ) {
        Load *load = domain->giveLoad( bodyLoadArray.at(i) );
        bcGeomType ltype = load->giveBCGeoType();
        if ( ltype == BodyLoadBGT && load->giveBCValType() == ForceLoadBVT ) {
            this->computeBodyLoadVectorAt(vec, load, atTime);
            answer.add(vec);
        }
    }
}

//copied from tr21stokes
void
TR1_2D_PFEM :: computeBodyLoadVectorAt(FloatArray &answer, Load *load, TimeStep *atTime)
{
    IntegrationRule *iRule = this->integrationRulesArray [ 0 ].get();
    FloatArray N, gVector;

    answer.resize(6);
    answer.zero();

    load->computeComponentArrayAt(gVector, atTime, VM_Total);
    if ( gVector.giveSize() ) {
        for ( auto &gp : *iRule ) {
            double rho = static_cast< FluidCrossSection * >( this->giveCrossSection() )->giveDensity( integrationRulesArray [ 0 ]->getIntegrationPoint(0) );
            double detJ = this->pressureInterpolation.giveTransformationJacobian( gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );
            double dA = detJ * gp->giveWeight();
            this->pressureInterpolation.evalN( N, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );
            for ( int j = 0; j < 3; j++ ) {
                answer(2 * j)     += N(j) * rho * gVector(0) * dA;
                answer(2 * j + 1) += N(j) * rho * gVector(1) * dA;
            }
        }
    }
}


// NOT IN USE
//copied from tr21stokes
void
TR1_2D_PFEM :: computeEdgeBCSubVectorAt(FloatArray &answer, Load *load, int iEdge, TimeStep *atTime)
{
    answer.resize(15);
    answer.zero();

    if ( load->giveType() == TransmissionBC ) { // Neumann boundary conditions (traction)
        BoundaryLoad *boundaryLoad = ( BoundaryLoad * ) load;

        int numberOfEdgeIPs = ( int ) ceil( ( boundaryLoad->giveApproxOrder() + 1. ) / 2. ) * 2;

        GaussIntegrationRule iRule(1, this, 1, 1);
        FloatArray N, t, f(6), f2(6);
        IntArray edge_mapping;

        f.zero();
        iRule.setUpIntegrationPoints(_Line, numberOfEdgeIPs, _Unknown);

        for ( auto &gp : iRule ) {
            const FloatArray &lcoords = gp->giveNaturalCoordinates();
            this->pressureInterpolation.edgeEvalN( N, iEdge, lcoords, FEIElementGeometryWrapper(this) );
            double dS = gp->giveWeight() * this->pressureInterpolation.edgeGiveTransformationJacobian( iEdge, lcoords, FEIElementGeometryWrapper(this) );

            if ( boundaryLoad->giveFormulationType() == Load :: FT_Entity ) {         // Edge load in xi-eta system
                boundaryLoad->computeValueAt(t, atTime, lcoords, VM_Total);
            } else {   // Edge load in x-y system
                FloatArray gcoords;
                this->pressureInterpolation.edgeLocal2global( gcoords, iEdge, lcoords, FEIElementGeometryWrapper(this) );
                boundaryLoad->computeValueAt(t, atTime, gcoords, VM_Total);
            }

            // Reshape the vector
            for ( int j = 0; j < 3; j++ ) {
                f(2 * j)   += N(j) * t(0) * dS;
                f(2 * j + 1) += N(j) * t(1) * dS;
            }
        }

        answer.assemble(f, this->edge_ordering [ iEdge - 1 ]);
    } else {
        OOFEM_ERROR("Strange boundary condition type");
    }
}

void
TR1_2D_PFEM :: computeDiagonalMassMtrx(FloatArray &answer, TimeStep *atTime)
{
    answer.resize(6);
    answer.zero();

    double rho = static_cast< FluidCrossSection * >( this->giveCrossSection() )->giveDensity( integrationRulesArray [ 0 ]->getIntegrationPoint(0) );
    double mm = rho * this->area / 3.0;
    for ( int i = 1; i <= 6; i++ ) {
        answer.at(i) = mm;
    }
}

void
TR1_2D_PFEM :: computeDiagonalMassMtrx(FloatMatrix &answer, TimeStep *atTime)
{
    answer.resize(6, 6);
    answer.zero();
    double rho = static_cast< FluidCrossSection * >( this->giveCrossSection() )->giveDensity( integrationRulesArray [ 0 ]->getIntegrationPoint(0) );
    double mm = rho * this->area / 3.0;
    for ( int i = 1; i <= 6; i++ ) {
        answer.at(i, i) = mm;
    }
}

Interface *
TR1_2D_PFEM :: giveInterface(InterfaceType /*interface*/)
{
    return NULL;
}

void
TR1_2D_PFEM :: computeDeviatoricStress(FloatArray &answer, GaussPoint *gp, TimeStep *tStep)
{
    /* one should call material driver instead */
    FloatArray u(6), eps(3);
    answer.resize(3);


    this->computeVectorOf(VM_Total, tStep, u);

    eps.at(1) = ( b [ 0 ] * u.at(1) + b [ 1 ] * u.at(3) + b [ 2 ] * u.at(5) );
    eps.at(2) = ( c [ 0 ] * u.at(2) + c [ 1 ] * u.at(4) + c [ 2 ] * u.at(6) );
    eps.at(3) = ( b [ 0 ] * u.at(2) + b [ 1 ] * u.at(4) + b [ 2 ] * u.at(6) + c [ 0 ] * u.at(1) + c [ 1 ] * u.at(3) + c [ 2 ] * u.at(5) );
    FluidDynamicMaterial *mat = static_cast< FluidCrossSection * >( this->giveCrossSection() )->giveFluidMaterial();
    mat->computeDeviatoricStressVector(answer, gp, eps, tStep);
}

void
TR1_2D_PFEM :: computeDeviatoricStressDivergence(FloatArray &answer, TimeStep *atTime)
{
    answer.resize(6);
    answer.zero();
    FloatArray stress;

    this->computeDeviatoricStress(stress, integrationRulesArray [ 0 ]->getIntegrationPoint(0), atTime);

    // \int dNu/dxj \Tau_ij
    for ( int i = 0; i < 3; i++ ) {
        answer.at( ( i ) * 2 + 1 ) = area * ( stress.at(1) * b [ i ] + stress.at(3) * c [ i ] );
        answer.at( ( i + 1 ) * 2 ) = area * ( stress.at(3) * b [ i ] + stress.at(2) * c [ i ] );
    }
}

int
TR1_2D_PFEM :: checkConsistency()
{
    Node *node1, *node2, *node3;
    double x1, x2, x3, y1, y2, y3;

    node1 = giveNode(1);
    node2 = giveNode(2);
    node3 = giveNode(3);

    // init geometry data
    x1 = node1->giveCoordinate(1);
    x2 = node2->giveCoordinate(1);
    x3 = node3->giveCoordinate(1);

    y1 = node1->giveCoordinate(2);
    y2 = node2->giveCoordinate(2);
    y3 = node3->giveCoordinate(2);

    this->area = 0.5 * ( x2 * y3 + x1 * y2 + y1 * x3 - x2 * y1 - x3 * y2 - x1 * y3 );

    b [ 0 ] = ( y2 - y3 ) / ( 2. * area );
    c [ 0 ] = ( x3 - x2 ) / ( 2. * area );
    b [ 1 ] = ( y3 - y1 ) / ( 2. * area );
    c [ 1 ] = ( x1 - x3 ) / ( 2. * area );
    b [ 2 ] = ( y1 - y2 ) / ( 2. * area );
    c [ 2 ] = ( x2 - x1 ) / ( 2. * area );

    return PFEMElement2d :: checkConsistency();
}

double
TR1_2D_PFEM :: computeCriticalTimeStep(TimeStep *tStep)
{
    double deltaT = domain->giveEngngModel()->giveSolutionStepWhenIcApply()->giveTimeIncrement();

    // assuming plane x-y !!!!!!!!!!!!!!!!!!!!!!

    FloatArray u;
    this->computeVectorOf(VM_Total, tStep, u);
    FloatArray u1(2), u2(2), u3(2);
    u1.at(1) = u.at(1);
    u1.at(2) = u.at(2);

    u2.at(1) = u.at(4);
    u2.at(2) = u.at(5);

    u3.at(1) = u.at(7);
    u3.at(2) = u.at(8);

    Node *node1 = giveNode(1);
    Node *node2 = giveNode(2);
    Node *node3 = giveNode(3);

    double x1, x2, x3, y1, y2, y3;

    x1 = node1->giveCoordinate(1);
    x2 = node2->giveCoordinate(1);
    x3 = node3->giveCoordinate(1);

    y1 = node1->giveCoordinate(2);
    y2 = node2->giveCoordinate(2);
    y3 = node3->giveCoordinate(2);

    FloatArray p1(2), p2(2), p3(2);
    p1.at(1) = x1;
    p1.at(2) = y1;

    p2.at(1) = x2;
    p2.at(2) = y2;

    p3.at(1) = x3;
    p3.at(2) = y3;

    FloatArray s12(p2), s23(p3), s31(p1);

    s12.subtract(p1);
    FloatArray s21(s12);
    s21.negated();

    s23.subtract(p2);
    FloatArray s32(s23);
    s32.negated();

    s31.subtract(p3);
    FloatArray s13(s31);
    s13.negated();

    FloatArray foot1(p2);
    FloatArray foot2(p3);
    FloatArray foot3(p1);

    double factor = s21.dotProduct(s23) / s23.dotProduct(s23);
    foot1.add(factor, s23);

    FloatArray altitude1(foot1);
    altitude1.subtract(p1);

    factor = s32.dotProduct(s31) / s31.dotProduct(s31);
    foot2.add(factor, s31);

    FloatArray altitude2(foot2);
    altitude2.subtract(p2);

    factor = s13.dotProduct(s12) / s12.dotProduct(s12);
    foot3.add(factor, s12);

    FloatArray altitude3(foot3);
    altitude3.subtract(p3);

    double dt1(deltaT), dt2(deltaT), dt3(deltaT);

    double u1_proj = u1.dotProduct(altitude1) / altitude1.computeNorm();
    if ( u1_proj > 1.e-6 ) {
        dt1 = altitude1.computeNorm() / u1_proj;
    }

    double u2_proj = u2.dotProduct(altitude2) / altitude2.computeNorm();
    if ( u2_proj > 1.e-6 ) {
        dt2 = altitude1.computeNorm() / u2_proj;
    }

    double u3_proj = u3.dotProduct(altitude3) / altitude3.computeNorm();
    if ( u3_proj > 1.e-6 ) {
        dt3 = altitude3.computeNorm() / u3_proj;
    }

    double dt_min = min( dt1, min(dt2, dt3) );

    return dt_min;
}




void
TR1_2D_PFEM :: updateYourself(TimeStep *tStep)
{
    PFEMElement :: updateYourself(tStep);
}


void
TR1_2D_PFEM :: printOutputAt(FILE *file, TimeStep *stepN)
// Performs end-of-step operations.
{
    PFEMElement :: printOutputAt(file, stepN);
}



contextIOResultType TR1_2D_PFEM :: saveContext(DataStream *stream, ContextMode mode, void *obj)
//
// saves full element context (saves state variables, that completely describe
// current state)
//
{
    contextIOResultType iores;

    if ( ( iores = PFEMElement :: saveContext(* stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    return CIO_OK;
}



contextIOResultType TR1_2D_PFEM :: restoreContext(DataStream *stream, ContextMode mode, void *obj)
//
// restores full element context (saves state variables, that completely describe
// current state)
//
{
    contextIOResultType iores;

    if ( ( iores = PFEMElement :: restoreContext(* stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }


    return CIO_OK;
}

double
TR1_2D_PFEM :: computeVolumeAround(GaussPoint *gp)
// Returns the portion of the receiver which is attached to gp.
{
    double determinant = fabs( this->velocityInterpolation.giveTransformationJacobian( gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) ) );
    return determinant * gp->giveWeight();
}

#ifdef __OOFEG
int
TR1_2D_PFEM :: giveInternalStateAtNode(FloatArray &answer, InternalStateType type, InternalStateMode mode,
                                       int node, TimeStep *atTime)
{
    //<RESTRICTED_SECTION>
    if ( type == IST_VOFFraction ) {
        answer.resize(1);
        //        answer.at(1) = this->giveTempVolumeFraction();
        return 1;
    } else
    //</RESTRICTED_SECTION>
    if ( type == IST_Density ) {
        answer.resize(1);
        answer.at(1) = static_cast< FluidCrossSection * >( this->giveCrossSection() )->giveDensity( integrationRulesArray [ 0 ]->getIntegrationPoint(0) );
        return 1;
    } else {
        return PFEMElement :: giveInternalStateAtNode(answer, type, mode, node, atTime);
    }
}

void
TR1_2D_PFEM :: drawRawGeometry(oofegGraphicContext &gc)
{
    WCRec p [ 3 ];
    GraphicObj *go;

    if ( !gc.testElementGraphicActivity(this) ) {
        return;
    }

    EASValsSetLineWidth(OOFEG_RAW_GEOMETRY_WIDTH);
    EASValsSetColor( gc.getElementColor() );
    EASValsSetEdgeColor( gc.getElementEdgeColor() );
    EASValsSetEdgeFlag(TRUE);
    EASValsSetLayer(OOFEG_RAW_GEOMETRY_LAYER);
    p [ 0 ].x = ( FPNum ) this->giveNode(1)->giveCoordinate(1);
    p [ 0 ].y = ( FPNum ) this->giveNode(1)->giveCoordinate(2);
    p [ 0 ].z = 0.;
    p [ 1 ].x = ( FPNum ) this->giveNode(2)->giveCoordinate(1);
    p [ 1 ].y = ( FPNum ) this->giveNode(2)->giveCoordinate(2);
    p [ 1 ].z = 0.;
    p [ 2 ].x = ( FPNum ) this->giveNode(3)->giveCoordinate(1);
    p [ 2 ].y = ( FPNum ) this->giveNode(3)->giveCoordinate(2);
    p [ 2 ].z = 0.;

    go =  CreateTriangle3D(p);
    EGWithMaskChangeAttributes(WIDTH_MASK | COLOR_MASK | EDGE_COLOR_MASK | EDGE_FLAG_MASK | LAYER_MASK, go);
    EGAttachObject(go, ( EObjectP ) this);
    EMAddGraphicsToModel(ESIModel(), go);
}

void TR1_2D_PFEM :: drawScalar(oofegGraphicContext &context)
{
    int i, indx, result = 0;
    WCRec p [ 3 ];
    GraphicObj *tr;
    TimeStep *tStep = this->giveDomain()->giveEngngModel()->giveCurrentStep();
    FloatArray v1, v2, v3;
    double s [ 3 ];
    IntArray map;

    if ( !context.testElementGraphicActivity(this) ) {
        return;
    }

    if ( context.giveIntVarMode() == ISM_recovered ) {
        result += this->giveInternalStateAtNode(v1, context.giveIntVarType(), context.giveIntVarMode(), 1, tStep);
        result += this->giveInternalStateAtNode(v2, context.giveIntVarType(), context.giveIntVarMode(), 2, tStep);
        result += this->giveInternalStateAtNode(v3, context.giveIntVarType(), context.giveIntVarMode(), 3, tStep);
    } else if ( context.giveIntVarMode() == ISM_local ) {
        GaussPoint *gp = integrationRulesArray [ 0 ]->getIntegrationPoint(0);
        result += giveIPValue(v1, gp, context.giveIntVarType(), tStep);
        v2 = v1;
        v3 = v1;
        result *= 3;
    }

    if ( result != 3 ) {
        return;
    }

    indx = context.giveIntVarIndx();

    s [ 0 ] = v1.at(indx);
    s [ 1 ] = v2.at(indx);
    s [ 2 ] = v3.at(indx);

    EASValsSetLayer(OOFEG_VARPLOT_PATTERN_LAYER);

    if ( context.getScalarAlgo() == SA_ISO_SURF ) {
        for ( i = 0; i < 3; i++ ) {
            p [ i ].x = ( FPNum ) this->giveNode(i + 1)->giveCoordinate(1);
            p [ i ].y = ( FPNum ) this->giveNode(i + 1)->giveCoordinate(2);
            p [ i ].z = 0.;
        }

        //EASValsSetColor(gc.getYieldPlotColor(ratio));
        context.updateFringeTableMinMax(s, 3);
        tr =  CreateTriangleWD3D(p, s [ 0 ], s [ 1 ], s [ 2 ]);
        EGWithMaskChangeAttributes(LAYER_MASK, tr);
        EMAddGraphicsToModel(ESIModel(), tr);
    } else if ( ( context.getScalarAlgo() == SA_ZPROFILE ) || ( context.getScalarAlgo() == SA_COLORZPROFILE ) ) {
        double landScale = context.getLandScale();

        for ( i = 0; i < 3; i++ ) {
            p [ i ].x = ( FPNum ) this->giveNode(i + 1)->giveCoordinate(1);
            p [ i ].y = ( FPNum ) this->giveNode(i + 1)->giveCoordinate(2);
            p [ i ].z = s [ i ] * landScale;
        }

        if ( context.getScalarAlgo() == SA_ZPROFILE ) {
            EASValsSetColor( context.getDeformedElementColor() );
            EASValsSetLineWidth(OOFEG_DEFORMED_GEOMETRY_WIDTH);
            EASValsSetFillStyle(FILL_SOLID);
            tr =  CreateTriangle3D(p);
            EGWithMaskChangeAttributes(WIDTH_MASK | COLOR_MASK | FILL_MASK | LAYER_MASK, tr);
        } else {
            context.updateFringeTableMinMax(s, 3);
            EASValsSetFillStyle(FILL_SOLID);
            tr =  CreateTriangleWD3D(p, s [ 0 ], s [ 1 ], s [ 2 ]);
            EGWithMaskChangeAttributes(FILL_MASK | LAYER_MASK, tr);
        }

        EMAddGraphicsToModel(ESIModel(), tr);
    }
}



#endif
} // end namespace oofem
