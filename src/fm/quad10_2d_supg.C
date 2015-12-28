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

#include "quad10_2d_supg.h"
#include "fei2dquadlin.h"
#include "fei2dquadconst.h"
#include "dofmanager.h"
#include "material.h"
#include "gausspoint.h"
#include "gaussintegrationrule.h"
#include "fluiddynamicmaterial.h"
#include "fluidcrosssection.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "intarray.h"
#include "domain.h"
#include "mathfem.h"
#include "engngm.h"
#include "timestep.h"
#include "materialinterface.h"
#include "contextioerr.h"
#include "crosssection.h"
#include "classfactory.h"

#ifdef __OOFEG
 #include "oofeggraphiccontext.h"
#endif

namespace oofem {
REGISTER_Element(Quad10_2D_SUPG);

FEI2dQuadLin Quad10_2D_SUPG :: velocityInterpolation(1, 2);
FEI2dQuadConst Quad10_2D_SUPG :: pressureInterpolation(1, 2);


Quad10_2D_SUPG :: Quad10_2D_SUPG(int n, Domain *aDomain) :
    SUPGElement2(n, aDomain), ZZNodalRecoveryModelInterface(this), pressureNode(1, aDomain, this)
{
    numberOfDofMans = 4;
}

Quad10_2D_SUPG :: ~Quad10_2D_SUPG()
{ }

FEInterpolation *
Quad10_2D_SUPG :: giveInterpolation() const
{
    return & this->velocityInterpolation;
}

FEInterpolation *
Quad10_2D_SUPG :: giveInterpolation(DofIDItem id) const
{
    if ( id == P_f ) {
        return & this->pressureInterpolation;
    } else {
        return & this->velocityInterpolation;
    }
}

DofManager *
Quad10_2D_SUPG :: giveInternalDofManager(int i) const
{
    return const_cast< ElementDofManager * >(& pressureNode);
}


int
Quad10_2D_SUPG :: computeNumberOfDofs()
{
    return 9;
}


void
Quad10_2D_SUPG :: giveDofManDofIDMask(int inode, IntArray &answer) const
{
    answer = {V_u, V_v};
}


void
Quad10_2D_SUPG :: giveInternalDofManDofIDMask(int i, IntArray &answer) const
{
    answer = {P_f};
}


IRResultType
Quad10_2D_SUPG :: initializeFrom(InputRecord *ir)
{
    this->pressureNode.initializeFrom(ir);

    return SUPGElement2 :: initializeFrom(ir);
}


void
Quad10_2D_SUPG :: giveInputRecord(DynamicInputRecord &input)
{
    SUPGElement2 :: giveInputRecord(input);
    this->pressureNode.giveInputRecord(input);
}


void
Quad10_2D_SUPG :: computeGaussPoints()
// Sets up the array containing the four Gauss points of the receiver.
{
    if ( integrationRulesArray.size() == 0 ) {
        integrationRulesArray.resize(3);


        integrationRulesArray [ 0 ].reset( new GaussIntegrationRule(1, this, 1, 3) );
        this->giveCrossSection()->setupIntegrationPoints(* integrationRulesArray [ 0 ], 4, this);

        //seven point Gauss integration
        integrationRulesArray [ 1 ].reset( new GaussIntegrationRule(2, this, 1, 3) );
        this->giveCrossSection()->setupIntegrationPoints(* integrationRulesArray [ 1 ], 4, this);

        integrationRulesArray [ 2 ].reset( new GaussIntegrationRule(3, this, 1, 3) );
        this->giveCrossSection()->setupIntegrationPoints(* integrationRulesArray [ 3 ], 4, this);
    }
}


void
Quad10_2D_SUPG :: computeNuMatrix(FloatMatrix &answer, GaussPoint *gp)
{
    FloatArray n;
    this->velocityInterpolation.evalN( n, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );
    answer.beNMatrixOf(n, 2);
 }


void
Quad10_2D_SUPG :: computeUDotGradUMatrix(FloatMatrix &answer, GaussPoint *gp, TimeStep *tStep)
{
    FloatMatrix n, dn;
    FloatArray u, un;
    this->velocityInterpolation.evaldNdx( dn, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );
    this->computeNuMatrix(n, gp);
    this->computeVectorOfVelocities(VM_Total, tStep, un);

    u.beProductOf(n, un);

    answer.resize(2, 8);
    answer.zero();
    for ( int i = 1; i <= 4; i++ ) {
        answer.at(1, 2 * i - 1) = dn.at(i, 1) * u.at(1) + dn.at(i, 2) * u.at(2);
        answer.at(2, 2 * i)     = dn.at(i, 1) * u.at(1) + dn.at(i, 2) * u.at(2);
    }
}

void
Quad10_2D_SUPG :: computeBMatrix(FloatMatrix &answer, GaussPoint *gp)
{
    FloatMatrix dn;
    this->velocityInterpolation.evaldNdx( dn, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );

    answer.resize(3, 8);
    answer.zero();

    for ( int i = 1; i <= 4; i++ ) {
        answer.at(1, 2 * i - 1) = dn.at(i, 1);
        answer.at(2, 2 * i)     = dn.at(i, 2);
        answer.at(3, 2 * i - 1) = dn.at(i, 2);
        answer.at(3, 2 * i)     = dn.at(i, 1);
    }
}

void
Quad10_2D_SUPG :: computeDivUMatrix(FloatMatrix &answer, GaussPoint *gp)
{
    FloatMatrix dn;
    velocityInterpolation.evaldNdx( dn, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );

    answer.resize(1, 8);
    answer.zero();

    for ( int i = 1; i <= 4; i++ ) {
        answer.at(1, 2 * i - 1) = dn.at(i, 1);
        answer.at(1, 2 * i) = dn.at(i, 2);
    }
}

void
Quad10_2D_SUPG :: computeNpMatrix(FloatMatrix &answer, GaussPoint *gp)
{
    FloatArray n;
    pressureInterpolation.evalN( n, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );

    answer.resize(1, 1);
    answer.zero();

    answer.at(1, 1) = n.at(1);
}


void
Quad10_2D_SUPG :: computeGradUMatrix(FloatMatrix &answer, GaussPoint *gp, TimeStep *tStep)
{
    FloatArray dnx(4), dny(4), u, u1(4), u2(4);
    FloatMatrix dn;

    answer.resize(2, 2);
    answer.zero();

    this->computeVectorOfVelocities(VM_Total, tStep, u);

    velocityInterpolation.evaldNdx( dn, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );
    for ( int i = 1; i <= 4; i++ ) {
        dnx.at(i) = dn.at(i, 1);
        dny.at(i) = dn.at(i, 2);

        u1.at(i) = u.at(2 * i - 1);
        u2.at(i) = u.at(2 * i);
    }


    answer.at(1, 1) =  u1.dotProduct(dnx);
    answer.at(1, 2) =  u1.dotProduct(dny);
    answer.at(2, 1) =  u2.dotProduct(dnx);
    answer.at(2, 2) =  u2.dotProduct(dny);
}

void
Quad10_2D_SUPG :: computeGradPMatrix(FloatMatrix &answer, GaussPoint *gp)
{
    FloatMatrix dn;
    pressureInterpolation.evaldNdx( dn, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );

    answer.beTranspositionOf(dn);
}


void
Quad10_2D_SUPG :: computeDivTauMatrix(FloatMatrix &answer, GaussPoint *gp, TimeStep *tStep)
{
    answer.resize(2, 8);
    answer.zero();
}


void
Quad10_2D_SUPG :: updateStabilizationCoeffs(TimeStep *tStep)
{
    double Re, norm_un, mu, mu_min, nu, norm_N, norm_N_d, norm_M_d, norm_LSIC, norm_G_c, norm_M_c, norm_N_c, t_p1, t_p2, t_p3, t_s1, t_s2, t_s3, rho;
    FloatMatrix N, N_d, M_d, LSIC, G_c, M_c, N_c;
    FloatArray u;

    mu_min = 1;
    rho = static_cast< FluidCrossSection * >( this->giveCrossSection() )->giveDensity( integrationRulesArray [ 0 ]->getIntegrationPoint(0) );
    for ( GaussPoint *gp: *integrationRulesArray [ 1 ] ) {
        mu = static_cast< FluidCrossSection * >( this->giveCrossSection() )->giveFluidMaterial()->giveEffectiveViscosity(gp, tStep);
        if ( mu_min > mu ) {
            mu_min = mu;
        }
    }

    nu = mu_min / rho;

    //this->computeVectorOfVelocities(VM_Total, tStep->givePreviousStep(), un);
    this->computeVectorOfVelocities(VM_Total, tStep, u);

    norm_un = u.computeNorm();

    this->computeAdvectionTerm(N, tStep);
    this->computeAdvectionDeltaTerm(N_d, tStep);
    this->computeMassDeltaTerm(M_d, tStep);
    this->computeLSICTerm(LSIC, tStep);
    this->computeLinearAdvectionTerm_MC(G_c, tStep);
    this->computeMassEpsilonTerm(M_c, tStep);
    this->computeAdvectionEpsilonTerm(N_c, tStep);


    norm_N = N.computeFrobeniusNorm();
    norm_N_d = N_d.computeFrobeniusNorm();
    norm_M_d = M_d.computeFrobeniusNorm();
    norm_LSIC = LSIC.computeFrobeniusNorm();
    norm_G_c = G_c.computeFrobeniusNorm();
    norm_M_c = M_c.computeFrobeniusNorm();
    norm_N_c = N_c.computeFrobeniusNorm();

    if ( ( norm_N == 0 ) || ( norm_N_d == 0 ) || ( norm_M_d == 0 ) ) {
        t_supg = 0;
    } else {
        Re = ( norm_un / nu ) * ( norm_N / norm_N_d );

        t_s1 = norm_N / norm_N_d;

        t_s2 = tStep->giveTimeIncrement() * ( norm_N / norm_M_d ) * 0.5;

        t_s3 = t_s1 * Re;

        t_supg =  1. / sqrt( 1. / ( t_s1 * t_s1 ) + 1. / ( t_s2 * t_s2 ) + 1. / ( t_s3 * t_s3 ) );
        // t_supg = 0;
    }

    if ( norm_LSIC == 0 ) {
        t_lsic = 0;
    } else {
        t_lsic = norm_N / norm_LSIC;

        // t_lsic = 0;
    }

    if ( ( norm_G_c == 0 ) || ( norm_N_c == 0 ) || ( norm_M_c == 0 ) ) {
        t_pspg = 0;
    } else {
        Re = ( norm_un / nu ) * ( norm_N / norm_N_d );

        t_p1 = norm_G_c / norm_N_c;

        t_p2 = tStep->giveTimeIncrement() * ( norm_G_c / norm_M_c ) * 0.5;

        t_p3 = t_p1 * Re;

        t_pspg =  1. / sqrt( 1. / ( t_p1 * t_p1 ) + 1. / ( t_p2 * t_p2 ) + 1. / ( t_p3 * t_p3 ) );
        // t_pspg = 0;
    }

    //t_pspg = 0;
}


void
Quad10_2D_SUPG :: computeAdvectionTerm(FloatMatrix &answer, TimeStep *tStep)
{
    FloatMatrix n, b;

    answer.clear();

    /* consistent part + supg stabilization term */
    for ( GaussPoint *gp: *this->integrationRulesArray [ 1 ] ) {
        this->computeNuMatrix(n, gp);
        this->computeUDotGradUMatrix(b, gp, tStep);
        double dV  = this->computeVolumeAround(gp);
        double rho = static_cast< FluidCrossSection * >( this->giveCrossSection() )->giveDensity(gp);
        answer.plusProductUnsym(n, b, rho * dV);
    }
}


void
Quad10_2D_SUPG :: computeAdvectionDeltaTerm(FloatMatrix &answer, TimeStep *tStep)
{
    FloatMatrix n, b;

    answer.clear();

    /* consistent part + supg stabilization term */
    for ( GaussPoint *gp: *this->integrationRulesArray [ 1 ] ) {
        this->computeNuMatrix(n, gp);
        this->computeUDotGradUMatrix(b, gp, tStep);
        double dV  = this->computeVolumeAround(gp);
        double rho = static_cast< FluidCrossSection * >( this->giveCrossSection() )->giveDensity(gp);

        answer.plusProductUnsym(b, b, rho * dV);
    }
}


void
Quad10_2D_SUPG :: computeMassDeltaTerm(FloatMatrix &answer, TimeStep *tStep)
{
    FloatMatrix n, b;

    answer.clear();

    /* mtrx for computing t_supg, norm of this mtrx is computed */
    for ( GaussPoint *gp: *this->integrationRulesArray [ 1 ] ) {
        this->computeNuMatrix(n, gp);
        this->computeUDotGradUMatrix(b, gp, tStep);
        double dV  = this->computeVolumeAround(gp);
        double rho = static_cast< FluidCrossSection * >( this->giveCrossSection() )->giveDensity(gp);

        answer.plusProductUnsym(b, n, rho * dV);
    }
}

void
Quad10_2D_SUPG :: computeLSICTerm(FloatMatrix &answer, TimeStep *tStep)
{
    FloatMatrix b;

    answer.clear();

    for ( GaussPoint *gp: *this->integrationRulesArray [ 1 ] ) {
        double dV  = this->computeVolumeAround(gp);
        double rho = static_cast< FluidCrossSection * >( this->giveCrossSection() )->giveDensity(gp);
        this->computeDivUMatrix(b, gp);

        answer.plusProductSymmUpper(b, b, dV * rho);
    }

    answer.symmetrized();
}


void
Quad10_2D_SUPG :: computeAdvectionEpsilonTerm(FloatMatrix &answer, TimeStep *tStep)
{
    //  to compute t_pspg
    FloatMatrix g, b;

    answer.clear();

    for ( GaussPoint *gp: *this->integrationRulesArray [ 1 ] ) {
        this->computeGradPMatrix(g, gp);
        this->computeUDotGradUMatrix(b, gp, tStep);
        double dV = this->computeVolumeAround(gp);

        answer.plusProductUnsym(g, b, dV);
    }
}

void
Quad10_2D_SUPG :: computeMassEpsilonTerm(FloatMatrix &answer, TimeStep *tStep)
{
    // to compute t_pspg
    FloatMatrix g, n;

    answer.clear();

    for ( GaussPoint *gp: *this->integrationRulesArray [ 1 ] ) {
        this->computeGradPMatrix(g, gp);
        this->computeNuMatrix(n, gp);
        double dV = this->computeVolumeAround(gp);

        answer.plusProductUnsym(g, n, dV);
    }
}

int
Quad10_2D_SUPG :: giveNumberOfSpatialDimensions()
{
    return 2;
}


double
Quad10_2D_SUPG :: LS_PCS_computeF(LevelSetPCS *ls, TimeStep *tStep)
{
#if 0
    FloatArray fi(3), un;

    this->computeVectorOfVelocities(VM_Total, tStep, un);
    for ( int i = 1; i <= 3; i++ ) {
        fi.at(i) = ls->giveLevelSetDofManValue( dofManArray.at(i) );
    }

    double fix = b [ 0 ] * fi.at(1) + b [ 1 ] * fi.at(2) + b [ 2 ] * fi.at(3);
    double fiy = c [ 0 ] * fi.at(1) + c [ 1 ] * fi.at(2) + c [ 2 ] * fi.at(3);
    double norm = sqrt(fix * fix + fiy * fiy);

    return ( 1. / 3. ) * ( fix * ( un.at(1) + un.at(3) + un.at(5) ) + fiy * ( un.at(2) + un.at(4) + un.at(6) ) ) / norm;

#endif
    return 0.0;
}


double
Quad10_2D_SUPG :: LS_PCS_computeS(LevelSetPCS *ls, TimeStep *tStep)
{
    return 0.0;
}


void
Quad10_2D_SUPG :: LS_PCS_computedN(FloatMatrix &answer)
{ }


void
Quad10_2D_SUPG :: LS_PCS_computeVOFFractions(FloatArray &answer, FloatArray &fi)
{ }


double
Quad10_2D_SUPG :: computeCriticalTimeStep(TimeStep *tStep)
{
    return 1.e6;
}


void
Quad10_2D_SUPG :: NodalAveragingRecoveryMI_computeNodalValue(FloatArray &answer, int node,
                                                             InternalStateType type, TimeStep *tStep)
{
    GaussPoint *gp = integrationRulesArray [ 0 ]->getIntegrationPoint(0);
    this->giveIPValue(answer, gp, type, tStep);
}


int
Quad10_2D_SUPG :: checkConsistency()
{
    return SUPGElement :: checkConsistency();
}


void
Quad10_2D_SUPG :: updateYourself(TimeStep *tStep)
{
    SUPGElement :: updateYourself(tStep);
}

int
Quad10_2D_SUPG :: giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep)
{
    if ( type == IST_VOFFraction ) {
        MaterialInterface *mi = domain->giveEngngModel()->giveMaterialInterface( domain->giveNumber() );
        if ( mi ) {
            FloatArray val;
            mi->giveElementMaterialMixture( val, gp->giveElement()->giveNumber() );
            answer = FloatArray{val.at(1)};
            return 1;
        } else {
            answer = FloatArray{1.0};
            return 1;
        }
    } else {
        return SUPGElement :: giveIPValue(answer, gp, type, tStep);
    }
}


contextIOResultType
Quad10_2D_SUPG :: saveContext(DataStream &stream, ContextMode mode, void *obj)
//
// saves full element context (saves state variables, that completely describe
// current state)
//
{
    contextIOResultType iores;

    if ( ( iores = SUPGElement :: saveContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    return CIO_OK;
}



contextIOResultType Quad10_2D_SUPG :: restoreContext(DataStream &stream, ContextMode mode, void *obj)
//
// restores full element context (saves state variables, that completely describe
// current state)
//
{
    contextIOResultType iores;

    if ( ( iores = SUPGElement :: restoreContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    return CIO_OK;
}


double
Quad10_2D_SUPG :: computeVolumeAround(GaussPoint *gp)
// Returns the portion of the receiver which is attached to gp.
{
    double determinant, weight, volume;

    determinant = fabs( this->velocityInterpolation.giveTransformationJacobian( gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) ) );


    weight = gp->giveWeight();
    volume = determinant * weight;

    return volume;
}

#if 0
double
Quad10_2D_SUPG :: computeVolumeAroundPressure(FEInterpolation2d &interpol, GaussPoint *gp)
// Returns the portion of the receiver which is attached to gp.
{
    double determinant, weight, volume;

    determinant = fabs( interpol.giveTransformationJacobian(domain, pressureDofManArray, gp->giveNaturalCoordinates(), 0.0) );

    weight = gp->giveWeight();
    volume = determinant * weight;

    return volume;
}
#endif

Interface *
Quad10_2D_SUPG :: giveInterface(InterfaceType interface)
{
    if ( interface == LevelSetPCSElementInterfaceType ) {
        return static_cast< LevelSetPCSElementInterface * >(this);
    } else if ( interface == ZZNodalRecoveryModelInterfaceType ) {
        return static_cast< ZZNodalRecoveryModelInterface * >(this);
    } else if ( interface == NodalAveragingRecoveryModelInterfaceType ) {
        return static_cast< NodalAveragingRecoveryModelInterface * >(this);
    }

    return NULL;
}


void
Quad10_2D_SUPG :: giveLocalVelocityDofMap(IntArray &map)
{
    map.enumerate(8);
}


void
Quad10_2D_SUPG :: giveLocalPressureDofMap(IntArray &map)
{
    map = {9};
}


#ifdef __OOFEG
int
Quad10_2D_SUPG :: giveInternalStateAtNode(FloatArray &answer, InternalStateType type, InternalStateMode mode,
                                          int node, TimeStep *tStep)
{
    return SUPGElement :: giveInternalStateAtNode(answer, type, mode, node, tStep);
}

void
Quad10_2D_SUPG :: drawRawGeometry(oofegGraphicContext &gc, TimeStep *tStep)
{
    WCRec p [ 3 ];
    GraphicObj *go;

    if ( !gc.testElementGraphicActivity(this) ) {
        return;
    }

    EASValsSetLineWidth(OOFEG_RAW_GEOMETRY_WIDTH);
    EASValsSetColor( gc.getElementColor() );
    EASValsSetEdgeColor( gc.getElementEdgeColor() );
    EASValsSetEdgeFlag(true);
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

void
Quad10_2D_SUPG :: drawScalar(oofegGraphicContext &gc, TimeStep *tStep)
{
    int i, indx, result = 0;
    WCRec p [ 3 ];
    GraphicObj *tr;
    FloatArray v1, v2, v3;
    double s [ 3 ];

    if ( !gc.testElementGraphicActivity(this) ) {
        return;
    }

    EASValsSetLayer(OOFEG_VARPLOT_PATTERN_LAYER);

    // if ((gc.giveIntVarMode() == ISM_local) && (gc.giveIntVarType() ==  IST_VOFFraction)) {
    if ( ( gc.giveIntVarType() ==  IST_VOFFraction ) && ( gc.giveIntVarMode() == ISM_local ) ) {
        Polygon matvolpoly;
        //this->formMaterialVolumePoly(matvolpoly, NULL, temp_normal, temp_p, false);
        EASValsSetColor( gc.getStandardSparseProfileColor() );
        //GraphicObj *go = matvolpoly.draw(gc,true,OOFEG_VARPLOT_PATTERN_LAYER);
        matvolpoly.draw(gc, true, OOFEG_VARPLOT_PATTERN_LAYER);
        return;
    }

    if ( gc.giveIntVarMode() == ISM_recovered ) {
        result += this->giveInternalStateAtNode(v1, gc.giveIntVarType(), gc.giveIntVarMode(), 1, tStep);
        result += this->giveInternalStateAtNode(v2, gc.giveIntVarType(), gc.giveIntVarMode(), 2, tStep);
        result += this->giveInternalStateAtNode(v3, gc.giveIntVarType(), gc.giveIntVarMode(), 3, tStep);
    } else if ( gc.giveIntVarMode() == ISM_local ) {
        GaussPoint *gp = integrationRulesArray [ 0 ]->getIntegrationPoint(0);
        result += giveIPValue(v1, gp, gc.giveIntVarType(), tStep);
        v2 = v1;
        v3 = v1;
        result *= 3;
    }

    if ( result != 3 ) {
        return;
    }

    indx = gc.giveIntVarIndx();

    s [ 0 ] = v1.at(indx);
    s [ 1 ] = v2.at(indx);
    s [ 2 ] = v3.at(indx);

    EASValsSetLayer(OOFEG_VARPLOT_PATTERN_LAYER);

    if ( gc.getScalarAlgo() == SA_ISO_SURF ) {
        for ( i = 0; i < 3; i++ ) {
            p [ i ].x = ( FPNum ) this->giveNode(i + 1)->giveCoordinate(1);
            p [ i ].y = ( FPNum ) this->giveNode(i + 1)->giveCoordinate(2);
            p [ i ].z = 0.;
        }

        //EASValsSetColor(gc.getYieldPlotColor(ratio));
        gc.updateFringeTableMinMax(s, 3);
        tr =  CreateTriangleWD3D(p, s [ 0 ], s [ 1 ], s [ 2 ]);
        EGWithMaskChangeAttributes(LAYER_MASK, tr);
        EMAddGraphicsToModel(ESIModel(), tr);
    } else if ( ( gc.getScalarAlgo() == SA_ZPROFILE ) || ( gc.getScalarAlgo() == SA_COLORZPROFILE ) ) {
        double landScale = gc.getLandScale();

        for ( i = 0; i < 3; i++ ) {
            p [ i ].x = ( FPNum ) this->giveNode(i + 1)->giveCoordinate(1);
            p [ i ].y = ( FPNum ) this->giveNode(i + 1)->giveCoordinate(2);
            p [ i ].z = s [ i ] * landScale;
        }

        if ( gc.getScalarAlgo() == SA_ZPROFILE ) {
            EASValsSetColor( gc.getDeformedElementColor() );
            EASValsSetLineWidth(OOFEG_DEFORMED_GEOMETRY_WIDTH);
            EASValsSetFillStyle(FILL_SOLID);
            tr =  CreateTriangle3D(p);
            EGWithMaskChangeAttributes(WIDTH_MASK | COLOR_MASK | FILL_MASK | LAYER_MASK, tr);
        } else {
            gc.updateFringeTableMinMax(s, 3);
            EASValsSetFillStyle(FILL_SOLID);
            tr =  CreateTriangleWD3D(p, s [ 0 ], s [ 1 ], s [ 2 ]);
            EGWithMaskChangeAttributes(FILL_MASK | LAYER_MASK, tr);
        }

        EMAddGraphicsToModel(ESIModel(), tr);
    }
}

#endif
} // end namespace oofem
