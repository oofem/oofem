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

#include "igaelements.h"
#include "flotarry.h"
#include "flotmtrx.h"
#include "domain.h"
#include "node.h"
#include "element.h"
#include "gausspnt.h"
#include "gaussintegrationrule.h"
#include "matresponsemode.h"
#include "crosssection.h"
#include "structuralcrosssection.h"
#include "mathfem.h"
#include "iga.h"

#ifdef __OOFEG
 #include "oofeggraphiccontext.h"
 #include "oofegutils.h"
#endif

namespace oofem {

BsplinePlaneStressElement :: BsplinePlaneStressElement(int n, Domain *aDomain) : IGAElement(n, aDomain), PlaneStressStructuralElementEvaluator(), interpolation(2) { }


IRResultType BsplinePlaneStressElement :: initializeFrom(InputRecord *ir) {
    BSplineInterpolation *interpol = ( BSplineInterpolation * ) this->giveInterpolation();
    IGAElement :: initializeFrom(ir);
    //PlaneStressStructuralElementEvaluator::initializeFrom(ir);

    // HUHU checkConsistency()
    if ( giveNumberOfDofManagers() != interpol->giveNumberOfControlPoints(1) * interpol->giveNumberOfControlPoints(2) ) {
        OOFEM_ERROR("BsplinePlaneStressElement::initializeFrom - number of control points mismatch");
    }

    return IRRT_OK;
}



NURBSPlaneStressElement :: NURBSPlaneStressElement(int n, Domain *aDomain) : IGAElement(n, aDomain), PlaneStressStructuralElementEvaluator(), interpolation(2) { }


IRResultType NURBSPlaneStressElement :: initializeFrom(InputRecord *ir) {
    NURBSInterpolation *interpol = ( NURBSInterpolation * ) this->giveInterpolation();
    IGAElement :: initializeFrom(ir);
    //PlaneStressStructuralElementEvaluator::initializeFrom(ir);

    // HUHU
    if ( giveNumberOfDofManagers() != interpol->giveNumberOfControlPoints(1) * interpol->giveNumberOfControlPoints(2) ) {
        OOFEM_ERROR("NURBSPlaneStressElement::initializeFrom - number of control points mismatch");
    }

    return IRRT_OK;
}


TSplinePlaneStressElement :: TSplinePlaneStressElement(int n, Domain *aDomain) : IGATSplineElement(n, aDomain), PlaneStressStructuralElementEvaluator(), interpolation(2) { }




NURBSSpace3dElement :: NURBSSpace3dElement(int n, Domain *aDomain) : IGAElement(n, aDomain), Space3dStructuralElementEvaluator(), interpolation(3) { }


IRResultType NURBSSpace3dElement :: initializeFrom(InputRecord *ir) {
    NURBSInterpolation *interpol = ( NURBSInterpolation * ) this->giveInterpolation();
    IGAElement :: initializeFrom(ir);
    //PlaneStressStructuralElementEvaluator::initializeFrom(ir);

    // HUHU
    if ( giveNumberOfDofManagers() != interpol->giveNumberOfControlPoints(1) * interpol->giveNumberOfControlPoints(2) * interpol->giveNumberOfControlPoints(3) ) {
        OOFEM_ERROR("NURBSSpace3dElement::initializeFrom - number of control points mismatch");
    }

    return IRRT_OK;
}

// HUHU should be implemented by IGA element (it is the same for Bspline NURBS and TSpline)
// however in such a case it should be generalized in terms of appropriately multiplying
// nseq for those integration elements which span more tham just a single knot span
// the reason is to ensure compatible division to quads over which scalar quantity is interpolated
// bilinearly !!!

#ifdef __OOFEG

//#define COMPUTE_STRESS
#define COMPUTE_STRAIN

void BsplinePlaneStressElement :: drawScalar(oofegGraphicContext &context)
{
    int indx;
    WCRec p [ 4 ];
    GraphicObj *go;
    TimeStep *tStep = this->giveDomain()->giveEngngModel()->giveCurrentStep();
    double s [ 4 ];
    IntArray map;
    FEInterpolation *interp = this->giveInterpolation();
    FloatArray val;

    if ( !context.testElementGraphicActivity(this) ) {
        return;
    }

    EASValsSetLayer(OOFEG_VARPLOT_PATTERN_LAYER);
    EASValsSetFillStyle(FILL_SOLID);
    numberOfIntegrationRules = this->giveNumberOfIntegrationRules();
    const double * const * knotVector = interp->giveKnotVector();
    const IntArray *span;
    IntegrationRule *iRule;
    int ir, j, nsd = this->giveNsd();
    FloatArray c [ 4 ], cg [ 4 ];
    IntArray sign[4];

    if ( nsd == 2 ) {
        for ( j = 0; j < 4; j++ ) {
            c [ j ].resize(2);
            cg [ j ].resize(2);
            sign[ j ].resize(2);
        }
        sign [ 0 ].at(1) = 1;
        sign [ 0 ].at(2) = 1;
        sign [ 1 ].at(1) = -1;
        sign [ 1 ].at(2) = 1;
        sign [ 2 ].at(1) = -1;
        sign [ 2 ].at(2) = -1;
        sign [ 3 ].at(1) = 1;
        sign [ 3 ].at(2) = -1;
    } else {
        OOFEM_ERROR2("drawRawGeometry: not implemented for nsd = %d", nsd);
    }

    this->giveIntVarCompFullIndx( map, context.giveIntVarType() );
    if ( ( indx = map.at( context.giveIntVarIndx() ) ) == 0 ) {
        return;
    }

    // loop over individual integration rules (i.e., knot spans)
    for ( ir = 0; ir < numberOfIntegrationRules; ir++ ) {
        iRule = this->giveIntegrationRule(ir);
        span = iRule->giveKnotSpan();
        if ( nsd == 2 ) {
            // divide span locally to get finer geometry rep.
            int i, j, k, nseg = 4;
            double du = ( knotVector [ 0 ] [ span->at(1) + 1 ] - knotVector [ 0 ] [ span->at(1) ] ) / nseg;
            double dv = ( knotVector [ 1 ] [ span->at(2) + 1 ] - knotVector [ 1 ] [ span->at(2) ] ) / nseg;
            for ( i = 1; i <= nseg; i++ ) {
                for ( j = 1; j <= nseg; j++ ) {
                    c [ 0 ].at(1) = knotVector [ 0 ] [ span->at(1) ] + du * ( i - 1 );
                    c [ 0 ].at(2) = knotVector [ 1 ] [ span->at(2) ] + dv * ( j - 1 );
                    c [ 1 ].at(1) = knotVector [ 0 ] [ span->at(1) ] + du * i;
                    c [ 1 ].at(2) = knotVector [ 1 ] [ span->at(2) ] + dv * ( j - 1 );
                    c [ 2 ].at(1) = knotVector [ 0 ] [ span->at(1) ] + du * i;
                    c [ 2 ].at(2) = knotVector [ 1 ] [ span->at(2) ] + dv * j;
                    c [ 3 ].at(1) = knotVector [ 0 ] [ span->at(1) ] + du * ( i - 1 );
                    c [ 3 ].at(2) = knotVector [ 1 ] [ span->at(2) ] + dv * j;

                    for ( k = 0; k < 4; k++ ) {
                        interp->local2global(cg [ k ], c [ k ], FEIIGAElementGeometryWrapper( this, iRule->giveKnotSpan() ));
                        p [ k ].x = ( FPNum ) cg [ k ].at(1);
                        p [ k ].y = ( FPNum ) cg [ k ].at(2);
                        p [ k ].z = 0.;

#ifdef QUARTER_PLATE_WITH_HOLE_SINGLE_PATCH
                        //move sampling gp out of boundary to overcome degeneracy on quarter plate with hole modelled by single patch
                        if(c[k].at(1) > 0.99999 && c[k].at(2) > 0.495 && c[k].at(2) < 0.505){
                            c[k].at(1) += sign[k].at(1)*du/10.0;
                            c[k].at(2) += sign[k].at(2)*dv/10.0;
                            c[k].at(3) += sign[k].at(3)*dw/10.0;
                        }
#endif

                        // create a dummy ip's
                        FloatArray *cc = new FloatArray (c[k]);         // constructor of gp does not make its own copy
                        GaussPoint gp(iRule, 999, cc, 1.0, _PlaneStress);
#ifdef COMPUTE_STRAIN
                        this->computeStrainVector (val, &gp, tStep);
#endif
#ifdef COMPUTE_STRESS
                        this->computeStressVector(val, & gp, tStep);
#endif
                        s [ k ] = val.at(indx);
                    }

                    if ( ( isnan(s [ 0 ]) ) || ( isnan(s [ 1 ]) ) || ( isnan(s [ 2 ]) ) || ( isnan(s [ 3 ]) ) ) {
                        continue;
                    }

                    if ( ( fabs(s [ 0 ]) > 1.e5 ) || ( fabs(s [ 1 ]) > 1.e5 ) || ( fabs(s [ 2 ]) > 1.e5 ) || ( fabs(s [ 3 ]) > 1.e5 ) ) {
                        continue;
                    }

                    //printf ("QWD: %e %e %e %e\n", s [ 0 ], s [ 1 ], s [ 2 ], s [ 3 ]);

                    go =  CreateQuadWD3D(p, s [ 0 ], s [ 1 ], s [ 2 ], s [ 3 ]);
                    EGWithMaskChangeAttributes(WIDTH_MASK | FILL_MASK | COLOR_MASK | EDGE_COLOR_MASK | EDGE_FLAG_MASK | LAYER_MASK, go);
                    EGAttachObject(go, ( EObjectP ) this);
                    EMAddGraphicsToModel(ESIModel(), go);
                }
            }
        }
    } // end loop over knot spans (irules)

}


void NURBSPlaneStressElement :: drawScalar(oofegGraphicContext &context) {
    int indx;
    WCRec p [ 4 ];
    GraphicObj *go;
    TimeStep *tStep = this->giveDomain()->giveEngngModel()->giveCurrentStep();
    double s [ 4 ];
    IntArray map;
    FEInterpolation *interp = this->giveInterpolation();
    FloatArray val;

    if ( !context.testElementGraphicActivity(this) ) {
        return;
    }

    EASValsSetLayer(OOFEG_VARPLOT_PATTERN_LAYER);
    EASValsSetFillStyle(FILL_SOLID);
    EASValsSetEdgeFlag(true);
    numberOfIntegrationRules = this->giveNumberOfIntegrationRules();
    const double * const * knotVector = interp->giveKnotVector();
    const IntArray *span;
    IntegrationRule *iRule;
    int ir, j, nsd = this->giveNsd();
    FloatArray c [ 4 ], cg [ 4 ];
    IntArray sign [ 4 ];

    if ( nsd == 2 ) {
        for ( j = 0; j < 4; j++ ) {
            c [ j ].resize(2);
            cg [ j ].resize(2);
            sign [ j ].resize(2);
        }
        sign [ 0 ].at(1) = 1;
        sign [ 0 ].at(2) = 1;
        sign [ 1 ].at(1) = -1;
        sign [ 1 ].at(2) = 1;
        sign [ 2 ].at(1) = -1;
        sign [ 2 ].at(2) = -1;
        sign [ 3 ].at(1) = 1;
        sign [ 3 ].at(2) = -1;
    } else {
        OOFEM_ERROR2("drawRawGeometry: not implemented for nsd = %d", nsd);
    }

    this->giveIntVarCompFullIndx( map, context.giveIntVarType() );
    if ( ( indx = map.at( context.giveIntVarIndx() ) ) == 0 ) {
        return;
    }

    //double maxs=-1.0e10, mins=1.0e10;

    // loop over individual integration rules (i.e., knot spans)
    for ( ir = 0; ir < numberOfIntegrationRules; ir++ ) {
        iRule = this->giveIntegrationRule(ir);
        span = iRule->giveKnotSpan();
        if ( nsd == 2 ) {
            // divide span locally to get finer geometry rep.
            int i, j, k, nseg = 8;
            double du = ( knotVector [ 0 ] [ span->at(1) + 1 ] - knotVector [ 0 ] [ span->at(1) ] ) / nseg;
            double dv = ( knotVector [ 1 ] [ span->at(2) + 1 ] - knotVector [ 1 ] [ span->at(2) ] ) / nseg;
            for ( i = 1; i <= nseg; i++ ) {
                for ( j = 1; j <= nseg; j++ ) {
                    c [ 0 ].at(1) = knotVector [ 0 ] [ span->at(1) ] + du * ( i - 1 );
                    c [ 0 ].at(2) = knotVector [ 1 ] [ span->at(2) ] + dv * ( j - 1 );
                    c [ 1 ].at(1) = knotVector [ 0 ] [ span->at(1) ] + du * i;
                    c [ 1 ].at(2) = knotVector [ 1 ] [ span->at(2) ] + dv * ( j - 1 );
                    c [ 2 ].at(1) = knotVector [ 0 ] [ span->at(1) ] + du * i;
                    c [ 2 ].at(2) = knotVector [ 1 ] [ span->at(2) ] + dv * j;
                    c [ 3 ].at(1) = knotVector [ 0 ] [ span->at(1) ] + du * ( i - 1 );
                    c [ 3 ].at(2) = knotVector [ 1 ] [ span->at(2) ] + dv * j;

                    for ( k = 0; k < 4; k++ ) {
                        interp->local2global(cg [ k ], c [ k ], FEIIGAElementGeometryWrapper( this, iRule->giveKnotSpan() ));
                        p [ k ].x = ( FPNum ) cg [ k ].at(1);
                        p [ k ].y = ( FPNum ) cg [ k ].at(2);
                        p [ k ].z = 0.;

                        //move sampling gp out of boundary to overcome degeneracy on quarter plate with hole modelled by single patch
                        if ( c [ k ].at(1) > 0.99999 && c [ k ].at(2) > 0.495 && c [ k ].at(2) < 0.505 ) {
                            c [ k ].at(1) += sign [ k ].at(1) * du / 10.0;
                            c [ k ].at(2) += sign [ k ].at(2) * dv / 10.0;
                        }

                        // create a dummy ip's
                        FloatArray *cc = new FloatArray (c[k]);         // constructor of gp does not make its own copy
                        GaussPoint gp(iRule, 999, cc, 1.0, _PlaneStress);
#ifdef COMPUTE_STRAIN
                        this->computeStrainVector (val, &gp, tStep);
#endif
#ifdef COMPUTE_STRESS
                        this->computeStressVector(val, & gp, tStep);
#endif
                        s [ k ] = val.at(indx);

#ifdef QUARTER_PLATE_WITH_HOLE
                        if( indx == 2 ) {
                            if(cg[k].at(1) <= 1.0 + 1.0e-10 && cg[k].at(2) <= 1.0e-10 ) {
                                fprintf(stderr, "A: syy = %e\n", s[k]);
                            }
                        }
                        if( indx == 1 ) {
                            if(cg[k].at(1) <= 1.0e-10 && cg[k].at(2) <= 1.0 + 1.0e-10 ) {
                                fprintf(stderr, "B: sxx = %e\n", s[k]);
                            }
                        }

                        double x, y, r, phi, rate, E, G, kap, ny;
                        x = cg[k].at(1);
                        y = cg[k].at(2);
                        if( x < 1.0e-10 ){
                            phi = M_PI / 2.0;
                            r = y;
                        }
                        else {
                            phi = atan(y/x);
                            r = x/cos(phi);
                        }

#ifdef STRESS
                        // exact stresses quarter plate with hole s0=1 a=1
                        rate=1.0/r;
                        rate *= rate;
                        if ( indx == 1 ) s[k]=0.5*(2.0+3.0*rate*rate*cos(4.0*phi)-rate*(3*cos(2.0*phi)+2.0*cos(4.0*phi)));
                        if ( indx == 2 ) s[k]=0.5*(-3.0*rate*rate*cos(4.0*phi)-rate*(cos(2.0*phi)-2.0*cos(4.0*phi)));
                        if ( indx == 3 ) s[k]=0.5*(3.0*rate*rate*sin(4.0*phi)-rate*(sin(2.0*phi)+2.0*sin(4.0*phi)));

                        if ( indx == 2 ) {
                            if(cg[k].at(1) <= 1.0 + 1.0e-10 && cg[k].at(2) <= 1.0e-10){
                                fprintf(stderr, "A: syy = %e\n", s[k]);
                            }
                        }
                        if ( indx == 1 ) {
                            if(cg[k].at(1) <= 1.0e-10 && cg[k].at(2) <= 1.0 + 1.0e-10){
                                fprintf(stderr, "B: sxx = %e\n", s[k]);
                            }
                        }
#endif
#ifdef DISPL
                        // exact displ quarter plate with hole s0=1 a=1
                        E = 15000.0;
                        ny = 0.25;
                        G = E/(2.0*(1.0+ny));
                        kap = (3.0 - ny)/(1.0+ny);
                        rate=1.0/r;
                        if ( indx == 1 ) s[k]=(r*(kap+1.0)*cos(phi)+2.0*rate*((1.0+kap)*cos(phi)+cos(3.0*phi))-2.0*rate*rate*rate*cos(3.0*phi))/(8.0*G);
                        if ( indx == 2 ) s[k]=(r*(kap-3.0)*sin(phi)+2.0*rate*((1.0-kap)*sin(phi)+sin(3.0*phi))-2.0*rate*rate*rate*sin(3.0*phi))/(8.0*G);

                        if( indx == 1 ) {
                            if(cg[k].at(1) <= 1.0 + 1.0e-10 && cg[k].at(2) <= 1.0e-10){
                                fprintf(stderr, "A: ux = %e\n", s[k]);
                            }
                        }
                        if (indx == 2 ) {
                            if(cg[k].at(1) <= 1.0e-10 && cg[k].at(2) <= 1.0 + 1.0e-10){
                                fprintf(stderr, "B: uy = %e\n", s[k]);
                            }
                        }
#endif

                        if(s[k] < mins)mins=s[k];
                        if(s[k] > maxs)maxs=s[k];
#endif

                    }

                    if ( ( isnan(s [ 0 ]) ) || ( isnan(s [ 1 ]) ) || ( isnan(s [ 2 ]) ) || ( isnan(s [ 3 ]) ) ) {
                        continue;
                    }

                    if ( ( fabs(s [ 0 ]) > 1.e5 ) || ( fabs(s [ 1 ]) > 1.e5 ) || ( fabs(s [ 2 ]) > 1.e5 ) || ( fabs(s [ 3 ]) > 1.e5 ) ) {
                        continue;
                    }

                    //printf ("QWD: %e %e %e %e\n", s [ 0 ], s [ 1 ], s [ 2 ], s [ 3 ]);

                    go =  CreateQuadWD3D(p, s [ 0 ], s [ 1 ], s [ 2 ], s [ 3 ]);
                    EGWithMaskChangeAttributes(WIDTH_MASK | FILL_MASK | COLOR_MASK | EDGE_COLOR_MASK | EDGE_FLAG_MASK | LAYER_MASK, go);
                    EGAttachObject(go, ( EObjectP ) this);
                    EMAddGraphicsToModel(ESIModel(), go);
                }
            }
        }
    } // end loop over knot spans (irules)

#ifdef QUARTER_PLATE_WITH_HOLE
    fprintf(stderr, "%d: %e %e %e %e\n", indx, mins, maxs, (10.0*mins+maxs)/11.0, (10.0*maxs+mins)/11.0);
#endif
}


// refinement of integration elements should be generalized in terms of appropriately multiplying
// nseq for those integration elements which span more tham just a single knot span
// the reason is to ensure compatible division to quads over which scalar quantity is interpolated
// bilinearly !!!

void TSplinePlaneStressElement :: drawScalar(oofegGraphicContext &context) {
    int indx;
    WCRec p [ 4 ];
    GraphicObj *go;
    TimeStep *tStep = this->giveDomain()->giveEngngModel()->giveCurrentStep();
    double s [ 4 ];
    IntArray map;
    FEInterpolation *interp = this->giveInterpolation();
    FloatArray val;

    if ( !context.testElementGraphicActivity(this) ) {
        return;
    }

    EASValsSetLayer(OOFEG_VARPLOT_PATTERN_LAYER);
    EASValsSetFillStyle(FILL_SOLID);
    EASValsSetEdgeFlag(true);
    numberOfIntegrationRules = this->giveNumberOfIntegrationRules();
    const double * const * knotVector = interp->giveKnotVector();
    const IntArray *span;
    IntegrationRule *iRule;
    int ir, j, nsd = this->giveNsd();
    FloatArray c [ 4 ], cg [ 4 ];
    IntArray sign [ 4 ];

    if ( nsd == 2 ) {
        for ( j = 0; j < 4; j++ ) {
            c [ j ].resize(2);
            cg [ j ].resize(2);
            sign [ j ].resize(2);
        }
        sign [ 0 ].at(1) = 1;
        sign [ 0 ].at(2) = 1;
        sign [ 1 ].at(1) = -1;
        sign [ 1 ].at(2) = 1;
        sign [ 2 ].at(1) = -1;
        sign [ 2 ].at(2) = -1;
        sign [ 3 ].at(1) = 1;
        sign [ 3 ].at(2) = -1;
    } else {
        OOFEM_ERROR2("drawRawGeometry: not implemented for nsd = %d", nsd);
    }

    this->giveIntVarCompFullIndx( map, context.giveIntVarType() );
    if ( ( indx = map.at( context.giveIntVarIndx() ) ) == 0 ) {
        return;
    }

    // loop over individual integration rules (i.e., knot spans)
    for ( ir = 0; ir < numberOfIntegrationRules; ir++ ) {
        iRule = this->giveIntegrationRule(ir);
        span = iRule->giveKnotSpan();
        if ( nsd == 2 ) {
            // divide span locally to get finer geometry rep.
            int i, j, k, nseg = 4;
            double du = ( knotVector [ 0 ] [ span->at(1) + 1 ] - knotVector [ 0 ] [ span->at(1) ] ) / nseg;
            double dv = ( knotVector [ 1 ] [ span->at(2) + 1 ] - knotVector [ 1 ] [ span->at(2) ] ) / nseg;
            for ( i = 1; i <= nseg; i++ ) {
                for ( j = 1; j <= nseg; j++ ) {
                    c [ 0 ].at(1) = knotVector [ 0 ] [ span->at(1) ] + du * ( i - 1 );
                    c [ 0 ].at(2) = knotVector [ 1 ] [ span->at(2) ] + dv * ( j - 1 );
                    c [ 1 ].at(1) = knotVector [ 0 ] [ span->at(1) ] + du * i;
                    c [ 1 ].at(2) = knotVector [ 1 ] [ span->at(2) ] + dv * ( j - 1 );
                    c [ 2 ].at(1) = knotVector [ 0 ] [ span->at(1) ] + du * i;
                    c [ 2 ].at(2) = knotVector [ 1 ] [ span->at(2) ] + dv * j;
                    c [ 3 ].at(1) = knotVector [ 0 ] [ span->at(1) ] + du * ( i - 1 );
                    c [ 3 ].at(2) = knotVector [ 1 ] [ span->at(2) ] + dv * j;

                    for ( k = 0; k < 4; k++ ) {
                        interp->local2global(cg [ k ], c [ k ], FEIIGAElementGeometryWrapper( this, iRule->giveKnotSpan() ));
                        p [ k ].x = ( FPNum ) cg [ k ].at(1);
                        p [ k ].y = ( FPNum ) cg [ k ].at(2);
                        p [ k ].z = 0.;

#ifdef QUARTER_PLATE_WITH_HOLE_SINGLE_PATCH
                        //move sampling gp out of boundary to overcome degeneracy on quarter plate with hole modelled by single patch
                        if(c[k].at(1) > 0.99999 && c[k].at(2) > 0.495 && c[k].at(2) < 0.505){
                            c[k].at(1) += sign[k].at(1)*du/10.0;
                            c[k].at(2) += sign[k].at(2)*dv/10.0;
                        }
#endif

                        // create a dummy ip's
                        FloatArray *cc = new FloatArray (c[k]);         // constructor of gp does not make its own copy
                        GaussPoint gp(iRule, 999, cc, 1.0, _PlaneStress);
#ifdef COMPUTE_STRAIN
                        this->computeStrainVector (val, &gp, tStep);
#endif
#ifdef COMPUTE_STRESS
                        this->computeStressVector(val, & gp, tStep);
#endif
                        s [ k ] = val.at(indx);
                    }

                    if ( ( isnan(s [ 0 ]) ) || ( isnan(s [ 1 ]) ) || ( isnan(s [ 2 ]) ) || ( isnan(s [ 3 ]) ) ) {
                        continue;
                    }

                    if ( ( fabs(s [ 0 ]) > 1.e5 ) || ( fabs(s [ 1 ]) > 1.e5 ) || ( fabs(s [ 2 ]) > 1.e5 ) || ( fabs(s [ 3 ]) > 1.e5 ) ) {
                        continue;
                    }

                    //printf ("QWD: %e %e %e %e\n", s [ 0 ], s [ 1 ], s [ 2 ], s [ 3 ]);

                    go =  CreateQuadWD3D(p, s [ 0 ], s [ 1 ], s [ 2 ], s [ 3 ]);
                    EGWithMaskChangeAttributes(WIDTH_MASK | FILL_MASK | COLOR_MASK | EDGE_COLOR_MASK | EDGE_FLAG_MASK | LAYER_MASK, go);
                    EGAttachObject(go, ( EObjectP ) this);
                    EMAddGraphicsToModel(ESIModel(), go);
                }
            }
        }
    } // end loop over knot spans (irules)

}



void NURBSSpace3dElement :: drawScalar(oofegGraphicContext &context) {
    int indx;
    WCRec p [ 8 ];
    GraphicObj *go;
    TimeStep *tStep = this->giveDomain()->giveEngngModel()->giveCurrentStep();
    double s [ 8 ];
    IntArray map;
    FEInterpolation *interp = this->giveInterpolation();
    FloatArray val;
    //int huhu = 0;

    if ( !context.testElementGraphicActivity(this) ) {
        return;
    }

    EASValsSetLayer(OOFEG_VARPLOT_PATTERN_LAYER);
    EASValsSetFillStyle(FILL_SOLID);
    EASValsSetEdgeFlag(true);
    numberOfIntegrationRules = this->giveNumberOfIntegrationRules();
    const double * const * knotVector = interp->giveKnotVector();
    const IntArray *span;
    IntegrationRule *iRule;
    int ir, j, nsd = this->giveNsd();
    FloatArray c [ 8 ], cg [ 8 ];
    IntArray sign [ 8 ];

    if ( nsd == 3 ) {
        for ( j = 0; j < 8; j++ ) {
            c [ j ].resize(3);
            cg [ j ].resize(3);
            sign [ j ].resize(3);
        }

        sign [ 0 ].at(1) = 1;
        sign [ 0 ].at(2) = 1;
        sign [ 0 ].at(3) = 1;
        sign [ 1 ].at(1) = -1;
        sign [ 1 ].at(2) = 1;
        sign [ 1 ].at(3) = 1;
        sign [ 2 ].at(1) = -1;
        sign [ 2 ].at(2) = -1;
        sign [ 2 ].at(3) = 1;
        sign [ 3 ].at(1) = 1;
        sign [ 3 ].at(2) = -1;
        sign [ 3 ].at(3) = 1;
        sign [ 4 ].at(1) = 1;
        sign [ 4 ].at(2) = 1;
        sign [ 4 ].at(3) = -1;
        sign [ 5 ].at(1) = -1;
        sign [ 5 ].at(2) = 1;
        sign [ 5 ].at(3) = -1;
        sign [ 6 ].at(1) = -1;
        sign [ 6 ].at(2) = -1;
        sign [ 6 ].at(3) = -1;
        sign [ 7 ].at(1) = 1;
        sign [ 7 ].at(2) = -1;
        sign [ 7 ].at(3) = -1;
    } else {
        OOFEM_ERROR2("drawRawGeometry: not implemented for nsd = %d", nsd);
    }

    this->giveIntVarCompFullIndx( map, context.giveIntVarType() );
    if ( ( indx = map.at( context.giveIntVarIndx() ) ) == 0 ) {
        return;
    }

#ifdef SPHERICAL_CS
    if(context.giveIntVarType() == IST_StrainTensor){
        huhu = 1;
    }
#endif
#ifdef MISSES_STRESS
    if(context.giveIntVarType() == IST_ErrorIndicatorLevel){
        huhu = 2;
    }
#endif

    //double maxs=-1.0e10, mins=1.0e10;

    // loop over individual integration rules (i.e., knot spans)
    for ( ir = 0; ir < numberOfIntegrationRules; ir++ ) {
        iRule = this->giveIntegrationRule(ir);
        span = iRule->giveKnotSpan();
        if ( nsd == 3 ) {
            // divide span locally to get finer geometry rep.
            int i, j, k, m, nseg = 8;
            double du = ( knotVector [ 0 ] [ span->at(1) + 1 ] - knotVector [ 0 ] [ span->at(1) ] ) / nseg;
            double dv = ( knotVector [ 1 ] [ span->at(2) + 1 ] - knotVector [ 1 ] [ span->at(2) ] ) / nseg;
            double dw = ( knotVector [ 2 ] [ span->at(3) + 1 ] - knotVector [ 2 ] [ span->at(3) ] ) / nseg;

#ifdef DRAW_VISIBLE_CONTOUR
            WCRec pp [ 4 ];
            double ss[4];
            int kk, ii, nd[6][4] = {{0, 1, 2, 3}, {4, 5, 6, 7}, {0, 1, 5, 4}, {1, 2, 6, 5}, {2, 3, 7, 6}, {3, 0, 4, 7}};
#endif

           for ( i = 1; i <= nseg; i++ ) {
                for ( j = 1; j <= nseg; j++ ) {
                    for ( m = 1; m <= nseg; m++ ) {
                        c [ 0 ].at(1) = knotVector [ 0 ] [ span->at(1) ] + du * ( i - 1 );
                        c [ 0 ].at(2) = knotVector [ 1 ] [ span->at(2) ] + dv * ( j - 1 );
                        c [ 0 ].at(3) = knotVector [ 2 ] [ span->at(3) ] + dw * ( m - 1 );
                        c [ 1 ].at(1) = knotVector [ 0 ] [ span->at(1) ] + du * i;
                        c [ 1 ].at(2) = knotVector [ 1 ] [ span->at(2) ] + dv * ( j - 1 );
                        c [ 1 ].at(3) = knotVector [ 2 ] [ span->at(3) ] + dw * ( m - 1 );
                        c [ 2 ].at(1) = knotVector [ 0 ] [ span->at(1) ] + du * i;
                        c [ 2 ].at(2) = knotVector [ 1 ] [ span->at(2) ] + dv * j;
                        c [ 2 ].at(3) = knotVector [ 2 ] [ span->at(3) ] + dw * ( m - 1 );
                        c [ 3 ].at(1) = knotVector [ 0 ] [ span->at(1) ] + du * ( i - 1 );
                        c [ 3 ].at(2) = knotVector [ 1 ] [ span->at(2) ] + dv * j;
                        c [ 3 ].at(3) = knotVector [ 2 ] [ span->at(3) ] + dw * ( m - 1 );
                        c [ 4 ].at(1) = knotVector [ 0 ] [ span->at(1) ] + du * ( i - 1 );
                        c [ 4 ].at(2) = knotVector [ 1 ] [ span->at(2) ] + dv * ( j - 1 );
                        c [ 4 ].at(3) = knotVector [ 2 ] [ span->at(3) ] + dw * m;
                        c [ 5 ].at(1) = knotVector [ 0 ] [ span->at(1) ] + du * i;
                        c [ 5 ].at(2) = knotVector [ 1 ] [ span->at(2) ] + dv * ( j - 1 );
                        c [ 5 ].at(3) = knotVector [ 2 ] [ span->at(3) ] + dw * m;
                        c [ 6 ].at(1) = knotVector [ 0 ] [ span->at(1) ] + du * i;
                        c [ 6 ].at(2) = knotVector [ 1 ] [ span->at(2) ] + dv * j;
                        c [ 6 ].at(3) = knotVector [ 2 ] [ span->at(3) ] + dw * m;
                        c [ 7 ].at(1) = knotVector [ 0 ] [ span->at(1) ] + du * ( i - 1 );
                        c [ 7 ].at(2) = knotVector [ 1 ] [ span->at(2) ] + dv * j;
                        c [ 7 ].at(3) = knotVector [ 2 ] [ span->at(3) ] + dw * m;

                        for ( k = 0; k < 8; k++ ) {
                            interp->local2global(cg [ k ], c [ k ], FEIIGAElementGeometryWrapper( this, iRule->giveKnotSpan() ));
                            p [ k ].x = ( FPNum ) cg [ k ].at(1);
                            p [ k ].y = ( FPNum ) cg [ k ].at(2);
                            p [ k ].z = ( FPNum ) cg [ k ].at(3);
                        }

#ifdef DRAW_VISIBLE_CONTOUR
#ifdef SPHERE_WITH_HOLE

                        // check whether drawing side visible in particular view !!!
                        int haha = 0;
                        for(kk=0;kk<6;kk++){
                            double xx = 0.0, yy = 0.0, zz = 0.0, rr, r;
                            for(ii=0;ii<4;ii++){
                                xx += (pp[ii].x=p[nd[kk][ii]].x);
                                yy += (pp[ii].y=p[nd[kk][ii]].y);
                                zz += (pp[ii].z=p[nd[kk][ii]].z);
                                ss[ii] = s[nd[kk][ii]];
                            }
                            xx /= 4.0;
                            yy /= 4.0;
                            zz /= 4.0;
                            rr = xx * xx + yy * yy;
                            r = rr + zz * zz;

                            if(zz < 2.0001 /* || xx < 0.0001 */|| yy < 0.0001 || rr < 1.0 || r < 25.0 || r > 5.99 * 5.99)haha = 1;
                        }
                        if(haha == 0)continue;
#endif
#endif

                        for (k=0;k<8; k++) {
                            // create a dummy ip's
                            FloatArray *cc = new FloatArray (c[k]);         // constructor of gp does not make its own copy
                            GaussPoint gp(iRule, 999, cc, 1.0, _3dMat);
#ifdef COMPUTE_STRAIN
                            this->computeStrainVector (val, &gp, tStep);
#endif
#ifdef COMPUTE_STRESS
                            this->computeStressVector(val, & gp, tStep);
#endif
                            s [ k ] = val.at(indx);

#ifdef MISSES_STRESS
                            if(huhu == 2){
                                double vonMisses;
                                vonMisses = sqrt(((val.at(1)-val.at(2))*(val.at(1)-val.at(2)) + (val.at(2)-val.at(3))*(val.at(2)-val.at(3)) + (val.at(1)-val.at(3))*(val.at(1)-val.at(3)) + 6.0*(val.at(4)*val.at(4)+val.at(5)*val.at(5)+val.at(6)*val.at(6)))/2.0);
                                s[k]=vonMisses;
                            }
#endif
#ifdef SPHERICAL_CS
                            if(huhu == 1){
                                double x, y, z, r, r2, rr;
                                FloatMatrix sigma(3,3), T(3,3), product(3,3), result(3,3);

                                sigma.at(1,1)=val.at(1);
                                sigma.at(1,2)=val.at(6);
                                sigma.at(1,3)=val.at(5);
                                sigma.at(2,1)=val.at(6);
                                sigma.at(2,2)=val.at(2);
                                sigma.at(2,3)=val.at(4);
                                sigma.at(3,1)=val.at(5);
                                sigma.at(3,2)=val.at(4);
                                sigma.at(3,3)=val.at(3);

                                x=cg[k].at(1);
                                y=cg[k].at(2);
                                z=cg[k].at(3);
                                r2=x*x+y*y+z*z;
                                r=sqrt(r2);
                                rr=sqrt(x*x+y*y);

                                T.at(1,1)=-z*z*y/rr/r2-y*rr/r2;
                                T.at(1,2)=z*z*x/rr/r2+x*rr/r2;
                                T.at(1,3)=0.0;
                                T.at(2,1)=-z*x/rr/r;
                                T.at(2,2)=-z*y/rr/r;
                                T.at(2,3)=rr/r;
                                T.at(3,1)=x/r;
                                T.at(3,2)=y/r;
                                T.at(3,3)=z/r;

                                product.beProductOf(T,sigma);
                                result.beProductTOf(product,T);

                                val.at(1)=result.at(1,1);
                                val.at(6)=result.at(1,2);
                                val.at(5)=result.at(1,3);
                                val.at(6)=result.at(2,1);
                                val.at(2)=result.at(2,2);
                                val.at(4)=result.at(2,3);
                                val.at(5)=result.at(3,1);
                                val.at(4)=result.at(3,2);
                                val.at(3)=result.at(3,3);

                                s[k]=val.at(indx);
                            }
#endif
#ifdef SPHERE_WITH_HOLE
                            if(s[k] < mins)mins=s[k];
                            if(s[k] > maxs)maxs=s[k];
#endif
                        }

                        if ( ( isnan(s [ 0 ]) ) || ( isnan(s [ 1 ]) ) || ( isnan(s [ 2 ]) ) || ( isnan(s [ 3 ]) ) ) {
                            continue;
                        }

                        if ( ( isnan(s [ 4 ]) ) || ( isnan(s [ 5 ]) ) || ( isnan(s [ 6 ]) ) || ( isnan(s [ 7 ]) ) ) {
                            continue;
                        }

                        if ( ( fabs(s [ 0 ]) > 1.e5 ) || ( fabs(s [ 1 ]) > 1.e5 ) || ( fabs(s [ 2 ]) > 1.e5 ) || ( fabs(s [ 3 ]) > 1.e5 ) ) {
                            continue;
                        }

                        if ( ( fabs(s [ 4 ]) > 1.e5 ) || ( fabs(s [ 5 ]) > 1.e5 ) || ( fabs(s [ 6 ]) > 1.e5 ) || ( fabs(s [ 7 ]) > 1.e5 ) ) {
                            continue;
                        }

                        //printf ("HWD: %e %e %e %e %e %e %e %e\n", s [ 0 ], s [ 1 ], s [ 2 ], s [ 3 ], s [ 4 ], s [ 5 ], s [ 6 ], s [ 7 ]);

#ifdef DRAW_VISIBLE_CONTOUR
#ifdef SPHERE_WITH_HOLE
                        for( kk=0; kk<6; kk++ ) {
                            double xx = 0.0, yy = 0.0, zz = 0.0, rr, r;
                            for(ii=0;ii<4;ii++){
                                xx += (pp[ii].x=p[nd[kk][ii]].x);
                                yy += (pp[ii].y=p[nd[kk][ii]].y);
                                zz += (pp[ii].z=p[nd[kk][ii]].z);
                                ss[ii] = s[nd[kk][ii]];
                            }
                            xx /= 4.0;
                            yy /= 4.0;
                            zz /= 4.0;
                            rr = xx * xx + yy * yy;
                            r = rr + zz * zz;

                            if ( zz < 2.0001 /* || xx < 0.0001 */|| yy < 0.0001 || rr < 1.0 || r < 25.0 || r > 5.99 * 5.99 ) {
                                go = CreateQuadWD3D(pp, ss[0], ss[1], ss[2], ss[3]);
                                EGWithMaskChangeAttributes(LAYER_MASK | EDGE_COLOR_MASK | EDGE_FLAG_MASK, go);
                                EMAddGraphicsToModel(ESIModel(), go);
                            }
                        }
#endif
#else
                        go = CreateHexahedronWD(p, s);
                        EGWithMaskChangeAttributes(WIDTH_MASK | FILL_MASK | COLOR_MASK | EDGE_COLOR_MASK | EDGE_FLAG_MASK | LAYER_MASK, go);
                        EGAttachObject(go, ( EObjectP ) this);
                        EMAddGraphicsToModel(ESIModel(), go);
#endif
                    }
                }
            }
        }
    } // end loop over knot spans (irules)

#ifdef SPHERE_WITH_HOLE
    fprintf(stderr, "%d %e %e %e %e\n", indx, mins, maxs, (10.0*mins+maxs)/11.0, (10.0*maxs+mins)/11.0);
#endif

}


#endif

} // end namespace oofem
