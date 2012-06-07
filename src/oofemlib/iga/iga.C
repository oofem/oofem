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

#include "inputrecord.h"
#include "flotarry.h"
#include "flotmtrx.h"
#include "mathfem.h"
#include "iga.h"
#include "feitspline.h"

#ifdef __OOFEG
 #include "oofeggraphiccontext.h"
 #include "oofegutils.h"
 #include "engngm.h"
 #include "structuralelementevaluator.h"
#endif


namespace oofem {

IRResultType IGAElement :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                 // Required by IR_GIVE_FIELD macro

    int indx = 0, ui, vi, wi, i, nsd, numberOfGaussPoints = 1, numberOfKnotSpans=1;
    double du, dv, dw;
    const FloatArray *gpcoords;
    FloatArray newgpcoords;
    IntArray knotSpan;

    Element :: initializeFrom(ir); // read nodes , material, cross section
    // set number of dofmanagers
    this->numberOfDofMans = dofManArray.giveSize();
    this->giveInterpolation()->initializeFrom(ir); // read geometry

    IR_GIVE_OPTIONAL_FIELD(ir, numberOfGaussPoints, IFT_IGAElement_NIP, "nip"); // Macro

    // generate individual IntegrationElements; one for each nonzero knot span
    nsd = this->giveNsd();
    if ( nsd == 1 ) {
        //HUHU
    } else if ( nsd == 2 ) {
        int numberOfKnotSpansU = this->giveInterpolation()->giveNumberOfKnotSpans(1);
        int numberOfKnotSpansV = this->giveInterpolation()->giveNumberOfKnotSpans(2);
        numberOfKnotSpans = numberOfKnotSpansU*numberOfKnotSpansV;
        const IntArray * knotMultiplicityU = this->giveInterpolation()->giveKnotMultiplicity(1);
        const IntArray * knotMultiplicityV = this->giveInterpolation()->giveKnotMultiplicity(2);
        const FloatArray * knotValuesU = this->giveInterpolation()->giveKnotValues(1);
        const FloatArray * knotValuesV = this->giveInterpolation()->giveKnotValues(2);

        newgpcoords.resize(2);
        knotSpan.resize(2);

        this->numberOfIntegrationRules = numberOfKnotSpansU * numberOfKnotSpansV;
        integrationRulesArray = new IntegrationRule * [ numberOfIntegrationRules ];

        knotSpan.at(2) = -1;
        for ( vi = 1; vi <= numberOfKnotSpansV; vi++ ) {
            dv = knotValuesV->at(vi + 1) - knotValuesV->at(vi);
            knotSpan.at(2) += knotMultiplicityV->at(vi);

            knotSpan.at(1) = -1;
            for ( ui = 1; ui <= numberOfKnotSpansU; ui++ ) {
                du = knotValuesU->at(ui + 1) - knotValuesU->at(ui);
                knotSpan.at(1) += knotMultiplicityU->at(ui);

                integrationRulesArray [ indx ] = new IGAIntegrationElement(indx, this, knotSpan);
                integrationRulesArray [ indx ]->setUpIntegrationPoints(_Square, numberOfGaussPoints, _PlaneStress); // HUHU _PlaneStress, rectangle

                // remap local subelement gp coordinates into knot span coordinates and update integration weight
                for ( i = 0; i < numberOfGaussPoints; i++ ) {
                    gpcoords = integrationRulesArray [ indx ]->getIntegrationPoint(i)->giveCoordinates();

                    newgpcoords.at(1) = knotValuesU->at(ui) + du * ( gpcoords->at(1) / 2.0 + 0.5 );
                    newgpcoords.at(2) = knotValuesV->at(vi) + dv * ( gpcoords->at(2) / 2.0 + 0.5 );
                    integrationRulesArray [ indx ]->getIntegrationPoint(i)->setCoordinates(newgpcoords);
                    integrationRulesArray [ indx ]->getIntegrationPoint(i)->setWeight(integrationRulesArray [ indx ]->getIntegrationPoint(i)->giveWeight() / 4.0 * du * dv);
                }

                indx++;
            }
        }
    } else if ( nsd == 3 ) {
        int numberOfKnotSpansU = this->giveInterpolation()->giveNumberOfKnotSpans(1);
        int numberOfKnotSpansV = this->giveInterpolation()->giveNumberOfKnotSpans(2);
        int numberOfKnotSpansW = this->giveInterpolation()->giveNumberOfKnotSpans(3);
        numberOfKnotSpans = numberOfKnotSpansU*numberOfKnotSpansV*numberOfKnotSpansW;
        const IntArray * knotMultiplicityU = this->giveInterpolation()->giveKnotMultiplicity(1);
        const IntArray * knotMultiplicityV = this->giveInterpolation()->giveKnotMultiplicity(2);
        const IntArray * knotMultiplicityW = this->giveInterpolation()->giveKnotMultiplicity(3);
        const FloatArray * knotValuesU = this->giveInterpolation()->giveKnotValues(1);
        const FloatArray * knotValuesV = this->giveInterpolation()->giveKnotValues(2);
        const FloatArray * knotValuesW = this->giveInterpolation()->giveKnotValues(3);

        newgpcoords.resize(3);
        knotSpan.resize(3);

        this->numberOfIntegrationRules = numberOfKnotSpansU * numberOfKnotSpansV * numberOfKnotSpansW;
        integrationRulesArray = new IntegrationRule * [ numberOfIntegrationRules ];

        knotSpan.at(3) = -1;
        for ( wi = 1; wi <= numberOfKnotSpansW; wi++ ) {
            dw = knotValuesW->at(wi + 1) - knotValuesW->at(wi);
            knotSpan.at(3) += knotMultiplicityW->at(wi);

            knotSpan.at(2) = -1;
            for ( vi = 1; vi <= numberOfKnotSpansV; vi++ ) {
                dv = knotValuesV->at(vi + 1) - knotValuesV->at(vi);
                knotSpan.at(2) += knotMultiplicityV->at(vi);

                knotSpan.at(1) = -1;
                for ( ui = 1; ui <= numberOfKnotSpansU; ui++ ) {
                    du = knotValuesU->at(ui + 1) - knotValuesU->at(ui);
                    knotSpan.at(1) += knotMultiplicityU->at(ui);

                    integrationRulesArray [ indx ] = new IGAIntegrationElement(indx, this, knotSpan);
                    integrationRulesArray [ indx ]->setUpIntegrationPoints(_Cube, numberOfGaussPoints, _3dMat);

                    // remap local subelement gp coordinates into knot span coordinates and update integration weight
                    for ( i = 0; i < numberOfGaussPoints; i++ ) {
                        gpcoords = integrationRulesArray [ indx ]->getIntegrationPoint(i)->giveCoordinates();

                        newgpcoords.at(1) = knotValuesU->at(ui) + du * ( gpcoords->at(1) / 2.0 + 0.5 );
                        newgpcoords.at(2) = knotValuesV->at(vi) + dv * ( gpcoords->at(2) / 2.0 + 0.5 );
                        newgpcoords.at(3) = knotValuesW->at(wi) + dw * ( gpcoords->at(3) / 2.0 + 0.5 );
                        integrationRulesArray [ indx ]->getIntegrationPoint(i)->setCoordinates(newgpcoords);
                        integrationRulesArray [ indx ]->getIntegrationPoint(i)->setWeight(integrationRulesArray [ indx ]->getIntegrationPoint(i)->giveWeight() / 8.0 * du * dv * dw);
                    }

                    indx++;
                }
            }
        }
    } else {
        OOFEM_ERROR2("unsupported number of spatial dimensions (nsd = %d)", nsd);
    }

#ifdef __PARALLEL_MODE
    // read optional knot span parallel mode
    int _i;
    this->knotSpanParallelMode.resize(numberOfKnotSpans);
    // set Element_local as default
    for (_i=1; _i<=numberOfKnotSpans; _i++) knotSpanParallelMode.at(_i)=Element_local;
    IR_GIVE_OPTIONAL_FIELD(ir, knotSpanParallelMode, IFT_IGAElement_KnotSpanParallelMode, "knotspanparmode"); // Macro
#endif


    return IRRT_OK;
}


#ifdef __PARALLEL_MODE
elementParallelMode
IGAElement::giveKnotSpanParallelMode(int knotSpanIndex) const
{
    elementParallelMode emode = this->giveParallelMode();
    if (emode == Element_remote) {
        return Element_remote;
    } else if (emode == Element_local) {
        return (elementParallelMode) this->knotSpanParallelMode.at(knotSpanIndex+1);
    } else {
        _error("Cannot determine elementParallelMode");
    }

    return Element_local;//to make compiler happy
}

#endif // __PARALLEL_MODE


// integration elements are setup in the same way as for IGAElement for now HUHU

IRResultType IGATSplineElement :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                 // Required by IR_GIVE_FIELD macro
    TSplineInterpolation *interpol = ( TSplineInterpolation * ) this->giveInterpolation();

    int indx = 0, ui, vi, i, nsd, numberOfGaussPoints = 1;
    double du, dv;
    const FloatArray *gpcoords;
    FloatArray newgpcoords;
    IntArray knotSpan;

    Element :: initializeFrom(ir); // read nodes , material, cross section
    // set number of dofmanagers
    this->numberOfDofMans = dofManArray.giveSize();
    // set number of control points before initialization HUHU HAHA
    interpol->setNumberOfControlPoints(this->numberOfDofMans);
    this->giveInterpolation()->initializeFrom(ir); // read geometry

    IR_GIVE_OPTIONAL_FIELD(ir, numberOfGaussPoints, IFT_IGAElement_NIP, "nip"); // Macro

    // generate individual IntegrationElements; one for each nonzero knot span
    nsd = giveNsd();
    if ( nsd == 2 ) {
        int numberOfKnotSpansU = this->giveInterpolation()->giveNumberOfKnotSpans(1);
        int numberOfKnotSpansV = this->giveInterpolation()->giveNumberOfKnotSpans(2);
        const IntArray * knotMultiplicityU = this->giveInterpolation()->giveKnotMultiplicity(1);
        const IntArray * knotMultiplicityV = this->giveInterpolation()->giveKnotMultiplicity(2);
        const FloatArray * knotValuesU = this->giveInterpolation()->giveKnotValues(1);
        const FloatArray * knotValuesV = this->giveInterpolation()->giveKnotValues(2);

        newgpcoords.resize(2);
        knotSpan.resize(2);

        this->numberOfIntegrationRules = numberOfKnotSpansU * numberOfKnotSpansV;
        integrationRulesArray = new IntegrationRule * [ numberOfIntegrationRules ];

        knotSpan.at(2) = -1;
        for ( vi = 1; vi <= numberOfKnotSpansV; vi++ ) {
            dv = knotValuesV->at(vi + 1) - knotValuesV->at(vi);
            knotSpan.at(2) += knotMultiplicityV->at(vi);

            knotSpan.at(1) = -1;
            for ( ui = 1; ui <= numberOfKnotSpansU; ui++ ) {
                du = knotValuesU->at(ui + 1) - knotValuesU->at(ui);
                knotSpan.at(1) += knotMultiplicityU->at(ui);

                integrationRulesArray [ indx ] = new IGAIntegrationElement(indx, this, knotSpan);
                integrationRulesArray [ indx ]->setUpIntegrationPoints(_Square, numberOfGaussPoints, _PlaneStress); // HUHU _PlaneStress, rectangle

                // remap local subelement gp coordinates into knot span coordinates and update integration weight
                for ( i = 0; i < numberOfGaussPoints; i++ ) {
                    gpcoords = integrationRulesArray [ indx ]->getIntegrationPoint(i)->giveCoordinates();

                    newgpcoords.at(1) = knotValuesU->at(ui) + du * ( gpcoords->at(1) / 2.0 + 0.5 );
                    newgpcoords.at(2) = knotValuesV->at(vi) + dv * ( gpcoords->at(2) / 2.0 + 0.5 );
                    integrationRulesArray [ indx ]->getIntegrationPoint(i)->setCoordinates(newgpcoords);
                    integrationRulesArray [ indx ]->getIntegrationPoint(i)->setWeight(integrationRulesArray [ indx ]->getIntegrationPoint(i)->giveWeight() / 4.0 * du * dv);
                }

                indx++;
            }
        }
    } else {
        OOFEM_ERROR2("unsupported number of spatial dimensions (nsd = %d)", nsd);
    }

    return IRRT_OK;
}





#ifdef __OOFEG

 #define DRAW_MESH

// if DRAW_MESH is defined only boundary of integration elements are drawn;
// currently mesh is not properly drawn for tsplines
// because integration elements (does not matter whether single span or multi span)
// are generaly finer than T-mesh;

void IGAElement :: drawRawGeometry(oofegGraphicContext &gc)
{
    WCRec p [ 8 ];
    GraphicObj *go;
    FEInterpolation *interp = this->giveInterpolation();
    int i, j, k, m, nseg;

 #ifdef DRAW_MESH
    WCRec pp [ 2 ];
 #endif

    if ( !gc.testElementGraphicActivity(this) ) {
        return;
    }

    EASValsSetLineWidth(OOFEG_RAW_GEOMETRY_WIDTH);
    EASValsSetColor( gc.getElementColor() );
    EASValsSetEdgeColor( gc.getElementEdgeColor() );
    EASValsSetEdgeFlag(true);
    EASValsSetLayer(OOFEG_RAW_GEOMETRY_LAYER);
    EASValsSetLineStyle(SOLID_STYLE);
    EASValsSetFillStyle(FILL_SOLID);
    //EASValsSetLineWidth(0);

 #ifdef DRAW_MESH
    nseg = 8;
 #else
    nseq = 4;
 #endif

    int numberOfIntegrationRules = this->giveNumberOfIntegrationRules();
    const double *const* knotVector = interp->giveKnotVector();
    const IntArray *span;
    IntegrationRule *iRule;
    int ir, nsd = this->giveNsd();

    if ( nsd == 1 ) {
        FloatArray c [ 2 ], cg [ 2 ];
        double du;

        for ( j = 0; j < 2; j++ ) {
            c [ j ].resize(1);
            cg [ j ].resize(1);
        }

        // loop over individual integration rules (i.e., knot spans)
        for ( ir = 0; ir < numberOfIntegrationRules; ir++ ) {
            iRule = this->giveIntegrationRule(ir);
            span = iRule->giveKnotSpan();
            // divide span locally to get finer geometry rep.
            du = ( knotVector [ 0 ] [ span->at(1) + 1 ] - knotVector [ 0 ] [ span->at(1) ] ) / nseg;
            for ( i = 1; i <= nseg; i++ ) {
                c [ 0 ].at(1) = knotVector [ 0 ] [ span->at(1) ] + du * ( i - 1 );
                c [ 1 ].at(1) = knotVector [ 0 ] [ span->at(1) ] + du * i;

                for ( k = 0; k < 2; k++ ) {
                    interp->local2global(cg [ k ], c [ k ], FEIIGAElementGeometryWrapper( this, iRule->giveKnotSpan() ));
                    p [ k ].x = ( FPNum ) cg [ k ].at(1);
                    p [ k ].y = 0.;
                    p [ k ].z = 0.;
                }

                go =  CreateLine3D(p);
                EGWithMaskChangeAttributes(WIDTH_MASK | STYLE_MASK | COLOR_MASK | LAYER_MASK, go);
                EGAttachObject(go, ( EObjectP ) this);
                EMAddGraphicsToModel(ESIModel(), go);
            }
        }                 // end loop over knot spans (irules)

    } else if ( nsd == 2 )      {
        FloatArray c [ 4 ], cg [ 4 ];
        double du, dv;

        for ( j = 0; j < 4; j++ ) {
            c [ j ].resize(2);
            cg [ j ].resize(2);
        }

        // loop over individual integration rules (i.e., knot spans)
        for ( ir = 0; ir < numberOfIntegrationRules; ir++ ) {
            iRule = this->giveIntegrationRule(ir);
            span = iRule->giveKnotSpan();
            // divide span locally to get finer geometry rep.
            du = ( knotVector [ 0 ] [ span->at(1) + 1 ] - knotVector [ 0 ] [ span->at(1) ] ) / nseg;
            dv = ( knotVector [ 1 ] [ span->at(2) + 1 ] - knotVector [ 1 ] [ span->at(2) ] ) / nseg;
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
                    }

 #ifdef DRAW_MESH
                    if ( i == 1 ) {
                        pp [ 0 ].x = p [ 0 ].x;
                        pp [ 0 ].y = p [ 0 ].y;
                        pp [ 0 ].z = p [ 0 ].z;

                        pp [ 1 ].x = p [ 3 ].x;
                        pp [ 1 ].y = p [ 3 ].y;
                        pp [ 1 ].z = p [ 3 ].z;

                        go =  CreateLine3D(pp);
                        EGWithMaskChangeAttributes(WIDTH_MASK | STYLE_MASK | COLOR_MASK | LAYER_MASK, go);
                        EGAttachObject(go, ( EObjectP ) this);
                        EMAddGraphicsToModel(ESIModel(), go);
                    }

                    if ( j == 1 ) {
                        pp [ 0 ].x = p [ 0 ].x;
                        pp [ 0 ].y = p [ 0 ].y;
                        pp [ 0 ].z = p [ 0 ].z;

                        pp [ 1 ].x = p [ 1 ].x;
                        pp [ 1 ].y = p [ 1 ].y;
                        pp [ 1 ].z = p [ 1 ].z;

                        go =  CreateLine3D(pp);
                        EGWithMaskChangeAttributes(WIDTH_MASK | STYLE_MASK | COLOR_MASK | LAYER_MASK, go);
                        EGAttachObject(go, ( EObjectP ) this);
                        EMAddGraphicsToModel(ESIModel(), go);
                    }

                    if ( i == nseg ) {
                        pp [ 0 ].x = p [ 1 ].x;
                        pp [ 0 ].y = p [ 1 ].y;
                        pp [ 0 ].z = p [ 1 ].z;

                        pp [ 1 ].x = p [ 2 ].x;
                        pp [ 1 ].y = p [ 2 ].y;
                        pp [ 1 ].z = p [ 2 ].z;

                        go =  CreateLine3D(pp);
                        EGWithMaskChangeAttributes(WIDTH_MASK | STYLE_MASK | COLOR_MASK | LAYER_MASK, go);
                        EGAttachObject(go, ( EObjectP ) this);
                        EMAddGraphicsToModel(ESIModel(), go);
                    }

                    if ( j == nseg ) {
                        pp [ 0 ].x = p [ 2 ].x;
                        pp [ 0 ].y = p [ 2 ].y;
                        pp [ 0 ].z = p [ 2 ].z;

                        pp [ 1 ].x = p [ 3 ].x;
                        pp [ 1 ].y = p [ 3 ].y;
                        pp [ 1 ].z = p [ 3 ].z;

                        go =  CreateLine3D(pp);
                        EGWithMaskChangeAttributes(WIDTH_MASK | STYLE_MASK | COLOR_MASK | LAYER_MASK, go);
                        EGAttachObject(go, ( EObjectP ) this);
                        EMAddGraphicsToModel(ESIModel(), go);
                    }

 #else
                    go =  CreateQuad3D(p);
                    EGWithMaskChangeAttributes(WIDTH_MASK | FILL_MASK | COLOR_MASK | EDGE_COLOR_MASK | EDGE_FLAG_MASK | LAYER_MASK, go);
                    EGAttachObject(go, ( EObjectP ) this);
                    EMAddGraphicsToModel(ESIModel(), go);
 #endif
                }
            }
        }                 // end loop over knot spans (irules)

    } else if ( nsd == 3 )      {
        FloatArray c [ 8 ], cg [ 8 ];
        double du, dv, dt;

        for ( j = 0; j < 8; j++ ) {
            c [ j ].resize(3);
            cg [ j ].resize(3);
        }

        // loop over individual integration rules (i.e., knot spans)
        for ( ir = 0; ir < numberOfIntegrationRules; ir++ ) {
            iRule = this->giveIntegrationRule(ir);
            span = iRule->giveKnotSpan();
            // divide span locally to get finer geometry rep.
            du = ( knotVector [ 0 ] [ span->at(1) + 1 ] - knotVector [ 0 ] [ span->at(1) ] ) / nseg;
            dv = ( knotVector [ 1 ] [ span->at(2) + 1 ] - knotVector [ 1 ] [ span->at(2) ] ) / nseg;
            dt = ( knotVector [ 2 ] [ span->at(3) + 1 ] - knotVector [ 2 ] [ span->at(3) ] ) / nseg;
            for ( i = 1; i <= nseg; i++ ) {
                for ( j = 1; j <= nseg; j++ ) {
                    for ( k = 1; k <= nseg; k++ ) {
                        c [ 0 ].at(1) = knotVector [ 0 ] [ span->at(1) ] + du * ( i - 1 );
                        c [ 0 ].at(2) = knotVector [ 1 ] [ span->at(2) ] + dv * ( j - 1 );
                        c [ 0 ].at(3) = knotVector [ 2 ] [ span->at(3) ] + dt * ( k - 1 );
                        c [ 1 ].at(1) = knotVector [ 0 ] [ span->at(1) ] + du * i;
                        c [ 1 ].at(2) = knotVector [ 1 ] [ span->at(2) ] + dv * ( j - 1 );
                        c [ 1 ].at(3) = knotVector [ 2 ] [ span->at(3) ] + dt * ( k - 1 );
                        c [ 2 ].at(1) = knotVector [ 0 ] [ span->at(1) ] + du * i;
                        c [ 2 ].at(2) = knotVector [ 1 ] [ span->at(2) ] + dv * j;
                        c [ 2 ].at(3) = knotVector [ 2 ] [ span->at(3) ] + dt * ( k - 1 );
                        c [ 3 ].at(1) = knotVector [ 0 ] [ span->at(1) ] + du * ( i - 1 );
                        c [ 3 ].at(2) = knotVector [ 1 ] [ span->at(2) ] + dv * j;
                        c [ 3 ].at(3) = knotVector [ 2 ] [ span->at(3) ] + dt * ( k - 1 );
                        c [ 4 ].at(1) = knotVector [ 0 ] [ span->at(1) ] + du * ( i - 1 );
                        c [ 4 ].at(2) = knotVector [ 1 ] [ span->at(2) ] + dv * ( j - 1 );
                        c [ 4 ].at(3) = knotVector [ 2 ] [ span->at(3) ] + dt * k;
                        c [ 5 ].at(1) = knotVector [ 0 ] [ span->at(1) ] + du * i;
                        c [ 5 ].at(2) = knotVector [ 1 ] [ span->at(2) ] + dv * ( j - 1 );
                        c [ 5 ].at(3) = knotVector [ 2 ] [ span->at(3) ] + dt * k;
                        c [ 6 ].at(1) = knotVector [ 0 ] [ span->at(1) ] + du * i;
                        c [ 6 ].at(2) = knotVector [ 1 ] [ span->at(2) ] + dv * j;
                        c [ 6 ].at(3) = knotVector [ 2 ] [ span->at(3) ] + dt * k;
                        c [ 7 ].at(1) = knotVector [ 0 ] [ span->at(1) ] + du * ( i - 1 );
                        c [ 7 ].at(2) = knotVector [ 1 ] [ span->at(2) ] + dv * j;
                        c [ 7 ].at(3) = knotVector [ 2 ] [ span->at(3) ] + dt * k;

                        for ( m = 0; m < 8; m++ ) {
                            interp->local2global(cg [ m ], c [ m ], FEIIGAElementGeometryWrapper( this, iRule->giveKnotSpan() ));
                            p [ m ].x = ( FPNum ) cg [ m ].at(1);
                            p [ m ].y = ( FPNum ) cg [ m ].at(2);
                            p [ m ].z = ( FPNum ) cg [ m ].at(3);
                        }

 #ifdef DRAW_MESH
                        if ( i == 1 && j == 1 ) {
                            pp [ 0 ].x = p [ 0 ].x;
                            pp [ 0 ].y = p [ 0 ].y;
                            pp [ 0 ].z = p [ 0 ].z;

                            pp [ 1 ].x = p [ 4 ].x;
                            pp [ 1 ].y = p [ 4 ].y;
                            pp [ 1 ].z = p [ 4 ].z;

#ifdef DRAW_VISIBLE_CONTOUR
#ifdef SPHERE_WITH_HOLE
                            double xx = 0.0, yy = 0.0, zz = 0.0, rr, r;
                            for( int ii=0; ii<2; ii++ ) {
                                xx += pp[ii].x;
                                yy += pp[ii].y;
                                zz += pp[ii].z;
                            }
                            xx /= 2.0;
                            yy /= 2.0;
                            zz /= 2.0;
                            rr = xx * xx + yy * yy;
                            r = rr + zz * zz;

                            if( zz < 2.0001 /* || xx < 0.0001 */|| yy < 0.0001 || rr < 1.001 * 1.001 || r < 25.0 || r > 5.98 * 5.98 ){
                                if( zz < 2.0001 || rr < 1.001 * 1.001 || yy < 0.0001 ){
                                    go = CreateLine3D(pp);
                                    EGWithMaskChangeAttributes(WIDTH_MASK | STYLE_MASK | COLOR_MASK | LAYER_MASK, go);
                                    EGAttachObject(go, ( EObjectP ) this);
                                    EMAddGraphicsToModel(ESIModel(), go);
                                }
                            }
#endif
#else
                            go =  CreateLine3D(pp);
                            EGWithMaskChangeAttributes(WIDTH_MASK | STYLE_MASK | COLOR_MASK | LAYER_MASK, go);
                            EGAttachObject(go, ( EObjectP ) this);
                            EMAddGraphicsToModel(ESIModel(), go);
#endif
                        }

                        if ( i == 1 && j == nseg ) {
                            pp [ 0 ].x = p [ 3 ].x;
                            pp [ 0 ].y = p [ 3 ].y;
                            pp [ 0 ].z = p [ 3 ].z;

                            pp [ 1 ].x = p [ 7 ].x;
                            pp [ 1 ].y = p [ 7 ].y;
                            pp [ 1 ].z = p [ 7 ].z;

#ifdef DRAW_VISIBLE_CONTOUR
#ifdef SPHERE_WITH_HOLE
                            double xx = 0.0, yy = 0.0, zz = 0.0, rr, r;
                            for(int ii=0;ii<2;ii++){
                                xx += pp[ii].x;
                                yy += pp[ii].y;
                                zz += pp[ii].z;
                            }
                            xx /= 2.0;
                            yy /= 2.0;
                            zz /= 2.0;
                            rr = xx * xx + yy * yy;
                            r = rr + zz * zz;

                            if(zz < 2.0001 /* || xx < 0.0001 */|| yy < 0.0001 || rr < 1.001 * 1.001 || r < 25.0 || r > 5.98 * 5.98){
                                if(zz < 2.0001 || rr < 1.001 * 1.001 || yy < 0.0001){
                                    go =  CreateLine3D(pp);
                                    EGWithMaskChangeAttributes(WIDTH_MASK | STYLE_MASK | COLOR_MASK | LAYER_MASK, go);
                                    EGAttachObject(go, ( EObjectP ) this);
                                    EMAddGraphicsToModel(ESIModel(), go);
                                }
                            }
#endif
#else
                            go =  CreateLine3D(pp);
                            EGWithMaskChangeAttributes(WIDTH_MASK | STYLE_MASK | COLOR_MASK | LAYER_MASK, go);
                            EGAttachObject(go, ( EObjectP ) this);
                            EMAddGraphicsToModel(ESIModel(), go);
#endif
                        }

                        if ( i == nseg && j == 1 ) {
                            pp [ 0 ].x = p [ 1 ].x;
                            pp [ 0 ].y = p [ 1 ].y;
                            pp [ 0 ].z = p [ 1 ].z;

                            pp [ 1 ].x = p [ 5 ].x;
                            pp [ 1 ].y = p [ 5 ].y;
                            pp [ 1 ].z = p [ 5 ].z;

#ifdef DRAW_VISIBLE_CONTOUR
#ifdef SPHERE_WITH_HOLE
                            double xx = 0.0, yy = 0.0, zz = 0.0, rr, r;
                            for(int ii=0;ii<2;ii++){
                                xx += pp[ii].x;
                                yy += pp[ii].y;
                                zz += pp[ii].z;
                            }
                            xx /= 2.0;
                            yy /= 2.0;
                            zz /= 2.0;
                            rr = xx * xx + yy * yy;
                            r = rr + zz * zz;

                            if(zz < 2.0001 /* || xx < 0.0001 */|| yy < 0.0001 || rr < 1.001 * 1.001 || r < 25.0 || r > 5.98 * 5.98){
                                if(zz < 2.0001 || rr < 1.001 * 1.001 || yy < 0.0001){
                                    go =  CreateLine3D(pp);
                                    EGWithMaskChangeAttributes(WIDTH_MASK | STYLE_MASK | COLOR_MASK | LAYER_MASK, go);
                                    EGAttachObject(go, ( EObjectP ) this);
                                    EMAddGraphicsToModel(ESIModel(), go);
                                }
                            }
#endif
#else
                            go =  CreateLine3D(pp);
                            EGWithMaskChangeAttributes(WIDTH_MASK | STYLE_MASK | COLOR_MASK | LAYER_MASK, go);
                            EGAttachObject(go, ( EObjectP ) this);
                            EMAddGraphicsToModel(ESIModel(), go);
#endif
                        }

                        if ( i == nseg && j == nseg ) {
                            pp [ 0 ].x = p [ 2 ].x;
                            pp [ 0 ].y = p [ 2 ].y;
                            pp [ 0 ].z = p [ 2 ].z;

                            pp [ 1 ].x = p [ 6 ].x;
                            pp [ 1 ].y = p [ 6 ].y;
                            pp [ 1 ].z = p [ 6 ].z;

#ifdef DRAW_VISIBLE_CONTOUR
#ifdef SPHERE_WITH_HOLE
                            double xx = 0.0, yy = 0.0, zz = 0.0, rr, r;
                            for(int ii=0;ii<2;ii++){
                                xx += pp[ii].x;
                                yy += pp[ii].y;
                                zz += pp[ii].z;
                            }
                            xx /= 2.0;
                            yy /= 2.0;
                            zz /= 2.0;
                            rr = xx * xx + yy * yy;
                            r = rr + zz * zz;

                            if(zz < 2.0001 /* || xx < 0.0001 */|| yy < 0.0001 || rr < 1.001 * 1.001 || r < 25.0 || r > 5.98 * 5.98){
                                if(zz < 2.0001 || rr < 1.001 * 1.001 || yy < 0.0001){
                                    go =  CreateLine3D(pp);
                                    EGWithMaskChangeAttributes(WIDTH_MASK | STYLE_MASK | COLOR_MASK | LAYER_MASK, go);
                                    EGAttachObject(go, ( EObjectP ) this);
                                    EMAddGraphicsToModel(ESIModel(), go);
                                }
                            }
#endif
#else
                            go =  CreateLine3D(pp);
                            EGWithMaskChangeAttributes(WIDTH_MASK | STYLE_MASK | COLOR_MASK | LAYER_MASK, go);
                            EGAttachObject(go, ( EObjectP ) this);
                            EMAddGraphicsToModel(ESIModel(), go);
#endif
                        }

                        if ( j == 1 && k == 1 ) {
                            pp [ 0 ].x = p [ 0 ].x;
                            pp [ 0 ].y = p [ 0 ].y;
                            pp [ 0 ].z = p [ 0 ].z;

                            pp [ 1 ].x = p [ 1 ].x;
                            pp [ 1 ].y = p [ 1 ].y;
                            pp [ 1 ].z = p [ 1 ].z;

#ifdef DRAW_VISIBLE_CONTOUR
#ifdef SPHERE_WITH_HOLE
                            double xx = 0.0, yy = 0.0, zz = 0.0, rr, r;
                            for(int ii=0;ii<2;ii++){
                                xx += pp[ii].x;
                                yy += pp[ii].y;
                                zz += pp[ii].z;
                            }
                            xx /= 2.0;
                            yy /= 2.0;
                            zz /= 2.0;
                            rr = xx * xx + yy * yy;
                            r = rr + zz * zz;

                            if(zz < 2.0001 /* || xx < 0.0001 */|| yy < 0.0001 || rr < 1.001 * 1.001 || r < 25.0 || r > 5.98 * 5.98){
                                go =  CreateLine3D(pp);
                                EGWithMaskChangeAttributes(WIDTH_MASK | STYLE_MASK | COLOR_MASK | LAYER_MASK, go);
                                EGAttachObject(go, ( EObjectP ) this);
                                EMAddGraphicsToModel(ESIModel(), go);
                            }
#endif
#else
                           go =  CreateLine3D(pp);
                            EGWithMaskChangeAttributes(WIDTH_MASK | STYLE_MASK | COLOR_MASK | LAYER_MASK, go);
                            EGAttachObject(go, ( EObjectP ) this);
                            EMAddGraphicsToModel(ESIModel(), go);
#endif
                        }

                        if ( j == 1 && k == nseg ) {
                            pp [ 0 ].x = p [ 4 ].x;
                            pp [ 0 ].y = p [ 4 ].y;
                            pp [ 0 ].z = p [ 4 ].z;

                            pp [ 1 ].x = p [ 5 ].x;
                            pp [ 1 ].y = p [ 5 ].y;
                            pp [ 1 ].z = p [ 5 ].z;

#ifdef DRAW_VISIBLE_CONTOUR
#ifdef SPHERE_WITH_HOLE
                            double xx = 0.0, yy = 0.0, zz = 0.0, rr, r;
                            for(int ii=0;ii<2;ii++){
                                xx += pp[ii].x;
                                yy += pp[ii].y;
                                zz += pp[ii].z;
                            }
                            xx /= 2.0;
                            yy /= 2.0;
                            zz /= 2.0;
                            rr = xx * xx + yy * yy;
                            r = rr + zz * zz;

                            if(zz < 2.0001 /* || xx < 0.0001 */|| yy < 0.0001 || rr < 1.001 * 1.001 || r < 25.0 || r > 5.98 * 5.98){
                                if(yy < 1.5 || zz < 2.0001){
                                    go =  CreateLine3D(pp);
                                    EGWithMaskChangeAttributes(WIDTH_MASK | STYLE_MASK | COLOR_MASK | LAYER_MASK, go);
                                    EGAttachObject(go, ( EObjectP ) this);
                                    EMAddGraphicsToModel(ESIModel(), go);
                                }
                            }
#endif
#else
                            go =  CreateLine3D(pp);
                            EGWithMaskChangeAttributes(WIDTH_MASK | STYLE_MASK | COLOR_MASK | LAYER_MASK, go);
                            EGAttachObject(go, ( EObjectP ) this);
                            EMAddGraphicsToModel(ESIModel(), go);
#endif
                        }

                        if ( j == nseg && k == 1 ) {
                            pp [ 0 ].x = p [ 3 ].x;
                            pp [ 0 ].y = p [ 3 ].y;
                            pp [ 0 ].z = p [ 3 ].z;

                            pp [ 1 ].x = p [ 2 ].x;
                            pp [ 1 ].y = p [ 2 ].y;
                            pp [ 1 ].z = p [ 2 ].z;

#ifdef DRAW_VISIBLE_CONTOUR
#ifdef SPHERE_WITH_HOLE
                            double xx = 0.0, yy = 0.0, zz = 0.0, rr, r;
                            for(int ii=0;ii<2;ii++){
                                xx += pp[ii].x;
                                yy += pp[ii].y;
                                zz += pp[ii].z;
                            }
                            xx /= 2.0;
                            yy /= 2.0;
                            zz /= 2.0;
                            rr = xx * xx + yy * yy;
                            r = rr + zz * zz;

                            if(zz < 2.0001 /* || xx < 0.0001 */|| yy < 0.0001 || rr < 1.001 * 1.001 || r < 25.0 || r > 5.98 * 5.98){
                                go =  CreateLine3D(pp);
                                EGWithMaskChangeAttributes(WIDTH_MASK | STYLE_MASK | COLOR_MASK | LAYER_MASK, go);
                                EGAttachObject(go, ( EObjectP ) this);
                                EMAddGraphicsToModel(ESIModel(), go);
                            }
#endif
#else
                            go =  CreateLine3D(pp);
                            EGWithMaskChangeAttributes(WIDTH_MASK | STYLE_MASK | COLOR_MASK | LAYER_MASK, go);
                            EGAttachObject(go, ( EObjectP ) this);
                            EMAddGraphicsToModel(ESIModel(), go);
#endif
                        }

                        if ( j == nseg && k == nseg ) {
                            pp [ 0 ].x = p [ 7 ].x;
                            pp [ 0 ].y = p [ 7 ].y;
                            pp [ 0 ].z = p [ 7 ].z;

                            pp [ 1 ].x = p [ 6 ].x;
                            pp [ 1 ].y = p [ 6 ].y;
                            pp [ 1 ].z = p [ 6 ].z;

#ifdef DRAW_VISIBLE_CONTOUR
#ifdef SPHERE_WITH_HOLE
                            double xx = 0.0, yy = 0.0, zz = 0.0, rr, r;
                            for(int ii=0;ii<2;ii++){
                                xx += pp[ii].x;
                                yy += pp[ii].y;
                                zz += pp[ii].z;
                            }
                            xx /= 2.0;
                            yy /= 2.0;
                            zz /= 2.0;
                            rr = xx * xx + yy * yy;
                            r = rr + zz * zz;

                            if(zz < 2.0001 /* || xx < 0.0001 */|| yy < 0.0001 || rr < 1.001 * 1.001 || r < 25.0 || r > 5.98 * 5.98){
                                if(yy < 1.5 || zz < 2.0001){
                                    go =  CreateLine3D(pp);
                                    EGWithMaskChangeAttributes(WIDTH_MASK | STYLE_MASK | COLOR_MASK | LAYER_MASK, go);
                                    EGAttachObject(go, ( EObjectP ) this);
                                    EMAddGraphicsToModel(ESIModel(), go);
                                }
                            }
#endif
#else
                            go =  CreateLine3D(pp);
                            EGWithMaskChangeAttributes(WIDTH_MASK | STYLE_MASK | COLOR_MASK | LAYER_MASK, go);
                            EGAttachObject(go, ( EObjectP ) this);
                            EMAddGraphicsToModel(ESIModel(), go);
#endif
                        }

                        if ( k == 1 && i == 1 ) {
                            pp [ 0 ].x = p [ 0 ].x;
                            pp [ 0 ].y = p [ 0 ].y;
                            pp [ 0 ].z = p [ 0 ].z;

                            pp [ 1 ].x = p [ 3 ].x;
                            pp [ 1 ].y = p [ 3 ].y;
                            pp [ 1 ].z = p [ 3 ].z;

#ifdef DRAW_VISIBLE_CONTOUR
#ifdef SPHERE_WITH_HOLE
                            double xx = 0.0, yy = 0.0, zz = 0.0, rr, r;
                            for(int ii=0;ii<2;ii++){
                                xx += pp[ii].x;
                                yy += pp[ii].y;
                                zz += pp[ii].z;
                            }
                            xx /= 2.0;
                            yy /= 2.0;
                            zz /= 2.0;
                            rr = xx * xx + yy * yy;
                            r = rr + zz * zz;

                            if(zz < 2.0001 /* || xx < 0.0001 */|| yy < 0.0001 || rr < 1.001 * 1.001 || r < 25.0 || r > 5.98 * 5.98){
                                go =  CreateLine3D(pp);
                                EGWithMaskChangeAttributes(WIDTH_MASK | STYLE_MASK | COLOR_MASK | LAYER_MASK, go);
                                EGAttachObject(go, ( EObjectP ) this);
                                EMAddGraphicsToModel(ESIModel(), go);
                            }
#endif
#else
                            go =  CreateLine3D(pp);
                            EGWithMaskChangeAttributes(WIDTH_MASK | STYLE_MASK | COLOR_MASK | LAYER_MASK, go);
                            EGAttachObject(go, ( EObjectP ) this);
                            EMAddGraphicsToModel(ESIModel(), go);
#endif
                        }

                        if ( k == 1 && i == nseg ) {
                            pp [ 0 ].x = p [ 1 ].x;
                            pp [ 0 ].y = p [ 1 ].y;
                            pp [ 0 ].z = p [ 1 ].z;

                            pp [ 1 ].x = p [ 2 ].x;
                            pp [ 1 ].y = p [ 2 ].y;
                            pp [ 1 ].z = p [ 2 ].z;

#ifdef DRAW_VISIBLE_CONTOUR
#ifdef SPHERE_WITH_HOLE
                            double xx = 0.0, yy = 0.0, zz = 0.0, rr, r;
                            for(int ii=0;ii<2;ii++){
                                xx += pp[ii].x;
                                yy += pp[ii].y;
                                zz += pp[ii].z;
                            }
                            xx /= 2.0;
                            yy /= 2.0;
                            zz /= 2.0;
                            rr = xx * xx + yy * yy;
                            r = rr + zz * zz;

                            if(zz < 2.0001 /* || xx < 0.0001 */|| yy < 0.0001 || rr < 1.001 * 1.001 || r < 25.0 || r > 5.98 * 5.98){
                                go =  CreateLine3D(pp);
                                EGWithMaskChangeAttributes(WIDTH_MASK | STYLE_MASK | COLOR_MASK | LAYER_MASK, go);
                                EGAttachObject(go, ( EObjectP ) this);
                                EMAddGraphicsToModel(ESIModel(), go);
                            }
#endif
#else
                            go =  CreateLine3D(pp);
                            EGWithMaskChangeAttributes(WIDTH_MASK | STYLE_MASK | COLOR_MASK | LAYER_MASK, go);
                            EGAttachObject(go, ( EObjectP ) this);
                            EMAddGraphicsToModel(ESIModel(), go);
#endif
                        }

                        if ( k == nseg && i == 1 ) {
                            pp [ 0 ].x = p [ 4 ].x;
                            pp [ 0 ].y = p [ 4 ].y;
                            pp [ 0 ].z = p [ 4 ].z;

                            pp [ 1 ].x = p [ 7 ].x;
                            pp [ 1 ].y = p [ 7 ].y;
                            pp [ 1 ].z = p [ 7 ].z;

#ifdef DRAW_VISIBLE_CONTOUR
#ifdef SPHERE_WITH_HOLE
                            double xx = 0.0, yy = 0.0, zz = 0.0, rr, r;
                            for(int ii=0;ii<2;ii++){
                                xx += pp[ii].x;
                                yy += pp[ii].y;
                                zz += pp[ii].z;
                            }
                            xx /= 2.0;
                            yy /= 2.0;
                            zz /= 2.0;
                            rr = xx * xx + yy * yy;
                            r = rr + zz * zz;

                            if(zz < 2.0001 /* || xx < 0.0001 */|| yy < 0.0001 || rr < 1.001 * 1.001 || r < 25.0 || r > 5.98 * 5.98){
                                if(yy < 0.0001){
                                    go =  CreateLine3D(pp);
                                    EGWithMaskChangeAttributes(WIDTH_MASK | STYLE_MASK | COLOR_MASK | LAYER_MASK, go);
                                    EGAttachObject(go, ( EObjectP ) this);
                                    EMAddGraphicsToModel(ESIModel(), go);
                                }
                            }
#endif
#else
                            go =  CreateLine3D(pp);
                            EGWithMaskChangeAttributes(WIDTH_MASK | STYLE_MASK | COLOR_MASK | LAYER_MASK, go);
                            EGAttachObject(go, ( EObjectP ) this);
                            EMAddGraphicsToModel(ESIModel(), go);
#endif
                        }

                        if ( k == nseg && i == nseg ) {
                            pp [ 0 ].x = p [ 5 ].x;
                            pp [ 0 ].y = p [ 5 ].y;
                            pp [ 0 ].z = p [ 5 ].z;

                            pp [ 1 ].x = p [ 6 ].x;
                            pp [ 1 ].y = p [ 6 ].y;
                            pp [ 1 ].z = p [ 6 ].z;

#ifdef DRAW_VISIBLE_CONTOUR
#ifdef SPHERE_WITH_HOLE
                            double xx = 0.0, yy = 0.0, zz = 0.0, rr, r;
                            for(int ii=0;ii<2;ii++){
                                xx += pp[ii].x;
                                yy += pp[ii].y;
                                zz += pp[ii].z;
                            }
                            xx /= 2.0;
                            yy /= 2.0;
                            zz /= 2.0;
                            rr = xx * xx + yy * yy;
                            r = rr + zz * zz;

                            if(zz < 2.0001 /* || xx < 0.0001 */|| yy < 0.0001 || rr < 1.001 * 1.001 || r < 25.0 || r > 5.98 * 5.98){
                                if(yy < 2.0001){
                                    go =  CreateLine3D(pp);
                                    EGWithMaskChangeAttributes(WIDTH_MASK | STYLE_MASK | COLOR_MASK | LAYER_MASK, go);
                                    EGAttachObject(go, ( EObjectP ) this);
                                    EMAddGraphicsToModel(ESIModel(), go);
                                }
                            }
#endif
#else
                            go =  CreateLine3D(pp);
                            EGWithMaskChangeAttributes(WIDTH_MASK | STYLE_MASK | COLOR_MASK | LAYER_MASK, go);
                            EGAttachObject(go, ( EObjectP ) this);
                            EMAddGraphicsToModel(ESIModel(), go);
#endif
                        }

 #else
                        go =  CreateHexahedron(p);
                        EGWithMaskChangeAttributes(WIDTH_MASK | FILL_MASK | COLOR_MASK | EDGE_COLOR_MASK | EDGE_FLAG_MASK | LAYER_MASK, go);
                        EGAttachObject(go, ( EObjectP ) this);
                        EMAddGraphicsToModel(ESIModel(), go);
 #endif
                    }
                }
            }
        }                 // end loop over knot spans (irules)

    } else   {
        OOFEM_ERROR2("drawRawGeometry: not implemented for nsd = %d", nsd);
    }
}


void drawIGAPatchDeformedGeometry(Element *elem, StructuralElementEvaluator *se, oofegGraphicContext &gc, UnknownType)
{
    WCRec p [ 8 ];
    GraphicObj *go;
    int i, j, k, m, n, nseg;
    FloatArray u;
    FloatMatrix N;
    FloatArray ur, d;
    IntArray lc;
    FEInterpolation *interp = elem->giveInterpolation();
    TimeStep *stepN = elem->giveDomain()->giveEngngModel()->giveCurrentStep();
    double defScale = gc.getDefScale();

    if ( !gc.testElementGraphicActivity(elem) ) {
        return;
    }

    EASValsSetLineWidth(OOFEG_DEFORMED_GEOMETRY_WIDTH);
    EASValsSetColor( gc.getDeformedElementColor() );
    EASValsSetEdgeColor( gc.getElementEdgeColor() );
    EASValsSetEdgeFlag(true);
    EASValsSetLayer(OOFEG_DEFORMED_GEOMETRY_LAYER);
    EASValsSetLineStyle(SOLID_STYLE);
    EASValsSetFillStyle(FILL_SOLID);
    EASValsSetLineWidth(0);

#ifdef DRAW_MESH
    //nseg = 8;
    nseg = 4;
#else
    nseg = 4;
#endif

    int numberOfIntegrationRules = elem->giveNumberOfIntegrationRules();
    const double *const* knotVector = interp->giveKnotVector();
    const IntArray *span;
    IntegrationRule *iRule;
    int ir, nsd = interp->giveNsd();

    se->computeVectorOf(EID_MomentumBalance, VM_Total, stepN, u);

    if ( nsd == 1 ) {
        FloatArray c [ 2 ], cg [ 2 ];
        double du;

        for ( j = 0; j < 2; j++ ) {
            c [ j ].resize(1);
            cg [ j ].resize(1);
        }

        // loop over individual integration rules (i.e., knot spans)
        for ( ir = 0; ir < numberOfIntegrationRules; ir++ ) {
            iRule = elem->giveIntegrationRule(ir);
            span = iRule->giveKnotSpan();
            // divide span locally to get finer geometry rep.
            du = ( knotVector [ 0 ] [ span->at(1) + 1 ] - knotVector [ 0 ] [ span->at(1) ] ) / nseg;
            for ( i = 1; i <= nseg; i++ ) {
                c [ 0 ].at(1) = knotVector [ 0 ] [ span->at(1) ] + du * ( i - 1 );
                c [ 1 ].at(1) = knotVector [ 0 ] [ span->at(1) ] + du * i;

                for ( k = 0; k < 2; k++ ) {
                    // create a dummy ip's
                    FloatArray *cc = new FloatArray (c[k]);   // constructor of gp does not make its own copy
                    GaussPoint gp(iRule, 999, cc, 1.0, _PlaneStress);

                    // compute displacements at gp
                    se->computeNMatrixAt(N, & gp);

                    // get local code numbers corresponding to ir
                    se->giveIntegrationElementLocalCodeNumbers(lc, elem, gp.giveIntegrationRule(), EID_MomentumBalance);
                    ur.resize( N.giveNumberOfColumns() );
                    for ( n = 1; n <= lc.giveSize(); n++ ) {
                        ur.at(n) = u.at( lc.at(n) );
                    }

                    // interpolate displacements
                    d.beProductOf(N, ur);

                    interp->local2global(cg [ k ], c [ k ], FEIIGAElementGeometryWrapper( elem, iRule->giveKnotSpan() ));
                    p [ k ].x = ( FPNum ) ( cg [ k ].at(1) + d.at(1) * defScale );
                    p [ k ].y = 0.;
                    p [ k ].z = 0.;
                }

                go =  CreateLine3D(p);
                EGWithMaskChangeAttributes(WIDTH_MASK | STYLE_MASK | COLOR_MASK | LAYER_MASK, go);
                EGAttachObject(go, ( EObjectP ) elem);
                EMAddGraphicsToModel(ESIModel(), go);
            }
        }                 // end loop over knot spans (irules)

    } else if ( nsd == 2 )      {
        FloatArray c [ 4 ], cg [ 4 ];
        double du, dv;

        for ( j = 0; j < 4; j++ ) {
            c [ j ].resize(2);
            cg [ j ].resize(2);
        }

        // loop over individual integration rules (i.e., knot spans)
        for ( ir = 0; ir < numberOfIntegrationRules; ir++ ) {
            iRule = elem->giveIntegrationRule(ir);
            span = iRule->giveKnotSpan();
            // divide span locally to get finer geometry rep.
            du = ( knotVector [ 0 ] [ span->at(1) + 1 ] - knotVector [ 0 ] [ span->at(1) ] ) / nseg;
            dv = ( knotVector [ 1 ] [ span->at(2) + 1 ] - knotVector [ 1 ] [ span->at(2) ] ) / nseg;
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
                        // create a dummy ip's
                        FloatArray *cc = new FloatArray (c[k]);   // constructor of gp does not make its own copy
                        GaussPoint gp(iRule, 999, cc, 1.0, _PlaneStress);

                        // compute displacements at gp
                        se->computeNMatrixAt(N, & gp);

                        // get local code numbers corresponding to ir
                        se->giveIntegrationElementLocalCodeNumbers(lc, elem, gp.giveIntegrationRule(), EID_MomentumBalance);
                        ur.resize( N.giveNumberOfColumns() );
                        for ( n = 1; n <= lc.giveSize(); n++ ) {
                            ur.at(n) = u.at( lc.at(n) );
                        }

                        // interpolate displacements
                        d.beProductOf(N, ur);

                        interp->local2global(cg [ k ], c [ k ], FEIIGAElementGeometryWrapper( elem, iRule->giveKnotSpan() ));
                        p [ k ].x = ( FPNum ) ( cg [ k ].at(1) + d.at(1) * defScale );
                        p [ k ].y = ( FPNum ) ( cg [ k ].at(2) + d.at(2) * defScale );
                        p [ k ].z = 0.;
                    }

                    go =  CreateQuad3D(p);
                    EGWithMaskChangeAttributes(WIDTH_MASK | FILL_MASK | COLOR_MASK | EDGE_COLOR_MASK | EDGE_FLAG_MASK | LAYER_MASK, go);
                    EGAttachObject(go, ( EObjectP ) elem);
                    EMAddGraphicsToModel(ESIModel(), go);
                }
            }
        }                 // end loop over knot spans (irules)

    } else if ( nsd == 3 )      {
        FloatArray c [ 8 ], cg [ 8 ];
        double du, dv, dt;

        for ( j = 0; j < 8; j++ ) {
            c [ j ].resize(3);
            cg [ j ].resize(3);
        }

        // loop over individual integration rules (i.e., knot spans)
        for ( ir = 0; ir < numberOfIntegrationRules; ir++ ) {
            iRule = elem->giveIntegrationRule(ir);
            span = iRule->giveKnotSpan();
            // divide span locally to get finer geometry rep.
            du = ( knotVector [ 0 ] [ span->at(1) + 1 ] - knotVector [ 0 ] [ span->at(1) ] ) / nseg;
            dv = ( knotVector [ 1 ] [ span->at(2) + 1 ] - knotVector [ 1 ] [ span->at(2) ] ) / nseg;
            dt = ( knotVector [ 2 ] [ span->at(3) + 1 ] - knotVector [ 2 ] [ span->at(3) ] ) / nseg;
            for ( i = 1; i <= nseg; i++ ) {
                for ( j = 1; j <= nseg; j++ ) {
                    for ( k = 1; k <= nseg; k++ ) {
                        c [ 0 ].at(1) = knotVector [ 0 ] [ span->at(1) ] + du * ( i - 1 );
                        c [ 0 ].at(2) = knotVector [ 1 ] [ span->at(2) ] + dv * ( j - 1 );
                        c [ 0 ].at(3) = knotVector [ 2 ] [ span->at(3) ] + dt * ( k - 1 );
                        c [ 1 ].at(1) = knotVector [ 0 ] [ span->at(1) ] + du * i;
                        c [ 1 ].at(2) = knotVector [ 1 ] [ span->at(2) ] + dv * ( j - 1 );
                        c [ 1 ].at(3) = knotVector [ 2 ] [ span->at(3) ] + dt * ( k - 1 );
                        c [ 2 ].at(1) = knotVector [ 0 ] [ span->at(1) ] + du * i;
                        c [ 2 ].at(2) = knotVector [ 1 ] [ span->at(2) ] + dv * j;
                        c [ 2 ].at(3) = knotVector [ 2 ] [ span->at(3) ] + dt * ( k - 1 );
                        c [ 3 ].at(1) = knotVector [ 0 ] [ span->at(1) ] + du * ( i - 1 );
                        c [ 3 ].at(2) = knotVector [ 1 ] [ span->at(2) ] + dv * j;
                        c [ 3 ].at(3) = knotVector [ 2 ] [ span->at(3) ] + dt * ( k - 1 );
                        c [ 4 ].at(1) = knotVector [ 0 ] [ span->at(1) ] + du * ( i - 1 );
                        c [ 4 ].at(2) = knotVector [ 1 ] [ span->at(2) ] + dv * ( j - 1 );
                        c [ 4 ].at(3) = knotVector [ 2 ] [ span->at(3) ] + dt * k;
                        c [ 5 ].at(1) = knotVector [ 0 ] [ span->at(1) ] + du * i;
                        c [ 5 ].at(2) = knotVector [ 1 ] [ span->at(2) ] + dv * ( j - 1 );
                        c [ 5 ].at(3) = knotVector [ 2 ] [ span->at(3) ] + dt * k;
                        c [ 6 ].at(1) = knotVector [ 0 ] [ span->at(1) ] + du * i;
                        c [ 6 ].at(2) = knotVector [ 1 ] [ span->at(2) ] + dv * j;
                        c [ 6 ].at(3) = knotVector [ 2 ] [ span->at(3) ] + dt * k;
                        c [ 7 ].at(1) = knotVector [ 0 ] [ span->at(1) ] + du * ( i - 1 );
                        c [ 7 ].at(2) = knotVector [ 1 ] [ span->at(2) ] + dv * j;
                        c [ 7 ].at(3) = knotVector [ 2 ] [ span->at(3) ] + dt * k;

                        for ( m = 0; m < 8; m++ ) {
                            // create a dummy ip's
                            FloatArray *cc = new FloatArray (c[m]);   // constructor of gp does not make its own copy
                            GaussPoint gp(iRule, 999, cc, 1.0, _3dMat);

                            // compute displacements at gp
                            se->computeNMatrixAt(N, & gp);

                            // get local code numbers corresponding to ir
                            se->giveIntegrationElementLocalCodeNumbers(lc, elem, gp.giveIntegrationRule(), EID_MomentumBalance);
                            ur.resize( N.giveNumberOfColumns() );
                            for ( n = 1; n <= lc.giveSize(); n++ ) {
                                ur.at(n) = u.at( lc.at(n) );
                            }

                            // interpolate displacements
                            d.beProductOf(N, ur);

                            interp->local2global(cg [ m ], c [ m ], FEIIGAElementGeometryWrapper( elem, iRule->giveKnotSpan() ));
                            p [ m ].x = ( FPNum ) ( cg [ m ].at(1) + d.at(1) * defScale );
                            p [ m ].y = ( FPNum ) ( cg [ m ].at(2) + d.at(2) * defScale );
                            p [ m ].z = ( FPNum ) ( cg [ m ].at(3) + d.at(3) * defScale );
                        }

                        go =  CreateHexahedron(p);
                        EGWithMaskChangeAttributes(WIDTH_MASK | FILL_MASK | COLOR_MASK | EDGE_COLOR_MASK | EDGE_FLAG_MASK | LAYER_MASK, go);
                        EGAttachObject(go, ( EObjectP ) elem);
                        EMAddGraphicsToModel(ESIModel(), go);
                    }
                }
            }
        }                 // end loop over knot spans (irules)

    } else   {
        OOFEM_ERROR2("drawDeformedGeometry: not implemented for nsd = %d", nsd);
    }
}




#endif
} // end namespace oofem
