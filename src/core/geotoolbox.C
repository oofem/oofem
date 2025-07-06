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

#include "geotoolbox.h"
#include "mathfem.h"

#ifdef __OOFEG
 #include "oofeggraphiccontext.h"
#endif

namespace oofem {
//#define GRAPH_DEBUG_PRINT
#define GRAPH_LENGTH_FRAC 1.e-4

int
Polygon :: testPoint(double x, double y) const
{
    bool oddNODES = false;

    Polygon :: PolygonEdgeIterator it(this);
    Vertex p1, p2;
    double x1, x2, y1, y2;
    while ( it.giveNext(p1, p2) ) {
        x1 = p1.coords(0);
        y1 = p1.coords(1);
        x2 = p2.coords(0);
        y2 = p2.coords(1);

        if ( ( ( y1 < y ) && ( y2 >= y ) ) ||
            ( ( y2 < y ) && ( y1 >= y ) ) ) {
            if ( x1 + ( y - y1 ) / ( y2 - y1 ) * ( x2 - x1 ) < x ) {
                oddNODES = !oddNODES;
            }
        }
    }

    return oddNODES;
}

double
Polygon :: computeVolume() const
{
    double area = 0.0;
    double x1, x2, x3, y1, y2, y3;
    Vertex p;
    Polygon :: PolygonVertexIterator it(this);

    it.giveNext(p);
    x1 = p.coords(0);
    y1 = p.coords(1);
    it.giveNext(p);
    x2 = p.coords(0);
    y2 = p.coords(1);
    while ( it.giveNext(p) ) {
        x3 = p.coords(0);
        y3 = p.coords(1);
        area += ( x2 * y3 + x1 * y2 + y1 * x3 - x2 * y1 - x3 * y2 - x1 * y3 );
        x2 = x3;
        y2 = y3;
    }

    return 0.5 * fabs(area);
}


double
Polygon :: pointDistance(double xp, double yp) const
{
    // computes distance of a given point from a closed polygon
    // distance is signed, for polygon with positive orientation
    // (anticlockwise vertex ordering) the positive distance is for a point
    // located outside

    Polygon :: PolygonEdgeIterator it(this);
    Vertex p1, p2;
    double x1, x2, y1, y2, d, l, tx, ty, t, nx, ny, dist = 0.0;
    bool init = true;

    while ( it.giveNext(p1, p2) ) {
        x1 = p1.coords(0);
        y1 = p1.coords(1);
        x2 = p2.coords(0);
        y2 = p2.coords(1);

        // first check start vertex (end vertex checked by next edge)
        d = sqrt( ( xp - x1 ) * ( xp - x1 ) + ( yp - y1 ) * ( yp - y1 ) );
        if ( init ) {
            dist = d;
            init = false;
        } else {
            dist = min(dist, d);
        }

        // edge unit tangent
        l = sqrt( ( x2 - x1 ) * ( x2 - x1 ) + ( y2 - y1 ) * ( y2 - y1 ) );
        tx = ( x2 - x1 ) / l;
        ty = ( y2 - y1 ) / l;
        // projection of position vector (xp-x1,yp-y1) onto edge
        t = ( xp - x1 ) * tx + ( yp - y1 ) * ty;
        // test if inside
        if ( ( t >= 0.0 ) && ( t <= l ) ) {
            // compute distance from edge
            nx = ty;
            ny = -tx;
            d = fabs( ( xp - x1 ) * nx + ( yp - y1 ) * ny );
            dist = min(dist, d);
        }
    }

    if ( this->testPoint(xp, yp) ) {
        return -1.0 * dist;
    } else {
        return dist;
    }
}



#ifdef __OOFEG
GraphicObj *
Polygon :: draw(oofegGraphicContext &gc, bool filled, int layer)
{
    GraphicObj *go;
    LIST ggroup = make_list();

    EASValsSetLayer(layer);
    //EASValsSetColor(gc.getElementColor());

    if ( filled ) {
        int count = 0;
        double xc = 0.0, yc = 0.0;
        Vertex p;
        WCRec r [ 3 ];
        Polygon :: PolygonVertexIterator it(this);

        while ( it.giveNext(p) ) {
            xc += p.coords(0);
            yc += p.coords(1);
            count++;
        }

        xc /= count;
        yc /= count;

        EASValsSetFillStyle(FILL_SOLID);
        r [ 0 ].x = xc;
        r [ 0 ].y = yc;
        r [ 0 ].z = 0.0;
        Polygon :: PolygonEdgeIterator it2(this);
        Vertex p1, p2;
        while ( it2.giveNext(p1, p2) ) {
            r [ 1 ].x = p1.coords(0);
            r [ 1 ].y = p1.coords(1);
            r [ 1 ].z = 0.0;
            r [ 2 ].x = p2.coords(0);
            r [ 2 ].y = p2.coords(1);
            r [ 2 ].z = 0.0;
            go =  CreateTriangle3D(r);
            add_to_tail(ggroup, go);
            //EMAddGraphicsToModel(ESIModel(), go);
            EGWithMaskChangeAttributes(COLOR_MASK | FILL_MASK | LAYER_MASK, go);
        }
    } else {
        WCRec p [ 2 ];
        Polygon :: PolygonEdgeIterator it(this);
        Vertex p1, p2;
        EASValsSetFillStyle(FILL_HOLLOW);
        while ( it.giveNext(p1, p2) ) {
            p [ 0 ].x = p1.coords(0);
            p [ 0 ].y = p1.coords(1);
            p [ 0 ].z = 0.0;
            p [ 1 ].x = p2.coords(0);
            p [ 1 ].y = p2.coords(1);
            p [ 1 ].z = 0.0;

            go = CreateLine3D(p);
            add_to_tail(ggroup, go);
            //EMAddGraphicsToModel(ESIModel(), go);
            EGWithMaskChangeAttributes(COLOR_MASK | FILL_MASK | LAYER_MASK, go);
        }
    }

    go = CreateGgroup(ggroup);
    EMAddGraphicsToModel(ESIModel(), go);
    //EGWithMaskChangeAttributes(FILL_MASK | LAYER_MASK, go);
    return go;
}
#endif




Graph :: ~Graph()
{
    this->clear();
}


void
Graph :: clear()
{
    node *next, *auxs = s, *auxc = c;
    if ( s ) {
        do {
            next = auxs->next;
            delete auxs;
            auxs = next;
        } while ( auxs != s );
    }

    if ( c ) {
        do {
            next = auxc->next;
            delete auxc;
            auxc = next;
        } while ( auxc != c );
    }

    s = c = NULL;
}


int
Graph :: testIfCoincident(node *p1, node *p2, node *q1, node *q2, double *alpha1, double *alpha2)
{
    /* return codes:
     * -1 ...  lines are not parallel
     * 0 .... lines do not coincide
     *
     * 1 .... lines are the same
     * 2 .... lines are the same (difffferent orientation)
     *
     * 10 ... p1,p2 inside q1,q2 (alpha1 is parametric coord of p1 on q1-q2, alpha2 is param.coord of p2 on q1-q2)
     * 20 ... q1,q2 inside p1,p2 (alpha1 is parametric coord of q1 on p1-p2, alpha2 is param.coord of q2 on p1-p2)
     *
     * 11 ... p1 inside q1,q2 (alpha1 is parametric coord of p1 on q1-q2)
     *        q1 inside p1,p1 (alpha2 is parametric coord of q1 on p1-p2)
     * 12 ... p1 inside q1,q2 (alpha1 is parametric coord of p1 on q1-q2)
     *        q2 inside p1,p1 (alpha2 is parametric coord of q2 on p1-p2)
     * 21 ... p2 inside q1,q2 (alpha1 is parametric coord of p2 on q1-q2)
     *        q1 inside p1,p1 (alpha2 is parametric coord of q1 on p1-p2)
     * 22 ... p2 inside q1,q2 (alpha1 is parametric coord of p2 on q1-q2)
     *        q2 inside p1,p1 (alpha2 is parametric coord of q2 on p1-p2)
     *
     *
     * 1102 ... common part is p1=q1, q2 (alpha1 is param coord of q2 on p1-p2)
     * 1120 ... common part is p1=q1, p2 (alpha1 is param coord of p2 on q1-q2)
     * 2102 ... common part is p2=q1, q2 (alpha1 is param coord of q2 on p1-p2)
     * 2110 ... common part is p2=q1, p1 (alpha1 is param coord of p1 on q1-q2)
     * 1201 ... common part is p1=q2, q1 (alpha1 is param coord of q1 on p1-p2)
     * 1220 ... common part is p1=q2, p2 (alpha1 is param coord of p2 on q1-q2)
     * 2201 ... common part is p2=q2, q1 (alpha1 is param coord of q1 on p1-p2)
     * 2210 ... common part is p2=q2, p1 (alpha1 is param coord of p1 on q1-q2)
     *
     * 110 ... common part is only p1=q1, lines are parallel
     * 210 ... common part is only p2=q1, lines are parallel
     * 120 ... common part is only p1=q2, lines are parallel
     * 220 ... common part is only p2=q2, lines are parallel
     */


    //par = ((p2->x - p1->x)*(q2->y - q1->y) -
    //       (p2->y - p1->y)*(q2->x - q1->x));

    double par = ( ( p1->x - p2->x ) * ( q2->y - q1->y ) -
           ( p1->y - p2->y ) * ( q2->x - q1->x ) );

#ifdef GRAPH_DEBUG_PRINT
    fprintf(stderr, "p1 [%e %e] p2 [%e %e]    q1 [%e %e] q2 [%e %e]\n", p1->x, p1->y, p2->x, p2->y, q1->x, q1->y, q2->x, q2->y);
    fprintf(stderr, "par %e\n", par);
#endif
    //if (!par) return 0;                               /* parallel lines */
    if ( fabs(par) < GT_TEPS ) { /* parallel lines */
        /* test if identical part exist */
        double plength = dist(p1->x, p1->y, p2->x, p2->y);
        double qlength = dist(q1->x, q1->y, q2->x, q2->y);
        double ptx = ( p2->x - p1->x ) / plength;
        double pty = ( p2->y - p1->y ) / plength;
        double pnx = -pty;
        double pny =  ptx;
        double d1 = ( pnx * ( q1->x - p1->x ) + pny * ( q1->y - p1->y ) );
        double d2 = ( pnx * ( q2->x - p1->x ) + pny * ( q2->y - p1->y ) );
#ifdef GRAPH_DEBUG_PRINT
        fprintf(stderr, "dist (q-line, p-line)=%e\n", d1);
#endif
        if ( ( fabs( d1 / min(plength, qlength) ) >= GT_EPS ) && ( fabs( d2 / min(plength, qlength) ) >= GT_EPS ) && ( sgn(d1) * sgn(d2) > 0. ) ) {
            return 0;                                                                                                              // lines are parallel with nonzero distance
        }

        double qt1 = ( ptx * ( q1->x - p1->x ) + pty * ( q1->y - p1->y ) ) / plength;
        double qt2 = ( ptx * ( q2->x - p1->x ) + pty * ( q2->y - p1->y ) ) / plength;
#ifdef GRAPH_DEBUG_PRINT
        fprintf(stderr, "param of q1 on pline t1=%e\n", qt1);
        fprintf(stderr, "param of q2 on pline t2=%e\n", qt2);
#endif
        if ( ( qt1 < 0. ) && ( qt2 < 0. ) ) {
            return 0;                 // lines do not have common part
        }

        if ( ( qt1 > 1. ) && ( qt2 > 1. ) ) {
            return 0;                    // lines do not have common part
        }

        if ( ( fabs(qt1) < GT_EPS ) && ( fabs(1. - qt2) < GT_EPS ) ) { // lines are the same
            return 1;
        }

        if ( ( fabs(qt2) < GT_EPS ) && ( fabs(1. - qt1) < GT_EPS ) ) { // lines are the same but orientation is different
            return 2;
        }

        double qtx = ( q2->x - q1->x ) / qlength;
        double qty = ( q2->y - q1->y ) / qlength;
        double pt1 = ( qtx * ( p1->x - q1->x ) + qty * ( p1->y - q1->y ) ) / qlength;
        double pt2 = ( qtx * ( p2->x - q1->x ) + qty * ( p2->y - q1->y ) ) / qlength;

        // remains - partial overlap
        if ( ( qt1 < -GT_EPS ) && ( qt2 > ( 1. + GT_EPS ) ) ) { // p1,p2 inside q1,q2
            * alpha1 = pt1;
            * alpha2 = pt2;
            return 10;
        } else if ( ( qt2 < -GT_EPS ) && ( qt1 > ( 1. + GT_EPS ) ) ) { // p1,p2 inside q1,q2
            * alpha1 = pt1;
            * alpha2 = pt2;
            return 10;
        } else if ( ( qt1 > GT_EPS ) && ( qt1 < ( 1. - GT_EPS ) ) && ( qt2 > GT_EPS ) && ( qt2 < ( 1. - GT_EPS ) ) ) { // q1,q2 inside p1,p2
            * alpha1 = qt1;
            * alpha2 = qt2;
            return 20;
        } else if ( fabs(qt1) < GT_EPS ) { // q1 concident with p1
            if ( ( qt2 > GT_EPS ) && ( qt2 < ( 1. - GT_EPS ) ) ) { // q2 inside p1,p2
                * alpha1 = qt2;
                return 1102;
            } else if ( ( pt2 > GT_EPS ) && ( pt2 < ( 1. - GT_EPS ) ) ) { // p2 inside q1,q2
                * alpha1 = pt2;
                return 1120;
            } else { // only q1=p1 but no common part
                return 110;
            }
        } else if ( fabs(1. - qt1) < GT_EPS ) { // q1 coincident with p2
            if ( ( qt2 > GT_EPS ) && ( qt2 < ( 1. - GT_EPS ) ) ) { // q2 inside p1,p2
                * alpha1 = qt2;
                return 2102;
            } else if ( ( pt1 > GT_EPS ) && ( pt1 < ( 1. - GT_EPS ) ) ) { // p1 inside q1,q2
                * alpha1 = pt1;
                return 2110;
            } else { // only q1=p2 but no common part
                return 210;
            }
        } else if ( fabs(qt2) < GT_EPS ) { // q2 coincident with p1
            if ( ( qt1 > GT_EPS ) && ( qt1 < ( 1. - GT_EPS ) ) ) { // q1 inside p1,p2
                * alpha1 = qt1;
                return 1201;
            } else if ( ( pt2 > GT_EPS ) && ( pt2 < ( 1. - GT_EPS ) ) ) { // p2 inside q1,q2
                * alpha1 = pt2;
                return 1220;
            } else { // only q2=p1 but no common part
                return 120;
            }
        } else if ( fabs(1. - qt2) < GT_EPS ) { // q2 coincident with p2
            if ( ( qt1 > GT_EPS ) && ( qt1 < ( 1. - GT_EPS ) ) ) { // q1 inside p1,p2
                * alpha1 = qt1;
                return 2201;
            } else if ( ( pt1 > GT_EPS ) && ( pt1 < ( 1. - GT_EPS ) ) ) { // p1 inside q1,q2
                * alpha1 = pt1;
                return 2210;
            } else { // only q2=p2 but no common part
                return 220;
            }
        } else if ( ( qt1 >= GT_EPS ) && ( qt1 <= ( 1. - GT_EPS ) ) ) { // q1 inside, q2 outside
            double qlength = dist(q1->x, q1->y, q2->x, q2->y);
            double qtx = ( q2->x - q1->x ) / qlength;
            double qty = ( q2->y - q1->y ) / qlength;
            double pt1 = ( qtx * ( p1->x - q1->x ) + qty * ( p1->y - q1->y ) ) / qlength;
            double pt2 = ( qtx * ( p2->x - q1->x ) + qty * ( p2->y - q1->y ) ) / qlength;
            if ( ( pt1 >= GT_EPS ) && ( pt1 <= ( 1. - GT_EPS ) ) ) { // p1 inside
                * alpha1 = pt1;
                * alpha2 = qt1;
                return 11;
            } else if ( ( pt2 >= GT_EPS ) && ( pt2 <= ( 1. - GT_EPS ) ) ) { // p2 inside
                * alpha1 = pt2;
                * alpha2 = qt1;
                return 21;
            } else if ( ( pt1 <= GT_EPS || ( 1. - pt1 ) <= GT_EPS ) && ( ( pt2 < 0. ) || ( pt2 > 1.0 ) ) ) { // p1 inside
                * alpha1 = pt1;
                * alpha2 = qt1;
                return 11;
            } else if ( ( pt2 <= GT_EPS || ( 1. - pt2 ) <= GT_EPS ) && ( ( pt1 < 0. ) || ( pt1 > 1.0 ) ) ) { // p2 inside
                * alpha1 = pt2;
                * alpha2 = qt1;
                return 21;
            } else {
                fprintf(stderr, "testIfConcident: unresolved case q1 inside, q2 outside, pt1(t=%e), pt22(t=%e)\n", pt1, pt2);
                return -1;
            }
        } else if ( ( qt2 >= GT_EPS ) && ( qt2 <= ( 1. - GT_EPS ) ) ) { // q2 inside, q1 outside
            double qlength = dist(q1->x, q1->y, q2->x, q2->y);
            double qtx = ( q2->x - q1->x ) / qlength;
            double qty = ( q2->y - q1->y ) / qlength;
            double pt1 = ( qtx * ( p1->x - q1->x ) + qty * ( p1->y - q1->y ) ) / qlength;
            double pt2 = ( qtx * ( p2->x - q1->x ) + qty * ( p2->y - q1->y ) ) / qlength;
            if ( ( pt1 >= GT_EPS ) && ( pt1 <= ( 1. - GT_EPS ) ) ) { // p1 inside
                * alpha1 = pt1;
                * alpha2 = qt2;
                return 12;
            } else if ( ( pt2 >= GT_EPS ) && ( pt2 <= ( 1. - GT_EPS ) ) ) { // p2 inside
                * alpha1 = pt2;
                * alpha2 = qt2;
                return 22;
            } else if ( ( pt1 <= GT_EPS || ( 1. - pt1 ) <= GT_EPS ) && ( ( pt2 < 0. ) || ( pt2 > 1.0 ) ) ) { // p1 inside
                * alpha1 = pt1;
                * alpha2 = qt2;
                return 12;
            } else if ( ( pt2 <= GT_EPS || ( 1. - pt2 ) <= GT_EPS ) && ( ( pt1 < 0. ) || ( pt1 > 1.0 ) ) ) { // p2 inside
                * alpha1 = pt2;
                * alpha2 = qt2;
                return 22;
            } else {
                fprintf(stderr, "testIfConcident: unresolved case q2 inside, q1 outside, pt1(t=%e), pt2(t=%e)\n", pt1, pt2);
                return -1;
            }
        } else {
            fprintf(stderr, "testIfConcident: unresolved case, q1(t=%e), q2(t=%e)\n", qt1, qt2);
        }
    }

    return -1;
}


int
Graph :: testIfIntersect(node *p1, node *p2, node *q1, node *q2,
                         double *alpha_p, double *alpha_q, double *xint, double *yint)
{
    /* return codes:
     * 0 - no intersection detected
     * 1 - regular intersection found (alpha_q, alpha_q, xint, yint returned)
     * -1 - unhandled case, internal error
     *
     * 11 - p1 q1 coincide, lines intersect there
     * 22 - p2 q2 coincide, lines intersect there
     * 12 - p1 q2 coincide, lines intersect there
     * 21 - p2 q1 coincide, lines intersect there
     *
     *
     * 110 - p1 being intersection on q1,q2 (alpha_q returned)
     * 120 - p2 being intersection on q1,q2 (alpha_q returned)
     * 101 - q1 being intersection on p1,p2 (alpha_p returned)
     * 102 - q2 being intersection on p1,p2 (alpha_p returned)
     */



    double x, y, tp, tq, par;

    par = ( ( p1->x - p2->x ) * ( q2->y - q1->y ) -
           ( p1->y - p2->y ) * ( q2->x - q1->x ) );

#ifdef GRAPH_DEBUG_PRINT
    fprintf(stderr, "p1 [%e %e] p2 [%e %e]    q1 [%e %e] q2 [%e %e]\n", p1->x, p1->y, p2->x, p2->y, q1->x, q1->y, q2->x, q2->y);
    fprintf(stderr, "par %e\n", par);
#endif
    //if (!par) return 0;                               /* parallel lines */
    if ( fabs(par) < GT_TEPS ) { /* parallel lines */
        fprintf(stderr, "testIfIntersect: Parallel lines detected....\n");
        fprintf(stderr, "p1 [%e %e] p2 [%e %e]    q1 [%e %e] q2 [%e %e]\n", p1->x, p1->y, p2->x, p2->y, q1->x, q1->y, q2->x, q2->y);
        fprintf(stderr, "par %e\n", par);
        return -1;
    } else {
        //tp = ((q1->x - p1->x)*(q2->y - q1->y) - (q1->y - p1->y)*(q2->x - q1->x))/par;
        //tq = ((p2->y - p1->y)*(q1->x - p1->x) - (p2->x - p1->x)*(q1->y - p1->y))/par;

        tp = ( ( p1->x - q1->x ) * ( q2->y - q1->y ) - ( q2->x - q1->x ) * ( p1->y - q1->y ) ) / par;
        tq = ( ( p1->x - p2->x ) * ( p1->y - q1->y ) - ( p1->x - q1->x ) * ( p1->y - p2->y ) ) / par;

        //if(tp<0.0 || tp>1.0 || tq<0.0 || tq>1.0) return 0;
        //if(tp<0.0 || tp>1.0 || tq<0.0 || tq>1.0) return 0;
        if ( ( tp < -GT_EPS ) || ( tp > ( 1. + GT_EPS ) ) || ( tq < -GT_EPS ) || ( tq > ( 1 + GT_EPS ) ) ) {
            return 0;
        } else if ( ( tp >= GT_EPS ) && ( tp <= ( 1. - GT_EPS ) ) && ( tq >= GT_EPS ) && ( tq <= ( 1. - GT_EPS ) ) ) {
            // regular intersection found
            x = p1->x + tp * ( p2->x - p1->x );
            y = p1->y + tp * ( p2->y - p1->y );

            * alpha_p = dist(p1->x, p1->y, x, y) / dist(p1->x, p1->y, p2->x, p2->y);
            * alpha_q = dist(q1->x, q1->y, x, y) / dist(q1->x, q1->y, q2->x, q2->y);
            * xint = x;
            * yint = y;

            return 1;
        } else if ( ( ( fabs(tp) <= GT_EPS ) || ( fabs(1. - tp) <= GT_EPS ) )  && ( ( tq >= GT_EPS ) && ( tq <= ( 1. - GT_EPS ) ) ) ) {
            // one of p1 or p2 is intersection inside q1,q2 line
            if ( fabs(tp) <= GT_EPS ) {
                // p1 is on q1, q2
                * alpha_q = dist(q1->x, q1->y, p1->x, p1->y) / dist(q1->x, q1->y, q2->x, q2->y);
                return 110;
            } else {
                // p2 is on q1,q2
                * alpha_q = dist(q1->x, q1->y, p2->x, p2->y) / dist(q1->x, q1->y, q2->x, q2->y);
                return 120;
            }
        } else if ( ( ( fabs(tq) <= GT_EPS ) || ( fabs(1. - tq) <= GT_EPS ) )  && ( ( tp >= GT_EPS ) && ( tp <= ( 1. - GT_EPS ) ) ) ) {
            // one of q1 or q2 is intersection inside p1,p2 line
            if ( fabs(tq) <= GT_EPS ) {
                // q1 is on q1, q2
                * alpha_p = dist(p1->x, p1->y, q1->x, q1->y) / dist(p1->x, p1->y, p2->x, p2->y);
                return 101;
            } else {
                // q2 is on q1,q2
                * alpha_p = dist(p1->x, p1->y, q2->x, q2->y) / dist(p1->x, p1->y, p2->x, p2->y);
                return 102;
            }
        } else if ( ( ( fabs(tq) <= GT_EPS ) || ( fabs(1. - tq) <= GT_EPS ) )  && ( ( fabs(tp) <= GT_EPS ) || ( fabs(1. - tp) <= GT_EPS ) ) ) {
            // lines intersect at some of its vertices
            if ( fabs(tq) <= GT_EPS ) {
                if ( ( fabs(tp) <= GT_EPS ) ) {
                    // q1 p1 coincide, lines intersect there
                    return 11;
                } else {
                    // q1 p2 coincide, lines intersect there
                    return 21;
                }
            } else {
                if ( ( fabs(tp) <= GT_EPS ) ) {
                    // q2 p1 coincide, lines intersect there
                    return 12;
                } else {
                    // q2 p2 coincide, lines intersect there
                    return 22;
                }
            }
        } else {
            fprintf(stderr, "testIfIntersect: unhandled case [tp %e, tq %e, par %e]\n", tp, tq, par);
            return -1;
        }

    }
}




void
Graph :: clip(Polygon &result, const Polygon &a, const Polygon &b)
{
    // create lists s,c from a and b
    Polygon :: PolygonVertexIterator it( &a);
    Vertex v;
    node *nw;
    int ret;
    // s,c are input polygons
    node *auxs = NULL, *auxc = NULL, *is, *ic;
    node *p1, *p2, *q1, *q2;
    double xi, yi;
    double alpha_s, alpha_c, alpha1, alpha2;
    node *crt;

    result.clear();

    it.init(& a);
    while ( it.giveNext(v) ) {
        nw = this->createNode(v.coords(0), v.coords(1), 0, 0, 0, 0, NS_Vertex, 0, 0, 0);
        if ( c ) {
            c->prev->next = nw;
            nw->next = c;
            nw->prev = c->prev;
            c->prev = nw;
        } else {
            c = nw;
            c->next = c->prev = c;
        }
    }

    it.init(& b);
    while ( it.giveNext(v) ) {
        nw = this->createNode(v.coords(0), v.coords(1), 0, 0, 0, 0, NS_Vertex, 0, 0, 0);
        if ( s ) {
            s->prev->next = nw;
            nw->next = s;
            nw->prev = s->prev;
            s->prev = nw;
        } else {
            s = nw;
            s->next = s->prev = s;
        }
    }

    if ( !s || !c ) {
        return;
    }

    //auxs = last_node(s);
    //createNode(s->x, s->y, 0, auxs, 0, 0, 0, 0, 0, 0.);
    //auxc = last_node(c);
    //createNode(c->x, c->y, 0, auxc, 0, 0, 0, 0, 0, 0.);

    // restart:

    auxs = s;
    do {
        //for(auxs = s; auxs == s; auxs = auxs->next)
        if ( auxs->status != NS_Intersection ) {
            p1 = auxs;
            p2 = next_node(auxs->next);
            if ( testCollapsedEdge(p1, p2) ) {
                continue;
            }

            auxc = c;
            do {
                // for(auxc = c; auxc == c; auxc = auxc->next)
                if ( auxc->status != NS_Intersection ) {
                    q1 = auxc;
                    q2 = next_node(auxc->next);
                    if ( testCollapsedEdge(q1, q2) ) {
                        continue;
                    }


                    ret = testIfCoincident(p1, p2, q1, q2, & alpha1, & alpha2);
                    /* return codes:
                     * 0 .... lines are parallel, but do not coincide
                     *
                     * 1 .... lines are the same
                     * 2 .... lines are the same (difffferent orientation)
                     *
                     * 10 ... p1,p2 inside q1,q2 (alpha1 is parametric coord of p1 on q1-q2, alpha2 is param.coord of p2 on q1-q2)
                     * 20 ... q1,q2 inside p1,p2 (alpha1 is parametric coord of q1 on p1-p2, alpha2 is param.coord of q2 on p1-p2)
                     *
                     * 11 ... p1 inside q1,q2 (alpha1 is parametric coord of p1 on q1-q2)
                     *        q1 inside p1,p1 (alpha2 is parametric coord of q1 on p1-p2)
                     * 12 ... p1 inside q1,q2 (alpha1 is parametric coord of p1 on q1-q2)
                     *        q2 inside p1,p1 (alpha2 is parametric coord of q2 on p1-p2)
                     * 21 ... p2 inside q1,q2 (alpha1 is parametric coord of p2 on q1-q2)
                     *        q1 inside p1,p1 (alpha2 is parametric coord of q1 on p1-p2)
                     * 22 ... p2 inside q1,q2 (alpha1 is parametric coord of p2 on q1-q2)
                     *        q2 inside p1,p1 (alpha2 is parametric coord of q2 on p1-p2)
                     *
                     * 1102 ... common part is p1=q1, q2 (alpha1 is param coord of q2 on p1-p2)
                     * 1120 ... common part is p1=q1, p2 (alpha1 is param coord of p2 on q1-q2)
                     * 2102 ... common part is p2=q1, q2 (alpha1 is param coord of q2 on p1-p2)
                     * 2110 ... common part is p2=q1, p1 (alpha1 is param coord of p1 on q1-q2)
                     * 1201 ... common part is p1=q2, q1 (alpha1 is param coord of q1 on p1-p2)
                     * 1220 ... common part is p1=q2, p2 (alpha1 is param coord of p2 on q1-q2)
                     * 2201 ... common part is p2=q2, q1 (alpha1 is param coord of q1 on p1-p2)
                     * 2210 ... common part is p2=q2, p1 (alpha1 is param coord of p1 on q1-q2)
                     *
                     * 110 ... common part is only p1=q1, lines are parallel
                     * 210 ... common part is only p2=q1, lines are parallel
                     * 120 ... common part is only p1=q2, lines are parallel
                     * 220 ... common part is only p2=q2, lines are parallel
                     */
                    if ( ret == 0 ) { /*linear are parallel, but do not coincide*/
                    } else if ( ret == 1 ) { /* lines are coincident */
                        /*
                         * q1->status=p1->status=q2->status=p2->status=NS_IntersectionVertex;
                         */
                        merge2vertex(p1, q1);
                        merge2vertex(p2, q2);
                        //vertex2IntersectionVertex(p1,q1,q2);
                        //vertex2IntersectionVertex(p2,q1,q2);
                        //vertex2IntersectionVertex(q1,p1,p2);
                        //vertex2IntersectionVertex(q2,p1,p2);

                        //p1->neighbor=q1;  q1->neighbor=p1;
                        //p2->neighbor=q2;  q2->neighbor=p2;
                        p1->entry = -1;
                        q1->entry = -1;
                    } else if ( ret == 2 ) { /*lines are the same (difffferent orientation)*/
                        //q1->status=p1->status=q2->status=p2->status=NS_IntersectionVertex;
                        merge2vertex(p1, q2);
                        merge2vertex(p2, q1);
                        //vertex2IntersectionVertex(p1,q1,q2);
                        //vertex2IntersectionVertex(p2,q1,q2);
                        //vertex2IntersectionVertex(q1,p1,p2);
                        //vertex2IntersectionVertex(q2,p1,p2);

                        //p1->neighbor=q2;  q2->neighbor=p1;
                        //p2->neighbor=q1;  q1->neighbor=p2;
                        p1->entry = -1;
                        q1->entry = -1;
                    } else if ( ret == 10 ) { /* 10 ... p1,p2 inside q1,q2 (alpha1 is param. of p1 on q1-q2, alpha2 is param of p2 on q1-q2)*/
                        //p1->status = NS_IntersectionVertex;
                        if ( vertex2IntersectionVertex(p1, q1, q2) ) {
                            is = createNode(p1->x, p1->y, 0, 0, 0, 0, NS_Intersection, 0, 0, alpha1);
                            p1->neighbor = is;
                            is->neighbor = p1;
                            insert( is, auxc, next_node(auxc->next) );
                            if ( alpha1 < alpha2 ) {
                                is->entry = -1;
                            }

                            testNewIntersectionVertexEdgeCollapse(p1, q1, q2);
                        }

                        p1->entry = -1;

                        //p2->status = NS_IntersectionVertex;
                        if ( vertex2IntersectionVertex(p2, q1, q2) ) {
                            ic = createNode(p2->x, p2->y, 0, 0, 0, 0, NS_Intersection, 0, 0, alpha2);
                            p2->neighbor = ic;
                            ic->neighbor = p2;
                            insert( ic, auxc, next_node(auxc->next) );
                            if ( alpha1 > alpha2 ) {
                                ic->entry = -1;
                            }

                            testNewIntersectionVertexEdgeCollapse(p2, q1, q2);
                        }
                    } else if ( ret == 20 ) { /* 20...q1,q2 inside p1,p2 (alpha1 is parametric of q1 on p1-p2, alpha2 is param of q2 on p1-p2) */
                        //q1->status = NS_IntersectionVertex;
                        if ( vertex2IntersectionVertex(q1, p1, p2) ) {
                            is = createNode(q1->x, q1->y, 0, 0, 0, 0, NS_Intersection, 0, 0, alpha1);
                            q1->neighbor = is;
                            is->neighbor = q1;
                            insert( is, auxs, next_node(auxs->next) );
                            if ( alpha1 < alpha2 ) {
                                is->entry = -1;
                            }

                            testNewIntersectionVertexEdgeCollapse(q1, p1, p2);
                        }

                        q1->entry = -1;

                        //q2->status = NS_IntersectionVertex;
                        if ( vertex2IntersectionVertex(q2, p1, p2) ) {
                            ic = createNode(q2->x, q2->y, 0, 0, 0, 0, NS_Intersection, 0, 0, alpha2);
                            q2->neighbor = ic;
                            ic->neighbor = q2;
                            insert( ic, auxs, next_node(auxs->next) );
                            if ( alpha1 > alpha2 ) {
                                ic->entry = -1;
                            }

                            testNewIntersectionVertexEdgeCollapse(q2, p1, p2);
                        }
                    } else if ( ret == 11 ) {
                        /*
                         * 11 ... p1 inside q1,q2 (alpha1 is parametric coord of p1 on q1-q2)
                         *       q1 inside p1,p1 (alpha2 is parametric coord of q1 on p1-p2)
                         */
                        //p1->status=NS_IntersectionVertex;
                        if ( vertex2IntersectionVertex(p1, q1, q2) ) {
                            is = createNode(p1->x, p1->y, 0, 0, 0, 0, NS_Intersection, 0, 0, alpha1);
                            p1->neighbor = is;
                            is->neighbor = p1;
                            insert( is, auxc, next_node(auxc->next) );
                            testNewIntersectionVertexEdgeCollapse(p1, q1, q2);
                        }

                        p1->entry = -1;
                        //q1->status = NS_IntersectionVertex;
                        if ( vertex2IntersectionVertex(q1, p1, p2) ) {
                            is = createNode(q1->x, q1->y, 0, 0, 0, 0, NS_Intersection, 0, 0, alpha2);
                            q1->neighbor = is;
                            is->neighbor = q1;
                            insert( is, auxs, next_node(auxs->next) );
                            testNewIntersectionVertexEdgeCollapse(q1, p1, p2);
                        }

                        q1->entry = -1;
                    } else if ( ret == 12 ) {
                        /*
                         * 12 ... p1 inside q1,q2 (alpha1 is parametric coord of p1 on q1-q2)
                         *       q2 inside p1,p1 (alpha2 is parametric coord of q2 on p1-p2)
                         */
                        //p1->status=NS_IntersectionVertex;
                        if ( vertex2IntersectionVertex(p1, q1, q2) ) {
                            is = createNode(p1->x, p1->y, 0, 0, 0, 0, NS_Intersection, 0, 0, alpha1);
                            p1->neighbor = is;
                            is->neighbor = p1;
                            insert( is, auxc, next_node(auxc->next) );
                            testNewIntersectionVertexEdgeCollapse(p1, q1, q2);
                        }

                        p1->entry = p1->neighbor->entry = -1;

                        //q2->status = NS_IntersectionVertex;
                        if ( vertex2IntersectionVertex(q2, p1, p2) ) {
                            is = createNode(q2->x, q2->y, 0, 0, 0, 0, NS_Intersection, 0, 0, alpha2);
                            q2->neighbor = is;
                            is->neighbor = q2;
                            insert( is, auxs, next_node(auxs->next) );
                            testNewIntersectionVertexEdgeCollapse(q2, p1, p2);
                        }
                    } else if ( ret == 21 ) {
                        /*
                         * 21 ... p2 inside q1,q2 (alpha1 is parametric coord of p2 on q1-q2)
                         *       q1 inside p1,p1 (alpha2 is parametric coord of q1 on p1-p2)
                         */
                        //p2->status=NS_IntersectionVertex;
                        if ( vertex2IntersectionVertex(p2, q1, q2) ) {
                            is = createNode(p2->x, p2->y, 0, 0, 0, 0, NS_Intersection, 0, 0, alpha1);
                            p2->neighbor = is;
                            is->neighbor = p2;
                            insert( is, auxc, next_node(auxc->next) );
                            testNewIntersectionVertexEdgeCollapse(p2, q1, q2);
                        }

                        //q1->status = NS_IntersectionVertex;
                        if ( vertex2IntersectionVertex(q1, p1, p2) ) {
                            is = createNode(q1->x, q1->y, 0, 0, 0, 0, NS_Intersection, 0, 0, alpha2);
                            q1->neighbor = is;
                            is->neighbor = q1;
                            insert( is, auxs, next_node(auxs->next) );
                            testNewIntersectionVertexEdgeCollapse(q1, p1, p2);
                        }

                        q1->entry = q1->neighbor->entry = -1;
                    } else if ( ret == 22 ) {
                        /*
                         * 22 ... p2 inside q1,q2 (alpha1 is parametric coord of p2 on q1-q2)
                         *       q2 inside p1,p1 (alpha2 is parametric coord of q2 on p1-p2)
                         */
                        //p2->status=NS_IntersectionVertex;
                        if ( vertex2IntersectionVertex(p2, q1, q2) ) {
                            is = createNode(p2->x, p2->y, 0, 0, 0, 0, NS_Intersection, 0, 0, alpha1);
                            p2->neighbor = is;
                            is->neighbor = p2;
                            insert( is, auxc, next_node(auxc->next) );
                            testNewIntersectionVertexEdgeCollapse(p2, q1, q2);
                        }

                        p2->neighbor->entry = -1;
                        //q2->status = NS_IntersectionVertex;
                        if ( vertex2IntersectionVertex(q2, p1, p2) ) {
                            is = createNode(q2->x, q2->y, 0, 0, 0, 0, NS_Intersection, 0, 0, alpha2);
                            q2->neighbor = is;
                            is->neighbor = q2;
                            insert( is, auxs, next_node(auxs->next) );
                            testNewIntersectionVertexEdgeCollapse(q2, p1, p2);
                        }

                        q2->neighbor->entry = -1;
                    } else if ( ret == 1102 ) { /*1102 ... common part is p1=q1, q2 (alpha1 is param coord of q2 on p1-p2)*/
                        //p1->status = q1->status = NS_IntersectionVertex;
                        merge2vertex(p1, q1);
                        //vertex2IntersectionVertex(p1,q1,q2);
                        //vertex2IntersectionVertex(q1,p1,p2);
                        //p1->neighbor = q1; q1->neighbor = p1;
                        p1->entry = q1->entry = -1;

                        if ( vertex2IntersectionVertex(q2, p1, p2) ) {
                            //q2->status = NS_IntersectionVertex;
                            is = createNode(q2->x, q2->y, 0, 0, 0, 0, NS_Intersection, 0, 0, alpha1);
                            q2->neighbor = is;
                            is->neighbor = q2;
                            insert( is, auxs, next_node(auxs->next) );
                            testNewIntersectionVertexEdgeCollapse(q2, p1, p2);
                        }
                    } else if ( ret == 1120 ) { /*1120 ... common part is p1=q1, p2 (alpha1 is param coord of p2 on q1-q2)*/
                        //p1->status = q1->status = NS_IntersectionVertex;
                        merge2vertex(p1, q1);
                        //vertex2IntersectionVertex(p1,q1,q2);
                        //vertex2IntersectionVertex(q1,p1,p2);
                        //p1->neighbor = q1; q1->neighbor = p1;
                        p1->entry = q1->entry = -1;

                        //p2->status = NS_IntersectionVertex;
                        if ( vertex2IntersectionVertex(p2, q1, q2) ) {
                            is = createNode(p2->x, p2->y, 0, 0, 0, 0, NS_Intersection, 0, 0, alpha1);
                            p2->neighbor = is;
                            is->neighbor = p2;
                            insert( is, auxc, next_node(auxc->next) );
                            testNewIntersectionVertexEdgeCollapse(p2, q1, q2);
                        }
                    } else if ( ret == 2102 ) { /*2102 ... common part is p2=q1, q2 (alpha1 is param coord of q2 on p1-p2)*/
                        //p2->status = q1->status = NS_IntersectionVertex;
                        merge2vertex(p2, q1);
                        //vertex2IntersectionVertex(p2,q1,q2);
                        //vertex2IntersectionVertex(q1,p1,p2);
                        //p2->neighbor = q1; q1->neighbor = p2;
                        q1->entry = -1;

                        //q2->status = NS_IntersectionVertex;
                        if ( vertex2IntersectionVertex(q2, p1, p2) ) {
                            is = createNode(q2->x, q2->y, 0, 0, 0, 0, NS_Intersection, 0, 0, alpha1);
                            q2->neighbor = is;
                            is->neighbor = q2;
                            insert( is, auxs, next_node(auxs->next) );
                            testNewIntersectionVertexEdgeCollapse(q2, p1, p2);
                        }

                        q2->neighbor->entry = -1;
                    } else if ( ret == 2110 ) { /*2110 ... common part is p2=q1, p1 (alpha1 is param coord of p1 on q1-q2*/
                        //p2->status = q1->status = NS_IntersectionVertex;
                        merge2vertex(p2, q1);
                        //vertex2IntersectionVertex(p2,q1,q2);
                        //vertex2IntersectionVertex(q1,p1,p2);
                        //p2->neighbor = q1; q1->neighbor = p2;
                        q1->entry = -1;
                        p1->entry = -1;

                        //p1->status = NS_IntersectionVertex;
                        if ( vertex2IntersectionVertex(p1, q1, q2) ) {
                            is = createNode(p1->x, p1->y, 0, 0, 0, 0, NS_Intersection, 0, 0, alpha1);
                            p1->neighbor = is;
                            is->neighbor = p1;
                            insert( is, auxc, next_node(auxc->next) );
                            testNewIntersectionVertexEdgeCollapse(p1, q1, q2);
                        }
                    } else if ( ret == 1201 ) { /*1201 ... common part is p1=q2, q1 (alpha1 is param coord of q1 on p1-p2)*/
                        //p1->status = q2->status = NS_IntersectionVertex;
                        merge2vertex(p1, q2);
                        //vertex2IntersectionVertex(p1,q1,q2);
                        //vertex2IntersectionVertex(q2,p1,p2);
                        //p1->neighbor = q2; q2->neighbor = p1;
                        p1->entry = q1->entry = -1;

                        //q1->status = NS_IntersectionVertex;
                        if ( vertex2IntersectionVertex(q1, p1, p2) ) {
                            is = createNode(q1->x, q1->y, 0, 0, 0, 0, NS_Intersection, 0, 0, alpha1);
                            q1->neighbor = is;
                            is->neighbor = q1;
                            insert( is, auxs, next_node(auxs->next) );
                            testNewIntersectionVertexEdgeCollapse(q1, p1, p2);
                        }
                    } else if ( ret == 1220 ) { /*1220 ... common part is p1=q2, p2 (alpha1 is param coord of p2 on q1-q2)*/
                        //p1->status = q2->status = NS_IntersectionVertex;
                        merge2vertex(p1, q2);
                        //vertex2IntersectionVertex(p1,q1,q2);
                        //vertex2IntersectionVertex(q2,p1,p2);
                        //p1->neighbor = q2; q2->neighbor = p1;
                        p1->entry = -1;

                        //p2->status = NS_IntersectionVertex;
                        if ( vertex2IntersectionVertex(p2, q1, q2) ) {
                            is = createNode(p2->x, p2->y, 0, 0, 0, 0, NS_Intersection, 0, 0, alpha1);
                            p2->neighbor = is;
                            is->neighbor = p2;
                            insert( is, auxc, next_node(auxc->next) );
                            testNewIntersectionVertexEdgeCollapse(p2, q1, q2);
                        }

                        p2->neighbor->entry = -1;
                    } else if ( ret == 2201 ) { /*2201 ... common part is p2=q2, q1 (alpha1 is param coord of q1 on p1-p2)*/
                        //p2->status = q2->status = NS_IntersectionVertex;
                        merge2vertex(p2, q2);
                        //vertex2IntersectionVertex(p2,q1,q2);
                        //vertex2IntersectionVertex(q2,p1,p2);
                        //p2->neighbor = q2; q2->neighbor = p2;
                        q1->entry = -1;

                        //q1->status = NS_IntersectionVertex;
                        if ( vertex2IntersectionVertex(q1, p1, p2) ) {
                            is = createNode(q1->x, q1->y, 0, 0, 0, 0, NS_Intersection, 0, 0, alpha1);
                            q1->neighbor = is;
                            is->neighbor = q1;
                            insert( is, auxs, next_node(auxs->next) );
                            testNewIntersectionVertexEdgeCollapse(q1, p1, p2);
                        }

                        q1->neighbor->entry = -1;
                    } else if ( ret == 2210 ) { /*2210 ... common part is p2=q2, p1 (alpha1 is param coord of p1 on q1-q2)*/
                        //p2->status = q2->status = NS_IntersectionVertex;
                        merge2vertex(p2, q2);
                        //vertex2IntersectionVertex(p2,q1,q2);
                        //vertex2IntersectionVertex(q2,p1,p2);
                        //p2->neighbor = q2; q2->neighbor = p2;
                        p1->entry = -1;

                        //p1->status = NS_IntersectionVertex;
                        if ( vertex2IntersectionVertex(p1, q1, q2) ) {
                            is = createNode(p1->x, p1->y, 0, 0, 0, 0, NS_Intersection, 0, 0, alpha1);
                            p1->neighbor = is;
                            is->neighbor = p1;
                            insert( is, auxc, next_node(auxc->next) );
                            testNewIntersectionVertexEdgeCollapse(p1, q1, q2);
                        }

                        p1->neighbor->entry = -1;
                    } else if ( ret == 110 ) { /*common part is only p1=q1, lines are parallel*/
                        //p1->status=q1->status = NS_IntersectionVertex;
                        merge2vertex(p1, q1);
                        //vertex2IntersectionVertex(p1,q1,q2);
                        //vertex2IntersectionVertex(q1,p1,p2);
                        //p1->neighbor = q1; q1->neighbor = p1;
                    } else if ( ret == 210 ) { /*common part is only p2=q1, lines are parallel*/
                        //p2->status=q1->status = NS_IntersectionVertex;
                        merge2vertex(p2, q1);
                        //vertex2IntersectionVertex(p2,q1,q2);
                        //vertex2IntersectionVertex(q1,p1,p2);
                        //p2->neighbor = q1; q1->neighbor = p2;
                    } else if ( ret == 120 ) { /*common part is only p1=q2, lines are parallel*/
                        //p1->status=q2->status = NS_IntersectionVertex;
                        merge2vertex(p1, q2);
                        //vertex2IntersectionVertex(p1,q1,q2);
                        //vertex2IntersectionVertex(q2,p1,p2);
                        //p1->neighbor = q2; q2->neighbor = p1;
                    } else if ( ret == 220 ) { /*common part is only p2=q2, lines are parallel*/
                        //p2->status=q2->status = NS_IntersectionVertex;
                        merge2vertex(p2, q2);
                        //vertex2IntersectionVertex(p2,q1,q2);
                        //vertex2IntersectionVertex(q2,p1,p2);
                        //p2->neighbor = q2; q2->neighbor = p2;
                    } else if ( ret == -1 ) {
                        ret = testIfIntersect(auxs, next_node(auxs->next), auxc, next_node(auxc->next),
                                              & alpha_s, & alpha_c, & xi, & yi);
                        /* return codes:
                         * 0 - no intersection detected
                         * 1 - regular intersection found (alpha_q, alpha_q, xint, yint returned)
                         * -1 - unhandled case, internal error
                         *
                         * 11 - p1 q1 coincide, lines intersect there
                         * 22 - p2 q2 coincide, lines intersect there
                         * 12 - p1 q2 coincide, lines intersect there
                         * 21 - p2 q1 coincide, lines intersect there
                         *
                         *
                         * 110 - p1 being intersection on q1,q2 (alpha_q returned)
                         * 120 - p2 being intersection on q1,q2 (alpha_q returned)
                         * 101 - q1 being intersection on p1,p2 (alpha_p returned)
                         * 102 - q2 being intersection on p1,p2 (alpha_p returned)
                         */
                        if ( ret == 0 ) { // no intersection
                        } else if ( ret == 1 ) {
                            bool skip = false;
                            // 1 .... regular intersection
                            if ( ( p1->status == NS_IntersectionVertex ) || ( p2->status == NS_IntersectionVertex ) ||
                                ( q1->status == NS_IntersectionVertex ) || ( q2->status == NS_IntersectionVertex ) ) {
                                bool _p1 = false, _p2 = false, _q1 = false, _q2 = false;
                                int c1 = 0, c2 = 0;
                                // check if IV not associated to same line
                                if ( p1->status == NS_IntersectionVertex ) {
                                    if ( ( _p1 = belongs(p1, q1, q2) ) ) {
                                        c1++;
                                    }
                                }

                                if ( p2->status == NS_IntersectionVertex ) {
                                    if ( ( _p2 = belongs(p2, q1, q2) ) ) {
                                        c1++;
                                    }
                                }

                                if ( q1->status == NS_IntersectionVertex ) {
                                    if ( ( _q1 = belongs(q1, p1, p2) ) ) {
                                        c2++;
                                    }
                                }

                                if ( q2->status == NS_IntersectionVertex ) {
                                    if ( ( _q2 = belongs(q2, p1, p2) ) ) {
                                        c2++;
                                    }
                                }

                                if ( ( c1 == 1 ) || ( c2 == 1 ) ) {
                                    skip = true;
                                } else if ( ( c1 > 1 ) || ( c2 > 1 ) ) {
                                    // q1 q2 is already boundary edge on p1 p2 (or vice versa)
                                    skip = true;
                                    //this->printYourself();
                                    //printf ("c1 %d, c2 %d, p1x %e p1y %e, q1x %e, q1y %e\n", c1,c2,p1->x,p1->y,q1->x, q1->y);
                                    //THROW_GT_EXCEPTIONM ("Graph::clip unhandled case");
                                }
                            }

                            if ( !skip ) {
                                is = createNode(xi, yi, 0, 0, 0, 0, NS_Intersection, 0, 0, alpha_s);
                                ic = createNode(xi, yi, 0, 0, 0, 0, NS_Intersection, 0, 0, alpha_c);
                                is->neighbor = ic;
                                ic->neighbor = is;
                                insert( is, auxs, next_node(auxs->next) );
                                insert( ic, auxc, next_node(auxc->next) );
                            }
                        } else if ( ret == 11 ) { /*11 - p1 q1 coincide, lines intersect there*/
                            //p1->status=q1->status=NS_IntersectionVertex;
                            merge2vertex(p1, q1);
                            //vertex2IntersectionVertex(p1,q1,q2);
                            //vertex2IntersectionVertex(q1,p1,p2);
                            //p1->neighbor=q1; q1->neighbor=p1;
                        } else if ( ret == 22 ) { /*22 - p2 q2 coincide, lines intersect there*/
                            //p2->status=q2->status=NS_IntersectionVertex;
                            merge2vertex(p2, q2);
                            //vertex2IntersectionVertex(p2,q1,q2);
                            //vertex2IntersectionVertex(q2,p1,p2);
                            //p2->neighbor=q2; q2->neighbor=p2;
                        } else if ( ret == 12 ) { /*12 - p1 q2 coincide, lines intersect there*/
                            //p1->status=q2->status=NS_IntersectionVertex;
                            merge2vertex(p1, q2);
                            //vertex2IntersectionVertex(p1,q1,q2);
                            //vertex2IntersectionVertex(q2,p1,p2);
                            //p1->neighbor=q2; q2->neighbor=p1;
                        } else if ( ret == 21 ) { /*21 - p2 q1 coincide, lines intersect there*/
                            //p2->status=q1->status=NS_IntersectionVertex;
                            merge2vertex(p2, q1);
                            //vertex2IntersectionVertex(p2,q1,q2);
                            //vertex2IntersectionVertex(q1,p1,p2);
                            //p2->neighbor=q1; q1->neighbor=p2;
                        } else if ( ret == 110 ) { /*110 - p1 being intersection on q1,q2 (alpha_q returned)*/
                            //p1->status=NS_IntersectionVertex;
                            if ( vertex2IntersectionVertex(p1, q1, q2) ) {
                                is = createNode(p1->x, p1->y, 0, 0, 0, 0, NS_Intersection, 0, 0, alpha_c);
                                p1->neighbor = is;
                                is->neighbor = p1;
                                insert( is, auxc, next_node(auxc->next) );
                                testNewIntersectionVertexEdgeCollapse(p1, q1, q2);
                            }
                        } else if ( ret == 120 ) { /*120 - p2 being intersection on q1,q2 (alpha_q returned)*/
                            //p2->status=NS_IntersectionVertex;
                            if ( vertex2IntersectionVertex(p2, q1, q2) ) {
                                is = createNode(p2->x, p2->y, 0, 0, 0, 0, NS_Intersection, 0, 0, alpha_c);
                                p2->neighbor = is;
                                is->neighbor = p2;
                                insert( is, auxc, next_node(auxc->next) );
                                testNewIntersectionVertexEdgeCollapse(p2, q1, q2);
                            }
                        } else if ( ret == 101 ) { /* 101 - q1 being intersection on p1,p2 (alpha_p returned)*/
                            //q1->status=NS_IntersectionVertex;
                            if ( vertex2IntersectionVertex(q1, p1, p2) ) {
                                is = createNode(q1->x, q1->y, 0, 0, 0, 0, NS_Intersection, 0, 0, alpha_s);
                                q1->neighbor = is;
                                is->neighbor = q1;
                                insert( is, auxs, next_node(auxs->next) );
                                testNewIntersectionVertexEdgeCollapse(q1, p1, p2);
                            }
                        } else if ( ret == 102 ) { /* 102 - q2 being intersection on p1,p2 (alpha_p returned)*/
                            //q2->status=NS_IntersectionVertex;
                            if ( vertex2IntersectionVertex(q2, p1, p2) ) {
                                is = createNode(q2->x, q2->y, 0, 0, 0, 0, NS_Intersection, 0, 0, alpha_s);
                                q2->neighbor = is;
                                is->neighbor = q2;
                                insert( is, auxs, next_node(auxs->next) );
                                testNewIntersectionVertexEdgeCollapse(q2, p1, p2);
                            }
                        } else {
                            char buff [ 100 ];
                            sprintf(buff, "Graph::clip: unhandled return code (%d) from testIfIntersect service", ret);
                            THROW_GT_EXCEPTIONM(buff);
                        }
                    } else {
                        char buff [ 100 ];
                        sprintf(buff, "Graph::clip: unhandled return code (%d) from testIfCoincident service", ret);
                        THROW_GT_EXCEPTIONM(buff);
                    }
                }
            } while ( ( auxc = auxc->next ) != c );
        }
    } while ( ( auxs = auxs->next ) != s );

#ifdef GRAPH_DEBUG_PRINT
    /* print graph */
    printYourself();
    /* end of debug printing */
#endif

    // ================END OF WORK ==========================================

    // assign status (outer/inner) to graph intersection vertices
    int cnt = 0;
    bool identical = true;
    double xx, yy;
    //insidea = e = testPoint(c, s->x, s->y);
    auxs = s;
    do {
        if ( auxs->entry != -1 ) {
            identical = false;
        }

        if ( auxs->status == NS_Intersection ) {
            if ( auxs->entry != -1 ) { // boundary
                double xn, yn;
                if ( auxs->next->status == NS_IntersectionVertex ) {
                    xn = 0.5 * ( auxs->next->x + auxs->next->neighbor->x );
                    yn = 0.5 * ( auxs->next->y + auxs->next->neighbor->y );
                } else {
                    xn = auxs->next->x;
                    yn = auxs->next->y;
                }

                xx = ( auxs->x + xn ) * 0.5;
                yy = ( auxs->y + yn ) * 0.5;
                auxs->entry = testPoint(c, xx, yy);
                cnt++;
            }
        } else if ( auxs->status == NS_IntersectionVertex ) {
            if ( auxs->entry != -1 ) { //boundary
                double xn, yn;
                if ( auxs->next->status == NS_IntersectionVertex ) {
                    xn = 0.5 * ( auxs->next->x + auxs->next->neighbor->x );
                    yn = 0.5 * ( auxs->next->y + auxs->next->neighbor->y );
                } else {
                    xn = auxs->next->x;
                    yn = auxs->next->y;
                }

                xx = ( 0.5 * ( auxs->x + auxs->neighbor->x ) + xn ) * 0.5;
                yy = ( 0.5 * ( auxs->y + auxs->neighbor->y ) + yn ) * 0.5;
                auxs->entry = testPoint(c, xx, yy);
                cnt++;
            }
        }

        auxs = auxs->next;
    } while ( auxs != s );


    //insideb = e = testPoint(s, c->x, c->y);
    do {
        if ( auxc->status == NS_Intersection ) {
            if ( auxc->entry != -1 ) { // boundary
                double xn, yn;
                if ( auxc->next->status == NS_IntersectionVertex ) {
                    xn = 0.5 * ( auxc->next->x + auxc->next->neighbor->x );
                    yn = 0.5 * ( auxc->next->y + auxc->next->neighbor->y );
                } else {
                    xn = auxc->next->x;
                    yn = auxc->next->y;
                }

                xx = ( auxc->x + xn ) * 0.5;
                yy = ( auxc->y + yn ) * 0.5;
                auxc->entry = testPoint(s, xx, yy);
            }
        } else if ( auxc->status == NS_IntersectionVertex ) {
            if ( auxc->entry != -1 ) { //boundary
                double xn, yn;
                if ( auxc->next->status == NS_IntersectionVertex ) {
                    xn = 0.5 * ( auxc->next->x + auxc->next->neighbor->x );
                    yn = 0.5 * ( auxc->next->y + auxc->next->neighbor->y );
                } else {
                    xn = auxc->next->x;
                    yn = auxc->next->y;
                }

                xx = ( 0.5 * ( auxc->x + auxc->neighbor->x ) + xn ) * 0.5;
                yy = ( 0.5 * ( auxc->y + auxc->neighbor->y ) + yn ) * 0.5;
                auxc->entry = testPoint(s, xx, yy);
            }
        }

        auxc = auxc->next;
    } while ( auxc != c );


    if ( identical ) {
        auxs = s;
        do {
            v.setCoords(auxs->x, auxs->y);
            result.addVertex(v);
            auxs = auxs->next;
        } while ( auxs != s );
    } else if ( cnt == 0 ) {
        xx = ( s->x + s->next->x ) * 0.5;
        yy = ( s->y + s->next->y ) * 0.5;
        if ( testPoint(c, xx, yy) ) {
            // s fully inside c
            auxs = s;
            do {
                v.setCoords(auxs->x, auxs->y);
                result.addVertex(v);
                auxs = auxs->next;
            } while ( auxs != s );
        } else if ( testPoint( s, 0.5 * ( c->x + c->next->x ), 0.5 * ( c->y + c->next->y ) ) ) {
            // c fully in s
            auxc = c;
            do {
                v.setCoords(auxc->x, auxc->y);
                result.addVertex(v);
                auxc = auxc->next;
            } while ( auxc != c );
        }
    } else {
        /* loop over all entry points (non visited intersection) */
        //while ((crt = first(s)) != s)
        while ( ( crt = first(s) ) ) {
            //old = 0;
            /* loop over not yet visited intersections */
            for ( ; !crt->visited; crt = crt->neighbor ) {
                if ( crt->entry == 0 ) {
                    crt->visited = 1;
                    crt = crt->neighbor; // if not entry intersection use neighbor
                    //if (crt->entry ==0 || crt->visited ) break;
                    if ( crt->entry == 0 ) {
                        break;
                    }
                } else if ( crt->entry == -1 ) {
                    crt->visited = 1;
                    crt->neighbor->visited = 1;
                }

                crt->visited = 1;
                crt->neighbor->visited = 1;


                for ( ; ; ) { // loop until next intersection
                    //nw = create(crt->x, crt->y, old, 0, 0, 0, 0, 0, 0, 0.);
                    v.setCoords(crt->x, crt->y);
                    result.addVertex(v);
                    //old = new;
                    crt->visited = 1;
                    crt = crt->next;
                    //if((crt->status==NS_Intersection) || ((crt->status==NS_IntersectionVertex) && (crt->entry != -1)))
                    if ( ( crt->status == NS_Intersection ) || ( ( crt->status == NS_IntersectionVertex ) ) ) {
                        crt->visited = 1;
                        break;
                    }
                }

                //old->nextPoly = root;
                //root = old;
            }
        }
    }

#ifdef GRAPH_DEBUG_PRINT
    /* print graph */
    printf("Resulting Graph:\n");
    it.init(& result);
    while ( it.giveNext(v) ) {
        printf("Node [xy %e %e %e]\n", v.coords(0), v.coords(1), 0.0);
    }

#endif

    //printf ("Area=%e\n", result.computeVolume());
}



Graph :: node *
Graph :: prev_node(node *p)
{
    node *aux = p;
    while ( aux && ( aux->status == NS_Intersection ) ) {
        aux = aux->prev;
    }

    return aux;
}


Graph :: node *
Graph :: next_node(node *p)
{
    node *aux = p;
    while ( aux && ( aux->status == NS_Intersection ) ) {
        aux = aux->next;
    }

    return aux;
}

Graph :: node *
Graph :: last_node(node *p)
{
    node *aux = p;
    if ( aux ) {
        while ( aux->next ) {
            aux = aux->next;
        }
    }

    return aux;
}

Graph :: node *
Graph :: first(node *p)
{
    node *aux = p;
    /*
     * if (aux)
     * do aux=aux->next;
     * while(aux!=p && (aux->status==NS_Vertex || aux->status==NS_Intersection && aux->visited || aux->status==NS_IntersectionVertex && ((aux->entry == -1) || (aux->visited)) ));
     * return aux;
     */
    if ( aux ) {
        do {
            if ( ( aux->status == NS_Intersection && !aux->visited ) || ( aux->status == NS_IntersectionVertex && ( aux->entry != -1 ) && ( !aux->visited ) ) ) {
                return aux;
            }

            aux = aux->next;
        } while ( aux != p );
    }

    return NULL;
}


bool
Graph :: belongs(node *n, node *v1, node *v2)
// tests if node n (assoc to graph1) belongs to line v1,v2 of graph2
{
    if ( n->status == NS_Vertex ) {
        return false;
    } else if ( n->neighbor->status == NS_Intersection ) {  //Intersection or intersection vertex
        if ( ( next_node(n->neighbor) == v2 ) && ( prev_node(n->neighbor) == v1 ) ) {
            return true;
        }

        if ( ( prev_node(n->neighbor) == v2 ) && ( next_node(n->neighbor) == v1 ) ) {
            return true;
        }
    } else { // Vertex to vertex
        if ( n->neighbor == v1 || n->neighbor == v2 ) {
            return true;
        }
    }

    return false;
}

double
Graph :: dist(double x1, double y1, double x2, double y2)
{
    return sqrt( ( x1 - x2 ) * ( x1 - x2 ) + ( y1 - y2 ) * ( y1 - y2 ) );
}

Graph :: node *
Graph :: createNode(double x, double y, node *next, node *prev, node *nextPoly,
                    node *neighbor, nodeStatus st, int entry, int visited, double alpha)
{
    node *nw = new node;
    nw->x = x;
    nw->y = y;
    nw->next = next;
    nw->prev = prev;
    if ( prev ) {
        nw->prev->next = nw;
    }

    if ( next ) {
        nw->next->prev = nw;
    }

    nw->nextPoly = nextPoly;
    nw->neighbor = neighbor;
    nw->status = st;
    nw->entry = entry;
    nw->visited = visited;
    nw->alpha = alpha;
    return nw;
}

void
Graph :: insert(node *ins, node *first, node *last)
{
    node *aux = first;
    while ( aux != last && aux->alpha <= ins->alpha ) {
        aux = aux->next;
    }

    ins->next = aux;
    ins->prev = aux->prev;
    ins->prev->next = ins;
    ins->next->prev = ins;
}

void
Graph :: remove(node *rem)
{
    if ( rem ) {
        rem->neighbor->entry = 0;
        rem->prev->next = rem->next;
        rem->next->prev = rem->prev;
        delete rem;
        rem = NULL;
    }
}

bool
Graph :: intersectionExist(node *p1, node *p2, node *q1, node *q2)
{
    node *aux = q1->next;
    while ( aux != q2 ) {
        if ( aux->status == NS_Intersection ) {
            if ( ( aux->neighbor->prev == p1 ) && ( aux->neighbor->next == p2 ) ) {
                return true;
            }
        }

        aux = aux->next;
    }

    return false;
}


bool
Graph :: testCollapsedEdge(node *p1, node *p2)
{
    if ( dist(p1->x, p1->y, p2->x, p2->y) <= GT_EPS ) {
        return true;
    } else {
        return false;
    }
}

void
Graph :: removeIntersectionIfExist(node *p1, node *p2, node *q1, node *q2)
{
    node *anext;
    node *aux = q1->next;
    while ( aux != q2 ) {
        anext = aux->next;
        if ( aux->status == NS_Intersection ) {
            if ( ( prev_node(aux->neighbor->prev) == p1 ) && ( next_node(aux->neighbor->next) == p2 ) ) {
                // remove aux and aux->neighbor
                node *neighbor = aux->neighbor;
                aux->prev->next = aux->next;
                aux->next->prev = aux->prev;
                delete aux;

                neighbor->prev->next = neighbor->next;
                neighbor->next->prev = neighbor->prev;
                delete neighbor;
            }
        }

        aux = anext;
    }
}

int
Graph :: vertex2IntersectionVertex(node *v, node *l1, node *l2)
{
    // transforms given vertex as IntersectionVertex on given line
    // checks for coorect topology, ie check whether previos intersections of lines from IV to l1,l2 are present
    // returns 1 if newly formed
    // returns 0 if already exist
    node *nv = next_node(v->next);
    node *pv = prev_node(v->prev);

    if ( v->status == NS_IntersectionVertex ) {
        if ( v->neighbor->status == NS_IntersectionVertex ) {
            if ( ( v->neighbor == l1 ) || ( v->neighbor == l2 ) ) {
                return 0;
            } else {
                THROW_GT_EXCEPTIONM("Graph::vertex2IntersectionVertex: topology error, neighbor is unrelated intersectionVertex");
            }
        } else if ( ( prev_node(v->neighbor->prev) == l1 ) && ( next_node(v->neighbor->next) == l2 ) ) {
            return 0;
        } else if ( ( prev_node(v->neighbor->prev) == l2 ) && ( next_node(v->neighbor->next) == l1 ) ) {
            return 0;
        } else {
            // intersection vertex associated to neiborhing edge -> move it to common vertex
            // merge vertex together, remove intersection now associated to IntersectionVertex
            if ( ( l1 == prev_node(v->neighbor->prev) ) || ( l1 == next_node(v->neighbor->next) ) ) {
                merge2vertex(v, l1);
            } else if ( ( l2 == prev_node(v->neighbor->prev) ) || ( l2 == next_node(v->neighbor->next) ) ) {
                merge2vertex(v, l2);
            } else {
                THROW_GT_EXCEPTIONM("Graph::vertex2IntersectionVertex: topology error");
            }

            return 0;
        }
    } else {
        removeIntersectionIfExist(pv, v, l1, l2);
        removeIntersectionIfExist(v, nv, l1, l2);


        v->status = NS_IntersectionVertex;
        return 1;
    }
}


void
Graph :: testNewIntersectionVertexEdgeCollapse(node *v, node *l1, node *l2)
{
    double v_par, vp_par, vn_par;

    if ( l1->status == NS_IntersectionVertex && l1->neighbor->status == NS_Intersection ) {
        if ( next_node(l1->neighbor->next) == v ) {
            l1->neighbor->entry = -1; // to v
            l1->entry = -1; // from l1
        } else if ( prev_node(l1->neighbor->prev) == v ) {
            l1->entry = -1; // from l1
            v->entry = -1; // from v
        }
    } else if ( l2->status == NS_IntersectionVertex && l2->neighbor->status == NS_Intersection ) {
        if ( next_node(l2->neighbor->next) == v ) {
            l2->neighbor->entry = -1; // to v
            v->neighbor->entry = -1; // to l2
        } else if ( prev_node(l2->neighbor->prev) == v ) {
            v->entry = -1; // from v
            v->neighbor->entry = -1; // to l2
        }
    } else if ( belongs(prev_node(v->prev), l1, l2) ) {
        node *vp = prev_node(v->prev);
        vp->entry = -1;

        if ( vp->neighbor == l1 ) {
            vp_par = 0.;
        } else if ( vp->neighbor == l2 ) {
            vp_par = 1.;
        } else {
            vp_par = vp->neighbor->alpha;
        }

        if ( v->neighbor == l1 ) {
            v_par = 0.;
        } else if ( v->neighbor == l2 ) {
            v_par = 1.;
        } else {
            v_par = v->neighbor->alpha;
        }

        if ( v_par < vp_par ) {
            v->neighbor->entry = -1;
        } else {
            vp->neighbor->entry = -1;
        }
    } else if ( belongs(next_node(v->next), l1, l2) ) {
        node *vn = next_node(v->next);
        v->entry = -1;

        if ( vn->neighbor == l1 ) {
            vn_par = 0.;
        } else if ( vn->neighbor == l2 ) {
            vn_par = 1.;
        } else {
            vn_par = vn->neighbor->alpha;
        }

        if ( v->neighbor == l1 ) {
            v_par = 0.;
        } else if ( v->neighbor == l2 ) {
            v_par = 1.;
        } else {
            v_par = v->neighbor->alpha;
        }

        if ( v_par < vn_par ) {
            v->neighbor->entry = -1;
        } else {
            vn->neighbor->entry = -1;
        }
    }
}



void
Graph :: merge2vertex(node *v1, node *v2)
{
    bool boundary = false;
    if ( v1->status == NS_IntersectionVertex && v2->status == NS_IntersectionVertex ) {
        if ( v1->neighbor->status == NS_IntersectionVertex && v2->neighbor->status == NS_IntersectionVertex &&
             v1->neighbor == v2 && v2->neighbor == v1 ) {
            return;
        } else if ( v1->neighbor->status == NS_Intersection && v2->neighbor->status == NS_Intersection ) {
            remove(v1->neighbor);
            remove(v2->neighbor);
            v1->neighbor = v2;
            v2->neighbor = v1;
            v1->status = v2->status = NS_IntersectionVertex;
        } else {
            THROW_GT_EXCEPTIONM("Graph::merge2vertex: consistency error (one of merged vertices (or both) linked to another vertex already");
        }
    } else if ( v1->status == NS_IntersectionVertex ) {
        if ( v1->neighbor->status == NS_Intersection ) {
            boundary = ( v1->neighbor->entry == -1 );
            remove(v1->neighbor);
        }

        v1->neighbor = v2;
        v2->neighbor = v1;
        v2->status = NS_IntersectionVertex;
        if ( boundary ) {
            v2->entry = -1;
        }
    } else if ( v2->status == NS_IntersectionVertex ) {
        if ( v2->neighbor->status == NS_Intersection ) {
            boundary = ( v2->neighbor->entry == -1 );
            remove(v2->neighbor);
        }

        v2->neighbor = v1;
        v1->neighbor = v2;
        v1->status = NS_IntersectionVertex;
        if ( boundary ) {
            v1->entry = -1;
        }
    } else {
        v1->neighbor = v2;
        v2->neighbor = v1;
        v1->status = v2->status = NS_IntersectionVertex;
    }
}

int
Graph :: testPoint(node *s, double x, double y) const
{
    bool oddNODES = false;
    double x1, x2, y1, y2;

    if ( s ) {
        node *auxs = s;
        do {
            x1 = auxs->x;
            y1 = auxs->y;
            x2 = auxs->next->x;
            y2 = auxs->next->y;

            if ( ( ( y1 < y ) && ( y2 >= y ) ) ||
                ( ( y2 < y ) && ( y1 >= y ) ) ) {
                if ( x1 + ( y - y1 ) / ( y2 - y1 ) * ( x2 - x1 ) < x ) {
                    oddNODES = !oddNODES;
                }
            }

            auxs = auxs->next;
        } while ( auxs != s );
    }

    return oddNODES;
}


void
Graph :: printYourself()
{
    node *auxs, *auxc;
    /* print graph */
    printf("Graph 1/2:\n");
    auxs = s;
    do {
        if ( auxs->status == NS_Vertex ) {
            printf("Vertex [xy %e %e %e] [e=%d]\n", auxs->x, auxs->y, 0.0, auxs->entry);
        } else if ( auxs->status == NS_IntersectionVertex ) {
            printf("IVertex [xy %e %e %e] [e=%d]\n", auxs->x, auxs->y, 0.0, auxs->entry);
        } else {
            printf("Intrsc [xy %e %e %e] [e=%d]\n", auxs->x, auxs->y, 0.0, auxs->entry);
        }

        auxs = auxs->next;
    } while ( auxs != s );

    /* print graph */
    printf("Graph 2/2:\n");
    auxc = c;
    do {
        if ( auxc->status == NS_Vertex ) {
            printf("Vertex [xy %e %e %e] [e=%d]\n", auxc->x, auxc->y, 0.0, auxc->entry);
        } else if ( auxc->status == NS_IntersectionVertex ) {
            printf("IVertex [xy %e %e %e] [e=%d]\n", auxc->x, auxc->y, 0.0, auxc->entry);
        } else {
            printf("Intrsc [xy %e %e %e] [e=%d]\n", auxc->x, auxc->y, 0.0, auxc->entry);
        }

        auxc = auxc->next;
    } while ( auxc != c );

    /* end of debug printing */
}


void
GT_Exception :: print()
{
    fprintf(stderr, "\nGT_Exception thrown in %s:%d\n", file, line);
    if ( msg ) {
        fprintf(stderr, "msg: %s\n", msg);
    }
}
} // end namespace oofem
