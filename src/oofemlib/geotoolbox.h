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

#ifndef geotoolbox_h
#define geotoolbox_h

#ifndef __MAKEDEPEND
 #include <list>
 #include <cstdlib>
#endif
#ifdef __OOFEG
 #include "oofeggraphiccontext.h"
#endif
#include "flotarry.h"


namespace oofem {
#define GT_EPS 1.e-12
// zero for parallel lines test
#define GT_TEPS 1.e-16


/**
 * Class representing vertex.
 * The attributes of a vertex are its coordinates.
 * This class is a part of geometry toolbox.
 */
class Vertex
{
public:
    FloatArray coords;
    Vertex(double x = 0, double y = 0) : coords(2) {
        coords(0) = x;
        coords(1) = y;
    }
    Vertex(double c [ 2 ]) : coords(2) {
        coords(0) = c [ 0 ];
        coords(1) = c [ 1 ];
    }
    Vertex(const Vertex &src) { coords = src.coords; }
    Vertex & operator=(const Vertex &src) {
        coords = src.coords;
        return * this;
    }
    void setCoords(double x, double y) {
        coords(0) = x;
        coords(1) = y;
    }
    const FloatArray *getCoords() const { return & coords; }
};


/**
 * Class representing 2D polygon.
 * Polygon is represented as a sequence of vertices. THe services for adding vertices, computing volume,
 * testing if point is inside polygon are provided. Also iterators over polygon vertices and edges are provided.
 * This class is a part of geometry toolbox.
 */
class Polygon
{
    std :: list< Vertex >vertices;
public:
    Polygon() { };
    void addVertex(Vertex v) { vertices.push_back(v); }
    double computeVolume() const;
    int testPoint(double x, double y) const;
    double pointDistance(double x, double y) const;
    void clear() { vertices.clear(); }
#ifdef __OOFEG
    GraphicObj *draw(oofegGraphicContext &, bool filled, int layer = OOFEG_DEBUG_LAYER);
#endif
public:
    class PolygonEdgeIterator
    {
        bool last;
        const Polygon *ptr;
        Vertex curr;
        std :: list< Vertex > :: const_iterator iter;
public:
        PolygonEdgeIterator(const Polygon *p) : iter() {
            iter = p->vertices.begin();
            if ( iter == p->vertices.end() ) { last = true; } else {
                curr = ( * iter );
                ptr = p;
                last = false;
            }
        }
        int giveNext(Vertex &p1, Vertex &p2) {
            Vertex next;
            if ( last ) { return 0; } else {
                //next=*(++iter);
                ++iter;
                if ( iter == ptr->vertices.end() ) {
                    next = * ( ptr->vertices.begin() );
                    last = true;
                } else {
                    next = * iter;
                }

                p1 = curr;
                p2 = next;
                curr = next;
                return 1;
            }
        }
    };

    class PolygonVertexIterator
    {
        const Polygon *ptr;
        std :: list< Vertex > :: const_iterator iter;
public:
        PolygonVertexIterator(const Polygon *p) : iter() {
            iter = p->vertices.begin();
            ptr = p;
        }
        void init(const Polygon *p) {
            iter = p->vertices.begin();
            ptr = p;
        }
        int giveNext(Vertex &p1) {
            if ( iter == ptr->vertices.end() ) { return 0; } else { p1 = ( * iter ); }

            ++iter;
            return 1;
        }
        int giveNext(const Vertex **p1) {
            if ( iter == ptr->vertices.end() ) { return 0; } else { * p1 = iter.operator->(); }

            ++iter;
            return 1;
        }
    };
};

/**
 * Class representing the special graph constructed from two polygons that is used to perform
 * boolean operation on polygons (polygon clipping in current implementation).
 * In short, the graph contains polygon vertices and intersections of polygon edges.
 * Then graph edges represent the edges between graph nodes (graph edge correspond to polygon edge if
 * this edge has no intersection with other edges of second polygon, or single polygon edge can be represented by two or more
 * graph edges if intersections are present). Each polygon vertex representing intersection contain link to its
 * twin representing the same intersection on second polygon. This allows to perform traversal walk over graph edges and vertices
 * that can represent various boolean operations.
 * This class is a part of geometry toolbox.
 */
class Graph
{
protected:
    enum nodeStatus { NS_Vertex, NS_IntersectionVertex, NS_Intersection };
    struct node
    {
        double x, y;
        struct node *next;
        struct node *prev;
        struct node *nextPoly;     /* pointer to the next polygon */
        struct node *neighbor;     /* the corresponding intersection point */
        //int intersect;            /* 1 if an intersection point, 0 otherwise */
        nodeStatus status;
        int entry;            /* 1 if entry point (edge starting at this node is inside), -1 boundary edge, 0 otherwise */
        int visited;          /* 1 if the node has been visited, 0 otherwise */
        double alpha;          /* intersection point placement */
    };

    node *s, *c;
public:
    Graph() { s = c = NULL; }
    ~Graph();

    void clip(Polygon &result, const Polygon &a, const Polygon &b);
protected:
    void insert(node *ins, node *first, node *last);     // insert node in graph
    /** Create new node struct */
    node *createNode(double x, double y, node *next, node *prev, node *nextPoly,
                     node *neighbor, nodeStatus st, int entry, int visited, double alpha);
    node *next_node(node *p);     // return next node in graph
    node *prev_node(node *p);     // return prev node in graph
    node *last_node(node *p);     // returns last node in graph
    node *first(node *p);
    void  remove(node *n);
    double dist(double x1, double y1, double x2, double y2);

    int testIfIntersect(node *p1, node *p2, node *q1, node *q2,
                        double *alpha_p, double *alpha_q, double *xint, double *yint);

    int testIfCoincident(node *p1, node *p2, node *q1, node *q2, double *alpha_1, double *alpha_2);
    int testPoint(node *poly, double x, double y) const;

    bool belongs(node *n, node *v1, node *v2);
    bool intersectionExist(node *p1, node *p2, node *q1, node *q2);
    void removeIntersectionIfExist(node *p1, node *p2, node *q1, node *q2);
    int  vertex2IntersectionVertex(node *v, node *l1, node *l2);
    void testNewIntersectionVertexEdgeCollapse(node *v, node *l1, node *l2);
    bool testCollapsedEdge(node *p1, node *p2);
    void merge2vertex(node *v1, node *v2);

    void printYourself();

    void clear();
};

class GT_Exception
{
    const char *msg, *file;
    int line;

public:

    GT_Exception(const char *file, int line) {
        this->file = file;
        this->line = line;
        this->msg = NULL;
    }
    GT_Exception(const char *msg, const char *file, int line) {
        this->file = file;
        this->line = line;
        this->msg = msg;
    }
    ~GT_Exception() { }

    void print();
};

#define THROW_GT_EXCEPTION() throw GT_Exception(__FILE__, __LINE__);
#define THROW_GT_EXCEPTIONM(m) throw GT_Exception(m, __FILE__, __LINE__);
} // end namespace oofem
#endif // geotoolbox_h
