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

#include "trianglemesherinterface.h"
#include "flotarry.h"
#include "intarray.h"
#include "mathfem.h"

#include <set>
#include <cstdio>

#ifdef __TRIANGLE_MODULE
extern "C" {
#include <triangle.h>

// Function prototypes from triangle.c
void parsecommandline(int argc, char **argv, struct behavior *b);
void transfernodes(struct mesh *m, struct behavior *b, REAL *pointlist, REAL *pointattriblist, int *pointmarkerlist, int numberofpoints, int numberofpointattribs);
void formskeleton(struct mesh *m, struct behavior *b, int *segmentlist, int *segmentmarkerlist, int numberofsegments);
void highorder(struct mesh *m, struct behavior *b);
void numbernodes(struct mesh *m, struct behavior *b);
void writeelements(struct mesh *m, struct behavior *b, int **trianglelist, REAL **triangleattriblist);
void writepoly(struct mesh *m, struct behavior *b, int **segmentlist, int **segmentmarkerlist);
void writeedges(struct mesh *m, struct behavior *b, int **edgelist, int **edgemarkerlist);
void statistics(struct mesh *m, struct behavior *b);
void writenodes(struct mesh *m, struct behavior *b, REAL **pointlist, REAL **pointattriblist, int **pointmarkerlist);
void checkmesh(struct mesh *m, struct behavior *b);
void checkdelaunay(struct mesh *m, struct behavior *b);
void writeneighbors(struct mesh *m, struct behavior *b, int **neighborlist);
void writevoronoi(struct mesh *m, struct behavior *b, REAL **vpointlist, REAL **vpointattriblist, int **vpointmarkerlist, int **vedgelist, int **vedgemarkerlist, REAL **vnormlist);
int reconstruct(struct mesh *m, struct behavior *b, int *trianglelist, REAL *triangleattriblist, REAL *trianglearealist, int elements, int corners, int attribs, int *segmentlist,int *segmentmarkerlist, int numberofsegments);
subseg *subsegtraverse(struct mesh *m);
void regionplague(struct mesh *m, struct behavior *b, REAL attribute, REAL area);
void plague(struct mesh *m, struct behavior *b);
void poolinit(struct memorypool *pool, int bytecount, int itemcount, int firstitemcount, int alignment);
void pooldeinit(struct memorypool *pool);

// Macros from triangle.c
#define setelemattribute(otri, attnum, value)                                 \
  ((REAL *) (otri).tri)[m->elemattribindex + (attnum)] = value
#define elemattribute(otri, attnum)                                           \
  ((REAL *) (otri).tri)[m->elemattribindex + (attnum)]

#define stpivot(osub, otri)                                                   \
  ptr = (triangle) (osub).ss[6 + (osub).ssorient];                            \
  decode(ptr, otri)
#define ssymself(osub)                                                        \
  (osub).ssorient = 1 - (osub).ssorient
#define mark(osub)  (* (int *) ((osub).ss + 8))
#define decode(ptr, otri)                                                     \
  (otri).orient = (int) ((long)ptr & (long)3l);         \
  (otri).tri = (triangle *) ((long)(ptr) ^ (long)(otri).orient)
#define infect(otri)                                                          \
  (otri).tri[6] = (triangle) ((long)(otri).tri[6] | (long)2l)
#define org(otri, vertexptr)                                                  \
  vertexptr = (vertex) (otri).tri[plus1mod3[(otri).orient] + 3]
#define dest(otri, vertexptr)                                                 \
  vertexptr = (vertex) (otri).tri[minus1mod3[(otri).orient] + 3]
#define apex(otri, vertexptr)                                                 \
  vertexptr = (vertex) (otri).tri[(otri).orient + 3]
#define sorg(osub, vertexptr)                                                 \
  vertexptr = (vertex) (osub).ss[2 + (osub).ssorient]
#define sdest(osub, vertexptr)                                                \
  vertexptr = (vertex) (osub).ss[3 - (osub).ssorient]


#define VIRUSPERBLOCK 1020   /* Number of virus triangles allocated at once. */

// This customized version carves holes defined by the regions.
// It'll remove any element which doesn't receive any region.
int custom_carveholes(struct mesh *m, struct behavior *b, const int *outside, const int *inside)
{
    struct otri neighbortri;
    struct otri triangleloop;
    struct osub subsegloop;
    triangle **temptri;
    triangle ptr;                         /* Temporary variable used by sym(). */
    double area;
    double attribute;
    double triangle_attribute;
    int sane_mesh;

    poolinit(&m->viri, sizeof(triangle *), VIRUSPERBLOCK, VIRUSPERBLOCK, 0);

    /* Assigns every triangle a regional attribute of -1 */
    traversalinit(&m->triangles);
    triangleloop.orient = 0;
    triangleloop.tri = triangletraverse(m);
    while (triangleloop.tri != (triangle *) NULL) {
      setelemattribute(triangleloop, m->eextras, -1.0);
      triangleloop.tri = triangletraverse(m);
    }

    sane_mesh = 1;

    /* Loop over all segments */
    traversalinit(&m->subsegs);
    subsegloop.ss = subsegtraverse(m);
    subsegloop.ssorient = 0;
    while (subsegloop.ss != (subseg *) NULL) {
        if (subsegloop.ss != m->dummysub) {
            /* First neighbor. */
            ssymself(subsegloop);
            stpivot(subsegloop, neighbortri);
            if ( neighbortri.tri != m->dummytri ) {
                area = 0.0; /// @todo Possible set this as the minimum as well.
                attribute = outside[mark(subsegloop)-1];
                triangle_attribute = elemattribute(neighbortri, m->eextras);

                /* The region no number yet. */
                if (triangle_attribute < 0 && attribute >= 0) {
                    infect(neighbortri);
                    temptri = (triangle **) poolalloc(&m->viri);
                    *temptri = neighbortri.tri;
                    /* Apply one region's attribute and/or area constraint. */
                    regionplague(m, b, attribute, area);
                    /* The virus pool should be empty now. */
                } else if (attribute >= 0 && triangle_attribute >= 0 && attribute != triangle_attribute) {
                    /* Check for problems. */
                    vertex v1, v2, v3;
                    org(neighbortri, v1);
                    dest(neighbortri, v2);
                    apex(neighbortri, v3);
                    fprintf(stdout, "Error: inconsistent region information (new %d, old %d from %d) (outside) at (%e, %e)\n",
                            (int)attribute, (int)triangle_attribute,(int)mark(subsegloop), (v1[0]+v2[0]+v3[0])/3, (v1[1]+v2[1]+v3[1])/3);
                    sane_mesh = 0;
                }
            }

            /* Second neighbor (same procedure). */
            ssymself(subsegloop);
            stpivot(subsegloop, neighbortri);
            if ( neighbortri.tri != m->dummytri) {
                area = 0.0;
                attribute = inside[mark(subsegloop)-1];
                triangle_attribute = elemattribute(neighbortri, m->eextras);

                if (triangle_attribute < 0 && attribute >= 0) {
                    infect(neighbortri);
                    temptri = (triangle **) poolalloc(&m->viri);
                    *temptri = neighbortri.tri;
                    regionplague(m, b, attribute, area);
                } else if (attribute >= 0 && triangle_attribute >= 0 && attribute != triangle_attribute) {
                    vertex v1, v2, v3;
                    org(neighbortri, v1);
                    dest(neighbortri, v2);
                    apex(neighbortri, v3);
                    fprintf(stdout, "Error: inconsistent region information (new %d, old %d from %d) (inside) at (%e, %e)\n",
                            (int)attribute, (int)triangle_attribute,(int)mark(subsegloop), (v1[0]+v2[0]+v3[0])/3, (v1[1]+v2[1]+v3[1])/3 );
                    sane_mesh = 0;
                }
            }

            subsegloop.ss = subsegtraverse(m);
        }
    }

    /* Remove all triangles with marker 0.0 */
    traversalinit(&m->triangles);
    triangleloop.tri = triangletraverse(m);
    int triangle_number=0;
    while (triangleloop.tri != (triangle *) NULL) {
        if ( triangleloop.tri != m->dummytri ) {
            triangle_number++;
            attribute = elemattribute(triangleloop, m->eextras);
            if (attribute == -1.0) {
                fprintf(stderr, "Broken mesh at triangle %d\n",triangle_number);
                sane_mesh = 0;
            } else if (attribute == 0.0) {
                infect(triangleloop);
                temptri = (triangle **) poolalloc(&m->viri);
                *temptri = triangleloop.tri;
            }
        }
        triangleloop.tri = triangletraverse(m);
    }
    /* Remove the marked elements */
    plague(m, b);

    if (b->regionattrib && !b->refine) {
        /* Note the fact that each triangle has an additional attribute. */
        m->eextras++;
    }

    /* Free up memory. */
    pooldeinit(&m->viri);

    return sane_mesh;
}

// Customized triangulation call.
void custom_triangulate(char *triswitches, struct triangulateio *in, struct triangulateio *out, struct triangulateio *vorout, const int *outside, const int *inside)
{
    struct mesh m;
    struct behavior b;
    REAL *holearray;
    REAL *regionarray;
    triangleinit(&m);
    parsecommandline(1, &triswitches, &b);
    //b.verbose=2;
    m.steinerleft = b.steiner;
    transfernodes(&m, &b, in->pointlist, in->pointattributelist,
            in->pointmarkerlist, in->numberofpoints,
            in->numberofpointattributes);
    if (b.refine) {
        m.hullsize = reconstruct(&m, &b, in->trianglelist,
                in->triangleattributelist, in->trianglearealist,
                in->numberoftriangles, in->numberofcorners,
                in->numberoftriangleattributes,
                in->segmentlist, in->segmentmarkerlist,
                in->numberofsegments);
    } else {
        m.hullsize = delaunay(&m, &b);
    }
    m.infvertex1 = (vertex) NULL;
    m.infvertex2 = (vertex) NULL;
    m.infvertex3 = (vertex) NULL;
    if (b.usesegments) {
        m.checksegments = 1;
        if (!b.refine) {
            formskeleton(&m, &b, in->segmentlist,
                    in->segmentmarkerlist, in->numberofsegments);
        }
    }
#if 0
    struct osub subsegloop;
    traversalinit(&m.subsegs);
    subsegloop.ss = subsegtraverse(&m);
    subsegloop.ssorient = 0;
    while (subsegloop.ss != (subseg *) NULL) {
        if (subsegloop.ss != m.dummysub) {
            REAL *p1, *p2;
            sorg(subsegloop,p1);
            sdest(subsegloop,p2);
            printf("  Connected (%f,%f) to (%f,%f)\n",p1[0],p1[1],p2[0],p2[1]);
            subsegloop.ss = subsegtraverse(&m);
        }
    }
#endif
    if (b.poly && (m.triangles.items > 0)) {
        holearray = in->holelist;
        m.holes = in->numberofholes;
        regionarray = in->regionlist;
        m.regions = in->numberofregions;
        if (!b.refine) {
            /* Only increase quality if the regions are properly defined. */
            int sane = custom_carveholes(&m, &b, outside, inside);
            b.quality *= sane;
            if (sane == 0) { printf("Probably bad PSLG\n"); exit(-1); }
            
        }
    } else {
        m.holes = 0;
        m.regions = 0;
    }
    if (b.quality && (m.triangles.items > 0)) {
        enforcequality(&m, &b);
    }
    m.edges = (3l * m.triangles.items + m.hullsize) / 2l;
    if (b.order > 1) {
        highorder(&m, &b);
    }
    if (!b.quiet) {
        printf("\n");
    }
    if (b.jettison) {
        out->numberofpoints = m.vertices.items - m.undeads;
    } else {
        out->numberofpoints = m.vertices.items;
    }
    out->numberofpointattributes = m.nextras;
    out->numberoftriangles = m.triangles.items;
    out->numberofcorners = (b.order + 1) * (b.order + 2) / 2;
    out->numberoftriangleattributes = m.eextras;
    out->numberofedges = m.edges;
    if (b.usesegments) {
        out->numberofsegments = m.subsegs.items;
    } else {
        out->numberofsegments = m.hullsize;
    }
    if (vorout != (struct triangulateio *) NULL) {
        vorout->numberofpoints = m.triangles.items;
        vorout->numberofpointattributes = m.nextras;
        vorout->numberofedges = m.edges;
    }
    if (b.nonodewritten || (b.noiterationnum && m.readnodefile)) {
        if (!b.quiet) {
            printf("NOT writing vertices.\n");
        }
        numbernodes(&m, &b);
    } else {
        writenodes(&m, &b, &out->pointlist, &out->pointattributelist,
                &out->pointmarkerlist);
    }

    // Simp. always write the triangles.
    writeelements(&m, &b, &out->trianglelist, &out->triangleattributelist);

    if (b.poly || b.convex) {
        writepoly(&m, &b, &out->segmentlist, &out->segmentmarkerlist);
        out->numberofholes = m.holes;
        out->numberofregions = m.regions;
        if (b.poly) {
            out->holelist = in->holelist;
            out->regionlist = in->regionlist;
        } else {
            out->holelist = (REAL *) NULL;
            out->regionlist = (REAL *) NULL;
        }
    }
    if (b.edgesout) {
        writeedges(&m, &b, &out->edgelist, &out->edgemarkerlist);
    }
    // Simp. no voronoi
    if (b.neighbors) {
        writeneighbors(&m, &b, &out->neighborlist);
    }
    // Simp. No statistics.
    if (b.docheck) {
        checkmesh(&m, &b);
        checkdelaunay(&m, &b);
    }
    triangledeinit(&m, &b);
}

}

// Convenience function.
void clearTriangulateIO(struct triangulateio &t)
{
    t.pointlist = NULL;
    t.pointattributelist = NULL;
    t.pointmarkerlist = NULL;
    t.numberofpoints = 0;
    t.numberofpointattributes = 0;

    t.trianglelist = NULL;
    t.triangleattributelist = NULL;
    t.trianglearealist = NULL;
    t.numberoftriangles = 0;
    t.numberofcorners = 0;
    t.numberoftriangleattributes = 0;

    t.segmentlist = NULL;
    t.segmentmarkerlist = NULL;
    t.numberofsegments = 0;

    t.holelist = NULL;
    t.numberofholes = 0;

    t.regionlist = NULL;
    t.numberofregions = 0;
}
#endif

namespace oofem {

bool TriangleMesherInterface :: meshPSLG(const Triangle_PSLG &pslg,
        const IntArray &outside, const IntArray &inside,
        AList<FloatArray> &nodes, AList<IntArray> &n_markers,
        AList<IntArray> &triangles, IntArray &t_markers,
        AList<IntArray> &segments, IntArray &s_markers) const
{
#ifdef __TRIANGLE_MODULE
    // 1. Fill the struct for triangle;
    struct triangulateio mesh;
    clearTriangulateIO(mesh);

    // 1.a. Copy over the node data.
    mesh.numberofpoints = pslg.nx.giveSize();
    mesh.pointlist = new REAL[mesh.numberofpoints*2];
    //mesh.pointmarkerlist = new REAL[mesh.numberofpoints];
    for (int i = 0; i < mesh.numberofpoints; ++i) {
        mesh.pointlist[i*2] = pslg.nx(i);
        mesh.pointlist[i*2+1] = pslg.ny(i);
        //mesh.pointmarkerlist[i] = pslg.n_marker(i);
    }

    // 1.b. Copy over the segment data
    printf("Copying segment data\n");
    mesh.numberofsegments = pslg.segment_a.giveSize();
    mesh.segmentlist = new int[mesh.numberofsegments*2];
    for (int i = 0; i < mesh.numberofsegments; ++i) {
        mesh.segmentlist[i*2] = pslg.segment_a(i);
        mesh.segmentlist[i*2+1] = pslg.segment_b(i);
    }

    if (pslg.segment_marker.giveSize() > 0) {
        mesh.segmentmarkerlist = new int[mesh.numberofsegments];
        for (int i = 0; i < mesh.numberofsegments; ++i) {
            mesh.segmentmarkerlist[i] = pslg.segment_marker(i);
        }
    }

    // 2. Triangulate
    char options[100];
    // Note: Not sure if -A is necessary when using the library interface.
    sprintf(options, "-p -q %f -a%f %s -A",this->minAngle, this->maxArea, this->quadratic ? "-o2":"");
    struct triangulateio output;
    clearTriangulateIO(output);
    custom_triangulate (options, &mesh, &output, NULL, outside.givePointer(), inside.givePointer());

    // 3. Copy back
    nodes.growTo(output.numberofpoints);
    //n_markers.resize(output.numberofpoints);
    for (int i = 0; i < output.numberofpoints; ++i) {
        FloatArray *node = new FloatArray(2);
        node->at(1) = output.pointlist[i*2];
        node->at(2) = output.pointlist[i*2+1];
        nodes.put(i+1, node);
        //n_markers(i) = output.pointmarkerlist[i]; // Not enough.
    }

    triangles.growTo(output.numberoftriangles);
    t_markers.resize(output.numberoftriangles);
    for (int i = 0; i < output.numberoftriangles; ++i) {
        IntArray *triangle = new IntArray(output.numberofcorners);
        for (int j = 0; j < 3; j++) { // First three
            triangle->at(j+1) = output.trianglelist[i*output.numberofcorners + j];
        }
        // Rearrange the strange ordering of the edge nodes.
        if (output.numberofcorners == 6) {
            triangle->at(4) = output.trianglelist[i*output.numberofcorners + 5];
            triangle->at(5) = output.trianglelist[i*output.numberofcorners + 3];
            triangle->at(6) = output.trianglelist[i*output.numberofcorners + 4];
        }
        triangles.put(i+1, triangle);
        t_markers.at(i+1) = round(output.triangleattributelist[i]);
    }

    // A somewhat annoying missing feature of triangle, it won't make the segments quadratic.
    std::set<int> *node_triangle;
    if (this->quadratic) {
        node_triangle = new std::set<int>[output.numberofpoints];
        for (int i = 1; i <= triangles.giveSize(); ++i) {
            IntArray *triangle = triangles.at(i);
            for (int j = 1; j <= 3; ++j) {
                node_triangle[triangle->at(j)-1].insert(i);
            }
        }
    }

    segments.growTo(output.numberofsegments);
    s_markers.resize(output.numberofsegments);
    for (int i = 0; i < output.numberofsegments; ++i) {
        IntArray *segment = this->quadratic ? new IntArray(3) : new IntArray(2);
        segment->at(1) = output.segmentlist[i*2 + 0];
        segment->at(2) = output.segmentlist[i*2 + 1];
        //segment->at(3) = output.segmentlist[i*3 + 2]; // Quadratic meshes only, not for segments.
        if (this->quadratic) {
            int a, b, c;
            std::set<int> tris = node_triangle[segment->at(1)-1];
            // Now look up any triangle with the other point included.
            for (std::set<int>::iterator it = tris.begin(); it != tris.end(); ++it) {
                IntArray *triangle = triangles.at(*it);
                if ( (b = triangle->findFirstIndexOf(segment->at(2))) > 0) {
                    a = triangle->findFirstIndexOf(segment->at(1));
                    if (a+b == 3) {
                        c = 4;
                    } else if (a+b == 5) {
                        c = 5;
                    } else {
                        c = 6;
                    }
                    segment->at(3) = triangle->at(c);
                    break;
                }
            }
        }
        segments.put(i+1, segment);
        s_markers.at(i+1) = output.segmentmarkerlist[i];
    }
    if (this->quadratic) {
        delete[] node_triangle;
    }

    this->fixNodeMarkers(nodes, n_markers, triangles, t_markers, segments, s_markers);

    // Deleting old memory
    delete[] mesh.pointlist;
    //delete[] mesh.pointattributelist;
    //delete[] mesh.pointmarkerlist;
    //delete[] mesh.trianglelist;
    //delete[] mesh.triangleattributelist;
    //delete[] mesh.trianglearealist;
    delete[] mesh.segmentlist;
    delete[] mesh.segmentmarkerlist;
    //if (mesh.holelist) { delete[] mesh.holelist; }
    //if (mesh.regionlist) { delete[] mesh.regionlist; }

    // Holes and regions are referenced in both.
    free(output.pointlist);
    //free(output.pointattributelist);
    free(output.pointmarkerlist);
    free(output.trianglelist);
    free(output.triangleattributelist);
    //free(output.trianglearealist);
    free(output.segmentlist);
    free(output.segmentmarkerlist);
    //free(output.holelist);
    //free(output.regionlist);
#else
    OOFEM_ERROR("TriangleMesherInterface :: meshPSLG - OOFEM is not compiled with support for triangle.");
#endif

    return true;
}

void TriangleMesherInterface :: fixNodeMarkers(const AList<FloatArray> &nodes, AList<IntArray> &n_markers,
        const AList<IntArray> &triangles, const IntArray &t_markers,
        const AList<IntArray> &segments, const IntArray &s_markers)
{
    n_markers.growTo(nodes.giveSize());
    for (int i = 1; i <= n_markers.giveSize(); ++i) {
        n_markers.put(i, new IntArray());
        //n_markers.at(i)->preallocate(nregions);
        ///@todo Maybe just a suboptimization.
    }

    for (int i = 1; i <= segments.giveSize(); ++i) {
        IntArray *sn = segments.at(i);
        for (int j = 1; j <= sn->giveSize(); ++j) {
            n_markers.at(sn->at(j))->insertSortedOnce(s_markers.at(i));
        }
    }
    for (int i = 1; i <= triangles.giveSize(); ++i) {
        IntArray *tn = triangles.at(i);
        for (int j = 1; j <= tn->giveSize(); ++j) {
            n_markers.at(tn->at(j))->insertSortedOnce(t_markers.at(i));
        }
    }
}

void TriangleMesherInterface :: simplifyPSLG(Triangle_PSLG &coarse, const Triangle_PSLG &pslg, double limit, double minlen)
{
    int segments = pslg.segment_a.giveSize();
    int nodes = pslg.nx.giveSize();

    // Calculate the inverted connection node->element
    std::set<int> *connectivity = new std::set<int>[nodes];
    for (int i = 0; i < segments; i++) {
        connectivity[pslg.segment_a(i)-1].insert(i+1);
        connectivity[pslg.segment_b(i)-1].insert(i+1);
    }

    // Some conservative error measure
    IntArray nodeRemoval(nodes);
    IntArray edgeRemoval(segments);
    edgeRemoval.zero();
    IntArray seg_a = pslg.segment_a;
    IntArray seg_b = pslg.segment_b;

    nodeRemoval.zero();
    FloatArray error_x(nodes), error_y(nodes), ab(2), ac(2), bc(2), err_vec(2);
    error_x.zero();
    error_y.zero();

    for (int j = 0; j < 2; ++j) {
        double allowed = j*limit/2;
        for (int i = 1; i <= nodes; i++) {
            std::set<int> &elems = connectivity[i-1];
            if (elems.size() == 2) {
                int e0 = *elems.begin();
                int e1 = *(++elems.begin());
                if (pslg.segment_marker.at(e0) == pslg.segment_marker.at(e1)) {
                    double abac, acac;
                    int n0, n1, n2;

                    n0 = seg_a.at(e0) == i ? seg_b.at(e0) : seg_a.at(e0);
                    n1 = i;
                    n2 = seg_a.at(e1) == i ? seg_b.at(e1) : seg_a.at(e1);

                    ab(0) = pslg.nx.at(n1) - pslg.nx.at(n0);
                    ab(1) = pslg.ny.at(n1) - pslg.ny.at(n0);
                    ac(0) = pslg.nx.at(n2) - pslg.nx.at(n0);
                    ac(1) = pslg.ny.at(n2) - pslg.ny.at(n0);

                    abac = ab.dotProduct(ac);
                    acac = ac.computeSquaredNorm();

                    // Find the error (how far the point would be moved to obtain the following line without it).
                    if (abac <= 0 || acac == 0) { // then -ab
                        err_vec(0) = -ab(0);
                        err_vec(1) = -ab(1);
                    } else if (abac >= acac) { // then bc
                        err_vec(0) = pslg.nx.at(n2) - pslg.nx.at(n1);
                        err_vec(1) = pslg.ny.at(n2) - pslg.ny.at(n1);
                    } else {
                        err_vec = ac;
                        err_vec.times(abac/acac);
                        err_vec.subtract(ab);
                    }

                    double ev_norm = err_vec.computeNorm();
                    // Max of new or old error;
                    double real_error = 0.0;
                    if (ev_norm == 0) {
                        real_error = 0.0;
                    } else {
                        error_x.at(n1) += err_vec(0);
                        error_y.at(n1) += err_vec(1);
                        real_error = sqrt(error_x.at(n1)*error_x.at(n1) + error_y.at(n1)*error_y.at(n1));
                    }

                    if (real_error <= allowed) {
                        // Mark node for removal, remove second edge, more first edge connection to next node;
                        nodeRemoval.at(i) = true;
                        edgeRemoval.at(e1) = true;
                        connectivity[n2-1].erase(e1);
                        connectivity[n2-1].insert(e0);
                        if (seg_a.at(e0) == n1) { // Doing the bothersome way to preserve direction of segments.
                            seg_a.at(e0) = n2;
                        } else {
                            seg_b.at(e0) = n2;
                        }
                        elems.clear();
                        // Accumulate the error vector
                        error_x.at(n0) += error_x.at(n1);
                        error_y.at(n0) += error_y.at(n1);
                        error_x.at(n2) += error_x.at(n1);
                        error_y.at(n2) += error_y.at(n1);
                    }
                }
            }
        }
    }

    // Deleting the elements which are too short.
    bool edgeRemoved = true;
    while (edgeRemoved) {
        edgeRemoved = false;
        for (int i = 1; i <= nodes; i++) {
            std::set<int> &elems = connectivity[i-1];
            if (elems.size() == 2) {
                int e0 = *elems.begin();
                int e1 = *(++elems.begin());
                if (pslg.segment_marker.at(e0) == pslg.segment_marker.at(e1)) {
                    int n0, n1, n2;

                    n0 = seg_a.at(e0) == i ? seg_b.at(e0) : seg_a.at(e0);
                    n1 = i;
                    n2 = seg_a.at(e1) == i ? seg_b.at(e1) : seg_a.at(e1);

                    ab(0) = pslg.nx.at(n1) - pslg.nx.at(n0);
                    ab(1) = pslg.ny.at(n1) - pslg.ny.at(n0);
                    bc(0) = pslg.nx.at(n2) - pslg.nx.at(n1);
                    bc(1) = pslg.ny.at(n2) - pslg.ny.at(n1);

                    if (ab.computeSquaredNorm() < minlen*minlen || bc.computeSquaredNorm() < minlen*minlen ) {
                        // Mark node for removal, remove second edge, more first edge connection to next node;
                        nodeRemoval.at(i) = true;
                        edgeRemoval.at(e1) = true;
                        connectivity[n2-1].erase(e1);
                        connectivity[n2-1].insert(e0);
                        if (seg_a.at(e0) == n1) { // Doing the bothersome way to preserve direction of segments.
                            seg_a.at(e0) = n2;
                        } else {
                            seg_b.at(e0) = n2;
                        }
                        elems.clear();
                        edgeRemoved = true;
                    }
                }
            }
        }
    }

    delete[] connectivity;

    // Cleanup
    int newNodes = 0;
    IntArray newNumber(nodes);
    newNumber.zero();
    coarse.nx.resize(nodes);
    coarse.ny.resize(nodes);
    for (int i = 1; i <= nodes; i++) {
        if (!nodeRemoval.at(i)) {
            newNodes++;
            coarse.nx.at(newNodes) = pslg.nx.at(i);
            coarse.ny.at(newNodes) = pslg.ny.at(i);
            newNumber.at(i) = newNodes;
        }
    }
    coarse.nx.resize(newNodes);
    coarse.ny.resize(newNodes);

    int newSegments = 0;
    coarse.segment_a.resize(segments);
    coarse.segment_b.resize(segments);
    coarse.segment_marker.resize(segments);

    for (int i = 1; i <= segments; i++) {
        if (!edgeRemoval.at(i)) {
            newSegments++;
            coarse.segment_a.at(newSegments) = newNumber.at(seg_a.at(i));
            coarse.segment_b.at(newSegments) = newNumber.at(seg_b.at(i));
            coarse.segment_marker.at(newSegments) = pslg.segment_marker.at(i);
        }
    }
    coarse.segment_a.resize(newSegments);
    coarse.segment_b.resize(newSegments);
    coarse.segment_marker.resize(newSegments);
}

};
