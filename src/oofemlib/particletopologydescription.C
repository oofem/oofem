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

#include "particletopologydescription.h"
#include "timestep.h"
#include "dofiditem.h"
#include "spatiallocalizer.h"
#include "mathfem.h"
#include "trianglemesherinterface.h"
#include "topologydescription.h"
#include "particlegrid.h"
#include "datareader.h"
#include "element.h"
#include "timer.h"
#include "classfactory.h"
#include "engngm.h"

#include "mixedgradientpressureneumann.h"
#include "node.h"
#include "activebc.h"
#include "activedof.h"
#include "masterdof.h"
#include "set.h"

#include <cstring>
#include <list>
#include <sstream>
#include <cstdio>

#ifdef __VTK_MODULE
 #include <vtkPoints.h>
 #include <vtkPointData.h>
 #include <vtkDoubleArray.h>
 #include <vtkVertex.h>
 #include <vtkCellArray.h>
 #include <vtkCellData.h>
 #include <vtkXMLPolyDataWriter.h>
 #include <vtkPolyData.h>
 #include <vtkSmartPointer.h>
#endif

namespace oofem
{
REGISTER_TopologyDescription(ParticleTopologyDescription);

ParticleTopologyDescription :: ParticleTopologyDescription(Domain *d) : TopologyDescription(d), resampled(false), maxdisp2(0)
{}

ParticleTopologyDescription :: ~ParticleTopologyDescription()
{
    if ( this->corners.size() ) {
        this->corners.clear();
    }
}

bool ParticleTopologyDescription :: instanciateYourself(DataReader *dr)
{
    // Required by IR_GIVE_FIELD macro
    IRResultType result;
    InputRecord *ir;
    std :: string name;
    int num;

    int res, nsd;
    double tubewidth;
    IntArray resolution;
    FloatArray bb0, bb1;
    // Input for geometry types
    FloatArray x0, x1, center;
    double radius, alpha0, alpha1;
    int id, nsegments;

    // Read topology description
    ir = dr->giveInputRecord(DataReader :: IR_gbpmRec, 1);
    IR_GIVE_FIELD(ir, nsd, _IFT_ParticleTopologyDescription_nsd);
    IR_GIVE_FIELD(ir, res, _IFT_ParticleTopologyDescription_baseResolution);
    IR_GIVE_FIELD(ir, bb0, _IFT_ParticleTopologyDescription_boundingBoxA);
    IR_GIVE_FIELD(ir, bb1, _IFT_ParticleTopologyDescription_boundingBoxB);
    IR_GIVE_FIELD(ir, nsegments, _IFT_ParticleTopologyDescription_numberOfSegments);
    tubewidth = 1.1;
    IR_GIVE_OPTIONAL_FIELD(ir, tubewidth, _IFT_ParticleTopologyDescription_tubeWidth);
    IR_GIVE_FIELD(ir, this->m, _IFT_ParticleTopologyDescription_neighbors);

    resolution.resize(nsd);
    resolution.zero();
    resolution.add(res);

    // Mappings
    IR_GIVE_OPTIONAL_FIELD(ir, this->regionInside, _IFT_ParticleTopologyDescription_regionInside);
    IR_GIVE_OPTIONAL_FIELD(ir, this->regionOutside, _IFT_ParticleTopologyDescription_regionOutside);
    IR_GIVE_OPTIONAL_FIELD(ir, this->regionElementType, _IFT_Meshing_elementType);
    IR_GIVE_OPTIONAL_FIELD(ir, this->regionSet, _IFT_Meshing_set);

    /// @todo Merging surfaces and such.
    int n = this->regionInside.giveSize();
    this->controlID.resize(n, n);
    this->mergeID.resize(n, n);

    this->tubeWidth = ( bb1(0) - bb0(0) ) / res * tubewidth;

    this->grid.reset( new ParticleGrid< ParticlePoint >(resolution, bb0, bb1) );

    for ( int i = 1; i <= nsegments; i++ ) {
        ir = dr->giveInputRecord(DataReader :: IR_geometryRec, i);
        IR_GIVE_RECORD_KEYWORD_FIELD(ir, name, num);
        IR_GIVE_FIELD(ir, id, _IFT_ParticleTopologyDescription_identification);
        if ( name.compare("line") == 0 ) {
            IR_GIVE_FIELD(ir, x0, _IFT_Line_start);
            IR_GIVE_FIELD(ir, x1, _IFT_Line_end);
            this->addLineSegment(id, x0, x1, *this->grid);
        } else if ( name.compare("circle") == 0 ) {
            IR_GIVE_FIELD(ir, center, _IFT_Circle_center);
            IR_GIVE_FIELD(ir, radius, _IFT_Circle_radius);
            alpha0 = -180;
            alpha1 = 180;
            IR_GIVE_OPTIONAL_FIELD(ir, alpha0, _IFT_Circle_start);
            IR_GIVE_OPTIONAL_FIELD(ir, alpha1, _IFT_Circle_end);
            alpha0 *= M_PI / 180;
            alpha1 *= M_PI / 180;
            this->addCircleSegment(id, center, radius, alpha0, alpha1, *this->grid);
        } else if ( name.compare("cornerpoint") == 0 ) {
            IR_GIVE_FIELD(ir, center, _IFT_Point_coords);
            this->addCorner(id, center, *this->grid);
        } else {
            OOFEM_ERROR( "Unknown geometry type '%s'", name.c_str() );
        }
        ir->finish();
    }

    this->checkOverlap();

#ifdef DEBUG
    //    this->writeVTKFile("initial_particle.vtp");
#endif
    /// @todo User input.
    this->writeVTK = false;
    return true;
}

void ParticleTopologyDescription :: addLineSegment(int id, const FloatArray &p0, const FloatArray &p1, ParticleGrid< ParticlePoint > &g) const
{
    ParticlePoint *point;
    double tubeWidth2 = this->tubeWidth * this->tubeWidth;
    double distance2, xi;
    double line_length;
    FloatArray v, new_normal(2), grid_p, new_foot, grid_p_v1;
    v.beDifferenceOf(p1, p0);
    line_length = v.computeNorm();
    v.times(1 / line_length);  // v is the directional vector along the line.
    new_normal(0) = -v(1);
    new_normal(1) = v(0);

    // Find the bounding box for the segment
    FloatArray x0, x1;
    x0.beMinOf(p0, p1);
    x0.add(-this->tubeWidth);
    x1.beMaxOf(p0, p1);
    x1.add(this->tubeWidth);

    for ( ParticleGrid< ParticlePoint > :: iterator it = this->grid->beginAt(x0, x1); !it.end(); ++it ) {
        it.getGridPoint(grid_p);
        // Find closest point on line, parameterized by xi
        grid_p_v1.beDifferenceOf(grid_p, p0);
        xi = grid_p_v1.dotProduct(v) / line_length;
        if ( xi < 0 || xi > 1 ) {
            continue;
        }

        // And it's coordinate
        new_foot = p0;
        new_foot.add(xi * line_length, v);

        distance2 = grid_p.distance_square(new_foot);

        if ( distance2 < tubeWidth2 ) {
            if ( ( point = it.getPoint() ) ) {
                if ( point->distance2 <= distance2 ) {
                    continue; // If distance already is smaller, then skip this.
                }
            }
            point = new ParticlePoint(new_foot, id, new_normal, distance2);
            it.setPoint(point);
        }
    }
}

void ParticleTopologyDescription :: addCircleSegment(int id, const FloatArray &c, double r, double v0, double v1,
                                                     ParticleGrid< ParticlePoint > &g) const
{
    ParticlePoint *point;
    FloatArray new_foot(2), new_normal(2), grid_p;
    double distance2, tubeWidth2 = this->tubeWidth * this->tubeWidth;

    // Find the bounding box for the segment
    FloatArray x0, x1;
    this->getBoundingBox(x0, x1, c, r + this->tubeWidth);

    for ( ParticleGrid< ParticlePoint > :: iterator it = this->grid->beginAt(x0, x1); !it.end(); ++it ) {
        it.getGridPoint(grid_p);

        double v = atan2( grid_p(1) - c(1), grid_p(0) - c(0) );
        if ( v < v0 || v > v1 ) {
            continue;
        }

        new_normal(0) = cos(v);
        new_normal(1) = sin(v);
        new_foot(0) = c(0) + cos(v) * r;
        new_foot(1) = c(1) + sin(v) * r;

        distance2 = grid_p.distance_square(new_foot);

        if ( distance2 < tubeWidth2 ) {
            if ( ( point = it.getPoint() ) ) {
                if ( point->distance2 <= distance2 ) {
                    continue; // If distance already is smaller, then skip this.
                }
            }
            point = new ParticlePoint(new_foot, id, new_normal, distance2);
            it.setPoint(point);
        }
    }
}

void ParticleTopologyDescription :: addCorner(int id, const FloatArray &c, ParticleGrid< ParticlePoint > &grid)
{
    /// @todo Remove this structure.
    ParticlePoint corner;
    corner.id = id;
    corner.foot = c;
    this->corners.push_front(corner);

    ParticlePoint *point;
    // Find the bounding box for the segment
    FloatArray x0, x1;
    this->getBoundingBox(x0, x1, c, this->tubeWidth);

    for ( ParticleGrid< ParticlePoint > :: iterator it = this->grid->beginAt(x0, x1); !it.end(); ++it ) {
        if ( ( point = it.getPoint() ) ) {
            if ( point->id == id ) {
                if ( point->corner.giveSize() < 0 || point->foot.distance_square(point->corner) > point->foot.distance_square(c) ) {
                    point->corner = c;
                }
            }
        }
    }
}

TopologyState ParticleTopologyDescription :: updateYourself(TimeStep *tStep)
{
    Timer timer;
    timer.startTimer();

    TopologyState ts = TS_OK;

    FloatArray displacement, grid_coord;
    ParticlePoint *p;
    this->maxdisp2 = 0; // Determines whether grid should be resampled.

    // Update all foot points
    this->d->giveSpatialLocalizer()->init(true);    // For SL we need to update every time step.

    this->resampled = false;
    for ( ParticleGrid< ParticlePoint > :: iterator it = this->grid->begin(); !it.end(); ++it ) {
        if ( ( p = it.getPoint() ) != NULL ) {
            bool state = this->findDisplacement(displacement, p->id, p->foot, tStep);
            if ( state ) {
                p->foot.add(displacement);
                p->total_displacement.add(displacement);
                this->maxdisp2 = max( this->maxdisp2, p->total_displacement.computeSquaredNorm() );

#if 0
                if ( p->corner.giveSize() > 0 ) {
                    this->findDisplacement(displacement, p->id, p->foot, tStep);
                    p->foot.add(displacement);
                    p->total_displacement.add(displacement);
                    this->maxdisp2 = max( this->maxdisp2, p->total_displacement.computeSquaredNorm() );
                }
#endif
            } else { // Can't find a suitable element to interpolate from, then just drop the point.
                grid->clearPosition( it.getIndex() );
            }
        }
    }

    // Remove this
#if 1
    for ( std :: list< ParticlePoint > :: iterator it = this->corners.begin(); it != this->corners.end(); ++it ) {
        bool state = this->findDisplacement(displacement, it->id, it->foot, tStep);
        if ( !state ) {
            OOFEM_ERROR("This always needs to find a displacement");
        }
        it->foot.add(displacement);
        it->total_displacement.add(displacement);
    }
#endif

    ts = this->checkOverlap();
    if ( ts != TS_OK || this->maxdisp2 > this->tubeWidth * this->tubeWidth ) {
        OOFEM_LOG_INFO("ParticleTopologyDescription :: updateYourself - Resampling\n");
        this->resample();
    }

    timer.stopTimer();
    OOFEM_LOG_INFO( "ParticleTopologyDescription: user time consumed by updating: %.2fs\n", timer.getUtime() );

    return ts;
}


TopologyState ParticleTopologyDescription :: checkOverlap()
{
    TopologyState ts = TS_OK;
    FloatArray x0, x1;
    std :: list< ParticlePoint * >points;
    ParticlePoint *p0;
    FloatArray p0p1;
    double normal_limit = 0.5, distance_limit = 0.0; // Very lenient requirement...
    distance_limit = this->tubeWidth / 2;

    for ( ParticleGrid< ParticlePoint > :: iterator it = this->grid->begin(); !it.end(); ++it ) {
        if ( ( p0 = it.getPoint() ) != NULL && p0->id > 0 ) {
            // Use either total distance or the local?
            this->getBoundingBox(x0, x1, p0->foot, this->tubeWidth + 2 * p0->total_displacement.computeNorm() + distance_limit);
            this->grid->getPointsWithin(points, x0, x1);

            std :: list< ParticlePoint * > :: iterator it_points;
            for ( it_points = points.begin(); it_points != points.end(); it_points++ ) {
                ParticlePoint *p1 = * it_points;
                if ( p1->id == p0->id ) { // For same surfaces, just merge them
                    if ( p1->normal.dotProduct(p0->normal) < normal_limit ) {
                        p0p1.beDifferenceOf(p1->foot, p0->foot);
                        // TODO: Is this a good criteria?
                        //if (p0->normal.dotProduct(p0p1) < distance_limit ) {
                        if ( p0p1.computeNorm() < distance_limit ) {
                            // Mark both for removal.
                            p0->removal = true;
                            p1->removal = true;
                            //ts = TS_NeedsRemeshing;
                        }
                    }
                }
                // This part turned out to be incredibly difficult; Merging of different surfaces (creating a new one, with end points accordingly).
#if 0
                else {
                    double distance = p0->foot.distance_square(p1->foot);
                    if ( distance < distance_limit ) {
                        if ( p0->corner.giveSize() > 0 && distance < p0->foot.distance_square(p0->corner) ) {
                            p0->removal = true;
                            p1->removal = true;
                        } else {
                            OOFEM_ERROR("Conditions for merging surface not implemented yet.");
                            // Not this simple, but something along these lines...
                            int newpid = ( int ) this->mergeID.at(p0->id, p1->id);
                            int control = ( int ) this->controlID.at(p0->id, p1->id);
                            if ( p0->id == control ) {
                                p0->id = newpid;
                                p1->removal = true;
                            } else {
                                p1->id = newpid;
                                p0->removal = true;
                            }
                            ts = TS_NeedsRemeshing;
                        }
                    }
                }
#endif
            }
        }
    }

    this->removePoints(*this->grid);
    return ts;
}

void ParticleTopologyDescription :: removePoints(ParticleGrid< ParticlePoint > &g) const
{
    ParticlePoint *p;
    FloatArray grid_coord;
    ParticleGrid< ParticlePoint > *sg;

    for ( int ind = 0; ind < g.getTotal(); ind++ ) {
        if ( g.getPoint(p, ind) ) {
            if ( p->removal ) {
                g.clearPosition(ind);
            }
        } else if ( this->grid->getSubGrid(sg, ind) ) {
            this->removePoints(*sg);
            // If all subpoint were removed, then remove the whole grid.
            if ( sg->isEmpty() ) {
                sg->clearPosition(ind);
            }
        }
    }
}

void ParticleTopologyDescription :: calculateShortestDistance(const ParticlePoint *p0, std :: list< ParticlePoint * > &neighbors, ParticleGrid< ParticlePoint > &g) const
// This routine is the most complicated one, and it's no unique in it's possible design.
// Current implementation does
// 1. Curve fit a second order polynomial
// 2. Calculate new shortest distances.
{
    int n = neighbors.size();
    FloatArray a; // Coefficients for the polynomial
    FloatArray tfooti; // Foot origin and mapped foot point.
    FloatArray temp;
    FloatArray normal;

    // The matrix representing the mapping to the normal space.
    FloatMatrix N;

    normal = p0->normal;
    for ( std :: list< ParticlePoint * > :: iterator it = neighbors.begin(); it != neighbors.end(); ++it ) {
        normal.add( ( * it )->normal );
    }
    normal.normalize();
    N.beLocalCoordSys(normal);

    // The limits of txi
    double txi_min = 0, txi_max = 0;

    FloatArray tx(n + 1), ty(n + 1);
    tx(0) = ty(0) = 0;       // Center point
    int i = 1;
    for ( std :: list< ParticlePoint * > :: iterator it = neighbors.begin(); it != neighbors.end(); ++it, ++i ) {
        // Calculate foot points in local coord.sys.
        temp.beDifferenceOf( ( * it )->foot, p0->foot );
        tfooti.beProductOf(N, temp);
        tx(i) = tfooti(0);
        ty(i) = tfooti(1);
        if ( tx(i) < txi_min ) {
            txi_min = tx(i);
        } else if ( tx(i) > txi_max ) {
            txi_max = tx(i);
        }
    }
    ls2fit(tx, ty, a);

    FloatArray x0, x1;
    IntArray ind0, ind1;
    double radius = this->tubeWidth + max(txi_max, -txi_min);
    this->getBoundingBox(x0, x1, p0->foot, radius);
    g.getBoundingBox(x0, x1, ind0, ind1);
    double tubeWidth2 = this->tubeWidth * this->tubeWidth;
    FloatArray new_foot, new_normal;

    double distance2;
    ParticlePoint *point;
    ParticleGrid< ParticlePoint > *subgrid;
    FloatArray grid_point;
    IntArray pos(2);
    /// @todo Adaptivity
    for ( pos(0) = ind0(0); pos(0) < ind1(0); pos(0)++ ) {
        for ( pos(1) = ind0(1); pos(1) < ind1(1); pos(1)++ ) {
            g.getGridCoord(grid_point, pos);
            if ( grid_point.distance_square(p0->foot) > radius * radius ) {  // Quick optimization (might give false negatives, but that is acceptable.)
                continue;
            }

            // Find the closest foot point
            distance2 = this->shortestDistanceFromCurve(a, txi_min, txi_max, normal, p0->foot,
                                                        grid_point, new_foot, new_normal);

            if ( distance2 < tubeWidth2 ) { // Otherwise we ignore it
                if ( g.getPoint(point, pos) ) {
                    if ( point->distance2 < distance2 ) {
                        continue;
                    }
                } else if ( g.getSubGrid(subgrid, pos) ) {
                    OOFEM_ERROR("Not implemented");
                }
                point = new ParticlePoint(new_foot, p0->id, new_normal, distance2);
                g.setPoint(point, pos);
            }
        }
    }
}

double ParticleTopologyDescription :: shortestDistanceFromCurve(const FloatArray &a, double txi_min, double txi_max,
                                                                const FloatArray &n0, const FloatArray &y0, const FloatArray &p, FloatArray &foot, FloatArray &normal) const
{
    double b3, b2, b1, b0; // Coefficients for the polynomial
    double temp1, temp2;
    FloatArray p_y0(2), t0(2);

    // Tangent
    t0(0) = n0(1);
    t0(1) = -n0(0);

    p_y0(0) = y0(0) - p(0);
    p_y0(1) = y0(1) - p(1);

    p_y0.beDifferenceOf(y0, p);
    temp1 = n0(0) * p_y0(0) + n0(1) * p_y0(1);           //n0.dotProduct(p_y0);
    temp2 = t0(0) * p_y0(0) + t0(1) * p_y0(1);           //t0.dotProduct(p_y0);

    // Same as [1]
    b3 = 2 * a(2) * a(2);
    b2 = 3 * a(1) * a(2);
    // Last terms with "tempX" missing in article [1]
    b1 = a(1) * a(1) + 2 * a(0) * a(2) + 1 + 2 * a(2) * temp1;
    b0 = a(0) * a(1) + a(1) * temp1 + temp2;

    // Find the roots
    int roots;
    double r [ 3 ] = {
        0, 0, 0
    };
    cubic(b3, b2, b1, b0, & r [ 0 ], & r [ 1 ], & r [ 2 ], & roots);
    // The potential points
    double limit = 0.6; // Found to be absolutely necessary for stability of surfaces when corner cases are included.
    double txi_list [ 5 ] = {
        txi_min *limit, txi_max * limit
    };
    int points = 0;
    for ( int i = 0; i < roots; i++ ) {
        if ( txi_min * limit < r [ i ] && r [ i ] < txi_max * limit ) {
            txi_list [ points ] = r [ i ];
            points++;
        }
    }

    // Normal and foot point in mapped domain (thus the t-prefix)
    double min_distance2 = 0.0, min_txi = 0.0;
    double dx, dy, distance2, f, fp, txi;
    for ( int i = 0; i < points; i++ ) {
        txi = txi_list [ i ];
        f = a(0) + a(1) * txi + a(2) * txi * txi;

        dx = n0(0) * f + t0(0) * txi + p_y0(0);
        dy = n0(1) * f + t0(1) * txi + p_y0(1);

        distance2 = dx * dx + dy * dy;

        if ( i == 0 || distance2 < min_distance2 ) {
            min_distance2 = distance2;
            min_txi = txi;
        }
    }

    if ( min_txi < txi_min || min_txi > txi_max ) {
        return 2 * this->tubeWidth * this->tubeWidth;
    }

    f = a(0) + a(1) * min_txi + a(2) * min_txi * min_txi;
    fp = a(1) + 2 * a(2) * min_txi;

    foot.resize(2);
    foot(0) = n0(0) * f + t0(0) * min_txi + y0(0);
    foot(1) = n0(1) * f + t0(1) * min_txi + y0(1);

    normal.resize(2);
    normal(0) = n0(0) + t0(0) * ( -fp );
    normal(1) = n0(1) + t0(1) * ( -fp );
    normal.normalize();

    return foot.distance_square(p);
}

void ParticleTopologyDescription :: resample()
// Copy structure from old grid for a new resampled grid.
{
    if ( this->resampled ) {
        return;
    }
    std :: list< ParticlePoint * >neighbours;
    ParticlePoint *origin;
    FloatArray grid0;
    FloatArray new_foot, new_normal;
    std :: unique_ptr< ParticleGrid< ParticlePoint > > new_grid( new ParticleGrid< ParticlePoint >(this->grid.get()) );

    double maxdisp = sqrt(this->maxdisp2);
    int total_points = 0;
    for ( ParticleGrid< ParticlePoint > :: iterator it = this->grid->begin(); !it.end(); ++it ) {
        total_points++;
        if ( ( origin = it.getPoint() ) != NULL ) {
            this->collectNeighbors(neighbours, origin, maxdisp);
            this->calculateShortestDistance(origin, neighbours, *new_grid);
        }
    }

    this->grid = std :: move( new_grid );
    this->resampled = true;
}

void ParticleTopologyDescription :: getBoundingBox(FloatArray &x0, FloatArray &x1, const FloatArray &c, double width)
{
    x0 = c;
    x0.add(-width);
    x1 = c;
    x1.add(width);
}

// Helper to sort the pairs
static bool compare_second(std :: pair< ParticlePoint *, double >a, std :: pair< ParticlePoint *, double >b)
{
    return fabs(a.second) < fabs(b.second);
}

void ParticleTopologyDescription :: collectNeighbors(std :: list< ParticlePoint * > &answer, const ParticlePoint *origin, double dist) const
{
    answer.clear();

    FloatArray temp;
    FloatArray lcoord;
    FloatMatrix toLocalCoord;

    toLocalCoord.beLocalCoordSys(origin->normal);

    // The algorithm uses a large bounding box, which should incorporate at least n points for well formed input.
    std :: list< ParticlePoint * >points;

    // Bounding box around grid point and collect all the potential points.
    FloatArray x0, x1;
    this->getBoundingBox(x0, x1, origin->foot, 1.1 * this->tubeWidth + dist);
    this->grid->getPointsWithin(points, x0, x1);

    typedef std :: pair< ParticlePoint *, double >dist_pair;
    std :: list< dist_pair >valid_points;

    /// @todo Change to the local corner information.
    for ( std :: list< ParticlePoint > :: const_iterator it_corners = corners.begin(); it_corners != this->corners.end(); ++it_corners ) {
        ParticlePoint *p = ( ParticlePoint * ) & ( * it_corners );
        if ( p->id != origin->id ) {
            continue;
        }

        if ( origin->foot.distance_square(p->foot) > 4 * this->tubeWidth * this->tubeWidth ) {
            continue;
        }
        temp.beDifferenceOf(p->foot, origin->foot);
        lcoord.beProductOf(toLocalCoord, temp);
        valid_points.push_front( dist_pair( p, lcoord(0) ) );
    }

    // Apply rejection criteria based on normals (and other things?)
    std :: list< ParticlePoint * > :: iterator it_points;
    for ( it_points = points.begin(); it_points != points.end(); ++it_points ) {
        ParticlePoint *p = * it_points;
        if ( p->id != origin->id ) {
            continue;
        }

        double projection = p->normal.dotProduct(origin->normal);
        if ( projection < 0.5 ) {
            continue;
        }
        temp.beDifferenceOf(p->foot, origin->foot);
        lcoord.beProductOf(toLocalCoord, temp);
        valid_points.push_front( dist_pair( p, lcoord(0) ) );
    }

    valid_points.sort(compare_second);
    //valid_points.unique(same_point);

    // Some special effort to take neighbors from both sides as well as
    int m_smaller = 0, m_larger = 0;
    double limit = 0.01 * this->grid->getGridStep(0);
    double smaller = -limit, larger = limit;
    // Copy them over to answer
    std :: list< dist_pair > :: iterator it = valid_points.begin();
    for ( int i = 0; i < ( int ) valid_points.size(); i++, it++ ) {
        if ( m_larger > 0 && m_smaller > 0 && m_larger + m_smaller >= this->m ) {
            break;
        }
        if ( it->second > larger ) {
            if ( m_larger >= this->m ) {
                continue;
            } else {
                larger = it->second + limit;
                m_larger++;
            }
        } else if ( it->second < smaller ) {
            if ( m_smaller >= this->m ) {
                continue;
            } else {
                smaller = it->second - limit;
                m_smaller++;
            }
        } else { // To close to center, skip it.
            continue;
        }
        answer.push_front(it->first);
    }
}

bool ParticleTopologyDescription :: findDisplacement(FloatArray &displacement, int id, const FloatArray &footpoint, TimeStep *tStep) const
{
    FloatArray lcoords, closest, offset;
    IntArray region(1);
    region(0) = id;
    // Finds a single element close to this point.
    Element *e = this->d->giveSpatialLocalizer()->giveElementClosestToPoint(lcoords, closest, footpoint, id);
    if ( !e ) { // No element at all for this domain. Then there is no displacement;
        OOFEM_WARNING("Couldn't find any element to interpolate from");
        displacement.resize( footpoint.giveSize() );
        displacement.zero();
        return false;
    }
    offset.beDifferenceOf(closest, footpoint);
    double dist = offset.computeSquaredNorm();
    if ( dist > 0.25 * this->grid->getGridStep(0) * this->grid->getGridStep(0) ) {    // TODO: Maybe also check normals?
        OOFEM_WARNING( "Foot point far from element (%e) (%e, %e)->(%e, %e).\n", sqrt(dist), footpoint.at(1), footpoint.at(2), closest.at(1), closest.at(2) );
        return false;
    }


    if ( this->useDisplacements ) {
        e->computeField(VM_Incremental, tStep, lcoords, displacement);    // Displacement
    } else {
        FloatArray fields;
        IntArray dofIds;
        displacement.resize( this->d->giveNumberOfSpatialDimensions() );
        displacement.zero();
        double dt = tStep->giveTimeIncrement();
        e->computeField(VM_Total, tStep, lcoords, fields);    // Velocities + pressures most likely.
        e->giveElementDofIDMask(dofIds);
        for ( int i = 0; i < dofIds.giveSize(); i++ ) {
            if ( dofIds(i) == V_u ) {
                displacement(0) = fields(i) * dt;
            } else if ( dofIds(i) == V_v ) {
                displacement(1) = fields(i) * dt;
            } else if ( dofIds(i) == V_w ) {
                displacement(2) = fields(i) * dt;
            }
        }
    }
    displacement.add(offset);
    return true;
}

void ParticleTopologyDescription :: doOutput(TimeStep *tStep)
{
    std :: ostringstream name;
    name << this->d->giveEngngModel()->giveOutputBaseFileName() << "." << tStep->giveNumber() << ".vtp";
    if ( writeVTK ) {
        this->writeVTKFile( name.str().c_str() );
    }
}


void ParticleTopologyDescription :: writeDataToFile(const char *name) const
{
    // Dumps everything in a ASCII table.
    FILE *fid = fopen(name, "w");

    int dims;
    FloatArray grid_coord;
    ParticlePoint *p;

    for ( ParticleGrid< ParticlePoint > :: iterator it = this->grid->begin(); !it.end(); ++it ) {
        if ( ( p = it.getPoint() ) != NULL ) {
            it.getGridPoint(grid_coord);

            dims = grid_coord.giveSize();
            // First write the grid coordinate
            for ( int i = 0; i < dims; i++ ) {
                fprintf( fid, "%e ", grid_coord(i) );
            }
            // Then the foot point
            for ( int i = 0; i < dims; i++ ) {
                fprintf( fid, "%e ", p->foot(i) );
            }
            // Auxiliary information, like normal
            for ( int i = 0; i < dims; i++ ) {
                fprintf( fid, "%e ", p->normal(i) );
            }
            fprintf(fid, "%d ", p->id);
            fprintf(fid, "\n");
        }
    }

    fclose(fid);
}

void ParticleTopologyDescription :: writeVTKFile(const char *name) const
{
#ifdef __VTK_MODULE
    ParticlePoint *p;
    int dims = this->grid->giveDimensions();
    // Create points and normals
    int npoints = this->grid->getNumberOfPoints();
    vtkSmartPointer< vtkPoints >points = vtkSmartPointer< vtkPoints > :: New();
    vtkSmartPointer< vtkVertex >vertex;
    vtkSmartPointer< vtkCellArray >vertices = vtkSmartPointer< vtkCellArray > :: New();
    vtkSmartPointer< vtkDoubleArray >pointNormalsArray = vtkSmartPointer< vtkDoubleArray > :: New();
    vtkSmartPointer< vtkIntArray >pointIDArray = vtkSmartPointer< vtkIntArray > :: New();
    vtkSmartPointer< vtkDoubleArray >pointEndPointsArray = vtkSmartPointer< vtkDoubleArray > :: New();

    pointIDArray->SetNumberOfComponents(1);
    pointNormalsArray->SetNumberOfComponents(3);
    pointEndPointsArray->SetNumberOfComponents(3);

    pointIDArray->SetName("ID");
    pointNormalsArray->SetName("Normals");
    pointEndPointsArray->SetName("End points");

    pointIDArray->SetNumberOfTuples(npoints);
    pointNormalsArray->SetNumberOfTuples(npoints);
    pointEndPointsArray->SetNumberOfTuples(npoints);

    int i = 0;
    for ( ParticleGrid< ParticlePoint > :: iterator it = this->grid->begin(); !it.end(); ++it ) {
        if ( ( p = it.getPoint() ) != NULL ) {
            points->InsertNextPoint(p->foot(0), p->foot(1), dims == 3 ? p->foot(2) : 0.0);
            vertex = vtkSmartPointer< vtkVertex > :: New();
            vertex->GetPointIds()->SetId(0, i);
            vertices->InsertNextCell(vertex);
            pointIDArray->SetTuple1(i,  p->id);
            pointNormalsArray->SetTuple3(i, p->normal(0), p->normal(1), dims == 3 ? p->normal(2) : 0.0);
            if ( p->corner.giveSize() > 0 ) {
                pointEndPointsArray->InsertNextTuple3(p->corner(0) - p->foot(0), p->corner(1) - p->foot(1), dims == 3 ? p->corner(2) - p->foot(2) : 0.0);
            } else {
                pointEndPointsArray->InsertNextTuple3(0.0, 0.0, 0.0);
            }
            ++i;
        }
    }

    // Create a polydata object and add the points to it.
    vtkSmartPointer< vtkPolyData >polydata = vtkSmartPointer< vtkPolyData > :: New();
    polydata->SetPoints(points);
    polydata->SetVerts(vertices);
    polydata->GetPointData()->SetNormals(pointNormalsArray);
    polydata->GetPointData()->SetVectors(pointEndPointsArray);
    polydata->GetPointData()->SetScalars(pointIDArray);

    // Write the file
    vtkSmartPointer< vtkXMLPolyDataWriter >writer = vtkSmartPointer< vtkXMLPolyDataWriter > :: New();
    writer->SetFileName(name);
    writer->SetInput(polydata);
    // Optional - set the mode. The default is binary.
    //writer->SetDataModeToBinary();
    //writer->SetDataModeToAscii();
    writer->Write();
#else
    OOFEM_WARNING("Not compiled with VTK support.");
#endif
}

class edge
{
public:
    int first;
    int second;
    int id;
};

class node
{
public:
    FloatArray c;
    int id;
};

static bool sort_edge(edge a, edge b)
{
    return ( a.first == b.first ) ? a.second > b.second : a.first > b.first;
}

static bool compare_edge(edge a, edge b)
{
    return a.first == b.first && a.second == b.second;
}

// There are much room for improvement here. In particular, it should be similar to
void ParticleTopologyDescription :: generatePSLG(Triangle_PSLG &pslg)
{
    this->resample();

    FloatArray lcoord, temp;
    FloatMatrix toLocalCoord;

    // First enumerate all the nodes, removing duplicates
    ParticlePoint *origin;
    std :: list< node >nodes;
    std :: list< ParticlePoint * >points;
    FloatArray x0, x1, grid_coord;

    // Nodal merge limit (for very close points, 1/10 of the grid step).
    double merge2 = this->grid->getGridStep(1) * this->grid->getGridStep(1) * 1e-2;
    // PSLG simplification limit, 1/10 of the grid step.
    double limit = this->grid->getGridStep(1) * 0.5;

    int total_nodes = 0;

    for ( ParticleGrid< ParticlePoint > :: iterator it = this->grid->begin(); !it.end(); ++it ) {
        if ( ( origin = it.getPoint() ) != NULL ) {
            origin->node = 0;
        }
    }

    for ( std :: list< ParticlePoint > :: iterator cit = this->corners.begin(); cit != this->corners.end(); ++cit ) {
        cit->node = 0;
    }

    // Enumerate the corner nodes first;
    for ( std :: list< ParticlePoint > :: iterator cit = this->corners.begin(); cit != this->corners.end(); ++cit ) {
        if ( cit->node != 0 ) {
            continue;
        }
        total_nodes++;
        cit->node = total_nodes;
        node n;
        n.id = cit->id;
        n.c = cit->foot;
        nodes.push_back(n);
        for ( std :: list< ParticlePoint > :: iterator cit_ = this->corners.begin(); cit_ != this->corners.end(); ++cit_ ) {
            if ( cit->foot.distance_square(cit_->foot) <= merge2 && cit_->node == 0 ) {
                cit_->node = cit->node;
                cit_->foot = cit->foot; // For consistency, merge the points which are merged in the mesh
            }
        }

        this->getBoundingBox(x0, x1, cit->foot, 2 * this->tubeWidth);
        this->grid->getPointsWithin(points, x0, x1);
        for ( auto pit = points.begin(); pit != points.end(); pit++ ) {
            if ( ( cit->foot.distance_square( ( * pit )->foot ) <= merge2 ) && ( * pit )->node == 0 ) {
                ( * pit )->node = cit->node;
                ( * pit )->foot = cit->foot;
            }
        }
    }

    for ( ParticleGrid< ParticlePoint > :: iterator it = this->grid->begin(); !it.end(); ++it ) {
        if ( ( origin = it.getPoint() ) != NULL ) {
            if ( origin->node != 0 ) {
                continue;
            }
            toLocalCoord.beLocalCoordSys(origin->normal);
            total_nodes++;
            origin->node = total_nodes;
            node n;
            n.id = origin->id;
            n.c = origin->foot;
            nodes.push_back(n);

            this->getBoundingBox(x0, x1, origin->foot, 2 * this->tubeWidth);
            this->grid->getPointsWithin(points, x0, x1);
            for ( auto pit = points.begin(); pit != points.end(); pit++ ) {
                if ( ( * pit )->node == 0 ) {
                    /*
                     * temp.beDifferenceOf((*pit)->foot, origin->foot);
                     * lcoord.beProductOf(toLocalCoord, temp);
                     * double ln = lcoord.at(lcoord.giveSize()); // Treating the buggy offset differently.
                     * lcoord.at(lcoord.giveSize()) = 0.0;
                     * if ( lcoord.computeSquaredNorm() <= limit*limit && (ln < this->tubeWidth && origin->normal.dotProduct((*pit)->normal) ) {
                     */
                    if ( origin->foot.distance_square( ( * pit )->foot ) <= merge2 ) {
                        ( * pit )->node = origin->node;
                        ( * pit )->foot = origin->foot;
                    }
                }
            }
        }
    }

    // Find segments
    std :: list< edge >edges;
    double max_value = 2 * this->tubeWidth; // No connections further away than this.
    double front_dist, back_dist, b_front_dist, b_back_dist;
    for ( ParticleGrid< ParticlePoint > :: iterator it = this->grid->begin(); !it.end(); ++it ) {
        if ( ( origin = it.getPoint() ) != NULL ) {
            toLocalCoord.beLocalCoordSys(origin->normal);

            edge e0, e1;
            e0.id = origin->id;
            e1.id = origin->id;
            e0.second = origin->node;
            e1.first = origin->node;
            e0.first = 0;
            e1.second = 0;
            front_dist = max_value;
            back_dist = -max_value;

            bool found_front = false, found_back = false;

            // Backup values
            b_front_dist = -max_value;
            b_back_dist = max_value;

            this->getBoundingBox(x0, x1, origin->foot, max_value);
            this->grid->getPointsWithin(points, x0, x1);

            // Loop over corner nodes as well
            for ( std :: list< ParticlePoint > :: iterator cit = this->corners.begin(); cit != this->corners.end(); ++cit ) {
                points.push_front(& * cit);
            }

            for ( auto pit = points.begin(); pit != points.end(); ++pit ) {
                ParticlePoint *p = * pit;

                if ( p->id != origin->id || origin->node == p->node ) {
                    continue;
                }

                temp.beDifferenceOf(p->foot, origin->foot);

                if ( temp.computeSquaredNorm() > max_value * max_value ) {
                    continue;
                }

                lcoord.beProductOf(toLocalCoord, temp);

                double projection = p->normal.giveSize() > 0 ? p->normal.dotProduct(origin->normal) : 1.0;
                if ( projection < 0.0 ) { // Then take the furthest point
#if 1
                    if ( lcoord(0) > b_front_dist && !found_front ) {
                        b_front_dist = lcoord(0);
                        e1.second = p->node;
                    }
                    if ( lcoord(0) < b_back_dist && !found_back ) {
                        b_back_dist = lcoord(0);
                        e0.first = p->node;
                    }
#endif
                    continue;
                }

                if ( lcoord(0) > 0 ) {
                    if ( lcoord(0) < front_dist ) {
                        found_front = true;
                        front_dist = lcoord(0);
                        e1.second = p->node;
                    }
                } else {
                    if ( lcoord(0) > back_dist ) {
                        found_back = true;
                        back_dist = lcoord(0);
                        e0.first = p->node;
                    }
                }
            }
            if ( e0.first != 0 ) {
                edges.push_front(e0);
            }
            if ( e1.second != 0 ) {
                edges.push_front(e1);
            }
        }
    }
    edges.sort(sort_edge);
    edges.unique(compare_edge);

    Triangle_PSLG fine;
    fine.nx.resize( nodes.size() );
    fine.ny.resize( nodes.size() );
    std :: list< node > :: iterator it_node;
    int n = 0;
    for ( it_node = nodes.begin(); it_node != nodes.end(); it_node++ ) {
        fine.nx(n) = it_node->c(0);
        fine.ny(n) = it_node->c(1);
        n++;
    }

    fine.segment_a.resize( edges.size() );
    fine.segment_b.resize( edges.size() );
    fine.segment_marker.resize( edges.size() );
    std :: list< edge > :: iterator it_edge;
    n = 0;
    for ( it_edge = edges.begin(); it_edge != edges.end(); it_edge++ ) {
        fine.segment_a(n) = it_edge->first;
        fine.segment_b(n) = it_edge->second;
        fine.segment_marker(n) = it_edge->id;
        n++;
    }

#if 0
    pslg.nx = fine.nx;
    pslg.ny = fine.ny;
    pslg.segment_a = fine.segment_a;
    pslg.segment_b = fine.segment_b;
    pslg.segment_marker = fine.segment_marker;
#else
    TriangleMesherInterface :: simplifyPSLG(pslg, fine, limit);
#endif

    // DEBUG: Print all to matlab file
#if 1
    this->writeVTKFile("initial_particle.vtp");

    FILE *file = fopen("pslg.m", "w");

    fprintf(file, "edges = [");
    for ( int i = 0; i < pslg.segment_a.giveSize(); i++ ) {
        fprintf( file, "%d, %d, %d;\n", pslg.segment_a(i), pslg.segment_b(i), pslg.segment_marker(i) );
    }
    fprintf(file, "];\n");

    fprintf(file, "nodes = [");
    for ( int i = 0; i < pslg.nx.giveSize(); i++ ) {
        fprintf( file, "%e, %e;\n", pslg.nx(i), pslg.ny(i) );
    }
    fprintf(file, "];\n");

    fprintf(file, "colors = 'grbkmcygrbkmcygrbkmcy';\n");
    fprintf(file, "hold on, axis equal\n");
    fprintf(file, "for i = 1:size(nodes,1)\n  plot(nodes(i,1),nodes(i,2),'rx'); %%text(nodes(i,1),nodes(i,2),num2str(i)); \nend\n");
    //fprintf(file, "for i = 1:size(edges,1)\n  plot(nodes(edges(i,1:2),1),nodes(edges(i,1:2),2),colors(edges(i,3)));\nend\n");
    fprintf(file, "for i = 1:size(edges,1)\n  quiver(nodes(edges(i,1),1),nodes(edges(i,1),2),nodes(edges(i,2),1)-nodes(edges(i,1),1),nodes(edges(i,2),2)-nodes(edges(i,1),2),1,colors(edges(i,3)));\nend\n");
    fprintf(file, "for i = 1:size(edges,1)\n  %%text(mean(nodes(edges(i,1:2),1)),mean(nodes(edges(i,1:2),2)),['\\color{red} ',num2str(i)]);\nend\n");

    fclose(file);
#endif
}

void ParticleTopologyDescription :: generateMesh(std :: vector< FloatArray > &nodes, std :: vector< IntArray > &elements, std :: vector< IntArray > &segments,
                                                 std :: vector< IntArray > &n_markers, IntArray &e_markers, IntArray &s_markers, IntArray &e_egt, IntArray &s_egt)
{
    double maxArea = 0.05;
    IntArray r_markers;

    TriangleMesherInterface tmi(30, maxArea, true);

    Triangle_PSLG pslg;
    this->generatePSLG(pslg);
    tmi.meshPSLG(pslg,
                 //holes, regions, r_markers,
                 this->regionOutside, this->regionInside,   // Input
                 nodes, n_markers, elements, e_markers, segments, s_markers);    // Output

    // Append the corner nodes manually. Didn't see any pretty way to do this.
#if 0 // Replace by this...
    for ( ParticleGrid< ParticlePoint > :: iterator it = this->grid->begin(); !it.end(); ++it ) {
        ParticlePoint *origin;
        if ( ( origin = it.getPoint() ) != NULL ) {
            double dist2, min_dist2 = 1e100;
            int min_i = 0;
            for ( int i = 1; i <= ( int ) nodes.size(); ++i ) {
                dist2 = nodes [ i - 1 ].distance_square(point->foot);
                if ( dist2 < min_dist2 ) {
                    min_dist2 = dist2;
                    min_i = i;
                }
            }
            n_markers.at(min_i)->insertSortedOnce(point->id);
        }
    }
#else
    for ( std :: list< ParticlePoint > :: iterator cit = this->corners.begin(); cit != this->corners.end(); ++cit ) {
        double dist2, min_dist2 = 1e100;
        int min_i = 0;
        for ( int i = 1; i <= ( int ) nodes.size(); ++i ) {
            dist2 = nodes [ i - 1 ].distance_square(cit->foot);
            if ( dist2 < min_dist2 ) {
                min_dist2 = dist2;
                min_i = i;
            }
        }
        n_markers [ min_i - 1 ].insertSortedOnce(cit->id);
    }
#endif

#if 1 // DEBUG: Print all to matlab file
    printf("Printing mesh.m\n");
    FILE *file = fopen("mesh.m", "w");
    fprintf(file, "nodes = [");
    for ( int i = 1; i <= ( int ) nodes.size(); i++ ) {
        FloatArray &x = nodes [ i - 1 ];
        for ( int j = 1; j <= 2; j++ ) {
            fprintf( file, "%e, ", x.at(j) );
        }
        fprintf( file, " %d", n_markers [ i - 1 ].at(1) );
        fprintf(file, ";\n");
    }
    fprintf(file, "];\n");

    fprintf(file, "triangles = [");
    for ( int i = 1; i <= ( int ) elements.size(); i++ ) {
        IntArray &x = elements [ i - 1 ];
        for ( int j = 1; j <= x.giveSize(); j++ ) {
            fprintf( file, "%d, ", x.at(j) );
        }
        fprintf( file, "%d", e_markers.at(i) );
        fprintf(file, ";\n");
    }
    fprintf(file, "];\n");

    fprintf(file, "segments = [");
    for ( int i = 1; i <= ( int ) segments.size(); i++ ) {
        IntArray &x = segments [ i - 1 ];
        for ( int j = 1; j <= x.giveSize(); j++ ) {
            fprintf( file, "%d, ", x.at(j) );
        }
        fprintf( file, "%d", s_markers.at(i) );
        fprintf(file, ";\n");
    }
    fprintf(file, "];\n");

    //fprintf(file, "clf\nhold on\nfor i = 1:size(triangles,1); plot(nodes(triangles(i,[1:3,1]),1), nodes(triangles(i,[1:3,1]),2), 'r'); end\n");
    fprintf(file, "clf\nhold on\n xx = nodes(:,1); yy = nodes(:,2); patch(xx(triangles(:,1:3))',yy(triangles(:,1:3))',triangles(:,end)'); \n");
    fprintf(file, "for i = 1:size(segments,1); plot(nodes(segments(i,1:2),1), nodes(segments(i,1:2),2), 'b'); end\n");

    fclose(file);
#endif
}

void ParticleTopologyDescription :: replaceFEMesh()
{
    std :: vector< FloatArray >nodes;
    std :: vector< IntArray >elements, segments, n_markers;
    ///@todo Boundary sets vs. bulk sets?
    std :: vector< IntArray >boundarySets( regionSet.maximum() );
    //std::vector< IntArray >buldSets( regionSet.maximum() );
    IntArray e_markers, s_markers, e_egt, s_egt;
    this->generateMesh(nodes, elements, segments, n_markers, e_markers, s_markers, e_egt, s_egt);

    for ( int i = 1; i <= ( int ) elements.size(); ++i ) {
        int r = e_markers.at(i);
        int set = regionSet.at(r);
        if ( set ) {
            boundarySets [ set - 1 ].followedBy(i);
            boundarySets [ set - 1 ].followedBy(0);
        }
    }
    for ( int i = 1; i <= ( int ) segments.size(); ++i ) {
        int r = s_markers.at(i);
        int set = regionSet.at(r);
        if ( set ) {
            boundarySets [ set - 1 ].followedBy( i + elements.size() );
            boundarySets [ set - 1 ].followedBy(0);
        }
    }

    // Replace all the elements with the new set.
    Domain *new_d;
    // Should i use a new domain and replace it instead?
    new_d = this->d;

    // Clear everything old;
    new_d->resizeElements( elements.size() + segments.size() );
    new_d->resizeDofManagers( nodes.size() );
    new_d->resizeSets( regionSet.giveSize() );

    for ( int i = 1; i <= ( int ) nodes.size(); ++i ) {
        //DynamicInputRecord *dir = new DynamicInputRecord(_IFT_Node_Name, i);
        //dir->setField(nodes[i-1], _IFT_Node_coords);
        //dr->insertInputRecord(IR_dofmanRec, dir);
        Node *node = new Node(i, new_d);
        node->setCoordinates(nodes [ i - 1 ]);
        new_d->setDofManager(i, node);
    }

    for ( int i = 1; i <= ( int ) elements.size(); ++i ) {
        int r = e_markers.at(i);
        //DynamicInputRecord *dir = new DynamicInputRecord(regionElementType[ r-1 ], i);
        //dir->setField(r, _IFT_Element_crosssect);
        //dir->setField(elements[i-1], _IFT_Element_nodes);
        //dr->insertInputRecord(IR_elemRec, dir);
        Element *element = classFactory.createElement(this->regionElementType [ r - 1 ].c_str(), i, new_d);
        element->setGlobalNumber(i);
        element->setCrossSection(r);
        element->setDofManagers(elements [ i - 1 ]);
        element->setParallelMode(Element_local);
        new_d->setElement(element->giveNumber(), element);
    }

    for ( int i = 1; i <= ( int ) segments.size(); ++i ) {
        int r = s_markers.at(i);
        //DynamicInputRecord *dir = new DynamicInputRecord(regionElementType[ r-1 ], segments.size() + i);
        //dir->setField(r, _IFT_Element_crosssect);
        //dir->setField(segments[i-1], _IFT_Element_nodes);
        //dr->insertInputRecord(IR_elemRec, dir);
        Element *element = classFactory.createElement(regionElementType [ r - 1 ].c_str(), elements.size() + i, new_d);
        if ( !element ) {
            printf( "Couldn't create: %s\n", regionElementType [ r - 1 ].c_str() );
        }
        element->setGlobalNumber( i + elements.size() );
        element->setCrossSection(r);
        element->setDofManagers(segments [ i - 1 ]);
        element->setParallelMode(Element_local);
        new_d->setElement(element->giveNumber(), element);
    }

    for ( int i = 1; i <= ( int ) boundarySets.size(); ++i ) {
        //DynamicInputRecord *dir = new DynamicInputRecord(_IFT_Set_Name, i);
        //dir->setField(boundarySets[i-1], _IFT_Set_elementBoundaries);
        //dr->insertInputRecord(IR_setRec, dir);
        Set *set = new Set(1, new_d);
        set->setBoundaryList(boundarySets [ i - 1 ]);
        new_d->setSet(i, set);
    }

    new_d->postInitialize();
    new_d->giveEngngModel()->checkConsistency();
    new_d->giveEngngModel()->forceEquationNumbering(); // TODO: This shouldn't be necessary
}
}
