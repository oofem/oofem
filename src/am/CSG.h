#pragma once

// Original CSG.JS library by Evan Wallace (http://madebyevan.com), under the MIT license.
// GitHub: https://github.com/evanw/csg.js/
//
// C++ port by Tomasz Dabrowski (http://28byteslater.com), under the MIT license.
// GitHub: https://github.com/dabroz/csgjscpp-cpp/
//
// Additional modifications by Jan Vorisek (https://github.com/janvorisek/)
// Added support for STL and VTK export, volume calculation
//
// Constructive Solid Geometry (CSG) is a modeling technique that uses Boolean
// operations like union and intersection to combine 3D solids. This library
// implements CSG operations on meshes elegantly and concisely using BSP trees,
// and is meant to serve as an easily understandable implementation of the
// algorithm. All edge cases involving overlapping coplanar polygons in both
// solids are correctly handled.

#include <algorithm>
#include <memory>
#include <math.h>
#include <vector>
#include <deque>
#include <map>

// `CSG.Plane.EPSILON` is the tolerance used by `splitPolygon()` to decide if a
// point is on the plane.
const double csgjs_EPSILON = 0.0001;

struct Vector {
    double x, y, z;

    Vector() :
        x( 0.0 ), y( 0.0 ), z( 0.0 )
    {
    }
    Vector( double x, double y, double z ) :
        x( x ), y( y ), z( z )
    {
    }
};

inline bool approxequal( double a, double b )
{
    return fabs( a - b ) < csgjs_EPSILON;
}

inline bool operator==( const Vector &a, const Vector &b )
{
    return approxequal( a.x, b.x ) && approxequal( a.y, b.y ) && approxequal( a.z, b.z );
}

inline bool operator!=( const Vector &a, const Vector &b )
{
    return !approxequal( a.x, b.x ) || !approxequal( a.y, b.y ) || !approxequal( a.z, b.z );
}

// Vector implementation

inline Vector operator+( const Vector &a, const Vector &b )
{
    return Vector( a.x + b.x, a.y + b.y, a.z + b.z );
}
inline Vector operator-( const Vector &a, const Vector &b )
{
    return Vector( a.x - b.x, a.y - b.y, a.z - b.z );
}
inline Vector operator*( const Vector &a, double b )
{
    return Vector( a.x * b, a.y * b, a.z * b );
}
inline Vector operator/( const Vector &a, double b )
{
    return a * ( (double)1.0 / b );
}
inline double dot( const Vector &a, const Vector &b )
{
    return a.x * b.x + a.y * b.y + a.z * b.z;
}
inline Vector lerp( const Vector &a, const Vector &b, double v )
{
    return a + ( b - a ) * v;
}
inline Vector negate( const Vector &a )
{
    return a * -(double)1.0;
}
inline double length( const Vector &a )
{
    return (double)sqrt( dot( a, a ) );
}

inline double lengthsquared( const Vector &a )
{
    return dot( a, a );
}

inline Vector unit( const Vector &a )
{
    return a / length( a );
}
inline Vector cross( const Vector &a, const Vector &b )
{
    return Vector( a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x );
}

inline Vector operator-( const Vector &a )
{
    return Vector( -a.x, -a.y, -a.z );
}

inline int lerp( int a, int b, double v )
{
    return a + (int)( ( b - a ) * v );
}

struct Vertex {
    Vector pos;
    Vector normal;
    int col;
};

inline bool operator==( const Vertex &a, const Vertex &b )
{
    return a.pos == b.pos && a.normal == b.normal && a.col == b.col;
}

inline bool operator!=( const Vertex &a, const Vertex &b )
{
    return a.pos != b.pos || a.normal != b.normal || a.col != b.col;
}

// interpolate vertex between two vertices, return t
inline double interpolate_between( const Vertex &a, const Vertex &b, const Vertex &v )
{
    return dot( v.pos - a.pos, b.pos - a.pos ) / lengthsquared( b.pos - a.pos );
}

struct Polygon;

// Represents a plane in 3D space.
struct Plane {
    Vector normal;
    double w;

    Plane();
    Plane( const Vector &a, const Vector &b, const Vector &c );

    inline bool ok() const
    {
        return length( this->normal ) > 0.0;
    }

    inline void flip()
    {
        this->normal = negate( this->normal );
        this->w *= -1.0;
    }

    void splitpolygon( const Polygon &poly, std::vector<Polygon> &coplanarFront,
        std::vector<Polygon> &coplanarBack, std::vector<Polygon> &front,
        std::vector<Polygon> &back ) const;

    enum Classification {
        COPLANAR = 0,
        FRONT    = 1,
        BACK     = 2,
        SPANNING = 3
    };
    inline Classification classify( const Vector &p ) const
    {
        double t         = dot( normal, p ) - this->w;
        Classification c = ( t < -csgjs_EPSILON ) ? BACK : ( ( t > csgjs_EPSILON ) ? FRONT : COPLANAR );
        return c;
    }
};

// Represents a convex polygon. The vertices used to initialize a polygon must
// be coplanar and form a convex loop. They do not have to be `CSG.Vertex`
// instances but they must behave similarly (duck typing can be used for
// customization).
//
// Each convex polygon has a `shared` property, which is shared between all
// polygons that are clones of each other or were split from the same polygon.
// This can be used to define per-polygon properties (such as surface color).
struct Polygon {
    std::vector<Vertex> vertices;
    Plane plane;

    Polygon();
    Polygon( const std::vector<Vertex> &list );

    inline void flip()
    {
        std::reverse( vertices.begin(), vertices.end() );
        for ( size_t i = 0; i < vertices.size(); i++ )
            vertices[i].normal = negate( vertices[i].normal );
        plane.flip();
    }
};

struct Model {

    using Index = int;

    std::vector<Vertex> vertices;
    std::vector<Index> indices;

    Index AddVertex( const Vertex &newv )
    {
        Index i = 0;
        for ( const auto &v : vertices ) {
            if ( v == newv ) {
                return i;
            }
            ++i;
        }
        vertices.push_back( newv );
        return i;
    }

    // Method to compute the volume of the mesh
    double volume() const
    {
        double volume = 0.0;
        for ( size_t i = 0; i < indices.size(); i += 3 ) {
            const Vertex &v1 = vertices[indices[i]];
            const Vertex &v2 = vertices[indices[i + 1]];
            const Vertex &v3 = vertices[indices[i + 2]];
            auto n           = cross( v2.pos - v1.pos, v3.pos - v1.pos );
            volume += dot( v1.pos, n );
        }
        return volume / 6.0;
    }
};

// Function to calculate the distance between a point and a line
double distanceToLine( const Vector &point, const Vector &linePoint1, const Vector &linePoint2 )
{
    // Calculate the components of the vector from linePoint1 to point
    double dx = point.x - linePoint1.x;
    double dy = point.y - linePoint1.y;
    double dz = point.z - linePoint1.z;

    // Calculate the components of the vector representing the line
    double line_dx = linePoint2.x - linePoint1.x;
    double line_dy = linePoint2.y - linePoint1.y;
    double line_dz = linePoint2.z - linePoint1.z;

    // Calculate the squared length of the line
    double line_length_squared = line_dx * line_dx + line_dy * line_dy + line_dz * line_dz;

    // Calculate the dot product of the two vectors
    double dot_product = dx * line_dx + dy * line_dy + dz * line_dz;

    // Calculate the parameter along the line where the point of interest lies
    double t = dot_product / line_length_squared;

    // Calculate the closest point on the line to the given point
    double closest_point_x = linePoint1.x + t * line_dx;
    double closest_point_y = linePoint1.y + t * line_dy;
    double closest_point_z = linePoint1.z + t * line_dz;

    // Calculate the distance between the point and its projection on the line
    double distance = sqrt( ( point.x - closest_point_x ) * ( point.x - closest_point_x ) + ( point.y - closest_point_y ) * ( point.y - closest_point_y ) + ( point.z - closest_point_z ) * ( point.z - closest_point_z ) );

    return distance;
}

void addMissingVertices( Model &mesh )
{
    std::vector<int> toDelete;
    // loop over triangles
    for ( int i = mesh.indices.size() - 3; i >= 0; i -= 3 ) {
        // get the three vertices of the triangle
        Vertex &v1 = mesh.vertices[mesh.indices[i]];
        Vertex &v2 = mesh.vertices[mesh.indices[i + 1]];
        Vertex &v3 = mesh.vertices[mesh.indices[i + 2]];

        // calculate the normal of the triangle
        Vector n = unit( cross( v2.pos - v1.pos, v3.pos - v1.pos ) );

        // Check other triangles that may touch triangle edges
        for ( int j = i - 3; j >= 0; j -= 3 ) {
            // get the three vertices of the triangle
            Vertex &v4 = mesh.vertices[mesh.indices[j]];
            Vertex &v5 = mesh.vertices[mesh.indices[j + 1]];
            Vertex &v6 = mesh.vertices[mesh.indices[j + 2]];

            // calculate the normal of the triangle
            Vector n2 = unit( cross( v5.pos - v4.pos, v6.pos - v4.pos ) );

            // TODO: maybe not needed
            if ( n.x > 0.0 && n2.x < 0.0 ) {
                n2.x *= -1;
                n2.y *= -1;
                n2.z *= -1;
            }

            // skip if normals are not the same
            if ( n != n2 ) {
                continue;
            }

            // std::cout << "dist " << distanceToLine(v4.pos, v1.pos, v2.pos) << std::endl;
            //  check if v4 lies on v1v2
            auto t = interpolate_between( v1, v2, v4 );
            std::cout << distanceToLine( v4.pos, v1.pos, v2.pos ) << std::endl;
            if ( approxequal( distanceToLine( v4.pos, v1.pos, v2.pos ), 0.0 ) && ( t > 1e-6 && t < ( 1 - 1e-6 ) ) ) {
                // Replace v1v2v3 with two triangles
                // toDelete.push_back(i);
                toDelete.push_back( i );

                // Create two new triangles v1v4v3 and v3v4v2
                mesh.indices.insert( mesh.indices.end(), { mesh.indices[i], mesh.indices[j], mesh.indices[i + 2] } );
                mesh.indices.insert( mesh.indices.end(), { mesh.indices[i + 2], mesh.indices[j], mesh.indices[i + 1] } );

                std::cout << "t=" << t << std::endl;
                std::cout << "added new triangles" << std::endl;
            }

            // check if v5 lies on v1v2
            t = interpolate_between( v1, v2, v5 );
            if ( approxequal( distanceToLine( v5.pos, v1.pos, v2.pos ), 0.0 ) && ( t > 1e-6 && t < ( 1 - 1e-6 ) ) ) {
                // Replace v1v2v3 with two triangles
                toDelete.push_back( i );

                // Create two new triangles v1v5v3 and v3v5v2
                mesh.indices.insert( mesh.indices.end(), { mesh.indices[i], mesh.indices[j + 1], mesh.indices[i + 2] } );
                mesh.indices.insert( mesh.indices.end(), { mesh.indices[i + 2], mesh.indices[j + 1], mesh.indices[i + 1] } );

                std::cout << "added new triangles" << std::endl;
            }

            // check if v6 lies on v1v2
            t = interpolate_between( v1, v2, v6 );
            if ( approxequal( distanceToLine( v6.pos, v1.pos, v2.pos ), 0.0 ) && ( t > 1e-6 && t < ( 1 - 1e-6 ) ) ) {
                // Replace v1v2v3 with two triangles
                toDelete.push_back( i );

                // Create two new triangles v1v6v3 and v3v6v2
                mesh.indices.insert( mesh.indices.end(), { mesh.indices[i], mesh.indices[j + 2], mesh.indices[i + 2] } );
                mesh.indices.insert( mesh.indices.end(), { mesh.indices[i + 2], mesh.indices[j + 2], mesh.indices[i + 1] } );

                std::cout << "added new triangles" << std::endl;
            }

            // check if v4 lies on v2v3
            t = interpolate_between( v2, v3, v4 );
            if ( approxequal( distanceToLine( v4.pos, v2.pos, v3.pos ), 0.0 ) && ( t > 1e-6 && t < ( 1 - 1e-6 ) ) ) {
                // Replace v1v2v3 with two triangles
                toDelete.push_back( i );

                // Create two new triangles v1v2v4 and v1v4v3
                mesh.indices.insert( mesh.indices.end(), { mesh.indices[i], mesh.indices[i + 1], mesh.indices[j] } );
                mesh.indices.insert( mesh.indices.end(), { mesh.indices[i], mesh.indices[j], mesh.indices[i + 2] } );

                std::cout << "added new triangles" << std::endl;
            }

            // check if v5 lies on v2v3
            t = interpolate_between( v2, v3, v5 );
            if ( approxequal( distanceToLine( v5.pos, v2.pos, v3.pos ), 0.0 ) && ( t > 1e-6 && t < ( 1 - 1e-6 ) ) ) {
                // Replace v1v2v3 with two triangles
                toDelete.push_back( i );

                // Create two new triangles v1v2v5 and v1v5v3
                mesh.indices.insert( mesh.indices.end(), { mesh.indices[i], mesh.indices[i + 1], mesh.indices[j + 1] } );
                mesh.indices.insert( mesh.indices.end(), { mesh.indices[i], mesh.indices[j + 1], mesh.indices[i + 2] } );

                std::cout << "added new triangles" << std::endl;
            }

            // check if v6 lies on v2v3
            t = interpolate_between( v2, v3, v6 );
            if ( approxequal( distanceToLine( v6.pos, v2.pos, v3.pos ), 0.0 ) && ( t > 1e-6 && t < ( 1 - 1e-6 ) ) ) {
                // Replace v1v2v3 with two triangles
                toDelete.push_back( i );

                // Create two new triangles v1v2v6 and v1v6v3
                mesh.indices.insert( mesh.indices.end(), { mesh.indices[i], mesh.indices[i + 1], mesh.indices[j + 2] } );
                mesh.indices.insert( mesh.indices.end(), { mesh.indices[i], mesh.indices[j + 2], mesh.indices[i + 2] } );

                std::cout << "added new triangles" << std::endl;
            }

            // check if v4 lies on v3v1
            t = interpolate_between( v3, v1, v4 );
            if ( approxequal( distanceToLine( v4.pos, v3.pos, v1.pos ), 0.0 ) && ( t > 1e-6 && t < ( 1 - 1e-6 ) ) ) {
                // Replace v1v2v3 with two triangles
                toDelete.push_back( i );

                // Create two new triangles v1v2v4 and v3v4v2
                mesh.indices.insert( mesh.indices.end(), { mesh.indices[i], mesh.indices[i + 1], mesh.indices[j] } );
                mesh.indices.insert( mesh.indices.end(), { mesh.indices[i + 2], mesh.indices[j], mesh.indices[i + 1] } );

                std::cout << "added new triangles" << std::endl;
            }

            // check if v5 lies on v2v1
            t = interpolate_between( v3, v1, v5 );
            if ( approxequal( distanceToLine( v5.pos, v3.pos, v1.pos ), 0.0 ) && ( t > 1e-6 && t < ( 1 - 1e-6 ) ) ) {
                // Replace v1v2v3 with two triangles
                toDelete.push_back( i );

                // Create two new triangles v1v2v5 and v3v5v2
                mesh.indices.insert( mesh.indices.end(), { mesh.indices[i], mesh.indices[i + 1], mesh.indices[j + 1] } );
                mesh.indices.insert( mesh.indices.end(), { mesh.indices[i + 2], mesh.indices[j + 1], mesh.indices[i + 1] } );

                std::cout << "added new triangles" << std::endl;
            }

            // check if v6 lies on v3v1
            t = interpolate_between( v3, v1, v6 );
            if ( approxequal( distanceToLine( v6.pos, v3.pos, v1.pos ), 0.0 ) && ( t > 1e-6 && t < ( 1 - 1e-6 ) ) ) {
                // Replace v1v2v3 with two triangles
                toDelete.push_back( i );

                // Create two new triangles v1v2v6 and v3v6v2
                mesh.indices.insert( mesh.indices.end(), { mesh.indices[i], mesh.indices[i + 1], mesh.indices[j + 2] } );
                mesh.indices.insert( mesh.indices.end(), { mesh.indices[i + 2], mesh.indices[j + 2], mesh.indices[i + 1] } );

                std::cout << "added new triangles" << std::endl;
            }
        }
    }
}

void writeSTL( const Model &model, std::string filename )
{
    std::ofstream file( filename );
    file << "solid" << std::endl;

    for ( size_t i = 0; i < model.indices.size(); i += 3 ) {
        const Vertex &v1 = model.vertices[model.indices[i]];
        const Vertex &v2 = model.vertices[model.indices[i + 1]];
        const Vertex &v3 = model.vertices[model.indices[i + 2]];
        auto n           = cross( v2.pos - v1.pos, v3.pos - v1.pos );
        n                = unit( n );
        file << "facet normal " << n.x << " " << n.y << " " << n.z << std::endl;
        file << "outer loop" << std::endl;
        file << "vertex " << v1.pos.x << " " << v1.pos.y << " " << v1.pos.z << std::endl;
        file << "vertex " << v2.pos.x << " " << v2.pos.y << " " << v2.pos.z << std::endl;
        file << "vertex " << v3.pos.x << " " << v3.pos.y << " " << v3.pos.z << std::endl;
        file << "endloop" << std::endl;
        file << "endfacet" << std::endl;
    }

    file << "endsolid" << std::endl;
}

void writeModelsSTL( const std::vector<Model> &models, std::string filename )
{
    std::ofstream file( filename );
    file << "solid" << std::endl;

    for ( auto &model : models ) {
        for ( size_t i = 0; i < model.indices.size(); i += 3 ) {
            const Vertex &v1 = model.vertices[model.indices[i]];
            const Vertex &v2 = model.vertices[model.indices[i + 1]];
            const Vertex &v3 = model.vertices[model.indices[i + 2]];
            auto n           = cross( v2.pos - v1.pos, v3.pos - v1.pos );
            n                = unit( n );
            file << "facet normal " << n.x << " " << n.y << " " << n.z << std::endl;
            file << "outer loop" << std::endl;
            file << "vertex " << v1.pos.x << " " << v1.pos.y << " " << v1.pos.z << std::endl;
            file << "vertex " << v2.pos.x << " " << v2.pos.y << " " << v2.pos.z << std::endl;
            file << "vertex " << v3.pos.x << " " << v3.pos.y << " " << v3.pos.z << std::endl;
            file << "endloop" << std::endl;
            file << "endfacet" << std::endl;
        }
    }

    file << "endsolid" << std::endl;
}

// public interface - not super efficient, if you use multiple CSG operations you should
// use BSP trees and convert them into model only once. Another optimization trick is
// replacing model with your own class.

Model csgunion( const Model &a, const Model &b );
Model csgintersection( const Model &a, const Model &b );
Model csgsubtract( const Model &a, const Model &b );

std::vector<Polygon> csgunion( const std::vector<Polygon> &a, const std::vector<Polygon> &b );
std::vector<Polygon> csgintersection( const std::vector<Polygon> &a, const std::vector<Polygon> &b );
std::vector<Polygon> csgsubtract( const std::vector<Polygon> &a, const std::vector<Polygon> &b );

/* API to build a set of polygons representning primatves. */
std::vector<Polygon> csgpolygon_cube( const Vector &center = { 0.0, 0.0, 0.0 },
    const Vector &dim = { 1.0, 1.0, 1.0 }, const int col = 0xFFFFFF );
std::vector<Polygon> csgpolygon_sphere( const Vector &center = { 0.0, 0.0, 0.0 }, double radius = 1.0,
    const int col = 0xFFFFFF, int slices = 16, int stacks = 8 );
std::vector<Polygon> csgpolygon_cylinder( const Vector &s = { 0.0, -1.0, 0.0 },
    const Vector &e = { 0.0, 1.0, 0.0 }, double radius = 1.0,
    const int col = 0xFFFFFF, int slices = 16 );

std::vector<Polygon> csgfixtjunc( const std::vector<Polygon> &polygons );

Model modelfrompolygons( const std::vector<Polygon> &polygons );

/* API to build models representing primatives */
Model csgmodel_cube( const Vector &center = { 0.0, 0.0, 0.0 }, const Vector &dim = { 1.0, 1.0, 1.0 },
    const int col = 0xFFFFFF );
Model csgmodel_sphere( const Vector &center = { 0.0, 0.0, 0.0 }, double radius = 1.0,
    const int col = 0xFFFFFF, int slices = 16, int stacks = 8 );

Model csgmodel_cylinder( const Vector &s = { 0.0, -1.0, 0.0 }, const Vector &e = { 0.0, 1.0, 0.0 },
    double radius = 1.0, const int col = 0xFFFFFF, int slices = 16 );

#include <assert.h>

// Holds a node in a BSP tree. A BSP tree is built from a collection of polygons
// by picking a polygon to split along. That polygon (and all other coplanar
// polygons) are added directly to that node and the other polygons are added to
// the front and/or back subtrees. This is not a leafy BSP tree since there is
// no distinction between internal and leaf nodes.
struct CSGNode {
    std::vector<Polygon> polygons;
    CSGNode *front;
    CSGNode *back;
    Plane plane;

    CSGNode();
    CSGNode( const std::vector<Polygon> &list );
    ~CSGNode();

    CSGNode *clone() const;
    void clipto( const CSGNode *other );
    void invert();
    void build( const std::vector<Polygon> &Polygon );
    std::vector<Polygon> clippolygons( const std::vector<Polygon> &list ) const;
    std::vector<Polygon> allpolygons() const;
};

// Vertex implementation

// Invert all orientation-specific data (e.g. Vertex normal). Called when the
// orientation of a polygon is flipped.
inline Vertex flip( Vertex v )
{
    v.normal = negate( v.normal );
    return v;
}

// Create a new Vertex between this Vertex and `other` by linearly
// interpolating all properties using a parameter of `t`. Subclasses should
// override this to interpolate additional properties.
inline Vertex interpolate( const Vertex &a, const Vertex &b, double t )
{
    Vertex ret;
    ret.pos    = lerp( a.pos, b.pos, t );
    ret.normal = lerp( a.normal, b.normal, t );
    ret.col    = lerp( a.col, b.col, t );
    return ret;
}

// Plane implementation

Plane::Plane() :
    normal(), w( 0.0 )
{
}

Plane::Plane( const Vector &a, const Vector &b, const Vector &c )
{
    this->normal = unit( cross( b - a, c - a ) );
    this->w      = dot( this->normal, a );
}

// Split `polygon` by this plane if needed, then put the polygon or polygon
// fragments in the appropriate lists. Coplanar polygons go into either
// `coplanarFront` or `coplanarBack` depending on their orientation with
// respect to this plane. Polygons in front or in back of this plane go into
// either `front` or `back`.
void Plane::splitpolygon( const Polygon &poly, std::vector<Polygon> &coplanarFront,
    std::vector<Polygon> &coplanarBack, std::vector<Polygon> &front,
    std::vector<Polygon> &back ) const
{

    // Classify each point as well as the entire polygon into one of the above
    // four classes.
    int polygonType = 0;
    for ( const auto &v : poly.vertices ) {
        polygonType |= classify( v.pos );
    }

    // Put the polygon in the correct list, splitting it when necessary.
    switch ( polygonType ) {
    case COPLANAR: {
        if ( dot( this->normal, poly.plane.normal ) > 0 )
            coplanarFront.push_back( poly );
        else
            coplanarBack.push_back( poly );
        break;
    }
    case FRONT: {
        front.push_back( poly );
        break;
    }
    case BACK: {
        back.push_back( poly );
        break;
    }
    case SPANNING: {
        std::vector<Vertex> f, b;

        for ( size_t i = 0; i < poly.vertices.size(); i++ ) {

            size_t j = ( i + 1 ) % poly.vertices.size();

            const Vertex &vi = poly.vertices[i];
            const Vertex &vj = poly.vertices[j];

            int ti = classify( vi.pos );
            int tj = classify( vj.pos );

            if ( ti != BACK )
                f.push_back( vi );
            if ( ti != FRONT )
                b.push_back( vi );
            if ( ( ti | tj ) == SPANNING ) {
                double t = ( this->w - dot( this->normal, vi.pos ) ) / dot( this->normal, vj.pos - vi.pos );
                Vertex v = interpolate( vi, vj, t );
                f.push_back( v );
                b.push_back( v );
            }
        }
        if ( f.size() >= 3 )
            front.push_back( Polygon( f ) );
        if ( b.size() >= 3 )
            back.push_back( Polygon( b ) );
        break;
    }
    }
}

// Polygon implementation

Polygon::Polygon()
{
}

Polygon::Polygon( const std::vector<Vertex> &list ) :
    vertices( list ), plane( vertices[0].pos, vertices[1].pos, vertices[2].pos )
{
}

// Node implementation

// Return a new CSG solid representing space in either this solid or in the
// solid `csg`. Neither this solid nor the solid `csg` are modified.
inline CSGNode *csg_union( const CSGNode *a1, const CSGNode *b1 )
{
    CSGNode *a = a1->clone();
    CSGNode *b = b1->clone();
    a->clipto( b );
    b->clipto( a );
    b->invert();
    b->clipto( a );
    b->invert();
    a->build( b->allpolygons() );
    CSGNode *ret = new CSGNode( a->allpolygons() );
    delete a;
    delete b;
    return ret;
}

// Return a new CSG solid representing space in this solid but not in the
// solid `csg`. Neither this solid nor the solid `csg` are modified.
inline CSGNode *csg_subtract( const CSGNode *a1, const CSGNode *b1 )
{
    CSGNode *a = a1->clone();
    CSGNode *b = b1->clone();
    a->invert();
    a->clipto( b );
    b->clipto( a );
    b->invert();
    b->clipto( a );
    b->invert();
    a->build( b->allpolygons() );
    a->invert();
    CSGNode *ret = new CSGNode( a->allpolygons() );
    delete a;
    delete b;
    return ret;
}

// Return a new CSG solid representing space both this solid and in the
// solid `csg`. Neither this solid nor the solid `csg` are modified.
inline CSGNode *csg_intersect( const CSGNode *a1, const CSGNode *b1 )
{
    CSGNode *a = a1->clone();
    CSGNode *b = b1->clone();
    a->invert();
    b->clipto( a );
    b->invert();
    a->clipto( b );
    b->clipto( a );
    a->build( b->allpolygons() );
    a->invert();
    CSGNode *ret = new CSGNode( a->allpolygons() );
    delete a;
    delete b;
    return ret;
}

// Convert solid space to empty space and empty space to solid space.
void CSGNode::invert()
{
    std::deque<CSGNode *> nodes;
    nodes.push_back( this );
    while ( nodes.size() ) {
        CSGNode *me = nodes.front();
        nodes.pop_front();

        for ( size_t i = 0; i < me->polygons.size(); i++ )
            me->polygons[i].flip();
        me->plane.flip();
        std::swap( me->front, me->back );
        if ( me->front )
            nodes.push_back( me->front );
        if ( me->back )
            nodes.push_back( me->back );
    }
}

// Recursively remove all polygons in `polygons` that are inside this BSP
// tree.
std::vector<Polygon> CSGNode::clippolygons( const std::vector<Polygon> &ilist ) const
{
    std::vector<Polygon> result;

    std::deque<std::pair<const CSGNode *const, std::vector<Polygon> > > clips;
    clips.push_back( std::make_pair( this, ilist ) );
    while ( clips.size() ) {
        const CSGNode *me                = clips.front().first;
        const std::vector<Polygon> &list = clips.front().second;

        if ( !me->plane.ok() ) {
            result.insert( result.end(), list.begin(), list.end() );
            clips.pop_front();
            continue;
        }

        std::vector<Polygon> list_front, list_back;
        for ( size_t i = 0; i < list.size(); i++ )
            me->plane.splitpolygon( list[i], list_front, list_back, list_front, list_back );

        if ( me->front )
            clips.push_back( std::make_pair( me->front, list_front ) );
        else
            result.insert( result.end(), list_front.begin(), list_front.end() );

        if ( me->back )
            clips.push_back( std::make_pair( me->back, list_back ) );

        clips.pop_front();
    }

    return result;
}

// Remove all polygons in this BSP tree that are inside the other BSP tree
// `bsp`.
void CSGNode::clipto( const CSGNode *other )
{
    std::deque<CSGNode *> nodes;
    nodes.push_back( this );
    while ( nodes.size() ) {
        CSGNode *me = nodes.front();
        nodes.pop_front();

        me->polygons = other->clippolygons( me->polygons );
        if ( me->front )
            nodes.push_back( me->front );
        if ( me->back )
            nodes.push_back( me->back );
    }
}

// Return a list of all polygons in this BSP tree.
std::vector<Polygon> CSGNode::allpolygons() const
{
    std::vector<Polygon> result;

    std::deque<const CSGNode *> nodes;
    nodes.push_back( this );
    while ( nodes.size() ) {
        const CSGNode *me = nodes.front();
        nodes.pop_front();

        result.insert( result.end(), me->polygons.begin(), me->polygons.end() );
        if ( me->front )
            nodes.push_back( me->front );
        if ( me->back )
            nodes.push_back( me->back );
    }

    return result;
}

CSGNode *CSGNode::clone() const
{
    CSGNode *ret = new CSGNode();

    std::deque<std::pair<const CSGNode *, CSGNode *> > nodes;
    nodes.push_back( std::make_pair( this, ret ) );
    while ( nodes.size() ) {
        const CSGNode *original = nodes.front().first;
        CSGNode *clone          = nodes.front().second;
        nodes.pop_front();

        clone->polygons = original->polygons;
        clone->plane    = original->plane;
        if ( original->front ) {
            clone->front = new CSGNode();
            nodes.push_back( std::make_pair( original->front, clone->front ) );
        }
        if ( original->back ) {
            clone->back = new CSGNode();
            nodes.push_back( std::make_pair( original->back, clone->back ) );
        }
    }

    return ret;
}

// Build a BSP tree out of `polygons`. When called on an existing tree, the
// new polygons are filtered down to the bottom of the tree and become new
// nodes there. Each set of polygons is partitioned using the first polygon
// (no heuristic is used to pick a good split).
void CSGNode::build( const std::vector<Polygon> &ilist )
{
    if ( !ilist.size() )
        return;

    std::deque<std::pair<CSGNode *, std::vector<Polygon> > > builds;
    builds.push_back( std::make_pair( this, ilist ) );

    while ( builds.size() ) {
        CSGNode *me                      = builds.front().first;
        const std::vector<Polygon> &list = builds.front().second;

        assert( list.size() > 0 && "logic error" );

        if ( !me->plane.ok() )
            me->plane = list[0].plane;
        std::vector<Polygon> list_front, list_back;

        // me->polygons.push_back(list[0]);
        for ( size_t i = 0; i < list.size(); i++ )
            me->plane.splitpolygon( list[i], me->polygons, me->polygons, list_front, list_back );

        if ( list_front.size() ) {
            if ( !me->front )
                me->front = new CSGNode;
            builds.push_back( std::make_pair( me->front, list_front ) );
        }
        if ( list_back.size() ) {
            if ( !me->back )
                me->back = new CSGNode;
            builds.push_back( std::make_pair( me->back, list_back ) );
        }

        builds.pop_front();
    }
}

CSGNode::CSGNode() :
    front( nullptr ), back( nullptr )
{
}

CSGNode::CSGNode( const std::vector<Polygon> &list ) :
    front( nullptr ), back( nullptr )
{
    build( list );
}

CSGNode::~CSGNode()
{

    std::deque<CSGNode *> nodes_to_delete;
    std::deque<CSGNode *> nodes_to_disassemble;

    nodes_to_disassemble.push_back( this );
    while ( nodes_to_disassemble.size() ) {
        CSGNode *me = nodes_to_disassemble.front();
        nodes_to_disassemble.pop_front();

        if ( me->front ) {
            nodes_to_disassemble.push_back( me->front );
            nodes_to_delete.push_back( me->front );
            me->front = NULL;
        }
        if ( me->back ) {
            nodes_to_disassemble.push_back( me->back );
            nodes_to_delete.push_back( me->back );
            me->back = NULL;
        }
    }

    for ( auto it = nodes_to_delete.begin(); it != nodes_to_delete.end(); ++it )
        delete *it;
}

// Public interface implementation

inline std::vector<Polygon> modeltopolygons( const Model &model )
{
    std::vector<Polygon> list;
    for ( size_t i = 0; i < model.indices.size(); i += 3 ) {
        std::vector<Vertex> triangle;
        for ( int j = 0; j < 3; j++ ) {
            Vertex v = model.vertices[model.indices[i + j]];
            triangle.push_back( v );
        }
        list.push_back( Polygon( triangle ) );
    }
    return list;
}

Model modelfrompolygons( const std::vector<Polygon> &polygons )
{
    Model model;

    for ( size_t i = 0; i < polygons.size(); i++ ) {
        const Polygon &poly = polygons[i];

        if ( poly.vertices.size() ) {

            Model::Index a = model.AddVertex( poly.vertices[0] );

            for ( size_t j = 2; j < poly.vertices.size(); j++ ) {

                Model::Index b = model.AddVertex( poly.vertices[j - 1] );
                Model::Index c = model.AddVertex( poly.vertices[j] );

                if ( a != b && b != c && c != a ) {
                    model.indices.push_back( a );
                    model.indices.push_back( b );
                    model.indices.push_back( c );
                }
            }
        }
    }
    return model;
}

typedef CSGNode *csg_function( const CSGNode *a1, const CSGNode *b1 );

std::vector<Polygon> csgjs_operation( const std::vector<Polygon> &apoly, const std::vector<Polygon> &bpoly,
    csg_function fun )
{

    CSGNode A( apoly );
    CSGNode B( bpoly );

    /* create a unique pointer here so we can delete AB on exit */
    std::unique_ptr<CSGNode> AB( fun( &A, &B ) );
    return AB->allpolygons();
}

inline std::vector<Polygon> csgjs_operation( const Model &a, const Model &b, csg_function fun )
{
    return csgjs_operation( modeltopolygons( a ), modeltopolygons( b ), fun );
}

std::vector<Polygon> csgpolygon_cube( const Vector &center, const Vector &dim, const int col )
{
    struct Quad {
        int indices[4];
        Vector normal;
    } quads[] = { { { 0, 4, 6, 2 }, { -1, 0, 0 } }, { { 1, 3, 7, 5 }, { +1, 0, 0 } }, { { 0, 1, 5, 4 }, { 0, -1, 0 } }, { { 2, 6, 7, 3 }, { 0, +1, 0 } }, { { 0, 2, 3, 1 }, { 0, 0, -1 } }, { { 4, 5, 7, 6 }, { 0, 0, +1 } } };

    std::vector<Polygon> polygons;
    for ( const auto &q : quads ) {

        std::vector<Vertex> verts;

        for ( auto i : q.indices ) {
            Vector pos( center.x + dim.x * ( 2.0 * !!( i & 1 ) - 1 ), center.y + dim.y * ( 2.0 * !!( i & 2 ) - 1 ),
                center.z + dim.z * ( 2.0 * !!( i & 4 ) - 1 ) );

            verts.push_back( { pos, q.normal, col } );
        }
        polygons.push_back( Polygon( verts ) );
    }
    return polygons;
}

Model csgmodel_cube( const Vector &center, const Vector &dim, int col )
{

    return modelfrompolygons( csgpolygon_cube( center, dim, col ) );
}

std::vector<Polygon> csgpolygon_sphere( const Vector &c, double r, int col, int slices, int stacks )
{
    std::vector<Polygon> polygons;

    auto mkvertex = [c, r, col]( double theta, double phi ) -> Vertex {
        theta *= (double)M_PI * 2;
        phi *= (double)M_PI;
        Vector dir( (double)cos( theta ) * (double)sin( phi ), (double)cos( phi ),
            (double)sin( theta ) * (double)sin( phi ) );

        return Vertex{ c + ( dir * r ), dir, col };
    };
    for ( double i = 0; i < slices; i++ ) {
        for ( double j = 0; j < stacks; j++ ) {

            std::vector<Vertex> vertices;

            vertices.push_back( mkvertex( i / slices, j / stacks ) );
            if ( j > 0 ) {
                vertices.push_back( mkvertex( ( i + 1 ) / slices, j / stacks ) );
            }
            if ( j < stacks - 1 ) {
                vertices.push_back( mkvertex( ( i + 1 ) / slices, ( j + 1 ) / stacks ) );
            }
            vertices.push_back( mkvertex( i / slices, ( j + 1 ) / stacks ) );
            polygons.push_back( Polygon( vertices ) );
        }
    }
    return polygons;
}

Model csgmodel_sphere( const Vector &c, double r, int col, int slices, int stacks )
{

    return modelfrompolygons( csgpolygon_sphere( c, r, col, slices, stacks ) );
}

std::vector<Polygon> csgpolygon_cylinder( const Vector &s, const Vector &e, double r, int col,
    int slices )
{
    Vector ray = e - s;

    Vector axisZ = unit( ray );
    bool isY     = fabs( axisZ.y ) > 0.5f;
    Vector axisX = unit( cross( Vector( isY, !isY, 0 ), axisZ ) );
    Vector axisY = unit( cross( axisX, axisZ ) );

    Vertex start{ s, -axisZ, col };
    Vertex end{ e, unit( axisZ ), col };

    std::vector<Polygon> polygons;

    auto point = [axisX, axisY, s, r, ray, axisZ, col]( double stack, double slice,
                     double normalBlend ) -> Vertex {
        double angle  = slice * (double)M_PI * 2;
        Vector out    = axisX * (double)cos( angle ) + axisY * (double)sin( angle );
        Vector pos    = s + ray * stack + out * r;
        Vector normal = out * ( 1.0 - fabs( normalBlend ) ) + axisZ * normalBlend;
        return Vertex{ pos, normal, col };
    };

    for ( double i = 0; i < slices; i++ ) {
        double t0 = i / slices;
        double t1 = ( i + 1 ) / slices;
        polygons.push_back( Polygon( { start, point( 0, t0, -1 ), point( 0, t1, -1 ) } ) );
        polygons.push_back( Polygon( { point( 0, t1, 0 ), point( 0, t0, 0 ), point( 1, t0, 0 ), point( 1, t1, 0 ) } ) );
        polygons.push_back( Polygon( { end, point( 1, t1, 1 ), point( 1, t0, 1 ) } ) );
    }
    return polygons;
}

Model csgmodel_cylinder( const Vector &s, const Vector &e, double r, int col, int slices )
{

    return modelfrompolygons( csgpolygon_cylinder( s, e, r, col, slices ) );
}

Model csgunion( const Model &a, const Model &b )
{
    return modelfrompolygons( csgjs_operation( a, b, csg_union ) );
}

Model csgintersection( const Model &a, const Model &b )
{
    return modelfrompolygons( csgjs_operation( a, b, csg_intersect ) );
}

Model csgsubtract( const Model &a, const Model &b )
{
    return modelfrompolygons( csgjs_operation( a, b, csg_subtract ) );
}

std::vector<Polygon> csgunion( const std::vector<Polygon> &a, const std::vector<Polygon> &b )
{
    return csgjs_operation( a, b, csg_union );
}

std::vector<Polygon> csgintersection( const std::vector<Polygon> &a, const std::vector<Polygon> &b )
{
    return csgjs_operation( a, b, csg_intersect );
}

std::vector<Polygon> csgsubtract( const std::vector<Polygon> &a, const std::vector<Polygon> &b )
{
    return csgjs_operation( a, b, csg_subtract );
}

template <class CONTTYPE, class T>
bool contains( CONTTYPE &cont, const T &t ) { return cont.find( t ) != cont.end(); }

template <class CONTTYPE, class T>
typename CONTTYPE::iterator find( CONTTYPE &cont, const T &t );

template <class T, typename... Ks>
typename std::map<Ks...>::iterator find( std::map<Ks...> &cont, const T &t ) { return cont.find( t ); }

template <typename T>
typename std::vector<T>::iterator find( std::vector<T> &cont, const T &t ) { return std::find( cont.begin(), cont.end(), t ); }

// template <typename CONTTYPE, typename T>
// void remove(CONTTYPE &cont, const T &t);

template <typename T>
void remove( std::vector<T> &cont, const T &t )
{
    auto i = find( cont, t );
    if ( i != cont.end() ) {
        cont.erase( i );
    }
}

template <typename T, typename... Ks>
void remove( std::map<Ks...> &cont, const T &t )
{
    auto i = cont.find( t );
    if ( i != cont.end() ) {
        cont.erase( i );
    }
}

const int kInvalidPolygonIndex = -1;

using Tag = size_t;
// Tag GetTag (const Vertex &v) {
//	return (Tag)&v;
// };

using SideTag = std::pair<Tag, Tag>;
bool operator==( const SideTag &a, const SideTag &b )
{
    return a.first == b.first && a.second == b.second;
}

// transcibed from
// https://github.com/jscad/csg.js/blob/6be72558e47355d59091d5684f3c4ed853476404/csg.js#L1090
/*
fixTJunctions:
Suppose we have two polygons ACDB and EDGF:
    A-----B
    |     |
    |     E--F
    |     |  |
    C-----D--G
Note that vertex E forms a T-junction on the side BD. In this case some STL slicers will complain
that the solid is not watertight. This is because the watertightness check is done by checking if
each side DE is matched by another side ED.
This function will return a new solid with ACDB replaced by ACDEB
Note that this can create polygons that are slightly non-convex (due to rounding errors). Therefore the result should
not be used for further CSG operations!
*/
std::vector<Polygon> csgfixtjunc( const std::vector<Polygon> &originalpolygons )
{

    // were going to need unique vertices so while we're at it
    // create a list of unique vertices and a list of polygons
    // with indexes into this list, we can use those indexes as
    // unique vertex id's in the core of the algo.
    struct IndexedVertex {
        Vertex vertex;
        size_t index; // index in to the vextor this vertex is in, also used as the unique tag.
    };

    struct Side {
        const IndexedVertex *vertex0;
        const IndexedVertex *vertex1;
        int polygonindex;
    };

    std::vector<IndexedVertex> uvertices;

    struct IndexPolygon {
        std::vector<size_t> vertexindex;
    };
    std::vector<IndexPolygon> polygons;
    for ( const auto &originalpoly : originalpolygons ) {
        IndexPolygon ipoly;
        for ( const auto &v : originalpoly.vertices ) {
            auto pos = std::find_if( uvertices.begin(), uvertices.end(), [v]( const IndexedVertex &a ) { return a.vertex == v; } );
            if ( pos != uvertices.end() ) {
                size_t index = pos - uvertices.begin();
                ipoly.vertexindex.push_back( index );
            } else {
                size_t index = uvertices.size();
                ipoly.vertexindex.push_back( index );
                uvertices.push_back( { v, index } );
            }
        }
        polygons.push_back( ipoly );
    }

    /* side map contains all sides that don't have a matching opposite
     * side AB is removed of a side BA is in the map for example.
     */
    std::map<SideTag, std::vector<Side> > sidemap = {};

    for ( int polygonindex = 0; polygonindex < (int)polygons.size(); polygonindex++ ) {
        auto &polygon      = polygons[polygonindex];
        size_t numvertices = polygon.vertexindex.size();
        if ( numvertices < 3 ) { // should be true
            continue;
        }

        IndexedVertex *vertex = &uvertices[polygon.vertexindex[0]];
        Tag vertextag         = vertex->index;
        for ( size_t vertexindex = 0; vertexindex < numvertices; vertexindex++ ) {
            size_t nextvertexindex = vertexindex + 1;
            if ( nextvertexindex == numvertices ) {
                nextvertexindex = 0;
            }

            IndexedVertex *nextvertex = &uvertices[polygon.vertexindex[nextvertexindex]];
            Tag nextvertextag         = nextvertex->index;
            auto sidetag              = std::make_pair( vertextag, nextvertextag );
            auto reversesidetag       = std::make_pair( nextvertextag, vertextag );

            auto sidemappos = sidemap.find( reversesidetag );
            if ( sidemappos != sidemap.end() ) {
                // this side matches the same side in another polygon. Remove from sidemap:
                // NOTE: dazza hmm not sre about just popping back here but the JS does splice(-1, 1) on it's array
                auto &ar = sidemappos->second;
                ar.pop_back();
                if ( ar.size() == 0 ) {
                    sidemap.erase( sidemappos );
                }
            } else {
                sidemap[sidetag].push_back( { vertex, nextvertex, polygonindex } );

                // var sideobj = {
                //     vertex0: vertex,
                //     vertex1 : nextvertex,
                //     polygonindex : polygonindex
                // };
                // if (!(sidetag in sidemap)) {
                //     sidemap[sidetag] = [sideobj];
                // } else {
                //     sidemap[sidetag].push(sideobj);
                // }
            }
            vertex    = nextvertex;
            vertextag = nextvertextag;
        }
    }

    // now sidemap contains 'unmatched' sides
    // i.e. side AB in one polygon does not have a matching side BA in another polygon

    // all sides that have vertex tag as it's start
    std::map<Tag, std::vector<SideTag> > vertextag2sidestart;
    // all sides that have vertex tag as it's end.
    std::map<Tag, std::vector<SideTag> > vertextag2sideend;
    std::deque<SideTag> sidestocheck;
    ;
    bool sidemapisempty = true;
    for ( const auto &iter : sidemap ) {
        const SideTag &sidetag            = iter.first;
        const std::vector<Side> &sideobjs = iter.second;

        sidemapisempty = false;
        sidestocheck.push_back( sidetag );

        for ( auto &sideobj : sideobjs ) {
            Tag starttag = sideobj.vertex0->index;
            Tag endtag   = sideobj.vertex1->index;
            vertextag2sidestart[starttag].push_back( sidetag );
            // if (starttag in vertextag2sidestart) {
            //     vertextag2sidestart[starttag].push(sidetag);
            // } else {
            //     vertextag2sidestart[starttag] = [sidetag];
            // }
            vertextag2sideend[endtag].push_back( sidetag );
            // if (endtag in vertextag2sideend) {
            //     vertextag2sideend[endtag].push(sidetag);
            // } else {
            //     vertextag2sideend[endtag] = [sidetag];
            // }
        }
    }

    // make a copy of the polygons array, since we are going to modify it:
    // auto polygons = inpolygons;
    if ( !sidemapisempty ) {

        auto deleteSide = [&sidemap, &vertextag2sidestart, &vertextag2sideend]( const IndexedVertex &vertex0, const IndexedVertex &vertex1, int polygonindex ) {
            auto starttag = vertex0.index;
            auto endtag   = vertex1.index;
            auto sidetag  = std::make_pair( starttag, endtag );
            // console.log("deleteSide("+sidetag+")");
            assert( contains( sidemap, sidetag ) && "logic error" );

            // todo this is better done with an iterator probably.
            int idx       = -1;
            auto sideobjs = sidemap.at( sidetag );
            for ( size_t i = 0; i < sideobjs.size(); i++ ) {
                auto &sideobj = sideobjs[i];
                if ( sideobj.vertex0->index != vertex0.index )
                    continue;
                if ( sideobj.vertex1->index != vertex1.index )
                    continue;
                if ( polygonindex != kInvalidPolygonIndex ) {
                    if ( sideobj.polygonindex != polygonindex )
                        continue;
                }
                idx = (int)i;
                break;
            }
            assert( idx >= 0 && "logic error" );
            sideobjs.erase( sideobjs.begin() + idx );
            if ( sideobjs.size() == 0 ) {
                remove( sidemap, sidetag );
            }
            auto siter = find( vertextag2sidestart[starttag], sidetag );
            assert( siter != vertextag2sidestart[starttag].end() && "logic error" );
            vertextag2sidestart[starttag].erase( siter );
            if ( vertextag2sidestart[starttag].size() == 0 ) {
                remove( vertextag2sidestart, starttag );
            }
            auto eiter = find( vertextag2sideend[endtag], sidetag );
            assert( eiter != vertextag2sideend[endtag].end() && "logic error" );
            vertextag2sideend[endtag].erase( eiter );
            if ( vertextag2sideend[endtag].size() == 0 ) {
                remove( vertextag2sideend, endtag );
            }
        };

        auto addSide = [&sidemap, &deleteSide, &vertextag2sidestart, &vertextag2sideend]( const IndexedVertex &vertex0, const IndexedVertex &vertex1, int polygonindex, SideTag &addedtag ) -> bool {
            auto starttag = vertex0.index;
            auto endtag   = vertex1.index;
            assert( starttag != endtag && "logic error" );
            auto newsidetag     = std::make_pair( starttag, endtag );
            auto reversesidetag = std::make_pair( endtag, starttag );
            if ( contains( sidemap, reversesidetag ) ) {
                // we have a matching reverse oriented side.
                // Instead of adding the new side, cancel out the reverse side:
                // console.log("addSide("+newsidetag+") has reverse side:");
                deleteSide( vertex1, vertex0, kInvalidPolygonIndex );
                return false;
            }
            //  console.log("addSide("+newsidetag+")");
            Side newsideobj{ &vertex0, &vertex1, polygonindex };
            sidemap[newsidetag].push_back( newsideobj );
            vertextag2sidestart[starttag].push_back( newsidetag );
            vertextag2sideend[endtag].push_back( newsidetag );
            addedtag = newsidetag;
            return true;
        };

        while ( true ) {
            // todo outerscope has a sidemapisempty so renamed this just incase for now
            bool sidemapisempty2 = true;
            for ( auto &iter : sidemap ) {
                const auto &sidetag = iter.first;
                sidemapisempty2     = false;
                sidestocheck.push_back( sidetag );
            }
            if ( sidemapisempty2 ) {
                break;
            }
            bool donesomething = false;
            while ( false == sidestocheck.empty() ) {

                SideTag sidetagtocheck = sidestocheck.front();
                sidestocheck.pop_front();
                // var sidetagtocheck = null;
                // for (var sidetag in sidestocheck) {
                //     sidetagtocheck = sidetag;
                //     break;
                // }
                // if (sidetagtocheck == = null) break; // sidestocheck is empty, we're done!
                bool donewithside = true;
                if ( contains( sidemap, sidetagtocheck ) ) {
                    auto &sideobjs = sidemap[sidetagtocheck];
                    assert( sideobjs.size() && "didn't expect an empty set of sides" );

                    Side &side = sideobjs[0];
                    for ( int directionindex = 0; directionindex < 2; directionindex++ ) {
                        auto startvertex    = ( directionindex == 0 ) ? side.vertex0 : side.vertex1;
                        auto endvertex      = ( directionindex == 0 ) ? side.vertex1 : side.vertex0;
                        auto startvertextag = startvertex->index;
                        auto endvertextag   = endvertex->index;
                        std::vector<SideTag> matchingsides; // TODO - this is gonna be copied into, ok with that?
                        if ( directionindex == 0 ) {
                            if ( contains( vertextag2sideend, startvertextag ) ) {
                                matchingsides = vertextag2sideend[startvertextag];
                            }
                        } else {
                            if ( contains( vertextag2sidestart, startvertextag ) ) {
                                matchingsides = vertextag2sidestart[startvertextag];
                            }
                        }
                        for ( const auto &matchingsidetag : matchingsides ) {

                            auto matchingside               = sidemap[matchingsidetag][0];
                            auto matchingsidestartvertex    = ( directionindex == 0 ) ? matchingside.vertex0 : matchingside.vertex1;
                            auto matchingsidestartvertextag = matchingsidestartvertex->index;
#if !defined( NDEBUG )
                            auto matchingsideendvertex    = ( directionindex == 0 ) ? matchingside.vertex1 : matchingside.vertex0;
                            auto matchingsideendvertextag = matchingsideendvertex->index;
                            assert( matchingsideendvertextag == startvertextag && "logic error" );
#endif
                            if ( matchingsidestartvertextag == endvertextag ) {
                                // matchingside cancels sidetagtocheck
                                deleteSide( *startvertex, *endvertex, kInvalidPolygonIndex );
                                deleteSide( *endvertex, *startvertex, kInvalidPolygonIndex );
                                donewithside   = false;
                                directionindex = 2; // skip reverse direction check
                                donesomething  = true;
                                break;
                            } else {
                                auto startpos  = startvertex->vertex.pos;
                                auto endpos    = endvertex->vertex.pos;
                                auto checkpos  = matchingsidestartvertex->vertex.pos;
                                auto direction = checkpos - startpos;
                                // Now we need to check if endpos is on the line startpos-checkpos:
                                double t = dot( ( endpos - startpos ), direction ) / dot( direction, direction );
                                if ( ( t > 0.0 ) && ( t < 1.0 ) ) {
                                    auto closestpoint    = startpos + direction * t;
                                    auto distancesquared = lengthsquared( closestpoint - endpos );
                                    if ( distancesquared < 1e-10 ) { // TODO - shouldn't this be epsilon constant?
                                        // Yes it's a t-junction! We need to split matchingside in two:
                                        auto polygonindex = matchingside.polygonindex;
                                        auto polygon      = polygons[polygonindex];
                                        // find the index of startvertextag in polygon:
                                        auto insertionvertextag     = matchingside.vertex1->index;
                                        int insertionvertextagindex = -1;
                                        for ( int i = 0; i < (int)polygon.vertexindex.size(); i++ ) {
                                            if ( polygon.vertexindex[i] == insertionvertextag ) {
                                                insertionvertextagindex = i;
                                                break;
                                            }
                                        }
                                        assert( insertionvertextagindex >= 0 && "logic error" );
                                        // split the side by inserting the vertex:
                                        auto newvertices = polygon.vertexindex; // deliberate copy TODO, is it needed, we're just inserting a vertex!
                                        newvertices.insert( newvertices.begin() + insertionvertextagindex, endvertex->index );
                                        polygons[polygonindex].vertexindex = newvertices;

                                        // remove the original sides from our maps:
                                        // deleteSide(sideobj.vertex0, sideobj.vertex1, null);
                                        deleteSide( *matchingside.vertex0, *matchingside.vertex1, polygonindex );
                                        SideTag newsidetag1, newsidetag2;
                                        if ( addSide( *matchingside.vertex0, *endvertex, polygonindex, newsidetag1 ) ) {
                                            sidestocheck.push_back( newsidetag1 );
                                        }
                                        if ( addSide( *endvertex, *matchingside.vertex1, polygonindex, newsidetag2 ) ) {
                                            sidestocheck.push_back( newsidetag2 );
                                        }
                                        donewithside   = false;
                                        directionindex = 2; // skip reverse direction check
                                        donesomething  = true;
                                        break;
                                    } // if(distancesquared < 1e-10)
                                } // if( (t > 0) && (t < 1) )
                            } // if(endingstidestartvertextag == endvertextag)
                        } // for matchingsideindex
                    } // for directionindex
                } // if(sidetagtocheck in sidemap)
                if ( donewithside ) {
                    // delete sidestocheck[sidetag];
                }
            }
            if ( !donesomething ) {
                break;
            }
        } // if(!sidemapisempty)
    }

    std::vector<Polygon> outpolys;
    for ( const auto &indexedpoly : polygons ) {
        Polygon p;
        for ( auto i : indexedpoly.vertexindex ) {
            p.vertices.push_back( uvertices[i].vertex );
        }
        assert( p.vertices.size() > 2 && "logic error" );
        p.plane = Plane( p.vertices[0].pos, p.vertices[1].pos, p.vertices[2].pos );
        outpolys.push_back( p );
    }
    return outpolys;
}
