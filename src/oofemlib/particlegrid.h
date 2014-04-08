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


#ifndef particlegrid_h
#define particlegrid_h

#include "floatarray.h"
#include "intarray.h"
#include "mathfem.h"

#include <list>

namespace oofem
{

template< class Point >
class ParticleGridIterator;

/**
 * Particle grid data structure for n-D grids.
 * Supports local refinements through recursively implementing sub grids into points.
 * @author Mikael Öhman
 */
template <class Point>
class OOFEM_EXPORT ParticleGrid
{
public:
    /// List iterator type
    typedef ParticleGridIterator<Point> iterator;

protected:
    /// Recursive data structure for
    struct RefinedParticlePoint {
        Point *point;
        ParticleGrid<Point> *subgrid;
    };

    /// Number of dimensions of the grid.
    int n;
    /// Resolution in each dimension.
    IntArray res;
    /// Total grid points.
    int total;
    /// Bounding boxes
    FloatArray bb0,bb1;
    /// Grid increments.
    FloatArray dx;
    /// Particle data structure.
    RefinedParticlePoint **data;
    /// Helper for index <-> position, the cumulative product of res.
    IntArray res_prod;

    /**
     * Translates from n-D position to total position.
     * @param pos n-D position.
     * @return Total position.
     */
    int giveIndex ( const IntArray &pos ) const;
    /**
     * Translates from total position to n-D position.
     * @param pos n-D position.
     * @param index Total position.
     */
    void givePosition ( IntArray &pos, int index ) const;

    /**
     * Initiation of structure, used by constructors.
     * @param res Resolution array.
     * @param bb0 Lower bounding box.
     * @param bb1 Upper bounding box.
     */
    void init ( const IntArray &res, const FloatArray &bb0, const FloatArray &bb1 );

public:
    /**
     * Creates a new empty particle grid
     * @param res Resolution in respective direction.
     * @param bb0 Lower corner of bounding box.
     * @param bb1 Upper corner of bounding box.
     */
    ParticleGrid ( const IntArray &res, const FloatArray &bb0, const FloatArray &bb1 );
    /**
     * Copies structure of old particle grid, but no content.
     * @param old ParticleGrid to copy.
     */
    ParticleGrid ( const ParticleGrid *old );
    /// Destructor.
    ~ParticleGrid();

    /// Returns Number of dimensions of grid.
    int giveDimensions() {
        return n;
    }
    /// Total number potential points.
    int getTotal() const {
        return this->total;
    }
    /// Number of active points.
    int getNumberOfPoints() const;
    /// Returns true if the entire grid is empty.
    bool isEmpty() const;
    /// Returns an iterator over the entire grid (with subgrids)
    iterator begin();
    /// Returns an iterator over specified region.
    iterator beginAt ( const FloatArray &x0, const FloatArray &x1 );

    /**
     * Returns the resolution in given dimension.
     * @param i Index of dimension.
     * @return Resolution in dimension i.
     */
    int getResolution ( int i ) {
        return this->res ( i );
    }
    /**
     * Returns the grid increment.
     * @param i Index of dimension.
     * @return Grid step in dimension i.
     */
    double getGridStep ( int i ) const {
        return this->dx ( i );
    }

    /**
     * Gives the real coordinates for the given identifier.
     * @param[out] answer The coordinate vector.
     * @param[in] pos Position for the grid point.
     */
    void getGridCoord ( FloatArray &answer, const IntArray &pos ) const;
    /**
     * Gives the real coordinates for the given identifier.
     * @param[out] answer The coordinate vector.
     * @param[in] ind Identifier for the grid point.
     */
    void getGridCoord ( FloatArray &answer, int ind ) const;

    /**
     * Gives the point at the position.
     * @param[out] answer The point as requested position
     * @param[in] pos The requested position
     * @return true if point is active, otherwise false
     */
    bool getPoint ( Point *&answer, const IntArray &pos ) const;
    bool getPoint ( Point *&answer, int ind ) const;

    /**
     * Gives the sub-grid at the position.
     * @param[out] subgrid Requested sub-grid.
     * @param[in] pos Position.
     * @return true if sub-grid is active, otherwise false.
     */
    bool getSubGrid ( ParticleGrid<Point> *&subgrid, const IntArray &pos ) const;
    bool getSubGrid ( ParticleGrid<Point> *&subgrid, int ind ) const;

    /**
     * Creates a 2 by 2 sub-grid at the given position.
     * Any existing point or sub-grid is deleted.
     * @param[out] subgrid Newly created sub-grid.
     * @param[in] res Resolution of new sub-grid.
     * @param[in] pos Position of sub-grid.
     */
    void createSubGrid ( ParticleGrid<Point> *&subgrid, const IntArray &res, const IntArray &pos );

    /**
     * Sets a point in the grid. Deletes any preexisting data.
     * @param[in] point Point to set at position pos.
     * @param[in] pos Position.
     */
    void setPoint ( Point *point, const IntArray &pos );
    void setPoint ( Point *point, int ind );

    /**
     * Deletes any preexisting data at position.
     * @param[in] pos Position to clear.
     */
    void clearPosition ( const IntArray &pos );
    void clearPosition ( int ind );

    void getPointsWithin ( std::list<Point*> &answer, const FloatArray &x0, const FloatArray &x1 );
    void appendAllPoints ( std::list<Point*> &answer );

    /**
     * Finds the indices [ind0, ind1) of the bounding box [x0,x1).
     * @param[in] x0 Lower limit of box.
     * @param[in] x1 Upper limit of box.
     * @param[out] ind0 Lower index corresponding to x0.
     * @param[out] ind1 Upper index
     */
    void getBoundingBox ( const FloatArray &x0, const FloatArray &x1, IntArray &ind0, IntArray &ind1 ) const;

    /**
     * Finds the position closest to coordinate x.
     * @param[in] x Coordinate.
     * @param[out] pos Position of closest grid point.
     */
    void getPosition ( const FloatArray &x, IntArray &pos ) const;

    std :: string errorInfo(const char *func) const { return std :: string("ParticleGrid :: ") + func; }

    friend class ParticleGridIterator<Point>;
};

template <class Point>
ParticleGrid<Point> :: ParticleGrid ( const IntArray &res, const FloatArray &bb0, const FloatArray &bb1 )
{
    this->init ( res, bb0, bb1 );
}

template <class Point>
ParticleGrid<Point> :: ParticleGrid ( const ParticleGrid *old )
{
    this->init ( old->res, old->bb0, old->bb1 );
}

template <class Point>
void ParticleGrid<Point> :: init ( const IntArray &res, const FloatArray &bb0, const FloatArray &bb1 )
{
    this->res = res;
    this->bb0 = bb0;
    this->bb1 = bb1;

    this->n = this->res.giveSize();
    this->dx.resize ( this->n );
    this->total = 1;
    this->res_prod.resize ( this->n );
    this->res_prod ( 0 ) = 1;

    for ( int i = 0; i < this->n; ++i ) {
        this->dx ( i ) = ( this->bb1 ( i ) - this->bb0 ( i ) ) / ( this->res ( i )-1 );

        this->total *= this->res ( i );

        this->res_prod ( i ) = 1;
        for ( int j = 0; j < i; ++j ) {
            this->res_prod ( i ) *= this->res ( j );
        }
    }

    this->data = new RefinedParticlePoint*[this->total];
    for ( int i = 0; i < this->total; ++i ) {
        this->data[i] = NULL;
    }
}

template <class Point>
ParticleGrid<Point> :: ~ParticleGrid()
{
    // Clear all points
    for ( int i = 0; i < this->total; ++i ) {
        this->clearPosition ( i );
    }
    delete[] this->data;
}

template <class Point>
int ParticleGrid<Point> :: getNumberOfPoints() const
{
    RefinedParticlePoint *p;
    int c = 0;
    for ( int i = 0; i < this->total; ++i ) {
        p = this->data[i];
        if ( p ) {
            if ( p->point ) {
                c += 1;
            } else if ( p->subgrid ) {
                c += p->subgrid->getNumberOfPoints();
            }
        }
    }
    return c;
}

template <class Point>
bool ParticleGrid<Point> :: isEmpty() const
{
    RefinedParticlePoint *p;
    for ( int i = 0; i < this->total; ++i ) {
        p = this->data[i];
        if ( p ) {
            if ( p->point ) {
                return false;
            } else if ( p->subgrid ) {
                if ( !p->subgrid->isEmpty() ) {
                    return false;
                }
            }
        }
    }
    return true;
}

template <class Point>
int ParticleGrid<Point> :: giveIndex ( const IntArray &pos ) const
{
    int index = 0;
    for ( int i = 0; i < this->n; ++i ) {
        index += this->res_prod ( i ) *pos ( i );
    }
    return index;
}

template <class Point>
void ParticleGrid<Point> :: givePosition ( IntArray &pos, int index ) const
{
    pos.resize ( this->n );
    for ( int i = 0; i < this->n; ++i ) {
        pos ( i ) = ( index/res_prod ( i ) ) % res ( i );
    }
}

template <class Point>
void ParticleGrid<Point> :: getGridCoord ( FloatArray &answer, const IntArray &pos ) const
{
    answer.resize ( this->n );
    for ( int i = 0; i < this->n; ++i ) {
        answer ( i ) = pos ( i ) *this->dx ( i ) + this->bb0 ( i );
    }
}

template <class Point>
void ParticleGrid<Point> :: getGridCoord ( FloatArray &answer, int ind ) const
{
    IntArray pos;
    this->givePosition ( pos, ind );
    this->getGridCoord ( answer, pos );
}

template <class Point>
void ParticleGrid<Point> :: getBoundingBox ( const FloatArray &x0, const FloatArray &x1, IntArray &ind0, IntArray &ind1 ) const
{
    ind0.resize ( this->n );
    ind1.resize ( this->n );
    for ( int i = 0; i < this->n; ++i ) {
        ind0 ( i ) = max ( ( int ) floor ( ( x0 ( i ) - this->bb0 ( i ) ) /this->dx ( i ) ), 0 );
        ind1 ( i ) = min ( ( int ) ceil ( ( x1 ( i ) - this->bb0 ( i ) ) /this->dx ( i ) ), this->res ( i ) );
    }
}

template <class Point>
void ParticleGrid<Point> :: getPosition ( const FloatArray &x, IntArray &pos ) const
{
    pos.resize ( this->n );
    for ( int i = 0; i < this->n; ++i ) {
        //pos(i) = (int)nearest((x(i) - this->bb0(i))/this->dx(i));
        pos ( i ) = ( int ) floor ( ( x ( i ) - this->bb0 ( i ) ) /this->dx ( i ) + 0.5 );
        pos ( i ) = min ( max ( pos ( i ),0 ),this->res ( i )-1 );
    }
}

template <class Point>
typename ParticleGrid<Point>::iterator ParticleGrid<Point> :: begin()
{
    iterator it ( this );
    return it;
}

template <class Point>
typename ParticleGrid<Point>::iterator ParticleGrid<Point> :: beginAt ( const FloatArray &x0, const FloatArray &x1 )
{
    IntArray ind0, ind1;
    this->getBoundingBox ( x0, x1, ind0, ind1 );
    iterator it ( this, ind0, ind1 );
    return it;
}

template <class Point>
void ParticleGrid<Point> :: getPointsWithin ( std::list<Point*> &answer, const FloatArray &x0, const FloatArray &x1 )
{
    answer.clear();

    // First find a suitable bounding box of positions (note: this isn't perfect in regards to sub grids)
    IntArray p0, p1;
    this->getBoundingBox ( x0, x1, p0, p1 );

    // Get all grid points in that region
    IntArray p ( this->n );
    if ( this->n == 2 ) {
        for ( int x = p0 ( 0 ); x < p1 ( 0 ); x++ ) {
            p ( 0 ) = x;
            for ( int y = p0 ( 1 ); y < p1 ( 1 ); y++ ) {
                p ( 1 ) = y;
                RefinedParticlePoint *rp = this->data[this->giveIndex ( p )];
                if ( rp ) {
                    if ( rp->subgrid ) {
                        rp->subgrid->appendAllPoints ( answer );
                    } else {
                        answer.push_front ( rp->point );
                    }
                }
            }
        }
    } else if ( this->n == 3 ) {
        OOFEM_ERROR ( "3D grids not implemented yet" );
    }
}

template <class Point>
void ParticleGrid<Point> :: appendAllPoints ( std::list<Point*> &answer )
{
    for ( int i = 0; i < this->total; ++i ) {
        RefinedParticlePoint *rp = this->data[i];
        if ( rp ) {
            if ( rp->subgrid ) {
                rp->subgrid->appendAllPoints ( answer );
            } else if ( rp->point ) {
                answer.push_front ( rp->point );
            } else {
                OOFEM_ERROR ( "Refined particle does not contain subgrid or point." );
            }
        }
    }
}

template <class Point>
bool ParticleGrid<Point> :: getPoint ( Point *&answer, const IntArray &pos ) const
{
    return this->getPoint ( answer, this->giveIndex ( pos ) );
}

template <class Point>
void ParticleGrid<Point> :: setPoint ( Point *point, const IntArray &pos )
{
    this->setPoint ( point, this->giveIndex ( pos ) );
}

template <class Point>
bool ParticleGrid<Point> :: getPoint ( Point *&answer, int ind ) const
{
    RefinedParticlePoint *rp = this->data[ind];
    if ( rp == NULL ) { // Not active
        answer = NULL;
        return false;
    } else {
        if ( rp->point ) { // Active and has a point
            answer = rp->point;
            return true;
        } else { // Active and has a sub grid
            answer = NULL;
            return false;
        }
    }
}

template <class Point>
void ParticleGrid<Point> :: setPoint ( Point *point, int ind )
{
    this->clearPosition ( ind );
    this->data[ind] = new RefinedParticlePoint;
    this->data[ind]->point = point;
    this->data[ind]->subgrid = NULL;
}

template <class Point>
void ParticleGrid<Point> :: clearPosition ( const IntArray &pos )
{
    return this->clearPosition ( this->giveIndex ( pos ) );
}

template <class Point>
void ParticleGrid<Point> :: clearPosition ( int ind )
{
    RefinedParticlePoint *rp = this->data[ind];
    if ( rp ) {
        if ( rp->point ) {
            delete rp->point;
        }
        if ( rp->subgrid ) {
            delete rp->subgrid;
        }
        delete rp;
    }
    this->data[ind] = NULL;
}

template <class Point>
bool ParticleGrid<Point> :: getSubGrid ( ParticleGrid<Point> *&subgrid, int ind ) const
{
    RefinedParticlePoint *rp = this->data[ind];
    if ( rp == NULL ) {
        subgrid = NULL;
        return false;
    } else {
        if ( rp->subgrid ) {
            subgrid = rp->subgrid;
            return true;
        } else {
            subgrid = NULL;
            return false;
        }
    }
}

template <class Point>
bool ParticleGrid<Point> :: getSubGrid ( ParticleGrid<Point> *&subgrid, const IntArray &pos ) const
{
    return this->getSubGrid ( subgrid, this->giveIndex ( pos ) );
}

template <class Point>
void ParticleGrid<Point> :: createSubGrid ( ParticleGrid<Point> *&subgrid, const IntArray &res, const IntArray &pos )
{
    // The bounding boxes for the sub domains.
    FloatArray grid_point;
    this->getGridCoord ( grid_point, pos );
    FloatArray sub_bb0, sub_bb1;
    sub_bb0 = grid_point;
    sub_bb0.add ( -0.5, dx );
    sub_bb1 = grid_point;
    sub_bb1.add ( 0.5, dx );
    subgrid = new ParticleGrid ( res, sub_bb0, sub_bb1 );

    int ind = this->giveIndex ( pos );
    this->clearPosition ( ind );
    this->data[ind] = new RefinedParticlePoint;
    this->data[ind]->point = NULL;
    this->data[ind]->subgrid = subgrid;
}

/**
 * A recursive iterator for a grid with refinements.
 */
template <class Point>
class ParticleGridIterator
{
protected:
    int index, endind;
    ParticleGridIterator<Point> *sub_it;
    ParticleGrid<Point> *grid;

    IntArray pos;
    const IntArray ind0, ind1;
    bool limited;

public:
    /// Constructor.
    ParticleGridIterator() :
        index ( 0 ), endind ( 0 ), sub_it ( NULL ), grid ( NULL ), ind0(), ind1(), limited ( false ) {}
    /// Constructor.
    ParticleGridIterator ( ParticleGrid<Point> *g ) :
        index ( 0 ), endind ( g->total ), sub_it ( NULL ), grid ( g ), ind0(), ind1(), limited ( false ) {
        if ( this->grid->data[this->index] && this->grid->data[this->index]->subgrid ) {
            this->sub_it = new ParticleGridIterator<Point> ( this->grid->data[this->index]->subgrid );
        }
    }
    /// Constructor.
    ParticleGridIterator ( ParticleGrid<Point> *g, IntArray ind0, IntArray ind1 ) :
        index ( 0 ), sub_it ( NULL ), grid ( g ), ind0 ( ind0 ), ind1 ( ind1 ), limited ( true ) {
        IntArray end = ind1;
        end.add ( -1 );
        this->endind = this->grid->giveIndex ( end ) + 1;
        this->pos = ind0;
        this->index = this->grid->giveIndex ( ind0 );
        if ( this->grid->data[this->index] && this->grid->data[this->index]->subgrid ) {
            this->sub_it = new ParticleGridIterator<Point> ( this->grid->data[this->index]->subgrid );
        }
    }
    /// Copy constructor.
    ParticleGridIterator ( const ParticleGridIterator<Point> &x ) :
        index ( x.index ), sub_it ( x.sub_it ), grid ( x.grid ), ind0 ( x.ind0 ), ind1 ( x.ind1 ), limited ( x.limited )  { }
    /// Destructor
    ~ParticleGridIterator() {
        if ( this->sub_it ) {
            delete sub_it;
        }
    }

    bool end() const {
        return this->index >= this->endind;
    }

    Point *getPoint() {
        if ( this->end() ) {
            return NULL;
        }
        if ( this->subGridActive() ) {
            return this->sub_it->getPoint();
        }
        if ( this->grid->data[this->index] != NULL ) {
            return this->grid->data[this->index]->point;
        }
        return NULL;
    }

    int getIndex() {
        return this->index;
    }

    void setPoint ( Point *p ) {
        if ( this->end() ) {
            OOFEM_ERROR ( "Can't set element in outside grid" );
        }
        if ( this->subGridActive() ) {
            this->sub_it->setPoint ( p );
        }
        this->grid->setPoint ( p, this->index );
    }

    bool subGridActive() const {
        return this->sub_it != NULL;
    }

    void getGridPoint ( FloatArray &x ) {
        if ( this->subGridActive() ) {
            this->sub_it->getGridPoint ( x );
        } else {
            this->grid->getGridCoord ( x, this->index );
        }
    }

    /// Lets the iterator step forward to the next element.
    void operator++() {
        if ( this->subGridActive() ) {
            this->sub_it++;
            if ( !this->sub_it->end() ) {
                return;
            } else {
                delete this->sub_it;
                this->sub_it = NULL;
            }
        }
        this->index++;
        if ( this->limited ) {
            this->pos ( 0 ) ++;
            for ( int i = 0; i < this->pos.giveSize() - 1; ++i ) {
                if ( this->pos ( i ) >= this->ind1 ( i ) ) {
                    this->pos ( i ) = this->ind0 ( i );
                    this->pos ( i+1 ) ++;
                } else {
                    break;
                }
            }
            index = this->grid->giveIndex ( pos );
        }
        if ( !this->end() && this->grid->data[this->index] && this->grid->data[this->index]->subgrid ) {
            this->sub_it = new ParticleGridIterator<Point> ( this->grid->data[this->index]->subgrid );
        }
    }

    std :: string errorInfo(const char *func) const { return std :: string("ParticleGridIterator :: ") + func; }

    friend class ParticleGrid<Point>;
};

} // end namespace oofem
#endif // particlegrid_h
