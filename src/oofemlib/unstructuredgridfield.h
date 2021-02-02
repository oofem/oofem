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

#ifndef unstructuredgridfield_h
#define unstructuredgridfield_h

#include "field.h"
#include "floatarray.h"
#include "intarray.h"
#include "elementgeometrytype.h"
#include "feinterpol.h"
#include "fei2dlinelin.h"
#include "fei2dlinequad.h"
#include "fei2dtrlin.h"
#include "fei2dtrquad.h"
#include "fei2dquadlin.h"
#include "fei2dquadquad.h"
#include "fei3dtetlin.h"
#include "fei3dhexalin.h"
#include "octreelocalizert.h"
#include "error.h"

namespace oofem {
/**
 * Field defined by values defined on unstructured grid.
 * This implementation does not require that underlying grid is oofem mesh (composed of DofMagars and Elements),
 * but rather simple objects (Vertices and Cells) which in turn does not have to provide complex FE services.
 * Typically used to map extermal field to mupif.
 * The performance is robust but quite slow due to excessive searching. For external fields, use rather DofManValueField.
 */
class OOFEM_EXPORT UnstructuredGridField : public Field
{
protected:
    class Vertex
    {
        FloatArray coords;
public:
        Vertex() : coords() {}
        Vertex(FloatArray &c) { coords = c; }
        Vertex &operator=(const Vertex &v) {
            coords = v.coords;
            return * this;
        }
        const FloatArray &getCoordinates() const { return coords; }
    };

    /* Auxiliary class to represent interpolation cell, basically holding vertices defining its shape and interpolation type */
    class Cell
    {
        Element_Geometry_Type itype;
        IntArray vertices;
        UnstructuredGridField *mesh;

        // array of interpolation instances used for supported ElementGeometryTypes
        static FEI2dLineLin i1;
        static FEI2dLineQuad i2;
        static FEI2dTrLin i3;
        static FEI2dTrQuad i4;
        static FEI2dQuadLin i5;
        static FEI2dQuadQuad i6;
        static FEI3dTetLin i7;
        static FEI3dHexaLin i8;

        static FEInterpolation *interpTable[];


        /**
         * Wrapper around element definition to provide FEICellGeometry interface.
         */
        class FEICellGeometryWrapper : public FEICellGeometry
        {
protected:
            const Cell *cell;
public:
            FEICellGeometryWrapper(const Cell *c) : FEICellGeometry() {
                this->cell = c;
            }
            virtual ~FEICellGeometryWrapper() { }

            int giveNumberOfVertices() const override {
                return this->cell->giveNumberOfVertices();
            }

            const FloatArray &giveVertexCoordinates(int i) const override
            {
                return ( ( cell->getVertex(i) )->getCoordinates() );
            }
        };

public:
        Cell() {
            itype = EGT_unknown;
            mesh = NULL;
        }
        Cell(Element_Geometry_Type t, IntArray &v, UnstructuredGridField *m) {
            itype = t;
            vertices = v;
            mesh = m;
        }

        Cell &operator=(const Cell &c) {
            itype = c.itype;
            vertices = c.vertices;
            mesh = c.mesh;
            return * this;
        }

        int giveNumberOfVertices() const { return vertices.giveSize(); }
        bool containsPoint(const FloatArray &coords) const {
            FloatArray tmp;
            BoundingBox b;
            this->giveBoundingBox(b);
            if ( b.contains(coords) ) {
                const FEICellGeometryWrapper cgw = FEICellGeometryWrapper(this);
                return ( this->getInterpolation()->global2local(tmp, coords, cgw) );
            }
            return false;
        }

        void giveBoundingBox(BoundingBox &bb) const {
            double size;
            FloatArray bb0, bb1;
            bb1 = bb0 = this->getVertex(1)->getCoordinates();

            for ( int i = 2; i <= this->giveNumberOfVertices(); i++ ) {
                const auto &coordinates = this->getVertex(i)->getCoordinates();
                bb0.beMinOf(bb0, coordinates);
                bb1.beMaxOf(bb1, coordinates);
            }
            bb1.subtract(bb0);
            int nsd = bb1.giveSize();
            IntArray mask(3);
            size = bb1.at(1);
            for ( int i = 1; i < nsd; i++ ) { size = max(size, bb1 [ i ]); }
            for ( int i = 0; i < nsd; i++ ) { mask [ i ] = 1; }
            bb.init(bb0, size, mask);
        }

        double giveClosestPoint(FloatArray &lcoords, FloatArray &closest, const FloatArray &gcoords) {
            FEInterpolation *interp = this->getInterpolation();

            if ( !interp->global2local(lcoords, gcoords, FEICellGeometryWrapper(this) ) ) { // Outside element
                interp->local2global(closest, lcoords, FEICellGeometryWrapper(this) );
                return distance(closest, gcoords);
            } else {
                closest = gcoords;
                return 0.0;
            }
        }

        FEInterpolation *getInterpolation() const {
            if ( this->itype == EGT_line_1 ) { return interpTable [ 0 ]; } else if ( this->itype == EGT_line_2 ) { return interpTable [ 1 ]; } else if ( this->itype == EGT_triangle_1 ) { return interpTable [ 2 ]; } else if ( this->itype == EGT_triangle_2 ) { return interpTable [ 3 ]; } else if ( this->itype == EGT_quad_1 ) { return interpTable [ 4 ]; } else if ( this->itype == EGT_quad_2 ) { return interpTable [ 5 ]; } else if ( this->itype == EGT_tetra_1 ) { return interpTable [ 6 ]; } else if ( this->itype == EGT_hexa_1 ) { return interpTable [ 7 ]; } else                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 {
                OOFEM_LOG_ERROR("UnstructuredGridField.Cell:: Unsupported cell type");
                return NULL;
            }
        }
        const Vertex *getVertex(int i) const {//1-based
            return mesh->getVertex(vertices [ i - 1 ]);
        }

        int interpolate(FloatArray &answer, const FloatArray &pos, FloatArray **vertexVals) {
            int i, j, size = vertexVals [ 0 ]->giveSize();
            FloatArray N, lcoords;

            FEInterpolation *it = this->getInterpolation();
            FEICellGeometryWrapper cw(this);
            it->global2local(lcoords, pos, cw);
            it->evalN(N, lcoords, cw);

            answer.resize(size);
            answer.zero();
            for ( i = 0; i < giveNumberOfVertices(); i++ ) {
                for ( j = 0; j < size; j++ ) {
                    answer(j) += N(i) * vertexVals [ i ]->at(j + 1);
                }
            }
            return 1;
        }

        int getVertexNum(int i) { //1-based
            return this->vertices [ i - 1 ];
        }
    };


    class CellInsertionFunctor : public SL_Insertion_Functor< Cell >
    {
private:
public:
        CellInsertionFunctor() {}
        bool evaluate(Cell &member, OctantRecT< Cell > *cell) override {
            BoundingBox b;
            member.giveBoundingBox(b);
            OctantRec::BoundingBoxStatus s = cell->testBoundingBox(b);
            if ( ( s == OctantRec::BBS_InsideCell ) || ( s == OctantRec::BBS_ContainsCell ) ) {
                return true;
            } else {
                return false;
            }
        }

        void registerInsertion(Cell &member, LocalInsertionData< Cell >lidata) override {}
        std::list< LocalInsertionData< Cell > > *giveInsertionList(Cell &m) override { return nullptr; }
    };

    class CellContainingPointFunctor : public SL_Evaluation_Functor< Cell >
    {
protected:
        FloatArray position;
        std::list< Cell >cells;
public:
        /**
         * Constructor
         * @param pos Starting position of the search
         * @param d Domain containing nodes the triangles are defined by
         */
        CellContainingPointFunctor(const FloatArray &pos) : position(pos) {}
        ~CellContainingPointFunctor() { }

        /**
         * Evaluates a triangle upon its circumscribed cricle.
         * @param DTptr Delaunay triangle. nodeNr Number of the node in the domain list
         * @returns true if the circumscribed circle of the Delaunay triangle contains the point, false otherwise
         */
        bool evaluate(Cell &c) override
        {
            if ( c.containsPoint(position) ) {
                cells.push_back(c);
                return true;
            } else {
                return false;
            }
        }

        /**
         * Gives the starting position of the search
         * @param position startingPosition
         */
        void giveStartingPosition(FloatArray &answer) override
        {
            answer = position;
        }

        /**
         * Gives the triangles containing the node
         * @param answer List containing Delaunay triangles
         */
        void giveResult(std::list< Cell > &answer) override
        { answer = cells; }

        bool isBBXStage1Defined(BoundingBox &BBXStage1) override { return true; }
        bool isBBXStage2Defined(BoundingBox &BBXStage2) override { return false; }
    };


    std::vector< Cell >cellList;
    std::vector< Vertex >vertexList;
    std::vector< FloatArray >valueList;
    /**
     * Spatial Localizer. It is build upon request.
     * Provides the spatial localization services.
     */
    OctreeSpatialLocalizerT< Cell >spatialLocalizer;
    /// octree build time stamp
    long int octreeTimeStamp;
    /// receiver timestamp
    long int timeStamp;
    /// octree origin shift
    double octreeOriginShift;
public:
    /**
     * Constructor. Creates a field, with unspecified field values.
     */
    UnstructuredGridField(int nvert, int ncells, double octreeOriginShift = 0.0) : Field(FieldType::FT_Unknown), spatialLocalizer()
    {
        this->timeStamp = this->octreeTimeStamp = 0;
        this->vertexList.resize(nvert);
        this->cellList.resize(ncells);
        this->valueList.resize(nvert);
        this->octreeOriginShift = octreeOriginShift;
    }
    virtual ~UnstructuredGridField() { }

    void setVertexValue(int num, const FloatArray &vv) {
        valueList [ num - 1 ] = vv;
    }

    void addVertex(int num, FloatArray &coords) {//1-based
        vertexList [ num - 1 ] = Vertex(coords);
        this->timeStamp++;
    }

    const FloatArray &getVertexCoordinates(int num) const
    {
        return ( this->vertexList [ num - 1 ].getCoordinates() );
    }

    Vertex *getVertex(int num) { //1-based
        return & this->vertexList [ num - 1 ];
    }

    void addCell(int num, Element_Geometry_Type type, IntArray &vertices) { //1-based
        cellList [ num - 1 ] = Cell(type, vertices, this);
        this->timeStamp++;
    }

    int evaluateAt(FloatArray &answer, const FloatArray &coords,
                   ValueModeType mode, TimeStep *tStep) override {
        std::list< Cell >elist;
        if ( ( mode == VM_Total ) || ( mode == VM_TotalIntrinsic ) ) {
            if ( this->cellList.size() > 0 ) {
                this->initOctree();
                CellContainingPointFunctor f(coords);
                this->spatialLocalizer.giveDataOnFilter(elist, f);
                if ( elist.size() ) {
                    Cell &c = elist.front(); // take first
                    // colect vertex values
                    int size = c.giveNumberOfVertices();
                    FloatArray **vertexValues = new FloatArray * [ size ];
                    for ( int i = 0; i < size; i++ ) {
                        vertexValues [ i ] = & ( this->valueList [ c.getVertexNum(i + 1) - 1 ] );
                    }
                    c.interpolate(answer, coords, vertexValues);
                    //answer.printYourself();
                    return 0;
                } else {
                    coords.printYourself();
                    OOFEM_ERROR("No cells defined");
                    return 1;
                }
            } else {
                if ( !this->vertexList.size() ) {
                    OOFEM_ERROR("No vertices defined");
                    return 1;
                }
                // alternative method - we use the value of the closest point
                double minDist = 0., dist = 0.;

                int idOfClosestPoint = -1;
                for ( int i = 0; i < ( int ) this->vertexList.size(); i++ ) {
                    const auto &pcoords = this->vertexList [ i ].getCoordinates();
                    dist = sqrt(pow(coords [ 0 ] - pcoords.at(1), 2) + pow(coords [ 1 ] - pcoords.at(2), 2) + pow(coords [ 2 ] - pcoords.at(3), 2) );
                    if ( ( dist < minDist ) || ( !i ) ) {
                        minDist = dist;
                        idOfClosestPoint = i;
                    }
                }
                answer = this->valueList [ idOfClosestPoint ];
                //printf("closest point is %d\n",idOfClosestPoint);
                return 0;
            }
        } else {
            OOFEM_ERROR("Unsupported ValueModeType");
        }
    }

    /**
     * Implementaton of Field::evaluateAt for DofManager.
     */
    int evaluateAt(FloatArray &answer, DofManager *dman,
                   ValueModeType mode, TimeStep *tStep) override {
        const auto &coords = dman->giveCoordinates();
        return this->evaluateAt(answer, coords, mode, tStep);
    }

    void saveContext(DataStream &stream) override { }
    void restoreContext(DataStream &stream) override { }

    /// @return Class name of the receiver.
    const char *giveClassName() const override { return "UnstructuredGridField"; }

protected:
    void initOctree() {
        if ( this->timeStamp != this->octreeTimeStamp ) {
            // rebuild octree
            this->spatialLocalizer.clear();
            // get octree bbox
            std::vector< Vertex >::iterator it = vertexList.begin();
            FloatArray cmax, cmin;
            cmax = cmin = ( * it ).getCoordinates();
            int nsd  = cmin.giveSize();
            ++it;
            for ( ; it != vertexList.end(); ++it ) {
                const auto &vc = ( * it ).getCoordinates();
                for ( int j = 0; j < nsd; j++ ) {
                    cmax(j) = max(cmax(j), vc(j) );
                    cmin(j) = min(cmin(j), vc(j) );
                }
            }
            // introduce origin shift (if set up)
            for ( int j = 0; j < nsd; j++ ) {
                cmin(j) -= octreeOriginShift;
            }
            double size = 0.0;
            IntArray mask(3);
            for ( int j = 0; j < nsd; j++ ) {
                size = max(size, cmax(j) - cmin(j) );
                mask [ j ] = 1;
            }
            BoundingBox bb;
            bb.init(cmin, size, mask);
            this->spatialLocalizer.init(bb);
            std::vector< Cell >::iterator cit;
            CellInsertionFunctor cf;
            //int cellcount = 1;
            for ( cit = this->cellList.begin(); cit != this->cellList.end(); ++cit ) {
                /*
                 * BoundingBox b; (*cit).giveBoundingBox(b);
                 * FloatArray o;
                 * b.giveOrigin(o);
                 * printf ("Inserting cell %d, Origin, Size:", cellcount++);
                 * o.printYourself();
                 * printf("%lf\n", b.giveSize());
                 */
                this->spatialLocalizer.insertMemberIntoOctree(* cit, cf);
            }

            // update timestamp
            this->octreeTimeStamp = this->timeStamp;
        }
    }
};
} // end namespace oofem
#endif // unstructuredgridfield_h
