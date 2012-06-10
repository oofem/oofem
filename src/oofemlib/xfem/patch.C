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

#include "flotarry.h"
#include "gausspnt.h"
#include "integrationrule.h"
#include "patch.h"
#include "gaussintegrationrule.h"
#include "delaunay.h"
#include "fei2dtrlin.h"
#include "crosssection.h"
#include "enrichmentitem.h"
#include "datastream.h"
#include "contextioerr.h"
#include "geometry.h"

namespace oofem {
Patch :: Patch(Element *parent) : BasicGeometry()
{
    this->parent = parent;
    this->material = -1;
}

Patch :: Patch(Element *parent, int material) : BasicGeometry()
{
    this->parent = parent;
    this->material = material;
}

Patch :: Patch(Element *parent, AList< FloatArray > *vertices) : BasicGeometry()
{
    this->parent = parent;
    this->vertices = vertices;
}

contextIOResultType
Patch :: saveContext(DataStream *stream, ContextMode mode, void *obj)
{
    contextIOResultType iores;

    if ( stream == NULL ) {
        OOFEM_ERROR("saveContex : can't write into NULL stream");
    }

    // save patch material id
    if ( !stream->write(& material, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    // save patch vertices
    int i, _nvert = vertices->giveSize();
    if ( !stream->write(& _nvert, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    for ( i = 1; i <= _nvert; i++ ) {
        if ( ( iores = vertices->at(i)->storeYourself(stream, mode) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }
    }

    return CIO_OK;
}

contextIOResultType
Patch :: restoreContext(DataStream *stream, ContextMode mode, void *obj)
{
    contextIOResultType iores;

    if ( stream == NULL ) {
        OOFEM_ERROR("restoreContex : can't write into NULL stream");
    }

    // read patch material id
    if ( !stream->read(& material, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }


    // restore patch vertices
    int i, _nvert;
    if ( !stream->read(& _nvert, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    this->vertices->growTo(_nvert);
    for ( i = 1; i <= _nvert; i++ ) {
        FloatArray *_arry = new FloatArray();
        if ( ( iores = _arry->restoreYourself(stream, mode) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        } else {
            this->vertices->put(i, _arry);
        }
    }

    return CIO_OK;
}

FEI2dTrLin TrianglePatch :: interpolation(1, 2);

void TrianglePatch :: convertGPIntoParental(GaussPoint *gp)
{
    FloatArray global;
    const FloatArray **coords = new const FloatArray * [ this->giveNrVertices() ];
    // this we should put into the function before
    for ( int i = 1; i <= this->giveNrVertices(); i++ ) {
        coords [ i - 1 ] = new FloatArray( *this->giveVertex(i) );
    }

    this->interpolation.local2global(global, * gp->giveCoordinates(),
                                     FEIVertexListGeometryWrapper(this->giveNrVertices(), coords));
    for ( int i = 1; i <= this->giveNrVertices(); i++ ) {
        delete coords [ i - 1 ];
    }

    delete [] coords;
    FloatArray local;
    parent->computeLocalCoordinates(local, global);
    gp->setCoordinates(local);
}


#ifdef __OOFEG
void
TrianglePatch :: draw(oofegGraphicContext &gc)
{
    WCRec p [ 3 ];
    GraphicObj *go;

    EASValsSetLineWidth(OOFEG_RAW_GEOMETRY_WIDTH);
    EASValsSetColor( gc.getElementColor() );
    EASValsSetEdgeColor( gc.getGeometryColor() );
    EASValsSetEdgeFlag(true);
    EASValsSetLayer(OOFEG_RAW_GEOMETRY_LAYER);
    p [ 0 ].x = ( FPNum ) vertices->at(1)->at(1);
    p [ 0 ].y = ( FPNum ) vertices->at(1)->at(2);
    p [ 0 ].z = 0.;
    p [ 1 ].x = ( FPNum ) vertices->at(2)->at(1);
    p [ 1 ].y = ( FPNum ) vertices->at(2)->at(2);
    p [ 1 ].z = 0.;
    p [ 2 ].x = ( FPNum ) vertices->at(3)->at(1);
    p [ 2 ].y = ( FPNum ) vertices->at(3)->at(2);
    p [ 2 ].z = 0.;

    go =  CreateTriangle3D(p);
    EGWithMaskChangeAttributes(WIDTH_MASK | COLOR_MASK | EDGE_COLOR_MASK | EDGE_FLAG_MASK | LAYER_MASK, go);
    EGAttachObject(go, ( EObjectP ) this);
    EMAddGraphicsToModel(ESIModel(), go);
}

void
TrianglePatch :: drawWD(oofegGraphicContext &gc, FloatArray &vd)
{
    WCRec p [ 3 ];
    double s [ 3 ];
    GraphicObj *go;

    EASValsSetLayer(OOFEG_VARPLOT_PATTERN_LAYER);
    EASValsSetColor( gc.getElementColor() );
    EASValsSetEdgeColor( gc.getGeometryColor() );
    EASValsSetEdgeFlag(true);
    p [ 0 ].x = ( FPNum ) vertices->at(1)->at(1);
    p [ 0 ].y = ( FPNum ) vertices->at(1)->at(2);
    p [ 0 ].z = 0.;
    p [ 1 ].x = ( FPNum ) vertices->at(2)->at(1);
    p [ 1 ].y = ( FPNum ) vertices->at(2)->at(2);
    p [ 1 ].z = 0.;
    p [ 2 ].x = ( FPNum ) vertices->at(3)->at(1);
    p [ 2 ].y = ( FPNum ) vertices->at(3)->at(2);
    p [ 2 ].z = 0.;

    if ( vd.giveSize() == 3 ) {
        s [ 0 ] = vd.at(1);
        s [ 1 ] = vd.at(2);
        s [ 2 ] = vd.at(3);
    }

    go =  CreateTriangleWD3D(p, s [ 0 ], s [ 1 ], s [ 2 ]);
    EGWithMaskChangeAttributes(WIDTH_MASK | COLOR_MASK | EDGE_COLOR_MASK | EDGE_FLAG_MASK | LAYER_MASK, go);
    EGAttachObject(go, ( EObjectP ) this);
    EMAddGraphicsToModel(ESIModel(), go);
}

#endif
} // end namespace oofem
