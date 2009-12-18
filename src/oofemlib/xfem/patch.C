#include "flotarry.h"
#include "gausspnt.h"
#include "integrationrule.h"
#include "patch.h"
#include "gaussintegrationrule.h"
#include "delaunay.h"
#include <math.h>
#include "fei2dtrlin.h"
#include "crosssection.h"
#include "enrichmentitem.h"

namespace oofem {

Patch :: Patch(Element *parent, int material) : BasicGeometry() {
    this->parent = parent;
    this->material = material;
}

Patch :: Patch(Element *parent, AList< FloatArray > *vertices) : BasicGeometry() {
    this->parent = parent;
    this->vertices = vertices;
}

FEI2dTrLin TrianglePatch :: interpolation(1, 2);

void TrianglePatch :: convertGPIntoParental(GaussPoint *gp) {
    FloatArray global;
    const FloatArray **coords = new const FloatArray * [ this->giveNrVertices() ];
    // this we should put into the function before
    for ( int i = 1; i <= this->giveNrVertices(); i++ ) {
        coords [ i - 1 ] = new FloatArray( * this->giveVertex(i) );
    }

    this->interpolation.local2global(global, coords, * gp->giveCoordinates(), 1.0);

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
    EASValsSetEdgeFlag(TRUE);
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
    EASValsSetEdgeFlag(TRUE);
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
