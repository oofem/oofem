#ifdef __OOFEG


#include "oofeggraphiccontext.h"
#include "mathfem.h"
#include "oofegutils.h"

namespace oofem {
void
oofeg_drawIsoLinesOnTriangle(WCRec coords [ 3 ], double s [ 3 ]) {
    /* Draws the iso-lines on triangle given by its coordinates (coords argument).
     * The values at vertices are stored in s argument.
     * the layer and width are changed, so they should be set before.
     */
    // isoline implementation
    double zcoord, zmin, zmax, minv, maxv;
    int indx, inode, jnode, iside, zlevel;
    int isc_color_scale_num_labels = COLOR_SCALE_NUM_LABELS;
    FPNum delta;
    EFringeTable ft;
    EPixel color;
    WCRec p [ 2 ];
    GraphicObj *tr;

    minv = min( s [ 0 ], min(s [ 1 ], s [ 2 ]) );
    maxv = max( s [ 0 ], max(s [ 1 ], s [ 2 ]) );

    ft = EMGetAssocFringeTable( ESIModel() );
    ColorFringesMinMax(ft, & minv, & maxv);
    delta = ( maxv - minv ) / isc_color_scale_num_labels;
    zcoord = minv + delta / 2;

    for ( zlevel = 1; zlevel <= isc_color_scale_num_labels; zlevel++ ) {
        indx = 0;
        if ( ( zcoord <= minv ) || ( zcoord >= maxv ) ) {
            return;
        }

        for ( iside = 1; iside <= 3; iside++ ) {
            inode = iside;
            jnode = ( iside == 3 ? 1 : iside + 1 );
            zmin = min(s [ inode - 1 ], s [ jnode - 1 ]);
            zmax = max(s [ inode - 1 ], s [ jnode - 1 ]);
            if ( ( zmax > zcoord ) && ( zmin < zcoord ) ) {
                // find intersection
                double ix, iy, iz, jx, jy, jz;
                double t, edgeLength;

                ix = coords [ inode - 1 ].x;
                iy = coords [ inode - 1 ].y;
                iz = 0.;
                jx = coords [ jnode - 1 ].x;
                jy = coords [ jnode - 1 ].y;
                jz = 0.;

                edgeLength = sqrt( ( ix - jx ) * ( ix - jx ) + ( iy - jy ) * ( iy - jy ) );
                t = ( zcoord - s [ inode - 1 ] ) * edgeLength / ( s [ jnode - 1 ] - s [ inode - 1 ] );
                t /= edgeLength;

                p [ indx ].x = ( FPNum ) ix + t * ( jx - ix );
                p [ indx ].y = ( FPNum ) iy + t * ( jy - iy );
                p [ indx ].z = 0.;
                indx++;
                if ( indx == 2 ) {
                    color = ColorFringeRangeToColor( ColorFringeValueToRange(ft, zcoord) );
                    EASValsSetColor(color);

                    tr =  CreateLine3D(p);
                    EGWithMaskChangeAttributes(LAYER_MASK | WIDTH_MASK | COLOR_MASK, tr);
                    EMAddGraphicsToModel(ESIModel(), tr);

                    break;
                }
            }
        }

        zcoord += delta;
    }
}

void
oofeg_drawIsoLinesOnQuad(WCRec coords [ 4 ], double s [ 4 ]) {
    /* Draws the iso-lines on quad given by its coordinates (coords argument).
     * The values at vertices are stored in s argument.
     * the layer and width are changed, so they should be set before.
     */
    // isoline implementation
    double zcoord, zmin, zmax, minv, maxv;
    int i, indx, inode, jnode, iside, zlevel;
    int isc_color_scale_num_labels = COLOR_SCALE_NUM_LABELS;
    FPNum delta;
    EFringeTable ft;
    EPixel color;
    WCRec p [ 4 ];
    GraphicObj *tr;

    minv = maxv = s [ 0 ];
    for ( i = 1; i < 4; i++ ) {
        minv = min(minv, s [ i ]);
        maxv = max(maxv, s [ i ]);
    }

    ft = EMGetAssocFringeTable( ESIModel() );
    ColorFringesMinMax(ft, & minv, & maxv);
    delta = ( maxv - minv ) / isc_color_scale_num_labels;
    zcoord = minv + delta / 2;

    for ( zlevel = 1; zlevel <= isc_color_scale_num_labels; zlevel++ ) {
        indx = 0;
        if ( ( zcoord <= minv ) || ( zcoord >= maxv ) ) {
            return;
        }

        for ( iside = 1; iside <= 4; iside++ ) {
            inode = iside;
            jnode = ( iside == 4 ? 1 : iside + 1 );
            zmin = min(s [ inode - 1 ], s [ jnode - 1 ]);
            zmax = max(s [ inode - 1 ], s [ jnode - 1 ]);
            if ( ( zmax > zcoord ) && ( zmin < zcoord ) ) {
                // find intersection
                double ix, iy, iz, jx, jy, jz;
                double t, edgeLength;

                ix = ( FPNum ) coords [ inode - 1 ].x;
                iy = ( FPNum ) coords [ inode - 1 ].y;
                iz = 0.;
                jx = ( FPNum ) coords [ jnode - 1 ].x;
                jy = ( FPNum ) coords [ jnode - 1 ].y;
                jz = 0.;


                edgeLength = sqrt( ( ix - jx ) * ( ix - jx ) + ( iy - jy ) * ( iy - jy ) );
                t = ( zcoord - s [ inode - 1 ] ) * edgeLength / ( s [ jnode - 1 ] - s [ inode - 1 ] );
                t /= edgeLength;

                p [ indx ].x = ( FPNum ) ix + t * ( jx - ix );
                p [ indx ].y = ( FPNum ) iy + t * ( jy - iy );
                p [ indx ].z = 0.;
                indx++;
                if ( indx == 4 ) {
                    break;
                }
            }
        }

        color = ColorFringeRangeToColor( ColorFringeValueToRange(ft, zcoord) );
        EASValsSetColor(color);

        if ( indx == 2 ) { // only two edges intersecting
            tr =  CreateLine3D(p);
            EGWithMaskChangeAttributes(LAYER_MASK | WIDTH_MASK | COLOR_MASK, tr);
            EMAddGraphicsToModel(ESIModel(), tr);
        } else if ( indx == 4 ) {
            // four edges intersecting
            if ( s [ 0 ] > zcoord ) { // s[2] > zcoord
                WCRec pp [ 2 ];
                pp [ 0 ].x = p [ 0 ].x;
                pp [ 0 ].y = p [ 0 ].y;
                pp [ 0 ].z = p [ 0 ].z;
                pp [ 1 ].x = p [ 3 ].x;
                pp [ 1 ].y = p [ 3 ].y;
                pp [ 1 ].z = p [ 3 ].z;
                tr =  CreateLine3D(pp);
                EGWithMaskChangeAttributes(LAYER_MASK | WIDTH_MASK | COLOR_MASK, tr);
                EMAddGraphicsToModel(ESIModel(), tr);

                pp [ 0 ].x = p [ 1 ].x;
                pp [ 0 ].y = p [ 1 ].y;
                pp [ 0 ].z = p [ 1 ].z;
                pp [ 1 ].x = p [ 2 ].x;
                pp [ 1 ].y = p [ 2 ].y;
                pp [ 1 ].z = p [ 2 ].z;
                tr =  CreateLine3D(pp);
                EGWithMaskChangeAttributes(LAYER_MASK | WIDTH_MASK | COLOR_MASK, tr);
                EMAddGraphicsToModel(ESIModel(), tr);
            } else {
                WCRec pp [ 2 ];
                pp [ 0 ].x = p [ 0 ].x;
                pp [ 0 ].y = p [ 0 ].y;
                pp [ 0 ].z = p [ 0 ].z;
                pp [ 1 ].x = p [ 1 ].x;
                pp [ 1 ].y = p [ 1 ].y;
                pp [ 1 ].z = p [ 1 ].z;
                tr =  CreateLine3D(pp);
                EGWithMaskChangeAttributes(LAYER_MASK | WIDTH_MASK | COLOR_MASK, tr);
                EMAddGraphicsToModel(ESIModel(), tr);

                pp [ 0 ].x = p [ 3 ].x;
                pp [ 0 ].y = p [ 3 ].y;
                pp [ 0 ].z = p [ 3 ].z;
                pp [ 1 ].x = p [ 2 ].x;
                pp [ 1 ].y = p [ 2 ].y;
                pp [ 1 ].z = p [ 2 ].z;
                tr =  CreateLine3D(pp);
                EGWithMaskChangeAttributes(LAYER_MASK | WIDTH_MASK | COLOR_MASK, tr);
                EMAddGraphicsToModel(ESIModel(), tr);
            }
        }

        zcoord += delta;
    }
}
} // end namespace oofem
#endif
