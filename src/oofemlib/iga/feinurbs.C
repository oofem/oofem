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

#include "feinurbs.h"
#include "flotarry.h"
#include "flotmtrx.h"
#include "iga.h"

namespace oofem {
// optimized version of A4.4 for d=1
#define OPTIMIZED_VERSION_A4dot4


NURBSInterpolation :: ~NURBSInterpolation() { }


void NURBSInterpolation :: evalN(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    FEIIGAElementGeometryWrapper *gw = ( FEIIGAElementGeometryWrapper * ) & cellgeo;
    IntArray span(nsd);
    double sum = 0.0, val;
    int count, c = 1, i, l, k, m, ind, indx, uind, vind, tind;
#ifdef HAVE_VARIABLE_ARRAY_SIZE
    FloatArray N [ nsd ];
#else
    FloatArray *N = new FloatArray [ nsd ];
#endif

    if ( gw->knotSpan ) {
        span = * gw->knotSpan;
    } else {
        for ( i = 0; i < nsd; i++ ) {
            span(i) = this->findSpan(numberOfControlPoints [ i ], degree [ i ], lcoords(i), knotVector [ i ]);
        }
    }

    for ( i = 0; i < nsd; i++ ) {
        this->basisFuns(N [ i ], span(i), lcoords(i), degree [ i ], knotVector [ i ]);
    }

    count = giveNumberOfKnotSpanBasisFunctions(span);
    answer.resize(count);

    if ( nsd == 1 ) {
        uind = span(0) - degree [ 0 ];
        ind = uind + 1;
        for ( k = 0; k <= degree [ 0 ]; k++ ) {
            answer.at(c++) = val = N [ 0 ](k) * cellgeo.giveVertexCoordinates(ind + k)->at(2);       // Nu*w
            sum += val;
        }
    } else if ( nsd == 2 ) {
        uind = span(0) - degree [ 0 ];
        vind = span(1) - degree [ 1 ];
        ind = vind * numberOfControlPoints [ 0 ] + uind + 1;
        for ( l = 0; l <= degree [ 1 ]; l++ ) {
            for ( k = 0; k <= degree [ 0 ]; k++ ) {
                answer.at(c++) = val = N [ 0 ](k) * N [ 1 ](l) * cellgeo.giveVertexCoordinates(ind + k)->at(3); // Nu*Nv*w
                sum += val;
            }

            ind += numberOfControlPoints [ 0 ];
        }
    } else if ( nsd == 3 ) {
        uind = span(0) - degree [ 0 ];
        vind = span(1) - degree [ 1 ];
        tind = span(2) - degree [ 2 ];
        ind = tind * numberOfControlPoints [ 0 ] * numberOfControlPoints [ 1 ] + vind * numberOfControlPoints [ 0 ] + uind + 1;
        for ( m = 0; m <= degree [ 2 ]; m++ ) {
            indx = ind;
            for ( l = 0; l <= degree [ 1 ]; l++ ) {
                for ( k = 0; k <= degree [ 0 ]; k++ ) {
                    answer.at(c++) = val = N [ 0 ](k) * N [ 1 ](l) * N [ 2 ](m) * cellgeo.giveVertexCoordinates(ind + k)->at(4);     // Nu*Nv*Nt*w
                    sum += val;
                }

                ind += numberOfControlPoints [ 0 ];
            }

            ind = indx + numberOfControlPoints [ 0 ] * numberOfControlPoints [ 1 ];
        }
    } else {
        OOFEM_ERROR2("evalN not implemented for nsd = %d", nsd);
    }

    while ( count ) {
        answer.at(count--) /= sum;
    }

#ifndef HAVE_VARIABLE_ARRAY_SIZE
    delete [] N;
#endif
}



void NURBSInterpolation :: evaldNdx(FloatMatrix &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    FEIIGAElementGeometryWrapper *gw = ( FEIIGAElementGeometryWrapper * ) & cellgeo;
    const FloatArray *vertexCoordsPtr;
    FloatMatrix jacobian(nsd, nsd);
    IntArray span(nsd);
    double Jacob, product, w, weight;
    int count, cnt, i, l, k, m, ind, indx, uind, vind, tind;
#ifdef HAVE_VARIABLE_ARRAY_SIZE
    FloatMatrix ders [ nsd ];
#else
    FloatMatrix *ders = new FloatMatrix [ nsd ];
#endif

    if ( gw->knotSpan ) {
        span = * gw->knotSpan;
    } else {
        for ( i = 0; i < nsd; i++ ) {
            span(i) = this->findSpan(numberOfControlPoints [ i ], degree [ i ], lcoords(i), knotVector [ i ]);
        }
    }

    for ( i = 0; i < nsd; i++ ) {
        this->dersBasisFuns(1, lcoords(i), span(i), degree [ i ], knotVector [ i ], ders [ i ]);
    }

    count = giveNumberOfKnotSpanBasisFunctions(span);
    answer.resize(count, nsd);

#if 0                       // code according NURBS book (too general allowing higher derivatives)
    if ( nsd == 2 ) {
        FloatArray tmp1(nsd + 1), tmp2(nsd + 1);    // allow for weight

        FloatMatrix Aders [ 2 ];      // derivatives in each coordinate direction on BSpline
 #ifndef OPTIMIZED_VERSION_A4dot4
        FloatMatrix Sders [ 2 ];      // derivatives in each coordinate direction on NURBS
 #endif
        FloatMatrix wders;              // derivatives in w direction on BSpline
        /*
         * IntArray Bin(2,2);      // binomial coefficients from 0 to d=1
         *                            // Bin(n,k)=(n above k)=n!/k!(n-k)! for n>=k>=0
         *                                                                                              // lower triangle corresponds to Pascal triangle
         *                                                                                              // according to A4.4 it seems that only coefficients in lower triangle except the first column are used
         */
        // resizing to (2,2) has nothing common with nsd
        // it is related to the fact that 0th and 1st derivatives are computed in each direction
        for ( i = 0; i < nsd; i++ ) {
            Aders [ i ].resize(2, 2);
            Aders [ i ].zero();
 #ifndef OPTIMIZED_VERSION_A4dot4
            Sders [ i ].resize(2, 2);
 #endif
        }

        wders.resize(2, 2);
        wders.zero();

        // calculation of jacobian matrix according to A4.4
        // calculate values and derivatives of nonrational Bspline surface with weights at first (Aders, wders)
        uind = span(0) - degree [ 0 ];
        vind = span(1) - degree [ 1 ];
        ind = vind * numberOfControlPoints [ 0 ] + uind + 1;
        for ( l = 0; l <= degree [ 1 ]; l++ ) {
            tmp1.zero();
            tmp2.zero();
            for ( k = 0; k <= degree [ 0 ]; k++ ) {
                vertexCoordsPtr = cellgeo.giveVertexCoordinates(ind + k);
                w = vertexCoordsPtr->at(3);

                tmp1(0) += ders [ 0 ](0, k) * vertexCoordsPtr->at(1) * w; // sum(Nu*x*w)
                tmp1(1) += ders [ 0 ](0, k) * vertexCoordsPtr->at(2) * w; // sum(Nu*y*w)
                tmp1(2) += ders [ 0 ](0, k) * w;           // sum(Nu*w)

                tmp2(0) += ders [ 0 ](1, k) * vertexCoordsPtr->at(1) * w; // sum(dNu/du*x*w)
                tmp2(1) += ders [ 0 ](1, k) * vertexCoordsPtr->at(2) * w; // sum(dNu/du*y*w)
                tmp2(2) += ders [ 0 ](1, k) * w;           // sum(dNu/du*w)
            }

            ind += numberOfControlPoints [ 0 ];

            Aders [ 0 ](0, 0) += ders [ 1 ](0, l) * tmp1(0); // xw=sum(Nv*sum(Nu*x*w))
            Aders [ 1 ](0, 0) += ders [ 1 ](0, l) * tmp1(1); // yw=sum(Nv*sum(Nu*y*w))
            wders(0, 0)    += ders [ 1 ](0, l) * tmp1(2); // w=sum(Nv*sum(Nu*w))

            Aders [ 0 ](0, 1) += ders [ 1 ](1, l) * tmp1(0); // dxw/dv=sum(dNv/dv*sum(Nu*x*w))
            Aders [ 1 ](0, 1) += ders [ 1 ](1, l) * tmp1(1); // dyw/dv=sum(dNv/dv*sum(Nu*y*w))
            wders(0, 1)    += ders [ 1 ](1, l) * tmp1(2); // dw/dv=sum(dNv/dv*sum(Nu*w))

            Aders [ 0 ](1, 0) += ders [ 1 ](0, l) * tmp2(0); // dxw/du=sum(Nv*sum(dNu/du*x*w))
            Aders [ 1 ](1, 0) += ders [ 1 ](0, l) * tmp2(1); // dyw/du=sum(Nv*sum(dNu/du*y*w))
            wders(1, 0)    += ders [ 1 ](0, l) * tmp2(2);       // dw/du=sum(Nv*sum(dNu/du*w))
        }

        weight = wders(0, 0);

 #ifndef OPTIMIZED_VERSION_A4dot4
        int j;
        const int d = 1;
        // calculate values and derivatives of NURBS surface (A4.4)
        // since all entries in Pascal triangle up to d=1 are 1, binomial coefficients are ignored
        for ( k = 0; k <= d; k++ ) {
            for ( l = 0; l <= d - k; l++ ) {
                tmp1(0) = Aders [ 0 ](k, l);
                tmp1(1) = Aders [ 1 ](k, l);
                for ( j = 1; j <= l; j++ ) {
                    tmp1(0) -= wders(0, j) * Sders [ 0 ](k, l - j);            // *Bin(l,j)
                    tmp1(1) -= wders(0, j) * Sders [ 1 ](k, l - j);            // *Bin(l,j)
                }

                for ( i = 1; i <= k; i++ ) {
                    tmp1(0) -= wders(i, 0) * Sders [ 0 ](k - i, l);            // *Bin(k,i)
                    tmp1(1) -= wders(i, 0) * Sders [ 1 ](k - i, l);            // *Bin(k,i)
                    tmp2.zero();
                    for ( j = 1; j <= l; j++ ) {
                        tmp2(0) += wders(i, j) * Sders [ 0 ](k - i, l - j);              // *Bin(l,j)
                        tmp2(1) += wders(i, j) * Sders [ 1 ](k - i, l - j);              // *Bin(l,j)
                    }

                    tmp1(0) -= tmp2(0);                     // *Bin(k,i)
                    tmp1(1) -= tmp2(1);                     // *Bin(k,i)
                }

                Sders [ 0 ](k, l) = tmp1(0) / weight;
                Sders [ 1 ](k, l) = tmp1(1) / weight;
            }
        }

        jacobian(0, 0) = Sders [ 0 ](1, 0);      // dx/du
        jacobian(0, 1) = Sders [ 1 ](1, 0);      // dy/du
        jacobian(1, 0) = Sders [ 0 ](0, 1);      // dx/dv
        jacobian(1, 1) = Sders [ 1 ](0, 1);      // dy/dv
 #else
        // optimized version of A4.4 for d=1, binomial coefficients ignored

        /*
         * // k=0 l=0 loop
         * Sders[0](0,0) = Aders[0](0,0) / weight;
         * Sders[1](0,0) = Aders[1](0,0) / weight;
         * // k=1 l=0 loop
         * Sders[0](1,0) = (Aders[0](1,0)-wders(1,0)*Sders[0](0,0)) / weight;
         * Sders[1](1,0) = (Aders[1](1,0)-wders(1,0)*Sders[1](0,0)) / weight;
         * // k=0 l=1 loop
         * Sders[0](0,1) = (Aders[0](0,1)-wders(0,1)*Sders[0](0,0)) / weight;
         * Sders[1](0,1) = (Aders[1](0,1)-wders(0,1)*Sders[1](0,0)) / weight;
         *
         * jacobian(0,0) = Sders[0](1,0);   // dx/du
         * jacobian(0,1) = Sders[1](1,0);   // dy/du
         * jacobian(1,0) = Sders[0](0,1);   // dx/dv
         * jacobian(1,1) = Sders[1](0,1);   // dy/dv
         */

        // k=0 l=0 loop
        tmp1(0) = Aders [ 0 ](0, 0) / weight;
        tmp1(1) = Aders [ 1 ](0, 0) / weight;
        // k=1 l=0 loop
        jacobian(0, 0) = ( Aders [ 0 ](1, 0) - wders(1, 0) * tmp1(0) ) / weight; // dx/du
        jacobian(0, 1) = ( Aders [ 1 ](1, 0) - wders(1, 0) * tmp1(1) ) / weight; // dy/du
        // k=0 l=1 loop
        jacobian(1, 0) = ( Aders [ 0 ](0, 1) - wders(0, 1) * tmp1(0) ) / weight; // dx/dv
        jacobian(1, 1) = ( Aders [ 1 ](0, 1) - wders(0, 1) * tmp1(1) ) / weight; // dy/dv
 #endif

        Jacob = jacobian.giveDeterminant();

        //calculation of derivatives of NURBS basis functions with respect to local parameters is not covered by NURBS book
        product = Jacob * weight * weight;
        cnt = 0;
        ind = vind * numberOfControlPoints [ 0 ] + uind + 1;
        for ( l = 0; l <= degree [ 1 ]; l++ ) {
            for ( k = 0; k <= degree [ 0 ]; k++ ) {
                w = cellgeo.giveVertexCoordinates(ind + k)->at(3);
                // dNu/du*Nv*w*sum(Nv*Nu*w) - Nu*Nv*w*sum(dNu/du*Nv*w)
                tmp1(0) = ders [ 0 ](1, k) * ders [ 1 ](0, l) * w * weight - ders [ 0 ](0, k) * ders [ 1 ](0, l) * w * wders(1, 0);
                // Nu*dNv/dv*w*sum(Nv*Nu*w) - Nu*Nv*w*sum(Nu*dNv/dv*w)
                tmp1(1) = ders [ 0 ](0, k) * ders [ 1 ](1, l) * w * weight - ders [ 0 ](0, k) * ders [ 1 ](0, l) * w * wders(0, 1);

                answer(cnt, 0) = ( +jacobian(1, 1) * tmp1(0) - jacobian(0, 1) * tmp1(1) ) / product;
                answer(cnt, 1) = ( -jacobian(1, 0) * tmp1(0) + jacobian(0, 0) * tmp1(1) ) / product;
                cnt++;
            }

            ind += numberOfControlPoints [ 0 ];
        }
    } else {
        OOFEM_ERROR2("evaldNdx not implemented for nsd = %d", nsd);
    }

#else
 #ifdef HAVE_VARIABLE_ARRAY_SIZE
    FloatArray Aders [ nsd ];  // 0th and 1st derivatives in each coordinate direction on BSpline
 #else
    FloatArray *Aders = new FloatArray [ nsd ];
 #endif
    FloatArray wders;          // 0th and 1st derivatives in w direction on BSpline

    for ( i = 0; i < nsd; i++ ) {
        Aders [ i ].resize(nsd + 1);
        Aders [ i ].zero();
    }

    wders.resize(nsd + 1);
    wders.zero();

    if ( nsd == 1 ) {
        // calculate values and derivatives of nonrational Bspline curve with weights at first (Aders, wders)
        uind = span(0) - degree [ 0 ];
        ind = uind + 1;
        for ( k = 0; k <= degree [ 0 ]; k++ ) {
            vertexCoordsPtr = cellgeo.giveVertexCoordinates(ind + k);
            w = vertexCoordsPtr->at(2);

            Aders [ 0 ](0) += ders [ 0 ](0, k) * vertexCoordsPtr->at(1) * w;   // xw=sum(Nu*x*w)
            wders(0)    += ders [ 0 ](0, k) * w;                               // w=sum(Nu*w)

            Aders [ 0 ](1) += ders [ 0 ](1, k) * vertexCoordsPtr->at(1) * w;   // dxw/du=sum(dNu/du*x*w)
            wders(1)    += ders [ 0 ](1, k) * w;                               // dw/du=sum(dNu/du*w)
        }

        weight = wders(0);

        // calculation of jacobian matrix according to Eq 4.7
        jacobian(0, 0) = ( Aders [ 0 ](1) - wders(1) * Aders [ 0 ](0) / weight ) / weight; // dx/du

        Jacob = jacobian.giveDeterminant();

        //calculation of derivatives of NURBS basis functions with respect to local parameters is not covered by NURBS book
        product = Jacob * weight * weight;
        cnt = 0;
        ind = uind + 1;
        for ( k = 0; k <= degree [ 0 ]; k++ ) {
            w = cellgeo.giveVertexCoordinates(ind + k)->at(2);
            // [dNu/du*w*sum(Nu*w) - Nu*w*sum(dNu/du*w)] / [J*sum(Nu*w)^2]
            answer(cnt, 0) = ders [ 0 ](1, k) * w * weight - ders [ 0 ](0, k) * w * wders(1) / product;
            cnt++;
        }
    } else if ( nsd == 2 ) {
        FloatArray tmp1(nsd + 1), tmp2(nsd + 1);    // allow for weight

        // calculate values and derivatives of nonrational Bspline surface with weights at first (Aders, wders)
        uind = span(0) - degree [ 0 ];
        vind = span(1) - degree [ 1 ];
        ind = vind * numberOfControlPoints [ 0 ] + uind + 1;
        for ( l = 0; l <= degree [ 1 ]; l++ ) {
            tmp1.zero();
            tmp2.zero();
            for ( k = 0; k <= degree [ 0 ]; k++ ) {
                vertexCoordsPtr = cellgeo.giveVertexCoordinates(ind + k);
                w = vertexCoordsPtr->at(3);

                tmp1(0) += ders [ 0 ](0, k) * vertexCoordsPtr->at(1) * w; // sum(Nu*x*w)
                tmp1(1) += ders [ 0 ](0, k) * vertexCoordsPtr->at(2) * w; // sum(Nu*y*w)
                tmp1(2) += ders [ 0 ](0, k) * w;           // sum(Nu*w)

                tmp2(0) += ders [ 0 ](1, k) * vertexCoordsPtr->at(1) * w; // sum(dNu/du*x*w)
                tmp2(1) += ders [ 0 ](1, k) * vertexCoordsPtr->at(2) * w; // sum(dNu/du*y*w)
                tmp2(2) += ders [ 0 ](1, k) * w;           // sum(dNu/du*w)
            }

            ind += numberOfControlPoints [ 0 ];

            Aders [ 0 ](0) += ders [ 1 ](0, l) * tmp1(0); // xw=sum(Nv*sum(Nu*x*w)
            Aders [ 1 ](0) += ders [ 1 ](0, l) * tmp1(1); // yw=sum(Nv*sum(Nu*y*w)
            wders(0)    += ders [ 1 ](0, l) * tmp1(2); // w=sum(Nv*sum(Nu*w)

            Aders [ 0 ](1) += ders [ 1 ](0, l) * tmp2(0); // dxw/du=sum(Nv*sum(dNu/du*x*w)
            Aders [ 1 ](1) += ders [ 1 ](0, l) * tmp2(1); // dyw/du=sum(Nv*sum(dNu/du*y*w)
            wders(1)    += ders [ 1 ](0, l) * tmp2(2);        // dw/du=sum(Nv*sum(dNu/du*w)

            Aders [ 0 ](2) += ders [ 1 ](1, l) * tmp1(0); // dxw/dv=sum(dNv/dv*sum(Nu*x*w)
            Aders [ 1 ](2) += ders [ 1 ](1, l) * tmp1(1); // dyw/dv=sum(dNv/dv*sum(Nu*y*w)
            wders(2)    += ders [ 1 ](1, l) * tmp1(2); // dw/dv=sum(dNv/dv*sum(Nu*w)
        }

        weight = wders(0);

        // calculation of jacobian matrix according to Eq 4.19
        tmp1(0) = Aders [ 0 ](0) / weight;
        tmp1(1) = Aders [ 1 ](0) / weight;
        jacobian(0, 0) = ( Aders [ 0 ](1) - wders(1) * tmp1(0) ) / weight; // dx/du
        jacobian(0, 1) = ( Aders [ 1 ](1) - wders(1) * tmp1(1) ) / weight; // dy/du
        jacobian(1, 0) = ( Aders [ 0 ](2) - wders(2) * tmp1(0) ) / weight; // dx/dv
        jacobian(1, 1) = ( Aders [ 1 ](2) - wders(2) * tmp1(1) ) / weight; // dy/dv

        Jacob = jacobian.giveDeterminant();

        //calculation of derivatives of NURBS basis functions with respect to local parameters is not covered by NURBS book
        product = Jacob * weight * weight;
        cnt = 0;
        ind = vind * numberOfControlPoints [ 0 ] + uind + 1;
        for ( l = 0; l <= degree [ 1 ]; l++ ) {
            for ( k = 0; k <= degree [ 0 ]; k++ ) {
                w = cellgeo.giveVertexCoordinates(ind + k)->at(3);
                // dNu/du*Nv*w*sum(Nu*Nv*w) - Nu*Nv*w*sum(dNu/du*Nv*w)
                tmp1(0) = ders [ 0 ](1, k) * ders [ 1 ](0, l) * w * weight - ders [ 0 ](0, k) * ders [ 1 ](0, l) * w * wders(1);
                // Nu*dNv/dv*w*sum(Nu*Nv*w) - Nu*Nv*w*sum(Nu*dNv/dv*w)
                tmp1(1) = ders [ 0 ](0, k) * ders [ 1 ](1, l) * w * weight - ders [ 0 ](0, k) * ders [ 1 ](0, l) * w * wders(2);

                answer(cnt, 0) = ( +jacobian(1, 1) * tmp1(0) - jacobian(0, 1) * tmp1(1) ) / product;
                answer(cnt, 1) = ( -jacobian(1, 0) * tmp1(0) + jacobian(0, 0) * tmp1(1) ) / product;
                cnt++;
            }

            ind += numberOfControlPoints [ 0 ];
        }
    } else if ( nsd == 3 ) {
        FloatArray tmp1(nsd + 1), tmp2(nsd + 1);    // allow for weight
        FloatArray temp1(nsd + 1), temp2(nsd + 1), temp3(nsd + 1); // allow for weight

        // calculate values and derivatives of nonrational Bspline solid with weights at first (Aders, wders)
        uind = span(0) - degree [ 0 ];
        vind = span(1) - degree [ 1 ];
        tind = span(2) - degree [ 2 ];
        ind = tind * numberOfControlPoints [ 0 ] * numberOfControlPoints [ 1 ] + vind * numberOfControlPoints [ 0 ] + uind + 1;
        for ( m = 0; m <= degree [ 2 ]; m++ ) {
            temp1.zero();
            temp2.zero();
            temp3.zero();
            indx = ind;
            for ( l = 0; l <= degree [ 1 ]; l++ ) {
                tmp1.zero();
                tmp2.zero();
                for ( k = 0; k <= degree [ 0 ]; k++ ) {
                    vertexCoordsPtr = cellgeo.giveVertexCoordinates(ind + k);
                    w = vertexCoordsPtr->at(4);

                    tmp1(0) += ders [ 0 ](0, k) * vertexCoordsPtr->at(1) * w;              // sum(Nu*x*w)
                    tmp1(1) += ders [ 0 ](0, k) * vertexCoordsPtr->at(2) * w;              // sum(Nu*y*w)
                    tmp1(2) += ders [ 0 ](0, k) * vertexCoordsPtr->at(3) * w;              // sum(Nu*z*w)
                    tmp1(3) += ders [ 0 ](0, k) * w;                                       // sum(Nu*w)

                    tmp2(0) += ders [ 0 ](1, k) * vertexCoordsPtr->at(1) * w;              // sum(dNu/du*x*w)
                    tmp2(1) += ders [ 0 ](1, k) * vertexCoordsPtr->at(2) * w;              // sum(dNu/du*y*w)
                    tmp2(2) += ders [ 0 ](1, k) * vertexCoordsPtr->at(3) * w;              // sum(dNu/du*z*w)
                    tmp2(3) += ders [ 0 ](1, k) * w;                                       // sum(dNu/du*w)
                }

                ind += numberOfControlPoints [ 0 ];

                temp1(0) += ders [ 1 ](0, l) * tmp1(0);            // sum(Nv*sum(Nu*x*w))
                temp1(1) += ders [ 1 ](0, l) * tmp1(1);            // sum(Nv*sum(Nu*y*w))
                temp1(2) += ders [ 1 ](0, l) * tmp1(2);            // sum(Nv*sum(Nu*z*w))
                temp1(3) += ders [ 1 ](0, l) * tmp1(3);            // sum(Nv*sum(Nu*w))

                temp2(0) += ders [ 1 ](0, l) * tmp2(0);            // sum(Nv*sum(dNu/du*x*w))
                temp2(1) += ders [ 1 ](0, l) * tmp2(1);            // sum(Nv*sum(dNu/du*y*w))
                temp2(2) += ders [ 1 ](0, l) * tmp2(2);            // sum(Nv*sum(dNu/du*z*w))
                temp2(3) += ders [ 1 ](0, l) * tmp2(3);            // sum(Nv*sum(dNu/du*w))

                temp3(0) += ders [ 1 ](1, l) * tmp1(0);            // sum(dNv/dv*sum(Nu*x*w))
                temp3(1) += ders [ 1 ](1, l) * tmp1(1);            // sum(dNv/dv*sum(Nu*y*w))
                temp3(2) += ders [ 1 ](1, l) * tmp1(2);            // sum(dNv/dv*sum(Nu*z*w))
                temp3(3) += ders [ 1 ](1, l) * tmp1(3);            // sum(dNv/dv*sum(Nu*w))
            }

            ind = indx + numberOfControlPoints [ 0 ] * numberOfControlPoints [ 1 ];

            Aders [ 0 ](0) += ders [ 2 ](0, m) * temp1(0);     // x=sum(Nt*sum(Nv*sum(Nu*x*w)))
            Aders [ 1 ](0) += ders [ 2 ](0, m) * temp1(1);     // y=sum(Nt*sum(Nv*sum(Nu*y*w)))
            Aders [ 2 ](0) += ders [ 2 ](0, m) * temp1(2);     // y=sum(Nt*sum(Nv*sum(Nu*y*w)))
            wders(0)    += ders [ 2 ](0, m) * temp1(3);        // w=sum(Nt*sum(Nv*sum(Nu*w)))

            Aders [ 0 ](1) += ders [ 2 ](0, m) * temp2(0);     // dx/du=sum(Nt*sum(Nv*sum(dNu/du*x*w)))
            Aders [ 1 ](1) += ders [ 2 ](0, m) * temp2(1);     // dy/du=sum(Nt*sum(Nv*sum(dNu/du*y*w)))
            Aders [ 2 ](1) += ders [ 2 ](0, m) * temp2(2);     // dy/du=sum(Nt*sum(Nv*sum(dNu/du*y*w)))
            wders(1)    += ders [ 2 ](0, m) * temp2(3);        // dw/du=sum(Nt*sum(Nv*sum(dNu/du*w)))

            Aders [ 0 ](2) += ders [ 2 ](0, m) * temp3(0);     // dx/dv=sum(Nt*sum(dNv/dv*sum(Nu*x*w)))
            Aders [ 1 ](2) += ders [ 2 ](0, m) * temp3(1);     // dy/dv=sum(Nt*sum(dNv/dv*sum(Nu*y*w)))
            Aders [ 2 ](2) += ders [ 2 ](0, m) * temp3(2);     // dy/dv=sum(Nt*sum(dNv/dv*sum(Nu*y*w)))
            wders(2)    += ders [ 2 ](0, m) * temp3(3);        // dw/dv=sum(Nt*sum(dNv/dv*sum(Nu*w)))

            Aders [ 0 ](3) += ders [ 2 ](1, m) * temp1(0);     // dx/dt=sum(dNt/dt*sum(Nv*sum(Nu*x*w)))
            Aders [ 1 ](3) += ders [ 2 ](1, m) * temp1(1);     // dy/dt=sum(dNt/dt*sum(Nv*sum(Nu*y*w)))
            Aders [ 2 ](3) += ders [ 2 ](1, m) * temp1(2);     // dy/dt=sum(dNt/dt*sum(Nv*sum(Nu*y*w)))
            wders(3)    += ders [ 2 ](1, m) * temp1(3);        // dw/dt=sum(dNt/dt*sum(Nv*sum(Nu*w)))
        }

        weight = wders(0);

        // calculation of jacobian matrix
        tmp1(0) = Aders [ 0 ](0) / weight;
        tmp1(1) = Aders [ 1 ](0) / weight;
        tmp1(2) = Aders [ 2 ](0) / weight;
        jacobian(0, 0) = ( Aders [ 0 ](1) - wders(1) * tmp1(0) ) / weight; // dx/du
        jacobian(0, 1) = ( Aders [ 1 ](1) - wders(1) * tmp1(1) ) / weight; // dy/du
        jacobian(0, 2) = ( Aders [ 2 ](1) - wders(1) * tmp1(2) ) / weight; // dz/du
        jacobian(1, 0) = ( Aders [ 0 ](2) - wders(2) * tmp1(0) ) / weight; // dx/dv
        jacobian(1, 1) = ( Aders [ 1 ](2) - wders(2) * tmp1(1) ) / weight; // dy/dv
        jacobian(1, 2) = ( Aders [ 2 ](2) - wders(2) * tmp1(2) ) / weight; // dz/dv
        jacobian(2, 0) = ( Aders [ 0 ](3) - wders(3) * tmp1(0) ) / weight; // dx/dt
        jacobian(2, 1) = ( Aders [ 1 ](3) - wders(3) * tmp1(1) ) / weight; // dy/dt
        jacobian(2, 2) = ( Aders [ 2 ](3) - wders(3) * tmp1(2) ) / weight; // dz/dt

        Jacob = jacobian.giveDeterminant();

        //calculation of derivatives of NURBS basis functions with respect to local parameters is not covered by NURBS book
        product = Jacob * weight * weight;
        cnt = 0;
        ind = tind * numberOfControlPoints [ 0 ] * numberOfControlPoints [ 1 ] + vind * numberOfControlPoints [ 0 ] + uind + 1;
        for ( m = 0; m <= degree [ 2 ]; m++ ) {
            indx = ind;
            for ( l = 0; l <= degree [ 1 ]; l++ ) {
                for ( k = 0; k <= degree [ 0 ]; k++ ) {
                    w = cellgeo.giveVertexCoordinates(ind + k)->at(4);
                    // dNu/du*Nv*Nt*w*sum(Nu*Nv*Nt*w) - Nu*Nv*Nt*w*sum(dNu/du*Nv*Nt*w)
                    tmp1(0) = ders [ 0 ](1, k) * ders [ 1 ](0, l) * ders [ 2 ](0, m) * w * weight - ders [ 0 ](0, k) * ders [ 1 ](0, l) * ders [ 2 ](0, m) * w * wders(1);
                    // Nu*dNv/dv*Nt*w*sum(Nu*Nv*Nt*w) - Nu*Nv*Nt*w*sum(Nu*dNv/dv*Nt*w)
                    tmp1(1) = ders [ 0 ](0, k) * ders [ 1 ](1, l) * ders [ 2 ](0, m) * w * weight - ders [ 0 ](0, k) * ders [ 1 ](0, l) * ders [ 2 ](0, m) * w * wders(2);
                    // Nu*Nv*dNt/dt*w*sum(Nu*Nv*Nt*w) - Nu*Nv*Nt*w*sum(Nu*Nv*dNt/dt*w)
                    tmp1(2) = ders [ 0 ](0, k) * ders [ 1 ](0, l) * ders [ 2 ](1, m) * w * weight - ders [ 0 ](0, k) * ders [ 1 ](0, l) * ders [ 2 ](0, m) * w * wders(3);

                    answer(cnt, 0) = ( ( jacobian(1, 1) * jacobian(2, 2) - jacobian(1, 2) * jacobian(2, 1) ) * tmp1(0) +
                                      ( jacobian(0, 2) * jacobian(2, 1) - jacobian(0, 1) * jacobian(2, 2) ) * tmp1(1) +
                                      ( jacobian(0, 1) * jacobian(1, 2) - jacobian(0, 2) * jacobian(1, 1) ) * tmp1(2) ) / product;                                                      // dN/dx
                    answer(cnt, 1) = ( ( jacobian(1, 2) * jacobian(2, 0) - jacobian(1, 0) * jacobian(2, 2) ) * tmp1(0) +
                                      ( jacobian(0, 0) * jacobian(2, 2) - jacobian(0, 2) * jacobian(2, 0) ) * tmp1(1) +
                                      ( jacobian(0, 2) * jacobian(1, 0) - jacobian(0, 0) * jacobian(1, 2) ) * tmp1(2) ) / product;                                                      // dN/dy
                    answer(cnt, 2) = ( ( jacobian(1, 0) * jacobian(2, 1) - jacobian(1, 1) * jacobian(2, 0) ) * tmp1(0) +
                                      ( jacobian(0, 1) * jacobian(2, 0) - jacobian(0, 0) * jacobian(2, 1) ) * tmp1(1) +
                                      ( jacobian(0, 0) * jacobian(1, 1) - jacobian(0, 1) * jacobian(1, 0) ) * tmp1(2) ) / product;                                                      // dN/dz
                    cnt++;
                }

                ind += numberOfControlPoints [ 0 ];
            }

            ind = indx + numberOfControlPoints [ 0 ] * numberOfControlPoints [ 1 ];
        }
    } else {
        OOFEM_ERROR2("evaldNdx not implemented for nsd = %d", nsd);
    }

 #ifndef HAVE_VARIABLE_ARRAY_SIZE
    delete [] Aders;
 #endif
#endif

#ifndef HAVE_VARIABLE_ARRAY_SIZE
    delete [] ders;
#endif
}


void NURBSInterpolation :: local2global(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    /* Based on SurfacePoint A4.3 implementation*/
    FEIIGAElementGeometryWrapper *gw = ( FEIIGAElementGeometryWrapper * ) & cellgeo;
    const FloatArray *vertexCoordsPtr;
    IntArray span(nsd);
    double w, weight = 0.0;
    int i, l, k, m, ind, indx, uind, vind, tind;
#ifdef HAVE_VARIABLE_ARRAY_SIZE
    FloatArray N [ nsd ];
#else
    FloatArray *N = new FloatArray [ nsd ];
#endif

    if ( gw->knotSpan ) {
        span = * gw->knotSpan;
    } else {
        for ( i = 0; i < nsd; i++ ) {
            span(i) = this->findSpan(numberOfControlPoints [ i ], degree [ i ], lcoords(i), knotVector [ i ]);
        }
    }

    for ( i = 0; i < nsd; i++ ) {
        this->basisFuns(N [ i ], span(i), lcoords(i), degree [ i ], knotVector [ i ]);
    }

    answer.resize(nsd);
    answer.zero();

    if ( nsd == 1 ) {
        uind = span(0) - degree [ 0 ];
        ind = uind + 1;
        for ( k = 0; k <= degree [ 0 ]; k++ ) {
            vertexCoordsPtr = cellgeo.giveVertexCoordinates(ind + k);
            w = vertexCoordsPtr->at(2);
            answer(0) += N [ 0 ](k) * vertexCoordsPtr->at(1) * w;       // xw=sum(Nu*x*w)
            weight    += N [ 0 ](k) * w;                                // w=sum(Nu*w)
        }
    } else if ( nsd == 2 ) {
        FloatArray tmp(nsd + 1);       // allow for weight

        uind = span(0) - degree [ 0 ];
        vind = span(1) - degree [ 1 ];
        ind = vind * numberOfControlPoints [ 0 ] + uind + 1;
        for ( l = 0; l <= degree [ 1 ]; l++ ) {
            tmp.zero();
            for ( k = 0; k <= degree [ 0 ]; k++ ) {
                vertexCoordsPtr = cellgeo.giveVertexCoordinates(ind + k);
                w = vertexCoordsPtr->at(3);
                tmp(0) += N [ 0 ](k) * vertexCoordsPtr->at(1) * w; // sum(Nu*x*w)
                tmp(1) += N [ 0 ](k) * vertexCoordsPtr->at(2) * w; // sum(Nu*y*w)
                tmp(2) += N [ 0 ](k) * w;            // sum(Nu*w)
            }

            ind += numberOfControlPoints [ 0 ];
            answer(0) += N [ 1 ](l) * tmp(0); // xw=sum(Nv*Nu*x*w)
            answer(1) += N [ 1 ](l) * tmp(1); // yw=sum(Nv*Nu*y*w)
            weight    += N [ 1 ](l) * tmp(2); // w=sum(Nv*Nu*w)
        }
    } else if ( nsd == 3 ) {
        FloatArray tmp(nsd + 1), temp(nsd + 1);    // allow for weight

        uind = span(0) - degree [ 0 ];
        vind = span(1) - degree [ 1 ];
        tind = span(2) - degree [ 2 ];
        ind = tind * numberOfControlPoints [ 0 ] * numberOfControlPoints [ 1 ] + vind * numberOfControlPoints [ 0 ] + uind + 1;
        for ( m = 0; m <= degree [ 2 ]; m++ ) {
            temp.zero();
            indx = ind;
            for ( l = 0; l <= degree [ 1 ]; l++ ) {
                tmp.zero();
                for ( k = 0; k <= degree [ 0 ]; k++ ) {
                    vertexCoordsPtr = cellgeo.giveVertexCoordinates(ind + k);
                    w = vertexCoordsPtr->at(4);
                    tmp(0) += N [ 0 ](k) * vertexCoordsPtr->at(1) * w;               // sum(Nu*x*w)
                    tmp(1) += N [ 0 ](k) * vertexCoordsPtr->at(2) * w;               // sum(Nu*y*w)
                    tmp(2) += N [ 0 ](k) * vertexCoordsPtr->at(3) * w;               // sum(Nu*z*w)
                    tmp(3) += N [ 0 ](k) * w;                                        // sum(Nu*w)
                }

                ind += numberOfControlPoints [ 0 ];

                temp(0) += N [ 1 ](l) * tmp(0);             // sum(Nv*Nu*x*w)
                temp(1) += N [ 1 ](l) * tmp(1);             // sum(Nv*Nu*y*w)
                temp(2) += N [ 1 ](l) * tmp(2);             // sum(Nv*Nu*z*w)
                temp(3) += N [ 1 ](l) * tmp(3);             // sum(Nv*Nu*w)
            }

            ind = indx + numberOfControlPoints [ 0 ] * numberOfControlPoints [ 1 ];

            answer(0) += N [ 2 ](m) * temp(0); // xw=sum(Nv*Nu*Nt*x*w)
            answer(1) += N [ 2 ](m) * temp(1); // yw=sum(Nv*Nu*Nt*y*w)
            answer(2) += N [ 2 ](m) * temp(2); // zw=sum(Nv*Nu*Nt*z*w)
            weight    += N [ 2 ](m) * temp(3); // w=sum(Nv*Nu*Nt*w)
        }
    } else {
        OOFEM_ERROR2("local2global not implemented for nsd = %d", nsd);
    }

    answer.times(1.0 / weight);

#ifndef HAVE_VARIABLE_ARRAY_SIZE
    delete [] N;
#endif
}


double NURBSInterpolation :: giveTransformationJacobian(const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    //
    // Based on Algorithm A4.4 (p. 137) for d=1
    //
    FEIIGAElementGeometryWrapper *gw = ( FEIIGAElementGeometryWrapper * ) & cellgeo;
    const FloatArray *vertexCoordsPtr;
    FloatMatrix jacobian(nsd, nsd);
    IntArray span(nsd);
    double Jacob, w, weight;
    int i, l, k, m, ind, indx, uind, vind, tind;
#ifdef HAVE_VARIABLE_ARRAY_SIZE
    FloatMatrix ders [ nsd ];
#else
    FloatMatrix *ders = new FloatMatrix [ nsd ];
#endif

    if ( gw->knotSpan ) {
        span = * gw->knotSpan;
    } else {
        for ( i = 0; i < nsd; i++ ) {
            span(i) = this->findSpan(numberOfControlPoints [ i ], degree [ i ], lcoords(i), knotVector [ i ]);
        }
    }

    for ( i = 0; i < nsd; i++ ) {
        this->dersBasisFuns(1, lcoords(i), span(i), degree [ i ], knotVector [ i ], ders [ i ]);
    }

#if 0                       // code according NURBS book (too general allowing higher derivatives)
    if ( nsd == 2 ) {
        FloatArray tmp1(nsd + 1), tmp2(nsd + 1);    // allow for weight

        FloatMatrix Aders [ 2 ];      // derivatives in each coordinate direction on BSpline
 #ifndef OPTIMIZED_VERSION_A4dot4
        FloatMatrix Sders [ 2 ];      // derivatives in each coordinate direction on NURBS
 #endif
        FloatMatrix wders;              // derivatives in w direction on BSpline
        /*
         * IntArray Bin(2,2);      // binomial coefficients from 0 to d=1
         *          // Bin(n,k)=(n above k)=n!/k!(n-k)! for n>=k>=0
         *                                                                                              // lower triangle corresponds to Pascal triangle
         *                                                                                              // according to A4.4 it seems that only coefficients in lower triangle except the first column are used
         */
        // resizing to (2,2) has nothing common with nsd
        // it is related to the fact that 0th and 1st derivatives are computed
        for ( i = 0; i < nsd; i++ ) {
            Aders [ i ].resize(2, 2);
            Aders [ i ].zero();
 #ifndef OPTIMIZED_VERSION_A4dot4
            Sders [ i ].resize(2, 2);
 #endif
        }

        wders.resize(2, 2);
        wders.zero();

        // calculation of jacobian matrix according to A4.4
        // calculate values and derivatives of nonrational Bspline surface with weights at first (Aders, wders)
        uind = span(0) - degree [ 0 ];
        vind = span(1) - degree [ 1 ];
        ind = vind * numberOfControlPoints [ 0 ] + uind + 1;
        for ( l = 0; l <= degree [ 1 ]; l++ ) {
            tmp1.zero();
            tmp2.zero();
            for ( k = 0; k <= degree [ 0 ]; k++ ) {
                vertexCoordsPtr = cellgeo.giveVertexCoordinates(ind + k);
                w = vertexCoordsPtr->at(3);

                tmp1(0) += ders [ 0 ](0, k) * vertexCoordsPtr->at(1) * w; // sum(Nu*x*w)
                tmp1(1) += ders [ 0 ](0, k) * vertexCoordsPtr->at(2) * w; // sum(Nu*y*w)
                tmp1(2) += ders [ 0 ](0, k) * w;           // sum(Nu*w)

                tmp2(0) += ders [ 0 ](1, k) * vertexCoordsPtr->at(1) * w; // sum(dNu/du*x*w)
                tmp2(1) += ders [ 0 ](1, k) * vertexCoordsPtr->at(2) * w; // sum(dNu/du*y*w)
                tmp2(2) += ders [ 0 ](1, k) * w;           // sum(dNu/du*w)
            }

            ind += numberOfControlPoints [ 0 ];

            Aders [ 0 ](0, 0) += ders [ 1 ](0, l) * tmp1(0); // xw=sum(Nv*sum(Nu*x*w))
            Aders [ 1 ](0, 0) += ders [ 1 ](0, l) * tmp1(1); // yw=sum(Nv*sum(Nu*y*w))
            wders(0, 0)    += ders [ 1 ](0, l) * tmp1(2); // w=sum(Nv*sum(Nu*w))

            Aders [ 0 ](0, 1) += ders [ 1 ](1, l) * tmp1(0); // dxw/dv=sum(dNv/dv*sum(Nu*x*w))
            Aders [ 1 ](0, 1) += ders [ 1 ](1, l) * tmp1(1); // dyw/dv=sum(dNv/dv*sum(Nu*y*w))
            wders(0, 1)    += ders [ 1 ](1, l) * tmp1(2); // dw/dv=sum(dNv/dv*sum(Nu*w))

            Aders [ 0 ](1, 0) += ders [ 1 ](0, l) * tmp2(0); // dxw/du=sum(Nv*sum(dNu/du*x*w))
            Aders [ 1 ](1, 0) += ders [ 1 ](0, l) * tmp2(1); // dyw/du=sum(Nv*sum(dNu/du*y*w))
            wders(1, 0)    += ders [ 1 ](0, l) * tmp2(2);       // dw/du=sum(Nv*sum(dNu/du*w))
        }

        weight = wders(0, 0);

 #ifndef OPTIMIZED_VERSION_A4dot4
        int j;
        const int d = 1;
        // calculate values and derivatives of NURBS surface (A4.4)
        // since all entries in Pascal triangle up to d=1 are 1, binomial coefficients are ignored
        for ( k = 0; k <= d; k++ ) {
            for ( l = 0; l <= d - k; l++ ) {
                tmp1(0) = Aders [ 0 ](k, l);
                tmp1(1) = Aders [ 1 ](k, l);
                for ( j = 1; j <= l; j++ ) {
                    tmp1(0) -= wders(0, j) * Sders [ 0 ](k, l - j);            // *Bin(l,j)
                    tmp1(1) -= wders(0, j) * Sders [ 1 ](k, l - j);            // *Bin(l,j)
                }

                for ( i = 1; i <= k; i++ ) {
                    tmp1(0) -= wders(i, 0) * Sders [ 0 ](k - i, l);            // *Bin(k,i)
                    tmp1(1) -= wders(i, 0) * Sders [ 1 ](k - i, l);            // *Bin(k,i)
                    tmp2.zero();
                    for ( j = 1; j <= l; j++ ) {
                        tmp2(0) += wders(i, j) * Sders [ 0 ](k - i, l - j);              // *Bin(l,j)
                        tmp2(1) += wders(i, j) * Sders [ 1 ](k - i, l - j);              // *Bin(l,j)
                    }

                    tmp1(0) -= tmp2(0);                     // *Bin(k,i)
                    tmp1(1) -= tmp2(1);                     // *Bin(k,i)
                }

                Sders [ 0 ](k, l) = tmp1(0) / weight;
                Sders [ 1 ](k, l) = tmp1(1) / weight;
            }
        }

        jacobian(0, 0) = Sders [ 0 ](1, 0);      // dx/du
        jacobian(0, 1) = Sders [ 1 ](1, 0);      // dy/du
        jacobian(1, 0) = Sders [ 0 ](0, 1);      // dx/dv
        jacobian(1, 1) = Sders [ 1 ](0, 1);      // dy/dv
 #else
        // optimized version of A4.4 for d=1, binomial coefficients ignored

        /*
         * // k=0 l=0 loop
         * Sders[0](0,0) = Aders[0](0,0) / weight;
         * Sders[1](0,0) = Aders[1](0,0) / weight;
         * // k=1 l=0 loop
         * Sders[0](1,0) = (Aders[0](1,0)-wders(1,0)*Sders[0](0,0)) / weight;
         * Sders[1](1,0) = (Aders[1](1,0)-wders(1,0)*Sders[1](0,0)) / weight;
         * // k=0 l=1 loop
         * Sders[0](0,1) = (Aders[0](0,1)-wders(0,1)*Sders[0](0,0)) / weight;
         * Sders[1](0,1) = (Aders[1](0,1)-wders(0,1)*Sders[1](0,0)) / weight;
         *
         * jacobian(0,0) = Sders[0](1,0);   // dx/du
         * jacobian(0,1) = Sders[1](1,0);   // dy/du
         * jacobian(1,0) = Sders[0](0,1);   // dx/dv
         * jacobian(1,1) = Sders[1](0,1);   // dy/dv
         */

        // k=0 l=0 loop
        tmp1(0) = Aders [ 0 ](0, 0) / weight;
        tmp1(1) = Aders [ 1 ](0, 0) / weight;
        // k=1 l=0 loop
        jacobian(0, 0) = ( Aders [ 0 ](1, 0) - wders(1, 0) * tmp1(0) ) / weight; // dx/du
        jacobian(0, 1) = ( Aders [ 1 ](1, 0) - wders(1, 0) * tmp1(1) ) / weight; // dy/du
        // k=0 l=1 loop
        jacobian(1, 0) = ( Aders [ 0 ](0, 1) - wders(0, 1) * tmp1(0) ) / weight; // dx/dv
        jacobian(1, 1) = ( Aders [ 1 ](0, 1) - wders(0, 1) * tmp1(1) ) / weight; // dy/dv
 #endif
    }     else {
        OOFEM_ERROR2("giveTransformationJacobian not implemented for nsd = %d", nsd);
    }

#else
 #ifdef HAVE_VARIABLE_ARRAY_SIZE
    FloatArray Aders [ nsd ];  // 0th and 1st derivatives in each coordinate direction on BSpline
 #else
    FloatArray *Aders = new FloatArray [ nsd ];
 #endif
    FloatArray wders;          // 0th and 1st derivatives in w direction on BSpline

    for ( i = 0; i < nsd; i++ ) {
        Aders [ i ].resize(nsd + 1);
        Aders [ i ].zero();
    }

    wders.resize(nsd + 1);
    wders.zero();

    if ( nsd == 1 ) {
        // calculate values and derivatives of nonrational Bspline curve with weights at first (Aders, wders)
        uind = span(0) - degree [ 0 ];
        ind = uind + 1;
        for ( k = 0; k <= degree [ 0 ]; k++ ) {
            vertexCoordsPtr = cellgeo.giveVertexCoordinates(ind + k);
            w = vertexCoordsPtr->at(2);

            Aders [ 0 ](0) += ders [ 0 ](0, k) * vertexCoordsPtr->at(1) * w;   // xw=sum(Nu*x*w)
            wders(0)    += ders [ 0 ](0, k) * w;                               // w=sum(Nu*w)

            Aders [ 0 ](1) += ders [ 0 ](1, k) * vertexCoordsPtr->at(1) * w;   // dxw/du=sum(dNu/du*x*w)
            wders(1)    += ders [ 0 ](1, k) * w;                               // dw/du=sum(dNu/du*w)
        }

        weight = wders(0);

        // calculation of jacobian matrix according to Eq 4.7
        jacobian(0, 0) = ( Aders [ 0 ](1) - wders(1) * Aders [ 0 ](0) / weight ) / weight; // dx/du
    } else if ( nsd == 2 ) {
        FloatArray tmp1(nsd + 1), tmp2(nsd + 1);    // allow for weight

        // calculate values and derivatives of nonrational Bspline surface with weights at first (Aders, wders)
        uind = span(0) - degree [ 0 ];
        vind = span(1) - degree [ 1 ];
        ind = vind * numberOfControlPoints [ 0 ] + uind + 1;
        for ( l = 0; l <= degree [ 1 ]; l++ ) {
            tmp1.zero();
            tmp2.zero();
            for ( k = 0; k <= degree [ 0 ]; k++ ) {
                vertexCoordsPtr = cellgeo.giveVertexCoordinates(ind + k);
                w = vertexCoordsPtr->at(3);

                tmp1(0) += ders [ 0 ](0, k) * vertexCoordsPtr->at(1) * w; // sum(Nu*x*w)
                tmp1(1) += ders [ 0 ](0, k) * vertexCoordsPtr->at(2) * w; // sum(Nu*y*w)
                tmp1(2) += ders [ 0 ](0, k) * w;           // sum(Nu*w)

                tmp2(0) += ders [ 0 ](1, k) * vertexCoordsPtr->at(1) * w; // sum(dNu/du*x*w)
                tmp2(1) += ders [ 0 ](1, k) * vertexCoordsPtr->at(2) * w; // sum(dNu/du*y*w)
                tmp2(2) += ders [ 0 ](1, k) * w;           // sum(dNu/du*w)
            }

            ind += numberOfControlPoints [ 0 ];

            Aders [ 0 ](0) += ders [ 1 ](0, l) * tmp1(0); // xw=sum(Nv*sum(Nu*x*w)
            Aders [ 1 ](0) += ders [ 1 ](0, l) * tmp1(1); // yw=sum(Nv*sum(Nu*y*w)
            wders(0)    += ders [ 1 ](0, l) * tmp1(2); // w=sum(Nv*sum(Nu*w)

            Aders [ 0 ](1) += ders [ 1 ](0, l) * tmp2(0); // dxw/du=sum(Nv*sum(dNu/du*x*w)
            Aders [ 1 ](1) += ders [ 1 ](0, l) * tmp2(1); // dyw/du=sum(Nv*sum(dNu/du*y*w)
            wders(1)    += ders [ 1 ](0, l) * tmp2(2);        // dw/du=sum(Nv*sum(dNu/du*w)

            Aders [ 0 ](2) += ders [ 1 ](1, l) * tmp1(0); // dxw/dv=sum(dNv/dv*sum(Nu*x*w)
            Aders [ 1 ](2) += ders [ 1 ](1, l) * tmp1(1); // dyw/dv=sum(dNv/dv*sum(Nu*y*w)
            wders(2)    += ders [ 1 ](1, l) * tmp1(2); // dw/dv=sum(dNv/dv*sum(Nu*w)
        }

        weight = wders(0);

        // calculation of jacobian matrix according to Eq 4.19
        tmp1(0) = Aders [ 0 ](0) / weight;
        tmp1(1) = Aders [ 1 ](0) / weight;
        jacobian(0, 0) = ( Aders [ 0 ](1) - wders(1) * tmp1(0) ) / weight; // dx/du
        jacobian(0, 1) = ( Aders [ 1 ](1) - wders(1) * tmp1(1) ) / weight; // dy/du
        jacobian(1, 0) = ( Aders [ 0 ](2) - wders(2) * tmp1(0) ) / weight; // dx/dv
        jacobian(1, 1) = ( Aders [ 1 ](2) - wders(2) * tmp1(1) ) / weight; // dy/dv
    } else if ( nsd == 3 ) {
        FloatArray tmp1(nsd + 1), tmp2(nsd + 1);    // allow for weight
        FloatArray temp1(nsd + 1), temp2(nsd + 1), temp3(nsd + 1); // allow for weight

        // calculate values and derivatives of nonrational Bspline solid with weights at first (Aders, wders)
        uind = span(0) - degree [ 0 ];
        vind = span(1) - degree [ 1 ];
        tind = span(2) - degree [ 2 ];
        ind = tind * numberOfControlPoints [ 0 ] * numberOfControlPoints [ 1 ] + vind * numberOfControlPoints [ 0 ] + uind + 1;
        for ( m = 0; m <= degree [ 2 ]; m++ ) {
            temp1.zero();
            temp2.zero();
            temp3.zero();
            indx = ind;
            for ( l = 0; l <= degree [ 1 ]; l++ ) {
                tmp1.zero();
                tmp2.zero();
                for ( k = 0; k <= degree [ 0 ]; k++ ) {
                    vertexCoordsPtr = cellgeo.giveVertexCoordinates(ind + k);
                    w = vertexCoordsPtr->at(4);

                    tmp1(0) += ders [ 0 ](0, k) * vertexCoordsPtr->at(1) * w;              // sum(Nu*x*w)
                    tmp1(1) += ders [ 0 ](0, k) * vertexCoordsPtr->at(2) * w;              // sum(Nu*y*w)
                    tmp1(2) += ders [ 0 ](0, k) * vertexCoordsPtr->at(3) * w;              // sum(Nu*z*w)
                    tmp1(3) += ders [ 0 ](0, k) * w;                                       // sum(Nu*w)

                    tmp2(0) += ders [ 0 ](1, k) * vertexCoordsPtr->at(1) * w;              // sum(dNu/du*x*w)
                    tmp2(1) += ders [ 0 ](1, k) * vertexCoordsPtr->at(2) * w;              // sum(dNu/du*y*w)
                    tmp2(2) += ders [ 0 ](1, k) * vertexCoordsPtr->at(3) * w;              // sum(dNu/du*z*w)
                    tmp2(3) += ders [ 0 ](1, k) * w;                                       // sum(dNu/du*w)
                }

                ind += numberOfControlPoints [ 0 ];

                temp1(0) += ders [ 1 ](0, l) * tmp1(0);            // sum(Nv*sum(Nu*x*w))
                temp1(1) += ders [ 1 ](0, l) * tmp1(1);            // sum(Nv*sum(Nu*y*w))
                temp1(2) += ders [ 1 ](0, l) * tmp1(2);            // sum(Nv*sum(Nu*z*w))
                temp1(3) += ders [ 1 ](0, l) * tmp1(3);            // sum(Nv*sum(Nu*w))

                temp2(0) += ders [ 1 ](0, l) * tmp2(0);            // sum(Nv*sum(dNu/du*x*w))
                temp2(1) += ders [ 1 ](0, l) * tmp2(1);            // sum(Nv*sum(dNu/du*y*w))
                temp2(2) += ders [ 1 ](0, l) * tmp2(2);            // sum(Nv*sum(dNu/du*z*w))
                temp2(3) += ders [ 1 ](0, l) * tmp2(3);            // sum(Nv*sum(dNu/du*w))

                temp3(0) += ders [ 1 ](1, l) * tmp1(0);            // sum(dNv/dv*sum(Nu*x*w))
                temp3(1) += ders [ 1 ](1, l) * tmp1(1);            // sum(dNv/dv*sum(Nu*y*w))
                temp3(2) += ders [ 1 ](1, l) * tmp1(2);            // sum(dNv/dv*sum(Nu*z*w))
                temp3(3) += ders [ 1 ](1, l) * tmp1(3);            // sum(dNv/dv*sum(Nu*w))
            }

            ind = indx + numberOfControlPoints [ 0 ] * numberOfControlPoints [ 1 ];

            Aders [ 0 ](0) += ders [ 2 ](0, m) * temp1(0);     // x=sum(Nt*sum(Nv*sum(Nu*x*w)))
            Aders [ 1 ](0) += ders [ 2 ](0, m) * temp1(1);     // y=sum(Nt*sum(Nv*sum(Nu*y*w)))
            Aders [ 2 ](0) += ders [ 2 ](0, m) * temp1(2);     // y=sum(Nt*sum(Nv*sum(Nu*y*w)))
            wders(0)    += ders [ 2 ](0, m) * temp1(3);        // w=sum(Nt*sum(Nv*sum(Nu*w)))

            Aders [ 0 ](1) += ders [ 2 ](0, m) * temp2(0);     // dx/du=sum(Nt*sum(Nv*sum(dNu/du*x*w)))
            Aders [ 1 ](1) += ders [ 2 ](0, m) * temp2(1);     // dy/du=sum(Nt*sum(Nv*sum(dNu/du*y*w)))
            Aders [ 2 ](1) += ders [ 2 ](0, m) * temp2(2);     // dy/du=sum(Nt*sum(Nv*sum(dNu/du*y*w)))
            wders(1)    += ders [ 2 ](0, m) * temp2(3);        // dw/du=sum(Nt*sum(Nv*sum(dNu/du*w)))

            Aders [ 0 ](2) += ders [ 2 ](0, m) * temp3(0);     // dx/dv=sum(Nt*sum(dNv/dv*sum(Nu*x*w)))
            Aders [ 1 ](2) += ders [ 2 ](0, m) * temp3(1);     // dy/dv=sum(Nt*sum(dNv/dv*sum(Nu*y*w)))
            Aders [ 2 ](2) += ders [ 2 ](0, m) * temp3(2);     // dy/dv=sum(Nt*sum(dNv/dv*sum(Nu*y*w)))
            wders(2)    += ders [ 2 ](0, m) * temp3(3);        // dw/dv=sum(Nt*sum(dNv/dv*sum(Nu*w)))

            Aders [ 0 ](3) += ders [ 2 ](1, m) * temp1(0);     // dx/dt=sum(dNt/dt*sum(Nv*sum(Nu*x*w)))
            Aders [ 1 ](3) += ders [ 2 ](1, m) * temp1(1);     // dy/dt=sum(dNt/dt*sum(Nv*sum(Nu*y*w)))
            Aders [ 2 ](3) += ders [ 2 ](1, m) * temp1(2);     // dy/dt=sum(dNt/dt*sum(Nv*sum(Nu*y*w)))
            wders(3)    += ders [ 2 ](1, m) * temp1(3);        // dw/dt=sum(dNt/dt*sum(Nv*sum(Nu*w)))
        }

        weight = wders(0);

        // calculation of jacobian matrix
        tmp1(0) = Aders [ 0 ](0) / weight;
        tmp1(1) = Aders [ 1 ](0) / weight;
        tmp1(2) = Aders [ 2 ](0) / weight;
        jacobian(0, 0) = ( Aders [ 0 ](1) - wders(1) * tmp1(0) ) / weight; // dx/du
        jacobian(0, 1) = ( Aders [ 1 ](1) - wders(1) * tmp1(1) ) / weight; // dy/du
        jacobian(0, 2) = ( Aders [ 2 ](1) - wders(1) * tmp1(2) ) / weight; // dz/du
        jacobian(1, 0) = ( Aders [ 0 ](2) - wders(2) * tmp1(0) ) / weight; // dx/dv
        jacobian(1, 1) = ( Aders [ 1 ](2) - wders(2) * tmp1(1) ) / weight; // dy/dv
        jacobian(1, 2) = ( Aders [ 2 ](2) - wders(2) * tmp1(2) ) / weight; // dz/dv
        jacobian(2, 0) = ( Aders [ 0 ](3) - wders(3) * tmp1(0) ) / weight; // dx/dt
        jacobian(2, 1) = ( Aders [ 1 ](3) - wders(3) * tmp1(1) ) / weight; // dy/dt
        jacobian(2, 2) = ( Aders [ 2 ](3) - wders(3) * tmp1(2) ) / weight; // dz/dt
    } else {
        OOFEM_ERROR2("giveTransformationJacobianMatrix not implemented for nsd = %d", nsd);
    }

 #ifndef HAVE_VARIABLE_ARRAY_SIZE
    delete [] Aders;
 #endif


#endif

    Jacob = jacobian.giveDeterminant();

    if ( fabs(Jacob) < 1.0e-10 ) {
        OOFEM_ERROR("giveTransformationJacobianMatrix - zero Jacobian");
    }

#ifndef HAVE_VARIABLE_ARRAY_SIZE
    delete [] ders;
#endif


    return Jacob;
}
} // end namespace oofem
