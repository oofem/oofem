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

#include "feitspline.h"
#include "floatarray.h"
#include "floatmatrix.h"
#include "mathfem.h"
#include "iga.h"

namespace oofem {
// optimized version of A4.4 for d=1
#define OPTIMIZED_VERSION_A4dot4


TSplineInterpolation :: ~TSplineInterpolation()
{
    for ( int i = 0; i <= numberOfControlPoints [ 0 ]; i++ ) {
        for ( int j = 0; j < nsd; j++ ) {
            delete [] localIndexKnotVector [ i ] [ j ];
        }

        delete [] localIndexKnotVector [ i ];
    }

    delete [] localIndexKnotVector;

    delete [] openLocalKnotVector;
}


IRResultType TSplineInterpolation :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                 // Required by IR_GIVE_FIELD macro

    BSplineInterpolation :: initializeFrom(ir);

    IntArray localIndexKnotVector_tmp;
    int *indexKnotVec, indexKnotVal;
    int pos, p;

    InputFieldType IFT_localIndexKnotVector [ 3 ] = {
        _IFT_TSplineInterpolation_localIndexKnotVectorU,
        _IFT_TSplineInterpolation_localIndexKnotVectorV,
        _IFT_TSplineInterpolation_localIndexKnotVectorW
    };
    int max_deg = 0;
    for ( int i = 0; i < nsd; i++ ) {
        if ( degree [ i ] > max_deg ) {
            max_deg = degree [ i ];
        }
    }

    openLocalKnotVector = new double [ 3 * max_deg + 2 ];

    localIndexKnotVector = new int ** [ totalNumberOfControlPoints ];
    for ( int i = 0; i < totalNumberOfControlPoints; i++ ) {
        localIndexKnotVector [ i ] = new int * [ nsd ];
    }

    for ( int n = 0; n < nsd; n++ ) {
        localIndexKnotVector_tmp.clear();
        IR_GIVE_FIELD(ir, localIndexKnotVector_tmp, IFT_localIndexKnotVector [ n ]);
        if ( localIndexKnotVector_tmp.giveSize() != totalNumberOfControlPoints * ( degree [ n ] + 2 ) ) {
            OOFEM_WARNING("invalid size of knot vector %s", IFT_localIndexKnotVector [ n ]);
            return IRRT_BAD_FORMAT;
        }

        pos = 0;
        for ( int i = 0; i < totalNumberOfControlPoints; i++ ) {
            indexKnotVec = localIndexKnotVector [ i ] [ n ] = new int [ degree [ n ] + 2 ];

            p = 0;
            for ( int j = 0; j < degree [ n ] + 2; j++ ) {
                indexKnotVec [ p++ ] = localIndexKnotVector_tmp(pos++);
            }

            // check for monotonicity of local index knot vector with multiplicity
            indexKnotVal = indexKnotVec [ 0 ];
            for ( int j = 1; j < degree [ n ] + 2; j++ ) {
                if ( indexKnotVal > indexKnotVec [ j ] ) {
                    OOFEM_WARNING("local index knot vector %s of control point %d is not monotonic",
                                 IFT_localIndexKnotVector [ n ], i + 1);
                    return IRRT_BAD_FORMAT;
                }

                /* this is only for the case when TSpline = NURBS
                 * if(indexKnotVal+1 < indexKnotVec[j])
                 *      OOFEM_ERROR("local index knot vector %s of control point %d is not continuous",
                 *                                                       IFT_localIndexKnotVectorString[n], i+1);
                 */
                indexKnotVal = indexKnotVec [ j ];
            }

            // check for nondegeneracy of local index knot vector
            if ( indexKnotVal == indexKnotVec [ 0 ] ) {
                OOFEM_WARNING("local index knot vector %s of control point %d is degenerated",
                             IFT_localIndexKnotVector [ n ], i + 1);
                return IRRT_BAD_FORMAT;
            }

            // check for range of local index knot vector
            if ( indexKnotVec [ 0 ] <= 0 || indexKnotVal > knotValues [ n ].giveSize() ) {
                OOFEM_WARNING("local index knot vector %s of control point %d out of range",
                             IFT_localIndexKnotVector [ n ], i + 1);
                return IRRT_BAD_FORMAT;
            }
        }
    }

    return IRRT_OK;
}


void TSplineInterpolation :: evalN(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    FEIIGAElementGeometryWrapper *gw = ( FEIIGAElementGeometryWrapper * ) & cellgeo;
    FloatArray N(nsd);
    IntArray span(nsd);
    IntArray mask;
    double sum = 0.0, val;
    int count;

    if ( nsd != 2 ) {
        OOFEM_ERROR("implemented for nsd = %d", nsd);
    }

    if ( gw->knotSpan ) {
        span = * gw->knotSpan;
    } else {
        for ( int i = 0; i < nsd; i++ ) {
            span(i) = this->findSpan(numberOfControlPoints [ i ], degree [ i ], lcoords(i), knotVector [ i ]);
        }
    }

    // identify which basis functions are nonzero
    giveKnotSpanBasisFuncMask(span, mask);
    count = mask.giveSize();
    answer.resize(count);

    if ( nsd == 2 ) {
        for ( int k = 0; k < count; k++ ) {
            for ( int i = 0; i < nsd; i++ ) {
                N(i) = this->basisFunction(lcoords(i), degree [ i ], * giveKnotValues(i + 1), localIndexKnotVector [ mask(k) - 1 ] [ i ]);
            }

            answer(k) = val = N(0) * N(1) * cellgeo.giveVertexCoordinates( mask(k) )->at(3);        // Nu*Nv*w
            sum += val;
        }
    }

    while ( count ) {
        answer.at(count--) /= sum;
    }
}


double TSplineInterpolation :: evaldNdx(FloatMatrix &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    FEIIGAElementGeometryWrapper *gw = ( FEIIGAElementGeometryWrapper * ) & cellgeo;
    const FloatArray *vertexCoordsPtr;
    FloatMatrix jacobian(nsd, nsd);
    FloatArray temp(nsd);
    IntArray span(nsd);
    IntArray mask;
    double Jacob = 0., product, w, xw, yw, weight;
    int count;
    std :: vector< FloatArray > tmp_ders(nsd);
    std :: vector< FloatMatrix > ders(nsd);

    /*
     * IntArray Bin(2,2);      // binomial coefficients from 0 to d=1
     *                // Bin(n,k)=(n above k)=n!/k!(n-k)! for n>=k>=0
     *                // lower triangle corresponds to Pascal triangle
     *                                                                                              // according to A4.4 it seems that only coefficients in lower triangle except the first column are used
     */
    if ( nsd != 2 ) {
        OOFEM_ERROR("not implemented for nsd = %d", nsd);
    }

    if ( gw->knotSpan ) {
        span = * gw->knotSpan;
    } else {
        for ( int i = 0; i < nsd; i++ ) {
            span(i) = this->findSpan(numberOfControlPoints [ i ], degree [ i ], lcoords(i), knotVector [ i ]);
        }
    }

    // identify which basis functions are nonzero
    giveKnotSpanBasisFuncMask(span, mask);
    count = mask.giveSize();
    answer.resize(count, nsd);

    for ( int i = 0; i < nsd; i++ ) {
        ders [ i ].resize(2, count);
    }

    if ( nsd == 2 ) {
        FloatMatrix Aders [ 2 ];      // derivatives in each coordinate direction on BSpline
        //FloatMatrix Sders[nsd]; // derivatives in each coordinate direction on TSpline
        FloatMatrix wders;              // derivatives in w direction on BSpline

        // resizing to (2,2) has nothing common with nsd
        // it is related to the fact that 0th and 1st derivatives are computed
        for ( int i = 0; i < nsd; i++ ) {
            Aders [ i ].resize(2, 2);
            Aders [ i ].zero();
            //Sders[i].resize(2,2);
        }

        wders.resize(2, 2);
        wders.zero();

        for ( int k = 0; k < count; k++ ) {
            for ( int i = 0; i < nsd; i++ ) {
                // it would be simpler if I could pass k-th column of ders[i] directly to dersBasisFunction HUHU array
                this->dersBasisFunction(1, lcoords(i), degree [ i ], * giveKnotValues(i + 1), localIndexKnotVector [ mask(k) - 1 ] [ i ], tmp_ders [ i ]);
                ders [ i ](0, k) = tmp_ders [ i ](0);
                ders [ i ](1, k) = tmp_ders [ i ](1);
            }

            // calculation of jacobian matrix in similar fashion as A4.4
            // calculate values and derivatives of nonrational Bspline surface with weights at first (Aders, wders)
            vertexCoordsPtr = cellgeo.giveVertexCoordinates( mask(k) );
            w = vertexCoordsPtr->at(3);
            xw = vertexCoordsPtr->at(1) * w;
            yw = vertexCoordsPtr->at(2) * w;

            product = tmp_ders [ 0 ](0) * tmp_ders [ 1 ](0);       // Nu*Nv
            Aders [ 0 ](0, 0) += product * xw;             // x=sum Nu*Nv*x*w
            Aders [ 1 ](0, 0) += product * yw;             // x=sum Nu*Nv*y*w
            wders(0, 0)    += product * w;                 // w=sum Nu*Nv*w

            product = tmp_ders [ 0 ](1) * tmp_ders [ 1 ](0);       // dNu/du*Nv
            Aders [ 0 ](1, 0) += product * xw; // dx/du=sum dNu/du*Nv*x*w
            Aders [ 1 ](1, 0) += product * yw; // dy/du=sum dNu/du*Nv*y*w
            wders(1, 0)    += product * w;                 // dw/du=sum dNu/du*Nv*w

            product = tmp_ders [ 0 ](0) * tmp_ders [ 1 ](1);       // Nu*dNv/dv
            Aders [ 0 ](0, 1) += product * xw; // dx/dv=sum Nu*dNv/dv*x*w
            Aders [ 1 ](0, 1) += product * yw; // dy/dv=sum Nu*dNv/dv*y*w
            wders(0, 1)    += product * w; // dw/dv=sum Nu*dNv/dv*w
        }

        weight = wders(0, 0);

        // optimized version of A4.4 for d=1, binomial coefficients ignored
        /*
         *      Sders[0](0,0) = Aders[0](0,0) / weight;
         *      Sders[1](0,0) = Aders[1](0,0) / weight;
         *      Sders[0](0,1) = (Aders[0](0,1)-wders(0,1)*Sders[0](0,0)) / weight;
         *      Sders[1](0,1) = (Aders[1](0,1)-wders(0,1)*Sders[1](0,0)) / weight;
         *      Sders[0](1,0) = (Aders[0](1,0)-wders(1,0)*Sders[0](0,0)) / weight;
         *      Sders[1](1,0) = (Aders[1](1,0)-wders(1,0)*Sders[1](0,0)) / weight;
         *
         *      jacobian(0,0) = Sders[0](1,0);   // dx/du
         *      jacobian(0,1) = Sders[1](1,0);   // dy/du
         *      jacobian(1,0) = Sders[0](0,1);   // dx/dv
         *      jacobian(1,1) = Sders[1](0,1);   // dy/dv
         */

        temp(0) = Aders [ 0 ](0, 0) / weight;
        temp(1) = Aders [ 1 ](0, 0) / weight;
        jacobian(1, 0) = ( Aders [ 0 ](0, 1) - wders(0, 1) * temp(0) ) / weight; // dx/dv
        jacobian(1, 1) = ( Aders [ 1 ](0, 1) - wders(0, 1) * temp(1) ) / weight; // dy/dv
        jacobian(0, 0) = ( Aders [ 0 ](1, 0) - wders(1, 0) * temp(0) ) / weight; // dx/du
        jacobian(0, 1) = ( Aders [ 1 ](1, 0) - wders(1, 0) * temp(1) ) / weight; // dy/du

        Jacob = jacobian.giveDeterminant();

        //calculation of derivatives of TSpline basis functions with respect to local parameters
        product = Jacob * weight * weight;

        for ( int k = 0; k < count; k++ ) {
            w = cellgeo.giveVertexCoordinates( mask(k) )->at(3);
            // dNu/du*Nv*w*sum(Nv*Nu*w) - Nu*Nv*w*sum(dNu/du*Nv*w)
            temp(0) = ders [ 0 ](1, k) * ders [ 1 ](0, k) * w * weight - ders [ 0 ](0, k) * ders [ 1 ](0, k) * w * wders(1, 0);
            // Nu*dNv/dv*w*sum(Nv*Nu*w) - Nu*Nv*w*sum(Nu*dNv/dv*w)
            temp(1) = ders [ 0 ](0, k) * ders [ 1 ](1, k) * w * weight - ders [ 0 ](0, k) * ders [ 1 ](0, k) * w * wders(0, 1);

            answer(k, 0) = ( jacobian(1, 1) * temp(0) - jacobian(0, 1) * temp(1) ) / product;
            answer(k, 1) = ( -jacobian(1, 0) * temp(0) + jacobian(0, 0) * temp(1) ) / product;
        }
    }

    return Jacob;
}


void TSplineInterpolation :: local2global(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    /* Based on SurfacePoint A4.3 implementation*/
    FEIIGAElementGeometryWrapper *gw = ( FEIIGAElementGeometryWrapper * ) & cellgeo;
    const FloatArray *vertexCoordsPtr;
    FloatArray N(nsd);
    IntArray span(nsd);
    IntArray mask;
    double w, xw, yw, product, weight = 0.0;
    int count;

    if ( nsd != 2 ) {
        OOFEM_ERROR("not implemented for nsd = %d", nsd);
    }

    if ( gw->knotSpan ) {
        span = * gw->knotSpan;
    } else {
        for ( int i = 0; i < nsd; i++ ) {
            span(i) = this->findSpan(numberOfControlPoints [ i ], degree [ i ], lcoords(i), knotVector [ i ]);
        }
    }

    // identify which basis functions are nonzero
    giveKnotSpanBasisFuncMask(span, mask);
    count = mask.giveSize();

    answer.resize(nsd);
    answer.zero();

    if ( nsd == 2 ) {
        for ( int k = 0; k < count; k++ ) {
            for ( int i = 0; i < nsd; i++ ) {
                N(i) = this->basisFunction(lcoords(i), degree [ i ], * giveKnotValues(i + 1), localIndexKnotVector [ mask(k) - 1 ] [ i ]);
            }

            vertexCoordsPtr = cellgeo.giveVertexCoordinates( mask(k) );
            w = vertexCoordsPtr->at(3);
            xw = vertexCoordsPtr->at(1) * w;
            yw = vertexCoordsPtr->at(2) * w;

            product = N(0) * N(1);                // Nu*Nv
            answer(0) += product * xw; // x=sum Nu*Nv*x*w
            answer(1) += product * yw; // y=sum Nu*Nv*y*w
            weight    += product * w; // w=sum Nu*Nv*w
        }
    }

    answer.times(1.0 / weight);
}


void TSplineInterpolation :: giveJacobianMatrixAt(FloatMatrix &jacobian, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    //
    // Based on Algorithm A4.4 (p. 137) for d=1
    //
    FEIIGAElementGeometryWrapper *gw = ( FEIIGAElementGeometryWrapper * ) & cellgeo;
    const FloatArray *vertexCoordsPtr;
    FloatArray temp(nsd);
    IntArray span(nsd);
    IntArray mask;
    double w, xw, yw, product, weight;
    int count;
    std :: vector< FloatArray > ders(nsd);
    jacobian.resize(nsd, nsd);

    /*
     * IntArray Bin(2,2);      // binomial coefficients from 0 to d=1
     *                // Bin(n,k)=(n above k)=n!/k!(n-k)! for n>=k>=0
     *                // lower triangle corresponds to Pascal triangle
     *                                                                                              // according to A4.4 it seems that only coefficients in lower triangle except the first column are used
     */
    if ( nsd != 2 ) {
        OOFEM_ERROR("not implemented for nsd = %d", nsd);
    }

    if ( gw->knotSpan ) {
        span = * gw->knotSpan;
    } else {
        for ( int i = 0; i < nsd; i++ ) {
            span(i) = this->findSpan(numberOfControlPoints [ i ], degree [ i ], lcoords(i), knotVector [ i ]);
        }
    }

    // identify which basis functions are nonzero
    giveKnotSpanBasisFuncMask(span, mask);
    count = mask.giveSize();

    if ( nsd == 2 ) {
        FloatMatrix Aders [ 2 ];      // derivatives in each coordinate direction on BSpline
        //FloatMatrix Sders[nsd]; // derivatives in each coordinate direction on TSpline
        FloatMatrix wders;              // derivatives in w direction on BSpline

        // resizing to (2,2) has nothing common with nsd
        // it is related to the fact that 0th and 1st derivatives are computed
        for ( int i = 0; i < nsd; i++ ) {
            Aders [ i ].resize(2, 2);
            Aders [ i ].zero();
            //Sders[i].resize(2,2);
        }

        wders.resize(2, 2);
        wders.zero();

        for ( int k = 0; k < count; k++ ) {
            for ( int i = 0; i < nsd; i++ ) {
                this->dersBasisFunction(1, lcoords(i), degree [ i ], * giveKnotValues(i + 1), localIndexKnotVector [ mask(k) - 1 ] [ i ], ders [ i ]);
            }

            // calculation of jacobian matrix in similar fashion as A4.4
            // calculate values and derivatives of nonrational Bspline surface with weights at first (Aders, wders)
            vertexCoordsPtr = cellgeo.giveVertexCoordinates( mask(k) );
            w = vertexCoordsPtr->at(3);
            xw = vertexCoordsPtr->at(1) * w;
            yw = vertexCoordsPtr->at(2) * w;

            product = ders [ 0 ](0) * ders [ 1 ](0);       // Nu*Nv
            Aders [ 0 ](0, 0) += product * xw;             // x=sum Nu*Nv*x*w
            Aders [ 1 ](0, 0) += product * yw;             // x=sum Nu*Nv*y*w
            wders(0, 0)    += product * w;                 // w=sum Nu*Nv*w

            product = ders [ 0 ](1) * ders [ 1 ](0);       // dNu/du*Nv
            Aders [ 0 ](1, 0) += product * xw; // dx/du=sum dNu/du*Nv*x*w
            Aders [ 1 ](1, 0) += product * yw; // dy/du=sum dNu/du*Nv*y*w
            wders(1, 0)    += product * w;                 // dw/du=sum dNu/du*Nv*w

            product = ders [ 0 ](0) * ders [ 1 ](1);       // Nu*dNv/dv
            Aders [ 0 ](0, 1) += product * xw; // dx/dv=sum Nu*dNv/dv*x*w
            Aders [ 1 ](0, 1) += product * yw; // dy/dv=sum Nu*dNv/dv*y*w
            wders(0, 1)    += product * w; // dw/dv=sum Nu*dNv/dv*w
        }

        weight = wders(0, 0);

        // optimized version of A4.4 for d=1, binomial coefficients ignored
#if 0
        Sders[0](0,0) = Aders[0](0,0) / weight;
        Sders[1](0,0) = Aders[1](0,0) / weight;
        Sders[0](0,1) = (Aders[0](0,1)-wders(0,1)*Sders[0](0,0)) / weight;
        Sders[1](0,1) = (Aders[1](0,1)-wders(0,1)*Sders[1](0,0)) / weight;
        Sders[0](1,0) = (Aders[0](1,0)-wders(1,0)*Sders[0](0,0)) / weight;
        Sders[1](1,0) = (Aders[1](1,0)-wders(1,0)*Sders[1](0,0)) / weight;

        jacobian(0,0) = Sders[0](1,0);   // dx/du
        jacobian(0,1) = Sders[1](1,0);   // dy/du
        jacobian(1,0) = Sders[0](0,1);   // dx/dv
        jacobian(1,1) = Sders[1](0,1);   // dy/dv
#endif

        temp(0) = Aders [ 0 ](0, 0) / weight;
        temp(1) = Aders [ 1 ](0, 0) / weight;
        jacobian(1, 0) = ( Aders [ 0 ](0, 1) - wders(0, 1) * temp(0) ) / weight; // dx/dv
        jacobian(1, 1) = ( Aders [ 1 ](0, 1) - wders(0, 1) * temp(1) ) / weight; // dy/dv
        jacobian(0, 0) = ( Aders [ 0 ](1, 0) - wders(1, 0) * temp(0) ) / weight; // dx/du
        jacobian(0, 1) = ( Aders [ 1 ](1, 0) - wders(1, 0) * temp(1) ) / weight; // dy/du
    }
}


// knotSpan corresponds to knot span in terms of BSpline;
// it should not matter which part of IGAIntegrationElement if covering more spans is addressed;
// this implementation relies on the fact that IGAIntegrationElements are those subsets
// of T-mesh cells on which there are fully (this means not only partially) nonzero relevant basis functions

int TSplineInterpolation :: giveKnotSpanBasisFuncMask(const IntArray &knotSpan, IntArray &mask)
{
    FloatArray knotStart(nsd), knotEnd(nsd);

    // resize the mask initially to the size corresponding to BSpline case
    // but there may be more nonzero basis functions
    if ( nsd == 2 ) {
        mask.preallocate( ( degree [ 0 ] + 1 ) * ( degree [ 1 ] + 1 ) );
    } else {
        OOFEM_ERROR("not implemented for nsd = %d", nsd);
    }

    // get starting and ending knots
    for ( int j = 0; j < nsd; j++ ) {
        knotStart(j) = knotVector [ j ] [ knotSpan(j) ];
        knotEnd(j) = knotVector [ j ] [ knotSpan(j) + 1 ];
    }

    // for each control point check
    for ( int i = 0; i < totalNumberOfControlPoints; i++ ) {
        // whether local knot vector overlaps the given knot span
        int nonzero = 1;
        for ( int j = 0; j < nsd; j++ ) {
            if ( ( knotEnd(j) <= knotValues [ j ].at(localIndexKnotVector [ i ] [ j ] [ 0 ]) ) ||
                ( knotStart(j) >= knotValues [ j ].at(localIndexKnotVector [ i ] [ j ] [ degree [ j ] + 1 ]) ) ) {
                nonzero = 0;
                break;
            }
        }

        if ( nonzero ) {
            mask.followedBy(i + 1, 4);
        }
    }

    return 1;
}



// knotSpan corresponds to knot span in terms of BSpline;
// it should not matter which part of IGAIntegrationElement if covering more spans is addressed;
// this implementation relies on the fact that IGAIntegrationElements are those subsets
// of T-mesh cells on which there are fully (this means not only partially) nonzero relevant basis functions

int TSplineInterpolation :: giveNumberOfKnotSpanBasisFunctions(const IntArray &knotSpan)
{
    int answer = 0;
    FloatArray knotStart(nsd), knotEnd(nsd);

    // get starting and ending knots
    for ( int j = 0; j < nsd; j++ ) {
        knotStart(j) = knotVector [ j ] [ knotSpan(j) ];
        knotEnd(j) = knotVector [ j ] [ knotSpan(j) + 1 ];
    }

    // for each control point check
    for ( int i = 0; i < totalNumberOfControlPoints; i++ ) {
        answer++;
        // whether local knot vector overlaps the given knot span
        for ( int j = 0; j < nsd; j++ ) {
            if ( ( knotEnd(j) <= knotValues [ j ].at(localIndexKnotVector [ i ] [ j ] [ 0 ]) ) ||
                ( knotStart(j) >= knotValues [ j ].at(localIndexKnotVector [ i ] [ j ] [ degree [ j ] + 1 ]) ) ) {
                answer--;
                break;
            }
        }
    }

    return answer;
}



// starKnotSpan and endKnotSpan correspond to knot span in terms of BSpline;
// should the number of non-zero basis function be calculated for single knot span
// starKnotSpan and endKnotSpan are equal;
// some of the basis function may not cover the whole knot span interval !!!

int TSplineInterpolation :: giveKnotSpanBasisFuncMask(const IntArray &startKnotSpan, const IntArray &endKnotSpan, IntArray &mask)
{
    int nonzero;
    FloatArray knotStart(nsd), knotEnd(nsd);

    // resize the mask initially to the size corresponding to BSpline case
    // but there may be more nonzero basis functions
    if ( nsd == 2 ) {
        mask.preallocate( ( degree [ 0 ] + 1 ) * ( degree [ 1 ] + 1 ) );
    } else {
        OOFEM_ERROR("not implemented for nsd = %d", nsd);
    }

    // get starting and ending knots
    for ( int j = 0; j < nsd; j++ ) {
        knotStart(j) = knotVector [ j ] [ startKnotSpan(j) ];
        knotEnd(j) = knotVector [ j ] [ endKnotSpan(j) + 1 ];
    }

    // for each control point check
    for ( int i = 0; i < totalNumberOfControlPoints; i++ ) {
        // whether local knot vector overlaps at least partially the knot span interval
        nonzero = 1;
        for ( int j = 0; j < nsd; j++ ) {
            if ( ( knotEnd(j) <= knotValues [ j ].at(localIndexKnotVector [ i ] [ j ] [ 0 ]) ) ||
                ( knotStart(j) >= knotValues [ j ].at(localIndexKnotVector [ i ] [ j ] [ degree [ j ] + 1 ]) ) ) {
                nonzero = 0;
                break;
            }
        }

        if ( nonzero ) {
            mask.followedBy(i + 1, 4);
        }
    }

    return 1;
}


// starKnotSpan and endKnotSpan correspond to knot apan in terms of BSpline;
// should the number of non-zero basis function be calculated for single knot span
// starKnotSpan and endKnotSpan are equal
// some of the basis function may not cover the whole knot span interval !!!

int TSplineInterpolation :: giveNumberOfKnotSpanBasisFunctions(const IntArray &startKnotSpan, const IntArray &endKnotSpan)
{
    int answer = 0;
    FloatArray knotStart(nsd), knotEnd(nsd);

    // get starting and ending knots
    for ( int j = 0; j < nsd; j++ ) {
        knotStart(j) = knotVector [ j ] [ startKnotSpan(j) ];
        knotEnd(j) = knotVector [ j ] [ endKnotSpan(j) + 1 ];
    }

    // for each control point check
    for ( int i = 0; i < totalNumberOfControlPoints; i++ ) {
        answer++;
        // whether local knot vector overlaps at least partially the knot span interval
        for ( int j = 0; j < nsd; j++ ) {
            if ( ( knotEnd(j) <= knotValues [ j ].at(localIndexKnotVector [ i ] [ j ] [ 0 ]) ) ||
                ( knotStart(j) >= knotValues [ j ].at(localIndexKnotVector [ i ] [ j ] [ degree [ j ] + 1 ]) ) ) {
                answer--;
                break;
            }
        }
    }

    return answer;
}



// call corresponding BSpline methods for open local knot vector

double TSplineInterpolation :: basisFunction(double u, int p, const FloatArray &U, const int *I)
{
    int span, prepend, append;
    FloatArray N;

    createLocalKnotVector(p, U, I, & prepend, & append);
    span = BSplineInterpolation :: findSpan(prepend + append, p, u, openLocalKnotVector);
    BSplineInterpolation :: basisFuns(N, span, u, p, openLocalKnotVector);

    // extract the middle basis function
    // this corresponds to index p-span, however prepended knotspans must be considered
    return N(p - span + prepend);
}


// call corresponding BSpline methods for open local knot vector

void TSplineInterpolation :: dersBasisFunction(int n, double u, int p, const FloatArray &U, const int *I, FloatArray &ders)
{
    int span, prepend, append;
    FloatMatrix Ders;

    createLocalKnotVector(p, U, I, & prepend, & append);
    span = BSplineInterpolation :: findSpan(prepend + append, p, u, openLocalKnotVector);
    BSplineInterpolation :: dersBasisFuns(n, u, span, p, openLocalKnotVector, Ders);

    // extract the middle basis function and its derivatives
    // this corresponds to index p-span, however prepended knotspans must be considered
    ders.resize(n + 1);
    for ( int i = 0; i <= n; i++ ) {
        ders(i) = Ders(i, p - span + prepend);
    }
}


void TSplineInterpolation :: createLocalKnotVector(int p, const FloatArray &U, const int *I, int *prepend, int *append)
{
    int j = 0, index_first = I [ 0 ], index_last = I [ p + 1 ], mult_first = 1, mult_last = 1;
    double first = U.at(index_first), last = U.at(index_last);

    for ( int i = 1; i < p + 1; i++ ) {
        if ( I [ i ] != index_first ) {
            break;
        }

        mult_first++;
    }

    for ( int i = p; i > 0; i-- ) {
        if ( I [ i ] != index_last ) {
            break;
        }

        mult_last++;
    }

    * prepend = p + 1 - mult_first;
    * append = p + 1 - mult_last;

    // prepend first knot (once more)
    for ( int i = 0; i <= * prepend; i++ ) {
        openLocalKnotVector [ j++ ] = first;
    }

    // copy middle of knot vector (without first and last)
    for ( int i = 1; i <= p; i++ ) {
        openLocalKnotVector [ j++ ] = U.at(I [ i ]);
    }

    // append last knot (once more)
    for ( int i = 0; i <= * append; i++ ) {
        openLocalKnotVector [ j++ ] = last;
    }
}
} // end namespace oofem
