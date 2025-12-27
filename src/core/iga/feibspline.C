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
 *               Copyright (C) 1993 - 2025   Borek Patzak
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

#include "floatarray.h"
#include "floatmatrix.h"
#include "iga.h"
#include "feibspline.h"
#include "mathfem.h"
#include "parametermanager.h"
#include "paramkey.h"

namespace oofem {

ParamKey BSplineInterpolation::IPK_BSplineInterpolation_degree("degree");
ParamKey BSplineInterpolation::IPK_BSplineInterpolation_knotVectorU("knotvectoru");
ParamKey BSplineInterpolation::IPK_BSplineInterpolation_knotVectorV("knotvectorv");
ParamKey BSplineInterpolation::IPK_BSplineInterpolation_knotVectorW("knotvectorw");
ParamKey BSplineInterpolation::IPK_BSplineInterpolation_knotMultiplicityU("knotmultiplicityu");
ParamKey BSplineInterpolation::IPK_BSplineInterpolation_knotMultiplicityV("knotmultiplicityv");
ParamKey BSplineInterpolation::IPK_BSplineInterpolation_knotMultiplicityW("knotmultiplicityw");



void
BSplineInterpolation :: initializeFrom(InputRecord &ir, ParameterManager&pm, int elnum, int priority)
{
    bool flag;
    IntArray degree_tmp;
    PM_UPDATE_PARAMETER_AND_REPORT(degree_tmp, pm, ir, elnum, IPK_BSplineInterpolation_degree, priority, flag);
    if (flag) {
        if ( degree_tmp.giveSize() != nsd ) {
            throw ValueInputException(ir, _IFT_BSplineInterpolation_degree, "degree size mismatch");
        }

        for ( int i = 0; i < nsd; i++ ) {
            degree [ i ] = degree_tmp.at(i + 1);
        }
    }
    PM_UPDATE_PARAMETER(knotValues [ 0 ], pm, ir, elnum, IPK_BSplineInterpolation_knotVectorU, priority);
    PM_UPDATE_PARAMETER(knotMultiplicity [ 0 ], pm, ir, elnum, IPK_BSplineInterpolation_knotMultiplicityU, priority);

    if (nsd>1) {
        PM_UPDATE_PARAMETER(knotValues [ 1 ], pm, ir, elnum, IPK_BSplineInterpolation_knotVectorV, priority);
        PM_UPDATE_PARAMETER(knotMultiplicity [ 1 ], pm, ir, elnum, IPK_BSplineInterpolation_knotMultiplicityV, priority);
    }

    if (nsd >2) {
        PM_UPDATE_PARAMETER(knotValues [ 2 ], pm, ir, elnum, IPK_BSplineInterpolation_knotVectorW, priority);
        PM_UPDATE_PARAMETER(knotMultiplicity [ 2 ], pm, ir, elnum, IPK_BSplineInterpolation_knotMultiplicityW, priority);
    }

}

void
BSplineInterpolation::postInitialize(ParameterManager&pm, int elnum)
{
    ParamKey* IPK_knotVector [ 3 ] = {
        &IPK_BSplineInterpolation_knotVectorU,
        &IPK_BSplineInterpolation_knotVectorV,
        &IPK_BSplineInterpolation_knotVectorW
    };

    ParamKey* IPK_knotMultiplicity [ 3 ] = {
        &IPK_BSplineInterpolation_knotMultiplicityU,
        &IPK_BSplineInterpolation_knotMultiplicityV,
        &IPK_BSplineInterpolation_knotMultiplicityW
    };

    for ( int n = 0; n < nsd; n++ ) {
        //IR_GIVE_FIELD(ir, knotValues [ n ], IFT_knotVector [ n ]);
        int size = knotValues [ n ].giveSize();
        if ( size < 2 ) {
            throw ComponentInputException(IPK_knotVector[n]->getName(), ComponentInputException::ComponentType::ctElement, elnum, "invalid size of knot vector");
        }

        // check for monotonicity of knot vector without multiplicity
        double knotVal = knotValues [ n ].at(1);
        for ( int i = 1; i < size; i++ ) {
            if ( knotValues [ n ].at(i + 1) <= knotVal ) {
                throw ComponentInputException(IPK_knotVector [ n ]->getName(), ComponentInputException::ComponentType::ctElement, elnum, "knot vector is not monotonic");
            }

            knotVal = knotValues [ n ].at(i + 1);
        }

        // transform knot vector to interval <0;1>
        double span = knotVal - knotValues [ n ].at(1);
        for ( int i = 1; i <= size; i++ ) {
            knotValues [ n ].at(i) /= span;
        }

        //IR_GIVE_OPTIONAL_FIELD(ir, knotMultiplicity [ n ], IFT_knotMultiplicity [ n ]);
        if ( knotMultiplicity [ n ].giveSize() == 0 ) {
            // default multiplicity
            knotMultiplicity [ n ].resize(size);
            // skip the first and last one
            for ( int i = 1; i < size - 1; i++ ) {
                knotMultiplicity [ n ].at(i + 1) = 1;
            }
        } else {
            if ( knotMultiplicity [ n ].giveSize() != size ) {
                throw ComponentInputException(IPK_knotMultiplicity [ n ]->getName(), ComponentInputException::ComponentType::ctElement, elnum, "knot multiplicity size mismatch");
            }

            // check for multiplicity range (skip the first and last one)
            for ( int i = 1; i < size - 1; i++ ) {
                if ( knotMultiplicity [ n ].at(i + 1) < 1 || knotMultiplicity [ n ].at(i + 1) > degree [ n ] ) {
                    throw ComponentInputException(IPK_knotMultiplicity [ n ]->getName(), ComponentInputException::ComponentType::ctElement, elnum, "knot multiplicity out of range");
                }
            }

            // check for multiplicity of the first and last one
            if ( knotMultiplicity [ n ].at(1) != degree [ n ] + 1 ) {
                OOFEM_LOG_RELEVANT("Multiplicity of the first knot in knot vector %s changed to %d\n", IPK_knotVector [ n ]->getNameCStr(), degree [ n ] + 1);
            }

            if ( knotMultiplicity [ n ].at(size) != degree [ n ] + 1 ) {
                OOFEM_LOG_RELEVANT("Multiplicity of the last knot in knot vector %s changed to %d\n", IPK_knotVector [ n ]->getNameCStr(), degree [ n ] + 1);
            }
        }

        // multiplicity of the 1st and last knot set to degree + 1
        knotMultiplicity [ n ].at(1) = knotMultiplicity [ n ].at(size) = degree [ n ] + 1;

        // sum the size of knot vector with multiplicity values
        int sum = 0;
        for ( int i = 0; i < size; i++ ) {
            sum += knotMultiplicity [ n ][i];
        }

        knotVector [ n ].resize( sum );

        // fill knot vector including multiplicity values
        int pos = 0;
        for ( int i = 0; i < size; i++ ) {
            for ( int j = 0; j < knotMultiplicity [ n ].at(i + 1); j++ ) {
                knotVector [ n ] [ pos++ ] = knotValues [ n ].at(i + 1);
            }
        }

        numberOfKnotSpans [ n ] = size - 1;
        numberOfControlPoints [ n ] = sum - degree [ n ] - 1;
    }
}
void BSplineInterpolation :: evalN(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const
{
    const FEIIGAElementGeometryWrapper &gw = static_cast< const FEIIGAElementGeometryWrapper& >(cellgeo);
    IntArray span(nsd);
    int c = 1;
    std :: vector< FloatArray > N(nsd);

    if ( gw.knotSpan ) {
        span = * gw.knotSpan;
    } else {
        for ( int i = 0; i < nsd; i++ ) {
            span[i] = this->findSpan(numberOfControlPoints [ i ], degree [ i ], lcoords[i], knotVector [ i ]);
        }
    }

    for ( int i = 0; i < nsd; i++ ) {
        this->basisFuns(N [ i ], span[i], lcoords[i], degree [ i ], knotVector [ i ]);
    }

    answer.resize(giveNumberOfKnotSpanBasisFunctions(span));

    if ( nsd == 1 ) {
        for ( int k = 0; k <= degree [ 0 ]; k++ ) {
            answer.at(c++) = N [ 0 ][k];
        }
    } else if ( nsd == 2 ) {
        for ( int l = 0; l <= degree [ 1 ]; l++ ) {
            for ( int k = 0; k <= degree [ 0 ]; k++ ) {
                answer.at(c++) = N [ 0 ][k] * N [ 1 ][l];
            }
        }
    } else if ( nsd == 3 ) {
        for ( int m = 0; m <= degree [ 2 ]; m++ ) {
            for ( int l = 0; l <= degree [ 1 ]; l++ ) {
                for ( int k = 0; k <= degree [ 0 ]; k++ ) {
                    answer.at(c++) = N [ 0 ][k] * N [ 1 ][l] * N [ 2 ][m];
                }
            }
        }
    } else {
        OOFEM_ERROR("evalN not implemented for nsd = %d", nsd);
    }
}


double BSplineInterpolation :: evaldNdx(FloatMatrix &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const
{
    const FEIIGAElementGeometryWrapper &gw = static_cast< const FEIIGAElementGeometryWrapper& >(cellgeo);
    FloatMatrix jacobian(nsd, nsd);
    IntArray span(nsd);
    double Jacob = 0.;
    std :: vector< FloatMatrix > ders(nsd);

    if ( gw.knotSpan ) {
        span = * gw.knotSpan;
    } else {
        for ( int i = 0; i < nsd; i++ ) {
            span[i] = this->findSpan(numberOfControlPoints [ i ], degree [ i ], lcoords[i], knotVector [ i ]);
        }
    }

    for ( int i = 0; i < nsd; i++ ) {
        this->dersBasisFuns(1, lcoords[i], span[i], degree [ i ], knotVector [ i ], ders [ i ]);
    }

    int count = giveNumberOfKnotSpanBasisFunctions(span);
    answer.resize(count, nsd);
    jacobian.zero();

    if ( nsd == 1 ) {
        int uind = span[0] - degree [ 0 ];
        int ind = uind + 1;
        for ( int k = 0; k <= degree [ 0 ]; k++ ) {
            const auto &vertexCoords = cellgeo.giveVertexCoordinates(ind + k);
            jacobian(0, 0) += ders [ 0 ](1, k) * vertexCoords[0];       // dx/du=sum(dNu/du*x)
        }

        Jacob = jacobian.giveDeterminant();

        if ( fabs(Jacob) < 1.0e-10 ) {
            OOFEM_ERROR("evaldNdx - zero Jacobian");
        }

        int cnt = 0;
        for ( int k = 0; k <= degree [ 0 ]; k++ ) {
            answer(cnt, 0) = ders [ 0 ](1, k) / Jacob;         // dN/dx=dN/du / dx/du
            cnt++;
        }
    } else if ( nsd == 2 ) {
        FloatArray tmp1(nsd), tmp2(nsd);

        int uind = span[0] - degree [ 0 ];
        int vind = span[1] - degree [ 1 ];
        int ind = vind * numberOfControlPoints [ 0 ] + uind + 1;
        for ( int l = 0; l <= degree [ 1 ]; l++ ) {
            tmp1.zero();
            tmp2.zero();
            for ( int k = 0; k <= degree [ 0 ]; k++ ) {
                const auto &vertexCoords = cellgeo.giveVertexCoordinates(ind + k);

                tmp1[0] += ders [ 0 ](1, k) * vertexCoords[0];            // sum(dNu/du*x)
                tmp1[1] += ders [ 0 ](1, k) * vertexCoords[1]; // sum(dNu/du*y)

                tmp2[0] += ders [ 0 ](0, k) * vertexCoords[0]; // sum(Nu*x)
                tmp2[1] += ders [ 0 ](0, k) * vertexCoords[1]; // sum(Nu*y)
            }

            ind += numberOfControlPoints [ 0 ];

            jacobian(0, 0) += ders [ 1 ](0, l) * tmp1[0]; // dx/du=sum(Nv*sum(dNu/du*x))
            jacobian(0, 1) += ders [ 1 ](0, l) * tmp1[1]; // dy/du=sum(Nv*sum(dNu/du*y))

            jacobian(1, 0) += ders [ 1 ](1, l) * tmp2[0]; // dx/dv=sum(dNv/dv*sum(Nu*x))
            jacobian(1, 1) += ders [ 1 ](1, l) * tmp2[1]; // dy/dv=sum(dNv/dv*sum(Nu*y))
        }

        Jacob = jacobian.giveDeterminant();

        if ( fabs(Jacob) < 1.0e-10 ) {
            OOFEM_ERROR("evaldNdx - zero Jacobian");
        }

        int cnt = 0;
        for ( int l = 0; l <= degree [ 1 ]; l++ ) {
            for ( int k = 0; k <= degree [ 0 ]; k++ ) {
                tmp1[0] = ders [ 0 ](1, k) * ders [ 1 ](0, l); // dN/du=dNu/du*Nv
                tmp1[1] = ders [ 0 ](0, k) * ders [ 1 ](1, l); // dN/dv=Nu*dNv/dv
                answer(cnt, 0) = ( +jacobian(1, 1) * tmp1[0] - jacobian(0, 1) * tmp1[1] ) / Jacob; // dN/dx
                answer(cnt, 1) = ( -jacobian(1, 0) * tmp1[0] + jacobian(0, 0) * tmp1[1] ) / Jacob; // dN/dy
                cnt++;
            }
        }
    } else if ( nsd == 3 ) {
        FloatArray tmp1(nsd), tmp2(nsd);
        FloatArray temp1(nsd), temp2(nsd), temp3(nsd);

        int uind = span[0] - degree [ 0 ];
        int vind = span[1] - degree [ 1 ];
        int tind = span[2] - degree [ 2 ];
        int ind = tind * numberOfControlPoints [ 0 ] * numberOfControlPoints [ 1 ] + vind * numberOfControlPoints [ 0 ] + uind + 1;
        for ( int m = 0; m <= degree [ 2 ]; m++ ) {
            temp1.zero();
            temp2.zero();
            temp3.zero();
            int indx = ind;
            for ( int l = 0; l <= degree [ 1 ]; l++ ) {
                tmp1.zero();
                tmp2.zero();
                for ( int k = 0; k <= degree [ 0 ]; k++ ) {
                    const auto &vertexCoords = cellgeo.giveVertexCoordinates(ind + k);

                    tmp1[0] += ders [ 0 ](1, k) * vertexCoords[0];                // sum(dNu/du*x)
                    tmp1[1] += ders [ 0 ](1, k) * vertexCoords[1];                // sum(dNu/du*y)
                    tmp1[2] += ders [ 0 ](1, k) * vertexCoords[2];                // sum(dNu/du*z)

                    tmp2[0] += ders [ 0 ](0, k) * vertexCoords[0];                // sum(Nu*x)
                    tmp2[1] += ders [ 0 ](0, k) * vertexCoords[1];                // sum(Nu*y)
                    tmp2[2] += ders [ 0 ](0, k) * vertexCoords[2];                // sum(Nu*y)
                }

                ind += numberOfControlPoints [ 0 ];

                temp1[0] += ders [ 1 ](0, l) * tmp1[0];            // sum(Nv*sum(dNu/du*x))
                temp1[1] += ders [ 1 ](0, l) * tmp1[1];            // sum(Nv*sum(dNu/du*y))
                temp1[2] += ders [ 1 ](0, l) * tmp1[2];            // sum(Nv*sum(dNu/du*z))

                temp2[0] += ders [ 1 ](1, l) * tmp2[0];            // sum(dNv/dv*sum(Nu*x))
                temp2[1] += ders [ 1 ](1, l) * tmp2[1];            // sum(dNv/dv*sum(Nu*y))
                temp2[2] += ders [ 1 ](1, l) * tmp2[1];            // sum(dNv/dv*sum(Nu*z))

                temp3[0] += ders [ 1 ](0, l) * tmp2[0];            // sum(Nv*sum(Nu*x))
                temp3[1] += ders [ 1 ](0, l) * tmp2[1];            // sum(Nv*sum(Nu*y))
                temp3[2] += ders [ 1 ](0, l) * tmp2[1];            // sum(Nv*sum(Nu*z))
            }

            ind = indx + numberOfControlPoints [ 0 ] * numberOfControlPoints [ 1 ];

            jacobian(0, 0) += ders [ 2 ](0, m) * temp1[0]; // dx/du=sum(Nt*sum(Nv*sum(dNu/du*x)))
            jacobian(0, 1) += ders [ 2 ](0, m) * temp1[1]; // dy/du=sum(Nt*sum(Nv*sum(dNu/du*y)))
            jacobian(0, 2) += ders [ 2 ](0, m) * temp1[2]; // dz/du=sum(Nt*sum(Nv*sum(dNu/du*z)))

            jacobian(1, 0) += ders [ 2 ](0, m) * temp2[0]; // dx/dv=sum(Nt*sum(dNv/dv*sum(Nu*x)))
            jacobian(1, 1) += ders [ 2 ](0, m) * temp2[1]; // dy/dv=sum(Nt*sum(dNv/dv*sum(Nu*y)))
            jacobian(1, 2) += ders [ 2 ](0, m) * temp2[2]; // dz/dv=sum(Nt*sum(dNv/dv*sum(Nu*z)))

            jacobian(2, 0) += ders [ 2 ](1, m) * temp3[0]; // dx/dt=sum(dNt/dt*sum(Nv*sum(Nu*x)))
            jacobian(2, 1) += ders [ 2 ](1, m) * temp3[1]; // dy/dt=sum(dNt/dt*sum(Nv*sum(Nu*y)))
            jacobian(2, 2) += ders [ 2 ](1, m) * temp3[2]; // dz/dt=sum(dNt/dt*sum(Nv*sum(Nu*z)))
        }

        Jacob = jacobian.giveDeterminant();

        if ( fabs(Jacob) < 1.0e-10 ) {
            OOFEM_ERROR("evaldNdx - zero Jacobian");
        }

        int cnt = 0;
        for ( int m = 0; m <= degree [ 2 ]; m++ ) {
            for ( int l = 0; l <= degree [ 1 ]; l++ ) {
                for ( int k = 0; k <= degree [ 0 ]; k++ ) {
                    tmp1[0] = ders [ 0 ](1, k) * ders [ 1 ](0, l) * ders [ 2 ](0, m);       // dN/du=dNu/du*Nv*Nt
                    tmp1[1] = ders [ 0 ](0, k) * ders [ 1 ](1, l) * ders [ 2 ](0, m);       // dN/dv=Nu*dNv/dv*Nt
                    tmp1[2] = ders [ 0 ](0, k) * ders [ 1 ](0, l) * ders [ 2 ](1, m);       // dN/dt=Nu*Nv*dNt/dt
                    answer(cnt, 0) = ( ( jacobian(1, 1) * jacobian(2, 2) - jacobian(1, 2) * jacobian(2, 1) ) * tmp1[0] +
                                      ( jacobian(0, 2) * jacobian(2, 1) - jacobian(0, 1) * jacobian(2, 2) ) * tmp1[1] +
                                      ( jacobian(0, 1) * jacobian(1, 2) - jacobian(0, 2) * jacobian(1, 1) ) * tmp1[2] ) / Jacob;  // dN/dx
                    answer(cnt, 1) = ( ( jacobian(1, 2) * jacobian(2, 0) - jacobian(1, 0) * jacobian(2, 2) ) * tmp1[0] +
                                      ( jacobian(0, 0) * jacobian(2, 2) - jacobian(0, 2) * jacobian(2, 0) ) * tmp1[1] +
                                      ( jacobian(0, 2) * jacobian(1, 0) - jacobian(0, 0) * jacobian(1, 2) ) * tmp1[2] ) / Jacob;  // dN/dy
                    answer(cnt, 2) = ( ( jacobian(1, 0) * jacobian(2, 1) - jacobian(1, 1) * jacobian(2, 0) ) * tmp1[0] +
                                      ( jacobian(0, 1) * jacobian(2, 0) - jacobian(0, 0) * jacobian(2, 1) ) * tmp1[1] +
                                      ( jacobian(0, 0) * jacobian(1, 1) - jacobian(0, 1) * jacobian(1, 0) ) * tmp1[2] ) / Jacob;  // dN/dz
                    cnt++;
                }
            }
        }
    } else {
        OOFEM_ERROR("evaldNdx not implemented for nsd = %d", nsd);
    }

    return Jacob;
}


void BSplineInterpolation :: local2global(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const
{
    /* Based on SurfacePoint A3.5 implementation*/
    const FEIIGAElementGeometryWrapper &gw = static_cast< const FEIIGAElementGeometryWrapper& >(cellgeo);
    IntArray span(nsd);
    std :: vector< FloatArray > N(nsd);


    if ( gw.knotSpan ) {
        span = * gw.knotSpan;
    } else {
        for ( int i = 0; i < nsd; i++ ) {
            span[i] = this->findSpan(numberOfControlPoints [ i ], degree [ i ], lcoords[i], knotVector [ i ]);
        }
    }

    for ( int i = 0; i < nsd; i++ ) {
        this->basisFuns(N [ i ], span[i], lcoords[i], degree [ i ], knotVector [ i ]);
    }

    answer.resize(nsd);
    answer.zero();

    if ( nsd == 1 ) {
        int uind = span[0] - degree [ 0 ];
        int ind = uind + 1;
        for ( int k = 0; k <= degree [ 0 ]; k++ ) {
            const auto &vertexCoords = cellgeo.giveVertexCoordinates(ind + k);
            answer[0] += N [ 0 ][k] * vertexCoords[0];
        }
    } else if ( nsd == 2 ) {
        FloatArray tmp(nsd);

        int uind = span[0] - degree [ 0 ];
        int vind = span[1] - degree [ 1 ];
        int ind = vind * numberOfControlPoints [ 0 ] + uind + 1;
        for ( int l = 0; l <= degree [ 1 ]; l++ ) {
            tmp.zero();
            for ( int k = 0; k <= degree [ 0 ]; k++ ) {
                const auto &vertexCoords = cellgeo.giveVertexCoordinates(ind + k);

                tmp[0] += N [ 0 ][k] * vertexCoords[0];
                tmp[1] += N [ 0 ][k] * vertexCoords[1];
            }

            ind += numberOfControlPoints [ 0 ];

            answer[0] += N [ 1 ][l] * tmp[0];
            answer[1] += N [ 1 ][l] * tmp[1];
        }
    } else if ( nsd == 3 ) {
        FloatArray tmp(nsd), temp(nsd);

        int uind = span[0] - degree [ 0 ];
        int vind = span[1] - degree [ 1 ];
        int tind = span[2] - degree [ 2 ];
        int ind = tind * numberOfControlPoints [ 0 ] * numberOfControlPoints [ 1 ] + vind * numberOfControlPoints [ 0 ] + uind + 1;
        for ( int m = 0; m <= degree [ 2 ]; m++ ) {
            temp.zero();
            int indx = ind;
            for ( int l = 0; l <= degree [ 1 ]; l++ ) {
                tmp.zero();
                for ( int k = 0; k <= degree [ 0 ]; k++ ) {
                    const auto &vertexCoords = cellgeo.giveVertexCoordinates(ind + k);
                    tmp.add(N [ 0 ][k], vertexCoords);
                }

                ind += numberOfControlPoints [ 0 ];

                temp.add(N [ 1 ][l], tmp);
            }

            ind = indx + numberOfControlPoints [ 0 ] * numberOfControlPoints [ 1 ];

            answer.add(N [ 2 ][m], temp);
        }
    } else {
        OOFEM_ERROR("local2global not implemented for nsd = %d", nsd);
    }
}


void BSplineInterpolation :: giveJacobianMatrixAt(FloatMatrix &jacobian, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const
{
    const FEIIGAElementGeometryWrapper &gw = static_cast< const FEIIGAElementGeometryWrapper& >(cellgeo);
    IntArray span(nsd);
    std :: vector< FloatMatrix > ders(nsd);
    jacobian.resize(nsd, nsd);

    if ( gw.knotSpan ) {
        span = * gw.knotSpan;
    } else {
        for ( int i = 0; i < nsd; i++ ) {
            span[i] = this->findSpan(numberOfControlPoints [ i ], degree [ i ], lcoords[i], knotVector [ i ]);
        }
    }

    for ( int i = 0; i < nsd; i++ ) {
        this->dersBasisFuns(1, lcoords[i], span[i], degree [ i ], knotVector [ i ], ders [ i ]);
    }

    jacobian.zero();

    if ( nsd == 1 ) {
        int uind = span[0] - degree [ 0 ];
        int ind = uind + 1;
        for ( int k = 0; k <= degree [ 0 ]; k++ ) {
            const auto &vertexCoords = cellgeo.giveVertexCoordinates(ind + k);
            jacobian(0, 0) += ders [ 0 ](1, k) * vertexCoords[0];       // dx/du=sum(dNu/du*x)
        }
    } else if ( nsd == 2 ) {
        FloatArray tmp1(nsd), tmp2(nsd);

        int uind = span[0] - degree [ 0 ];
        int vind = span[1] - degree [ 1 ];
        int ind = vind * numberOfControlPoints [ 0 ] + uind + 1;
        for ( int l = 0; l <= degree [ 1 ]; l++ ) {
            tmp1.zero();
            tmp2.zero();
            for ( int k = 0; k <= degree [ 0 ]; k++ ) {
                const auto &vertexCoords = cellgeo.giveVertexCoordinates(ind + k);

                tmp1[0] += ders [ 0 ](1, k) * vertexCoords[0]; // sum(dNu/du*x)
                tmp1[1] += ders [ 0 ](1, k) * vertexCoords[1]; // sum(dNu/du*y)

                tmp2[0] += ders [ 0 ](0, k) * vertexCoords[0]; // sum(Nu*x)
                tmp2[1] += ders [ 0 ](0, k) * vertexCoords[1]; // sum(Nu*y)
            }

            ind += numberOfControlPoints [ 0 ];

            jacobian(0, 0) += ders [ 1 ](0, l) * tmp1[0]; // dx/du=sum(Nv*sum(dNu/du*x))
            jacobian(0, 1) += ders [ 1 ](0, l) * tmp1[1]; // dy/du=sum(Nv*sum(dNu/du*y))

            jacobian(1, 0) += ders [ 1 ](1, l) * tmp2[0]; // dx/dv=sum(dNv/dv*sum(Nu*x))
            jacobian(1, 1) += ders [ 1 ](1, l) * tmp2[1]; // dy/dv=sum(dNv/dv*sum(Nu*y))
        }
    } else if ( nsd == 3 ) {
        FloatArray tmp1(nsd), tmp2(nsd);
        FloatArray temp1(nsd), temp2(nsd), temp3(nsd);

        int uind = span[0] - degree [ 0 ];
        int vind = span[1] - degree [ 1 ];
        int tind = span[2] - degree [ 2 ];
        int ind = tind * numberOfControlPoints [ 0 ] * numberOfControlPoints [ 1 ] + vind * numberOfControlPoints [ 0 ] + uind + 1;
        for ( int m = 0; m <= degree [ 2 ]; m++ ) {
            temp1.zero();
            temp2.zero();
            temp3.zero();
            int indx = ind;
            for ( int l = 0; l <= degree [ 1 ]; l++ ) {
                tmp1.zero();
                tmp2.zero();
                for ( int k = 0; k <= degree [ 0 ]; k++ ) {
                    const auto &vertexCoords = cellgeo.giveVertexCoordinates(ind + k);

                    tmp1.add(ders [ 0 ](1, k), vertexCoords); // sum(dNu/du*x)
                    tmp2.add(ders [ 0 ](0, k), vertexCoords); // sum(Nu*x)
                }

                ind += numberOfControlPoints [ 0 ];

                temp1.add(ders [ 1 ](0, l), tmp1); // sum(Nv*sum(dNu/du*x)
                temp2.add(ders [ 1 ](1, l), tmp2); // sum(dNv/dv*sum(Nu*x)
                temp3.add(ders [ 1 ](0, l), tmp2); // sum(Nv*sum(Nu*x)
            }

            ind = indx + numberOfControlPoints [ 0 ] * numberOfControlPoints [ 1 ];

            jacobian(0, 0) += ders [ 2 ](0, m) * temp1[0]; // dx/du=sum(Nt*sum(Nv*sum(dNu/du*x)))
            jacobian(0, 1) += ders [ 2 ](0, m) * temp1[1]; // dy/du=sum(Nt*sum(Nv*sum(dNu/du*y)))
            jacobian(0, 2) += ders [ 2 ](0, m) * temp1[2]; // dz/du=sum(Nt*sum(Nv*sum(dNu/du*z)))

            jacobian(1, 0) += ders [ 2 ](0, m) * temp2[0]; // dx/dv=sum(Nt*sum(dNv/dv*sum(Nu*x)))
            jacobian(1, 1) += ders [ 2 ](0, m) * temp2[1]; // dy/dv=sum(Nt*sum(dNv/dv*sum(Nu*y)))
            jacobian(1, 2) += ders [ 2 ](0, m) * temp2[2]; // dz/dv=sum(Nt*sum(dNv/dv*sum(Nu*z)))

            jacobian(2, 0) += ders [ 2 ](1, m) * temp3[0]; // dx/dt=sum(dNt/dt*sum(Nv*sum(Nu*x)))
            jacobian(2, 1) += ders [ 2 ](1, m) * temp3[1]; // dy/dt=sum(dNt/dt*sum(Nv*sum(Nu*y)))
            jacobian(2, 2) += ders [ 2 ](1, m) * temp3[2]; // dz/dt=sum(dNt/dt*sum(Nv*sum(Nu*z)))
        }
    } else {
        OOFEM_ERROR("giveTransformationJacobian not implemented for nsd = %d", nsd);
    }
}


int BSplineInterpolation :: giveKnotSpanBasisFuncMask(const IntArray &knotSpan, IntArray &mask) const
{
    int c = 1;
    int size = giveNumberOfKnotSpanBasisFunctions(knotSpan);
    mask.resize(size);

    if ( nsd == 1 ) {
        for ( int i = 0; i <= degree [ 0 ]; i++ ) {
            int iindx = ( i + knotSpan[0] - degree [ 0 ] );
            mask.at(c++) = iindx + 1;
        }
    } else if ( nsd == 2 ) {
        for ( int j = 0; j <= degree [ 1 ]; j++ ) {
            int jindx = ( j + knotSpan[1] - degree [ 1 ] );
            for ( int i = 0; i <= degree [ 0 ]; i++ ) {
                int iindx = ( i + knotSpan[0] - degree [ 0 ] );
                mask.at(c++) = jindx * numberOfControlPoints [ 0 ] + iindx + 1;
            }
        }
    } else if ( nsd == 3 ) {
        for ( int k = 0; k <= degree [ 2 ]; k++ ) {
            int kindx = ( k + knotSpan[2] - degree [ 2 ] );
            for ( int j = 0; j <= degree [ 1 ]; j++ ) {
                int jindx = ( j + knotSpan[1] - degree [ 1 ] );
                for ( int i = 0; i <= degree [ 0 ]; i++ ) {
                    int iindx = ( i + knotSpan[0] - degree [ 0 ] );
                    mask.at(c++) = kindx * numberOfControlPoints [ 0 ] * numberOfControlPoints [ 1 ] + jindx * numberOfControlPoints [ 0 ] + iindx + 1;
                }
            }
        }
    } else {
        OOFEM_ERROR("not implemented for nsd = %d", nsd);
    }
    return 1;
}


// for pure Bspline the number of nonzero basis functions is the same for each knot span
int BSplineInterpolation :: giveNumberOfKnotSpanBasisFunctions(const IntArray &knotSpan) const
{
    int answer = 1;
    // there are always degree+1 nonzero basis functions on each knot span
    ///@todo This loop seems meaningless. It just returns degree[nsd-1]+1 in the end ?
    for ( int i = 0; i < nsd; i++ ) {
        answer *= ( degree [ i ] + 1 );
    }

    return answer;
}


// generally it is redundant to pass p and U as these data are part of BSplineInterpolation
// and can be retrieved for given spatial dimension;
// however in such a case this function could not be used for calculation on local knot vector of TSpline;
// it is also redundant to pass the span which can be calculated
// but we want to profit from knowing the span a priori

void BSplineInterpolation :: basisFuns(FloatArray &N, int span, double u, int p, const FloatArray &U) const
{
    //
    // Based on Algorithm A2.2 (p. 70)
    //
    FloatArray right(p + 1);
    FloatArray left(p + 1);

    N.resize(p + 1);
    N[0] = 1.0;
    for ( int j = 1; j <= p; j++ ) {
        left[j] = u - U [ span + 1 - j ];
        right[j] = U [ span + j ] - u;
        double saved = 0.0;
        for ( int r = 0; r < j; r++ ) {
            double temp = N[r] / ( right(r + 1) + left(j - r) );
            N[r] = saved + right(r + 1) * temp;
            saved = left(j - r) * temp;
        }

        N[j] = saved;
    }
}


// generally it is redundant to pass p and U as these data are part of BSplineInterpolation
// and can be retrieved for given spatial dimension;
// however in such a case this function could not be used for calculation on local knot vector of TSpline;
// it is also redundant to pass the span which can be calculated
// but we want to profit from knowing the span a priori

void BSplineInterpolation :: dersBasisFuns(int n, double u, int span, int p, const FloatArray &U, FloatMatrix &ders) const
{
    //
    // Based on Algorithm A2.3 (p. 72)
    //
    FloatArray left(p + 1);
    FloatArray right(p + 1);
    FloatMatrix ndu(p + 1, p + 1);

    ders.resize(n + 1, p + 1);

    ndu(0, 0) = 1.0;
    for ( int j = 1; j <= p; j++ ) {
        left[j] = u - U [ span + 1 - j ];
        right[j] = U [ span + j ] - u;
        double saved = 0.0;

        for ( int r = 0; r < j; r++ ) {
            // Lower triangle
            ndu(j, r) = right(r + 1) + left(j - r);
            double temp = ndu(r, j - 1) / ndu(j, r);
            // Upper triangle
            ndu(r, j) = saved + right(r + 1) * temp;
            saved = left(j - r) * temp;
        }

        ndu(j, j) = saved;
    }

    for ( int j = 0; j <= p; j++ ) {
        ders(0, j) = ndu(j, p);
    }

    // Compute the derivatives
    FloatMatrix a(2, p + 1);
    for ( int r = 0; r <= p; r++ ) {
        int s1, s2;
        s1 = 0;
        s2 = 1;       // alternate rows in array a
        a(0, 0) = 1.0;
        // Compute the kth derivative
        for ( int k = 1; k <= n; k++ ) {
            double d;
            int rk, pk, j1, j2;
            d = 0.0;
            rk = r - k;
            pk = p - k;

            if ( r >= k ) {
                a(s2, 0) = a(s1, 0) / ndu(pk + 1, rk);
                d = a(s2, 0) * ndu(rk, pk);
            }

            if ( rk >= -1 ) {
                j1 = 1;
            } else {
                j1 = -rk;
            }

            if ( r - 1 <= pk ) {
                j2 = k - 1;
            } else {
                j2 = p - r;
            }

            for ( int j = j1; j <= j2; j++ ) {
                a(s2, j) = ( a(s1, j) - a(s1, j - 1) ) / ndu(pk + 1, rk + j);
                d += a(s2, j) * ndu(rk + j, pk);
            }

            if ( r <= pk ) {
                a(s2, k) = -a(s1, k - 1) / ndu(pk + 1, r);
                d += a(s2, k) * ndu(r, pk);
            }

            ders(k, r) = d;
            std :: swap(s1, s2); // Switch rows
        }
    }

    // Multiply through by the correct factors
    int r = p;
    for ( int k = 1; k <= n; k++ ) {
        for ( int j = 0; j <= p; j++ ) {
            ders(k, j) *= r;
        }

        r *= p - k;
    }
}



// generally it is redundant to pass n, p and U as these data are part of BSplineInterpolation
// and can be retrieved for given spatial dimension;
// however in such a case this function could not be used for span localization in local knot vector of TSpline

// jaky ma vyznam const = ve funkci se objekt nesmi zmenit

int BSplineInterpolation :: findSpan(int n, int p, double u, const FloatArray &U) const
{
    if ( u == U [ n + 1 ] ) {
        return n;
    }

    int low  = p;
    int high = n + 1;
    int mid = ( low + high ) / 2;

    while ( u < U [ mid ] || u >= U [ mid + 1 ] ) {
        if ( u < U [ mid ] ) {
            high = mid;
        } else {
            low = mid;
        }

        mid = ( low + high ) / 2;
    }

    return mid;
}
} // end namespace oofem
