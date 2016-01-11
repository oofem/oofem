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

#include "mmaleastsquareprojection.h"
#include "mathfem.h"
#include "gausspoint.h"
#include "element.h"
#include "domain.h"
#include "spatiallocalizer.h"
#include "integrationrule.h"
#include "connectivitytable.h"
#include "dynamicinputrecord.h"

namespace oofem {
MMALeastSquareProjection :: MMALeastSquareProjection() : MaterialMappingAlgorithm()
{
    this->stateFilter = 0;
    this->regionFilter = 1;
}

MMALeastSquareProjection :: ~MMALeastSquareProjection() { }

void
MMALeastSquareProjection :: __init(Domain *dold, IntArray &type, FloatArray &coords, Set &elemSet, TimeStep *tStep, bool iCohesiveZoneGP)
//(Domain* dold, IntArray& varTypes, GaussPoint* gp, TimeStep* tStep)
{
    GaussPoint *sourceIp;
    Element *sourceElement;
    SpatialLocalizer *sl = dold->giveSpatialLocalizer();
    IntegrationRule *iRule;

    IntArray patchList;

    this->patchDomain = dold;
    // find the closest IP on old mesh
    sourceElement = sl->giveElementContainingPoint(coords, elemSet);

    if ( !sourceElement ) {
        OOFEM_ERROR("no suitable source element found");
    }

    // determine the type of patch
    Element_Geometry_Type egt = sourceElement->giveGeometryType();
    if ( egt == EGT_line_1 ) {
        this->patchType = MMALSPPatchType_1dq;
    } else if ( ( egt == EGT_triangle_1 ) || ( egt == EGT_quad_1 ) ) {
        this->patchType = MMALSPPatchType_2dq;
    } else {
        OOFEM_ERROR("unsupported material mode");
    }

    /* Determine the state of closest point.
     * Only IP in the neighbourhood with same state can be used
     * to interpolate the values.
     */
    FloatArray dam;
    int state = 0;
    if ( this->stateFilter ) {
        iRule = sourceElement->giveDefaultIntegrationRulePtr();
        for ( GaussPoint *gp: *iRule ) {
            sourceElement->giveIPValue(dam, gp, IST_PrincipalDamageTensor, tStep);
            if ( dam.computeNorm() > 1.e-3 ) {
                state = 1; // damaged
            }
        }
    }

    // from source neighbours the patch will be constructed
    Element *element;
    IntArray neighborList;
    patchList.resize(1);
    patchList.at(1) = sourceElement->giveNumber();
    int minNumberOfPoints = this->giveNumberOfUnknownPolynomialCoefficients(this->patchType);
    int actualNumberOfPoints = sourceElement->giveDefaultIntegrationRulePtr()->giveNumberOfIntegrationPoints();
    int nite = 0;
    int elemFlag;
    // check if number of IP in patchList is sufficient
    // some recursion control would be appropriate
    while ( ( actualNumberOfPoints < minNumberOfPoints ) && ( nite <= 2 ) ) {
        //if not,  construct the neighborhood
        dold->giveConnectivityTable()->giveElementNeighbourList(neighborList, patchList);
        // count number of available points
        patchList.clear();
        actualNumberOfPoints = 0;
        for ( int i = 1; i <= neighborList.giveSize(); i++ ) {
            if ( this->stateFilter ) {
                element = patchDomain->giveElement( neighborList.at(i) );
                // exclude elements in different regions
                if ( !elemSet.hasElement( element->giveNumber() ) ) {
                    continue;
                }

                iRule = element->giveDefaultIntegrationRulePtr();
                elemFlag = 0;
                for ( GaussPoint *gp: *iRule ) {
                    element->giveIPValue(dam, gp, IST_PrincipalDamageTensor, tStep);
                    if ( state && ( dam.computeNorm() > 1.e-3 ) ) {
                        actualNumberOfPoints++;
                        elemFlag = 1;
                    } else if ( ( state == 0 ) && ( dam.computeNorm() < 1.e-3 ) ) {
                        actualNumberOfPoints++;
                        elemFlag = 1;
                    }
                }

                if ( elemFlag ) {
                    // include this element with corresponding state in neighbor search.
                    patchList.followedBy(neighborList.at(i), 10);
                }
            } else { // if (! yhis->stateFilter)
                element = patchDomain->giveElement( neighborList.at(i) );
                // exclude elements in different regions
                if ( !elemSet.hasElement( element->giveNumber() ) ) {
                    continue;
                }

                actualNumberOfPoints += element->giveDefaultIntegrationRulePtr()->giveNumberOfIntegrationPoints();

                patchList.followedBy(neighborList.at(i), 10);
            }
        } // end loop over neighbor list

        nite++;
    }

    if ( nite > 2 ) {
        // not enough points -> take closest point projection
        patchGPList.clear();
        sourceIp = sl->giveClosestIP(coords, elemSet);
        patchGPList.push_front(sourceIp);
        //fprintf(stderr, "MMALeastSquareProjection: too many neighbor search iterations\n");
        //exit (1);
        return;
    }

#ifdef MMALSP_ONLY_CLOSEST_POINTS
    // select only the nval closest IP points
    GaussPoint **gpList = ( GaussPoint ** ) malloc(sizeof( GaussPoint * ) * actualNumberOfPoints);
    FloatArray dist(actualNumberOfPoints), srcgpcoords;
    int npoints = 0;
    // check allocation of gpList
    if ( gpList == NULL ) {
        OOFEM_FATAL("memory allocation error");
    }

    for ( int ielem = 1; ielem <= patchList.giveSize(); ielem++ ) {
        element = patchDomain->giveElement( patchList.at(ielem) );
        iRule = element->giveDefaultIntegrationRulePtr();
        for ( GaussPoint *srcgp: *iRule ) {
            if ( element->computeGlobalCoordinates( srcgpcoords, * ( srcgp->giveNaturalCoordinates() ) ) ) {
                element->giveIPValue(dam, srcgp, IST_PrincipalDamageTensor, tStep);
                if ( this->stateFilter ) {
                    // consider only points with same state
                    if ( ( ( state == 1 ) && ( norm(dam) > 1.e-3 ) ) || ( ( ( state == 0 ) && norm(dam) < 1.e-3 ) ) ) {
                        npoints++;
                        dist.at(npoints) = coords.distance(srcgpcoords);
                        gpList [ npoints - 1 ] = srcgp;
                    }
                } else {
                    // take all points into account
                    npoints++;
                    dist.at(npoints) = coords.distance(srcgpcoords);
                    gpList [ npoints - 1 ] = srcgp;
                }
            } else {
                OOFEM_ERROR("computeGlobalCoordinates failed");
            }
        }
    }

    if ( npoints != actualNumberOfPoints ) {
        OOFEM_ERROR("internal error");
    }

    //minNumberOfPoints = min (actualNumberOfPoints, minNumberOfPoints+2);

    patchGPList.clear();
    // now find the minNumberOfPoints with smallest distance
    // from point of interest
    double swap, minDist;
    int minDistIndx = 0;
    // loop over all points
    for ( int i = 1; i <= minNumberOfPoints; i++ ) {
        minDist = dist.at(i);
        minDistIndx = i;
        // search for point with i-th smallest distance
        for ( j = i + 1; j <= actualNumberOfPoints; j++ ) {
            if ( dist.at(j) < minDist ) {
                minDist = dist.at(j);
                minDistIndx = j;
            }
        }

        // remember this ip
        patchGPList.push_front(gpList [ minDistIndx - 1 ]);
        swap = dist.at(i);
        dist.at(i) = dist.at(minDistIndx);
        dist.at(minDistIndx) = swap;
        srcgp = gpList [ i - 1 ];
        gpList [ i - 1 ] = gpList [ minDistIndx - 1 ];
        gpList [ minDistIndx - 1 ] = srcgp;
    }

    if ( patchGPList.size() != minNumberOfPoints ) {
        OOFEM_ERROR("internal error 2");
        exit(1);
    }

    free(gpList);

#else

    // take all neighbors
    patchGPList.clear();
    for ( int ielem = 1; ielem <= patchList.giveSize(); ielem++ ) {
        element = patchDomain->giveElement( patchList.at(ielem) );
        iRule = element->giveDefaultIntegrationRulePtr();
        for ( GaussPoint *gp: *iRule ) {
            patchGPList.push_front( gp );
        }
    }

#endif
}


void
MMALeastSquareProjection :: finish(TimeStep *tStep)
{ }


int
MMALeastSquareProjection :: __mapVariable(FloatArray &answer, FloatArray &targetCoords,
                                          InternalStateType type, TimeStep *tStep)
{
    //int nelem, ielem,
    int neq = this->giveNumberOfUnknownPolynomialCoefficients(this->patchType);
    int nval = 0;
    FloatArray ipVal, coords, P;
    FloatMatrix a, rhs, x;
    Element *element;
    //IntegrationRule* iRule;

    a.resize(neq, neq);
    a.zero();

    // determine the value from patch
    int size = patchGPList.size();
    if ( size == 1 ) {
        GaussPoint *srcgp  = *patchGPList.begin();
        srcgp->giveElement()->giveIPValue(answer, srcgp, type, tStep);
    } else if ( size < neq ) {
        OOFEM_ERROR("internal error");
    } else {
        for ( auto &srcgp: patchGPList ) {
            element = srcgp->giveElement();
            element->giveIPValue(ipVal, srcgp, type, tStep);
            if ( nval == 0 ) {
                nval = ipVal.giveSize();
                rhs.resize(neq, nval);
                rhs.zero();
            }
            if ( element->computeGlobalCoordinates( coords, srcgp->giveNaturalCoordinates() ) ) {
                coords.subtract(targetCoords);
                // compute ip contribution
                this->computePolynomialTerms(P, coords, patchType);
                for ( int j = 1; j <= neq; j++ ) {
                    for ( int k = 1; k <= nval; k++ ) {
                        rhs.at(j, k) += P.at(j) * ipVal.at(k);
                    }

                    for ( int k = 1; k <= neq; k++ ) {
                        a.at(j, k) += P.at(j) * P.at(k);
                    }
                }
            } else {
                OOFEM_ERROR("computeGlobalCoordinates failed");
            }
        }

#if 0
        // check for correct solution
        FloatMatrix aa = a;
#endif
        a.solveForRhs(rhs, x);
#if 0
        // check for correct solution
        FloatMatrix tmp;
        tmp.beProductOf(aa, x);
        for ( j = 1; j <= neq; j++ ) {
            for ( k = 1; k <= nval; k++ ) {
                if ( fabs( tmp.at(j, k) - rhs.at(j, k) ) > 1.e-3 ) {
                    printf("(SE)");
                }
            }
        }

#endif
        // determine the value from patch
        targetCoords.zero();
        //gp->giveElement()->computeGlobalCoordinates (coords, gp->giveNaturalCoordinates());
        this->computePolynomialTerms(P, targetCoords, patchType);

        answer.resize(nval);
        answer.zero();
        for ( int i = 1; i <= nval; i++ ) {
            for ( int j = 1; j <= neq; j++ ) {
                answer.at(i) += P.at(j) * x.at(j, i);
            }
        }
    }

    /*
     * double ee;
     * aa.subtract (answer);
     * ee = sqrt(dotProduct(aa,aa,aa.giveSize()));
     * if ((aa.giveSize() != answer.giveSize()) || (ee > 1.e-6)) {
     * printf("Diference @@@@@!\n");
     * exit (1);
     * }
     */
    return 1;
}

int
MMALeastSquareProjection :: mapStatus(MaterialStatus &oStatus) const
{
    OOFEM_ERROR("not implemented yet.")

    return 0;
}

void
MMALeastSquareProjection :: computePolynomialTerms(FloatArray &P, FloatArray &coords, MMALeastSquareProjectionPatchType type)
{
    if ( type == MMALSPPatchType_2dq ) {
        /*
         * P.resize (1);
         * P.at(1) = 1.0;
         */
        /*
         * P.resize (3);
         * P.at(1) = coords.at(1);
         * P.at(2) = coords.at(2);
         * P.at(3) = 1.0;
         */

        P.resize(6);
        P.at(1) = 1.0;
        P.at(2) = coords.at(1);
        P.at(3) = coords.at(2);
        P.at(4) = coords.at(1) * coords.at(2);
        P.at(5) = coords.at(1) * coords.at(1);
        P.at(6) = coords.at(2) * coords.at(2);
    } else if ( type == MMALSPPatchType_1dq ) {
        P.resize(3);
        P.at(1) = coords.at(1) * coords.at(1);
        P.at(2) = coords.at(1);
        P.at(3) = 1.0;
    } else {
        OOFEM_ERROR("unknown regionType");
    }
}

int
MMALeastSquareProjection :: giveNumberOfUnknownPolynomialCoefficients(MMALeastSquareProjectionPatchType regType)
{
    if ( regType == MMALSPPatchType_2dq ) {
        return 6;
    } else if ( regType == MMALSPPatchType_1dq ) {
        return 3;
    } else {
        return 0;
    }
}


IRResultType
MMALeastSquareProjection :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    this->stateFilter = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, this->stateFilter, _IFT_MMALeastSquareProjection_statefilter);

    this->regionFilter = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, this->regionFilter, _IFT_MMALeastSquareProjection_regionfilter);

    return IRRT_OK;
}


void MMALeastSquareProjection :: giveInputRecord(DynamicInputRecord &input)
{
    if ( this->stateFilter ) {
        input.setField(this->stateFilter, _IFT_MMALeastSquareProjection_statefilter);
    }
    if ( this->regionFilter ) {
        input.setField(this->stateFilter, _IFT_MMALeastSquareProjection_regionfilter);
    }
}
} // end namespace oofem
