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

#include "fastmarchingmethod.h"
#include "domain.h"
#include "mathfem.h"
#include "node.h"
#include "element.h"
#include "conTable.h"

#ifndef __MAKEDEPEND
 #include <cstdlib>
#endif

namespace oofem {
void
FastMarchingMethod :: solve(FloatArray &dmanValues,
                            const std :: list< int > &bcDofMans,
                            double F)
{
    int i, j, jn, candidate;
    const IntArray *neighborElemList;
    ConnectivityTable *ct = domain->giveConnectivityTable();
    Element *ie;


    this->dmanValuesPtr = & dmanValues;
    // tag points with boundary value as known
    // then tag as trial all points that are one grid point away
    // finally tag as far all other grid points
    this->initialize(dmanValues, bcDofMans, F);

    // let candidate be the thrial point with smallest T value
    while ( ( candidate = this->getSmallestTrialDofMan() ) ) {
        // add the candidate to known, remove it from trial
        dmanRecords.at(candidate - 1).status = FMM_Status_KNOWN;
        // tag as trial all neighbors of candidate that are not known
        // if the neighbor is in far, remove and add it to the trial set
        // and recompute the values of T at all trial neighbors of candidate
        neighborElemList = ct->giveDofManConnectivityArray(candidate);
        for ( i = 1; i <= neighborElemList->giveSize(); i++ ) {
            ie = domain->giveElement( neighborElemList->at(i) );
            for ( j = 1; j <= ie->giveNumberOfDofManagers(); j++ ) {
                jn = ie->giveDofManagerNumber(j);
                if ( dmanRecords.at(jn - 1).status != FMM_Status_KNOWN ) {
                    // recompute the value of T at candidate trial neighbor
                    this->updateTrialValue(dmanValues, jn, F);
                }
            }
        }
    }
}



void
FastMarchingMethod :: initialize(FloatArray &dmanValues,
                                 const std :: list< int > &bcDofMans,
                                 double F)
{
    // tag points with boundary value as known
    // then tag as trial all points that are one grid point away
    // finally tag as far all other grid points
    int i, k, l;
    int jnode, neighborNode, nnode = domain->giveNumberOfDofManagers();
    Element *neighborElem;
    const IntArray *neighborElemList;
    ConnectivityTable *ct = domain->giveConnectivityTable();

    // all points are far by default
    dmanRecords.resize(nnode);
    for ( i = 0; i < nnode; i++ ) {
        dmanRecords [ i ].status = FMM_Status_FAR;
    }

    // first tag all boundary points
    std :: list< int > :: const_iterator it;

    for ( it = bcDofMans.begin(); it != bcDofMans.end(); ++it ) {
        if ( ( jnode = * it ) > 0 ) {
            dmanRecords.at(jnode - 1).status = FMM_Status_KNOWN;
        } else {
            dmanRecords.at(-jnode - 1).status = FMM_Status_KNOWN_BOUNDARY;
        }
    }

    // tag as trial all points that are one grid point away
    for ( it = bcDofMans.begin(); it != bcDofMans.end(); ++it ) {
        jnode = abs(* it);

        neighborElemList = ct->giveDofManConnectivityArray(jnode);
        for ( k = 1; k <= neighborElemList->giveSize(); k++ ) {
            if ( neighborElemList->at(k) == i ) {
                continue;
            }

            neighborElem = domain->giveElement( neighborElemList->at(k) );
            for ( l = 1; l <= neighborElem->giveNumberOfDofManagers(); l++ ) {
                neighborNode = neighborElem->giveDofManagerNumber(l);
                if ( ( dmanRecords.at(neighborNode - 1).status != FMM_Status_KNOWN ) &&
                    ( dmanRecords.at(neighborNode - 1).status != FMM_Status_KNOWN_BOUNDARY ) &&
                    ( dmanRecords.at(neighborNode - 1).status != FMM_Status_TRIAL ) ) {
                    this->updateTrialValue(dmanValues, neighborNode, F);
                }
            }
        }
    }
}


void
FastMarchingMethod :: updateTrialValue(FloatArray &dmanValues, int id, double F)
{
    int i, j, ai, bi, ci, h, nroot, _ind = 0;
    double at, bt, ht, a, b, u, cos_fi, sin_fi, _a, _b, _c, r1, r2, r3, t = 0.0, _h;
    bool reg_upd_flag;
    Element *ie;
    FloatArray *ac, *bc, *cc, cb, ca;
    const IntArray *neighborElemList;
    ConnectivityTable *ct = domain->giveConnectivityTable();


    // first look for suitable element that can produce admissible value
    // algorithm limited to non-obtuse 2d triangulations
    neighborElemList = ct->giveDofManConnectivityArray(id);
    for ( i = 1; i <= neighborElemList->giveSize(); i++ ) {
        // test element if admissible
        ie = domain->giveElement( neighborElemList->at(i) );
        if ( ie->giveGeometryType() == EGT_triangle_1 ) {
            for ( j = 1; j <= 3; j++ ) {
                if ( ie->giveDofManagerNumber(j) == id ) {
                    _ind = j;
                }
            }

            ci = ie->giveDofManagerNumber(_ind);
            ai = ie->giveDofManagerNumber(1 + ( _ind ) % 3);
            bi = ie->giveDofManagerNumber(1 + ( _ind + 1 ) % 3);

            if ( ( dmanRecords.at(ai - 1).status == FMM_Status_KNOWN ) &&
                ( dmanRecords.at(bi - 1).status == FMM_Status_KNOWN ) ) {
                at = dmanValues.at(ai);
                bt = dmanValues.at(bi);
                if ( fabs(at) > fabs(bt) ) {
                    h = ai;
                    ai = bi;
                    bi = h;
                    ht = at;
                    at = bt;
                    bt = ht;
                }

                // get nodal coordinates
                ac = domain->giveNode(ai)->giveCoordinates();
                bc = domain->giveNode(bi)->giveCoordinates();
                cc = domain->giveNode(ci)->giveCoordinates();

                // a = distance of BC
                a = cc->distance(bc);
                // b = distance of AC
                b = cc->distance(ac);
                // compute fi angle
                cb.beDifferenceOf(*bc, *cc);
                cb.normalize();
                ca.beDifferenceOf(*ac, *cc);
                ca.normalize();
                cos_fi = cb.dotProduct(ca);
                sin_fi = sqrt(1.0 - cos_fi * cos_fi);
                u = fabs(bt);
                // compute quadratic equation coefficients for t
                _a = ( a * a + b * b - 2.0 * a * b * cos_fi );
                _b = 2.0 * b * u * ( a * cos_fi - b );
                _c = b * b * ( u * u - F * F * a * a * sin_fi * sin_fi );
                cubic3r(0.0, _a, _b, _c, & r1, & r2, & r3, & nroot);
                reg_upd_flag = true;

                if ( nroot == 0 ) {
                    reg_upd_flag = false;
                } else if ( nroot == 1 ) {
                    t = r1;
                } else if ( r1 >= 0.0 ) {
                    t = r1;
                } else if ( r2 >= 0.0 ) {
                    t = r2;
                } else {
                    reg_upd_flag = false;
                }

                if ( reg_upd_flag ) {
                    _h = b * ( t - u ) / t;
                    if ( ( t > u ) && ( _h > a * cos_fi ) && ( _h < a / cos_fi ) ) {
                        if ( dmanRecords.at(ci - 1).status == FMM_Status_FAR ) {
                            dmanValues.at(ci) = sgn(F) * t + at;
                        } else if ( F > 0. ) {
                            dmanValues.at(ci) = min(dmanValues.at(ci), sgn(F) * t + at);
                        } else {
                            dmanValues.at(ci) = max(dmanValues.at(ci), sgn(F) * t + at);
                        }
                    } else {
                        reg_upd_flag = false;
                    }
                }

                if ( !reg_upd_flag ) {
                    if ( F > 0. ) {
                        _h = min(b * F + at, a * F + bt);
                    } else {
                        _h = max(b * F + at, a * F + bt);
                    }

                    if ( dmanRecords.at(ci - 1).status == FMM_Status_FAR ) {
                        dmanValues.at(ci) = _h;
                    } else if ( F > 0. ) {
                        dmanValues.at(ci) = min(dmanValues.at(ci), _h);
                    } else {
                        dmanValues.at(ci) = max(dmanValues.at(ci), _h);
                    }
                }

                // if not yet in queue (trial for the first time), put it there
                if ( dmanRecords.at(ci - 1).status != FMM_Status_TRIAL ) {
                    dmanTrialQueue.push(ci);
                }

                dmanRecords.at(ci - 1).status = FMM_Status_TRIAL;
            } // admissible triangle

        } // end EGT_triangle_1 element type

    }
}

int
FastMarchingMethod :: getSmallestTrialDofMan()
{
    int answer;
    if ( dmanTrialQueue.empty() ) {
        return 0;
    }

    answer = dmanTrialQueue.top();
    dmanTrialQueue.pop();
    return answer;
}
} // end namespace oofem
