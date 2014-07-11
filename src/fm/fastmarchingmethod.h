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

#ifndef fastmarchingmethod_h
#define fastmarchingmethod_h

#include "floatarray.h"
#include "mathfem.h"

#include <vector>
#include <list>
#include <queue>

namespace oofem {
class Domain;

/**
 * Fast Marching Method for unstructured grids.
 * Used to solve Eikonal equation and especially to construct
 * signed distance function.
 */
class FastMarchingMethod
{
protected:

    /// Type describing node status for fast marching method.
    enum FNM_Status_Type {
        FMM_Status_FAR,            ///< Nodes not yet visited.
        FMM_Status_TRIAL,          ///< Trial nodes, candidates for known (accepted).
        FMM_Status_KNOWN,          ///< Accepted nodes.
        FMM_Status_KNOWN_BOUNDARY, ///< Boundary nodes, from which the front will not propagate.
    };
    /// DofManager Fast Marching data record.
    class FMM_DofmanRecord
    {
public:
        FNM_Status_Type status;
    };

    /// Array of DofManager records.
    std :: vector< FMM_DofmanRecord >dmanRecords;
    /// Pointer to working set of dmanValues.
    const FloatArray *dmanValuesPtr;
    // Delegate of FMM_DofmanRecord; stored in priority queue.
    //class FMM_DofmanRecordDelegate {
    //public:FMM_DofmanRecordDelegaFMM_DofmanRecordDelegatete
    //  int id;
    //};

    struct FMM_DofmanRecordDelegate_greater {
        const FloatArray **dmanValuesPtrRef;
        FMM_DofmanRecordDelegate_greater(const FloatArray **dv) : dmanValuesPtrRef(dv)
        { }
        bool operator() (const int & p, const int & q) const
        { return ( fabs( ( * dmanValuesPtrRef )->at(p) ) > fabs( ( * dmanValuesPtrRef )->at(q) ) ); }
    };

    /// Domain.
    Domain *domain;

    /// Priority queue for trial T values.
    std :: priority_queue< int, std :: vector< int > , FMM_DofmanRecordDelegate_greater >dmanTrialQueue;

public:
    /**
     * Constructor. Takes two two arguments. Creates
     * FastMarchingMethod material interface instance with given number and belonging to given domain.
     * @param d Domain to which component belongs to.
     */
    FastMarchingMethod(Domain * d) : dmanTrialQueue( FMM_DofmanRecordDelegate_greater(& this->dmanValuesPtr) ) {
        domain = d;
    }
    ~FastMarchingMethod() { }

    /**
     * Solution of problem.
     * @param[in,out] dmanValues On input should contain boundary
     * values for those dofnam that are known; on output will contain solution.
     * @param bcDofMans A list containing IDs (numbers) of those
     * dofmans, for which boundary value is known. If this number is positive,
     * then the front will propagate from this dofman, if negative, then the front
     * will not propagate from this dofman (usefull, when one needs to construct
     * "one sided" solution).
     * @param F is the front propagation speed.
     */
    void solve(FloatArray &dmanValues, const std :: list< int > &bcDofMans, double F);

    // identification
    const char *giveClassName() const { return "FastMarchingMethod"; }

protected:
    /// Initialize receiver.
    void initialize(FloatArray &dmanValues, const std :: list< int > &bcDofMans, double F);

    /// Updates the distance of trial node with given id).
    void updateTrialValue(FloatArray &dmanValues, int id, double F);

    /// Get the trial point with smallest T; zero if empty.
    int  getSmallestTrialDofMan();
};
} // end namespace oofem
#endif // fastmarchingmethod_h
