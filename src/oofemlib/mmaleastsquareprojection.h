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

#ifndef mmaleastsquareprojection_h
#define mmaleastsquareprojection_h

#include "materialmappingalgorithm.h"
#include "interface.h"

#include <list>

///@name Input fields for MMALeastSquareProjection
//@{
#define _IFT_MMALeastSquareProjection_statefilter "mmalsp_statefilter"
#define _IFT_MMALeastSquareProjection_regionfilter "mmalsp_regionfilter"
//@}

namespace oofem {
class Domain;
class Element;
class TimeStep;
class DynamicInputRecord;

enum MMALeastSquareProjectionPatchType { MMALSPPatchType_1dq, MMALSPPatchType_2dq };
/*
 * Defines, whether only the necessary number of closest points will be used
 * to fit a polynomial. If not defined, more points from gp neighborhood will
 * be used, based on element connectivity.
 */
//#define MMALSP_ONLY_CLOSEST_POINTS

/**
 * The class implements the transfer of state variables based on
 * Least square fit over old mesh integration points (in the neighborhood of point of interest).
 *
 * The algorithm of projecting internal vars (q) can be summarized as follows:
 * -# The "source" element on old mesh containing point of interest is found
 * -# The element patch is constructed from "source" element neighbourhood
 * -# The least square fit is done
 * -# Value in point of interest is evaluated.
 *
 * It is obvious, that this mapper operates locally and therefore there is no need to declared this
 * mapper as static material member.
 *
 */
class OOFEM_EXPORT MMALeastSquareProjection : public MaterialMappingAlgorithm
{
protected:
    /// If set, then only IP in the neighbourhood with same state can be used to interpolate the values.
    int stateFilter;
    /// If set, then only IP in the same region are taken into account.
    int regionFilter;
    /// List of Gp participating in patch.
    std :: list< GaussPoint * >patchGPList;
    /// Patch domain.
    Domain *patchDomain;
    /// Type of patch.
    MMALeastSquareProjectionPatchType patchType;
public:
    /// Constructor
    MMALeastSquareProjection();
    /// Destructor
    virtual ~MMALeastSquareProjection();

    virtual void __init(Domain *dold, IntArray &type, FloatArray &coords, Set &sourceElemSet, TimeStep *tStep, bool iCohesiveZoneGP = false);

    virtual void finish(TimeStep *tStep);

    virtual int __mapVariable(FloatArray &answer, FloatArray &coords, InternalStateType type, TimeStep *tStep);

    virtual int mapStatus(MaterialStatus &oStatus) const;

    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual void giveInputRecord(DynamicInputRecord &input);

    virtual const char *giveClassName() const { return "MMALeastSquareProjectionPatchType"; }

protected:
    void computePolynomialTerms(FloatArray &P, FloatArray &coords, MMALeastSquareProjectionPatchType type);
    int giveNumberOfUnknownPolynomialCoefficients(MMALeastSquareProjectionPatchType regType);
};
} // end namespace oofem
#endif // mmaleastsquareprojection_h
