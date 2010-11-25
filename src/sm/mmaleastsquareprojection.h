/* $Header: /home/cvs/bp/oofem/sm/src/mmaleastsquareprojection.h,v 1.6 2003/05/19 13:04:00 bp Exp $ */
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
 *               Copyright (C) 1993 - 2008   Borek Patzak
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



//   *********************************************************
//   *** CLASS SHAPE FUNCTION PROJECTION MAPPING ALGORITHM ***
//   *********************************************************

#ifndef mmaleastsquareprojection_h
#define mmaleastsquareprojection_h

#include "alist.h"
#include "materialmappingalgorithm.h"
#include "interface.h"
#include "compiler.h"

#include "dynalist.h"

namespace oofem {
class Domain;
class Element;
class TimeStep;

enum MMALeastSquareProjectionPatchType { MMALSPPatchType_1dq, MMALSPPatchType_2dq };
/**
 * Defines, whether only the necessary number of colosest points will be used
 * to fit a polynomial. If not defined, more points from gp neighborhood will
 * be used, based on element conectivity.
 */
//#define MMALSP_ONLY_CLOSEST_POINTS

/**
 * The class implements the transfer of state variables based on
 * Least square fit over old mesh integration points (in the neighborhood of point of interest).
 *
 * The algorithm of projecting internal vars (q) can be summarized as follows:
 * 1) The "source" element on old mesh containing point of interest is found
 * 2) The element patch is constructed from "source" element neighbourhood
 * 3) The least square fit is done
 * 4) Value in point of interest is evaluated.
 *
 * It is obvious, that this mapper operates locally and therefore there is no need to declared this
 * mapper as static material member.
 *
 */
class MMALeastSquareProjection : public MaterialMappingAlgorithm
{
protected:
    /// if set, then only IP in the neighbourhood with same state can be used to interpolate the values.
    int stateFilter;
    /// if set, then only IP in the same region are taken into account.
    int regionFilter;
    /// List of Gp participating in patch
    dynaList< GaussPoint * >patchGPList;
    /// patch domain
    Domain *patchDomain;
    /// Type of patch
    MMALeastSquareProjectionPatchType patchType;
public:
    /** Constructor
     * @param oldd old mesh reference
     */
    MMALeastSquareProjection();
    /// Destructor
    ~MMALeastSquareProjection();
    /**
     * Initializes the receiver state before mapping. The idea is to place some
     * comon global oprerations before mapping particular IP's if necessary.
     * Stores Times stamp of last initialization, so multiple calls for same step
     * do not initialize receiver again.
     * @param dold old domain
     * @param varTypes array of InternalStateType values, identifiing all vars to be mapped
     * @param coords coordinates of the receiver point
     * @param region if > 0 region id of receiver point,, if < 0 ignore regions.
     * @param tStep time step
     */
    virtual void __init(Domain *dold, IntArray &type, FloatArray &coords, int region, TimeStep *tStep);
    /**
     * Finishes the mapping for given time step. Used to perform cleanup.
     * Typically some mappers reguire to compute some global mesh data related to
     * current step, which are valid for all IPs - so they are computed only once for
     * all IPs.
     */
    virtual void finish(TimeStep *tStep);
    /** Maps and update the unknown of given type from
     * old mesh oldd to new mesh to which gp belongs to. The result is stored in answer array.
     * The closest IP point to specified one is found and its state variable is
     * returned as a result.
     * @param answer contains result
     * @param type determines the type of internal variable
     * @param coords coordinates of receiver to which mapping occur
     * @param tStep time step
     * @return nonzero if o.k.
     */
    virtual int __mapVariable(FloatArray &answer, FloatArray &coords, InternalStateType type, TimeStep *tStep);
    /** Initializes receiver acording to object description stored in input record.
     * InitString can be imagined as data record in component database
     * belonging to receiver. Receiver may use value-name extracting functions
     * to extract particular field from record.
     * @see readInteger, readDouble and similar functions */
    virtual IRResultType initializeFrom(InputRecord *ir);
    /** Setups the input record string of receiver
     * @param str string to be filled by input record
     * @param keyword print record keyword (default true)
     */
    virtual int giveInputRecordString(std :: string &str, bool keyword = true);
    /** Returns class name of the receiver */
    const char *giveClassName() const { return "MMALeastSquareProjectionPatchType"; }



protected:
    void computePolynomialTerms(FloatArray &P, FloatArray &coords, MMALeastSquareProjectionPatchType type);
    int  giveNumberOfUnknownPolynomialCoefficients(MMALeastSquareProjectionPatchType regType);
};
} // end namespace oofem
#endif // mmaleastsquareprojection_h
