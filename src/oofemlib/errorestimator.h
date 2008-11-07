/* $Header: /home/cvs/bp/oofem/oofemlib/src/errorestimator.h,v 1.10.4.1 2004/04/05 15:19:43 bp Exp $ */
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


//   *****************************************
//   *** CLASS ERROR ESTIMATOR (INDICATOR) ***
//   *****************************************

#ifndef errorestimator_h
#define errorestimator_h

#include "femcmpnn.h"
#include "compiler.h"

#include "interface.h"
#include "remeshingcrit.h"
#include "classtype.h"
#include "errorestimatortype.h"

class Domain;
class Element;
class TimeStep;

/** Type characterizing different type of errors*/
enum EE_ValueType { globalNormEEV, globalErrorEEV, globalWeightedErrorEEV };
/** Type characterizcing different type of element errors*/
enum EE_ErrorType { unknownET, indicatorET, internalStressET, primaryUnknownET };
/** Type determining whether temp or nontemp variables are used for error evaluation */
enum EE_ErrorMode { equilibratedEM, temporaryEM };


/**
 * The base class for all error estimation or error indicator algorithms.
 * The basic task is to evaluate the error on associated domain. If this task requires
 * the special element algorithms, these should be included using interface concept.
 *
 * This estimator should also provide the compatible Remeshing Criteria class, which
 * based on warious error measures will evaluate the required mesh density of a new domain.
 */
class ErrorEstimator : public FEMComponent
{
protected:

    ErrorEstimatorType eeType;
    RemeshingCriteria *rc;
    /** map indicating regions to skip (region - cross section model).
     * Do not access this variable directly, since this variable is read from input and could have size different
     * from actual number of regions - use alvays the skipRegion method, since it performs size check.  and handles
     * all regions correctly.
     */
    IntArray regionSkipMap;
    ///number of skipped elements
    int skippedNelems;

public:
    /// Constructor
    ErrorEstimator(int n, Domain *d) : FEMComponent(n, d) { rc = NULL;
                                                            skippedNelems = 0;
                                                            regionSkipMap.resize(0); }
    /// Destructor
    ~ErrorEstimator() { if ( rc ) { delete rc; } }
    /// Sets Domain; should also re-initialize attributes if necessary.
    void setDomain(Domain *d);
    /** Returns the element error. The estimateError service should be called before.
     * @param type error type
     * @param elem element for which error requested
     * @param tStep time step
     */
    virtual double giveElementError(EE_ErrorType type, Element *elem, TimeStep *tStep) = 0;
    /** Returns the characteristic value of given type. The estimateError service should be called before.
     * This method is supposed to be used by associated remeshingCriteria to acces some characteristic
     * values already computed or known at error estinator level.
     * @param type error type
     * @param tStep time step
     */
    virtual double giveValue(EE_ValueType type, TimeStep *tStep) = 0;
    /** Returns number of elements skipped in error estimation
     */
    int giveNumberOfSkippedElements() { return skippedNelems; }
    /**
     * Estimates the error on associated domain at given timeSte. The estimated values can be
     * requested using giveElementError and giveValue methods. The type of errors provided
     * depends on error estimator type implementing the service.
     * @param mode erorr mode
     * @param tStep time step
     */
    virtual int estimateError(EE_ErrorMode mode, TimeStep *tStep) = 0;
    /** Returns reference to associated remeshing criteria.
     */
    virtual RemeshingCriteria *giveRemeshingCrit() = 0;
    /** Initializes receiver acording to object description stored in input reader.
     * This function is called immediately after creating object using
     * constructor. InitString can be imagined as data record in component database
     * belonging to receiver. Receiver may use value-name extracting functions
     * to extract particular field from record.
     * @see readInteger, readDouble and similar functions */
    virtual IRResultType initializeFrom(InputRecord *ir);
    /// Returns class name of the receiver.
    const char *giveClassName() const { return "ErrorEstimator"; }
    /** Returns classType id of receiver.
     * @see FEMComponent::giveClassID
     */
    classType                giveClassID() const { return ErrorEstimatorClass; }
    /** Returns error estimation type of receiver.
     */
    const ErrorEstimatorType giveErrorEstimatorType() { return eeType; }

    /** Returns nonzero if region has been skipped in error estimation (user option)
     * It is strongly recomended to use this function, instead of direc access to
     * regionSkipMap variable by derivad classes, since the size check is done here.
     */
    int skipRegion(int reg) { if ( reg <= regionSkipMap.giveSize() ) { return regionSkipMap.at(reg); } else { return 0; } }
    virtual void reinitialize() {this->rc->reinitialize();}

protected:
};


#endif // errorestimator_h






