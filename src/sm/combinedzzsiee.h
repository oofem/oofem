/* $Header: /home/cvs/bp/oofem/sm/src/combinedzzsiee.h,v 1.5 2003/04/06 14:08:30 bp Exp $ */
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

//   *************************************************************************************
//   *** CLASS COMBINED ZIENKIEWICZ AND ZHU ERROR ESTIMATOR AND SCALAR ERROR INDICATOR ***
//   *************************************************************************************

#ifndef combinedzzsiee_h
#define combinedzzsiee_h

#include "zzerrorestimator.h"
#include "scalarerrorindicator.h"
#include "directerrorindicatorrc.h"
/**
 * The implementation of combined criteria: Zienkiewicz Zhu Error Estimator for elastic regime and
 * scalar error indicator in non-linear regime.
 * The basic task is to evaluate the stress error on associated domain.
 * The algorithm is written in general way, so it is possible to to evaluate
 * different errors (for example temperature error). Then corresponding
 * member attribute identifying the type of quantity used should be declared and initialized
 * (for example in instanciateYourself() service). Then the modification is required
 * only when requesting element contributions.
 *
 */
class CombinedZZSIErrorEstimator : public ErrorEstimator
{
protected:
    ZZErrorEstimator zzee;
    ScalarErrorIndicator siee;
public:
    /// Constructor
    CombinedZZSIErrorEstimator(int n, Domain *d) : ErrorEstimator(n, d), zzee(n, d), siee(n, d) { eeType = EET_CZZSI; }
    /// Destructor
    ~CombinedZZSIErrorEstimator() { }
    /** Returns the element error of requested type. The estimateError service should be called before.
     * @param type error type
     * @param elem element for which error requested
     * @param tStep time step
     */
    virtual double giveElementError(EE_ErrorType type, Element *elem, TimeStep *tStep);
    /** Returns the characteristic value of given type.
     * The estimateError service should be called before. Intended to be used by remeshingCriterias to querry
     * various values provided by specific error estimator.
     * @param type value type
     * @param tStep time step
     */
    virtual double giveValue(EE_ValueType type, TimeStep *tStep);

    /**
     * Estimates the error on associated domain at given timeSte.
     * @param tStep time step
     */
    virtual int estimateError(EE_ErrorMode mode, TimeStep *tStep);
    /** Returns reference to associated remeshing criteria.
     */
    virtual RemeshingCriteria *giveRemeshingCrit();
    /** Initializes receiver acording to object description stored in input record.
     * This function is called immediately after creating object using
     * constructor. InitString can be imagined as data record in component database
     * belonging to receiver. Receiver may use value-name extracting functions
     * to extract particular field from record.
     * @see readInteger, readDouble and similar functions */
    virtual IRResultType initializeFrom(InputRecord *ir);
    /// Returns class name of the receiver.
    const char *giveClassName() const { return "CombinedZZSIErrorEstimator"; }
    /** Returns classType id of receiver.
     * @see FEMComponent::giveClassID
     */
    classType                giveClassID() const { return CombinedZZSIErrorEstimatorClass; }
    /// Sets Domain
    void setDomain(Domain *d);
protected:
};

/**
 * The class represent the corresponding remeshing criteria to CombinedZZSIErrorEstimator.
 * In regions, where the error indicator is larger than given treshold, the
 * mesh is refined according this indicator (currently linear interpolation).
 * Otherwise, the mesh size is determined from Zienkiewicz-Zhu remeshing criteria.
 * (Assumes that error is equally distributed between elements, then the requirement for max. permissible error
 * can be translated into placing a limit on the error on each element.)
 * The basic task is to evaluate the required mesh density (at nodes) on given domain,
 * based on informations provided by the compatible error ertimator.
 *
 *
 * The remeshing criteria is maintained by the corresponding error estimator. This is mainly due to fact, that is
 * necessary for given EE to create compatible RC. In our concept, the EE is responsible.
 *
 */
class CombinedZZSIRemeshingCriteria : public RemeshingCriteria
{
protected:

    ZZRemeshingCriteria zzrc;
    DirectErrorIndicatorRC dirc;
public:
    /// Constructor
    CombinedZZSIRemeshingCriteria(int n, ErrorEstimator *e);
    /// Destructor
    ~CombinedZZSIRemeshingCriteria() { }

    /** Returns the required mesh size n given dof manager.
     * The mesh density is defined as a required element size
     * (in 1D the element length, in 2D the square from element area).
     * @param num dofman  number
     * @param tStep time step
     * @param relative if zero, then actual density is returned, otherwise the relative density to current is returned.
     */
    virtual double giveRequiredDofManDensity(int num, TimeStep *tStep, int relative = 0);
    /**
     * Returns existing mesh size for given dof manager.
     * @param num dofMan number
     */
    virtual double giveDofManDensity(int num);
    /**
     * Determines, if the remeshing is needed, and if nedded, the type of strategy used
     */
    virtual RemeshingStrategy giveRemeshingStrategy(TimeStep *tStep);
    /**
     * Estimates the nodal densities.
     * @param tStep time step
     */
    virtual int estimateMeshDensities(TimeStep *tStep);
    /** Initializes receiver acording to object description stored in input string.
     * This function is called immediately after creating object using
     * constructor. InitString can be imagined as data record in component database
     * belonging to receiver. Receiver use value-name extracting functions
     * to extract particular field from record.
     * @see readInteger, readDouble and similar functions */
    virtual IRResultType initializeFrom(InputRecord *ir);

    /// Returns "ZZErrorEstimator" - class name of the receiver.
    const char *giveClassName() const { return "CombinedZZSIRemeshingCriteria"; }
    /** Returns ZZRemeshingCriteriaClass - classType id of receiver.
     * @see FEMComponent::giveClassID
     */
    classType                giveClassID() const { return CombinedZZSIRemeshingCriteriaClass; }
    /// Sets Domain
    void setDomain(Domain *d);

protected:
};

#endif // combinedzzsiee_h

