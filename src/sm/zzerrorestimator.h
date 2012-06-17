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

#ifndef zzerrorestimator_h
#define zzerrorestimator_h

#include "errorestimator.h"
#include "interface.h"
#include "internalstatetype.h"
#include "flotarry.h"
#include "statecountertype.h"

namespace oofem {
#define ZZErrorEstimator_ElementResultCashed

class Element;
class GaussPoint;

/**
 * The implementation of Zienkiewicz Zhu Error Estimator (Zienkiewicz and Zhu: A simple error
 * estimator and adaptive procedure for practical engineering analysis, International Journal
 * for Numerical Methods in Engineering, vol. 24, 337-357, 1987).
 * The basic task is to evaluate the stress error on associated domain.
 * The algorithm is written in general way, so it is possible to to evaluate
 * different errors (for example temperature error). Then corresponding
 * member attribute identifying the type of quantity used should be declared and initialized
 * (for example in instanciateYourself() service). Then the modification is required
 * only when requesting element contributions.
 *
 * This task requires the special element algorithms, which are supported at element level
 * using interface concept.
 * This estimator also provides the compatible Remeshing Criteria, which
 * based on error measure will evaluate the required mesh density of a new domain.
 */
class ZZErrorEstimator : public ErrorEstimator
{
public:
    /// Type of norm used.
    enum NormType { L2Norm, EnergyNorm };
    /// Nodal recovery type.
    enum NodalRecoveryType { ZZRecovery, SPRRecovery };

protected:
    /// Global error norm.
    double globalENorm;
    /// Global norm of quantity which error is evaluated.
    double globalSNorm;
#ifdef ZZErrorEstimator_ElementResultCashed
    /// Cache storing element norms.
    FloatArray eNorms;
#endif
    /// Type of norm used.
    NormType normType;
    /// Nodal recovery type.
    NodalRecoveryType nodalRecoveryType;

    /// Actual state counter.
    StateCounterType stateCounter;

public:
    /// Constructor
    ZZErrorEstimator(int n, Domain *d) : ErrorEstimator(n, d) {
        eeType = EET_ZZEE;
        stateCounter = 0;
        normType = L2Norm;
        nodalRecoveryType = ZZRecovery;
    }
    /// Destructor
    virtual ~ZZErrorEstimator() { }

    virtual double giveElementError(EE_ErrorType type, Element *elem, TimeStep *tStep);
    virtual double giveValue(EE_ValueType type, TimeStep *tStep);

    virtual int estimateError(EE_ErrorMode mode, TimeStep *tStep);
    virtual RemeshingCriteria *giveRemeshingCrit();

    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual const char *giveClassName() const { return "ZZErrorEstimator"; }
    virtual classType giveClassID() const { return ZZErrorEstimatorClass; }
};


/**
 * The element interface corresponding to ZZErrorEstimator.
 * It declares necessary services provided by element to be compatible with ZZErrorEstimator.
 */
class ZZErrorEstimatorInterface : public Interface
{
public:
    /// Constructor
    ZZErrorEstimatorInterface() { }

    /// Returns reference to corresponding element
    virtual Element *ZZErrorEstimatorI_giveElement() = 0;
    /// Computes the interpolation matrix used for recovered values
    virtual void ZZErrorEstimatorI_computeEstimatedStressInterpolationMtrx(FloatArray &answer, GaussPoint *gp,
                                                                           InternalStateType type) = 0;
    /**
     * Computes the element contributions to global norms.
     * @param eNorm Element contribution to error norm.
     * @param sNorm Element contribution to norm of quantity which error is evaluated.
     * @param norm Determines the type of norm to evaluate.
     * @param type Determines the type of internal state to compute for.
     * @param tStep Time step.
     */
    virtual void ZZErrorEstimatorI_computeElementContributions(double &eNorm, double &sNorm, ZZErrorEstimator :: NormType norm,
                                                               InternalStateType type, TimeStep *tStep);
};


/**
 * The class representing Zienkiewicz-Zhu remeshing criteria.
 * (Assumes that error is equally distributed between elements, then the requirement for max. permissible error
 * can be translated into placing a limit on the error on each element.)
 * The basic task is to evaluate the required mesh density (at nodes) on given domain,
 * based on informations provided by the compatible error estimator.
 *
 * The remeshing criteria is maintained by the corresponding error estimator. This is mainly due to fact, that is
 * necessary for given EE to create compatible RC. In our concept, the EE is responsible.
 */
class ZZRemeshingCriteria : public RemeshingCriteria
{
public:
    /// Mode of receiver, allows to use it in more general situations.
    enum ZZRemeshingCriteriaModeType { stressBased };

protected:
    /// Array of nodal mesh densities.
    FloatArray nodalDensities;
    /// Remeshing strategy proposed.
    RemeshingStrategy remeshingStrategy;
    /// Actual values (densities) state counter.
    StateCounterType stateCounter;
    /// Mode of receiver.
    ZZRemeshingCriteriaModeType mode;
    /// Required error to obtain.
    double requiredError;
    /// Minimum element size allowed.
    double minElemSize;

public:
    /// Constructor
    ZZRemeshingCriteria(int n, ErrorEstimator *e);
    /// Destructor
    virtual ~ZZRemeshingCriteria() { }


    virtual double giveRequiredDofManDensity(int num, TimeStep *tStep, int relative = 0);
    virtual double giveDofManDensity(int num);
    virtual RemeshingStrategy giveRemeshingStrategy(TimeStep *tStep);
    virtual int estimateMeshDensities(TimeStep *tStep);
    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual const char *giveClassName() const { return "ZZErrorEstimator"; }
    virtual classType giveClassID() const { return ZZRemeshingCriteriaClass; }
};

/**
 * The corresponding element interface to ZZRemeshingCriteria class.
 * Declares the necessary services, which have to be provided by particular elements.
 */

class ZZRemeshingCriteriaInterface : public Interface
{
public:
    /// Constructor
    ZZRemeshingCriteriaInterface() : Interface() { }
    /**
     * Determines the characteristic size of element. This quantity is defined as follows:
     * For 1D it is the element length, for 2D it is the square root of element area.
     */
    virtual double ZZRemeshingCriteriaI_giveCharacteristicSize() = 0;
    /**
     * Returns the polynomial order of receiver trial functions.
     */
    virtual int ZZRemeshingCriteriaI_givePolynOrder() = 0;
};
} // end namespace oofem
#endif // zzerrorestimator_h
