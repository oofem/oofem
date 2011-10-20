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
 *               Copyright (C) 1993 - 2011   Borek Patzak
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

#ifndef directerrorindicatorrc_h
#define directerrorindicatorrc_h

#include "compiler.h"

#include "interface.h"
#include "remeshingcrit.h"

namespace oofem {
class Domain;
class Element;
class TimeStep;

/**
 * The corresponding element interface to DirectErrorIndicatorRC class.
 * Declares the necessary services, which have to be provided by particular elements.
 */
class DirectErrorIndicatorRCInterface : public Interface
{
public:
    /// Constructor
    DirectErrorIndicatorRCInterface() : Interface() { }
    /**
     * Determines the characteristic size of element. This quantity is defined as follows:
     * For 1D it is the element length, for 2D it is the square root of element area.
     */
    virtual double DirectErrorIndicatorRCI_giveCharacteristicSize() = 0;
};


/**
 * The class is an implementation of "direct" remeshing criteria, which
 * maps the error indication, which is usually the value of observed internal variable
 * to the corresponding required element size. The error estimate and global error for
 * error indicator could not be obtained. The error indication value only tells, that
 * the remeshing is necessary, but no information about the actual error are provided.
 * The task of this class is to simply transform the indicator value to corresponding mesh size.
 */
class DirectErrorIndicatorRC : public RemeshingCriteria
{
protected:
    double minIndicatorLimit, maxIndicatorLimit;
    double minIndicatorDensity, maxIndicatorDensity;
    /// Default mesh density for Indicator value < minIndicatorLimit
    double zeroIndicatorDensity;
    /**
     * Ratio between proposedDensity and currDensity.
     * The remeshing is forced, whenever the actual ratio is smaller than
     * this value. Default value is 0.80
     */
    double remeshingDensityRatioToggle;
    FloatArray nodalDensities;
    /// Actual values (densities) state counter.
    StateCounterType stateCounter;
    RemeshingStrategy currStrategy;
#ifdef __PARALLEL_MODE
    std :: map< int, double >sharedDofManDensities;
    std :: map< int, double >sharedDofManIndicatorVals;
    bool dofManDensityExchangeFlag;
#endif

public:
    /// Constructor
    DirectErrorIndicatorRC(int n, ErrorEstimator *e);
    virtual ~DirectErrorIndicatorRC();

    virtual double giveRequiredDofManDensity(int num, TimeStep *tStep, int relative = 0);
    virtual double giveDofManDensity(int num);

    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual int estimateMeshDensities(TimeStep *tStep);
    virtual RemeshingStrategy giveRemeshingStrategy(TimeStep *tStep);
    /// returns the minIndicatorDensity
    double giveMinIndicatorLimit() { return minIndicatorLimit; }
    double giveMinIndicatorDensity() { return minIndicatorDensity; }

    void giveNodeChar(int inode, TimeStep *tStep, double &indicatorVal, double &currDensity);
    double giveZeroIndicatorDensity() { return zeroIndicatorDensity; }
    void reinitialize();

    virtual void setDomain(Domain *d);

protected:
    double giveLocalDofManDensity(int num);
    /**
     * Returns dof man indicator values.
     * @param num DofMan number.
     * @param tStep Time step.
     */
    double giveDofManIndicator(int num, TimeStep *tStep);
    double giveLocalDofManIndicator(int num, TimeStep *tStep);
#ifdef __PARALLEL_MODE
    void exchangeDofManDensities();
    int packSharedDofManLocalDensities(ProcessCommunicator &processComm);
    int unpackSharedDofManLocalDensities(ProcessCommunicator &processComm);

    void exchangeDofManIndicatorVals(TimeStep *tStep);
    int packSharedDofManLocalIndicatorVals(ProcessCommunicator &processComm);
    int unpackSharedDofManLocalIndicatorVals(ProcessCommunicator &processComm);
#endif
};
} // end namespace oofem
#endif // directerrorindicatorrc_h
