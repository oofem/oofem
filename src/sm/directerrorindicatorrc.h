/* $Header: /home/cvs/bp/oofem/sm/src/directerrorindicatorrc.h,v 1.4 2003/04/06 14:08:30 bp Exp $ */
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

//   *******************************************************
//   *** CLASS DIRECT ERROR INDICATOR REMESHING CRITERIA ***
//   *******************************************************

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
protected:
public:
    /// Constructor
    DirectErrorIndicatorRCInterface() : Interface() { }
    /*
     * Determines the characteristic size of element. This quantity is defined as follows:
     * For 1D it is the element length, for 2D it is the square root of element area.
     */
    virtual double DirectErrorIndicatorRCI_giveCharacteristicSize() = 0;
};




/**
 * The class is an implementation of "direct" remeshing criteria, which
 * maps the error indication, which is ususally the value of observed internal variable
 * to the corresponding required element size. The error estimate and global error for
 * error undicator could not be obtained. The error indication value only tells, that
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
    /// actual values (densities) state counter.
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
    /** Returns the required mesh size n given dof manager.
     * The mesh density is defined as a required element size
     * (in 1D the element length, in 2D the square from element area).
     * @param num dofman  number
     * @param tStep time step
     */
    virtual double giveRequiredDofManDensity(int num, TimeStep *tStep, int relative = 0);
    /**
     * Returns existing mesh size for given dof manager.
     * @param num dofMan number
     */
    virtual double giveDofManDensity(int num);


    /** Initializes receiver acording to object description stored in input record.
     * This function is called immediately after creating object using
     * constructor. InitString can be imagined as data record in component database
     * belonging to receiver. Receiver may use value-name extracting functions
     * to extract particular field from record.
     * @see readInteger, readDouble and similar functions */
    virtual IRResultType initializeFrom(InputRecord *ir);
    /**
     * Estimates the nodal densities.
     * @param tStep time step
     */
    virtual int estimateMeshDensities(TimeStep *tStep);
    /**
     * Determines, if the remeshing is needed, and if nedded, the type of strategy used
     */
    virtual RemeshingStrategy giveRemeshingStrategy(TimeStep *tStep);
    /// returns the minIndicatorDensity
    double giveMinIndicatorLimit() { return minIndicatorLimit; }
    double giveMinIndicatorDensity() { return minIndicatorDensity; }

    //protected:
    void giveNodeChar(int inode, TimeStep *tStep, double &indicatorVal, double &currDensity);
    double giveZeroIndicatorDensity() { return zeroIndicatorDensity; }
    void reinitialize();
    /// sets associated Domain
    virtual void         setDomain(Domain *d);

protected:
    double giveLocalDofManDensity(int num);
    /**
     * Returns dof man indicator values.
     * @param num dofMan number
     */

    double giveDofManIndicator(int num, TimeStep *);
    double giveLocalDofManIndicator(int num, TimeStep *);
#ifdef __PARALLEL_MODE
    void exchangeDofManDensities();
    int packSharedDofManLocalDensities(ProcessCommunicator &processComm);
    int unpackSharedDofManLocalDensities(ProcessCommunicator &processComm);

    void exchangeDofManIndicatorVals(TimeStep *);
    int packSharedDofManLocalIndicatorVals(ProcessCommunicator &processComm);
    int unpackSharedDofManLocalIndicatorVals(ProcessCommunicator &processComm);
#endif
};
} // end namespace oofem
#endif // directerrorindicatorrc_h
