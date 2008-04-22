/* $Header: /home/cvs/bp/oofem/sm/src/scalarerrorindicator.h,v 1.5 2003/04/06 14:08:31 bp Exp $ */
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

//   ************************************
//   *** CLASS SCALAR ERROR INDICATOR ***
//   ************************************

#ifndef scalarerrorindicator_h
#define scalarerrorindicator_h

#include "compiler.h"
#include "cltypes.h"
#include "errorestimator.h"

class Domain;
class Element;
class TimeStep;
class RemeshingCriteria;


/**
 * The class representing scalar error indicator.
 * It indicates element error based on the value of some suitable scalar value obtained from the
 * element integration points and corresponding material model.
 */
class ScalarErrorIndicator : public ErrorEstimator
{
protected:
    /// type of internal variable to be indicator (type for temp and nontemp varsion)
    int indicatorType;
    /// corresponding internal state type
    InternalStateType varType;
public:
    /// Constructor
    ScalarErrorIndicator(int n, Domain *d) : ErrorEstimator(n, d) { eeType = EET_SEI; }
    /// Destructor
    virtual ~ScalarErrorIndicator() { }
    /** Returns the element error of requested type.
     * @param type error type
     * @param elem element for which error requested
     * @param tStep time step
     */
    virtual double giveElementError(EE_ErrorType type, Element *elem, TimeStep *tStep);
    /** Returns the global domain value of given type.
     * @param type error type
     * @param tStep time step
     */
    virtual double giveValue(EE_ValueType type, TimeStep *tStep) { return 0.0; }
    /**
     * Estimates the error on associated domain at given timeSte.
     * Empty implementation, the indicator value is assumed to be directly the scalar internal variable
     * obtained from elements andcorresponding IP.
     * @param tStep time step
     */
    virtual int estimateError(EE_ErrorMode mode, TimeStep *tStep);
    /** Returns reference to associated remeshing criteria.
     */
    virtual RemeshingCriteria *giveRemeshingCrit();
    /* Initalizes the receiver from input record */
    virtual IRResultType initializeFrom(InputRecord *ir);
    /// Returns class name of the receiver.
    const char *giveClassName() const { return "ScalarErrorIndicator"; }
    /** Returns classType id of receiver.
     * @see FEMComponent::giveClassID
     */
    classType                giveClassID() const { return ScalarErrorIndicatorClass; }

protected:
};

#endif // scalarerrorindicator_h
