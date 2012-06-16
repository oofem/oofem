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

#ifndef enrichmentfunction_h
#define enrichmentfunction_h

#include "femcmpnn.h"
#include "intarray.h"

namespace oofem {
class EnrichmentItem;
class BasicGeometry;
class GaussPoint;

/**
 * Abstract class representing global shape function
 * Base class declares abstract interface common to all enrichment functions.
 * Particularly, evaluateAt() and evaluateDerivativeAt() services are declared
 * to evaluate the value and spatial derivatives at a given point of the receiver.
 * @author chamrova
 */
class EnrichmentFunction : public FEMComponent
{
public:
    /**
     * Constructor.
     * @param n Number associated with receiver.
     * @param aDomain Reference to domain.
     */
    EnrichmentFunction(int n, Domain *aDomain) : FEMComponent(n, aDomain) { }
    /// Destructor
    virtual ~EnrichmentFunction() { };
    /// Evaluates a function at a particular point
    virtual double evaluateFunctionAt(FloatArray *point, EnrichmentItem *ei) = 0;
    /// Evaluates a function derivative at a particular point
    virtual void evaluateDerivativeAt(FloatArray &answer, FloatArray *point, EnrichmentItem *ei) = 0;
    /// Evaluates a function at a particular point
    virtual double evaluateFunctionAt(GaussPoint *gp, EnrichmentItem *ei);
    /// Evaluates a function derivative at a particular point
    virtual void evaluateDerivativeAt(FloatArray &answer, GaussPoint *gp, EnrichmentItem *ei);
    // Inserts EnrichmentItem into associatedEnrItem array
    // void insertEnrichmentItem(EnrichmentItem *er);
    // Sets a particular EnrichmentItem active
    // void setActive(EnrichmentItem *er);
    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual const char *giveClassName() const { return "EnrichmentFunction"; }
    /// Accessor.
    int giveNumberOfDofs() { return numberOfDofs; }

protected:
    /// EnrichmentItems associated with this EnrichmentFunction.
    IntArray assocEnrItemArray;
    /// Active EnrichmentItem.
    EnrichmentItem *activeEnrItem;
    /// Number of dofs to enrich.
    int numberOfDofs;
};

/** Class representing Heaviside EnrichmentFunction. */
class DiscontinuousFunction : public EnrichmentFunction
{
public:

    DiscontinuousFunction(int n, Domain *aDomain) : EnrichmentFunction(n, aDomain) {
        this->numberOfDofs = 2;
    }
    double evaluateFunctionAt(FloatArray *point, EnrichmentItem *ei);
    void evaluateDerivativeAt(FloatArray &answer, FloatArray *point, EnrichmentItem *ei);
};

/** Class representing Branch EnrichmentFunction. */
class BranchFunction : public EnrichmentFunction
{
public:

    BranchFunction(int n, Domain *aDomain) : EnrichmentFunction(n, aDomain) {
        this->numberOfDofs = 2;
    }
    double evaluateFunctionAt(FloatArray *point, EnrichmentItem *ei);
    void evaluateDerivativeAt(FloatArray &answer, FloatArray *point, EnrichmentItem *ei);
};

/** Class representing bimaterial interface. */
class RampFunction : public EnrichmentFunction
{
public:

    RampFunction(int n, Domain *aDomain) : EnrichmentFunction(n, aDomain) {
        this->numberOfDofs = 2;
    }
    double evaluateFunctionAt(FloatArray *point, EnrichmentItem *ei);
    void evaluateDerivativeAt(FloatArray &answer, FloatArray *point, EnrichmentItem *ei);
    double evaluateFunctionAt(GaussPoint *gp, EnrichmentItem *ei);
    void evaluateDerivativeAt(FloatArray &answer, GaussPoint *gp, EnrichmentItem *ei);
};
} // end namespace oofem
#endif  // enrichmentfunction_h
