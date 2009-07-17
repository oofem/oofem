/*
 * File:   enrichmentfunction.h
 * Author: chamrova
 *
 * Created on October 26, 2008, 3:15 PM
 */

#ifndef _ENRICHMENTFUNCTION_H
#define _ENRICHMENTFUNCTION_H

#include "femcmpnn.h"
#include "gausspnt.h"

class EnrichmentItem;
class BasicGeometry;

/** Abstract class representing global shape function
 *  Base class declares abstract interface common to all enrichment functions.
 *  Particularly, evaluateAt() and evaluateDerivativeAt() services are declared
 *  to evaluate the value and spatial derivatives at a given point of the receiver.
 */
class EnrichmentFunction : public FEMComponent
{
public:

    /**
     * Constructor.
     * @param n number associated with receiver
     * @param aDomain reference to domain.
     */
    EnrichmentFunction(int n, Domain *aDomain) : FEMComponent(n, aDomain) { }
    /// Destructor
    ~EnrichmentFunction() { };
    /// Evaluates a function at a particular point
    virtual double evaluateFunctionAt(FloatArray *point) = 0;
    /// Evaluates a function derivative at a particular point
    virtual void evaluateDerivativeAt(FloatArray &answer, FloatArray *point) = 0;
    /// Evaluates a function at a particular point
    virtual double evaluateFunctionAt(GaussPoint *gp);
    /// Evaluates a function derivative at a particular point
    virtual void evaluateDerivativeAt(FloatArray &answer, GaussPoint *gp);
    /// Accessor
    BasicGeometry *giveGeometry();
    /// Inserts EnrichmentItem into associatedEnrItem array
    void insertEnrichmentItem(EnrichmentItem *er);
    /// Sets a particular EnrichmentItem active
    void setActive(EnrichmentItem *er);
    /// Initializes EnrichmentItem from InputRecord
    IRResultType initializeFrom(InputRecord *ir);

    const char *giveClassName() const {
        return "EnrichmentFunction";
    }
    /// Accessor
    int giveNumberOfDofs() { return numberOfDofs; }

protected:
    /// EnrichmentItems associated with this EnrichmentFunction
    IntArray assocEnrItemArray;
    /// active EnrichmentItem
    EnrichmentItem *activeEnrItem;
    // number of dofs to enrich
    int numberOfDofs;
};

/** Class representing Heaviside EnrichmentFunction */
class DiscontinuousFunction : public EnrichmentFunction
{
public:

    DiscontinuousFunction(int n, Domain *aDomain) : EnrichmentFunction(n, aDomain) {
        this->numberOfDofs = 2;
    }
    double evaluateFunctionAt(FloatArray *point);
    void evaluateDerivativeAt(FloatArray &answer, FloatArray *point);
};

/** Class representing Branch EnrichmentFunction */
class BranchFunction : public EnrichmentFunction
{
public:

    BranchFunction(int n, Domain *aDomain) : EnrichmentFunction(n, aDomain) {
        this->numberOfDofs = 2;
    }
    double evaluateFunctionAt(FloatArray *point);
    void evaluateDerivativeAt(FloatArray &answer, FloatArray *point);
};

/** Class representing bimaterial interface */
class RampFunction : public EnrichmentFunction
{
public:

    RampFunction(int n, Domain *aDomain) : EnrichmentFunction(n, aDomain) {
        this->numberOfDofs = 2;
    }
    double evaluateFunctionAt(FloatArray *point);
    void evaluateDerivativeAt(FloatArray &answer, FloatArray *point);
    double evaluateFunctionAt(GaussPoint *gp);
    void evaluateDerivativeAt(FloatArray &answer, GaussPoint *gp);
};

#endif  /* _ENRICHMENTFUNCTION_H */


