/*
 * File:   enrichmentfunction.h
 * Author: chamrova
 *
 * Created on October 26, 2008, 3:15 PM
 */

#ifndef _ENRICHMENTFUNCTION_H
#define _ENRICHMENTFUNCTION_H

#include "femcmpnn.h"

class EnrichmentItem;
class BasicGeometry;

/** Abstract class representing EnrichmentItem */
class EnrichmentFunction : public FEMComponent
{
public:

    /**
     * Constructor.
     * @param n number associated with receiver
     * @param aDomain reference to domain.
     */
    EnrichmentFunction(int n, Domain *aDomain) : FEMComponent(n, aDomain) {}
    /// Destructor
    ~EnrichmentFunction() {};
    /// Evaluates a function at a particular point
    virtual void evaluateFunctionAt(FloatArray &answer, FloatArray *point) = 0;
    /// Evaluates a function derivative at a particular point
    virtual void evaluateDerivativeAt(FloatMatrix &answer, FloatArray *point) = 0;
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
    virtual void evaluateFunctionAt(FloatArray &answer, FloatArray *point);
    virtual void evaluateDerivativeAt(FloatMatrix &answer, FloatArray *point);
};

/** Class representing Branch EnrichmentFunction */
class BranchFunction : public EnrichmentFunction
{
public:

    BranchFunction(int n, Domain *aDomain) : EnrichmentFunction(n, aDomain) {
        this->numberOfDofs = 2;
    }
    virtual void evaluateFunctionAt(FloatArray &answer, FloatArray *point);
    virtual void evaluateDerivativeAt(FloatMatrix &answer, FloatArray *point);
};

#endif  /* _ENRICHMENTFUNCTION_H */


