/* 
 * File:   enrichmentfunction.h
 * Author: chamrova
 *
 * Created on October 26, 2008, 3:15 PM
 */

#ifndef _ENRICHMENTFUNCTION_H
#define	_ENRICHMENTFUNCTION_H

#include "femcmpnn.h"
class EnrichmentItem; class Geometry;

class EnrichmentFunction : public FEMComponent {
    public:
         EnrichmentFunction (int n,Domain* aDomain) : FEMComponent (n, aDomain){}   // constructor
         ~EnrichmentFunction () {};
         virtual double evaluateFunctionAt(FloatArray* point) = 0;
         virtual void evaluateDerivativeAt(FloatArray &answer, FloatArray *point) = 0;
         Geometry*  giveGeometry ();
         void insertEnrichmentItem(EnrichmentItem *er);
         void setActive(EnrichmentItem *er);
         IRResultType initializeFrom (InputRecord* ir);
         const char* giveClassName () const { return "EnrichmentFunction" ; }
         void setNumber (int number) { this->number = number;}
         
    protected:
         IntArray assocEnrItemArray;
         EnrichmentItem* activeEnrItem;
         // for branch functions
         int number;
    
};


class DiscontinuousFunction : public EnrichmentFunction {
    public:
        DiscontinuousFunction(int n,Domain* aDomain): EnrichmentFunction(n, aDomain) {}
        double evaluateFunctionAt(FloatArray* point);
        void evaluateDerivativeAt(FloatArray &answer, FloatArray *point);
};

class BranchFunction : public EnrichmentFunction {
    public:
        BranchFunction(int n,Domain* aDomain): EnrichmentFunction(n, aDomain) {}
        double evaluateFunctionAt(FloatArray* point);
        void evaluateDerivativeAt(FloatArray &answer, FloatArray *point);
};

#endif	/* _ENRICHMENTFUNCTION_H */


