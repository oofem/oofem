# include "enrichmentfunction.h"
# include "enrichmentitem.h"
# include "domain.h"
# include <math.h>

// temporary function
void EnrichmentFunction::insertEnrichmentItem(EnrichmentItem *er){
    int sz = assocEnrItemArray.giveSize();
    assocEnrItemArray.resize(sz+1);
    this->assocEnrItemArray.at(sz+1) = er->giveNumber();
}

void EnrichmentFunction::setActive(EnrichmentItem *er){
    this->activeEnrItem = er->giveNumber();
}

IRResultType EnrichmentFunction::initializeFrom (InputRecord* ir)
{
    return IRRT_OK;
} 

double DiscontinuousFunction::evaluateFunctionAt(FloatArray *point){
    double dist = domain->giveEnrichmentItem(activeEnrItem)->giveGeometry()->giveDistanceTo(point);
    //  here instead of zero some tolerance should be inserted
    if(dist > 0){
        return 1.0;
    }
    else return -1.0;      
}

