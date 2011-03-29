#ifndef _PATCHINTEGRATIONRULE_H
#define _PATCHINTEGRATIONRULE_H

#include "integrationrule.h"
#include "patch.h"
#include "alist.h"
#include "gaussintegrationrule.h"

namespace oofem {
/*
 * class PatchIntegrationRule : public IntegrationRule
 * {
 * protected:
 *  /// list of patches
 *  AList< Patch > *patchSet;
 * public:
 *  /// Constructor
 *  PatchIntegrationRule(int n, Element *e, AList< Triangle > *patchSet);
 *  /// Destructor
 *  ~PatchIntegrationRule();
 *  /// Performes integration over the whole patchset
 *  void computeGps(MaterialMode matMode);
 *  /// Accessor
 *  Patch *givePatch(int n) { return patchSet->at(n); }
 *  /// Computes gps on a single patch
 *  void computeGpsForPatch(Patch *patch, AList< GaussPoint > *answer, MaterialMode matMode);
 *  /// Wrap up function for computeGps
 *  int SetUpPointsOnSquare(int, MaterialMode matMode, GaussPoint ***);
 *  int giveNumberOfPatches() { return patchSet->giveSize(); }
 * };*/

/**
 * Represents an IntegrationRule for an interacted element
 * the standard integration is replaced by an integration over a
 * patchset
 */
class PatchIntegrationRule : public GaussIntegrationRule
{
protected:
    /// patch
    Patch *patch;

public:
    /// Constructor
    PatchIntegrationRule(int n, Element *e, Patch *p);
    /// Destructor
    ~PatchIntegrationRule();
    int SetUpPointsOnTriagle(int, MaterialMode, GaussPoint * * *);
    int giveMaterial() { return this->patch->giveMaterial(); }
    Patch *givePatch() { return this->patch; }

    virtual contextIOResultType saveContext(DataStream *stream, ContextMode mode, void *obj);
    virtual contextIOResultType restoreContext(DataStream *stream, ContextMode mode, void *obj);
    virtual classType giveClassID() const { return PatchIntegrationRuleClass; }
};
} // end namespace oofem
#endif
