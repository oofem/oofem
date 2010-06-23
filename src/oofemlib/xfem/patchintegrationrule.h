#ifndef _PATCHINTEGRATIONRULE_H
#define _PATCHINTEGRATIONRULE_H

#include "integrationrule.h"
#include "patch.h"
#include "alist.h"
#include "gaussintegrationrule.h"

namespace oofem {

/** represents an IntegrationRule for an interacted element
 * the standard integration is replaced by an integration over a
 * patchset
 */
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

class PatchIntegrationRule : public GaussIntegrationRule
{
protected:
    /// patch
    Patch *patch;

public:
    /// Constructor
    PatchIntegrationRule(int n, Element *e, Patch *p); // HUHU
    /// Destructor
    ~PatchIntegrationRule();
    int SetUpPointsOnTriagle(int, MaterialMode, GaussPoint * * *);
    int giveMaterial() { return this->patch->giveMaterial(); } // HUHU
    Patch *givePatch() { return this->patch; }
    /**
     * Saves receiver's context to stream.
     * Calls saveContext service for all receiver's integration points.
     * Note: does not call the FEMComponent::saveContext service, in order not
     * to write class id info for each integration rule.
     * @exception throws an ContextIOERR exception if error encountered.
     */
    virtual contextIOResultType saveContext(DataStream *stream, ContextMode mode, void *obj);
    /**
     * Restores receiver's context to stream.
     * Calls restoreContext service for all receiver's integration points.
     * Note: does not call the FEMComponent::restoreContext service, in order not
     * to write class id info for each integration rule.
     * @param obj should be a pointer to invoking element, ie., to which the receiver will belong to.
     * @exception throws an ContextIOERR exception if error encountered.
     */
    virtual contextIOResultType restoreContext(DataStream *stream, ContextMode mode, void *obj);
    ///Returns classType id of receiver.
    virtual classType giveClassID() const { return PatchIntegrationRuleClass; }

};

} // end namespace oofem
#endif
