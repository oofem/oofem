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
 *               Copyright (C) 1993 - 2025   Borek Patzak
 *
 *
 *
 *       Czech Technical University, Faculty of Civil Engineering,
 *   Department of Structural Mechanics, 166 29 Prague, Czech Republic
 *
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with this library; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */

#ifndef contactsearch_h
#define contactsearch_h

#include "aabb.h"
#include "intarray.h"
#include <memory>
#include <vector>

namespace oofem {
class FloatArray;
class ContactPair;
class FEContactSurface;
class Domain;
class TimeStep;
class ContactElement;
class DataStream;

/**
 * @brief Abstract base class for contact search algorithms.
 *
 * The ContactSearch class defines a common interface for algorithms responsible
 * for detecting potential contact interactions between contact surfaces.
 * It provides the functionality required to identify candidate contact pairs
 * based on geometric proximity, bounding-box intersection, or other search
 * strategies.
 *
 * Concrete implementations may employ different search techniques, such as
 * brute-force search, spatial partitioning, hierarchical bounding volumes,
 * or grid-based methods. The output of a contact search is typically a set of
 * potential contact pairs that can be further processed by contact enforcement
 * algorithms.
 *
 * This class focuses solely on contact detection and does not prescribe any
 * particular contact formulation or constraint enforcement method.
 */
  
class ContactSearchAlgorithm
{
protected:
  Domain *domain;
  std::vector<std::unique_ptr<ContactPair>> contactPairs;

public:
  ContactSearchAlgorithm(Domain *d) {domain = d;}
  virtual ~ContactSearchAlgorithm() = default;
  /**
   * @brief Creates initial contact pairs based on the current configuration.
   *
   * Performs broad-phase (and possibly narrow-phase) detection to produce a set
   * of candidate contact pairs used by the contact boundary condition.
   */
  virtual void createContactPairs() = 0;
  /**
   * @brief Updates previously created contact pairs for the current time step.
   *
   * Typical responsibilities include re-projection, pair filtering, and updating
   * kinematic measures (gap, normals, tangents) as the configuration changes.
   *
   * @param tStep Current time step.
   */
  virtual void updateContactPairs(TimeStep *tStep) = 0;
  /** Commits broad-phase history after an accepted solution step. */
  virtual void commitSearchState(TimeStep *tStep);
  /** Stores committed search and pair state. */
  virtual void saveContext(DataStream &stream) const;
  /** Restores committed search and pair state. */
  virtual void restoreContext(DataStream &stream);
  /**
   * @brief Returns the internally stored list of contact pairs.
   *
   */
  std::vector<std::unique_ptr<ContactPair>>& getContactPairs() { return contactPairs; }
  /**
   * @brief Returns the internally stored list of contact pairs (const overload).
   */
  const std::vector<std::unique_ptr<ContactPair>>& getContactPairs() const { return contactPairs; }
};


  class ContactSearchAlgorithm_Surface2FESurface : public ContactSearchAlgorithm
{
protected: 
  FEContactSurface *slaveContactSurface;
  FEContactSurface *masterContactSurface;
  int surface_dimension;
  double searchPadding = 1.e-12;
  /// Negative selects the geometry-based default; nonnegative is an absolute distance.
  double configuredSearchPadding = -1.0;
  bool generalizedFeatures = false;
  bool directionalProjection = false;
  /**
   * Relative dead-band a rival master facet must beat the pair's currently
   * owned facet by (scaled by the candidates' own distance magnitude) before
   * it is allowed to take over ownership. Zero reproduces the previous
   * machine-epsilon-only tie-break, which lets ownership flip on an
   * infinitesimal, noise-level distance difference. This is needed at
   * near-degenerate shared edges/ridges between adjacent facets (e.g. the
   * saddle point of two crossing curved surfaces approximated by flat
   * quads), where the true closest facet is genuinely ambiguous and the tight
   * tie-break otherwise reassigns ownership on essentially every residual
   * evaluation, preventing Newton convergence. See
   * doc/contact-improvement-handoff.md, 2026-07-28 entries.
   */
  double facetOwnershipHysteresis = 0.0;
  std::vector<AABB> masterSearchAABBs;
  std::vector<AABB> committedMasterSearchAABBs;
  void updateMasterSearchAABBs(TimeStep *tStep);
  AABB computeSweptMasterSearchAABB(int masterIndex, TimeStep *tStep) const;
  std::vector<int> findCandidateMasterContactElements(ContactPair *cp, TimeStep *tStep);
  void appendCurrentActiveMasterCandidate(ContactPair *cp, std::vector<int> &candidates) const;
  void updatePairFrom3dCandidates(ContactPair *cp, const std::vector<int> &masterCandidateIndices, TimeStep *tStep);
  void updatePairFrom2dCandidates(ContactPair *cp, const std::vector<int> &masterCandidateIndices, TimeStep *tStep);
  double computeSearchPadding(TimeStep *tStep);
public:
  ContactSearchAlgorithm_Surface2FESurface(FEContactSurface *scs, FEContactSurface *mcs, Domain *d, int sd);
  ~ContactSearchAlgorithm_Surface2FESurface() override = default;
  void createContactPairs() override;
  void updateContactPairs(TimeStep *tStep) override{;}
  void setSearchPadding(double padding) { configuredSearchPadding = padding; }
  void setGeneralizedFeatures(bool enabled) { generalizedFeatures = enabled; }
  void setDirectionalProjection(bool enabled) { directionalProjection = enabled; }
  void setFacetOwnershipHysteresis(double relativeMargin) { facetOwnershipHysteresis = relativeMargin; }
  void commitSearchState(TimeStep *tStep) override;
  void saveContext(DataStream &stream) const override;
  void restoreContext(DataStream &stream) override;
};



  class ContactSearchAlgorithm_Surface2FESurface_3d : public ContactSearchAlgorithm_Surface2FESurface
{
protected:
public:
  ContactSearchAlgorithm_Surface2FESurface_3d(FEContactSurface *scs, FEContactSurface *mcs, Domain *d);
  ~ContactSearchAlgorithm_Surface2FESurface_3d() override = default;
  void updateContactPairs(TimeStep *tStep) override;
};
 
  class ContactSearchAlgorithm_Surface2FESurface_2d : public ContactSearchAlgorithm_Surface2FESurface
{
protected:
public:
  ContactSearchAlgorithm_Surface2FESurface_2d(FEContactSurface *scs, FEContactSurface *mcs, Domain *d);
  ~ContactSearchAlgorithm_Surface2FESurface_2d() override = default;
  void updateContactPairs(TimeStep *tStep) override;
};  



  


  
} // end namespace oofem
#endif //contactsearch_h
