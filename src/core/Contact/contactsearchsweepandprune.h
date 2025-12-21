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

#ifndef contactsearchsweepandprune_h
#define contactsearchsweepandprune_h

#include "contactsearch.h"
#include "contactelement.h"
#include "contactpair.h"
#include "aabb.h"
#include <memory>
#include <map>
#include <set>
#include <tuple>
#include <array>
#include <utility>
namespace oofem {

struct Bound {
  public:
    double value;
    unsigned int id;
    bool isSlave;
    bool isMin;
};


 /**
 * @brief Contact search algorithm based on the sweep-and-prune strategy.
 *
 * The ContactSearchSweepAndPrune class implements a contact detection algorithm
 * using the sweep-and-prune method. The algorithm projects bounding volumes of
 * contact surfaces or contact entities onto one or more coordinate axes and
 * efficiently identifies potentially intersecting pairs by sorting and
 * sweeping these projections.
 *
 * This approach is particularly efficient for problems with a large number of
 * contact entities and relatively small changes in their spatial configuration
 * between consecutive time steps, as it exploits temporal coherence of the
 * motion.
 *
 * The class is intended to be used as a concrete implementation of the
 * ContactSearch interface and provides candidate contact pairs that can be
 * further processed by contact evaluation and enforcement algorithms.
 */ 

class ContactSearchAlgorithm_Surface2FESurface_3d_SweepAndPrune : public ContactSearchAlgorithm_Surface2FESurface_3d
{
public:
  using IJ = std::pair<unsigned int, unsigned int>;
  using SetIJ = std::set<IJ>;
protected:
  SetIJ potentialPairs;
  std::array<SetIJ,3> potentialPairsPerAxes;
  std::array<std::vector<Bound>,3> boundss;
  std::vector<AABB> aabbs;
  bool isInitialized;
public:
  /**
   * @brief Constructs a 3D surface-to-surface contact search algorithm using sweep-and-prune.
   *
   * Initializes the search algorithm for a given slave and master FE contact surface
   * within the provided computational domain. The algorithm will use axis-aligned
   * bounding boxes (AABBs) and per-axis endpoint sorting to efficiently identify
   * potential contact candidate pairs.
   *
   * @param scs Pointer to the slave FE contact surface.
   * @param mcs Pointer to the master FE contact surface.
   * @param d   Pointer to the associated domain.
   */
  ContactSearchAlgorithm_Surface2FESurface_3d_SweepAndPrune(FEContactSurface *scs, FEContactSurface *mcs, Domain *d);
  ~ContactSearchAlgorithm_Surface2FESurface_3d_SweepAndPrune(){;}
  /**
   * @brief Updates the list of contact pairs using sweep-and-prune broad-phase detection.
   *
   * Updates internal geometric data (AABBs and axis bounds) for the current time step,
   * performs the sweep-and-prune procedure to obtain candidate master–slave pairs,
   * and updates the internally stored set of potential pairs used for contact evaluation.
   *
   * @param tStep Current time step.
   */
  void updateContactPairs(TimeStep *tStep) override;
protected:
  /**
   * @brief Initializes internal data structures for sweep-and-prune.
   *
   * Builds initial AABBs and per-axis bound lists, and initializes the sets
   * tracking potential pairs. Typically called on the first update before
   * incremental sorting is applied.
   */
    void initialize();
  /**
   * @brief Recomputes AABBs used by the search algorithm.
   *
   * Updates the axis-aligned bounding boxes of the contact entities participating
   * in the search (typically the slave entities, and/or both sides depending on
   * implementation details). These AABBs are the input to the sweep-and-prune
   * endpoint generation.
   */
    void updateAABBs();
  /**
   * @brief Updates per-axis bound (endpoint) arrays from current AABBs.
   *
   * Rebuilds the arrays of min/max endpoints (bounds) for each coordinate axis
   * from the current AABBs. These bound arrays are then processed by sorting
   * (initially or incrementally) to detect overlaps efficiently.
   */
    void updateBoundss();
  /**
   * @brief Performs incremental sorting of bounds using insertion sort.
   *
   * Applies insertion sort to the per-axis bound arrays to exploit temporal coherence
   * (small changes between time steps). During sorting, inversions (swapped endpoints)
   * are detected and handled to update the set of potential overlapping pairs.
   */
    void insertionSort();
  /**
   * @brief Adds a candidate pair to the set of potential contact pairs.
   *
   * Inserts the pair identifier into the global potential-pair set and/or
   * updates auxiliary per-axis sets. A pair is typically promoted to a global
   * potential pair once overlap has been detected consistently across axes.
   *
   * @param ij Pair identifier (e.g. indices of slave/master entities).
   */
    void addPotentialPair(IJ);
  /**
   * @brief Removes a candidate pair from the set of potential contact pairs.
   *
   * Deletes the pair identifier from the global potential-pair set and/or
   * updates auxiliary per-axis sets when overlap is no longer present.
   *
   * @param ij Pair identifier (e.g. indices of slave/master entities).
   */
    void deletePotentialPair(IJ);
  /**
   * @brief Handles a bound inversion event detected during sorting on a given axis.
   *
   * When two neighboring bounds swap order during incremental sorting, the overlap
   * status of the corresponding entities along the given axis may change. This method
   * updates the per-axis potential-pair bookkeeping (and possibly the global set)
   * based on whether the inversion creates or destroys an overlap interval.
   *
   * @param axis Axis index on which the inversion occurred (0=x, 1=y, 2=z).
   * @param b1   First bound involved in the inversion.
   * @param b2   Second bound involved in the inversion.
   */
    void solveInversion(int axis, const Bound& b1, const Bound& b2);
  /**
   * @brief Utility: checks whether a pair identifier is contained in a set.
   * @param set Set of pair identifiers.
   * @param ij  Pair identifier to check.
   * @return True if @p ij is present in @p set, false otherwise.
   */
    static bool contains(const SetIJ& set, IJ ij) {
      // return set.contains(ij); // C++20
      return set.find(ij) != set.end();
    }
};
  
} // end namespace oofem
#endif //contactsearchsweepandprune_h
