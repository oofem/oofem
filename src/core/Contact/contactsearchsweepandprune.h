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
 *               Copyright (C) 1993 - 2013   Borek Patzak
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
  ContactSearchAlgorithm_Surface2FESurface_3d_SweepAndPrune(FEContactSurface *scs, FEContactSurface *mcs, Domain *d);
  ~ContactSearchAlgorithm_Surface2FESurface_3d_SweepAndPrune(){;}
  void updateContactPairs(TimeStep *tStep) override;
protected:
    void initialize();
    void updateAABBs();
    void updateBoundss();
    void insertionSort();
    void addPotentialPair(IJ);
    void deletePotentialPair(IJ);
    void solveInversion(int axis, const Bound& b1, const Bound& b2);
    static bool contains(const SetIJ& set, IJ ij) {
      // return set.contains(ij); // C++20
      return set.find(ij) != set.end();
    }
};
  
} // end namespace oofem
#endif //contactsearchsweepandprune_h
