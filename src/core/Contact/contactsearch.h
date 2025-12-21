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

#include "intarray.h"
#include <memory>

namespace oofem {
class FloatArray;
class ContactPair;
class FEContactSurface;
class Domain;
class TimeStep;

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
  ~ContactSearchAlgorithm(){;}
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
public:
  ContactSearchAlgorithm_Surface2FESurface(FEContactSurface *scs, FEContactSurface *mcs, Domain *d, int sd);
  ~ContactSearchAlgorithm_Surface2FESurface(){;}
  void createContactPairs() override;
  void updateContactPairs(TimeStep *tStep) override{;}
};



  class ContactSearchAlgorithm_Surface2FESurface_3d : public ContactSearchAlgorithm_Surface2FESurface
{
protected:
public:
  ContactSearchAlgorithm_Surface2FESurface_3d(FEContactSurface *scs, FEContactSurface *mcs, Domain *d);
  ~ContactSearchAlgorithm_Surface2FESurface_3d(){;}
  void updateContactPairs(TimeStep *tStep) override;
};
 
  class ContactSearchAlgorithm_Surface2FESurface_2d : public ContactSearchAlgorithm_Surface2FESurface
{
protected:
public:
  ContactSearchAlgorithm_Surface2FESurface_2d(FEContactSurface *scs, FEContactSurface *mcs, Domain *d);
  ~ContactSearchAlgorithm_Surface2FESurface_2d(){;}
  void updateContactPairs(TimeStep *tStep) override;
};  



  


  
} // end namespace oofem
#endif //contactsearch_h
