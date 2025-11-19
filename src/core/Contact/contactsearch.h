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
 * Abstract base class for all contact finite elements. Derived classes should be  base
 * classes for specific analysis type (for example base class for structural analysis,
 * thermal analysis or magnetostatics one). These derived classes then declare
 * analysis-specific part of interface and they provide default implementation
 * for these methods.
 * This abstract class declares (and possibly implements) general data and methods
 * common to all element types. General methods for obtaining characteristic vectors,
 * matrices and values are introduced and should be used instead of calling directly
 * specific member functions (these must be overloaded by derived analysis-specific
 * classes in order to invoke proper method according to type of component requested).
 */

class ContactSearchAlgorithm
{
protected:
  Domain *domain;
  std::vector<std::unique_ptr<ContactPair>> contactPairs;

public:
  ContactSearchAlgorithm(Domain *d) {domain = d;}
  ~ContactSearchAlgorithm(){;}
  virtual void createContactPairs() = 0;
  virtual void updateContactPairs(TimeStep *tStep) = 0;
  std::vector<std::unique_ptr<ContactPair>>& getContactPairs() { return contactPairs; }
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
