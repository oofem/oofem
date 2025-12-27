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

#ifndef entityrenumberingscheme_h
#define entityrenumberingscheme_h

#include <map>
#include <unordered_map>

#include "oofemenv.h"
#include "error.h"

namespace oofem {
/**
 * Type allowing to specify the required renumbering scheme;
 * One can have a renumbering scheme for dof managers
 * and another one for elements;
 */
enum EntityRenumberingScheme {
    ERS_DofManager,
    ERS_Element
};


class OOFEM_EXPORT EntityRenumberingFunctor
{
public:
    virtual ~EntityRenumberingFunctor() { }
    // two possible functions to call member function. virtual cause derived
    // classes will use a pointer to an object and a pointer to a member function
    // to make the function call
    /** 
     * Call using operator.
     * @return New number of component. Throws std::out_of_range if not found.
     */
    virtual int operator() (int, EntityRenumberingScheme) = 0;
    /// Call using function.
    virtual int call(int, EntityRenumberingScheme) = 0;
};

/// Derived template class
template< class TClass > class SpecificEntityRenumberingFunctor : public EntityRenumberingFunctor
{
private:
    /// Pointer to member function
    int ( TClass :: *fpt )(int, EntityRenumberingScheme);
    /// Pointer to object
    TClass *pt2Object;

public:

    // constructor - takes pointer to an object and pointer to a member and stores
    // them in two private variables
    SpecificEntityRenumberingFunctor( TClass * _pt2Object, int ( TClass :: *_fpt )(int, EntityRenumberingScheme) )
    {
        pt2Object = _pt2Object;
        fpt = _fpt;
    };

    virtual int operator() (int n, EntityRenumberingScheme ers)
    {
        return ( * pt2Object.*fpt )(n, ers);
    };

    virtual int call(int n, EntityRenumberingScheme ers)
    { return ( * pt2Object.*fpt )(n, ers); };
};


/// Renumbering functor based on provided maps.
class OOFEM_EXPORT MapBasedEntityRenumberingFunctor : public EntityRenumberingFunctor
{
private:
    std :: unordered_map< int, int > &dofmanMap, &elemMap;

public:
    MapBasedEntityRenumberingFunctor(std :: unordered_map< int, int > & _dofmanMap, std :: unordered_map< int, int > & _elemMap) :
        dofmanMap(_dofmanMap), elemMap(_elemMap)
    { }

    int operator() (int n, EntityRenumberingScheme ers) override
    {
        if ( ers == ERS_DofManager ) {
            auto pos = dofmanMap.find(n);
            if (pos != dofmanMap.end()) {
                return pos->second;
            }
        } else if ( ers == ERS_Element ) {
            auto pos = elemMap.find(n);
            if ( pos != elemMap.end() ) {
                return pos->second;
            }
        } else {
            OOFEM_ERROR("unsupported EntityRenumberingScheme");
        }
        throw std::out_of_range("entry not found");
        //OOFEM_ERROR("component label %d not found", n);
    }

    int call(int n, EntityRenumberingScheme ers) override
    { return this->operator() (n, ers); };
};
} // end namespace oofem
#endif // entityrenumberingscheme_h
