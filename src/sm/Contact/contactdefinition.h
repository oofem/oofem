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

#ifndef contactdefinition_h
#define contactdefinition_h

#include "oofemcfg.h"
#include "datareader.h"
#include "inputrecord.h"


#include <unordered_map>
#include <list>
#include <vector>
#include <memory>

///@name Input fields for _IFT_ContactManager
//@{
//#define _IFT_ContactManager_Name "contactmanager"

//@}

namespace oofem {
class Domain;
class ContactManager;
class ContactObject;
class ContactElement;



/**
 * This class manages a particular contact definition. 
 * This keeps track of the discretization, how the contact constraints are enforced 
 *
 * @author Jim Brouzoulis
 */
class OOFEM_EXPORT ContactDefinition
{
protected:
    ContactManager *cMan;

    std :: vector< ContactElement *> masterElementList;
    
public:

    /// Constructor.
    ContactDefinition(ContactManager *cMan);
    /// Destructor.
    virtual ~ContactDefinition();

    virtual IRResultType initializeFrom(InputRecord *ir);;

    virtual int instanciateYourself(DataReader *dr);
    virtual const char *giveClassName() const { return "ContactDefinition"; }
    //virtual const char *giveInputRecordName() const { return _IFT_ContactDefinition_Name; }

    
    virtual void computeContactForces(FloatArray &answer, TimeStep *tStep, CharType type, ValueModeType mode,
                                const UnknownNumberingScheme &s, Domain *domain, FloatArray *eNorms);
    
    virtual void computeContactTangent(SparseMtrx *answer, TimeStep *tStep,
                      CharType type, const UnknownNumberingScheme &r_s, const UnknownNumberingScheme &c_s); 
    
    ContactElement *giveContactElemnt(const int num) { return this->masterElementList[num-1]; };
    int giveNumbertOfContactElements() { return this->masterElementList.size(); };
    
    // objects can be of different kinds
    // nodes, segments, surfaces, analytical functions
    
    /* Each contact definition has several master objects that each keep track of generally several slave objects
     * Ex:
     * -Master can be a node and only keep track of another node -> node2node
     * -Master can be a node and keep track of several other nodes -> for larger displacements
     * -Master can be a segment and only keep track of another segment -> segment2segment 
     * -etc.
     *
     * how should they be stored? 
     */
    //MasterObjects
    //SlaveObjects
    // assembleVectorOf...
    // assembleTangentOf...
    
    
    
    //ConstrainType - How should the contact constraints be fulfilled
    // Penalty, Lagrange multiplier, Augmented (mix), Mortar (weakly fulfilled)
    
    // InterfaceModel - a constitutive model for the contact - depends on the physics
    // linear/nonlinear for stresses, thermal conductance, etc.
    // Normal/tangential
    
    
    
};








} // end namespace oofem
#endif // contactdefinition_h
