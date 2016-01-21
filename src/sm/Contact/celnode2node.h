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

#ifndef celnode2node_h
#define celnode2node_h

#include "contact/contactelement.h"

///@name Input fields for _IFT_ContactElement
//@{
//#define _IFT_ContactManager_Name "contactmanager"

//@}

namespace oofem {
class Domain;
class ContactManager;
class ContactDefinition;
class SparseMtrx;
class TimeStep;
class DofManager;
class GaussPoint;
class UnknownNumberingScheme;
class FloatMatrix;
class IntegrationRule;
class ContactElement;

class OOFEM_EXPORT Node2NodeContact : public ContactElement
{
protected:
    ContactDefinition *cDef;

private:
    DofManager *masterNode;
    DofManager *slaveNode;
    
    // should be set by input:
    double area; // The area associated with the node (default = 1)- in order to represent some physical dimension.  
    double epsN;  // penalty stiffness 
    
    FloatArray normal;
public:

    /// Constructor.
    Node2NodeContact(DofManager *master, DofManager *slave);
    /// Destructor.
    virtual ~Node2NodeContact(){};
    virtual int instanciateYourself(DataReader *dr);
    virtual void setupIntegrationPoints();
    
    virtual void computeGap(FloatArray &answer, TimeStep *tStep);
    virtual void computeContactTractionAt(GaussPoint *gp, FloatArray &t, FloatArray &gap, TimeStep *tStep);
    virtual void computeCmatrixAt(GaussPoint *gp, FloatArray &answer, TimeStep *TimeStep);
    FloatArray &giveNormal() { return this->normal; }
    
    
    // Necessary methods - pure virtual in base class
    virtual void computeContactForces(FloatArray &answer, TimeStep *tStep, CharType type, ValueModeType mode,
                                const UnknownNumberingScheme &s, Domain *domain, FloatArray *eNorms);    
    
    virtual void computeContactTangent(FloatMatrix &answer, CharType type, TimeStep *tStep);
    
    virtual void giveLocationArray(IntArray &answer, const UnknownNumberingScheme &s);
};




class OOFEM_EXPORT Node2NodeContactL : public Node2NodeContact
{
protected:
    ContactDefinition *cDef;

private:
    DofManager *masterNode;
    DofManager *slaveNode;
    int lagrangeId; // dof Id associated with the Lagrange multiplier
    
    // should be set by input:
    double area; // The area associated with the node (default = 1)- in order to represent some physical dimension.  
  
    
public:

    /// Constructor.
    Node2NodeContactL(DofManager *master, DofManager *slave);
    /// Destructor.
    virtual ~Node2NodeContactL(){};
    virtual void giveDofManagersToAppendTo(IntArray &answer); 
    virtual void computeContactTractionAt(GaussPoint *gp, FloatArray &t, FloatArray &gap, TimeStep *tStep);
    
    // Necessary methods - pure virtual in base class
    virtual void computeContactForces(FloatArray &answer, TimeStep *tStep, CharType type, ValueModeType mode,
                                const UnknownNumberingScheme &s, Domain *domain, FloatArray *eNorms);    
    
    virtual void computeContactTangent(FloatMatrix &answer, CharType type, TimeStep *tStep);
    
    virtual void giveLocationArray(IntArray &answer, const UnknownNumberingScheme &s);
};



} // end namespace oofem
#endif // celnode2node_h
