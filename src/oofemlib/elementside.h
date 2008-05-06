/* $Header: /home/cvs/bp/oofem/oofemlib/src/elementside.h,v 1.8 2003/04/06 14:08:24 bp Exp $ */
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
 *               Copyright (C) 1993 - 2008   Borek Patzak
 *
 *
 *
 *       Czech Technical University, Faculty of Civil Engineering,
 *   Department of Structural Mechanics, 166 29 Prague, Czech Republic
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */


//   **************************
//   *** CLASS Element Side ***
//   **************************


#ifndef elementside_h
#define elementside_h

#include "dofmanager.h"
#include "domain.h"

#ifndef __MAKEDEPEND
#include <stdio.h>
#endif

class Dof;
class NodalLoad;
class TimeStep;
class FloatArray;
class IntArray;

/**
 * Class implementing element side having some DOFs in finite element mesh.
 * ElementSide posses degrees of freedom (see base class DofManager).
 * ElementSide id usually attribute of few elements and is managed by domain.
 */
class ElementSide : public DofManager
{
    /*
     * This class implements an element side in a finite element mesh. A side is an attri-
     * bute of a domain. It is usually also attribute of a few elements.
     * DESCRIPTION
     * The side possesses 'numberOfDofs' degrees of freedom, stored in 'dofArray'.
     * In 'loadArray' it stores the number of every nodal load it is subjected to
     * (typically, concentrated forces and moments).
     * In 'locationArray' the node stores the equation number of each of its dofs.
     * This location array is used by the side for assembling its load vector to
     * the right-hand side of the linear system ; it is also used by elements for
     * calculating their own location arrays.
     * TASKS
     * - managing its degrees of freedom (method 'giveDof') ;
     * - calculating its  load vector;
     * - printing and updating at end of step ;
     * - managing its swapping to and from disk.
     */

private:
public:
    /**
     * Constructor. Creates a element side belonging to domain.
     * @param n side number in domain aDomain
     * @param aDomain domain to which side belongs
     */
    ElementSide(int n, Domain *aDomain);            // constructor
    /// Destructor.
    ~ElementSide();                                 // destructor

    // miscellaneous
    /// Returns class name of the receiver.
    const char *giveClassName() const { return "Side"; }
    /** Returns classType id of receiver.
     * @see FEMComponent::giveClassID
     */
    classType    giveClassID() const { return ElementSideClass; }
    ///Initializes receiver acording to object description stored in input record.
    IRResultType initializeFrom(InputRecord *ir);
    //virtual IntArray* ResolveDofIDArray (char* initString);
    /// prints receiver state on stdout. Usefull for debuging.
    void         printYourself();

    /** Computes receiver transformation matrix from global cs. to dofManager specific
     * coordinate system (in which governing equations are assembled, for example the
     * local coordinate system in node).
     * @param answer computed transformation matrix. It has generally dofIDArry.size rows and
     * if loc is obtained using giveLocationArray(dofIDArry, loc) call, loc.giveSize() columns.
     * This is because this transformation should generally include not only transformation to
     * dof manager local coordinate system, but receiver dofs can be expressed using
     * dofs of another dofManager (In this case, squre answer is produced anly if all
     * dof transformation is required).
     * @param dofIDArry array containing DofIDItem-type values (this is enumeration
     * identifying physical meaning of particular DOF, see cltypes.h) for which transfromation mtrx is
     * assembled. if dofIDArry is NULL, then all receiver dofs are assumed.
     */
    virtual void computeTransformation(FloatMatrix &answer, const IntArray *dofIDArry);
    /**
     * Indicates, whether dofManager requires the transformation from global c.s. to
     * dof manager specific coordinate system.
     * @return nonzero if transformation is necessary, even for single dof.
     */
    virtual int requiresTransformation() { return 0; }
    /// Returns true if dof of given type is allowed to be associated to receiver
    virtual bool isDofTypeCompatible(dofType type) const { return ( type == DT_master || type == DT_simpleSlave ); }
};


#endif // elementside_h






