/* $Header: /home/cvs/bp/oofem/oofemlib/src/load.h,v 1.10 2003/04/18 10:37:07 bp Exp $ */
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


//   ******************
//   *** CLASS LOAD ***
//   ******************


#ifndef load_h
#define load_h

#include "generalbc.h"
#include "domain.h"
#include "flotarry.h"
#include "dictionr.h"

/**
 * Load is base abstract class for all loads.
 * Load is an aribute of the domain (it belongs to).
 * Load is also atrribute of several elements, nodes,
 * which are subjected to loading type boundary condition.
 *
 * The value (or the components) of a load will be
 * the product of its value (stored in componentArray) by the value of
 * the associated load time function at given time step.
 */
class Load : public GeneralBoundaryCondition
{
    /*
     * This abstract class is the superclass of the classes that implement loads
     * (body load, nodal load, boundary conditions, etc). A load is an attribute
     * of the domain. It is usually also attribute of several elements, nodes or
     * dofs.
     * DESCRIPTION
     * The load stores its values in 'componentArray'. The components of a load
     * at a given time step is the product of 'componentArray' by the value of
     * the function 'loadTimeFunction' at that time step.
     * TASK
     * Returning its components and its load-time function ;
     */
protected:
    /// Components of boundary condition
    FloatArray componentArray;
    /** The load is pecified for all dofs of object to which is associated.
     * For some types of boundary conditions the zero value of load does not mean
     * that the load is not applied (newton's type of bc, for example). Then
     * some mask, which allows to exclude specific dofs is necessary.
     * The dofMask attribute is introduced to alow this.
     * By default it is of the same size as componentArray, filled with zeroes.
     * If some value of dofExcludeMask is set to nonzero, then the corresponding componentArray
     * is set to zero. */
    IntArray dofExcludeMask;
public:

    /**
     * Constructor. Creates boundary condition with given number, belonging to given domain.
     * @param n boundary condition number
     * @param d domain to which new object will belongs.
     */
    Load(int, Domain *);                           // constructor
    /// Destructor.
    virtual ~Load()  { }   // destructor

    /**
     * Computes boundary condition value - its components values at given time.
     * Default implementation returns as the answer its component array multiplied
     * with load time function value (load response mode is taken in to account)
     * @param answer computed boundary conditions components
     * @param stepN time step, for which components are computed.
     * @param mode determines response mode.
     */
    virtual void  computeComponentArrayAt(FloatArray &answer, TimeStep *, ValueModeType mode);
    /**
     * Computes components values of load at given point - global coordinates (coordinates given).
     * Default implementation computes product of aproximation matrix (computeNArray service) and
     * with "vertex" value array attribute and the result is then multiplied by
     * corresponding load time function value respecting load response mode.
     * @param answer component values at given point and time
     * @param stepN time step representing time
     * @param coords global (or local) problem coordinates, which are used to
     * evaluate components values.
     * @param mode determines response mode.
     */
    virtual void         computeValueAt(FloatArray &answer, TimeStep *atTime, FloatArray &coords, ValueModeType mode) = 0;
    /**
     * Returns the value of dofExcludeMask corresponding to given indx.
     * See the description of dofExcludeMask attribute for more details.
     */
    int isDofExcluded(int indx);
    // definition of a load
    /// Initializes receiver acording to object description stored in input record.
    IRResultType initializeFrom(InputRecord *ir);
    /** Setups the input record string of receiver
     * @param str string to be filled by input record
     * @param keyword print record keyword (default true)
     */
    virtual int giveInputRecordString(std :: string &str, bool keyword = true);

    /// Returns classType id of receiver.
    classType    giveClassID() const { return LoadClass; }
    /// Returns class name of the receiver.
    const char *giveClassName() const { return "Load"; }


protected:
    /**
     * Returns pointer to receiver component array, where component values of boundary condition are stored.
     */
    FloatArray &giveComponentArray();
};

#endif // load_h

