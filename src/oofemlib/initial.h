/* $Header: /home/cvs/bp/oofem/oofemlib/src/initial.h,v 1.8 2003/04/06 14:08:24 bp Exp $ */
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


//   *******************************
//   *** CLASS INITIAL CONDITION ***
//   *******************************

#ifndef initial_h
#define initial_h

#include "generalbc.h"
#include "dictionr.h"

#include "classtype.h"
#include "bcvaltype.h"
#include "valuemodetype.h"

/**
 * Class implementing general initial condition. Initial condition is usually attribute of
 * one or more degrees of freedom (DOFs).
 *
 * //The inherited attribute 'componentArray' is not used. It is replaced with
 * //the more adequate dictionary 'initialValueDictionary', which entries are
 * //referenced by names rather than by indices.
 *
 * One particular DOF (with its physical meaning - for example displacement)
 * can have associated only single initial condition.
 * Initial condition trefore must represent several initial conditions for particular DOF
 * (for example velocity and acceleration of unknown can be prescribed using single
 * initial condition instance). These multiple entries are distinguished by their ValueModeType value.
 * The ValueModeType value is also used as key in initialValueDictionary.
 *
 * Initial conditions apply and should be taken into account
 * only in one particular time step, which number is determined from engineering model
 * giveNumberOfTimeStepWhenIcApply service.
 * If in this time step both boundary condition on unknown and also initial condition
 * on value of this unknown (TotalMode ValueModeType) are prescribed, then
 * always value reported by boundary condition is taken into account.
 * @see EngngModel::giveNumberOfTimeStepWhenIcApply
 */
class InitialCondition : public FEMComponent
{
    /*
     * This class implements an initial condition. An initial condition is usu-
     * ally attribute of one or more degrees of freedom.
     */
private:
    /// Dictionary of initial values.
    Dictionary initialValueDictionary;
    /// Physical maening of bc value
    bcValType valType;


public:
    /**
     * Constructor. Creates initial condition with given number, belonging to given domain.
     * @param n initial condition number
     * @param d domain to which new object will belongs.
     */
    InitialCondition(int i, Domain *d) : FEMComponent(i, d), initialValueDictionary()
    { }
    /// Destructor.
    ~InitialCondition()                { }

    //double  give (char) ;
    /**
     * Returns value of initial condition for given unknown mode (determines whether total or velocity or acceleration
     * mode of unknown is requested).
     * @param mode charcteristic mode of unknown, characteristic type depends on DOF (represent physical meaning).
     * @return value of initial condition for given mode.
     */
    double give(ValueModeType);
    /**
     * Returns receiver load type. It distinguish particular boundary conditions according to
     * their "physical" meaning (like StructuralTemperatureLoadLT, StructuralLoadLT).
     * Derived classes should always overload, default implementation returns UnknownLT value.
     * See cltypes.h file for details.
     */
    bcValType  giveICValType() const { return this->valType; }
    /**
     * Tests if receiver has initial condition for specific unknown-mode.
     * @param u ValueModeType value of unknown
     * @return nonzero if given mode has initial condition, zero otherwise
     */
    int     hasConditionOn(int u);
    /**
     * Tests if receiver has initial condition for specific unknown-mode.
     * @param u ValueModeType value of unknown
     * @return nonzero if given mode has initial condition, zero otherwise
     */
    int     hasConditionOn(ValueModeType type);
    /**
     * Scales the receiver value (determined by ValueModeType) by given value.
     * Typically used in nondimensional analysis to scale down BCs and ICs.
     */
    virtual void scale(ValueModeType, double);



    /// Initializes receiver acording to object description stored in input record.
    IRResultType initializeFrom(InputRecord *ir);
    /// Print receiver state and contents. Usefull for debuging.
    void    printYourself();

    /// Returns classType id of receiver.
    classType    giveClassID() const { return InitialConditionClass; }
    /// Returns class name of the receiver.
    const char *giveClassName() const { return "InitialCondition"; }
};

#endif // initial_h
