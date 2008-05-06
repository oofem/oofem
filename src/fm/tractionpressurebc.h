/* $Header: /home/cvs/bp/oofem/oofemlib/src/boundary.h,v 1.11 2003/04/06 14:08:23 bp Exp $ */
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

//   ********************************
//   *** CLASS BOUNDARY CONDITION ***
//   ********************************


#ifndef tractionpressurebc_h
#define tractionpressurebc_h

#include "boundary.h"

#ifndef __MAKEDEPEND
#include <string.h>
#endif

/**
 * Class implementing prescribed pressure bc due to prescribed tractions (Dirichlet boundary condition on DOF).
 * This boundary condition is usually attribute of one or more degrees of freedom (DOF).
 */
class TractionPressureBC : public BoundaryCondition
{
private:
public:
    /**
     * Constructor. Creates boundary condition with given number, belonging to given domain.
     * @param n boundary condition number
     * @param d domain to which new object will belongs.
     */
    TractionPressureBC(int i, Domain *d) : BoundaryCondition(i, d)
    { }
    /// Destructor
    ~TractionPressureBC()            { }


    /**
     * Returns the value of a prescribed unknown, respecting requested mode for given time.
     * Its physical meaning is determined by corresponding DOF.
     * @param dof determines the dof subjected to receiver bc
     * @param mode unknown char type (if total or incremental value is returned)
     * @return prescribed value of unknown or zero if not prescribed
     */
    virtual double give(Dof *, ValueModeType, TimeStep *);

    /// Initializes receiver acording to object description stored in input record.
    IRResultType initializeFrom(InputRecord *ir);

    /** Setups the input record string of receiver
     *  @param str string to be filled by input record
     *  @param keyword print record keyword (default true)
     */
    virtual int giveInputRecordString(std :: string &str, bool keyword = true);
    /**
     * Scales the receiver according to given value. Typically used in nondimensional analysis to scale down BCs and ICs.
     */
    virtual void scale(double s) { }
};


#endif // tractionpressurebc_h

