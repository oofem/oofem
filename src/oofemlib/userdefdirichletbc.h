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

/*
 * userdefdirichletbc.h
 *
 *  Created on: Aug 7, 2013
 *
 *  Class representing user defined Dirichlet boundary conditions. The boundary
 *  condition is specified in a Python function with the syntax
 *
 *  def giveUserDefBC(iX, iY, iZ, iDofNum):
 *
 *  where (iX, iY, iZ) are the node coordinates and
 *  iDofNum is the dof number, i.e. 1 for x-direction,
 *  2 for y-direction and so on.
 *
 *  The Python function should return the prescribed value in the node.
 *
 *	When the boundary condition is created in the input file, the
 *	file name of the Python function needs to be specified under
 *	the filename keyword. Use lower case letters in the file name!
 *
 *	A separate PythonInitializer class is used to initialize and finalize
 *	the Python interpreter, see below.
 *
 *	Status: experimental.
 *
 * 	@author: Erik Svenning
 */

#ifndef USERDEFDIRICHLETBC_H_
#define USERDEFDIRICHLETBC_H_

#include "boundarycondition.h"

#define _IFT_UserDefDirichletBC_Name "userdefdirichletbc"

#define _IFT_UserDefDirichletBC_filename "filename"

#ifdef USERDEFDIRICHLETBC

// forward declare PyObject
// as suggested on the python mailing list
// http://mail.python.org/pipermail/python-dev/2003-August/037601.html
#ifndef PyObject_HEAD
struct _object;
typedef _object PyObject;
#endif

//class PyObject;
#endif

namespace oofem {

class PythonInitializer;

class UserDefDirichletBC : public BoundaryCondition
{
protected:
    /// Prescribed values for each resp. dof
    FloatArray values;

    /*
     * 	To run a Python script, the functions Py_Initialize() and Py_Finalize() must be called.
     * 	Calling these functions every time UserDefDirichletBC::give() is called becomes very
     * 	expensive. A solution to the problem is to put Py_Initialize() in the contructor of a
     * 	separate class PythonInitializer and to put Py_Finalize() in the destructor of the same
     * 	class. Creating a static PythonInitializer object ensures that Py_Initialize() and
     * 	Py_Finalize() are called exactly once, which is what we want.
     */
#ifdef USERDEFDIRICHLETBC
    static PythonInitializer mPythonInitializer;

    PyObject *mpName;
    PyObject *mpModule;
    PyObject *mpFunc;
#endif

    std::string mFileName;

public:
    /**
     * Constructor. Creates boundary condition with given number, belonging to given domain.
     * @param i Boundary condition number.
     * @param d Domain to which new object will belongs.
     */
    UserDefDirichletBC(int i, Domain *d);

    /// Destructor
    virtual ~UserDefDirichletBC();

    /**
     * Returns the value of a prescribed unknown, respecting requested mode for given time.
     * Its physical meaning is determined by corresponding DOF.
     * This function should only be used if the BC is imposed.
     * @see isImposed
     * @param dof Determines the dof subjected to receiver BC.
     * @param mode Unknown char type (if total or incremental value is returned).
     * @param tStep Time step to give value for.
     * @return Prescribed value of unknown or zero if not prescribed.
     */
    virtual double give(Dof *dof, ValueModeType mode, TimeStep *tStep);

    /**
     * Set prescribed value at the input record string of receiver
     * @param s prescribed value
     * @todo This function isn't as meaningful anymore. Possibly keep it if we change it to a vector. No inheriting b.c.s can overload this in a meaningful way.
     */
    virtual void setPrescribedValue(double s);

    // Overloaded methods:
    virtual bcType giveType() const { return DirichletBT; }
    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual void giveInputRecord(DynamicInputRecord &input);
    virtual void scale(double s);
    virtual const char *giveClassName() const { return "UserDefDirichletBC"; }
    virtual const char *giveInputRecordName() const { return _IFT_UserDefDirichletBC_Name; }
    virtual classType giveClassID() const { return UserDefDirichletBCClass; }
};

class PythonInitializer {
public:
	PythonInitializer();
	~PythonInitializer();
};

} // end namespace oofem

#endif /* USERDEFDIRICHLETBC_H_ */
