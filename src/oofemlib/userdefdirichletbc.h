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

#ifndef USERDEFDIRICHLETBC_H_
#define USERDEFDIRICHLETBC_H_

#include "boundarycondition.h"

#define _IFT_UserDefDirichletBC_Name "userdefdirichletbc"
#define _IFT_UserDefDirichletBC_filename "filename"

// forward declare PyObject
// as suggested on the python mailing list
// http://mail.python.org/pipermail/python-dev/2003-August/037601.html
#ifndef PyObject_HEAD
struct _object;
typedef _object PyObject;
#endif

//class PyObject;

/**
 * Class representing user defined Dirichlet boundary conditions. The boundary
 * condition is specified in a Python function with the syntax
 * @code{.py}
 * def giveUserDefBC(coord, iDofNum, time):
 * @endcode
 * where coord is the node coordinate, iDofNum is the dof number, i.e. 1 for D_u, 2 D_v and so on.
 * The time argument is the target time (the time at the end of the time step).
 *
 * The Python function should return the prescribed value in the node.
 *
 * When the boundary condition is created in the input file, the
 * file name of the Python function needs to be specified under
 * the filename keyword. Use lower case letters in the file name!
 *
 * Status: experimental.
 *
 * @date Aug 7, 2013
 * @author: Erik Svenning
 */
namespace oofem {

class PythonInitializer;

class UserDefDirichletBC : public BoundaryCondition
{
protected:
    /// Prescribed values for each resp. dof
    FloatArray values;

    PyObject *mpName;
    PyObject *mpModule;
    PyObject *mpFunc;

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

    virtual double give(Dof *dof, ValueModeType mode, TimeStep *tStep);

    virtual void setPrescribedValue(double s);

    // Overloaded methods:
    virtual bcType giveType() const { return DirichletBT; }
    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual void giveInputRecord(DynamicInputRecord &input);
    virtual void scale(double s);
    virtual const char *giveClassName() const { return "UserDefDirichletBC"; }
    virtual const char *giveInputRecordName() const { return _IFT_UserDefDirichletBC_Name; }
};

} // end namespace oofem

#endif /* USERDEFDIRICHLETBC_H_ */
