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

#ifndef linearconstraintbc_h
#define linearconstraintbc_h

#include "activebc.h"
#include "intarray.h"
#include "equationid.h"
#include "chartype.h"
#include "valuemodetype.h"
#include "dofmanager.h"
#include "error.h"

#define _IFT_LinearConstraintBC_Name   "linearconstraintbc"

///@name Input fields for active boundary condition
//@{
#define _IFT_LinearConstraintBC_weights "weights"
#define _IFT_LinearConstraintBC_weightsltf "weightsltf"
#define _IFT_LinearConstraintBC_rhs "rhs"
#define _IFT_LinearConstraintBC_rhsltf "rhsltf"
#define _IFT_LinearConstraintBC_dofmans "dofmans"
#define _IFT_LinearConstraintBC_dofs "dofs"
//@}

namespace oofem {
/**
 * Class implementing linear constraint on selected DOFs in the form $\sum_i w_ir_i = c$,
 * where r_i is i-th degree of freedom entering contraint, w_i is corresponding weight,
 * and c is given value. The weights and constant can have associated loadtime function,
 * so they can evolve in time. By default, the loadtimefunction value for all weights and 
 * constant is set to 1.0. 
 * This boundary condition is introduced as additional stationary condition in energy 
 * functional using Lagrange multiplier, which is an additional degree of freedom introduced 
 * by this boundary condition. 
 */
class LinearConstraintBC : public ActiveBoundaryCondition
{
protected:
    FloatArray weights;
    IntArray weightsLtf;
    double rhs;
    int rhsLtf;
    IntArray dofmans;
    IntArray dofs;
    DofManager *md;

public:
    LinearConstraintBC(int n, Domain *d);
    /// Destructor.
    virtual ~LinearConstraintBC() { delete this->md; }

    IRResultType initializeFrom(InputRecord *ir);
    virtual const char *giveInputRecordName() const { return _IFT_LinearConstraintBC_Name; }
    virtual classType giveClassID() const { return LinearConstraintClass; }
    virtual void assemble(SparseMtrx *answer, TimeStep *tStep, EquationID eid,
                          CharType type, const UnknownNumberingScheme &r_s, const UnknownNumberingScheme &c_s);
    virtual void assembleVector(FloatArray &answer, TimeStep *tStep, EquationID eid,
                                CharType type, ValueModeType mode,
                                const UnknownNumberingScheme &s, FloatArray *eNorms = NULL);

    virtual void giveLocationArrays(std :: vector< IntArray > &rows, std :: vector< IntArray > &cols, 
				    EquationID eid, CharType type, const UnknownNumberingScheme &r_s, 
				    const UnknownNumberingScheme &c_s);

    /// Gives the number of internal dof managers.
    virtual int giveNumberOfInternalDofManagers() { return 1; }
    /// Gives an internal dof manager from receiver.
    virtual DofManager *giveInternalDofManager(int i) { return this->md; }

    virtual contextIOResultType saveContext(DataStream *stream, ContextMode mode, void *obj = NULL);
    virtual contextIOResultType restoreContext(DataStream *stream, ContextMode mode, void *obj = NULL);
protected:
    void giveLocArray(const UnknownNumberingScheme &r_s,  IntArray &locr, int& lambdaeq);
};
} //end of oofem namespace
#endif // LinearConstraintBC_h
