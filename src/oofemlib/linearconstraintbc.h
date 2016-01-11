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

#ifndef linearconstraintbc_h
#define linearconstraintbc_h

#include "activebc.h"
#include "floatarray.h"
#include "intarray.h"
#include "chartype.h"
#include "valuemodetype.h"
#include "dofmanager.h"
#include "error.h"

#include <memory>

#define _IFT_LinearConstraintBC_Name   "linearconstraintbc"

///@name Input fields for active boundary condition
//@{
#define _IFT_LinearConstraintBC_weights "weights"
#define _IFT_LinearConstraintBC_weightsfuncs "weightsltf"
#define _IFT_LinearConstraintBC_rhs "rhs"
#define _IFT_LinearConstraintBC_rhsfuncs "rhsltf"
#define _IFT_LinearConstraintBC_dofmans "dofmans"
#define _IFT_LinearConstraintBC_dofs "dofs"
#define _IFT_LinearConstraintBC_lhstype "lhstype"
#define _IFT_LinearConstraintBC_rhstype "rhstype"

//@}

namespace oofem {
/**
 * Class implementing linear constraint on selected DOFs in the form @f$ \sum_i w_i r_i = c @f$,
 * where @f$ r_i @f$ is i-th degree of freedom entering contraint, w_i is corresponding weight,
 * and @f$ c @f$ is given value. The weights and constant can have associated loadtime function,
 * so they can evolve in time. By default, the Function value for all weights and
 * constant is set to 1.0.
 * This boundary condition is introduced as additional stationary condition in energy
 * functional using Lagrange multiplier, which is an additional degree of freedom introduced
 * by this boundary condition.
 */
class OOFEM_EXPORT LinearConstraintBC : public ActiveBoundaryCondition
{
protected:
    FloatArray weights;
    IntArray weightsTf;
    double rhs;
    int rhsTf;
    IntArray dofmans;
    IntArray dofs;
    std :: unique_ptr< DofManager > md;


    // characteristicType of LHS and RHS contributions (this makes this bc trully general, as one can customize, to which
    // characteristic component the contibution will be assembled)
    IntArray lhsType;
    IntArray rhsType;

public:
    LinearConstraintBC(int n, Domain * d);
    /// Destructor.
    virtual ~LinearConstraintBC();

    IRResultType initializeFrom(InputRecord *ir);
    virtual const char *giveInputRecordName() const { return _IFT_LinearConstraintBC_Name; }
    virtual void assemble(SparseMtrx &answer, TimeStep *tStep,
                          CharType type, const UnknownNumberingScheme &r_s, const UnknownNumberingScheme &c_s);
    virtual void assembleVector(FloatArray &answer, TimeStep *tStep,
                                CharType type, ValueModeType mode,
                                const UnknownNumberingScheme &s, FloatArray *eNorms = NULL);

    virtual void giveLocationArrays(std :: vector< IntArray > &rows, std :: vector< IntArray > &cols,
                                    CharType type, const UnknownNumberingScheme &r_s,
                                    const UnknownNumberingScheme &c_s);

    /// Gives the number of internal dof managers.
    virtual int giveNumberOfInternalDofManagers() { return 1; }
    /// Gives an internal dof manager from receiver.
    virtual DofManager *giveInternalDofManager(int i) { return this->md.get(); }

    virtual contextIOResultType saveContext(DataStream &stream, ContextMode mode, void *obj = NULL);
    virtual contextIOResultType restoreContext(DataStream &stream, ContextMode mode, void *obj = NULL);

    virtual const char *giveClassName() const { return "LinearConstraintBC"; }

protected:
    void giveLocArray(const UnknownNumberingScheme &r_s,  IntArray &locr, int &lambdaeq);
};
} //end of oofem namespace
#endif // LinearConstraintBC_h
