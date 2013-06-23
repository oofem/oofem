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

#ifndef linearconstrainbc_h
#define linearconstrainbc_h

#include "activebc.h"
#include "intarray.h"
#include "equationid.h"
#include "chartype.h"
#include "valuemodetype.h"
#include "dofmanager.h"
#include "error.h"

#define _IFT_LinearConstrainBC_Name   "linearconstrainbc"

///@name Input fields for active boundary condition
//@{
#define _IFT_LinearConstrainBC_weights "weights"
#define _IFT_LinearConstrainBC_dofmans "dofmans"
#define _IFT_LinearConstrainBC_dofs "dofs"
#define _IFT_LinearConstrainBC_rhs "rhs"
//@}

namespace oofem {
/**
 * Abstract base class for all active boundary conditions.
 * Design of active boundary conditions are subject to change.
 */
class LinearConstrainBC : public ActiveBoundaryCondition
{
 protected:
  FloatArray weights;
  FloatArray rhs;
  IntArray dofmans;
  IntArray dofs;
  DofManager* md;

 public:
  LinearConstrainBC (int n, Domain *d);
  /// Destructor.
  virtual ~LinearConstrainBC() { delete this->md; }

  IRResultType initializeFrom(InputRecord *ir);
  virtual const char *giveInputRecordName() const { return _IFT_LinearConstrainBC_Name; }
  virtual void assemble(SparseMtrx *answer, TimeStep *tStep, EquationID eid,
			CharType type, const UnknownNumberingScheme &r_s, const UnknownNumberingScheme &c_s);
  virtual void assembleVector(FloatArray &answer, TimeStep *tStep, EquationID eid,
			      CharType type, ValueModeType mode,
			      const UnknownNumberingScheme &s, FloatArray *eNorms = NULL);

  /// Gives the number of internal dof managers.
  virtual int giveNumberOfInternalDofManagers() { return 1; }
  /// Gives an internal dof manager from receiver.
  virtual DofManager *giveInternalDofManager(int i) { return this->md; }
 protected:
  void giveLocArray (const UnknownNumberingScheme& r_s,  IntArray &locr);

  
};

}//end of oofem namespace
#endif // linearconstrainbc_h
