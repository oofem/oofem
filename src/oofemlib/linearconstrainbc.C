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

#include "linearconstrainbc.h"
#include "classfactory.h"
#include "masterdof.h"
#include "floatmatrix.h"
#include "sparsemtrx.h"
#include "unknownnumberingscheme.h"

namespace oofem {

REGISTER_BoundaryCondition( LinearConstrainBC );


LinearConstrainBC::LinearConstrainBC (int n, Domain *d) : ActiveBoundaryCondition (n, d) 
{
  this->md = new Node(0, domain);
  // this is internal lagrange multiplier uses to enforce the receiver constrain
  this->md->appendDof(new MasterDof(0, this->md, (DofIDItem)(d->giveNextFreeDofID())));
}



IRResultType LinearConstrainBC::initializeFrom(InputRecord *ir)
    {
        ActiveBoundaryCondition :: initializeFrom(ir);
        const char *__proc = "initializeFrom";
        IRResultType result;

        IR_GIVE_FIELD(ir, weights, _IFT_LinearConstrainBC_weights);
        IR_GIVE_FIELD(ir, dofmans, _IFT_LinearConstrainBC_dofmans);
        IR_GIVE_FIELD(ir, dofs, _IFT_LinearConstrainBC_dofs);
	IR_GIVE_FIELD(ir, rhs, _IFT_LinearConstrainBC_rhs);

        return IRRT_OK;
    }
 

void LinearConstrainBC::giveLocArray (const UnknownNumberingScheme& r_s,  IntArray &locr)
{
  int size = this->weights.giveSize();
  Dof* idof;

  locr.resize(size+1);
  // aseemble location array
  for (int _i=1; _i<=size; _i++) {
    idof = this->domain->giveDofManager(this->dofmans.at(_i))->giveDof(this->dofs.at(_i));
    locr.at(_i) = r_s.giveDofEquationNumber(idof);
  }
  locr.at(size+1) = r_s.giveDofEquationNumber(md->giveDof(1));
}


void LinearConstrainBC::assemble(SparseMtrx *answer, TimeStep *tStep, EquationID eid,
				 CharType type, const UnknownNumberingScheme &r_s, 
				 const UnknownNumberingScheme &c_s) 
{
  int size = this->weights.giveSize();
  FloatMatrix contrib(size+1, size+1);
  IntArray locr (size+1), locc (size+1);

  for (int _i=1; _i<=size; _i++) { // loop over dofs
    contrib.at(_i, size+1) = this->weights.at(_i);
    contrib.at(size+1, _i) = this->weights.at(_i);
  }
  contrib.at(size+1, size+1) = 1.0; // internal multiplier DOF

  this->giveLocArray (r_s,locr);
  this->giveLocArray (c_s,locc);
  answer->assemble(locr, contrib);
}

void LinearConstrainBC::assembleVector(FloatArray &answer, TimeStep *tStep, EquationID eid,
				       CharType type, ValueModeType mode,
				       const UnknownNumberingScheme &s, FloatArray *eNorms) 
{
  IntArray loc;

  this->giveLocArray (s,loc);
  answer.assemble(rhs, loc);
}


}//end of oofem namespace
