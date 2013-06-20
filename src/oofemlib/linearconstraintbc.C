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

#include "linearconstraintbc.h"
#include "classfactory.h"
#include "masterdof.h"
#include "floatmatrix.h"
#include "sparsemtrx.h"
#include "unknownnumberingscheme.h"
#include "loadtimefunction.h"
#include "timestep.h"

namespace oofem {
REGISTER_BoundaryCondition(LinearConstraintBC);


LinearConstraintBC :: LinearConstraintBC(int n, Domain *d) : ActiveBoundaryCondition(n, d)
{
    this->md = new Node(0, domain);
    // this is internal lagrange multiplier uses to enforce the receiver constrain
    this->md->appendDof( new MasterDof( 0, this->md, ( DofIDItem ) ( d->giveNextFreeDofID() ) ) );
}


IRResultType LinearConstraintBC :: initializeFrom(InputRecord *ir)
{
    ActiveBoundaryCondition :: initializeFrom(ir);
    const char *__proc = "initializeFrom";
    IRResultType result;
    rhsLtf = 0.;
    
    IR_GIVE_FIELD(ir, weights, _IFT_LinearConstraintBC_weights);
    IR_GIVE_FIELD(ir, rhs, _IFT_LinearConstraintBC_rhs);
    IR_GIVE_FIELD(ir, dofmans, _IFT_LinearConstraintBC_dofmans);
    IR_GIVE_FIELD(ir, dofs, _IFT_LinearConstraintBC_dofs);
    if(weights.giveSize() != dofmans.giveSize()){
        OOFEM_ERROR3("Size mismatch, weights %d and dofmans %d", weights.giveSize(), dofmans.giveSize());
    }
    IR_GIVE_OPTIONAL_FIELD(ir, weightsLtf, _IFT_LinearConstraintBC_weightsltf);
    IR_GIVE_OPTIONAL_FIELD(ir, rhsLtf, _IFT_LinearConstraintBC_rhsltf);
    return IRRT_OK;
}


void LinearConstraintBC :: giveLocArray(const UnknownNumberingScheme &r_s,  IntArray &locr)
{
    int size = this->weights.giveSize();
    Dof *idof;

    locr.resize(size + 1);
    // assemble location array
    for ( int _i = 1; _i <= size; _i++ ) {
        idof = this->domain->giveDofManager( this->dofmans.at(_i) )->giveDof( this->dofs.at(_i) );
        locr.at(_i) = r_s.giveDofEquationNumber(idof);
    }

    locr.at(size + 1) = r_s.giveDofEquationNumber( md->giveDof(1) );
}


void LinearConstraintBC :: assemble(SparseMtrx *answer, TimeStep *tStep, EquationID eid,
                                    CharType type, const UnknownNumberingScheme &r_s,
                                    const UnknownNumberingScheme &c_s)
{
    int size = this->weights.giveSize();
    FloatMatrix contrib(size + 1, size + 1);
    IntArray locr(size + 1), locc(size + 1);

    for ( int _i = 1; _i <= size; _i++ ) { // loop over dofs
        double factor=1.;
        if(weightsLtf.giveSize()){
            domain->giveLoadTimeFunction(weightsLtf.at(_i))->__at(tStep->giveIntrinsicTime());
        }
        contrib.at(_i, size + 1) = this->weights.at(_i)*factor;
        contrib.at(size + 1, _i) = this->weights.at(_i)*factor;
    }

    contrib.at(size + 1, size + 1) = 0.0; // internal Lagrange multiplier

    this->giveLocArray(r_s, locr);
    this->giveLocArray(c_s, locc);
    answer->assemble(locr, contrib);
}

void LinearConstraintBC :: assembleVector(FloatArray &answer, TimeStep *tStep, EquationID eid,
                                          CharType type, ValueModeType mode,
                                          const UnknownNumberingScheme &s, FloatArray *eNorms)
{
    IntArray loc;
    FloatArray vec;
    double factor=1.;
    if(rhsLtf){
        domain->giveLoadTimeFunction(rhsLtf)->__at(tStep->giveIntrinsicTime());
    }
    this->giveLocArray(s, loc);
    vec.resize( loc.giveSize() );
    vec.at( vec.giveSize() ) = rhs*factor;
    answer.assemble(vec, loc);
}

void LinearConstraintBC :: giveLocationArrays(std :: vector< IntArray > &rows, std :: vector< IntArray > &cols, EquationID eid, CharType type, const UnknownNumberingScheme &r_s, const UnknownNumberingScheme &c_s) {
    rows.resize(1);
    IntArray loc;
    this->giveLocArray(r_s, loc);
    rows [ 0 ] = loc;
}
} //end of oofem namespace
