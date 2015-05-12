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

#include "linearconstraintbc.h"
#include "classfactory.h"
#include "masterdof.h"
#include "floatmatrix.h"
#include "sparsemtrx.h"
#include "unknownnumberingscheme.h"
#include "function.h"
#include "timestep.h"
#include "datastream.h"
#include "contextioerr.h"
#include "node.h"
#include "domain.h"

namespace oofem {
REGISTER_BoundaryCondition(LinearConstraintBC);

LinearConstraintBC :: LinearConstraintBC(int n, Domain *d) : ActiveBoundaryCondition(n, d),
    md( new Node(0, domain) )
{
    // this is internal lagrange multiplier used to enforce the receiver constrain
    // this allocates a new equation related to this constraint
    this->md->appendDof( new MasterDof( this->md.get(), ( DofIDItem ) ( d->giveNextFreeDofID() ) ) );
    this->lhsType.clear();
    this->rhsType.clear();
}


LinearConstraintBC :: ~LinearConstraintBC()
{
}


IRResultType LinearConstraintBC :: initializeFrom(InputRecord *ir)
{
    IRResultType result;
    rhsTf = 0;

    IR_GIVE_FIELD(ir, weights, _IFT_LinearConstraintBC_weights);
    IR_GIVE_FIELD(ir, rhs, _IFT_LinearConstraintBC_rhs);
    IR_GIVE_FIELD(ir, dofmans, _IFT_LinearConstraintBC_dofmans);
    IR_GIVE_FIELD(ir, dofs, _IFT_LinearConstraintBC_dofs);
    if ( weights.giveSize() != dofmans.giveSize() ) {
        OOFEM_WARNING("Size mismatch, weights %d and dofmans %d", weights.giveSize(), dofmans.giveSize());
        return IRRT_BAD_FORMAT;
    }
    IR_GIVE_OPTIONAL_FIELD(ir, weightsTf, _IFT_LinearConstraintBC_weightsfuncs);
    IR_GIVE_OPTIONAL_FIELD(ir, rhsTf, _IFT_LinearConstraintBC_rhsfuncs);

    IR_GIVE_FIELD(ir, lhsType, _IFT_LinearConstraintBC_lhstype);
    IR_GIVE_FIELD(ir, rhsType, _IFT_LinearConstraintBC_rhstype);

    return ActiveBoundaryCondition :: initializeFrom(ir);
}


void LinearConstraintBC :: giveLocArray(const UnknownNumberingScheme &r_s,  IntArray &locr, int &lambda_eq)
{
    int size = this->weights.giveSize();

    locr.resize(size);
    // assemble location array
    for ( int _i = 1; _i <= size; _i++ ) {
        Dof *idof = this->domain->giveDofManager( this->dofmans.at(_i) )->giveDofWithID( this->dofs.at(_i) );
        locr.at(_i) = r_s.giveDofEquationNumber(idof);
    }

    lambda_eq = r_s.giveDofEquationNumber( *md->begin() );
}


void LinearConstraintBC :: assemble(SparseMtrx &answer, TimeStep *tStep,
                                    CharType type, const UnknownNumberingScheme &r_s,
                                    const UnknownNumberingScheme &c_s)
{
    int size = this->weights.giveSize();
    IntArray lambdaeq(1);
    FloatMatrix contrib(size, 1), contribt;
    IntArray locr(size);

    if ( !this->lhsType.contains( ( int ) type ) ) {
        return;
    }
    this->giveLocArray( r_s, locr, lambdaeq.at(1) );

    if ( this->isImposed(tStep) ) {
        for ( int _i = 1; _i <= size; _i++ ) { // loop over dofs
            double factor = 1.;
            if ( weightsTf.giveSize() ) {
                factor = domain->giveFunction( weightsTf.at(_i) )->evaluateAtTime( tStep->giveIntrinsicTime() );
            }
            contrib.at(_i, 1) = this->weights.at(_i) * factor;
        }
        contribt.beTranspositionOf(contrib);

        answer.assemble(lambdaeq, locr, contribt);
        answer.assemble(locr, lambdaeq, contrib);
    } else {
        // the bc is not imposed at specific time step, however in order to make the equation system regular
        // we initialize the allocated equation to the following form 1*labmda = 0, forcing lagrange multiplier
        // of inactive condition to be zero.
        FloatMatrix help(1, 1);
        help.at(1, 1) = 1.0;
        answer.assemble(lambdaeq, lambdaeq, help);
    }
}

void LinearConstraintBC :: assembleVector(FloatArray &answer, TimeStep *tStep,
                                          CharType type, ValueModeType mode,
                                          const UnknownNumberingScheme &s, FloatArray *eNorms)
{
    IntArray loc, lambdaeq(1);
    FloatArray vec(1);
    double factor = 1.;

    if ( !this->rhsType.contains( ( int ) type ) ) {
        return;
    }

    if ( type == InternalForcesVector ) {
        // compute true residual
        int size = this->weights.giveSize();
        Dof *mdof = *md->begin();
        Dof *idof;

        // assemble location array
        for ( int _i = 1; _i <= size; _i++ ) {
            factor = 1.;
            if ( weightsTf.giveSize() ) {
                factor = domain->giveFunction( weightsTf.at(_i) )->evaluateAtTime( tStep->giveIntrinsicTime() );
            }
            idof = this->domain->giveDofManager( this->dofmans.at(_i) )->giveDofWithID( this->dofs.at(_i) );
            if ( s.giveDofEquationNumber(idof) ) {
                answer.at( s.giveDofEquationNumber(idof) ) += mdof->giveUnknown(mode, tStep) * this->weights.at(_i) * factor;
            }
            if ( s.giveDofEquationNumber( mdof ) ) {
                answer.at( s.giveDofEquationNumber( mdof ) ) += idof->giveUnknown(mode, tStep) * this->weights.at(_i) * factor;
            }
        }
    } else if ( type == ExternalForcesVector ) {
        // use rhs value

        if ( rhsTf ) {
            factor = domain->giveFunction(rhsTf)->evaluateAtTime( tStep->giveIntrinsicTime() );
        }
        this->giveLocArray( s, loc, lambdaeq.at(1) );
        vec.at(1) = rhs * factor;
        answer.assemble(vec, lambdaeq);
    }
}

void LinearConstraintBC :: giveLocationArrays(std :: vector< IntArray > &rows, std :: vector< IntArray > &cols, CharType type, const UnknownNumberingScheme &r_s, const UnknownNumberingScheme &c_s)
{
    rows.resize(3);
    cols.resize(3);

    IntArray loc, lambdaeq(1);
    this->giveLocArray( r_s, loc, lambdaeq.at(1) );
    // column block
    rows [ 0 ] = loc;
    cols [ 0 ] = lambdaeq;
    // row block
    cols [ 1 ] = loc;
    rows [ 1 ] = lambdaeq;
    // diagonal enry (some sparse mtrx implementation requaire this)
    rows [ 2 ] = lambdaeq;
    cols [ 2 ] = lambdaeq;
}


contextIOResultType
LinearConstraintBC :: saveContext(DataStream &stream, ContextMode mode, void *obj)
{
    contextIOResultType iores;
    if ( mode & CM_Definition ) {
        if ( ( iores = weights.storeYourself(stream) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }
        if ( ( iores = weightsTf.storeYourself(stream) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }
        if ( ( iores = dofmans.storeYourself(stream) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }
        if ( ( iores = dofs.storeYourself(stream) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }
        if ( !stream.write(rhs) ) {
            THROW_CIOERR(CIO_IOERR);
        }
        if ( !stream.write(rhsTf) ) {
            THROW_CIOERR(CIO_IOERR);
        }
        if ( ( iores = lhsType.storeYourself(stream) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }
        if ( ( iores = rhsType.storeYourself(stream) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }
    }

    if ( ( iores = md->saveContext(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    return CIO_OK;
}


contextIOResultType
LinearConstraintBC :: restoreContext(DataStream &stream, ContextMode mode, void *obj)
{
    contextIOResultType iores;
    if ( mode & CM_Definition ) {
        if ( ( iores = weights.restoreYourself(stream) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }
        if ( ( iores = weightsTf.restoreYourself(stream) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }
        if ( ( iores = dofmans.restoreYourself(stream) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }
        if ( ( iores = dofs.restoreYourself(stream) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }
        if ( !stream.read(rhs) ) {
            THROW_CIOERR(CIO_IOERR);
        }
        if ( !stream.read(rhsTf) ) {
            THROW_CIOERR(CIO_IOERR);
        }
        if ( ( iores = lhsType.restoreYourself(stream) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }
        if ( ( iores = rhsType.restoreYourself(stream) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }
    }

    if ( ( iores = md->restoreContext(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    return CIO_OK;
}
} //end of oofem namespace
