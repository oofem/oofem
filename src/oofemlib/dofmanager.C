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

#include "dofmanager.h"
#include "masterdof.h"
#include "simpleslavedof.h"
#include "timestep.h"
#include "load.h"
#include "floatarray.h"
#include "floatmatrix.h"
#include "intarray.h"
#include "classfactory.h"
#include "datastream.h"
#include "contextioerr.h"
#include "mathfem.h"
#include "dynamicinputrecord.h"
#include "domain.h"
#include "unknownnumberingscheme.h"
#include "entityrenumberingscheme.h"
#include "engngm.h"

namespace oofem {
DofManager :: DofManager(int n, Domain *aDomain) :
    FEMComponent(n, aDomain), dofArray(), loadArray(), partitions()
{
    isBoundaryFlag = false;
    hasSlaveDofs  = false;
    dofidmask = NULL;
    dofTypemap = NULL;
    dofMastermap = NULL;
    dofBCmap = NULL;
    dofICmap = NULL;
    parallel_mode = DofManager_local;
}


DofManager :: ~DofManager()
{
    for ( Dof *dof: dofArray ) {
        delete dof;
    }

    delete dofidmask;
    delete dofTypemap;
    delete dofMastermap;
    delete dofBCmap;
    delete dofICmap;
}


IntArray *DofManager :: giveLoadArray()
// Returns the list containing the number of every nodal loads that act on
// the receiver. If this list does not exist yet, constructs it. This list
// is not to be confused with the load vector.
{
    return & loadArray;
}


void DofManager :: setLoadArray(IntArray &la)
{
    this->loadArray = la;
}


void DofManager :: computeLoadVector(FloatArray &answer, Load *load, CharType type, TimeStep *tStep, ValueModeType mode)
{
    if ( load->giveBCGeoType() != NodalLoadBGT ) {
        OOFEM_ERROR("incompatible load type applied");
    }

    answer.clear();
    if ( type != ExternalForcesVector ) {
        return;
    }

    if ( load->giveDofIDs().giveSize() == 0 ) {
        load->computeComponentArrayAt(answer, tStep, mode);
    } else {
        answer.resize(this->giveNumberOfDofs());
        FloatArray tmp;
        load->computeComponentArrayAt(tmp, tStep, mode);
        answer.assemble(tmp, load->giveDofIDs());
    }
}


Dof *DofManager :: giveDofWithID(int dofID) const
// Returns the degree of freedom of the receiver with 'dofID'.
{
    auto pos = this->findDofWithDofId( ( DofIDItem ) dofID );

#ifdef DEBUG
    if ( pos == this->end() ) {
        OOFEM_ERROR("dof with given DofID does not exists");
    }
#endif

    return *pos;
}


void DofManager :: appendDof(Dof *dof)
// appends a Dof to the end position in dofArray
// because dofArray is not resizeable, the data need
// to be copied out of the dofArray, dofArray deleted
// and again constructed with new Dof
{
#ifdef DEBUG
    // check if dofID of DOF to be added not already present
    if ( this->findDofWithDofId( dof->giveDofID() ) != this->end() ) {
        OOFEM_ERROR("DOF with dofID %d already present (dofman %d)", dof->giveDofID(), this->number);
    }
#endif

    this->dofArray.push_back(dof);
}


void DofManager :: removeDof(DofIDItem id)
{
    int i = 0;
    for ( Dof *dof: *this ) {
        if ( dof->giveDofID() == id ) {
            delete dof;
            this->dofArray.erase( i + this->begin() );
            return;
        }
        i++;
    }
    OOFEM_WARNING("no DOF with dofID %d found", id);
}


bool DofManager :: hasDofID(DofIDItem id)
{
    for ( Dof *dof: *this ) {
        if ( dof->giveDofID() == id ) {
            return true;
        }
    }
    return false;
}


void DofManager :: giveLocationArray(const IntArray &dofIDArry, IntArray &locationArray, const UnknownNumberingScheme &s) const
{
    if ( !hasSlaveDofs ) {
        int size;
        // prevents some size problem when connecting different elements with
        // different number of dofs
        size = dofIDArry.giveSize();
        locationArray.resize(size);
        for ( int i = 1; i <= size; i++ ) {
            auto pos = this->findDofWithDofId( ( DofIDItem ) dofIDArry.at(i) );
            if ( pos == this->end() ) {
                OOFEM_ERROR("incompatible dof (%d) requested", dofIDArry.at(i));
            }

            locationArray.at(i) = s.giveDofEquationNumber( *pos );
        }
    } else {
        IntArray mstrEqNmbrs;
        locationArray.clear();

        for ( int dofid: dofIDArry ) {
            auto pos = this->findDofWithDofId( ( DofIDItem ) dofid );
            if ( pos == this->end() ) {
                OOFEM_ERROR("incompatible dof (%d) requested", dofid);
            }
            (*pos)->giveEquationNumbers(mstrEqNmbrs, s);
            locationArray.followedBy(mstrEqNmbrs);
        }
    }
}


void DofManager :: giveMasterDofIDArray(const IntArray &dofIDArry, IntArray &masterDofIDs) const
{
    if ( !hasSlaveDofs ) {
        masterDofIDs = dofIDArry;
    } else {
        IntArray temp;
        masterDofIDs.clear();

        for ( int dofid: dofIDArry ) {
            auto pos = this->findDofWithDofId( ( DofIDItem ) dofid );
            if ( pos == this->end() ) {
                OOFEM_ERROR("incompatible dof (%d) requested", dofid);
            }
            (*pos)->giveDofIDs(temp);
            masterDofIDs.followedBy(temp);
        }
    }
}


void DofManager :: giveCompleteLocationArray(IntArray &locationArray, const UnknownNumberingScheme &s) const
{
    if ( !hasSlaveDofs ) {
        // prevents some size problem when connecting different elements with
        // different number of dofs
        locationArray.resizeWithValues(0, this->giveNumberOfDofs());
        for ( Dof *dof: *this ) {
            locationArray.followedBy( s.giveDofEquationNumber( dof ) );
        }
    } else {
        IntArray temp;
        locationArray.resize(0);
        for ( Dof *dof: *this ) {
            dof->giveEquationNumbers(temp, s);
            locationArray.followedBy(temp);
        }
    }
}


void DofManager :: giveCompleteMasterDofIDArray(IntArray &dofIDArray) const
{
    if ( !hasSlaveDofs ) {
        dofIDArray.resizeWithValues(0, this->giveNumberOfDofs());
        for ( Dof *dof: *this ) {
            dofIDArray.followedBy( dof->giveDofID() );
        }
    } else {
        IntArray temp;
        for ( Dof *dof: *this ) {
            dof->giveDofIDs(temp);
            dofIDArray.followedBy(temp);
        }
    }
}


std :: vector< Dof* > :: const_iterator  DofManager :: findDofWithDofId(DofIDItem dofID) const
{
    int i = 0;
    for ( Dof *dof: *this ) {
        if ( dof->giveDofID() == dofID ) {
            return this->begin() + i;
        }
        i++;
    }
    return this->end();
}


int DofManager :: giveNumberOfDofs() const
// Returns the number of degrees of freedom of the receiver.
{
    return (int)dofArray.size();
}


void DofManager :: setNumberOfDofs(int _ndofs)
{
    for ( Dof *dof: *this ) {
        delete dof;
    }
    this->dofArray.assign(_ndofs, NULL);
}


void DofManager :: askNewEquationNumbers(TimeStep *tStep)
{
    for ( Dof *dof: *this ) {
        dof->askNewEquationNumber(tStep);
    }
}


int DofManager :: giveNumberOfPrimaryMasterDofs(const IntArray &dofIDArray) const
{
    if ( !hasSlaveDofs ) {
        return dofIDArray.giveSize();
    }

    int answer = 0;

    for ( int dofid: dofIDArray ) {
        auto pos = this->findDofWithDofId((DofIDItem)dofid);
#ifdef DEBUG
        if ( pos == this->end() ) {
            OOFEM_ERROR("Dof with ID %d doesn't exist", dofid);
        }
#endif
        answer += (*pos)->giveNumberOfPrimaryMasterDofs();
    }

    return answer;
}


IRResultType
DofManager :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                 // Required by IR_GIVE_FIELD macro

    delete dofidmask;
    dofidmask = NULL;
    delete dofTypemap;
    dofTypemap = NULL;
    delete dofMastermap;
    dofMastermap = NULL;
    delete dofBCmap;
    dofBCmap = NULL;
    delete dofICmap;
    dofICmap = NULL;

    IntArray dofIDArry;
    IntArray ic, masterMask, dofTypeMask;

    loadArray.clear();
    IR_GIVE_OPTIONAL_FIELD(ir, loadArray, _IFT_DofManager_load);

    ///@todo This is unnecessary, we should just check if user has supplied a dofidmask field or not and just drop "numberOfDofs"). It is left for now because it will give lots of warnings otherwise, but it is effectively ignored.
    int dummy;
    IR_GIVE_OPTIONAL_FIELD(ir, dummy, "ndofs");

    if ( ir->hasField(_IFT_DofManager_dofidmask) ) {
        IR_GIVE_FIELD(ir, dofIDArry, _IFT_DofManager_dofidmask);
        this->dofidmask = new IntArray(dofIDArry);
    } else {
        dofIDArry = domain->giveDefaultNodeDofIDArry();
    }

    mBC.clear();
    IR_GIVE_OPTIONAL_FIELD(ir, mBC, _IFT_DofManager_bc);

    ic.clear();
    IR_GIVE_OPTIONAL_FIELD(ir, ic, _IFT_DofManager_ic);

    // reads master mask - in this array are numbers of master dofManagers
    // to which are connected dofs in receiver.
    // if master mask index is zero then dof is created as master (i.e., having own equation number)
    // othervise slave dof connected to master DofManager is created.
    // by default if masterMask is not specifyed, all dofs are created as masters.
    dofTypeMask.clear(); // termitovo
    IR_GIVE_OPTIONAL_FIELD(ir, dofTypeMask, _IFT_DofManager_doftypemask);

    // read boundary flag
    if ( ir->hasField(_IFT_DofManager_boundaryflag) ) {
        isBoundaryFlag = true;
    }


    partitions.clear();
    IR_GIVE_OPTIONAL_FIELD(ir, partitions, _IFT_DofManager_partitions);

    if ( ir->hasField(_IFT_DofManager_sharedflag) ) {
        parallel_mode = DofManager_shared;
    } else if ( ir->hasField(_IFT_DofManager_remoteflag) ) {
        parallel_mode = DofManager_remote;
    } else if ( ir->hasField(_IFT_DofManager_nullflag) ) {
        parallel_mode = DofManager_null;
    } else {
        parallel_mode = DofManager_local;
    }

    // in parallel mode,  slaves are allowed, because ((Dr. Rypl promissed)
    // masters have to be in same partition as slaves. They can be again Remote copies.


    bool hasIc = ic.giveSize() != 0;
    bool hasBc = mBC.giveSize() != 0;
    bool hasTypeinfo = dofTypeMask.giveSize() != 0;

    ///@todo This should eventually be removed, still here to preserve backwards compatibility:
    if ( ( hasIc || hasBc || hasTypeinfo ) && !this->dofidmask ) {
        this->dofidmask = new IntArray(dofIDArry);
    }

    // check sizes
    if ( hasBc ) {
        if ( mBC.giveSize() != dofIDArry.giveSize() ) {
            OOFEM_ERROR("bc size mismatch. Size is %d and need %d", mBC.giveSize(), dofIDArry.giveSize());
        }
        this->dofBCmap = new std :: map< int, int >();
        for ( int i = 1; i <= mBC.giveSize(); ++i ) {
            if ( mBC.at(i) > 0 ) {
                ( * this->dofBCmap ) [ dofIDArry.at(i) ] = mBC.at(i);
            }
        }
    }

    if ( hasIc ) {
        if ( ic.giveSize() != dofIDArry.giveSize() ) {
            OOFEM_ERROR("ic size mismatch. Size is %d and need %d", ic.giveSize(), dofIDArry.giveSize());
        }
        this->dofICmap = new std :: map< int, int >();
        for ( int i = 1; i <= ic.giveSize(); ++i ) {
            if ( ic.at(i) > 0 ) {
                ( * this->dofICmap ) [ dofIDArry.at(i) ] = ic.at(i);
            }
        }
    }

    if ( hasTypeinfo ) {
        if ( dofTypeMask.giveSize() != dofIDArry.giveSize() ) {
            OOFEM_ERROR("dofTypeMask size mismatch. Size is %d and need %d", dofTypeMask.giveSize(), dofIDArry.giveSize());
        }
        this->dofTypemap = new std :: map< int, int >();
        for ( int i = 1; i <= dofTypeMask.giveSize(); ++i ) {
            if ( dofTypeMask.at(i) != DT_master ) {
                ( * this->dofTypemap ) [ dofIDArry.at(i) ] = dofTypeMask.at(i);
            }
        }
        // For simple slave dofs:
        if ( dofTypeMask.contains(DT_simpleSlave) ) {
            IR_GIVE_FIELD(ir, masterMask, _IFT_DofManager_mastermask);
            if ( masterMask.giveSize() != dofIDArry.giveSize() ) {
                OOFEM_ERROR("mastermask size mismatch");
            }
            this->dofMastermap = new std :: map< int, int >();
            for ( int i = 1; i <= masterMask.giveSize(); ++i ) {
                if ( masterMask.at(i) > 0 ) {
                    ( * this->dofMastermap ) [ dofIDArry.at(i) ] = masterMask.at(i);
                }
            }
        }
    }

    return IRRT_OK;
}


void DofManager :: giveInputRecord(DynamicInputRecord &input)
{
    FEMComponent :: giveInputRecord(input);

    IntArray mbc, dofids;
    // Ignore here mBC and dofidmask as they may not correspod to actual state.
    // They just decribe the state at init, but after some dofs may have been
    // added dynamically (xfem, subdivision, etc).
    for ( Dof *dof: *this ) {
      dofids.followedBy(dof->giveDofID(),3);
      if (dof->giveBcId()) mbc.followedBy(dof->giveBcId(),3);
      else mbc.followedBy(0,3);
    }
    input.setField(mbc, _IFT_DofManager_bc);
    input.setField(dofids, _IFT_DofManager_dofidmask);


    if ( this->dofTypemap ) {
        IntArray typeMask( this->dofidmask->giveSize() );
        for ( int i = 1; i <= dofidmask->giveSize(); ++i ) {
            typeMask.at(i) = ( * this->dofTypemap ) [ dofidmask->at(i) ];
        }
        input.setField(typeMask, _IFT_DofManager_doftypemask);
    }

    if ( this->dofMastermap ) {
        IntArray masterMask( this->dofidmask->giveSize() );
        for ( int i = 1; i <= dofidmask->giveSize(); ++i ) {
            masterMask.at(i) = ( * this->dofMastermap ) [ dofidmask->at(i) ];
        }
        input.setField(masterMask, _IFT_DofManager_mastermask);
    }

    if ( isBoundaryFlag ) {
        input.setField(_IFT_DofManager_boundaryflag);
    }


    if ( this->partitions.giveSize() > 0 ) {
        input.setField(this->partitions, _IFT_DofManager_partitions);
    }

    if ( parallel_mode == DofManager_shared ) {
        input.setField(_IFT_DofManager_sharedflag);
    } else if ( parallel_mode == DofManager_remote ) {
        input.setField(_IFT_DofManager_remoteflag);
    } else if ( parallel_mode == DofManager_null ) {
        input.setField(_IFT_DofManager_nullflag);
    }
}


void DofManager :: printOutputAt(FILE *stream, TimeStep *tStep)
{
    EngngModel *emodel = this->giveDomain()->giveEngngModel();

    fprintf( stream, "%-8s%8d (%8d):\n", this->giveClassName(), this->giveLabel(), this->giveNumber() );
    for ( Dof *dof: *this ) {
        emodel->printDofOutputAt(stream, dof, tStep);
    }
}


void DofManager :: printYourself()
// Prints the receiver on screen.
{
    printf("DofManager %d\n", number);
    for ( Dof *dof: *this ) {
        dof->printYourself();
    }

    loadArray.printYourself();
    printf("\n");
}


void DofManager :: updateYourself(TimeStep *tStep)
// Updates the receiver at end of step.
{
    for ( Dof *dof: *this ) {
        dof->updateYourself(tStep);
    }
}


contextIOResultType DofManager :: saveContext(DataStream &stream, ContextMode mode, void *obj)
//
// saves full node context (saves state variables, that completely describe
// current state)
//
{
    int _val;
    contextIOResultType iores;

    if ( ( iores = FEMComponent :: saveContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }


    int numberOfDofs = this->giveNumberOfDofs();
    if ( !stream.write(numberOfDofs) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    // store dof types
    for ( Dof *dof: *this ) {
        _val = dof->giveDofType();
        if ( !stream.write(_val) ) {
            THROW_CIOERR(CIO_IOERR);
        }
    }

    // store dof types
    for ( Dof *dof: *this ) {
        _val = dof->giveDofID();
        if ( !stream.write(_val) ) {
            THROW_CIOERR(CIO_IOERR);
        }
    }

    if ( mode & CM_Definition ) {
        if ( ( iores = loadArray.storeYourself(stream) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }

        if ( !stream.write(isBoundaryFlag) ) {
            THROW_CIOERR(CIO_IOERR);
        }

        if ( !stream.write(hasSlaveDofs) ) {
            THROW_CIOERR(CIO_IOERR);
        }

        if ( !stream.write(globalNumber) ) {
            THROW_CIOERR(CIO_IOERR);
        }

        _val = ( int ) parallel_mode;
        if ( !stream.write(_val) ) {
            THROW_CIOERR(CIO_IOERR);
        }

        if ( ( iores = partitions.storeYourself(stream) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }
    }

    for ( Dof *dof: *this ) {
        if ( ( iores = dof->saveContext(stream, mode, obj) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }
    }

    return CIO_OK;
}


contextIOResultType DofManager :: restoreContext(DataStream &stream, ContextMode mode, void *obj)
//
// restores full node context (saves state variables, that completely describe
// current state)
//
{
    contextIOResultType iores;

    if ( ( iores = FEMComponent :: restoreContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    int _numberOfDofs;
    if ( !stream.read(_numberOfDofs) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    IntArray dtypes(_numberOfDofs);
    // restore dof types
    for ( int i = 1; i <= _numberOfDofs; i++ ) {
        if ( !stream.read(dtypes.at(i)) ) {
            THROW_CIOERR(CIO_IOERR);
        }
    }
    
    IntArray dofids(_numberOfDofs);
    // restore dof ids
    for ( int i = 1; i <= _numberOfDofs; i++ ) {
        if ( !stream.read(dofids.at(i)) ) {
            THROW_CIOERR(CIO_IOERR);
        }
    }

    // allocate new ones
    for ( auto &d: dofArray) { delete d; } ///@todo Smart pointers would be nicer here
    dofArray.clear();
    for ( int i = 0; i < _numberOfDofs; i++ ) {
        Dof *dof = classFactory.createDof( ( dofType ) dtypes(i), (DofIDItem)dofids(i), this );
        this->appendDof(dof);
    }

    if ( mode & CM_Definition ) {
        if ( ( iores = loadArray.restoreYourself(stream) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }

        if ( !stream.read(isBoundaryFlag) ) {
            THROW_CIOERR(CIO_IOERR);
        }

        if ( !stream.read(hasSlaveDofs) ) {
            THROW_CIOERR(CIO_IOERR);
        }

        if ( !stream.read(globalNumber) ) {
            THROW_CIOERR(CIO_IOERR);
        }

        int _val;
        if ( !stream.read(_val) ) {
            THROW_CIOERR(CIO_IOERR);
        }

        parallel_mode = ( dofManagerParallelMode ) _val;
        if ( ( iores = partitions.restoreYourself(stream) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }
    }

    for ( Dof *dof: *this ) {
        if ( ( iores = dof->restoreContext(stream, mode, obj) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }
    }

    return CIO_OK;
}


void DofManager :: giveUnknownVector(FloatArray &answer, const IntArray &dofIDArry, ValueModeType mode, TimeStep *tStep, bool padding)
{
    answer.resize( dofIDArry.giveSize() );
    if ( dofIDArry.giveSize() == 0 ) return;

    int k = 0;
    for ( auto &dofid: dofIDArry ) {
        auto pos = this->findDofWithDofId( ( DofIDItem ) dofid );
        if ( pos == this->end() ) {
            if ( padding ) {
                answer.at(++k) = 0.;
                continue;
            } else {
                continue;
            }
        }
        answer.at(++k) = (*pos)->giveUnknown(mode, tStep);
    }
    answer.resizeWithValues(k);

    // Transform to global c.s.
    FloatMatrix L2G;
    if ( this->computeL2GTransformation(L2G, dofIDArry) ) {
        answer.rotatedWith(L2G, 'n');
    }
}


void DofManager :: giveUnknownVector(FloatArray &answer, const IntArray &dofIDArry,
                                     PrimaryField &field, ValueModeType mode, TimeStep *tStep, bool padding)
{
    answer.resize( dofIDArry.giveSize() );

    int k = 0;
    for ( auto &dofid: dofIDArry ) {
        auto pos = this->findDofWithDofId( ( DofIDItem ) dofid );
        if ( pos == this->end() ) {
            if ( padding ) {
                answer.at(++k) = 0.;
                continue;
            } else {
                continue;
            }
        }
        answer.at(++k) = (*pos)->giveUnknown(field, mode, tStep);
    }
    answer.resizeWithValues(k);

    // Transform to global c.s.
    FloatMatrix L2G;
    if ( this->computeL2GTransformation(L2G, dofIDArry) ) {
        answer.rotatedWith(L2G, 'n');
    }
}


void DofManager :: giveCompleteUnknownVector(FloatArray &answer, ValueModeType mode, TimeStep *tStep)
{
    int i = 1;
    answer.resize(this->giveNumberOfDofs());
    for ( Dof *dof: *this ) {
        answer.at(i) = dof->giveUnknown(mode, tStep);
        i++;
    }
}


void DofManager :: givePrescribedUnknownVector(FloatArray &answer, const IntArray &dofIDArry,
                                               ValueModeType mode, TimeStep *tStep)
{
    answer.resize(dofIDArry.giveSize());

    int j = 1;
    for ( int dofid: dofIDArry ) {
        answer.at(j++) = this->giveDofWithID( dofid )->giveBcValue(mode, tStep);
    }

    // Transform to global c.s.
    FloatMatrix L2G;
    if ( this->computeL2GTransformation(L2G, dofIDArry) ) {
        answer.rotatedWith(L2G, 'n');
    }
}


void DofManager :: giveUnknownVectorOfType(FloatArray &answer, UnknownType ut, ValueModeType mode, TimeStep *tStep)
{
    int k = 1;
    FloatArray localVector(3);
    IntArray dofIDArry(3);

    // This is a bit cumbersome. I first construct the local vector, which might have a odd order, e.g [D_w, D_u], which is later added to the global vector "answer"
    // I also store the dof id's, so that I can construct the local 2 global transformation afterwards (if its necessary). / Mikael
    for ( Dof *d: *this ) {
        double val = d->giveUnknown(mode, tStep);
        if ( ut == DisplacementVector || ut == EigenVector ) { // Just treat eigenvectors as displacement vectors (they are redundant)
            if      ( d->giveDofID() == D_u ) {
                dofIDArry.at(k) = D_u;
                localVector.at(k) = val;
                k++;
            } else if ( d->giveDofID() == D_v ) {
                dofIDArry.at(k) = D_v;
                localVector.at(k) = val;
                k++;
            } else if ( d->giveDofID() == D_w ) {
                dofIDArry.at(k) = D_w;
                localVector.at(k) = val;
                k++;
            }
        } else if ( ut == VelocityVector ) {
            if      ( d->giveDofID() == V_u ) {
                dofIDArry.at(k) = V_u;
                localVector.at(k) = val;
                k++;
            } else if ( d->giveDofID() == V_v ) {
                dofIDArry.at(k) = V_v;
                localVector.at(k) = val;
                k++;
            } else if ( d->giveDofID() == V_w ) {
                dofIDArry.at(k) = V_w;
                localVector.at(k) = val;
                k++;
            }
        } else {
            OOFEM_ERROR("Can't produce vector for unknown type: %d", ut);
        }
    }

    FloatMatrix L2G;
    if ( this->computeL2GTransformation(L2G, dofIDArry) ) {
        // Transform to global c.s.
        answer.beProductOf(L2G, localVector);
    } else {
        // No local c.s, just copy the values to respective index;
        answer.resize(3);
        answer.zero();
        for ( int i = 1; i <= k; i++ ) {
            if ( dofIDArry.at(i) == D_u || dofIDArry.at(i) == V_u ) {
                answer.at(1) = localVector.at(i);
            } else if ( dofIDArry.at(i) == D_v || dofIDArry.at(i) == V_v ) {
                answer.at(2) = localVector.at(i);
            } else if ( dofIDArry.at(i) == D_w || dofIDArry.at(i) == V_w ) {
                answer.at(3) = localVector.at(i);
            }
        }
    }
}


bool DofManager :: hasAnySlaveDofs()
{
    for ( Dof *dof: *this ) {
        if ( !dof->isPrimaryDof() ) {
            return true;
        }
    }

    return false;
}


bool DofManager :: giveMasterDofMans(IntArray &masters)
{
    IntArray _dof_masters;
    bool answer = false;

    masters.clear();
    for ( Dof *dof: *this ) {
        if ( !dof->isPrimaryDof() ) {
            answer = true;
            dof->giveMasterDofManArray(_dof_masters);
            for ( int j = 1; j <= _dof_masters.giveSize(); j++ ) {
                masters.insertSortedOnce(_dof_masters.at(j), 2);
            }
        }
    }

    return answer;
}


void DofManager :: postInitialize()
{
    hasSlaveDofs = false;
    for ( Dof *dof: *this ) {
        if ( !dof->isPrimaryDof() ) {
            hasSlaveDofs = true;
            continue;
        }
    }
}


bool DofManager :: computeM2GTransformation(FloatMatrix &answer, const IntArray &dofMask)
// computes transformation matrix of receiver.
// transformation should include transformation from global cs to nodal cs,
// as well as further necessary transformations (for example in case
// rigid arms this must include transformation to master dofs).
{
    FloatMatrix L2G, M2L;

    bool hasL2G = computeL2GTransformation(L2G, dofMask);
    bool hasM2L = computeM2LTransformation(M2L, dofMask);

    if ( !hasL2G && !hasM2L ) {
        answer.clear();
        return false;
    } else if ( hasL2G && hasM2L ) {
        answer.beProductOf(L2G, M2L);
    } else if ( hasL2G ) {
        answer = L2G;
    } else {
        answer = M2L;
    }
    return true;
}


bool DofManager :: computeL2GTransformation(FloatMatrix &answer, const IntArray &dofIDArry)
{
    return false;
}


bool DofManager :: computeM2LTransformation(FloatMatrix &answer, const IntArray &dofIDArry)
{
    if ( !this->hasAnySlaveDofs() ) {
        return false;
    }

    FloatArray mstrContrs;

    if ( dofIDArry.isEmpty() ) {
        ///@todo I don't think this should be called like this since it relies on the order of the dofs.
        int cols = 0;
        for ( Dof *dof: *this ) {
            cols += dof->giveNumberOfPrimaryMasterDofs();
        }
        answer.resize( this->giveNumberOfDofs(), cols );
        answer.zero();

        int k = 1, i = 1;
        for ( Dof *dof: *this ) {
            dof->computeDofTransformation(mstrContrs);
            answer.copySubVectorRow(mstrContrs, i, k);
            k += mstrContrs.giveSize();
            i++;
        }
    } else {
        answer.resize( dofIDArry.giveSize(), giveNumberOfPrimaryMasterDofs(dofIDArry) );
        answer.zero();

        int k = 1;
        for ( int i = 1; i <= dofIDArry.giveSize(); i++ ) {
            this->giveDofWithID(dofIDArry.at(i))->computeDofTransformation(mstrContrs);
            answer.copySubVectorRow(mstrContrs, i, k);
            k += mstrContrs.giveSize();
        }
    }
    return true;
}


bool DofManager :: requiresTransformation()
{
    return this->hasAnySlaveDofs();
}


void DofManager :: updateLocalNumbering(EntityRenumberingFunctor &f)
{
    //update masterNode numbering
    if ( this->dofMastermap ) {
        for ( auto & mapper: *this->dofMastermap ) {
            mapper.second = f( mapper.second, ERS_DofManager );
        }
    }

    for ( Dof *dof: *this ) {
        dof->updateLocalNumbering(f);
    }
}


void DofManager :: mergePartitionList(IntArray &_p)
{
    // more optimized version can be made requiring sorted partition list of receiver
    int size = _p.giveSize();
    for ( int i = 1; i <= size; i++ ) {
        partitions.insertOnce( _p.at(i) );
    }
}


int DofManager :: givePartitionsConnectivitySize()
{
    int n = partitions.giveSize();
    int myrank = this->giveDomain()->giveEngngModel()->giveRank();
    if ( partitions.findFirstIndexOf(myrank) ) {
        return n;
    } else {
        return n + 1;
    }
}


bool DofManager :: isLocal()
{
    if ( parallel_mode == DofManager_local ) {
        return true;
    }

    if ( parallel_mode == DofManager_shared ) {
        // determine if problem is the lowest one sharing the dofman; if yes the receiver is responsible to
        // deliver number
        int n = partitions.giveSize();
        int myrank = this->giveDomain()->giveEngngModel()->giveRank();
        int minrank = myrank;

        for ( int j = 1; j <= n; j++ ) {
            minrank = min( minrank, partitions.at(j) );
        }

        if ( minrank == myrank ) {
            return true;
        }
    }

    return false;
}

} // end namespace oofem
