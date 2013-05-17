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

#include "dofmanager.h"
#include "masterdof.h"
#include "slavedof.h"
#include "simpleslavedof.h"
#include "activedof.h"
#include "timestep.h"
#include "load.h"
#include "floatarray.h"
#include "floatmatrix.h"
#include "intarray.h"
#include "classfactory.h"
#include "datastream.h"
#include "contextioerr.h"
#include "mathfem.h"

namespace oofem {
DofManager :: DofManager(int n, Domain *aDomain) :
    FEMComponent(n, aDomain), loadArray()
    // Constructor. Creates a node with number n, belonging to aDomain.
{
    numberOfDofs  = 0;
    dofArray      = NULL;
    isBoundaryFlag = false;
    hasSlaveDofs  = false;
    dofidmask = NULL;
    dofTypemap = NULL;
    dofMastermap = NULL;
    dofBCmap = NULL;
    dofICmap = NULL;
#ifdef __PARALLEL_MODE
    partitions.resize(0);
    parallel_mode = DofManager_local;
#endif
}


DofManager :: ~DofManager()
// Destructor.
{
    int i = numberOfDofs;

    if ( numberOfDofs ) {
        while ( i-- ) {
            delete dofArray [ i ];
        }

        delete[] dofArray;
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


void DofManager :: computeLoadVectorAt(FloatArray &answer, TimeStep *stepN, ValueModeType mode)
// Computes the vector of the nodal loads of the receiver.
{
    FloatArray contribution;

    if ( this->giveLoadArray()->isEmpty() ) {
        answer.resize(0);
        return;
    } else {
        answer.resize(0);
        int nLoads = loadArray.giveSize();     // the node may be subjected
        for ( int i = 1; i <= nLoads; i++ ) {   // to more than one load
            int n = loadArray.at(i);
            Load *loadN = domain->giveLoad(n);
            computeLoadVector(contribution, loadN, ExternalForcesVector, stepN, mode);
            answer.add(contribution);
        }
    }
}


void DofManager :: computeLoadVector(FloatArray &answer, Load *load, CharType type, TimeStep *stepN, ValueModeType mode)
{
    if ( load->giveBCGeoType() != NodalLoadBGT ) {
        _error("computeLoadVector: incompatible load type applied");
    }
    load->computeComponentArrayAt(answer, stepN, mode);
}


Dof *DofManager :: giveDof(int i) const
// Returns the i-th degree of freedom of the receiver.
{
#ifdef DEBUG
    if ( !dofArray ) {
        _error("giveDof: dof is not defined");
    }
#endif

    return dofArray [ i - 1 ];
}


Dof *DofManager :: giveDofWithID(int dofID) const
// Returns the degree of freedom of the receiver with 'dofID'.
{
    int indx = this->findDofWithDofId( ( DofIDItem ) dofID );

#ifdef DEBUG
    if ( !indx ) {
        _error("giveDofWithID: dof with given DofID does not exists");
    }
#endif

    return dofArray [ indx - 1 ];
}


void DofManager :: setDof(int i, Dof *dof)
{
    if ( i <= numberOfDofs ) {
        if ( dofArray [ i - 1 ] ) {
            delete dofArray [ i - 1 ];
        }

        dofArray [ i - 1 ] = dof;
    } else {
        _error("setDof: DOF index out of range");
    }
}


void DofManager :: appendDof(Dof *dof)
// appends a Dof to the end position in dofArray
// because dofArray is not resizeable, the data need
// to be copied out of the dofArray, dofArray deleted
// and again constructed with new Dof
{
#ifdef DEBUG
    // check if dofID of DOF to be added not already present
    if (this->findDofWithDofId(dof->giveDofID()) != 0) {
        _error3("DofManager::addDof: DOF with dofID %d already present (dofman %d)", dof->giveDofID(), this->number);
    }
#endif

    Dof **dofArray_new = new Dof * [ numberOfDofs + 1 ];
    for ( int j = 0; j < numberOfDofs; j++ ) {
        dofArray_new[j] = dofArray[j];
    }
    dofArray_new[numberOfDofs] = dof;
    delete [] dofArray;
    this->dofArray = dofArray_new;
    numberOfDofs++;
}


void DofManager :: removeDof(DofIDItem id)
{
    int position = this->findDofWithDofId(id);
    if (position) { // dof id found
        Dof **dofArray_new = new Dof * [ numberOfDofs - 1 ];
        for ( int j = 0; j < position-1; j++ ) {
            dofArray_new[j] = dofArray[j];
        }

        for ( int j = position; j < numberOfDofs; j++ ) {
            dofArray_new[j-1] = dofArray[j];
        }
        // delete dof
        delete dofArray[position-1];
        delete []dofArray;
        dofArray = dofArray_new;
        numberOfDofs--;
    } else {
        _error2("removeDof::no DOF with dofID %d found", id);
    }
}


bool DofManager :: hasDofID(DofIDItem id)
{
    for ( int l = 1; l <= giveNumberOfDofs(); l++ ) {
        if ( dofArray [ l - 1 ]->giveDofID() == id ) {
            return true;
        }
    }
    return false;
}


void DofManager :: giveLocationArray(const IntArray &dofIDArry, IntArray &locationArray, const UnknownNumberingScheme &s) const
{
    if ( !hasSlaveDofs ) {
        int size, indx;
        // prevents some size problem when connecting different elements with
        // different number of dofs
        size = dofIDArry.giveSize();
        locationArray.resize(size);
        for ( int i = 1; i <= size; i++ ) {
            if ( ( indx = this->findDofWithDofId( ( DofIDItem ) dofIDArry.at(i) ) ) == 0 ) {
                _error("giveLocationArray: incompatible dof requested");
            }

            locationArray.at(i) = s.giveDofEquationNumber( this->giveDof(indx) );
        }

    } else {
        IntArray dofArray, mstrEqNmbrs;
        int masterDofs = giveNumberOfPrimaryMasterDofs(dofArray);

        this->giveDofArray(dofIDArry, dofArray);
        locationArray.resize( masterDofs );

        for ( int k = 1, i = 1; i <= dofArray.giveSize(); i++ ) {
            this->giveDof(dofArray.at(i))->giveEquationNumbers(mstrEqNmbrs, s);
            locationArray.copySubVector(mstrEqNmbrs, k);
            k += mstrEqNmbrs.giveSize();
        }
    }
}


void DofManager :: giveMasterDofIDArray(const IntArray &dofIDArry, IntArray &masterDofIDs) const
{
    if ( !hasSlaveDofs ) {
        masterDofIDs = dofIDArry;
    } else {
        IntArray dofArray, temp;
        int masterDofs = giveNumberOfPrimaryMasterDofs(dofArray);

        this->giveDofArray(dofIDArry, dofArray);
        masterDofIDs.resize( masterDofs );

        for ( int k = 1, i = 1; i <= numberOfDofs; i++ ) {
            this->giveDof(i)->giveDofIDs(temp);
            masterDofIDs.copySubVector(temp, k);
            k += temp.giveSize();
        }
    }
}


void DofManager :: giveCompleteLocationArray(IntArray &locationArray, const UnknownNumberingScheme &s) const
{
    if ( !hasSlaveDofs ) {
        // prevents some size problem when connecting different elements with
        // different number of dofs
        locationArray.resize(numberOfDofs);
        for ( int i = 1; i <= numberOfDofs; i++ ) {
            locationArray.at(i) = s.giveDofEquationNumber( this->giveDof(i) );
        }
    } else {
        IntArray temp;
        int nMasterDofs = 0;
        for ( int i = 1; i <= numberOfDofs; i++ ) {
            nMasterDofs += this->giveDof(i)->giveNumberOfPrimaryMasterDofs();
        }
        locationArray.resize(nMasterDofs);
        for ( int k = 1, i = 1; i <= numberOfDofs; i++ ) {
            this->giveDof(i)->giveEquationNumbers(temp, s);
            locationArray.copySubVector(temp, k);
            k += temp.giveSize();
        }
    }
}


void DofManager :: giveCompleteMasterDofIDArray(IntArray &dofIDArray) const
{
    if ( !hasSlaveDofs ) {
        dofIDArray.resize(numberOfDofs);
        for ( int i = 1; i <= numberOfDofs; i++ ) {
            dofIDArray.at(i) = this->giveDof(i)->giveDofID();
        }
    } else {
        IntArray temp;
        int nMasterDofs = 0;
        for ( int i = 1; i <= numberOfDofs; i++ ) {
            nMasterDofs += this->giveDof(i)->giveNumberOfPrimaryMasterDofs();
        }
        dofIDArray.resize(nMasterDofs);
        for ( int k = 1, i = 1; i <= numberOfDofs; i++ ) {
            this->giveDof(i)->giveDofIDs(temp);
            dofIDArray.copySubVector(temp, k);
            k += temp.giveSize();
        }
    }
}


void DofManager :: giveDofArray(const IntArray &dofIDArry, IntArray &answer) const
// Returns the dof index array of the receiver.
// The location array contains the indexes of particular requsted DOFs
// In dofIDArray are stored DofID's of requsted DOFs in receiver.
// The DofID's are determining the physical meaning of particular DOFs
{
    int size;
    // prevents some size problem when connecting different elements with
    // different number of dofs
    size = dofIDArry.giveSize();
    answer.resize(size);
    for ( int i = 1; i <= size; i++ ) {
        if ( ( answer.at(i) = this->findDofWithDofId( ( DofIDItem ) dofIDArry.at(i) ) ) == 0 ) {
            _error("giveDofArray : incompatible dof requested");
        }
    }
}


int DofManager :: findDofWithDofId(DofIDItem dofID) const
{
    // finds index of DOF in receivers node with dofID
    // if such DOF does not exists, returns zero value
    for ( int i = 1; i <= numberOfDofs; i++ ) {
        if ( this->giveDof(i)->giveDofID() == dofID ) {
            return i;
        }
    }

    // nothing found
    return 0;
}


int DofManager :: giveNumberOfDofs() const
// Returns the number of degrees of freedom of the receiver.
{
    return numberOfDofs;
}


void DofManager :: setNumberOfDofs(int _ndofs)
{
    if ( _ndofs != this->giveNumberOfDofs() ) {
        if ( dofArray ) {
            int i = numberOfDofs;
            if ( numberOfDofs ) {
                while ( i-- ) {
                    delete dofArray [ i ];
                }

                delete[] dofArray;
            }
        }

        dofArray = new Dof * [ _ndofs ];
        for ( int i = 0; i < _ndofs; i++ ) {
            dofArray [ i ] = NULL;
        }

        this->numberOfDofs = _ndofs;
    }
}


void DofManager :: askNewEquationNumbers(TimeStep *tStep)
{
    for ( int i = 1; i <= this->numberOfDofs; i++ ) {
        this->giveDof(i)->askNewEquationNumber(tStep);
    }
}


int DofManager :: giveNumberOfPrimaryMasterDofs(IntArray &dofArray) const
{
    if ( !hasSlaveDofs ) {
        return dofArray.giveSize();
    }

    int answer = 0;

    for ( int i = 1; i <= dofArray.giveSize(); i++ ) {
        answer += this->giveDof( dofArray.at(i) )->giveNumberOfPrimaryMasterDofs();
    }

    return answer;
}


IRResultType DofManager ::  resolveDofIDArray(InputRecord *ir, IntArray &dofIDArry)
{
    const char *__proc = "resolveDofIDArray";
    IRResultType result;

    numberOfDofs = -1;
    IR_GIVE_OPTIONAL_FIELD(ir, numberOfDofs, _IFT_DofManager_ndofs);

    // returns nonzero if succes
    if ( numberOfDofs == -1 ) {
        dofIDArry = domain->giveDefaultNodeDofIDArry();
        numberOfDofs = dofIDArry.giveSize();
    } else {
        // if ndofs is prescribed, read the physical meaning of particular dofs
        // for detailed values of DofMask array see cltypes.h file
        // for exaple 1 is for D_u (displacemet in u dir), 2 for D_v, 3 for D_w, ...
        IR_GIVE_FIELD(ir, dofIDArry, _IFT_DofManager_dofidmask);

        if ( dofIDArry.giveSize() != numberOfDofs ) {
            _error("resolveDofIDArray : DofIDMask size mismatch");
        }

        this->dofidmask = new IntArray(dofIDArry);
    }

    return IRRT_OK;
}

IRResultType
DofManager :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                 // Required by IR_GIVE_FIELD macro

    delete dofidmask; dofidmask = NULL;
    delete dofTypemap; dofTypemap = NULL;
    delete dofMastermap; dofMastermap = NULL;
    delete dofBCmap; dofBCmap = NULL;
    delete dofICmap; dofICmap = NULL;

    IntArray dofIDArry;
    IntArray bc, ic, masterMask, dofTypeMask;

    loadArray.resize(0);
    IR_GIVE_OPTIONAL_FIELD(ir, loadArray, _IFT_DofManager_load);

#if 1
    ///@todo This method is unnecessary, we should just check if user has supplied a dofidmask field or not and just drop "numberOfDofs")
    this->resolveDofIDArray(ir, dofIDArry);
#else
    if ( ir->hasField(_IFT_DofManager_dofidmask) ) {
        IR_GIVE_FIELD(ir, dofIDArry, _IFT_DofManager_dofidmask);
        this->dofidmask = new IntArray(dofIDArry);
    } else {
        dofIDArry = domain->giveDefaultNodeDofIDArry();
    }
    numberOfDofs = dofIDArry.giveSize();
#endif

    bc.resize(0);
    IR_GIVE_OPTIONAL_FIELD(ir, bc, _IFT_DofManager_bc);

    ic.resize(0);
    IR_GIVE_OPTIONAL_FIELD(ir, ic, _IFT_DofManager_ic);

    // reads master mask - in this array are numbers of master dofManagers
    // to which are connected dofs in receiver.
    // if master mask index is zero then dof is created as master (i.e., having own equation number)
    // othervise slave dof connected to master DofManager is created.
    // by default if masterMask is not specifyed, all dofs are created as masters.
    dofTypeMask.resize(0); // termitovo
    IR_GIVE_OPTIONAL_FIELD(ir, dofTypeMask, _IFT_DofManager_doftypemask);

    // read boundary flag
    if ( ir->hasField(_IFT_DofManager_boundaryflag) ) {
        isBoundaryFlag = true;
    }


#ifdef __PARALLEL_MODE
    partitions.resize(0);
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
#endif


    int dofIc = 0, dofBc = 0;

    bool hasIc = !( ic.giveSize() == 0 );
    bool hasBc = !( bc.giveSize() == 0 );
    bool hasTypeinfo = !( dofTypeMask.giveSize() == 0 );

    ///@todo This should eventually be removed, still here to preserve backwards compatibility:
    if ( ( hasIc || hasBc || hasTypeinfo ) && !this->dofidmask )  this->dofidmask = new IntArray(dofIDArry);

    // check sizes
    if ( hasBc ) {
        if ( bc.giveSize() != this->giveNumberOfDofs() ) {
            _error3( "initializeFrom: bc size mismatch. Size is %d and need %d", bc.giveSize(), this->giveNumberOfDofs() );
        }
        this->dofBCmap = new std::map< int, int >();
        for (int i = 1; i <= bc.giveSize(); ++i) {
            if ( bc.at(i) > 0 ) {
                (*this->dofBCmap)[dofIDArry.at(i)] = bc.at(i);
            }
        }
    }

    if ( hasIc ) {
        if ( ic.giveSize() != this->giveNumberOfDofs() ) {
            _error3( "initializeFrom: ic size mismatch. Size is %d and need %d", ic.giveSize(), this->giveNumberOfDofs() );
        }
        this->dofICmap = new std::map< int, int >();
        for (int i = 1; i <= ic.giveSize(); ++i) {
            if ( ic.at(i) > 0 ) {
                (*this->dofICmap)[dofIDArry.at(i)] = ic.at(i);
            }
        }
    }

    if ( hasTypeinfo ) {
        if ( dofTypeMask.giveSize() != this->giveNumberOfDofs() ) {
            _error3( "initializeFrom: dofTypeMask size mismatch. Size is %d and need %d", dofTypeMask.giveSize(), this->giveNumberOfDofs() );
        }
        this->dofTypemap = new std::map< int, int >();
        for (int i = 1; i <= dofTypeMask.giveSize(); ++i) {
            if ( dofTypeMask.at(i) != DT_master ) {
                (*this->dofTypemap)[dofIDArry.at(i)] = dofTypeMask.at(i);
            }
        }
        // For simple slave dofs:
        if ( dofTypeMask.contains(DT_simpleSlave) ) {
            IR_GIVE_FIELD(ir, masterMask, _IFT_DofManager_mastermask);
            if ( masterMask.giveSize() != numberOfDofs ) {
                _error("initializeFrom: mastermask size mismatch");
            }
            this->dofMastermap = new std::map< int, int >();
            for (int i = 1; i <= masterMask.giveSize(); ++i) {
                if ( masterMask.at(i) > 0 ) {
                    (*this->dofMastermap)[dofIDArry.at(i)] = masterMask.at(i);
                }
            }
        }
    }

    dofArray = new Dof * [ this->giveNumberOfDofs() ];
    for ( int j = 0; j < numberOfDofs; j++ ) {
        dofType dtype;
        if ( hasTypeinfo ) {
            dtype = ( dofType ) dofTypeMask.at(j + 1);
        } else {
            dtype = DT_master;
        }

        if ( this->isDofTypeCompatible(dtype) ) {
            if ( hasIc ) {
                dofIc = ic.at(j + 1);
            }
            if ( hasBc ) {
                dofBc = bc.at(j + 1);
            }

            if ( dtype == DT_master ) {
                dofArray [ j ] = new MasterDof( j + 1, this, dofBc, dofIc, ( DofIDItem ) dofIDArry.at(j + 1) );
            } else if ( dtype == DT_active ) {
                dofArray [ j ] = new ActiveDof( j + 1, this, dofBc, ( DofIDItem ) dofIDArry.at(j + 1) );
            } else if ( dtype == DT_simpleSlave ) { // Simple slave dof
                dofArray [ j ] = new SimpleSlaveDof( j + 1, this, masterMask.at(j + 1), ( DofIDItem ) dofIDArry.at(j + 1) );
            } else if ( dtype == DT_slave ) { // Slave dof
                dofArray [ j ] = new SlaveDof( j + 1, this, ( DofIDItem ) dofIDArry.at(j + 1) );
            } else {
                _error2( "initializeFrom: unknown dof type (%s)",  __dofTypeToString(dtype) );
            }
        } else {
            _error("initializeFrom: incompatible dof type");
        }
    }

    return IRRT_OK;
}


void DofManager :: printOutputAt(FILE *stream, TimeStep *stepN)
{
    EngngModel *emodel = this->giveDomain()->giveEngngModel();

    fprintf( stream, "%-8s%8d (%8d):\n", this->giveClassName(), this->giveLabel(), this->giveNumber() );
    for ( int i = 1; i <= numberOfDofs; i++ ) {
        emodel->printDofOutputAt(stream, this->giveDof(i), stepN);
    }
}


void DofManager :: printYourself()
// Prints the receiver on screen.
{
    printf("DofManager %d\n", number);
    for ( int i = 0; i < numberOfDofs; i++ ) {
        if ( dofArray [ i ] ) {
            dofArray [ i ]->printYourself();
        } else {
            printf("dof %d is nil \n", i + 1);
        }
    }

    loadArray.printYourself();
    printf("\n");
}


void DofManager :: updateYourself(TimeStep *tStep)
// Updates the receiver at end of step.
{
    for ( int i = 1; i <= numberOfDofs; i++ ) {
        this->giveDof(i)->updateYourself(tStep);
    }
}


contextIOResultType DofManager :: saveContext(DataStream *stream, ContextMode mode, void *obj)
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


    if ( !stream->write(& numberOfDofs, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    // store dof types
    for ( int i = 1; i <= numberOfDofs; i++ ) {
        _val = this->giveDof(i)->giveDofType();
        if ( !stream->write(& _val, 1) ) {
            THROW_CIOERR(CIO_IOERR);
        }
    }

    if ( mode & CM_Definition ) {
        if ( ( iores = loadArray.storeYourself(stream, mode) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }

        _val = ( int ) isBoundaryFlag;
        if ( !stream->write(& _val, 1) ) {
            THROW_CIOERR(CIO_IOERR);
        }

        _val = ( int ) hasSlaveDofs;
        if ( !stream->write(& _val, 1) ) {
            THROW_CIOERR(CIO_IOERR);
        }

#ifdef __PARALLEL_MODE
        if ( !stream->write(& globalNumber, 1) ) {
            THROW_CIOERR(CIO_IOERR);
        }

        _val = ( int ) parallel_mode;
        if ( !stream->write(& _val, 1) ) {
            THROW_CIOERR(CIO_IOERR);
        }

        if ( ( iores = partitions.storeYourself(stream, mode) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }

#endif
    }

    for ( int i = 1; i <= numberOfDofs; i++ ) {
        if ( ( iores = this->giveDof(i)->saveContext(stream, mode, obj) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }
    }

    return CIO_OK;
}


contextIOResultType DofManager :: restoreContext(DataStream *stream, ContextMode mode, void *obj)
//
// restores full node context (saves state variables, that completely describe
// current state)
//
{
    contextIOResultType iores;
    int _val;

    if ( ( iores = FEMComponent :: restoreContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    int _numberOfDofs;
    if ( !stream->read(& _numberOfDofs, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    IntArray dtypes(_numberOfDofs);
    // restore dof types
    for ( int i = 1; i <= _numberOfDofs; i++ ) {
        if ( !stream->read(& dtypes.at(i), 1) ) {
            THROW_CIOERR(CIO_IOERR);
        }
    }

    // create new dofs if necessary
    bool samedofs = ( numberOfDofs == _numberOfDofs );
    if ( samedofs ) {
        // if size match, check types
        for ( int i = 1; i <= _numberOfDofs; i++ ) {
            if ( this->giveDof(i)->giveDofType() != dtypes.at(i) ) {
                samedofs = false;
                break;
            }
        }
    }

    if ( !samedofs ) {
        // delete old dofs
        if ( numberOfDofs ) {
            int i = numberOfDofs;
            while ( i-- ) {
                delete dofArray [ i ];
            }

            delete[] dofArray;
        }

        // allocate new ones
        dofArray = new Dof * [ _numberOfDofs ];
        for ( int i = 0; i < _numberOfDofs; i++ ) {
            dofArray [ i ] = classFactory.createDof( ( dofType ) dtypes(i), i + 1, this );
        }

        numberOfDofs = _numberOfDofs;
    }

    if ( mode & CM_Definition ) {
        if ( ( iores = loadArray.restoreYourself(stream, mode) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }

        if ( !stream->read(& _val, 1) ) {
            THROW_CIOERR(CIO_IOERR);
        }

        isBoundaryFlag = ( bool ) _val;
        if ( !stream->read(& _val, 1) ) {
            THROW_CIOERR(CIO_IOERR);
        }

        hasSlaveDofs = ( bool ) _val;
#ifdef __PARALLEL_MODE
        if ( !stream->read(& globalNumber, 1) ) {
            THROW_CIOERR(CIO_IOERR);
        }

        if ( !stream->read(& _val, 1) ) {
            THROW_CIOERR(CIO_IOERR);
        }

        parallel_mode = ( dofManagerParallelMode ) _val;
        if ( ( iores = partitions.restoreYourself(stream, mode) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }

#endif
    }

    for ( int i = 1; i <= numberOfDofs; i++ ) {
        if ( ( iores = this->giveDof(i)->restoreContext(stream, mode, obj) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }
    }

    return CIO_OK;
}


void DofManager :: giveUnknownVector(FloatArray &answer, const IntArray &dofIDArry, ValueModeType mode, TimeStep *stepN)
{
    int size;
    IntArray dofArray;

    answer.resize( size = dofIDArry.giveSize() );

    for ( int i = 1; i <= size; i++ ) {
        int pos = this->findDofWithDofId( ( DofIDItem ) dofIDArry.at(i) );
#ifdef DEBUG
        if ( pos == 0 ) {
            OOFEM_ERROR2("DofManager :: giveUnknownVector - Couldn't find dof with Dof ID %d", dofIDArry.at(i));
        }
#endif
        answer.at(i) = this->giveDof(pos)->giveUnknown(mode, stepN);
    }

    // Transform to global c.s.
    FloatMatrix L2G;
    if (this->computeL2GTransformation(L2G, dofIDArry)) {
        answer.rotatedWith(L2G, 'n');
    }
}


void DofManager :: giveUnknownVector(FloatArray &answer, const IntArray &dofIDArry,
                                PrimaryField &field, ValueModeType mode, TimeStep *stepN)
{
    int size;
    IntArray dofArray;

    answer.resize( size = dofIDArry.giveSize() );

    for ( int i = 1; i <= size; i++ ) {
        int pos = this->findDofWithDofId( ( DofIDItem ) dofIDArry.at(i) );
#ifdef DEBUG
        if ( pos == 0 ) {
            OOFEM_ERROR2("DofManager :: giveUnknownVector - Couldn't find dof with Dof ID %d", dofIDArry.at(i));
        }
#endif
        answer.at(i) = this->giveDof(pos)->giveUnknown(field, mode, stepN);
    }

    // Transform to global c.s.
    FloatMatrix L2G;
    if (this->computeL2GTransformation(L2G, dofIDArry)) {
        answer.rotatedWith(L2G, 'n');
    }
}


void DofManager :: giveCompleteUnknownVector(FloatArray &answer, ValueModeType mode, TimeStep *stepN)
{
    answer.resize( this->numberOfDofs );
    for ( int i = 1; i <= this->numberOfDofs; i++ ) {
        answer.at(i) = this->giveDof(i)->giveUnknown(mode, stepN);
    }
}


void DofManager :: givePrescribedUnknownVector(FloatArray &answer, const IntArray &dofIDArry,
                                          ValueModeType mode, TimeStep *stepN)
{
    int size;
    IntArray dofArray;

    answer.resize( size = dofIDArry.giveSize() );
    this->giveDofArray(dofIDArry, dofArray);

    for ( int j = 1; j <= size; j++ ) {
        answer.at(j) = this->giveDof( dofArray.at(j) )->giveBcValue(mode, stepN);
    }

    // Transform to global c.s.
    FloatMatrix L2G;
    if (this->computeL2GTransformation(L2G, dofIDArry)) {
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
    for ( int i = 1; i <= this->numberOfDofs; i++ ) {
        Dof *d = this->giveDof(i);
        double val = d->giveUnknown(mode, tStep);
        if (ut == DisplacementVector || ut == EigenVector) { // Just treat eigenvectors as displacement vectors (they are redundant)
            if      (d->giveDofID() == D_u) { dofIDArry.at(k) = D_u; localVector.at(k) = val; k++; }
            else if (d->giveDofID() == D_v) { dofIDArry.at(k) = D_v; localVector.at(k) = val; k++; }
            else if (d->giveDofID() == D_w) { dofIDArry.at(k) = D_w; localVector.at(k) = val; k++; }
        } else if (ut == VelocityVector) {
            if      (d->giveDofID() == V_u) { dofIDArry.at(k) = V_u; localVector.at(k) = val; k++; }
            else if (d->giveDofID() == V_v) { dofIDArry.at(k) = V_v; localVector.at(k) = val; k++; }
            else if (d->giveDofID() == V_w) { dofIDArry.at(k) = V_w; localVector.at(k) = val; k++; }
        } else {
            OOFEM_ERROR2("DofManager :: giveUnknownVectorOfType - Can't produce vector for unknown type: %d", ut);
        }
    }

    FloatMatrix L2G;
    if (this->computeL2GTransformation(L2G, dofIDArry)) {
        // Transform to global c.s.
        answer.beProductOf(L2G, localVector);
    } else {
        // No local c.s, just copy the values to respective index;
        answer.resize(3);
        answer.zero();
        for ( int i = 1; i <= k; i++ ) {
            if ( dofIDArry.at(i) == D_u || dofIDArry.at(i) == V_u )
                answer.at(1) = localVector.at(i);
            else if ( dofIDArry.at(i) == D_v || dofIDArry.at(i) == V_v )
                answer.at(2) = localVector.at(i);
            else if ( dofIDArry.at(i) == D_w || dofIDArry.at(i) == V_w )
                answer.at(3) = localVector.at(i);
        }
    }
}


int DofManager :: hasAnySlaveDofs()
{
    for ( int i = 1; i <= numberOfDofs; i++ ) {
        if ( !this->giveDof(i)->isPrimaryDof() ) {
            return 1;
        }
    }

    return 0;
}


bool DofManager :: giveMasterDofMans(IntArray &masters)
{
    IntArray _dof_masters;
    bool answer = false;

    masters.resize(0);
    for ( int i = 1; i <= numberOfDofs; i++ ) {
        if ( !this->giveDof(i)->isPrimaryDof() ) {
            answer = true;
            this->giveDof(i)->giveMasterDofManArray(_dof_masters);
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
    for ( int i = 1; i <= numberOfDofs; i++ ) {
        if ( !this->giveDof(i)->isPrimaryDof() ) {
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

    if (!hasL2G && !hasM2L) {
        answer.beEmptyMtrx();
        return false;
    } else if (hasL2G && hasM2L) {
        answer.beProductOf(L2G, M2L);
    } else if (hasL2G) {
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
    if (!this->hasAnySlaveDofs()) {
        return false;
    }

    IntArray dofArray;
    FloatArray mstrContrs;

    if ( dofIDArry.isEmpty() ) {
        dofArray.resize(numberOfDofs);
        for ( int i = 1; i <= numberOfDofs; i++ ) {
            dofArray.at(i) = i;
        }
    } else {
        this->giveDofArray( dofIDArry, dofArray);
    }

    answer.resize( dofArray.giveSize(), giveNumberOfPrimaryMasterDofs(dofArray) );
    answer.zero();

    int indx;
    for ( int k = 1, i = 1; i <= dofArray.giveSize(); i++ ) {
        indx = dofArray.at(i);
        this->giveDof(indx)->computeDofTransformation(mstrContrs);
        answer.copySubVectorRow(mstrContrs, i, k);
        k += mstrContrs.giveSize();
    }
    return true;
}


bool DofManager :: requiresTransformation()
{
    return this->hasAnySlaveDofs();
}


void DofManager :: updateLocalNumbering(EntityRenumberingFunctor &f)
{
    for ( int i = 1; i <= numberOfDofs; i++ ) {
        this->giveDof(i)->updateLocalNumbering(f);
    }
}


#ifdef __PARALLEL_MODE
void DofManager :: mergePartitionList(IntArray &_p)
{
    // more optimized version can be made requiring sorted partition list of receiver
    int size = _p.giveSize();
    for ( int i = 1; i <= size; i++ ) {
        partitions.insertOnce( _p.at(i) );
    }
}


int
DofManager :: packDOFsUnknowns(CommunicationBuffer &buff,
                               ValueModeType mode, TimeStep *stepN)
{
    int result = 1;
    for ( int i = 1; i <= numberOfDofs; i++ ) {
        result &= this->giveDof(i)->packUnknowns(buff, mode, stepN);
    }

    return result;
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

#endif
} // end namespace oofem
