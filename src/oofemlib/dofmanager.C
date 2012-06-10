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
 *               Copyright (C) 1993 - 2012   Borek Patzak
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
#include "flotarry.h"
#include "flotmtrx.h"
#include "intarray.h"
#include "usrdefsub.h"
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
    int i, n, nLoads;
    Load *loadN;
    FloatArray contribution;

    if ( this->giveLoadArray()->isEmpty() ) {
        answer.resize(0);
        return;
    } else {
        answer.resize(0);
        nLoads = loadArray.giveSize();     // the node may be subjected
        for ( i = 1; i <= nLoads; i++ ) {   // to more than one load
            n            = loadArray.at(i);
            loadN        = domain->giveLoad(n);
            if ( loadN->giveBCGeoType() != NodalLoadBGT ) {
                _error("computeLoadVectorAt: incompatible load type applied");
            }

            loadN->computeComponentArrayAt(contribution, stepN, mode); // can be NULL
            answer.add(contribution);
        }
    }
}


Dof *DofManager :: giveDof(int i) const
// Returns the i-th degree of freedom of the receiver.
{
    if ( !dofArray ) {
        _error("giveDof: dof is not defined");
    }

    return dofArray [ i - 1 ];
}


Dof *DofManager :: giveDofWithID(int dofID) const
// Returns the degree of freedom of the receiver with 'dofID'.
{
    int indx = this->findDofWithDofId( ( DofIDItem ) dofID );

    if ( !indx ) {
        _error("giveDofWithID: dof with given DofID doesnot exists");
    }

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
// Returns the location array of the receiver. Creates this array if it
// does not exist yet. The location array contains the equation number of
// every requested degree of freedom of the receiver.
// In dofIDArray are stored DofID's of requsted DOFs in receiver.
// The DofID's are determining the physical meaning of particular DOFs
{
    if ( !hasSlaveDofs ) {
        int i, size, indx;
        // prevents some size problem when connecting different elements with
        // different number of dofs
        size = dofIDArry.giveSize();
        locationArray.resize(size);
        for ( i = 1; i <= size; i++ ) {
            if ( ( indx = this->findDofWithDofId( ( DofIDItem ) dofIDArry.at(i) ) ) == 0 ) {
                _error("giveLocationArray: incompatible dof requested");
            }

            locationArray.at(i) = s.giveDofEquationNumber( this->giveDof(indx) );
        }
    } else {
        int i, k, indx;
        IntArray dofArray, mstrEqNmbrs;

        this->giveDofArray(dofIDArry, dofArray);
        locationArray.resize( giveNumberOfPrimaryMasterDofs(dofArray) );

        for ( k = 1, i = 1; i <= dofArray.giveSize(); i++ ) {
            indx = dofArray.at(i);
            this->giveDof(indx)->giveEquationNumbers(mstrEqNmbrs, s);
            locationArray.copySubVector(mstrEqNmbrs, k);
            k += mstrEqNmbrs.giveSize();
        }
    }
}


void DofManager :: giveCompleteLocationArray(IntArray &locationArray, const UnknownNumberingScheme &s) const
// Returns the complete location array of the receiver.
// including all available dofs
{
    if ( !hasSlaveDofs ) {
        int i;
        // prevents some size problem when connecting different elements with
        // different number of dofs
        locationArray.resize(numberOfDofs);
        for ( i = 1; i <= numberOfDofs; i++ ) {
            locationArray.at(i) = s.giveDofEquationNumber( this->giveDof(i) );
        }
    } else {
        giveLocationArray(* giveCompleteGlobalDofIDArray(), locationArray, s);
    }
}


void DofManager :: giveDofArray(const IntArray &dofIDArry, IntArray &answer) const
// Returns the dof index array of the receiver.
// The location array contains the indexes of particular requsted DOFs
// In dofIDArray are stored DofID's of requsted DOFs in receiver.
// The DofID's are determining the physical meaning of particular DOFs
{
    int i, size;
    // prevents some size problem when connecting different elements with
    // different number of dofs
    size = dofIDArry.giveSize();
    answer.resize(size);
    for ( i = 1; i <= size; i++ ) {
        if ( ( answer.at(i) = this->findDofWithDofId( ( DofIDItem ) dofIDArry.at(i) ) ) == 0 ) {
            _error("giveDofArray : incompatible dof requested");
        }
    }
}


int DofManager :: findDofWithDofId(DofIDItem dofID) const
{
    // finds index of DOF in receivers node with dofID
    // if such DOF does not exists, returns zero value
    int i;
    for ( i = 1; i <= numberOfDofs; i++ ) {
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
    int i;
    if ( _ndofs != this->giveNumberOfDofs() ) {
        if ( dofArray ) {
            i = numberOfDofs;
            if ( numberOfDofs ) {
                while ( i-- ) {
                    delete dofArray [ i ];
                }

                delete[] dofArray;
            }
        }

        dofArray = new Dof * [ _ndofs ];
        for ( i = 0; i < _ndofs; i++ ) {
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

    int i, answer = 0;

    for ( i = 1; i <= dofArray.giveSize(); i++ ) {
        answer += this->giveDof( dofArray.at(i) )->giveNumberOfPrimaryMasterDofs();
    }

    return answer;
}


IRResultType DofManager ::  resolveDofIDArray(InputRecord *ir, IntArray &dofIDArry)
{
    const char *__keyword, *__proc = "resolveDofIDArray";
    IRResultType result;

    numberOfDofs = -1;
    __keyword = "ndofs";
    result = ir->giveOptionalField(numberOfDofs, IFT_DofManager_ndofs, __keyword);
    if ( result != IRRT_OK ) {
        IR_IOERR(giveClassName(), __proc, IFT_DofManager_ndofs, __keyword, ir, result);
    }

    // returns nonzero if succes
    if ( numberOfDofs == -1 ) {
        dofIDArry = domain->giveDefaultNodeDofIDArry();
        numberOfDofs = dofIDArry.giveSize();
    } else {
        // if ndofs is prescribed, read the physical meaning of particular dofs
        // for detailed values of DofMask array see cltypes.h file
        // for exaple 1 is for D_u (displacemet in u dir), 2 for D_v, 3 for D_w, ...
        __keyword = "dofidmask";
        result = ir->giveField(dofIDArry, IFT_DofManager_dofidmask, __keyword);
        if ( result != IRRT_OK ) {
            IR_IOERR(giveClassName(), __proc, IFT_DofManager_dofidmask, __keyword, ir, result);
        }

        if ( dofIDArry.giveSize() != numberOfDofs ) {
            _error("resolveDofIDArray : DofIDMask size mismatch");
        }
    }

    return IRRT_OK;
}

IRResultType
DofManager :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                 // Required by IR_GIVE_FIELD macro

    int j;
    IntArray dofIDArry;
    IntArray bc, ic, masterMask, dofTypeMask; // termitovo

    loadArray.resize(0);
    IR_GIVE_OPTIONAL_FIELD(ir, loadArray, IFT_DofManager_load, "load"); // Macro

    if ( this->resolveDofIDArray(ir, dofIDArry) != IRRT_OK ) {
        IR_IOERR(giveClassName(), __proc,  IFT_Unknown, "", ir, result);
    }

    // numberOfDofs = domain->giveNumberOfDofs () ;
    bc.resize(0);
    IR_GIVE_OPTIONAL_FIELD(ir, bc, IFT_DofManager_bc, "bc"); // Macro

    ic.resize(0);
    IR_GIVE_OPTIONAL_FIELD(ir, ic, IFT_DofManager_ic, "ic"); // Macro
    // reads master mask - in this array are numbers of master dofManagers
    // to which are connected dofs in receiver.
    // if master mask index is zero then dof is created as master (i.e., having own equation number)
    // othervise slave dof connected to master DofManager is created.
    // by default if masterMask is not specifyed, all dofs are created as masters.
    dofTypeMask.resize(0); // termitovo
    IR_GIVE_OPTIONAL_FIELD(ir, dofTypeMask, IFT_DofManager_doftypemask, "doftype"); // Macro

    // read boundary flag
    if ( ir->hasField(IFT_DofManager_boundaryflag, "boundary") ) {
        isBoundaryFlag = true;
    }


#ifdef __PARALLEL_MODE
 #ifndef __ENABLE_COMPONENT_LABELS
    globalNumber = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, globalNumber, IFT_DofManager_globnum, "globnum"); // Macro
 #endif

    partitions.resize(0);
    IR_GIVE_OPTIONAL_FIELD(ir, partitions, IFT_DofManager_partitions, "partitions"); // Macro

    if ( ir->hasField(IFT_DofManager_sharedflag, "shared") ) {
        parallel_mode = DofManager_shared;
    } else if ( ir->hasField(IFT_DofManager_remoteflag, "remote") ) {
        parallel_mode = DofManager_remote;
    } else if ( ir->hasField(IFT_DofManager_nullflag, "null") ) {
        parallel_mode = DofManager_null;
    } else {
        parallel_mode = DofManager_local;
    }

    // in parallel mode,  slaves are allowed, because ((Dr. Rypl promissed)
    // masters have to be in same partition as slaves. They can be again Remote copies.
#endif



    int hasIc, hasBc, dofIc = 0, dofBc = 0, hasTypeinfo = 0;
    dofType dtype;

    hasIc = !( ic.giveSize() == 0 );
    hasBc = !( bc.giveSize() == 0 );
    hasTypeinfo = !( dofTypeMask.giveSize() == 0 );

    // check sizes
    if ( hasBc ) {
        if ( bc.giveSize() != this->giveNumberOfDofs() ) {
            _error3( "initializeFrom: bc size mismatch. Size is %d and need %d", bc.giveSize(), this->giveNumberOfDofs() );
        }
    }

    if ( hasIc ) {
        if ( ic.giveSize() != this->giveNumberOfDofs() ) {
            _error3( "initializeFrom: ic size mismatch. Size is %d and need %d", ic.giveSize(), this->giveNumberOfDofs() );
        }
    }

    if ( hasTypeinfo ) {
        if ( dofTypeMask.giveSize() != this->giveNumberOfDofs() ) {
            _error3( "initializeFrom: dofTypeMask size mismatch. Size is %d and need %d", dofTypeMask.giveSize(), this->giveNumberOfDofs() );
        }
    }

    dofArray = new Dof * [ this->giveNumberOfDofs() ];
    for ( j = 0; j < numberOfDofs; j++ ) {
        if ( hasTypeinfo ) {
            dtype = ( dofType ) dofTypeMask.at(j + 1);
        } else {
            dtype = DT_master;
        }

        if ( this->isDofTypeCompatible(dtype) ) {
            if ( dtype == DT_master ) {
                if ( hasIc ) {
                    dofIc = ic.at(j + 1);
                }

                if ( hasBc ) {
                    dofBc = bc.at(j + 1);
                }

                dofArray [ j ] = new MasterDof( j + 1, this, dofBc, dofIc, ( DofIDItem ) dofIDArry.at(j + 1) );
            } else if ( dtype == DT_active ) {
                if ( hasBc ) {
                    dofBc = bc.at(j + 1);
                }
                dofArray [ j ] = new ActiveDof( j + 1, this, dofBc, ( DofIDItem ) dofIDArry.at(j + 1) );
            } else if ( dtype == DT_simpleSlave ) { // Simple slave dof
                if ( masterMask.giveSize() == 0 ) {
                    IR_GIVE_FIELD(ir, masterMask, IFT_DofManager_mastermask, "mastermask"); // Macro
                    if ( masterMask.giveSize() != numberOfDofs ) {
                        _error("initializeFrom: mastermask size mismatch");
                    }
                }

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
    int i;

#if defined( __PARALLEL_MODE ) || defined( __ENABLE_COMPONENT_LABELS )
    fprintf( stream, "%-8s%8d (%8d):\n", this->giveClassName(), this->giveLabel(), this->giveNumber() );
#else
    fprintf( stream, "%-8s%8d:\n", this->giveClassName(), this->giveNumber() );
#endif
    for ( i = 1; i <= numberOfDofs; i++ ) {
        emodel->printDofOutputAt(stream, this->giveDof(i), stepN);
    }
}


void DofManager :: printYourself()
// Prints the receiver on screen.
{
    int i;

    printf("DofManager %d\n", number);
    for ( i = 0; i < numberOfDofs; i++ ) {
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
    int i;

    for ( i = 1; i <= numberOfDofs; i++ ) {
        this->giveDof(i)->updateYourself(tStep);
    }
}


contextIOResultType DofManager :: saveContext(DataStream *stream, ContextMode mode, void *obj)
//
// saves full node context (saves state variables, that completely describe
// current state)
//
{
    int i, _val;
    contextIOResultType iores;

    if ( ( iores = FEMComponent :: saveContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }


    if ( !stream->write(& numberOfDofs, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    // store dof types
    for ( i = 1; i <= numberOfDofs; i++ ) {
        _val =  this->giveDof(i)->giveClassID();
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

    for ( i = 1; i <= numberOfDofs; i++ ) {
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
    int i, _val;

    if ( ( iores = FEMComponent :: restoreContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    int _numberOfDofs;
    if ( !stream->read(& _numberOfDofs, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    IntArray dtypes(_numberOfDofs);
    // restore dof types
    for ( i = 1; i <= _numberOfDofs; i++ ) {
        if ( !stream->read(& dtypes.at(i), 1) ) {
            THROW_CIOERR(CIO_IOERR);
        }
    }

    // create new dofs if necessary
    bool samedofs = ( numberOfDofs == _numberOfDofs );
    if ( samedofs ) {
        // if size match, check types
        for ( i = 1; i <= _numberOfDofs; i++ ) {
            if ( this->giveDof(i)->giveClassID() != dtypes.at(i) ) {
                samedofs = false;
                break;
            }
        }
    }

    if ( !samedofs ) {
        // delete old dofs
        if ( numberOfDofs ) {
            i = numberOfDofs;
            while ( i-- ) {
                delete dofArray [ i ];
            }

            delete[] dofArray;
        }

        // allocate new ones
        dofArray = new Dof * [ _numberOfDofs ];
        for ( i = 0; i < _numberOfDofs; i++ ) {
            dofArray [ i ] = CreateUsrDefDofOfType( ( classType ) dtypes(i), i + 1, this );
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

    for ( i = 1; i <= numberOfDofs; i++ ) {
        if ( ( iores = this->giveDof(i)->restoreContext(stream, mode, obj) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }
    }

    return CIO_OK;
}


void DofManager :: giveUnknownVector(FloatArray &answer, const IntArray &dofIDArry,
                                EquationID type, ValueModeType mode, TimeStep *stepN)
{
    int j, size;
    IntArray dofArray;

    answer.resize( size = dofIDArry.giveSize() );
    this->giveDofArray(dofIDArry, dofArray);

    for ( j = 1; j <= size; j++ ) {
        answer.at(j) = this->giveDof( dofArray.at(j) )->giveUnknown(type, mode, stepN);
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
    int j, size;
    IntArray dofArray;

    answer.resize( size = dofIDArry.giveSize() );
    this->giveDofArray(dofIDArry, dofArray);

    for ( j = 1; j <= size; j++ ) {
        answer.at(j) = this->giveDof( dofArray.at(j) )->giveUnknown(field, mode, stepN);
    }

    // Transform to global c.s.
    FloatMatrix L2G;
    if (this->computeL2GTransformation(L2G, dofIDArry)) {
        answer.rotatedWith(L2G, 'n');
    }
}


void DofManager :: givePrescribedUnknownVector(FloatArray &answer, const IntArray &dofIDArry,
                                          ValueModeType mode, TimeStep *stepN)
{
    int j, size;
    IntArray dofArray;

    answer.resize( size = dofIDArry.giveSize() );
    this->giveDofArray(dofIDArry, dofArray);

    for ( j = 1; j <= size; j++ ) {
        answer.at(j) = this->giveDof( dofArray.at(j) )->giveBcValue(mode, stepN);
    }

    // Transform to global c.s.
    FloatMatrix L2G;
    if (this->computeL2GTransformation(L2G, dofIDArry)) {
        answer.rotatedWith(L2G, 'n');
    }
}


int DofManager :: hasAnySlaveDofs()
{
    int i;
    for ( i = 1; i <= numberOfDofs; i++ ) {
        if ( !this->giveDof(i)->isPrimaryDof() ) {
            return 1;
        }
    }

    return 0;
}


bool DofManager :: giveMasterDofMans(IntArray &masters)
{
    int i, j;
    IntArray _dof_masters;
    bool answer = false;

    masters.resize(0);
    for ( i = 1; i <= numberOfDofs; i++ ) {
        if ( !this->giveDof(i)->isPrimaryDof() ) {
            answer = true;
            this->giveDof(i)->giveMasterDofManArray(_dof_masters);
            for ( j = 1; j <= _dof_masters.giveSize(); j++ ) {
                masters.insertSortedOnce(_dof_masters.at(j), 2);
            }
        }
    }

    return answer;
}


int DofManager :: checkConsistency()
// Checks internal data consistency in node.
// Current implementation checks (when receiver has simple slave dofs) if receiver
// has the same coordinate system as master dofManager of slave dof.
{
    hasSlaveDofs = false;
    for (int i = 1; i <= numberOfDofs; i++ ) {
        if ( !this->giveDof(i)->isPrimaryDof() ) {
            hasSlaveDofs = true;
            continue;
        }
    }

    return 1;
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
    for (int k = 1, i = 1; i <= dofArray.giveSize(); i++ ) {
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


IntArray *DofManager :: giveCompleteGlobalDofIDArray() const
{
    IntArray *answer = new IntArray(numberOfDofs);

    for ( int i = 1; i <= numberOfDofs; i++ ) {
        answer->at(i) = ( int ) this->giveDof(i)->giveDofID();
    }

    return answer;
}


void DofManager :: updateLocalNumbering(EntityRenumberingFunctor &f)
{
    int i;
    for ( i = 1; i <= numberOfDofs; i++ ) {
        this->giveDof(i)->updateLocalNumbering(f);
    }
}


#ifdef __PARALLEL_MODE
void DofManager :: mergePartitionList(IntArray &_p)
{
    // more optimized version can be made requiring sorted partition list of receiver
    int i, size = _p.giveSize();
    for ( i = 1; i <= size; i++ ) {
        partitions.insertOnce( _p.at(i) );
    }
}


int
DofManager :: packDOFsUnknowns(CommunicationBuffer &buff, EquationID type,
                               ValueModeType mode, TimeStep *stepN)
{
    int i, result = 1;
    for ( i = 1; i <= numberOfDofs; i++ ) {
        result &= this->giveDof(i)->packUnknowns(buff, type, mode, stepN);
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
