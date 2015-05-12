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

#include "AbaqusUserElement.h"
#include "gausspoint.h"
#include "classfactory.h"
#include "dynamicinputrecord.h"
#include "dofmanager.h"
#include "node.h"

#ifdef _WIN32 //_MSC_VER and __MINGW32__ included
 #include <Windows.h>
#else
 #include <dlfcn.h>
#endif

#include <cstring>

namespace oofem {
REGISTER_Element(AbaqusUserElement);


AbaqusUserElement :: AbaqusUserElement(int n, Domain *d) :
    NLStructuralElement(n, d), uelobj(NULL), hasTangentFlag(false), uel(NULL)
{}

AbaqusUserElement :: ~AbaqusUserElement()
{
#ifdef _WIN32
    if ( this->uelobj ) {
        FreeLibrary( ( HMODULE ) this->uelobj );
    }
#else
    if ( this->uelobj ) {
        dlclose(this->uelobj);
    }
#endif
}


IRResultType AbaqusUserElement :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                                        // Required by IR_GIVE_FIELD macro

    result = StructuralElement :: initializeFrom(ir);
    if ( result != IRRT_OK ) {
        return result;
    }

    this->numberOfDofMans = dofManArray.giveSize();

    // necessary to prevent an array dimension error in Init
    IR_GIVE_FIELD(ir, nCoords, _IFT_AbaqusUserElement_numcoords);

    IR_GIVE_OPTIONAL_FIELD(ir, this->dofs, _IFT_AbaqusUserElement_dofs);
    IR_GIVE_FIELD(ir, this->numSvars, _IFT_AbaqusUserElement_numsvars);
    if ( this->numSvars < 0 ) {
        OOFEM_ERROR("'numsvars' field has an invalid value");
    }
    IR_GIVE_FIELD(ir, this->props, _IFT_AbaqusUserElement_properties);
    IR_GIVE_FIELD(ir, this->jtype, _IFT_AbaqusUserElement_type);
    if ( this->jtype < 0 ) {
        OOFEM_ERROR("'type' has an invalid value");
    }
    IR_GIVE_FIELD(ir, this->filename, _IFT_AbaqusUserElement_userElement);

#if 0
    uelname = "uel";
    IR_GIVE_OPTIONAL_FIELD(ir, uelname, _IFT_AbaqusUserElement_name);
#endif

#ifdef _WIN32
    this->uelobj = ( void * ) LoadLibrary( filename.c_str() );
    if ( !this->uelobj ) {
        OOFEM_ERROR( "couldn't load \"%s\",\ndlerror: %s", filename.c_str() );
    }

    * ( FARPROC * ) ( & this->uel ) = GetProcAddress( ( HMODULE ) this->uelobj, "uel_" );  //works for MinGW 32bit
    if ( !this->uel ) {
        // char *dlresult = GetLastError();
        DWORD dlresult = GetLastError();                 //works for MinGW 32bit
        OOFEM_ERROR("couldn't load symbol uel,\nerror: %s\n", dlresult);
    }
#else
    this->uelobj = dlopen(filename.c_str(), RTLD_NOW);
    if ( !this->uelobj ) {
        OOFEM_ERROR( "couldn't load \"%s\",\ndlerror: %s", filename.c_str(), dlerror() );
    }

    * ( void ** ) ( & this->uel ) = dlsym(this->uelobj, "uel_");
    char *dlresult = dlerror();
    if ( dlresult ) {
        OOFEM_ERROR("couldn't load symbol uel,\ndlerror: %s\n", dlresult);
    }
#endif

    return IRRT_OK;
}


void AbaqusUserElement :: postInitialize()
{
    NLStructuralElement :: postInitialize();

    this->ndofel = this->numberOfDofMans * this->nCoords;
    this->mlvarx = this->ndofel;
    this->nrhs = 2;
    this->rhs.resize(this->ndofel, this->nrhs);
    this->amatrx.resize(this->ndofel, this->ndofel);
    this->svars.resize(this->numSvars);
    this->lFlags.resize(5);
    this->predef.resize( this->npredef * this->numberOfDofMans * 2 );
    this->energy.resize(8);
    this->U.resize(this->ndofel);
    this->V.resize(this->ndofel);
    this->A.resize(this->ndofel);
    this->DU.resize(this->ndofel, this->nrhs);

    if ( !this->coords.isNotEmpty() ) {
        this->coords.resize(this->numberOfDofMans, this->mcrd);
        for ( int j = 1; j <= numberOfDofMans; j++ ) {
            Node *dm = this->giveNode(j);
            for ( int i = 1; i <= mcrd; i++ ) {
                this->coords.at(i, j) = dm->giveCoordinate(i);
            }
        }
    }
}


void AbaqusUserElement :: giveInputRecord(DynamicInputRecord &input)
{
    StructuralElement :: giveInputRecord(input);

    input.setField(this->coords, _IFT_AbaqusUserElement_numcoords);
    input.setField(this->dofs, _IFT_AbaqusUserElement_dofs);
    input.setField(this->numSvars, _IFT_AbaqusUserElement_numsvars);
    input.setField(this->props, _IFT_AbaqusUserElement_properties);
    input.setField(this->jtype, _IFT_AbaqusUserElement_type);
    input.setField(this->filename, _IFT_AbaqusUserElement_userElement);
}

Interface *AbaqusUserElement :: giveInterface(InterfaceType it)
{
    return NULL;
}

void AbaqusUserElement :: giveDofManDofIDMask(int inode, IntArray &answer) const
{
    answer = this->dofs;
}

void AbaqusUserElement :: computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
{
    if ( !hasTangent() ) {
        // use uel to calculate the tangent
        FloatArray forces;
        giveInternalForcesVector(forces, tStep, U, DU, 0);
    }
    // give tangent
    answer = giveTempTangent();
    // add stuff to behave differently if mUseNumericalTangent is set?
}

void AbaqusUserElement :: updateYourself(TimeStep *tStep)
{
    StructuralElement :: updateYourself(tStep);
    svars = tempSvars;
    amatrx = tempAmatrx;
    rhs = tempRHS;
    hasTangentFlag = false;
}

void AbaqusUserElement :: updateInternalState(TimeStep *tStep)
{
    FloatArray tmp;
    this->giveInternalForcesVector(tmp, tStep, 0);
}

void AbaqusUserElement :: giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord)
{
    // init U vector
    this->computeVectorOf(VM_Total, tStep, U);
    FloatArray tempIntVect;
    // init DU vector
    computeVectorOf(VM_Incremental, tStep, tempIntVect);
    this->giveDomain()->giveClassName();
    DU.zero();
    DU.setColumn(tempIntVect, 1);
    //this->computeVectorOf(VM_Total, tStep, DU);
    this->giveInternalForcesVector(answer, tStep, U, DU, useUpdatedGpRecord);
}

void AbaqusUserElement :: giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, 
                                            FloatArray &U, FloatMatrix &DU, int useUpdatedGpRecord)
{
    if ( useUpdatedGpRecord ) {
        this->rhs.copyColumn(answer, 1);
    } else {
        this->lFlags.at(1) = 1;                 // 1 based access
        this->lFlags.at(3) = 1;                 // 1 based access
        this->lFlags.at(4) = 0;                 // 1 based access

        int nprops = props.giveSize();
        int njprops = jprops.giveSize();

        FloatMatrix loc_rhs(this->ndofel, this->nrhs);
        FloatMatrix loc_amatrx(this->ndofel, this->ndofel);
        FloatArray loc_svars = this->giveStateVector();

        //this->getSvars();
        double period = 0., pnewdt = 0.;
        double dtime = tStep->giveTimeIncrement();
        double time[] = {tStep->giveTargetTime() - dtime, tStep->giveTargetTime()};
        this->uel(
            loc_rhs.givePointer(),
            loc_amatrx.givePointer(),
            loc_svars.givePointer(),
            energy.givePointer(),
            & ndofel,
            & nrhs,
            & numSvars,
            props.givePointer(),
            & nprops,
            coords.givePointer(),
            & mcrd,
            & this->numberOfDofMans,
            U.givePointer(),
            DU.givePointer(),
            V.givePointer(),
            A.givePointer(),
            & jtype,
            time,
            & dtime,
            & kstep,
            & kinc,
            & ( this->number ),
            params,
            & ndLoad,
            jdltype,
            adlmag.givePointer(),
            predef.givePointer(),
            & npredef,
            lFlags.givePointer(),
            & mlvarx,
            ddlmag.givePointer(),
            & mdLoad,
            & pnewdt,
            jprops.givePointer(),
            & njprops,
            & period);
        //FloatArray vout;
        //vout.resize(12);
        //for (int i = 1; i <= 3; i++)
        //{
        //	vout.at(i) = rhs.at(i, 1);
        //	vout.at(i+6) = rhs.at(i+3, 1);
        //}
        //answer = vout;
        //this->rhs.copyColumn(answer, 1);
        //answer.negated();
        loc_rhs.negated();                      //really needed???
        loc_rhs.copyColumn(answer, 1);
        letTempRhsBe(loc_rhs);
        letTempTangentBe(loc_amatrx);
        letTempSvarsBe(loc_svars);
    }
}


void
AbaqusUserElement :: computeConsistentMassMatrix(FloatMatrix &answer, TimeStep *tStep, double &mass, const double *ipDensity)
{
    answer.resize(ndofel, ndofel);
    answer.zero();
}
}       // namespace oofem
