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


#include "structuralinterfacematerial.h"
#include "contextioerr.h"

#include "structuralmaterial.h"
#include "structuralinterfacematerialstatus.h"
#include "gausspoint.h"
namespace oofem {

StructuralInterfaceMaterialStatus :: StructuralInterfaceMaterialStatus(int n, Domain *d, GaussPoint *g) :
    MaterialStatus(n, d, g), jump(), traction(), firstPKTraction(),
    tempTraction(), tempJump(), F(), tempF(), tempFirstPKTraction()
{
    int size = StructuralMaterial :: giveSizeOfVoigtSymVector( gp->giveMaterialMode() ); ///@todo how to best get rid of matMode?

    this->jump.resize(size);
    this->traction.resize(size);
    this->firstPKTraction.resize(size);
    this->F.resize(size,size);

    // reset temp vars.
    this->tempJump            = this->jump;
    this->tempTraction        = this->traction;
    this->tempFirstPKTraction = this->firstPKTraction;
    this->tempF               = this->F;

}


StructuralInterfaceMaterialStatus :: ~StructuralInterfaceMaterialStatus() { }


void StructuralInterfaceMaterialStatus :: printOutputAt(FILE *File, TimeStep *tNow)
// Prints the strains and stresses on the data file.
{
    //FloatArray helpVec;
    //int n;

    //MaterialStatus :: printOutputAt(File, tNow);

    //fprintf(File, "  strains ");
    //StructuralMaterial :: giveFullSymVectorForm(helpVec, strainVector, gp->giveMaterialMode());
    //n = helpVec.giveSize();
    //for ( int i = 1; i <= n; i++ ) {
    //    fprintf( File, " % .4e", helpVec.at(i) );
    //}

    //fprintf(File, "\n              stresses");
    //StructuralMaterial :: giveFullSymVectorForm(helpVec, stressVector, gp->giveMaterialMode());

    //n = helpVec.giveSize();
    //for ( int i = 1; i <= n; i++ ) {
    //    fprintf( File, " % .4e", helpVec.at(i) );
    //}
    //fprintf(File, "\n");
}


void StructuralInterfaceMaterialStatus :: updateYourself(TimeStep *tStep)
// Performs end-of-step updates.
{
    MaterialStatus :: updateYourself(tStep);

    this->jump            = this->tempJump;
    this->traction        = this->tempTraction; 
    this->firstPKTraction = this->tempFirstPKTraction; 
    this->F               = this->tempF;
}


void StructuralInterfaceMaterialStatus :: initTempStatus()
//
// initialize record at the begining of new load step
//
{
    MaterialStatus :: initTempStatus();

    // see if vectors describing reached equilibrium are defined
    if ( this->giveJump().giveSize() == 0 ) {
        this->jump.resize( StructuralMaterial :: giveSizeOfVoigtSymVector( gp->giveMaterialMode() ) );
    }

    if ( this->giveTraction().giveSize() == 0 ) {
        this->traction.resize( StructuralMaterial :: giveSizeOfVoigtSymVector( gp->giveMaterialMode() ) );
    }

    // reset temp vars.
    this->tempJump            = this->jump;
    this->tempTraction        = this->traction;
    this->tempFirstPKTraction = this->firstPKTraction;
    this->tempF               = this->F;
}


contextIOResultType
StructuralInterfaceMaterialStatus :: saveContext(DataStream *stream, ContextMode mode, void *obj)
//
// saves full ms context (saves state variables, that completely describe
// current state)
{
    //contextIOResultType iores;
    //if ( stream == NULL ) {
    //    _error("saveContex : can't write into NULL stream");
    //}

    //if ( ( iores = MaterialStatus :: saveContext(stream, mode, obj) ) != CIO_OK ) {
    //    THROW_CIOERR(iores);
    //}

    //if ( ( iores = strainVector.storeYourself(stream, mode) ) != CIO_OK ) {
    //    THROW_CIOERR(iores);
    //}

    //if ( ( iores = stressVector.storeYourself(stream, mode) ) != CIO_OK ) {
    //    THROW_CIOERR(iores);
    //}

    return CIO_OK;
}


contextIOResultType
StructuralInterfaceMaterialStatus :: restoreContext(DataStream *stream, ContextMode mode, void *obj)
//
// restores full material context (saves state variables, that completely describe
// current state)
//
{
    //contextIOResultType iores;
    //if ( stream == NULL ) {
    //    _error("saveContex : can't write into NULL stream");
    //}

    //if ( ( iores = MaterialStatus :: restoreContext(stream, mode, obj) ) != CIO_OK ) {
    //    THROW_CIOERR(iores);
    //}

    //if ( ( iores = strainVector.restoreYourself(stream, mode) ) != CIO_OK ) {
    //    THROW_CIOERR(iores);
    //}

    //if ( ( iores = stressVector.restoreYourself(stream, mode) ) != CIO_OK ) {
    //    THROW_CIOERR(iores);
    //}

    return CIO_OK;
}
} // end namespace oofem
