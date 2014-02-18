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

#include "structuralms.h"
#include "structuralmaterial.h"
#include "contextioerr.h"
#include "nlstructuralelement.h"
#include "gausspoint.h"

namespace oofem {

StructuralMaterialStatus :: StructuralMaterialStatus(int n, Domain *d, GaussPoint *g) :
    MaterialStatus(n, d, g), strainVector(), stressVector(),
    tempStressVector(), tempStrainVector(), FVector(), tempFVector()
{
    int rsize = StructuralMaterial :: giveSizeOfVoigtSymVector( gp->giveMaterialMode() );
    strainVector.resize(rsize);
    stressVector.resize(rsize);

    // reset temp vars.
    tempStressVector = stressVector;
    tempStrainVector = strainVector;

    if ( NLStructuralElement * el = dynamic_cast< NLStructuralElement * >( gp->giveElement() ) ) {
        if ( el->giveGeometryMode() == 1  ) { // if large def, initialize F and P
            PVector.resize(9);
            PVector.zero();
            FVector.resize(9);
            FVector.zero();
            FVector.at(1) = FVector.at(2) = FVector.at(3) = 1.;
        }
        tempPVector = PVector;
        tempFVector = FVector;
    }
}


StructuralMaterialStatus :: ~StructuralMaterialStatus() { }


void StructuralMaterialStatus :: printOutputAt(FILE *File, TimeStep *tNow)
// Prints the strains and stresses on the data file.
{
    FloatArray helpVec;
    int n;

    MaterialStatus :: printOutputAt(File, tNow);

    fprintf(File, "  strains ");
    StructuralMaterial :: giveFullSymVectorForm( helpVec, strainVector, gp->giveMaterialMode() );
    n = helpVec.giveSize();
    for ( int i = 1; i <= n; i++ ) {
        fprintf( File, " % .4e", helpVec.at(i) );
    }

    fprintf(File, "\n              stresses");
    StructuralMaterial :: giveFullSymVectorForm( helpVec, stressVector, gp->giveMaterialMode() );

    n = helpVec.giveSize();
    for ( int i = 1; i <= n; i++ ) {
        fprintf( File, " % .4e", helpVec.at(i) );
    }
    fprintf(File, "\n");
}


void StructuralMaterialStatus :: updateYourself(TimeStep *tStep)
// Performs end-of-step updates.
{
    MaterialStatus :: updateYourself(tStep);

    stressVector = tempStressVector;
    strainVector = tempStrainVector;
    PVector      = tempPVector;
    FVector      = tempFVector;
}


void StructuralMaterialStatus :: initTempStatus()
//
// initialize record at the begining of new load step
//
{
    MaterialStatus :: initTempStatus();

    // see if vectors describing reached equilibrium are defined
    if ( this->giveStrainVector().giveSize() == 0 ) {
        strainVector.resize( StructuralMaterial :: giveSizeOfVoigtSymVector( gp->giveMaterialMode() ) );
    }

    if ( this->giveStressVector().giveSize() == 0 ) {
        stressVector.resize( StructuralMaterial :: giveSizeOfVoigtSymVector( gp->giveMaterialMode() ) );
    }

    // reset temp vars.
    tempStressVector = stressVector;
    tempStrainVector = strainVector;
    tempPVector      = PVector;
    tempFVector      = FVector;
}


contextIOResultType
StructuralMaterialStatus :: saveContext(DataStream *stream, ContextMode mode, void *obj)
//
// saves full ms context (saves state variables, that completely describe
// current state)
{
    contextIOResultType iores;
    if ( stream == NULL ) {
        _error("saveContex : can't write into NULL stream");
    }

    if ( ( iores = MaterialStatus :: saveContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = strainVector.storeYourself(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = stressVector.storeYourself(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    return CIO_OK;
}


contextIOResultType
StructuralMaterialStatus :: restoreContext(DataStream *stream, ContextMode mode, void *obj)
//
// restores full material context (saves state variables, that completely describe
// current state)
//
{
    contextIOResultType iores;
    if ( stream == NULL ) {
        _error("saveContex : can't write into NULL stream");
    }

    if ( ( iores = MaterialStatus :: restoreContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = strainVector.restoreYourself(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = stressVector.restoreYourself(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    return CIO_OK;
}

void StructuralMaterialStatus :: copyStateVariables(const MaterialStatus &iStatus)
{
    MaterialStatus &tmpStat = const_cast< MaterialStatus & >(iStatus);
    const StructuralMaterialStatus &structStatus = dynamic_cast< StructuralMaterialStatus & >(tmpStat);

    strainVector = structStatus.giveStrainVector();
    stressVector = structStatus.giveStressVector();
    tempStressVector = structStatus.giveTempStressVector();
    tempStrainVector = structStatus.giveTempStrainVector();

    PVector = structStatus.givePVector();
    tempPVector = structStatus.giveTempPVector();
    CVector = structStatus.giveCVector();
    tempCVector = structStatus.giveTempCVector();
    FVector = structStatus.giveFVector();
    tempFVector = structStatus.giveTempFVector();
}

void StructuralMaterialStatus :: addStateVariables(const MaterialStatus &iStatus)
{
    printf("Entering StructuralMaterialStatus :: addStateVariables().\n");
}
} // end namespace oofem
