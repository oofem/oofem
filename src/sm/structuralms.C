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

#include "structuralms.h"
#include "structuralcrosssection.h"
#include "structuralmaterial.h"
#include "contextioerr.h"
#include "nlstructuralelement.h"
namespace oofem {
StructuralMaterialStatus :: StructuralMaterialStatus(int n, Domain *d, GaussPoint *g) :
    MaterialStatus(n, d, g), strainVector(), stressVector(), FVector(),
    tempStressVector(), tempStrainVector(), tempFVector()
{
    StructuralMaterial *mat = static_cast< StructuralMaterial * >( gp->giveMaterial() );
    int rsize = mat->giveSizeOfReducedStressStrainVector( gp->giveMaterialMode() );
    strainVector.resize(rsize);
    stressVector.resize(rsize);

    // reset temp vars.
    tempStressVector = stressVector;
    tempStrainVector = strainVector;

    if ( NLStructuralElement *el = dynamic_cast< NLStructuralElement * > ( gp->giveElement() ) ) {
        if ( el->giveGeometryMode() == 1  ) { // if large def, initialize F and P
            mat->giveIdentityVector(FVector, gp->giveMaterialMode());
            PVector.resize(FVector.giveSize());                        
            tempPVector = PVector;
            tempFVector = FVector;
        }
    }

}


StructuralMaterialStatus :: ~StructuralMaterialStatus() { }


void StructuralMaterialStatus :: printOutputAt(FILE *File, TimeStep *tNow)
// Prints the strains and stresses on the data file.
{
    FloatArray helpVec;
    StructuralCrossSection *cs = static_cast< StructuralCrossSection * >( gp->giveCrossSection() );
    int n;

    MaterialStatus :: printOutputAt(File, tNow);

    fprintf(File, "  strains ");
    cs->giveFullCharacteristicVector(helpVec, gp, strainVector);
    n = helpVec.giveSize();
    for ( int i = 1; i <= n; i++ ) {
        fprintf( File, " % .4e", helpVec.at(i) );
    }

    fprintf(File, "\n              stresses");
    cs->giveFullCharacteristicVector(helpVec, gp, stressVector);

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
        strainVector.resize( static_cast< StructuralMaterial * >( gp->giveMaterial() )->
                            giveSizeOfReducedStressStrainVector( gp->giveMaterialMode() ) );
    }

    if ( this->giveStressVector().giveSize() == 0 ) {
        stressVector.resize( static_cast< StructuralMaterial * >( gp->giveMaterial() )->
                            giveSizeOfReducedStressStrainVector( gp->giveMaterialMode() ) );
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
} // end namespace oofem
