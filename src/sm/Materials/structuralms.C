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

#include "sm/Materials/structuralms.h"
#include "sm/Materials/structuralmaterial.h"
#include "contextioerr.h"
#include "sm/Elements/nlstructuralelement.h"
#include "gausspoint.h"

namespace oofem {
StructuralMaterialStatus :: StructuralMaterialStatus(GaussPoint *g) :
    MaterialStatus(g), strainVector(), stressVector(),
    tempStressVector(), tempStrainVector(), FVector(), tempFVector()
{
    int rsize = StructuralMaterial :: giveSizeOfVoigtSymVector( gp->giveMaterialMode() );
    strainVector.resize(rsize);
    stressVector.resize(rsize);

    // reset temp vars.
    tempStressVector = stressVector;
    tempStrainVector = strainVector;

    /// Hack to prevent crashing when there are only "loose" gausspoints
    if ( gp->giveIntegrationRule() == NULL ) {
        return;
    }
    if ( NLStructuralElement * el = dynamic_cast< NLStructuralElement * >( gp->giveElement() ) ) {
      if ( el->giveGeometryMode() == 1  ) { // if large def, initialize F and P
	PVector.resize(9);
	FVector.resize(9);
	FVector.at(1) = FVector.at(2) = FVector.at(3) = 1.;
	tempPVector = PVector;
	tempFVector = FVector;
	
      }
    }

}


void StructuralMaterialStatus :: printOutputAt(FILE *File, TimeStep *tStep) const
// Prints the strains and stresses on the data file.
{
    FloatArray helpVec;

    MaterialStatus :: printOutputAt(File, tStep);
    NLStructuralElement * el = static_cast< NLStructuralElement * >( gp->giveElement());
    if ( el->giveGeometryMode() == 1) {
      fprintf(File, "  F ");
      StructuralMaterial :: giveFullVectorFormF( helpVec, FVector, gp->giveMaterialMode() );
      for ( auto &var : helpVec ) {
        fprintf( File, " %.4e", var );
      }

      fprintf(File, "\n  P");
      StructuralMaterial :: giveFullVectorForm( helpVec, PVector, gp->giveMaterialMode() );
      for ( auto &var : helpVec ) {
        fprintf( File, " %.4e", var );
      }
    } else {
      fprintf(File, "  strains ");
      StructuralMaterial :: giveFullSymVectorForm( helpVec, strainVector, gp->giveMaterialMode() );
      for ( auto &var : helpVec ) {
        fprintf( File, " %.4e", var );
      }
      
      fprintf(File, "\n              stresses");
      StructuralMaterial :: giveFullSymVectorForm( helpVec, stressVector, gp->giveMaterialMode() );
      
      for ( auto &var : helpVec ) {
        fprintf( File, " %.4e", var );
      }
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


void
StructuralMaterialStatus :: saveContext(DataStream &stream, ContextMode mode)
{
    MaterialStatus :: saveContext(stream, mode);

    contextIOResultType iores;
    if ( ( iores = strainVector.storeYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = stressVector.storeYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }
}


void
StructuralMaterialStatus :: restoreContext(DataStream &stream, ContextMode mode)
{
    MaterialStatus :: restoreContext(stream, mode);

    contextIOResultType iores;
    if ( ( iores = strainVector.restoreYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = stressVector.restoreYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }
}

void StructuralMaterialStatus :: copyStateVariables(const MaterialStatus &iStatus)
{
    const StructuralMaterialStatus &structStatus = static_cast< const StructuralMaterialStatus & >(iStatus);

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
