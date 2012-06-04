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

#include "fluiddynamicmaterial.h"
#include "domain.h"
#include "verbose.h"
#include "gausspnt.h"
#include "flotmtrx.h"
#include "flotarry.h"

#include "mathfem.h"
#include "fieldmanager.h"

#ifndef __MAKEDEPEND
 #include <stdlib.h>
#endif
#include "contextioerr.h"

namespace oofem {
void
FluidDynamicMaterial :: updateInternalState(const FloatArray &vec, GaussPoint *gp, TimeStep *tStep)
{ }


void
FluidDynamicMaterial :: computeDeviatoricStressVector(FloatArray &stress_dev, double &epsp_vol, GaussPoint *gp, const FloatArray &eps, double pressure, TimeStep *tStep)
{
    if ( gp->giveMaterialMode() == _2dFlow ) {
        epsp_vol = -(eps.at(1) + eps.at(2));
    } else {
        epsp_vol = -(eps.at(1) + eps.at(2) + eps.at(3));
    }
    this->computeDeviatoricStressVector(stress_dev, gp, eps, tStep);
}


void
FluidDynamicMaterial :: giveDeviatoricPressureStiffness(FloatArray &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
    int size = ((FluidDynamicMaterialStatus*)gp->giveMaterialStatus())->giveDeviatoricStressVector().giveSize();
    answer.resize(size);
    answer.zero();
}


void
FluidDynamicMaterial :: giveVolumetricDeviatoricStiffness(FloatArray &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
    int size = ((FluidDynamicMaterialStatus*)gp->giveMaterialStatus())->giveDeviatoricStressVector().giveSize();
    answer.resize(size);
    answer.zero();
}


void
FluidDynamicMaterial :: giveVolumetricPressureStiffness(double &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
    answer = 0.0;
}


FluidDynamicMaterialStatus :: FluidDynamicMaterialStatus(int n, Domain *d, GaussPoint *g) :
    MaterialStatus(n, d, g), deviatoricStressVector()
{ }

void
FluidDynamicMaterialStatus :: printOutputAt(FILE *File, TimeStep *tNow)
// Prints the strains and stresses on the data file.
{
    int i, n;

    fprintf(File, "\n deviatoric stresses");
    n = deviatoricStressVector.giveSize();
    for ( i = 1; i <= n; i++ ) {
        fprintf( File, " % .4e", deviatoricStressVector.at(i) );
    }

    fprintf(File, "\n");
}

void
FluidDynamicMaterialStatus :: updateYourself(TimeStep *tStep)
// Performs end-of-step updates.
{
    MaterialStatus :: updateYourself(tStep);
}


void
FluidDynamicMaterialStatus :: initTempStatus()
//
// initialize record at the beginning of new load step
//
{
    MaterialStatus :: initTempStatus();
}


int
FluidDynamicMaterial :: giveIPValue(FloatArray &answer, GaussPoint *aGaussPoint, InternalStateType type, TimeStep *atTime)
{
  /*
    FluidDynamicMaterialStatus *status = ( FluidDynamicMaterialStatus * ) this->giveStatus(aGaussPoint);
    if ( type == IST_DeviatoricStrain ) {
        answer = status->giveDeviatoricStressVector();
        return 1;
    }
  */
    if (type == IST_Viscosity ) {
        answer.resize(1); answer.at(1) = this->giveCharacteristicValue (MRM_Viscosity, aGaussPoint, atTime);
        return 1;
    } else {
        return Material :: giveIPValue(answer, aGaussPoint, type, atTime);
    }
}


InternalStateValueType
FluidDynamicMaterial :: giveIPValueType(InternalStateType type)
{
  /*
  if ( ( type == IST_DeviatoricStrain ) ) {
    return ISVT_TENSOR_S3;
  } */
    if (type == IST_Viscosity ) {
        return ISVT_SCALAR;
    } else {
        return Material :: giveIPValueType(type);
    }
}


int
FluidDynamicMaterial :: giveIntVarCompFullIndx(IntArray &answer, InternalStateType type, MaterialMode mmode)
{
  /*
    if ( ( type == IST_DeviatoricStrain ) ) {
      if ( mmode == _2dFlow ) {
	answer.resize(6);
	answer.at(1) = 1;
	answer.at(2) = 2;
	answer.at(6) = 3;
	return 1;
      } else if (mmode == _2dAxiFlow ) {
	answer.resize(6);
	answer.at(1) = 1;
	answer.at(2) = 2;
	answer.at(3) = 3;
	answer.at(6) = 4;
	return 1;
      } else if (mmode == _3dFlow ) {
	answer.resize(6);
	int i;
	for (i=1; i<=6; i++) answer.at(i) = i;
	return 1;
      } else {
	OOFEM_ERROR ("FluidDynamicMaterial :: giveIntVarCompFullIndx: material mode not supported");
	return 0;
      }
    }
  */
    if (type == IST_Viscosity ) {
        answer.resize(1); answer.at(1) = 1;
        return 1;
    }  else {
        return Material :: giveIntVarCompFullIndx(answer, type, mmode);
    }
}


int
FluidDynamicMaterial :: giveIPValueSize(InternalStateType type, GaussPoint *aGaussPoint)
{
  /*
    MaterialMode mmode = aGaussPoint->giveMaterialMode();
    if ( ( type == IST_DeviatoricStrain ) ) {
      if ( mmode == _2dFlow ) {
	return 3;
      } else if (mmode == _2dAxiFlow ) {
	return 4;
      } else if (mmode == _3dFlow ) {
	return 6;
      } else {
	OOFEM_ERROR ("FluidDynamicMaterial :: giveIPValueSize: material mode not supported");
	return 0;
      }
   }
  */
    if (type == IST_Viscosity ) {
        return 1;
    } else {
        return Material :: giveIPValueSize(type, aGaussPoint);
    }
}


contextIOResultType
FluidDynamicMaterialStatus :: saveContext(DataStream *stream, ContextMode mode, void *obj)
//
// saves full ms context (saves state variables, that completely describe
// current state)
// saving the data in  TDictionary is left to material (yield crit. level).
{
    contextIOResultType iores;
    if ( stream == NULL ) {
        _error("saveContex : can't write into NULL stream");
    }

    if ( ( iores = MaterialStatus :: saveContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = deviatoricStressVector.storeYourself(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    return CIO_OK;
}


contextIOResultType
FluidDynamicMaterialStatus :: restoreContext(DataStream *stream, ContextMode mode, void *obj)
//
// restores full material context (saves state variables, that completely describe
// current state)
//
{
    // FloatArray *s;
    contextIOResultType iores;
    if ( stream == NULL ) {
        _error("saveContex : can't write into NULL stream");
    }

    if ( ( iores = MaterialStatus :: restoreContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = deviatoricStressVector.restoreYourself(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    return CIO_OK;
}
} // end namespace oofem
