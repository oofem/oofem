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

#ifndef lsmastermatgrad_h
#define lsmastermatgrad_h

#include "lsmastermat.h"
#include "structuralmaterial.h"
#include "structuralms.h"
#include "linearelasticmaterial.h"
#include "dictionary.h"
#include "floatarray.h"
#include "floatmatrix.h"

namespace oofem {
class GaussPoint;
class Domain;

/**
 * This class implements an isotropic elastoplastic material
 * with Mises yield condition, associated flow rule
 * and linear isotropic hardening.
 *
 * It differs from other similar materials (such as J2Mat)
 * by implementation - here we use the radial return, which
 * is the most efficient algorithm for this specific model.
 * Also, an extension to large strain will be available.
 *
 */
class LsMasterMatGrad : public LsMasterMat
{
public:
    LsMasterMatGrad(int n, Domain *d);
    virtual ~LsMasterMatGrad();

    virtual int hasMaterialModeCapability(MaterialMode mode);
    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual void giveCharacteristicMatrix(FloatMatrix &answer, MatResponseForm form, MatResponseMode rMode, GaussPoint *gp, TimeStep *atTime);
    void give3dKappaMatrix(FloatMatrix &answer, MatResponseForm form, MatResponseMode mode, GaussPoint *gp, TimeStep *atTime);
    void give3dGprime(FloatMatrix &answer, MatResponseForm form, MatResponseMode mode, GaussPoint *gp, TimeStep *atTime);
    void giveInternalLength(FloatMatrix &answer, MatResponseForm form, MatResponseMode mode, GaussPoint *gp, TimeStep *atTime);
    virtual const char *giveClassName() const { return "LsMasterMatGrad"; }
    virtual classType giveClassID() const { return LsMasterMatClass; }

    MaterialStatus *CreateStatus(GaussPoint *gp) const;
    virtual void give3dMaterialStiffnessMatrix(FloatMatrix & answer,
                                       MatResponseForm, MatResponseMode,
                                       GaussPoint * gp,
                                       TimeStep * atTime);

    virtual void giveRealStressVector(FloatArray & answer,  MatResponseForm, GaussPoint *,
                              const FloatArray &, TimeStep *);
};

//=============================================================================


class LsMasterMatGradStatus : public LsMasterMatStatus
{
public:
    LsMasterMatGradStatus(int n, Domain *d, GaussPoint *g, int s);
    virtual ~LsMasterMatGradStatus();

    virtual  void printOutputAt(FILE *file, TimeStep *tStep);

    virtual void initTempStatus();

    virtual void updateYourself(TimeStep *);

    virtual contextIOResultType saveContext(DataStream *stream, ContextMode mode, void *obj = NULL);
    virtual contextIOResultType restoreContext(DataStream *stream, ContextMode mode, void *obj = NULL);
    virtual const char *giveClassName() const { return "LsMasterMatGradStatus"; }

    classType giveClassID() const { return LsMasterMatStatusClass; }
};
} // end namespace oofem
#endif // misesmat_h
