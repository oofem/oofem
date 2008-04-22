/* $Header: /home/cvs/bp/oofem/tm/src/isoheatmat.h,v 1.1 2003/04/14 16:01:39 bp Exp $ */
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
 *               Copyright (C) 1993 - 2008   Borek Patzak
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


//   *****************************************************************
//   *** CLASS Heat and Moisture transport MATERIAL WITH HYDRATION ***
//   *****************************************************************

#ifndef hydratinghemomat_h
#define hydratinghemomat_h

#include "hemotkmat.h"
#include "hydratingisoheatmat.h"
#include "../sm/hydram.h"
#include "dictionr.h"
#include "flotarry.h"
#include "flotmtrx.h"
#include "cltypes.h"

class GaussPoint;

class HydratingHeMoMaterial : public HeMoTKMaterial, public HydrationModelInterface
{
    /*
     * This class implements a isotropic linear heat  material in a finite element problem.
     * A material is an attribute of a domain. It is usually also attribute of many elements.
     *
     * DESCRIPTION
     * Isotropic Linear Heat Material with interface to the Hydration Model
     *
     * TASK
     *
     */

protected:
    int hydration, hydrationHeat, hydrationLHS, teplotaOut;

public:

    HydratingHeMoMaterial(int n, Domain *d) : HeMoTKMaterial(n, d), HydrationModelInterface() { }
    ~HydratingHeMoMaterial() { }

    void setMixture(MixtureType mix);

    virtual int hasInternalSource(); // return true if hydration heat source is present
    virtual void computeInternalSourceVector(FloatArray &val, GaussPoint *gp, TimeStep *atTime, ValueModeType mode);
    /**
     * Updates internal state of material according to new state vector.
     * @param vec new state vector
     * @param gp integration point
     * @param tStep solution step
     */
    virtual void updateInternalState(const FloatArray &vec, GaussPoint *gp, TimeStep *);

    /*
     * void  giveCharacteristicMatrix (FloatMatrix& answer,
     *                               MatResponseForm form,
     *                               MatResponseMode mode,
     *                               GaussPoint* gp,
     *                               TimeStep* atTime);
     */

    virtual double  giveCharacteristicValue(MatResponseMode mode,
                                            GaussPoint *gp,
                                            TimeStep *atTime);

    // saves current context(state) into stream
    contextIOResultType saveContext(DataStream *stream, ContextMode mode, void *obj = NULL);
    contextIOResultType restoreContext(DataStream *stream, ContextMode mode, void *obj = NULL);

    // identification and auxiliary functions
    const char *giveClassName() const { return "HydratingHeMoMaterial"; }
    classType giveClassID()         const { return HydratingHeMoMaterialClass; }

    IRResultType initializeFrom(InputRecord *ir);

    /*
     * // non-standard - returns time independent material constant
     * double  give (int) ;
     */

    // post-processing
    virtual int giveIPValue(FloatArray &answer, GaussPoint *aGaussPoint, InternalStateType type, TimeStep *atTime);
    virtual InternalStateValueType giveIPValueType(InternalStateType type);
    virtual int giveIntVarCompFullIndx(IntArray &answer, InternalStateType type, MaterialMode mmode);
    virtual int giveIPValueSize(InternalStateType type, GaussPoint *aGaussPoint);

protected:
    virtual MaterialStatus *CreateStatus(GaussPoint *gp) const;
};


#endif // hydratinghemomat_h
