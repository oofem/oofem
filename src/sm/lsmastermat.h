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

#ifndef lsmastermat_h
#define lsmastermat_h

#include "structuralmaterial.h"
#include "structuralms.h"
#include "linearelasticmaterial.h"
#include "dictionary.h"
#include "floatarray.h"
#include "floatmatrix.h"

///@name Input fields for LsMasterMat
//@{
#define _IFT_LsMasterMat_m "m"
#define _IFT_LsMasterMat_slaveMat "slavemat"
//@}

namespace oofem {
class GaussPoint;
class Domain;

/**
 * Large strain master material.
 * Stress and stiffness are computed from small strain(slaveMat) material model 
 * using a strain tensor from the Seth-Hill strain tensors family (depends on parameter m,
 * m = 0 logarithmic strain,m = 1 Green-Lagrange strain ...)
 * then stress and stiffness are transformed 2.PK stress and appropriate stiffness 
 */
class LsMasterMat : public StructuralMaterial
{
protected:
    /// Reference to the basic elastic material.
    LinearElasticMaterial *linearElasticMaterial;

    /// 'slave' material model number.
    int slaveMat;
    /// Specifies the strain tensor.
    double m;

 
public:
    LsMasterMat(int n, Domain *d);
    virtual ~LsMasterMat();

    virtual int hasMaterialModeCapability(MaterialMode mode);
    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual int hasNonLinearBehaviour() { return 1; }
    virtual const char *giveClassName() const { return "LsMasterMat"; }
    virtual classType giveClassID() const { return LsMasterMatClass; }

    LinearElasticMaterial *giveLinearElasticMaterial() { return linearElasticMaterial; }

    virtual bool isCharacteristicMtrxSymmetric(MatResponseMode rMode) { return false; }

    virtual MaterialStatus *CreateStatus(GaussPoint *gp) const;

    virtual void give3dMaterialStiffnessMatrix(FloatMatrix & answer,
                                       MatResponseForm, MatResponseMode,
                                       GaussPoint * gp,
                                       TimeStep * atTime);

    virtual void giveRealStressVector(FloatArray & answer,  MatResponseForm, GaussPoint *,
                              const FloatArray &, TimeStep *);

    /// transformation matrix
    void constructTransformationMatrix(FloatMatrix F, GaussPoint *gp);

    
protected:
    virtual int giveIPValue(FloatArray &answer, GaussPoint *aGaussPoint, InternalStateType type, TimeStep *atTime);
    virtual int giveIntVarCompFullIndx(IntArray &answer, InternalStateType type, MaterialMode mmode);
    virtual InternalStateValueType giveIPValueType(InternalStateType type);
    virtual int giveIPValueSize(InternalStateType type, GaussPoint *aGaussPoint);
};

//=============================================================================


class LsMasterMatStatus : public StructuralMaterialStatus
{
protected:
    FloatMatrix Pmatrix,TLmatrix, transformationMatrix;
    int slaveMat;

public:
    LsMasterMatStatus(int n, Domain *d, GaussPoint *g, int s);
    virtual ~LsMasterMatStatus();


    void givePmatrix(FloatMatrix &answer)
    { answer = Pmatrix; }
    void giveTLmatrix(FloatMatrix &answer)
    { answer = TLmatrix; }
    void giveTransformationMatrix(FloatMatrix &answer)
    { answer = transformationMatrix;}

    void setPmatrix(FloatMatrix values) { Pmatrix = values; }
    void setTLmatrix(FloatMatrix values) { TLmatrix = values; }
    void setTransformationMatrix(FloatMatrix values){transformationMatrix = values;}

    virtual void printOutputAt(FILE *file, TimeStep *tStep);

    virtual void initTempStatus();

    virtual void updateYourself(TimeStep *);

    virtual contextIOResultType saveContext(DataStream *stream, ContextMode mode, void *obj = NULL);

    virtual contextIOResultType restoreContext(DataStream *stream, ContextMode mode, void *obj = NULL);

    virtual const char *giveClassName() const { return "LsMasterMatStatus"; }

    virtual classType giveClassID() const { return LsMasterMatStatusClass; }
};
} // end namespace oofem
#endif // misesmat_h
