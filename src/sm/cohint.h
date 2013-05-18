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

#ifndef cohint_h
#define cohint_h

#include "material.h"
#include "linearelasticmaterial.h"
#include "structuralmaterial.h"
#include "structuralms.h"
#include "cltypes.h"

namespace oofem {

///@name Input fields for RankineMat
//@{
#define _IFT_CohesiveInterfaceMaterial_Name "cohint"
#define _IFT_IsoInterfaceDamageMaterial_kn "kn"
#define _IFT_IsoInterfaceDamageMaterial_ks "ks"
//@}

/**
 * This class implements associated Material Status to CohesiveInterfaceMaterial.
 * It is attribute of matStatusDictionary at every GaussPoint, for which this material
 * is active.
 *
 * @author Milan Jirasek
 */
class CohesiveInterfaceMaterialStatus : public StructuralMaterialStatus
{
public:
    /// Constructor
    CohesiveInterfaceMaterialStatus(int n, Domain *d, GaussPoint *g);
    /// Destructor
    virtual ~CohesiveInterfaceMaterialStatus();

    virtual void printOutputAt(FILE *file, TimeStep *tStep);

    // definition
    virtual const char *giveClassName() const { return "CohesiveInterfaceMaterialStatus"; }
    virtual classType giveClassID() const { return CohesiveInterfaceMaterialStatusClass; }

    virtual void initTempStatus();
    virtual void updateYourself(TimeStep *tStep);

    virtual contextIOResultType saveContext(DataStream *stream, ContextMode mode, void *obj = NULL);
    virtual contextIOResultType restoreContext(DataStream *stream, ContextMode mode, void *obj = NULL);
};


/**
 * Class representing cohesive interface material model.
 */
class CohesiveInterfaceMaterial : public StructuralMaterial
{
protected:
    /// Coefficient of thermal dilatation.
    double tempDillatCoeff;
    /// Elastic properties (normal and shear moduli).
    double kn, ks;

public:
    /// Constructor
    CohesiveInterfaceMaterial(int n, Domain *d);
    /// Destructor
    virtual ~CohesiveInterfaceMaterial() {};

    virtual int hasNonLinearBehaviour() { return 0; }
    virtual int hasMaterialModeCapability(MaterialMode mode) { return ( mode == _3dInterface ); }

    virtual const char *giveClassName() const { return "CohesiveInterfaceMaterial"; }
    virtual classType giveClassID() const { return CohesiveInterfaceMaterialClass; }

    virtual void give3dMaterialStiffnessMatrix(FloatMatrix & answer,
                                               MatResponseForm, MatResponseMode,
                                               GaussPoint *gp,
                                               TimeStep *atTime);

    double computeVolumetricStrain(GaussPoint *gp, TimeStep *atTime);

    virtual void giveRealStressVector(FloatArray & answer,  MatResponseForm, GaussPoint *gp,
                              const FloatArray &reducedStrain, TimeStep *tStep);

    virtual void  giveCharacteristicMatrix(FloatMatrix &answer,
                                           MatResponseForm form,
                                           MatResponseMode mode,
                                           GaussPoint *gp,
                                           TimeStep *atTime);

    virtual int giveStressStrainComponentIndOf(MatResponseForm, MaterialMode mmode, int);
    virtual void giveStressStrainMask(IntArray & answer, MatResponseForm, MaterialMode mmode) const;
    virtual int giveSizeOfReducedStressStrainVector(MaterialMode);
    void giveReducedCharacteristicVector(FloatArray &answer, GaussPoint *gp,
                                         const FloatArray &charVector3d);
    void giveFullCharacteristicVector(FloatArray &answer,  GaussPoint *gp,
                                      const FloatArray &strainVector);

    virtual int giveIPValue(FloatArray &answer, GaussPoint *aGaussPoint, InternalStateType type, TimeStep *atTime);
    virtual int giveIntVarCompFullIndx(IntArray &answer, InternalStateType type, MaterialMode mmode);

    virtual InternalStateValueType giveIPValueType(InternalStateType type);
    virtual int giveIPValueSize(InternalStateType type, GaussPoint *aGaussPoint);
    virtual void giveThermalDilatationVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep);

    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual int giveInputRecordString(std :: string &str, bool keyword = true);
    virtual void giveInputRecord(DynamicInputRecord &input);

    virtual MaterialStatus *CreateStatus(GaussPoint *gp) const { return new CohesiveInterfaceMaterialStatus(1, FEMComponent :: domain, gp); }

protected:
    // Overloaded to use specialized versions of these services possibly implemented by linearElastic member
    void give3dInterfaceMaterialStiffnessMatrix(FloatMatrix &answer, MatResponseForm form, MatResponseMode rMode,
                                                GaussPoint *gp, TimeStep *atTime);
};
} // namespace oofem
#endif

