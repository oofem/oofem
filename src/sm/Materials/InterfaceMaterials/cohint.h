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

#ifndef cohint_h
#define cohint_h

#include "material.h"
#include "Materials/linearelasticmaterial.h"
#include "../sm/Materials/structuralmaterial.h"
#include "../sm/Materials/structuralms.h"
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
 * @author Milan Jirasek
 */
class CohesiveInterfaceMaterialStatus : public StructuralMaterialStatus
{
public:
    /// Constructor
    CohesiveInterfaceMaterialStatus(int n, Domain * d, GaussPoint * g);
    /// Destructor
    virtual ~CohesiveInterfaceMaterialStatus();

    virtual void printOutputAt(FILE *file, TimeStep *tStep);

    // definition
    virtual const char *giveClassName() const { return "CohesiveInterfaceMaterialStatus"; }

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
    CohesiveInterfaceMaterial(int n, Domain * d);
    /// Destructor
    virtual ~CohesiveInterfaceMaterial() { };

    virtual int hasNonLinearBehaviour() { return 0; }
    virtual int hasMaterialModeCapability(MaterialMode mode) { return ( mode == _3dInterface ); }

    virtual const char *giveClassName() const { return "CohesiveInterfaceMaterial"; }
    virtual const char *giveInputRecordName() const { return _IFT_CohesiveInterfaceMaterial_Name; }

    virtual void give3dMaterialStiffnessMatrix(FloatMatrix &answer,
                                               MatResponseMode,
                                               GaussPoint *gp,
                                               TimeStep *tStep);

    double computeVolumetricStrain(GaussPoint *gp, TimeStep *tStep);

    virtual void giveRealStressVector(FloatArray &answer, GaussPoint *gp,
                                      const FloatArray &reducedStrain, TimeStep *tStep);

    virtual void  giveStiffnessMatrix(FloatMatrix &answer,
                                      MatResponseMode mode,
                                      GaussPoint *gp,
                                      TimeStep *tStep);

    virtual int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep);

    virtual void giveThermalDilatationVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep);

    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual void giveInputRecord(DynamicInputRecord &input);

    virtual MaterialStatus *CreateStatus(GaussPoint *gp) const { return new CohesiveInterfaceMaterialStatus(1, FEMComponent :: domain, gp); }

protected:
    // Overloaded to use specialized versions of these services possibly implemented by linearElastic member
    void give3dInterfaceMaterialStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode,
                                                GaussPoint *gp, TimeStep *tStep);
};
} // namespace oofem
#endif
