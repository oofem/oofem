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

#ifndef simpleinterfacemat_h
#define simpleinterfacemat_h

#include "sm/Materials/InterfaceMaterials/structuralinterfacematerial.h"
#include "sm/Materials/InterfaceMaterials/structuralinterfacematerialstatus.h"

///@name Input fields for SimpleInterfaceMaterial
//@{
#define _IFT_SimpleInterfaceMaterial_Name "simpleintermat"
#define _IFT_SimpleInterfaceMaterial_kn "kn"
#define _IFT_SimpleInterfaceMaterial_ks "ks"
#define _IFT_SimpleInterfaceMaterial_frictCoeff "fc"
#define _IFT_SimpleInterfaceMaterial_stiffCoeff "stiffcoeff"
#define _IFT_SimpleInterfaceMaterial_normalClearance "normalclearance"
//@}

namespace oofem {
/**
 * This class implements associated Material Status to SimpleInterfaceMaterial.
 */
class SimpleInterfaceMaterialStatus : public StructuralInterfaceMaterialStatus
{
protected:
    bool shearYieldingFlag;
    FloatArray shearStressShift, tempShearStressShift;

public:
    /// Constructor
    SimpleInterfaceMaterialStatus(int n, Domain * d, GaussPoint * g);
    /// Destructor
    virtual ~SimpleInterfaceMaterialStatus();

    void printOutputAt(FILE *file, TimeStep *tStep) override;

    const char *giveClassName() const override { return "SimpleInterfaceMaterialStatus"; }

    void initTempStatus() override;
    void updateYourself(TimeStep *tStep) override;

    const FloatArray &giveShearStressShift();
    void setTempShearStressShift(FloatArray newShearStressShift) { tempShearStressShift = newShearStressShift; }
    bool giveShearYieldingFlag(){return shearYieldingFlag;}
    void setShearYieldingFlag(bool sY){ shearYieldingFlag = sY;}

    void saveContext(DataStream &stream, ContextMode mode) override;
    void restoreContext(DataStream &stream, ContextMode mode) override;
};


/**
 * Base class representing general isotropic damage model.
 * It is based on isotropic damage concept, assuming that damage evolution law
 * is postulated in explicit form, relatin damage parameter (omega) to scalar measure
 * of the largest strain level ever reached in material (kappa).
 */
class SimpleInterfaceMaterial : public StructuralInterfaceMaterial
{
protected:
  double kn, ks;
    double stiffCoeff;
    double frictCoeff;
    /// Normal distance which needs to be closed when interface element should act in compression (distance is 0 by default).
    double normalClearance;

public:
    /// Constructor
    SimpleInterfaceMaterial(int n, Domain * d);
    /// Destructor
    virtual ~SimpleInterfaceMaterial();

    bool hasAnalyticalTangentStiffness() const override { return true; }

    const char *giveInputRecordName() const override { return _IFT_SimpleInterfaceMaterial_Name; }
    const char *giveClassName() const override { return "SimpleInterfaceMaterial"; }

    void giveEngTraction_3d(FloatArray &answer, GaussPoint *gp, const FloatArray &jump, TimeStep *tStep) override;
    void give3dStiffnessMatrix_Eng(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) override;

    int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep) override;

    IRResultType initializeFrom(InputRecord *ir) override;
    void giveInputRecord(DynamicInputRecord &input) override;

    MaterialStatus *CreateStatus(GaussPoint *gp) const override { return new SimpleInterfaceMaterialStatus(1, domain, gp); }
};
} // end namespace oofem
#endif // simpleinterfacemat_h
