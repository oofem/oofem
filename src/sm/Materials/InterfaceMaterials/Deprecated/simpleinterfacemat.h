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

#include "../structuralinterfacematerial.h"
#include "../structuralinterfacematerialstatus.h"

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

    virtual void printOutputAt(FILE *file, TimeStep *tStep);

    // definition
    virtual const char *giveClassName() const { return "SimpleInterfaceMaterialStatus"; }

    virtual void initTempStatus();
    virtual void updateYourself(TimeStep *tStep);

    const FloatArray &giveShearStressShift();
    void setTempShearStressShift(FloatArray newShearStressShift) { tempShearStressShift = newShearStressShift; }
    bool giveShearYieldingFlag(){return shearYieldingFlag;}
    void setShearYieldingFlag(bool sY){ shearYieldingFlag = sY;}

    virtual contextIOResultType saveContext(DataStream &stream, ContextMode mode, void *obj = NULL);
    virtual contextIOResultType restoreContext(DataStream &stream, ContextMode mode, void *obj = NULL);
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

    virtual int hasNonLinearBehaviour() { return 1; }
    virtual bool hasAnalyticalTangentStiffness() const { return true; }

    virtual const char *giveInputRecordName() const { return _IFT_SimpleInterfaceMaterial_Name; }
    virtual const char *giveClassName() const { return "SimpleInterfaceMaterial"; }

    virtual void giveEngTraction_3d(FloatArray &answer, GaussPoint *gp, const FloatArray &jump, TimeStep *tStep);
    virtual void give3dStiffnessMatrix_Eng(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep);

    virtual int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep);

    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual void giveInputRecord(DynamicInputRecord &input);

    virtual MaterialStatus *CreateStatus(GaussPoint *gp) const { return new SimpleInterfaceMaterialStatus(1, domain, gp); }
};
} // end namespace oofem
#endif // simpleinterfacemat_h
