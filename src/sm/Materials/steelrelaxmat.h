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

#ifndef steelrelaxmat_h
#define steelrealxmat_h

#include "../sm/Materials/structuralmaterial.h"
#include "../sm/Materials/structuralms.h"
#include "dictionary.h"
#include "floatarray.h"
#include "floatmatrix.h"

///@name Input fields for SteelRelaxMat
//@{
#define _IFT_SteelRelaxMat_Name "steelrelaxmat"
#define _IFT_SteelRelaxMat_E "e"
#define _IFT_SteelRelaxMat_reinfClass "reinfclass"
#define _IFT_SteelRelaxMat_k1 "k1"
#define _IFT_SteelRelaxMat_k2 "k2"
#define _IFT_SteelRelaxMat_rho1000 "rho1000"
#define _IFT_SteelRelaxMat_stiffnessFactor "stiffnessfactor"
#define _IFT_SteelRelaxMat_timeFactor "timefactor"
//#define _IFT_SteelRelaxMat_prestress "prestress"
#define _IFT_SteelRelaxMat_charStrength "charstrength"
#define _IFT_SteelRelaxMat_approach "approach"
#define _IFT_SteelRelaxMat_tolerance "tolerance"
#define _IFT_SteelRelaxMat_relRelaxBound "relrelaxbound"
//@}

namespace oofem {
class GaussPoint;
class Domain;

/**
 *
 */
class SteelRelaxMat : public StructuralMaterial
{
protected:

  double E;
  double k1;
  double k2;
  double rho1000;
  double mu;
  double timeFactor;
  double stiffnessFactor;
  //  double prestress;
  double charStrength;
  double tolerance;
  double relRelaxBound;

  enum approachType { Bazant_EC2, EquivTime_EC2} Approach;

  

public:
    SteelRelaxMat(int n, Domain * d);
    virtual ~SteelRelaxMat();

    virtual void giveRealStressVector_1d(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedE, TimeStep *tStep)
    { this->giveRealStressVector(answer, gp, reducedE, tStep); }

    virtual void giveRealStressVector(FloatArray &answer, GaussPoint *gp,
                                      const FloatArray &reducedStrain, TimeStep *tStep);

    virtual void give1dStressStiffMtrx(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep); 

    virtual void giveStressDependentPartOfStrainVector(FloatArray &answer, GaussPoint *gp, const FloatArray &totalStrain,
                                               TimeStep *tStep, ValueModeType mode);

    void computeStressRelaxationStrainVector(FloatArray &answer, GaussPoint *gp, const FloatArray &totalStrain, 
						     TimeStep *tStep, ValueModeType mode);

    void evalStressRelaxationAtConstStrain(double &answer, GaussPoint *gp, double dt);

    void evalIncrOfPrestressLossAtVarStrain(double &answer, GaussPoint *gp, double stress);

    void computeIncrOfPrestressLossAtVarStrain(double &answer, GaussPoint *gp, TimeStep *tStep, double stress);

    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual int hasNonLinearBehaviour() { return 1; }

    virtual int hasMaterialModeCapability(MaterialMode mode);

    virtual bool isCharacteristicMtrxSymmetric(MatResponseMode rMode) { return true; }

    virtual const char *giveInputRecordName() const { return _IFT_SteelRelaxMat_Name; }
    virtual const char *giveClassName() const { return "SteelRelaxMat"; }

    virtual MaterialStatus *CreateStatus(GaussPoint *gp) const;

    virtual int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep);
};

//=============================================================================


class SteelRelaxMatStatus : public StructuralMaterialStatus
{
protected:

  // for Bazant this internal variable is a cumulative viscous strain
  // while for eurocode it is a cumulative prestress loss

  double relaxIntVariable;
  double tempRelaxIntVariable;

  double prestress;

public:
    SteelRelaxMatStatus(int n, Domain * d, GaussPoint * g);
    virtual ~SteelRelaxMatStatus();



    void setTempRelaxIntVariable(double src) { tempRelaxIntVariable = src; }
    double giveTempRelaxIntVariable(void) { return tempRelaxIntVariable; }
    double giveRelaxIntVariable(void) { return relaxIntVariable; }

    void setPrestress(double src) { prestress = src; }
    double givePrestress(void) { return prestress; }


    virtual void printOutputAt(FILE *file, TimeStep *tStep);

    virtual void initTempStatus();

    virtual void updateYourself(TimeStep *tStep);

    virtual contextIOResultType saveContext(DataStream &stream, ContextMode mode, void *obj = NULL);
    virtual contextIOResultType restoreContext(DataStream &stream, ContextMode mode, void *obj = NULL);

    virtual const char *giveClassName() const { return "SteelRelaxMatStatus"; }
};
} // end namespace oofem
#endif // steelrelaxmat_h
