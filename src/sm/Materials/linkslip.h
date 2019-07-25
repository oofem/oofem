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

#ifndef linkslip_h
#define linkslip_h

#include "latticelinearelastic.h"
#include "latticematstatus.h"

///@name Input fields for LatticeSlip
//@{
#define _IFT_LinkSlip_Name "linkslip"
#define _IFT_LinkSlip_talpha "talpha"
#define _IFT_LinkSlip_eNormal "e"
#define _IFT_LinkSlip_alphaOne "a1"
#define _IFT_LinkSlip_alphaTwo "a2"
#define _IFT_LinkSlip_tauZero "t0"
#define _IFT_LinkSlip_localrandomtype "randomtype"
#define _IFT_LinkSlip_coefficientOfVariation "cov"
#define _IFT_LinkSlip_calpha "calpha"
//@}

namespace oofem {



class LinkSlipStatus : public StructuralMaterialStatus
{
protected:

  double plasticStrain;
  
  double tempPlasticStrain;
  


public:

    /// Constructor
    LinkSlipStatus(GaussPoint *g);
    /// Destructor
    ~LinkSlipStatus() {}
    
    double  giveTempPlasticStrain(){ return this->tempPlasticStrain; }
    double  givePlasticStrain(){ return this->plasticStrain; }
    void  letTempPlasticStrainBe(const double &v)
    { this->tempPlasticStrain = v; }

    void   printOutputAt(FILE *file, TimeStep *tStep);

    const char *giveClassName() const { return "LinkSlipStatus"; }

    virtual void initTempStatus();

    virtual void updateYourself(TimeStep *); // update after new equilibrium state reached

    contextIOResultType    saveContext(DataStream &stream, ContextMode mode, void *obj = NULL);

    contextIOResultType    restoreContext(DataStream &stream, ContextMode mode, void *obj = NULL);
};



/**
 * This class implements a slip model for concrete for lattice elements.
 */
class LinkSlip : public LinearElasticMaterial
    //
{
protected:

    ///Normal modulus
    double eNormalMean;

    ///Ratio of shear and normal modulus
    double alphaOne;

    ///Ratio of torsion and normal modulus
    double alphaTwo;

    ///Strength for slip component
    double tauZero;
      
    /// coefficient variation of the Gaussian distribution
    double coefficientOfVariation;

    /// flag which chooses between no distribution (0) and Gaussian distribution (1)
    double localRandomType;

    double cAlpha;

    double tAlphaMean;

public:

    /// Constructor
    LinkSlip(int n, Domain *d) : LinearElasticMaterial(n, d) { };


    LinkSlip(int n, Domain *d, double eNormalMean, double alphaOne, double alphaTwo);

    /// Destructor
    virtual ~LinkSlip();

    virtual const char *giveInputRecordName() const { return _IFT_LinkSlip_Name; }
    virtual const char *giveClassName() const { return "LinkSlip"; }

    virtual IRResultType initializeFrom(InputRecord *ir);

    //  virtual void computeStressIndependentStrainVector(FloatArray &answer, GaussPoint *gp, TimeStep *stepN, ValueModeType mode);

    virtual void  giveThermalDilatationVector(FloatArray &answer,  GaussPoint *gp,  TimeStep *tStep);


    virtual bool isCharacteristicMtrxSymmetric(MatResponseMode rMode) { return false; }


    virtual void give3dMaterialStiffnessMatrix(FloatMatrix &answer,
				 MatResponseMode rmode,
				 GaussPoint *gp,
				 TimeStep *atTime);

    
    virtual int hasMaterialModeCapability(MaterialMode mode);


    virtual Interface *giveInterface(InterfaceType);

    virtual void giveRealStressVector_3d(FloatArray &answer, GaussPoint *,
                                      const FloatArray &, TimeStep *);

    virtual MaterialStatus *CreateStatus(GaussPoint *gp) const;


    virtual double  give(int aProperty, GaussPoint *gp);



protected:

    virtual int giveIPValue(FloatArray &answer,
                            GaussPoint *gp,
                            InternalStateType type,
                            TimeStep *atTime);
};
} // end namespace oofem

#endif
