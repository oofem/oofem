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
#define _IFT_LinkSlip_kn "kn"
#define _IFT_LinkSlip_a1 "a1"
#define _IFT_LinkSlip_a2 "a2"
#define _IFT_LinkSlip_type "type"
#define _IFT_LinkSlip_t0 "t0"
#define _IFT_LinkSlip_s1 "s1"
#define _IFT_LinkSlip_s2 "s2"
#define _IFT_LinkSlip_s3 "s3"
#define _IFT_LinkSlip_tf "tf"
#define _IFT_LinkSlip_alpha "alpha"
//@}

namespace oofem {



class LinkSlipStatus : public StructuralMaterialStatus
{
protected:

  double plasticStrain;
  
  double tempPlasticStrain;

  double kappa;

  double tempKappa;

public:

    /// Constructor
    LinkSlipStatus(GaussPoint *g);
    /// Destructor
    ~LinkSlipStatus() {}
    
    double  giveTempPlasticStrain(){ return this->tempPlasticStrain; }

    double  givePlasticStrain(){ return this->plasticStrain; }

    double  giveKappa(){ return this->kappa; }

    void  letTempPlasticStrainBe(const double &v)
    { this->tempPlasticStrain = v; }

    void  letTempKappaBe(const double &v)
    { this->tempKappa = v; }

    void   printOutputAt(FILE *file, TimeStep *tStep);

    const char *giveClassName() const { return "LinkSlipStatus"; }

    virtual void initTempStatus();

    virtual void updateYourself(TimeStep *); // update after new equilibrium state reached

    void saveContext(DataStream &stream, ContextMode mode);

    void restoreContext(DataStream &stream, ContextMode mode);
};

/**
 * This class implements a slip model for concrete for lattice elements.
 */
class LinkSlip : public LinearElasticMaterial
    //
{
protected:

    ///Normal modulus
    double kNormal;

    ///Ratio of shear and normal modulus
    double alphaOne;

    int type;
    
    ///Strength for slip component
    double tauMax,tauFinal;
    
    double s1,s2,s3;

    double alpha;

public:

    LinkSlip(int n, Domain *d);

    /// Destructor
    virtual ~LinkSlip();

    virtual const char *giveInputRecordName() const { return _IFT_LinkSlip_Name; }
    virtual const char *giveClassName() const { return "LinkSlip"; }

    virtual IRResultType initializeFrom(InputRecord *ir);

    //  virtual void computeStressIndependentStrainVector(FloatArray &answer, GaussPoint *gp, TimeStep *stepN, ValueModeType mode);

    virtual void  giveThermalDilatationVector(FloatArray &answer,  GaussPoint *gp,  TimeStep *tStep);


    virtual bool isCharacteristicMtrxSymmetric(MatResponseMode rMode) { return false; }
    
    double evaluateBondStress(const double kappa) const;

    virtual void give3dMaterialStiffnessMatrix(FloatMatrix &answer,
				 MatResponseMode rmode,
				 GaussPoint *gp,
				 TimeStep *atTime);

    
    virtual int hasMaterialModeCapability(MaterialMode mode);


    virtual Interface *giveInterface(InterfaceType);

    virtual void giveRealStressVector_3d(FloatArray &answer, GaussPoint *,
                                      const FloatArray &, TimeStep *);

    virtual MaterialStatus *CreateStatus(GaussPoint *gp) const;


protected:

    virtual int giveIPValue(FloatArray &answer,
                            GaussPoint *gp,
                            InternalStateType type,
                            TimeStep *atTime);
};
} // end namespace oofem

#endif
