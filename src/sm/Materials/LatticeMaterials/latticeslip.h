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
 *               Copyright (C) 1993 - 2019   Borek Patzak
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

#ifndef latticeslip_h
#define latticeslip_h

#include "latticelinearelastic.h"
#include "latticematstatus.h"

///@name Input fields for LatticeSlip
//@{
#define _IFT_LatticeSlip_Name "latticeslip"
#define _IFT_LatticeSlip_talpha "talpha"
#define _IFT_LatticeSlip_e "e"
#define _IFT_LatticeSlip_a1 "a1"
#define _IFT_LatticeSlip_a2 "a2"
#define _IFT_LatticeSlip_t0 "t0"
//@}

namespace oofem {
class LatticeSlipStatus : public LatticeMaterialStatus
{
protected:


public:

    /// Constructor
    LatticeSlipStatus(GaussPoint *g);
    /// Destructor
    ~LatticeSlipStatus() {}

    void  letTempPlasticStrainBe(const FloatArray &v)
    { tempPlasticStrain = v; }

    void   printOutputAt(FILE *file, TimeStep *tStep);

    const char *giveClassName() const override { return "LatticeSlipStatus"; }

    virtual void initTempStatus() override;

    virtual void updateYourself(TimeStep *) override;

    void saveContext(DataStream &stream, ContextMode mode) override;

    void restoreContext(DataStream &stream, ContextMode mode) override;
};



/**
 * This class implements a slip model for interface elements between lattice and beam elements.
 */

class LatticeSlip : public LatticeLinearElastic
    //
{
protected:

    ///Normal modulus
    double eNormal;

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
    LatticeSlip(int n, Domain *d);


    /// Destructor
    virtual ~LatticeSlip();

    virtual const char *giveInputRecordName() const override { return _IFT_LatticeSlip_Name; }
    virtual const char *giveClassName() const override { return "LatticeSlip"; }

    virtual void initializeFrom(InputRecord &ir) override;

    virtual void  giveThermalDilatationVector(FloatArray &answer,  GaussPoint *gp,  TimeStep *tStep) override;


    virtual bool isCharacteristicMtrxSymmetric(MatResponseMode rMode) { return false; }


    virtual void give1dLatticeStiffMtrx(FloatMatrix &answer,
                                        MatResponseMode rmode,
                                        GaussPoint *gp,
                                        TimeStep *atTime) override;

    virtual void give2dLatticeStiffMtrx(FloatMatrix &answer,
                                        MatResponseMode rmode,
                                        GaussPoint *gp,
                                        TimeStep *atTime) override;

    virtual void give3dLatticeStiffMtrx(FloatMatrix &answer,
                                        MatResponseMode rmode,
                                        GaussPoint *gp,
                                        TimeStep *atTime) override;

    double computeDeltaDissipation(GaussPoint *gp,
                                   TimeStep *atTime);

    virtual int hasMaterialModeCapability(MaterialMode mode);


    virtual Interface *giveInterface(InterfaceType) override;

    virtual void giveRealStressVector(FloatArray &answer, GaussPoint *,
                                      const FloatArray &, TimeStep *) override;

    virtual void giveRealStressVector_Lattice2d(FloatArray &answer, GaussPoint *gp, const FloatArray &totalStrain, TimeStep *atTime) override { this->giveRealStressVector(answer, gp, totalStrain, atTime); }

    virtual void giveRealStressVector_Lattice3d(FloatArray &answer, GaussPoint *gp, const FloatArray &totalStrain, TimeStep *atTime) override { this->giveRealStressVector(answer, gp, totalStrain, atTime); }

    virtual void giveRandomParameters(FloatArray &param) override;


    virtual MaterialStatus *CreateStatus(GaussPoint *gp) const override;


    virtual double  give(int aProperty, GaussPoint *gp);



protected:

    virtual int giveIPValue(FloatArray &answer,
                            GaussPoint *gp,
                            InternalStateType type,
                            TimeStep *atTime) override;
};
} // end namespace oofem

#endif
