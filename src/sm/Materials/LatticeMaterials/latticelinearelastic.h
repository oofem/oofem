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

#ifndef latticelinearelastic_h
#define latticelinearelastic_h

#include "latticestructuralmaterial.h"
#include "cltypes.h"
#include "randommaterialext.h"
#include "strainvector.h"
#include "stressvector.h"
#include "latticematstatus.h"

///@name Input fields for LatticeLinearElastic
//@{
#define _IFT_LatticeLinearElastic_Name "latticelinearelastic"
#define _IFT_LatticeLinearElastic_talpha "talpha"
#define _IFT_LatticeLinearElastic_eNormal "e"
#define _IFT_LatticeLinearElastic_nu "n"
#define _IFT_LatticeLinearElastic_alphaOne "a1"
#define _IFT_LatticeLinearElastic_alphaTwo "a2"
#define _IFT_LatticeLinearElastic_localrandomtype "randomtype"
#define _IFT_LatticeLinearElastic_coefficientOfVariation "cov"
#define _IFT_LatticeLinearElastic_calpha "calpha"
//@}

namespace oofem {
/**
 * This class implements a local random linear elastic model for lattice elements.
 */
class LatticeLinearElastic : public LatticeStructuralMaterial, public RandomMaterialExtensionInterface
    //
{
protected:

    ///Normal modulus
    double eNormalMean;
    ///Ratio of shear and normal modulus
    double alphaOne;
    ///Ratio of torsion and normal modulus
    double alphaTwo;

    double nu;

    /// coefficient variation of the Gaussian distribution
    double coefficientOfVariation;

    /// flag which chooses between no distribution (0) and Gaussian distribution (1)
    double localRandomType;

    double cAlpha;

    double tAlphaMean;

public:

    /// Constructor
    LatticeLinearElastic(int n, Domain *d) : LatticeStructuralMaterial(n, d), RandomMaterialExtensionInterface() { };


    LatticeLinearElastic(int n, Domain *d, double eNormalMean, double alphaOne, double alphaTwo);

    /// Destructor
    virtual ~LatticeLinearElastic();

    const char *giveInputRecordName() const override { return _IFT_LatticeLinearElastic_Name; }
    const char *giveClassName() const override { return "LatticeLinearElastic"; }

    IRResultType initializeFrom(InputRecord *ir) override;

  
    void  giveThermalDilatationVector(FloatArray &answer,  GaussPoint *gp,  TimeStep *tStep) override;


    bool isCharacteristicMtrxSymmetric(MatResponseMode rMode) { return false; }

    
    virtual FloatArrayF<3> giveLatticeStress2d(const FloatArrayF<3> &strain, GaussPoint *gp, TimeStep *tStep) override;
    
    virtual FloatArrayF<6> giveLatticeStress3d(const FloatArrayF<6> &strain, GaussPoint *gp, TimeStep *tStep) override;

    virtual FloatMatrixF<1,1> give1dLatticeStiffnessMatrix(MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) const override;
    
    virtual FloatMatrixF<3,3> give2dLatticeStiffnessMatrix(MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) const override;
    
    virtual FloatMatrixF<6,6> give3dLatticeStiffnessMatrix(MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) const override;

    
    int hasMaterialModeCapability(MaterialMode mode);


    Interface *giveInterface(InterfaceType) override;

    virtual void giveRandomParameters(FloatArray &param);

    MaterialStatus *CreateStatus(GaussPoint *gp) const override;

    MaterialStatus *giveStatus(GaussPoint *gp) const override;

    double  give(int aProperty, GaussPoint *gp) const;

protected:

    int giveIPValue(FloatArray &answer,
                    GaussPoint *gp,
                    InternalStateType type,
                    TimeStep *atTime) override;
};


class LatticeLinearElasticMaterialStatus : public LatticeMaterialStatus
{
protected:

public:

    /// Constructor
    LatticeLinearElasticMaterialStatus(GaussPoint *g);
    /// Destructor
    ~LatticeLinearElasticMaterialStatus() {}

    void   printOutputAt(FILE *file, TimeStep *tStep);

    const char *giveClassName() const { return "LatticeLinearElasticMaterialStatus"; }
};
} // end namespace oofem

#endif
