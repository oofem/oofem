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

#ifndef hyperelasticmaterial_h
#define hyperelasticmaterial_h

#include "structuralmaterial.h"
#include "structuralms.h"

///@name Input fields for HyperElasticMaterial
//@{
#define _IFT_HyperElasticMaterial_Name "hyperelmat"
#define _IFT_HyperElasticMaterial_k "k"
#define _IFT_HyperElasticMaterial_g "g"
//@}

namespace oofem {
/**
 * This class implements associated MaterialStatus for HyperElasticMaterial.
 */
class HyperElasticMaterialStatus : public StructuralMaterialStatus
{
public:
    /// Constructor
    HyperElasticMaterialStatus(int n, Domain *d, GaussPoint *g);
    /// Destructor
    virtual ~HyperElasticMaterialStatus();

    virtual void printOutputAt(FILE *file, TimeStep *tStep);

    // definition
    virtual const char *giveClassName() const { return "HyperElasticMaterialStatus"; }
    virtual classType giveClassID() const { return HyperElasticMaterialStatusClass; }

    virtual void initTempStatus();
    virtual void updateYourself(TimeStep *tStep);
};


class HyperElasticMaterial : public StructuralMaterial
{
protected:
    double K; ///< Bulk modulus.
    double G; ///< Shear modulus.

public:
    HyperElasticMaterial(int n, Domain *d);

    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual void give3dMaterialStiffnessMatrix(FloatMatrix &answer,
                                               MatResponseMode mode, GaussPoint *gp,
                                               TimeStep *tStep);


    virtual void giveRealStressVector_3d(FloatArray &answer, GaussPoint *gp,
                                      const FloatArray &reducedStrain, TimeStep *tStep);
    
    virtual MaterialStatus *CreateStatus(GaussPoint *gp) const;

    virtual int hasMaterialModeCapability(MaterialMode);
    virtual const char *giveInputRecordName() const { return _IFT_HyperElasticMaterial_Name; }
    virtual const char *giveClassName() const { return "HyperElasticMaterial"; }
    virtual classType giveClassID() const { return HyperElasticMaterialClass; }
};
} // end namespace oofem
#endif
