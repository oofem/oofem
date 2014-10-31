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
#ifndef structmatsettable_h
#define structmatsettable_h

#include "floatarray.h"
#include "floatmatrix.h"

#include "../sm/Materials/structuralms.h"
#include "../sm/Materials/structuralmaterial.h"

#include "Materials/isolinearelasticmaterial.h"

///@name Input fields for StructuralMaterialSettable
//@{
#define _IFT_StructuralMaterialSettable_Name "structmatsettable"
#define _IFT_StructuralMaterialSettable_e "e"
#define _IFT_StructuralMaterialSettable_nu "nu"
//@}

namespace oofem {

/**
 * This class implements TODO
 */
class StructuralMaterialSettable : public StructuralMaterial
{
private:
    IsotropicLinearElasticMaterial *isoLE;

public:
    /// Constructor
    StructuralMaterialSettable(int n, Domain *d);
    /// Destructor
    virtual ~StructuralMaterialSettable();

    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual const char *giveClassName() const { return "StructuralMaterialSettable"; }
    virtual const char *giveInputRecordName() const { return _IFT_StructuralMaterialSettable_Name; }

    virtual void giveRealStressVector_3d(FloatArray &answer,
                              GaussPoint *gp,
                              const FloatArray &strainVector,
                              TimeStep *atTime);

    virtual void  give3dMaterialStiffnessMatrix(FloatMatrix &answer,
                                                MatResponseMode mode,
                                                GaussPoint *gp,
                                                TimeStep *atTime);

    virtual MaterialStatus *CreateStatus(GaussPoint *gp) const;

};
} // end namespace oofem
#endif // structmatsettable_h
