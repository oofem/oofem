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

#ifndef anisolinearelasticmaterial_h
#define anisolinearelasticmaterial_h

#include "linearelasticmaterial.h"
#include "dictionary.h"
#include "floatarray.h"
#include "floatmatrix.h"
#include "matconst.h"
#include "element.h"

///@name Input fields for AnisotropicLinearElasticMaterial
//@{
#define _IFT_AnisotropicLinearElasticMaterial_Name "anisole"
#define _IFT_AnisotropicLinearElasticMaterial_stiff "stiff"
#define _IFT_AnisotropicLinearElasticMaterial_talpha "talpha"
//@}

namespace oofem {
class GaussPoint;

/**
 * This class implements a general anisotropic linear elastic material in a finite
 * element problem.
 *
 * Tasks:
 * - Returning standard material stiffness matrix for 3d-case.
 *   according to current state determined by using data stored
 *   in Gausspoint, and local coordinate system defined in gp.
 * - Returning real stress state vector(tensor) at gauss point for 3d - case.
 */
class AnisotropicLinearElasticMaterial : public LinearElasticMaterial
{
protected:
    FloatMatrix stiffmat;
    FloatArray alpha;

public:

    AnisotropicLinearElasticMaterial(int n, Domain *d) : LinearElasticMaterial(n, d)
    {
        stiffmat.resize(6, 6);
        stiffmat.zero();
        alpha.resize(3);
        alpha.zero();
    }
    virtual ~AnisotropicLinearElasticMaterial()
    {}


    // identification and auxiliary functions
    virtual const char *giveInputRecordName() const { return _IFT_AnisotropicLinearElasticMaterial_Name; }
    virtual const char *giveClassName() const { return "AnisotropicLinearElasticMaterial"; }
    virtual IRResultType initializeFrom(InputRecord *ir);
    void giveInputRecord(DynamicInputRecord &input);

    // important functions
    virtual void give3dMaterialStiffnessMatrix(FloatMatrix &answer,
                                               MatResponseMode mode, GaussPoint *gp,
                                               TimeStep *tStep);
    virtual void giveThermalDilatationVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep);

    virtual MaterialStatus *CreateStatus(GaussPoint *gp) const;

protected:

    friend class CrossSection;
};
} // end namespace oofem
#endif // anisolinearelasticmaterial_h
