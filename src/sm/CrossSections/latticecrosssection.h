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

#ifndef latticecrosssection_h
#define latticecrosssection_h

#include "sm/Materials/LatticeMaterials/latticestructuralmaterial.h"
#include "crosssection.h"
#include "classfactory.h"

///@name Input fields for LatticeCrossSection
//@{
#define _IFT_LatticeCrossSection_Name "latticecs"
#define _IFT_LatticeCrossSection_Material "material"
#define _IFT_LatticeCrossSection_thickness "thickness"
//@}

namespace oofem {
class GaussPoint;
class Element;
class FloatArray;
class FloatMatrix;
typedef GaussPoint IntegrationPoint;

/**
 * Base class for all structural lattice cross section models.
 * Keeps track of the lattice material, the geometric thickness (for 2d elements)
 *
 */
class LatticeCrossSection : public CrossSection
{
public:
    /**
     * Constructor. Creates cross section with given number, belonging to given domain.
     * @param n Cross section number.
     * @param d Domain to which new cross section will belong.
     */
    LatticeCrossSection(int n, Domain *d) : CrossSection(n, d)
    {
        materialNum = 0;
        crossSectionType = CS_LatticeStructuralCapability;
    }
    /// Destructor.
    virtual ~LatticeCrossSection() { }

    void initializeFrom(InputRecord &ir) override;

    int testCrossSectionExtension(CrossSectExtension ext) override { return this->crossSectionType; }

    /**
     * Computes the real stress vector for given strain and integration point.
     * The service should use previously reached equilibrium history variables. Also
     * it should update temporary history variables in status according to newly reached state.
     * The temporary history variables are moved into equilibrium ones after global structure
     * equilibrium has been reached by iteration process.
     * Elements should always pass their requests to their cross section model, which
     * performs necessary integration over its volume and invokes necessary material
     * services for corresponding material model defined for given integration point.
     * @param answer Contains result.
     * @param gp Integration point.
     * @param reduthiscedF Deformation gradient in reduced form.
     * @param tStep Current time step (most models are able to respond only when tStep is current time step).
     */
    //@{
    // Pass all calls to the material

    double giveLatticeStress1d(double strain, GaussPoint *gp, TimeStep *tStep) const;

    FloatArrayF< 3 >giveLatticeStress2d(const FloatArrayF< 3 > &strain, GaussPoint *gp, TimeStep *tStep) const;

    FloatArrayF< 6 >giveLatticeStress3d(const FloatArrayF< 6 > &strain, GaussPoint *gp, TimeStep *tStep) const;

    FloatMatrixF< 1, 1 >give1dStiffnessMatrix(MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) const;

    FloatMatrixF< 3, 3 >give2dStiffnessMatrix(MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) const;

    FloatMatrixF< 6, 6 >give3dStiffnessMatrix(MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) const;
    //@}

    LatticeStructuralMaterial *giveLatticeMaterial() const;

    int checkConsistency() override;

    int giveIPValue(FloatArray &answer, GaussPoint *ip, InternalStateType type, TimeStep *tStep) override;

    Material *giveMaterial(IntegrationPoint *ip) const override;
    int giveMaterialNumber() const { return this->materialNum; }
    void setMaterialNumber(int matNum) { this->materialNum = matNum; }

    int packUnknowns(DataStream &buff, TimeStep *tStep, GaussPoint *gp) override;
    int unpackAndUpdateUnknowns(DataStream &buff, TimeStep *tStep, GaussPoint *gp) override;
    int estimatePackSize(DataStream &buff, GaussPoint *gp) override;

    // identification and auxiliary functions
    const char *giveClassName() const override { return "LatticeCrossSection"; }
    const char *giveInputRecordName() const override { return _IFT_LatticeCrossSection_Name; }

    virtual double give(int aProperty, GaussPoint *gp) const override;
    
    //    void giveCharMaterialStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep);


    CrossSectExtension crossSectionType;
private:
    int materialNum;
};
} // end namespace oofem
#endif // structuralcrosssection_h
