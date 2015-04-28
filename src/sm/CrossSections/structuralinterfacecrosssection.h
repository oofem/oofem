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

#ifndef structuralinterfacecrosssection_h
#define structuralinterfacecrosssection_h

#include "../sm/Materials/InterfaceMaterials/structuralinterfacematerial.h"
#include "crosssection.h"
#include "classfactory.h"

///@name Input fields for StructuralInterfaceCrossSection
//@{
#define _IFT_StructuralInterfaceCrossSection_Name "interfacecs"
#define _IFT_StructuralInterfaceCrossSection_Material "material"
#define _IFT_StructuralInterfaceCrossSection_thickness "thickness"
//@}

namespace oofem {
class GaussPoint;
class Element;
class FloatArray;
class FloatMatrix;
typedef GaussPoint IntegrationPoint;

/**
 * Base class for all structural interface cross section models.
 * Keeps track of the interface material, the geometric thickness (for 2d elements)
 * and possibly (in the future) the integration rule (Gauss, Lobatto etc)
 */
class StructuralInterfaceCrossSection : public CrossSection
{
public:
    /**
     * Constructor. Creates cross section with given number, belonging to given domain.
     * @param n Cross section number.
     * @param d Domain to which new cross section will belong.
     */
    StructuralInterfaceCrossSection(int n, Domain * d) : CrossSection(n, d)
    {
        materialNum = 0;
        crossSectionType = CS_StructuralInterfaceCapability;
    }
    /// Destructor.
    virtual ~StructuralInterfaceCrossSection() { }

    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual int testCrossSectionExtension(CrossSectExtension ext) { return this->crossSectionType; }

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
     * @param reducedF Deformation gradient in reduced form.
     * @param tStep Current time step (most models are able to respond only when tStep is current time step).
     */
    //@{
    // Pass all calls to the material
    void giveFirstPKTraction_1d( FloatArray &answer, GaussPoint *gp, const FloatArray &jump, const FloatMatrix &F, TimeStep *tStep )
    { this->giveInterfaceMaterial()->giveFirstPKTraction_1d(answer, gp, jump, F, tStep); }
    void giveFirstPKTraction_2d( FloatArray &answer, GaussPoint *gp, const FloatArray &jump, const FloatMatrix &F, TimeStep *tStep )
    { this->giveInterfaceMaterial()->giveFirstPKTraction_2d(answer, gp, jump, F, tStep); }
    void giveFirstPKTraction_3d(FloatArray &answer, GaussPoint *gp, const FloatArray &jump, const FloatMatrix &F, TimeStep *tStep)
    { this->giveInterfaceMaterial()->giveFirstPKTraction_3d(answer, gp, jump, F, tStep); }


    void giveEngTraction_1d(FloatArray &answer, GaussPoint *gp, const FloatArray &jump, TimeStep *tStep)
    {
        this->giveInterfaceMaterial()->giveEngTraction_1d(answer, gp, jump, tStep);
    }
    void giveEngTraction_2d(FloatArray &answer, GaussPoint *gp, const FloatArray &jump, TimeStep *tStep)
    {
        this->giveInterfaceMaterial()->giveEngTraction_2d(answer, gp, jump, tStep);   
    }
    void giveEngTraction_3d(FloatArray &answer, GaussPoint *gp, const FloatArray &jump, TimeStep *tStep)
    {
        this->giveInterfaceMaterial()->giveEngTraction_3d(answer, gp, jump, tStep);
    }

    void give1dStiffnessMatrix_dTdj( FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep );
    void give2dStiffnessMatrix_dTdj( FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep );
    void give3dStiffnessMatrix_dTdj( FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep );


    void give1dStiffnessMatrix_Eng( FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep );
    void give2dStiffnessMatrix_Eng( FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep );
    void give3dStiffnessMatrix_Eng( FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep );
    //@}

    StructuralInterfaceMaterial *giveInterfaceMaterial();
    const FloatArray &giveTraction(IntegrationPoint *ip);

    virtual int checkConsistency();

    virtual int giveIPValue(FloatArray &answer, GaussPoint *ip, InternalStateType type, TimeStep *tStep);

    virtual int packUnknowns(DataStream &buff, TimeStep *tStep, GaussPoint *gp);
    virtual int unpackAndUpdateUnknowns(DataStream &buff, TimeStep *tStep, GaussPoint *gp);
    virtual int estimatePackSize(DataStream &buff, GaussPoint *gp);

    // identification and auxiliary functions
    virtual const char *giveClassName() const { return "StructuralInterfaceCrossSection"; }
    virtual const char *giveInputRecordName() const { return _IFT_StructuralInterfaceCrossSection_Name; }

    CrossSectExtension crossSectionType;
private:
    int materialNum;
};
} // end namespace oofem
#endif // structuralcrosssection_h
