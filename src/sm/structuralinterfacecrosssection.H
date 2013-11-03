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

#ifndef structuralcrosssection_h
#define structuralcrosssection_h

#include "crosssection.h"
#include "crosssectextension.h"
#include "structuralinterfacematerial.h"

namespace oofem {
class GaussPoint;
class Element;
class FloatArray;
class FloatMatrix;

/**
 * Abstract base class for all structural cross section models. It declares commons services provided by all
 * structural cross section models. The implementation of this services is left on derived classes,
 * which will implement cross section model dependent part. However, some general services are
 * implemented here.
 * For information, how to introduce integration points in cross section volume for
 * macro integration point, see @ref CrossSection reference manual.
 *
 * At structural level of cross section or constitutive models are introduced several stress/strain modes.
 * Full and reduced formats of stress/strain vectors are also introduced for convenience.
 * The full format includes all components, even if they are zero due to stress/strain mode nature,
 * but in the reduced format, only generally nonzero components are stored.
 * (full format must used only if absolutely necessary, to avoid wasting of space. It is used
 * by output routines to print results in general form). Methods for converting vectors between
 * full and reduced format are provided.
 * General full strain vector has one of the following forms:
 * -# strainVector3d {eps_x,eps_y,eps_z,gamma_yz,gamma_zx,gamma_xy}
 * -# For integrated cross section models (2d and 3d beams, plates and general shells)
 *    strainVectorShell {eps_x,eps_y,gamma_xy, kappa_x, kappa_y, kappa_xy, gamma_zx, gamma_zy}
 */
class StructuralInterfaceCrossSection : public CrossSection
{
public:
    /**
     * Constructor. Creates cross section with given number, belonging to given domain.
     * @param n Cross section number.
     * @param d Domain to which new cross section will belong.
     */
    StructuralInterfaceCrossSection(int n, Domain *d) : CrossSection(n, d) 
    { 
        materialNum = 0; 
        crossSectionType = CS_StructuralInterfaceCapability;
    }
    /// Destructor.
    virtual ~StructuralInterfaceCrossSection() { }

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
     * @param reducedStrain Strain vector in reduced form.
     * @param tStep Current time step (most models are able to respond only when tStep is current time step).
     */
    //@{

    // Pass all calls to the material
    void giveFirstPKTraction_1d(FloatArray &answer, GaussPoint *gp, const FloatArray &jump, const FloatArray &reducedF, TimeStep *tStep)
        { this->giveInterfaceMaterial()->giveFirstPKTraction_1d(answer, gp, jump, reducedF, tStep); }
    void giveFirstPKTraction_2d(FloatArray &answer, GaussPoint *gp, const FloatArray &jump, const FloatArray &reducedF, TimeStep *tStep) 
        { this->giveInterfaceMaterial()->giveFirstPKTraction_2d(answer, gp, jump, reducedF, tStep); }
    void giveFirstPKTraction_3d(FloatArray &answer, GaussPoint *gp, const FloatArray &jump, const FloatMatrix &F, TimeStep *tStep) 
        { this->giveInterfaceMaterial()->giveFirstPKTraction_3d(answer, gp, jump, F, tStep); }

    void giveEngTraction_1d(FloatArray &answer, GaussPoint *gp, const FloatArray &jump, TimeStep *tStep)
        { this->giveInterfaceMaterial()->giveEngTraction_1d(answer, gp, jump, tStep); }
    void giveEngTraction_2d(FloatArray &answer, GaussPoint *gp, const FloatArray &jump, TimeStep *tStep)
        { this->giveInterfaceMaterial()->giveEngTraction_2d(answer, gp, jump, tStep); }
    void giveEngTraction_3d(FloatArray &answer, GaussPoint *gp, const FloatArray &jump, TimeStep *tStep)
        { this->giveInterfaceMaterial()->giveEngTraction_3d(answer, gp, jump, tStep); }
    
    void give1dStiffnessMatrix_dTdj(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
        { this->giveInterfaceMaterial()->give1dStiffnessMatrix_dTdj(answer, rMode, gp, tStep); }
    void give2dStiffnessMatrix_dTdj(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
        { this->giveInterfaceMaterial()->give2dStiffnessMatrix_dTdj(answer, rMode, gp, tStep); }
    void give3dStiffnessMatrix_dTdj(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
        { this->giveInterfaceMaterial()->give3dStiffnessMatrix_dTdj(answer, rMode, gp, tStep); }

    void give1dStiffnessMatrix_Eng(FloatMatrix &answer,  MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
        { this->giveInterfaceMaterial()->give1dStiffnessMatrix_Eng(answer, rMode, gp, tStep); }
    void give2dStiffnessMatrix_Eng(FloatMatrix &answer,  MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
        { this->giveInterfaceMaterial()->give2dStiffnessMatrix_Eng(answer, rMode, gp, tStep); }
    void give3dStiffnessMatrix_Eng(FloatMatrix &answer,  MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
        { this->giveInterfaceMaterial()->give3dStiffnessMatrix_Eng(answer, rMode, gp, tStep); }
    //@}
    
     StructuralInterfaceMaterial *giveInterfaceMaterial();
     const FloatArray &giveTraction(IntegrationPoint *ip);

     virtual int StructuralInterfaceCrossSection :: checkConsistency();


    // identification and auxiliary functions
    virtual const char *giveClassName() const { return "StructuralInterfaceCrossSection"; }
    CrossSectExtension crossSectionType;
    

private:
    int materialNum;

};
} // end namespace oofem
#endif // structuralcrosssection_h

