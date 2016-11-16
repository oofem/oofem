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

#ifndef simplecrosssection_h
#define simplecrosssection_h

#include "../sm/CrossSections/structuralcrosssection.h"
#include "../sm/Materials/structuralmaterial.h"
#include "floatarray.h"
#include "floatmatrix.h"

///@name Input fields for SimpleCrossSection
//@{
#define _IFT_SimpleCrossSection_Name "simplecs"
#define _IFT_SimpleCrossSection_thick "thick"
#define _IFT_SimpleCrossSection_width "width"
#define _IFT_SimpleCrossSection_area "area"
#define _IFT_SimpleCrossSection_iy "iy" ///< Inertia moment y
#define _IFT_SimpleCrossSection_iz "iz" ///< Inertia moment z
#define _IFT_SimpleCrossSection_ik "ik" ///< Torsion moment x
#define _IFT_SimpleCrossSection_shearcoeff "beamshearcoeff"
#define _IFT_SimpleCrossSection_shearareay "shearareay" ///< Shear area y direction
#define _IFT_SimpleCrossSection_shearareaz "shearareaz" ///< Shear area z direction
#define _IFT_SimpleCrossSection_drillStiffness "drillstiffness" ///< Penalty term for drilling stiffness.
#define _IFT_SimpleCrossSection_MaterialNumber "material" ///< Material number for the bulk material
#define _IFT_SimpleCrossSection_directorx "directorx"
#define _IFT_SimpleCrossSection_directory "directory"
#define _IFT_SimpleCrossSection_directorz "directorz"
//@}

namespace oofem {
/**
 * Class implementing "simple" cross section model in finite element problem.
 * A cross section  is an attribute of a domain. It is usually also attribute of many
 * elements.
 *
 * The simple cross section implementation does not perform any integration over cross-section volume,
 * it represents a cross section model, where the whole cross section model is represented by single integration point.
 * and therefore all requests for characteristic contributions (stiffness) and for real stress computations are simply
 * passed to parent StructuralCrossSection class, which invokes corresponding material mode services.
 * Please note, that it is assumed that material model will support these material modes and provide
 * corresponding services for characteristic components and stress evaluation.
 * For description, how to incorporate more elaborate models of cross section, please read
 * base CrossSection documentation.
 *
 * The overloaded methods giveFullCharacteristicVector and giveFullCharacteristicVector add some additional support
 * for integrated cross section models - _3dShell, _3dBeam, _2dPlate and _2dBeam.
 *
 * This class also reads into its property dictionary necessary geometric cross section characteristics,
 * which are requested by particular material models. These parameters can be requested using get service and
 * include those defined by CrossSectionProperty.
 */
class OOFEM_EXPORT SimpleCrossSection : public StructuralCrossSection
{
public:
    /**
     * Constructor.
     * @param n Cross section number.
     * @param d Associated domain.
     */
    SimpleCrossSection(int n, Domain *d) : StructuralCrossSection(n, d) {
        materialNumber = 0;
        czMaterialNumber = 0;
    }

    virtual void giveRealStress_3d(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedStrain, TimeStep *tStep);
    virtual void giveRealStress_3dDegeneratedShell(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedStrain, TimeStep *tStep);
    virtual void giveRealStress_PlaneStrain(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedStrain, TimeStep *tStep);
    virtual void giveRealStress_PlaneStress(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedStrain, TimeStep *tStep);
    virtual void giveRealStress_1d(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedStrain, TimeStep *tStep);
    virtual void giveRealStress_Warping(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedStrain, TimeStep *tStep);

    virtual void giveStiffnessMatrix_3d(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);
    virtual void giveStiffnessMatrix_PlaneStress(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);
    virtual void giveStiffnessMatrix_PlaneStrain(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);
    virtual void giveStiffnessMatrix_1d(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);


    virtual void giveGeneralizedStress_Beam2d(FloatArray &answer, GaussPoint *gp, const FloatArray &generalizedStrain, TimeStep *tStep);
    virtual void giveGeneralizedStress_Beam3d(FloatArray &answer, GaussPoint *gp, const FloatArray &generalizedStrain, TimeStep *tStep);
    virtual void giveGeneralizedStress_Plate(FloatArray &answer, GaussPoint *gp, const FloatArray &generalizedStrain, TimeStep *tStep);
    virtual void giveGeneralizedStress_Shell(FloatArray &answer, GaussPoint *gp, const FloatArray &generalizedStrain, TimeStep *tStep);
    virtual void giveGeneralizedStress_MembraneRot(FloatArray &answer, GaussPoint *gp, const FloatArray &generalizedStrain, TimeStep *tStep);
    virtual void giveGeneralizedStress_PlateSubSoil(FloatArray &answer, GaussPoint *gp, const FloatArray &generalizedStrain, TimeStep *tStep);

    virtual void giveCharMaterialStiffnessMatrix(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);
    virtual bool isCharacteristicMtrxSymmetric(MatResponseMode mode);


    virtual void give3dDegeneratedShellStiffMtrx(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep);
    virtual void give2dBeamStiffMtrx(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);
    virtual void give3dBeamStiffMtrx(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);
    virtual void give2dPlateStiffMtrx(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);
    virtual void give3dShellStiffMtrx(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);
    virtual void giveMembraneRotStiffMtrx(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);
    virtual void give2dPlateSubSoilStiffMtrx(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);


    /**
     * Initializes receiver acording to object description stored in input record.
     * Calls CrossSection initializeFrom service and reads the values of
     * - 'thick' thickness
     * - 'width' width
     * - 'area' area
     * - 'iy' Moment of inertia around y
     * - 'iz' Moment of inertia around z
     * - 'ik' Torsion moment around x
     * - 'beamshearcoeff' Beam shear coefficient
     * @param ir Record to read off.
     */
    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual void giveInputRecord(DynamicInputRecord &input);

    virtual void createMaterialStatus(GaussPoint &iGP); // ES


    // identification and auxiliary functions
    virtual const char *giveClassName() const { return "SimpleCrossSection"; }
    virtual const char *giveInputRecordName() const { return _IFT_SimpleCrossSection_Name; }

    virtual double give(int aProperty, GaussPoint *gp);
    virtual double give(CrossSectionProperty a, GaussPoint *gp) { return CrossSection :: give(a, gp); }
    virtual double give(CrossSectionProperty a, const FloatArray &coords, Element *elem, bool local) { return CrossSection :: give(a, coords, elem, local); }
    virtual int giveIPValue(FloatArray &answer, GaussPoint *ip, InternalStateType type, TimeStep *tStep);
    virtual Material *giveMaterial(IntegrationPoint *ip);

    int giveMaterialNumber() const { return this->materialNumber; }
    void setMaterialNumber(int matNum) { this->materialNumber = matNum; }
    virtual int checkConsistency();
    virtual Interface *giveMaterialInterface(InterfaceType t, IntegrationPoint *ip);




    virtual void giveFirstPKStresses(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedFIncrement, TimeStep *tStep);
    virtual void giveCauchyStresses(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedFIncrement, TimeStep *tStep);
    virtual void giveEshelbyStresses(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedvF, TimeStep *tStep);
    virtual void giveStiffnessMatrix_dPdF(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);
    virtual void giveStiffnessMatrix_dCde(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);

    virtual void giveTemperatureVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep);

    virtual int packUnknowns(DataStream &buff, TimeStep *tStep, GaussPoint *gp);
    virtual int unpackAndUpdateUnknowns(DataStream &buff, TimeStep *tStep, GaussPoint *gp);
    virtual int estimatePackSize(DataStream &buff, GaussPoint *gp);

protected:
    int materialNumber;   // material number
    int czMaterialNumber; // cohesive zone material number
};
} // end namespace oofem
#endif // simplecrosssection_h
