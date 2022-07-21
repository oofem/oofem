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

#include "sm/CrossSections/structuralcrosssection.h"
#include "sm/Materials/structuralmaterial.h"
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
#define _IFT_SimpleCrossSection_ik "ik" ///< Saint-Venant torsional constant
#define _IFT_SimpleCrossSection_shearcoeff "beamshearcoeff"
#define _IFT_SimpleCrossSection_shearareay "shearareay" ///< Shear area y direction
#define _IFT_SimpleCrossSection_shearareaz "shearareaz" ///< Shear area z direction
#define _IFT_SimpleCrossSection_drillStiffness "drillstiffness" ///< Penalty term for drilling stiffness.
#define _IFT_SimpleCrossSection_relDrillStiffness "reldrillstiffness" ///< Relative penalty term for drilling stiffness.
#define _IFT_SimpleCrossSection_drillType "drilltype" ///< Type of artificially added stiffnes for drilling DOFs.
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
protected:
    int materialNumber = 0;   ///< Material number
    int czMaterialNumber = 0; ///< Cohesive zone material number

public:
    /**
     * Constructor.
     * @param n Cross section number.
     * @param d Associated domain.
     */
    SimpleCrossSection(int n, Domain *d) : StructuralCrossSection(n, d) { }

    FloatArrayF< 6 >giveRealStress_3d(const FloatArrayF< 6 > &reducedStrain, GaussPoint *gp, TimeStep *tStep) const override;
    FloatArrayF< 6 >giveRealStress_3dDegeneratedShell(const FloatArrayF< 6 > &reducedStrain, GaussPoint *gp, TimeStep *tStep) const override;
    FloatArrayF< 4 >giveRealStress_PlaneStrain(const FloatArrayF< 4 > &reducedStrain, GaussPoint *gp, TimeStep *tStep) const override;
    FloatArrayF< 3 >giveRealStress_PlaneStress(const FloatArrayF< 3 > &reducedStrain, GaussPoint *gp, TimeStep *tStep) const override;
    FloatArrayF< 1 >giveRealStress_1d(const FloatArrayF< 1 > &reducedStrain, GaussPoint *gp, TimeStep *tStep) const override;
    FloatArrayF< 2 >giveRealStress_Warping(const FloatArrayF< 2 > &reducedStrain, GaussPoint *gp, TimeStep *tStep) const override;

    FloatMatrixF< 6, 6 >giveStiffnessMatrix_3d(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const override;
    FloatMatrixF< 3, 3 >giveStiffnessMatrix_PlaneStress(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const override;
    FloatMatrixF< 4, 4 >giveStiffnessMatrix_PlaneStrain(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const override;
    FloatMatrixF< 1, 1 >giveStiffnessMatrix_1d(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const override;

    FloatArrayF< 3 >giveGeneralizedStress_Beam2d(const FloatArrayF< 3 > &generalizedStrain, GaussPoint *gp, TimeStep *tStep) const override;
    FloatArrayF< 6 >giveGeneralizedStress_Beam3d(const FloatArrayF< 6 > &generalizedStrain, GaussPoint *gp, TimeStep *tStep) const override;
    FloatArrayF< 5 >giveGeneralizedStress_Plate(const FloatArrayF< 5 > &generalizedStrain, GaussPoint *gp, TimeStep *tStep) const override;
    FloatArrayF< 8 >giveGeneralizedStress_Shell(const FloatArrayF< 8 > &generalizedStrain, GaussPoint *gp, TimeStep *tStep) const override;
    FloatArrayF< 9 >giveGeneralizedStress_ShellRot(const FloatArrayF< 9 > &generalizedStrain, GaussPoint *gp, TimeStep *tStep) const override;
    FloatArrayF< 4 >giveGeneralizedStress_MembraneRot(const FloatArrayF< 4 > &generalizedStrain, GaussPoint *gp, TimeStep *tStep) const override;
    FloatArrayF< 3 >giveGeneralizedStress_PlateSubSoil(const FloatArrayF< 3 > &generalizedStrain, GaussPoint *gp, TimeStep *tStep) const override;

    void giveCharMaterialStiffnessMatrix(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) override;
    bool isCharacteristicMtrxSymmetric(MatResponseMode mode) const override;

    FloatMatrixF< 6, 6 >give3dDegeneratedShellStiffMtrx(MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) const override;
    FloatMatrixF< 3, 3 >give2dBeamStiffMtrx(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const override;
    FloatMatrixF< 6, 6 >give3dBeamStiffMtrx(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const override;
    FloatMatrixF< 5, 5 >give2dPlateStiffMtrx(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const override;
    FloatMatrixF< 8, 8 >give3dShellStiffMtrx(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const override;
    FloatMatrixF< 9, 9 >give3dShellRotStiffMtrx(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const override;
    FloatMatrixF< 4, 4 >giveMembraneRotStiffMtrx(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const override;
    FloatMatrixF< 3, 3 >give2dPlateSubSoilStiffMtrx(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const override;


    //@{
    virtual FloatArrayF< 9 >giveFirstPKStress_3d(const FloatArrayF< 9 > &reducedF, GaussPoint *gp, TimeStep *tStep) const override;
    virtual FloatArrayF< 5 >giveFirstPKStress_PlaneStrain(const FloatArrayF< 5 > &reducedF, GaussPoint *gp, TimeStep *tStep) const override;
    virtual FloatArrayF< 4 >giveFirstPKStress_PlaneStress(const FloatArrayF< 4 > &reducedF, GaussPoint *gp, TimeStep *tStep) const override;
    virtual FloatArrayF< 1 >giveFirstPKStress_1d(const FloatArrayF< 1 > &reducedF, GaussPoint *gp, TimeStep *tStep) const override;

    void giveCharMaterialStiffnessMatrix_dPdF(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) override;
    FloatMatrixF< 9, 9 >giveStiffnessMatrix_dPdF_3d(MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) const override;
    FloatMatrixF< 5, 5 >giveStiffnessMatrix_dPdF_PlaneStrain(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const override;
    FloatMatrixF< 4, 4 >giveStiffnessMatrix_dPdF_PlaneStress(MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) const override;
    FloatMatrixF< 1, 1 >giveStiffnessMatrix_dPdF_1d(MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) const override;



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
    void initializeFrom(InputRecord &ir) override;
    void giveInputRecord(DynamicInputRecord &input) override;

    void createMaterialStatus(GaussPoint &iGP) override;


    const char *giveClassName() const override { return "SimpleCrossSection"; }
    const char *giveInputRecordName() const override { return _IFT_SimpleCrossSection_Name; }

    double give(int aProperty, GaussPoint *gp) const override;
    double give(CrossSectionProperty a, GaussPoint *gp) const override { return CrossSection::give(a, gp); }
    double give(CrossSectionProperty a, const FloatArray &coords, Element *elem, bool local) const override { return CrossSection::give(a, coords, elem, local); }
    int giveIPValue(FloatArray &answer, GaussPoint *ip, InternalStateType type, TimeStep *tStep) override;
    Material *giveMaterial(IntegrationPoint *ip) const override;

    int giveMaterialNumber() const { return this->materialNumber; }
    void setMaterialNumber(int matNum) { this->materialNumber = matNum; }
    int checkConsistency() override;
    Interface *giveMaterialInterface(InterfaceType t, IntegrationPoint *ip) override;


    void giveCauchyStresses(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedFIncrement, TimeStep *tStep) override;
    void giveEshelbyStresses(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedvF, TimeStep *tStep) override;
    void giveStiffnessMatrix_dCde(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) override;

    void giveTemperatureVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep) const;

    int packUnknowns(DataStream &buff, TimeStep *tStep, GaussPoint *gp) override;
    int unpackAndUpdateUnknowns(DataStream &buff, TimeStep *tStep, GaussPoint *gp) override;
    int estimatePackSize(DataStream &buff, GaussPoint *gp) override;

    void saveContext(DataStream &stream, ContextMode mode) override;
    void restoreContext(DataStream &stream, ContextMode mode) override;
};
} // end namespace oofem
#endif // simplecrosssection_h
