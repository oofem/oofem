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

#ifndef fiberedcs_h
#define fiberedcs_h

#include "sm/CrossSections/structuralcrosssection.h"
#include "sm/Materials/structuralmaterial.h"
#include "element.h"
#include "dictionary.h"
#include "floatarray.h"
#include "floatmatrix.h"
#include "interface.h"

///@name Input fields for FiberedCrossSection
//@{
#define _IFT_FiberedCrossSection_Name "fiberedcs"
#define _IFT_FiberedCrossSection_nfibers "nfibers"
#define _IFT_FiberedCrossSection_fibermaterials "fibermaterials"
#define _IFT_FiberedCrossSection_thicks "thicks"
#define _IFT_FiberedCrossSection_widths "widths"
#define _IFT_FiberedCrossSection_fiberycentrecoords "fiberycentrecoords"
#define _IFT_FiberedCrossSection_fiberzcentrecoords "fiberzcentrecoords"
#define _IFT_FiberedCrossSection_thick "thick"
#define _IFT_FiberedCrossSection_width "width"
//@}

namespace oofem {
class GaussPoint;
class FiberedCrossSectionModelInterface;

/**
 * This class implements a fibered cross section in a finite element problem. A cross
 * section  is an attribute of a domain. It is usually also attribute of many
 * elements.
 *
 * The attribute 'propertyDictionary' contains all the properties of a
 * layered cross section, like thickness and width of each layer.
 * The attribute 'layerMaterials' contains an array of Materials corresponding
 * to each layer.
 *
 * it uses master - slave GaussPoint approach, where master gp has more slaves gp.
 * slave gp represent for each fiber material point. It's coordinate sections
 * contains y,z-coordinates from mid-section. the slaves are manageg completely
 * ( created, saved their context.,,,) from this class. Master gp only deletes
 * slaves in destructor.
 *
 * Tasks:
 * - Returning standard material stiffness marices (like 3d stress-strain, 2d plane ,
 *   plate, 3dbeam, 2d beam ..) according to current state determined by parameter
 *   StressMode by calling gp->material->GiveMaterialStiffnessMatrix (....) and by
 *   possible modifying returned matrix. (for example in layered mode approach
 *   each layer  is asked for 3dMatrialStiffnes and this is integrated for example
 *   over thickness for plate bending problems)
 * - Returning RealStress state in gauss point and for given Stress mode.
 * - Returning a properties of cross section like thickness or area.
 */
class FiberedCrossSection : public StructuralCrossSection
{
protected:
    IntArray fiberMaterials; ///< Material of each fiber.
    FloatArray fiberThicks; ///< Thickness for each fiber.
    FloatArray fiberWidths; ///< Width for each fiber.
    int numberOfFibers = 0;     ///< Number of fibers.
    double thick = 0.; ///< Total thickness.
    double width = 0.; ///< Total width.
    double area = 0.;  ///< Total area.
    FloatArray fiberYcoords, fiberZcoords;

public:
    FiberedCrossSection(int n, Domain * d) : StructuralCrossSection(n, d)
    { }

    FloatArrayF<6> giveRealStress_3d(const FloatArrayF<6> &reducedStrain, GaussPoint *gp, TimeStep *tStep) const override;
    FloatArrayF<4> giveRealStress_PlaneStrain(const FloatArrayF<4> &reducedStrain, GaussPoint *gp, TimeStep *tStep) const override;
    FloatArrayF<3> giveRealStress_PlaneStress(const FloatArrayF<3> &reducedStrain, GaussPoint *gp, TimeStep *tStep) const override;
    FloatArrayF<1> giveRealStress_1d(const FloatArrayF<1> &reducedStrain, GaussPoint *gp, TimeStep *tStep) const override;
    FloatArrayF<2> giveRealStress_Warping(const FloatArrayF<2> &reducedStrain, GaussPoint *gp, TimeStep *tStep) const override;


    FloatMatrixF<6,6> giveStiffnessMatrix_3d(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const override;
    FloatMatrixF<3,3> giveStiffnessMatrix_PlaneStress(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const override;
    FloatMatrixF<4,4> giveStiffnessMatrix_PlaneStrain(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const override;
    FloatMatrixF<1,1> giveStiffnessMatrix_1d(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const override;

    FloatArrayF<3> giveGeneralizedStress_Beam2d(const FloatArrayF<3> &generalizedStrain, GaussPoint *gp, TimeStep *tStep) const override;
    FloatArrayF<6> giveGeneralizedStress_Beam3d(const FloatArrayF<6> &generalizedStrain, GaussPoint *gp, TimeStep *tStep) const override;
    FloatArrayF<5> giveGeneralizedStress_Plate(const FloatArrayF<5> &generalizedStrain, GaussPoint *gp, TimeStep *tStep) const override;
    FloatArrayF<8> giveGeneralizedStress_Shell(const FloatArrayF<8> &generalizedStrain, GaussPoint *gp, TimeStep *tStep) const override;
    FloatArrayF<4> giveGeneralizedStress_MembraneRot(const FloatArrayF<4> &generalizedStrain, GaussPoint *gp, TimeStep *tStep) const override;
    FloatArrayF<3> giveGeneralizedStress_PlateSubSoil(const FloatArrayF<3> &generalizedStrain, GaussPoint *gp, TimeStep *tStep) const override;

    void giveCharMaterialStiffnessMatrix(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) override;

    FloatMatrixF<3,3> give2dBeamStiffMtrx(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const override;
    FloatMatrixF<6,6> give3dBeamStiffMtrx(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const override;
    FloatMatrixF<5,5> give2dPlateStiffMtrx(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const override;
    FloatMatrixF<8,8> give3dShellStiffMtrx(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const override;
    FloatMatrixF<4,4> giveMembraneRotStiffMtrx(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const override;
    FloatMatrixF<3,3> give2dPlateSubSoilStiffMtrx(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const override;

    bool isCharacteristicMtrxSymmetric(MatResponseMode mode) const override;
    double give(int aProperty, GaussPoint *gp) const override
    {
        OOFEM_ERROR("not implemented yet");
        return 0.0;
    }
    FloatArray *imposeStressConstrainsOnGradient(GaussPoint *gp, FloatArray *gradientStressVector3d) override;
    FloatArray *imposeStrainConstrainsOnGradient(GaussPoint *gp, FloatArray *gradientStrainVector3d) override;

    int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep) override;

    double give(CrossSectionProperty a, GaussPoint *gp) const override;

    // identification and auxiliary functions
    const char *giveInputRecordName() const override { return _IFT_FiberedCrossSection_Name; }
    const char *giveClassName() const override { return "FiberedCrossSection"; }
    void initializeFrom(InputRecord &ir) override;

    void createMaterialStatus(GaussPoint &iGP) override; // ES

    void printYourself() override;
    double computeIntegralThickWidth();
    static MaterialMode giveCorrespondingSlaveMaterialMode(MaterialMode);
    GaussPoint *giveSlaveGaussPoint(GaussPoint *gp, int) const;

    void saveIPContext(DataStream &stream, ContextMode mode, GaussPoint *gp) override;
    void restoreIPContext(DataStream &stream, ContextMode mode, GaussPoint *gp) override;

    int checkConsistency() override;

    int packUnknowns(DataStream &buff, TimeStep *tStep, GaussPoint *ip) override
    {
        OOFEM_ERROR("not implemented");
        return 0;
    }

    int unpackAndUpdateUnknowns(DataStream &buff, TimeStep *tStep, GaussPoint *ip) override
    {
        OOFEM_ERROR("not implemented");
        return 0;
    }

    int estimatePackSize(DataStream &buff, GaussPoint *ip) override
    {
        OOFEM_ERROR("not implemented");
        return 0;
    }

    void giveFirstPKStresses(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedFIncrement, TimeStep *tStep) override
    { OOFEM_ERROR("not implemented"); }
    void giveCauchyStresses(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedFIncrement, TimeStep *tStep) override
    { OOFEM_ERROR("not implemented"); }
    void giveStiffnessMatrix_dPdF(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) override
    { OOFEM_ERROR("not implemented"); }
    void giveStiffnessMatrix_dCde(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) override
    { OOFEM_ERROR("not implemented"); }
    
    Material *giveMaterial(IntegrationPoint *ip) const override;
};

/**
 * The element interface required by FiberedCrossSection.
 */
class FiberedCrossSectionInterface : public Interface
{
public:
    FiberedCrossSectionInterface() { }

    /**
     * Computes full 3d strain vector in element fiber. This function is necessary
     * if layered cross section is specified.
     * @param answer Full fiber strain vector.
     * @param masterGpStrain Strain vector at master gauss point.
     * @param slaveGp Slave integration point representing particular fiber.
     * @param tStep Time step.
     */
    virtual void FiberedCrossSectionInterface_computeStrainVectorInFiber(FloatArray &answer, const FloatArray &masterGpStrain,
                                                                         GaussPoint *slaveGp, TimeStep *tStep) = 0;
};
} // end namespace oofem
#endif // fiberedcs_h
