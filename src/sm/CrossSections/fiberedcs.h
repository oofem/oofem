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
    int numberOfFibers;     ///< Number of fibers.
    double thick; ///< Total thickness.
    double width; ///< Total width.
    double area;  ///< Total area.
    FloatArray fiberYcoords, fiberZcoords;

public:
    FiberedCrossSection(int n, Domain * d) : StructuralCrossSection(n, d),
        thick(0.), width(0.), area(-1.0), fiberYcoords(), fiberZcoords()
    { }

    virtual ~FiberedCrossSection()  { }

    void giveRealStress_3d(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedStrain, TimeStep *tStep) override;
    void giveRealStress_PlaneStrain(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedStrain, TimeStep *tStep) override;
    void giveRealStress_PlaneStress(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedStrain, TimeStep *tStep) override;
    void giveRealStress_1d(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedStrain, TimeStep *tStep) override;
    void giveRealStress_Warping(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedStrain, TimeStep *tStep) override;

    void giveStiffnessMatrix_3d(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) override;
    void giveStiffnessMatrix_PlaneStress(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) override;
    void giveStiffnessMatrix_PlaneStrain(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) override;
    void giveStiffnessMatrix_1d(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) override;

    void giveGeneralizedStress_Beam2d(FloatArray &answer, GaussPoint *gp, const FloatArray &generalizedStrain, TimeStep *tStep) override;
    void giveGeneralizedStress_Beam3d(FloatArray &answer, GaussPoint *gp, const FloatArray &generalizedStrain, TimeStep *tStep) override;
    void giveGeneralizedStress_Plate(FloatArray &answer, GaussPoint *gp, const FloatArray &generalizedStrain, TimeStep *tStep) override;
    void giveGeneralizedStress_Shell(FloatArray &answer, GaussPoint *gp, const FloatArray &generalizedStrain, TimeStep *tStep) override;
    void giveGeneralizedStress_MembraneRot(FloatArray &answer, GaussPoint *gp, const FloatArray &generalizedStrain, TimeStep *tStep) override;
    void giveGeneralizedStress_PlateSubSoil(FloatArray &answer, GaussPoint *gp, const FloatArray &generalizedStrain, TimeStep *tStep) override;

    void giveCharMaterialStiffnessMatrix(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) override;

    void give2dBeamStiffMtrx(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) override;
    void give3dBeamStiffMtrx(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) override;
    void give2dPlateStiffMtrx(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) override;
    void give3dShellStiffMtrx(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) override;
    void giveMembraneRotStiffMtrx(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) override;
    void give2dPlateSubSoilStiffMtrx(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) override;

    bool isCharacteristicMtrxSymmetric(MatResponseMode mode) override;
    double give(int aProperty, GaussPoint *gp) override
    {
        OOFEM_ERROR("not implemented yet");
        return 0.0;
    }
    FloatArray *imposeStressConstrainsOnGradient(GaussPoint *gp, FloatArray *gradientStressVector3d) override;
    FloatArray *imposeStrainConstrainsOnGradient(GaussPoint *gp, FloatArray *gradientStrainVector3d) override;

    int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep) override;

    /**
     * Method for computing 1d fiber stiffness matrix of receiver.
     * Default implementation computes 3d stiffness matrix using give3dMaterialStiffnessMatrix and
     * reduces it to 1d fiber stiffness using reduce method described above.
     * However, this reduction is quite time consuming and if it is possible,
     * it is recommended to overload this method and provide direct method for computing
     * particular stiffness matrix.
     * @param fiberMatrix Stiffness matrix.
     * @param mode Material response mode.
     * @param layerGp Integration point.
     * @param tStep Time step (most models are able to respond only when tStep is current time step).
     */
    void giveFiberMaterialStiffnessMatrix(FloatMatrix &fiberMatrix, MatResponseMode mode, GaussPoint *layerGp, TimeStep *tStep);

    double give(CrossSectionProperty a, GaussPoint *gp) override;

    // identification and auxiliary functions
    const char *giveInputRecordName() const override { return _IFT_FiberedCrossSection_Name; }
    const char *giveClassName() const override { return "FiberedCrossSection"; }
    IRResultType initializeFrom(InputRecord *ir) override;

    void createMaterialStatus(GaussPoint &iGP) override; // ES

    void printYourself() override;
    double computeIntegralThickWidth();
    MaterialMode giveCorrespondingSlaveMaterialMode(MaterialMode);
    GaussPoint *giveSlaveGaussPoint(GaussPoint *gp, int);

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
