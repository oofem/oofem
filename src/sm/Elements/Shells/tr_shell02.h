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

#ifndef tr_shell02_h
#define tr_shell02_h

#include "sm/Elements/structuralelement.h"
#include "zznodalrecoverymodel.h"
#include "sm/ErrorEstimators/zzerrorestimator.h"
#include "sm/Elements/Plates/dkt3d.h"
#include "sm/Elements/PlaneStress/trplanestressrotallman3d.h"
#include "spatiallocalizer.h"

#include <memory>

#define _IFT_TR_SHELL02_Name "tr_shell02"

namespace oofem {
/**
 * This class implements an triangular three-node shell finite element, composed of
 * dkt3d and trplanestressrotallman3d elements.
 * Each node has 6 degrees of freedom.
 *
 */
class TR_SHELL02 : public StructuralElement, public ZZNodalRecoveryModelInterface, public NodalAveragingRecoveryModelInterface, public ZZErrorEstimatorInterface, public SpatialLocalizerInterface
{
protected:
    /// Pointer to plate element.
    std :: unique_ptr< DKTPlate3d > plate;
    /// Pointer to membrane (plane stress) element.
    std :: unique_ptr< TrPlanestressRotAllman3d > membrane;
    /**
     * Element integraton rule (plate and membrane parts have their own integration rules)
     * this one used to integrate element error and perhaps can be (re)used for other putrposes.
     * Created on demand.
     */
    std :: unique_ptr< IntegrationRule > compositeIR;

    static IntArray loc_plate;
    static IntArray loc_membrane;

public:
    /// Constructor
    TR_SHELL02(int n, Domain * d);
    /// Destructor
    virtual ~TR_SHELL02() {}

    FEInterpolation *giveInterpolation() const override { return plate->giveInterpolation(); }

    int computeNumberOfDofs() override { return 18; }
    void giveDofManDofIDMask(int inode, IntArray &answer) const override
    { plate->giveDofManDofIDMask(inode, answer); }
    // definition & identification
    const char *giveInputRecordName() const override { return _IFT_TR_SHELL02_Name; }
    const char *giveClassName() const override { return "TR_SHELL02"; }
    void initializeFrom(InputRecord &ir) override;

    void giveCharacteristicVector(FloatArray &answer, CharType mtrx, ValueModeType mode, TimeStep *tStep) override;
    void giveCharacteristicMatrix(FloatMatrix &answer, CharType mtrx, TimeStep *tStep) override;
    double computeVolumeAround(GaussPoint *gp) override;
    bool giveRotationMatrix(FloatMatrix &answer) override;

    void updateYourself(TimeStep *tStep) override;
    void updateInternalState(TimeStep *tStep) override;
    void printOutputAt(FILE *file, TimeStep *tStep) override;
    void saveContext(DataStream &stream, ContextMode mode) override;
    void restoreContext(DataStream &stream, ContextMode mode) override;
    void postInitialize() override;
    void updateLocalNumbering(EntityRenumberingFunctor &f) override;
    void setCrossSection(int csIndx) override;
#ifdef __OOFEG
    void drawRawGeometry(oofegGraphicContext &gc, TimeStep *tStep) override;
    void drawDeformedGeometry(oofegGraphicContext &gc, TimeStep *tStep, UnknownType type) override;
    void drawScalar(oofegGraphicContext &gc, TimeStep *tStep) override;
#endif
    // the membrane and plate irules are same (chacked in initializeFrom)
    int giveDefaultIntegrationRule() const override { return plate->giveDefaultIntegrationRule(); }
    IntegrationRule *giveDefaultIntegrationRulePtr() override { return plate->giveDefaultIntegrationRulePtr(); }
    IntegrationRule *giveIntegrationRule(int i) override { return plate->giveIntegrationRule(i); }
    Element_Geometry_Type giveGeometryType() const override { return EGT_triangle_1; }
    integrationDomain giveIntegrationDomain() const override { return _Triangle; }
    MaterialMode giveMaterialMode() override { return _Unknown; }

    Interface *giveInterface(InterfaceType it) override;
    int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep) override;

    void NodalAveragingRecoveryMI_computeNodalValue(FloatArray &answer, int node,
                                                    InternalStateType type, TimeStep *tStep) override;

    IntegrationRule *ZZErrorEstimatorI_giveIntegrationRule() override;
    void ZZErrorEstimatorI_computeLocalStress(FloatArray &answer, FloatArray &sig) override;

    // SpatialLocalizerI
    void SpatialLocalizerI_giveBBox(FloatArray &bb0, FloatArray &bb1) override;

    int computeGlobalCoordinates(FloatArray &answer, const FloatArray &lcoords) override {
        return this->plate->computeGlobalCoordinates(answer, lcoords);
    }

    bool computeLocalCoordinates(FloatArray &answer, const FloatArray &gcoords) override {
        return this->plate->computeLocalCoordinates(answer, gcoords);
    }

protected:
    void computeBmatrixAt(GaussPoint *, FloatMatrix &, int = 1, int = ALL_STRAINS) override
    { OOFEM_ERROR("calling of this function is not allowed"); }
    void computeNmatrixAt(const FloatArray &iLocCoord, FloatMatrix &) override
    { OOFEM_ERROR("calling of this function is not allowed"); }

    /// @todo In time delete
protected:
    void computeGaussPoints() override
    {
        this->membrane->computeGaussPoints();
        this->plate->computeGaussPoints();
    }
    void computeStressVector(FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep) override
    { OOFEM_ERROR("calling of this function is not allowed"); }
    void computeConstitutiveMatrixAt(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) override
    { OOFEM_ERROR("calling of this funciton is not allowed"); }
    void computeBodyLoadVectorAt(FloatArray &answer, Load *forLoad, TimeStep *tStep, ValueModeType mode) override;

public:
    void computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep) override
    { OOFEM_ERROR("calling of this function is not allowed"); }
    void computeMassMatrix(FloatMatrix &answer, TimeStep *tStep) override
    { OOFEM_ERROR("calling of this function is not allowed"); }
    void giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord) override
    { OOFEM_ERROR("calling of this function is not allowed"); }
};
} // end namespace oofem
#endif
