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

#ifndef TRPLANSTRSSXFEM_H_
#define TRPLANSTRSSXFEM_H_

#include "sm/Elements/PlaneStress/trplanstrss.h"
#include "sm/xfem/xfemstructuralelementinterface.h"
#include "vtkxmlexportmodule.h"


#define _IFT_TrPlaneStress2dXFEM_Name "trplanestress2dxfem"

#define _IFT_TrPlaneStress2dXFEM_RegCoeff "reg_coeff"
#define _IFT_TrPlaneStress2dXFEM_RegCoeffTol "reg_coeff_tol"

namespace oofem {
/**
 * 3-node triangle with XFEM kinematics
 * @author Erik Svenning
 */
class TrPlaneStress2dXFEM : public TrPlaneStress2d, public XfemStructuralElementInterface, public VTKXMLExportModuleElementInterface
{
protected:
    void updateYourself(TimeStep *tStep) override;
    void postInitialize() override;

    double mRegCoeff, mRegCoeffTol;

public:
    TrPlaneStress2dXFEM(int n, Domain * d) : TrPlaneStress2d(n, d), XfemStructuralElementInterface(this), VTKXMLExportModuleElementInterface() { numberOfDofMans = 3; mRegCoeff = 1.0e-6; mRegCoeffTol = 1.0e-6;}

    virtual ~TrPlaneStress2dXFEM();

    int checkConsistency() override;

    int testElementExtension(ElementExtension ext) override { return ( ( ext == Element_EdgeLoadSupport ) ? 1 : 0 ); }

    Interface *giveInterface(InterfaceType it) override;

    const char *giveInputRecordName() const override { return _IFT_TrPlaneStress2dXFEM_Name; }
    const char *giveClassName() const override { return "TrPlaneStress2dXFEM"; }

    int computeNumberOfDofs() override;
    void computeGaussPoints() override;
    void computeNmatrixAt(const FloatArray &iLocCoord, FloatMatrix &answer) override;
    void computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int lowerIndx = 1, int upperIndx = ALL_STRAINS) override;
    void computeBHmatrixAt(GaussPoint *gp, FloatMatrix &answer) override;
    void giveDofManDofIDMask(int inode, IntArray &answer) const override;
    void computeConstitutiveMatrixAt(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) override;
    void computeStressVector(FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep) override;
    void computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep) override;

    void computeDeformationGradientVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep) override;

    void giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord) override;

    void computeConsistentMassMatrix(FloatMatrix &answer, TimeStep *tStep, double &mass, const double *ipDensity = NULL) override { XfemStructuralElementInterface :: XfemElementInterface_computeConsistentMassMatrix(answer, tStep, mass, ipDensity); }

    Element_Geometry_Type giveGeometryType() const override;

#ifdef __OOFEG
    // TODO: Implement OOFEG functions
    void drawRawGeometry(oofegGraphicContext &gc, TimeStep *tStep) override;
    //void drawDeformedGeometry(oofegGraphicContext &gc, TimeStep *tStep, UnknownType);
    void drawScalar(oofegGraphicContext &gc, TimeStep *tStep) override;
    //void drawSpecial(oofegGraphicContext &gc, TimeStep *tStep) override;
#endif

    void initializeFrom(InputRecord &ir) override;
    MaterialMode giveMaterialMode() override;
    void giveInputRecord(DynamicInputRecord &input) override;

    void computeField(ValueModeType mode, TimeStep *tStep, const FloatArray &lcoords, FloatArray &answer) override;
    void giveElementDofIDMask(IntArray &answer) const override;

    /// VTK Interface
    void giveCompositeExportData(std::vector< VTKPiece > &vtkPieces, IntArray &primaryVarsToExport, IntArray &internalVarsToExport, IntArray cellVarsToExport, TimeStep *tStep) override;

};
} /* namespace oofem */
#endif /* TRPLANSTRSSXFEM_H_ */
