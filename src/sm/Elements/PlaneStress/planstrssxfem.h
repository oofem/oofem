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

#ifndef planstrssxfem_h
#define planstrssxfem_h

#include "sm/Elements/PlaneStress/planstrss.h"
#include "sm/xfem/xfemstructuralelementinterface.h"
#include "vtkxmlexportmodule.h"

#define _IFT_PlaneStress2dXfem_Name "planestress2dxfem"

namespace oofem {
/**
 * Temporary class for testing
 * in the usual case instead of PlaneStress2dXfem
 * there will be the standard PlaneStress2d
 */
class PlaneStress2dXfem : public PlaneStress2d, public XfemStructuralElementInterface, public VTKXMLExportModuleElementInterface
{
protected:
    void updateYourself(TimeStep *tStep) override;
    void postInitialize() override;

public:
    /// Constructor
    PlaneStress2dXfem(int n, Domain * d) : PlaneStress2d(n, d), XfemStructuralElementInterface(this), VTKXMLExportModuleElementInterface() { numberOfDofMans = 4; }
    /// Destructor
    virtual ~PlaneStress2dXfem() { }

    Interface *giveInterface(InterfaceType it) override;

    const char *giveInputRecordName() const override { return _IFT_PlaneStress2dXfem_Name; }
    const char *giveClassName() const override { return "PlaneStress2dXfem"; }
    int computeNumberOfDofs() override;
    void computeGaussPoints() override;

    void increaseNumGP();

    void computeNmatrixAt(const FloatArray &iLocCoord, FloatMatrix &answer) override;
    void computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int lowerIndx = 1, int upperIndx = ALL_STRAINS) override;
    void computeBHmatrixAt(GaussPoint *gp, FloatMatrix &answer) override;
    void giveDofManDofIDMask(int inode, IntArray &answer) const override;
    void computeConstitutiveMatrixAt(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *, TimeStep *tStep) override;
    void computeStressVector(FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep) override;
    void computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep) override;

    void computeDeformationGradientVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep) override;

    void giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord) override;

    void computeConsistentMassMatrix(FloatMatrix &answer, TimeStep *tStep, double &mass, const double *ipDensity = NULL) override
    { XfemStructuralElementInterface :: XfemElementInterface_computeConsistentMassMatrix(answer, tStep, mass, ipDensity); }

    Element_Geometry_Type giveGeometryType() const override;

#ifdef __OOFEG
    void drawRawGeometry(oofegGraphicContext &gc, TimeStep *tStep) override;
    //void drawDeformedGeometry(oofegGraphicContext &gc, TimeStep *tStep, UnknownType);
    void drawScalar(oofegGraphicContext &gc, TimeStep *tStep) override;
    //void drawSpecial(oofegGraphicContext &gc, TimeStep *tStep) override;
#endif

    void initializeFrom(InputRecord &ir) override;
    MaterialMode giveMaterialMode() override;
    void giveInputRecord(DynamicInputRecord &input) override;

    /// VTK Interface
    void giveCompositeExportData(std::vector< VTKPiece > &vtkPieces, IntArray &primaryVarsToExport, IntArray &internalVarsToExport, IntArray cellVarsToExport, TimeStep *tStep) override;

};
} // end namespace oofem
#endif
