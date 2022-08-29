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

#ifndef libeam3d2_h
#define libeam3d2_h

#include "sm/Elements/nlstructuralelement.h"
#include "sm/CrossSections/fiberedcs.h"

///@name Input fields for LIBeam3d2
//@{
#define _IFT_LIBeam3d2_Name "libeam3d2"
#define _IFT_LIBeam3d2_refnode "refnode"
//@}

namespace oofem {
/**
 * This class implements a 3-dimensional Linear Isoparametric
 * Mindlin theory beam element, with reduced integration.
 * VERY SIMPLE support for geometric nonlinearity is taken into account.
 * The central triad (which defines local coordinate system - corotational system)
 * is updated using average increment of rotational pseudovector (same as in Simo -Qu Voc element).
 * Nonlinearity in geometric equations is not taken into account,
 * the only nonlinearity accounted is due to rotation of local (corotational) frame.
 * Increments of element's strains are computed using increment of displacements in
 * updated local system and total strain by adding this increment to stored total deformation.
 * The local stiffness or nodal forces are transformed using updated centre triad to global system.
 */
class LIBeam3d2 : public NLStructuralElement, public FiberedCrossSectionInterface
{
private:
    double length;
    int referenceNode;

    /// Last equilibrium triad at the centre.
    FloatMatrix tc;
    /// Temporary triad at the centre.
    FloatMatrix tempTc;
    /// Time stamp of temporary centre triad.
    StateCounterType tempTcCounter;

public:
    LIBeam3d2(int n, Domain *d);
    virtual ~LIBeam3d2() { }

    void initializeFrom(InputRecord &ir) override;

    void updateYourself(TimeStep *tStep) override;

    void computeLumpedMassMatrix(FloatMatrix &answer, TimeStep *tStep) override;
    void computeMassMatrix(FloatMatrix &answer, TimeStep *tStep) override
    { computeLumpedMassMatrix(answer, tStep); }
    void computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep) override;
    bool computeGtoLRotationMatrix(FloatMatrix &answer) override;

    int testElementExtension(ElementExtension ext) override;

    int computeNumberOfDofs() override { return 12; }
    void giveDofManDofIDMask(int inode, IntArray &) const override;
    double computeVolumeAround(GaussPoint *gp) override;
    int giveLocalCoordinateSystem(FloatMatrix &answer) override;
    void computeStrainVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep) override;

    int computeGlobalCoordinates(FloatArray &answer, const FloatArray &lcoords) override;

    void computeStressVector(FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep) override;
    void computeConstitutiveMatrixAt(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) override;
    void computeConstitutiveMatrix_dPdF_At(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) override;


    // Fibered cross section support functions
    void FiberedCrossSectionInterface_computeStrainVectorInFiber(FloatArray &answer, const FloatArray &masterGpStrain,
                                                                 GaussPoint *slaveGp, TimeStep *tStep) override;

    Interface *giveInterface(InterfaceType it) override;

    // definition & identification
    const char *giveInputRecordName() const override { return _IFT_LIBeam3d2_Name; }
    const char *giveClassName() const override { return "LIBeam3d2"; }
    Element_Geometry_Type giveGeometryType() const override { return EGT_line_1; }

    void saveContext(DataStream &stream, ContextMode mode) override;
    void restoreContext(DataStream &stream, ContextMode mode) override;

#ifdef __OOFEG
    void drawRawGeometry(oofegGraphicContext &gc, TimeStep *tStep) override;
    void drawDeformedGeometry(oofegGraphicContext &gc, TimeStep *tStep, UnknownType) override;
    void drawScalar(oofegGraphicContext &gc, TimeStep *tStep) override;
#endif

    integrationDomain giveIntegrationDomain() const override { return _Line; }
    MaterialMode giveMaterialMode() override { return _3dBeam; }

protected:
    // edge load support
    void giveEdgeDofMapping(IntArray &answer, int iEdge) const override;
    double computeEdgeVolumeAround(GaussPoint *gp, int iEdge) override;
    int computeLoadLEToLRotationMatrix(FloatMatrix &answer, int iEdge, GaussPoint *gp) override;
    int computeLoadGToLRotationMtrx(FloatMatrix &answer) override;
    void computeBodyLoadVectorAt(FloatArray &answer, Load *load, TimeStep *tStep, ValueModeType mode) override;

    void computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int = 1, int = ALL_STRAINS) override;
    void computeNmatrixAt(const FloatArray &iLocCoord, FloatMatrix &answer) override;
    void computeGaussPoints() override;
    double computeLength() override;

    // nonlinearity
    void updateTempTriad(TimeStep *tStep);
    void computeSMtrx(FloatMatrix &answer, FloatArray &vec);
    void computeRotMtrx(FloatMatrix &answer, FloatArray &psi);
    double giveCurrentLength(TimeStep *tStep);
    void initForNewStep() override;
};
} // end namespace oofem
#endif // libeam3d2_h
