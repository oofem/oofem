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

#ifndef beam2d_h
#define beam2d_h

#include "sm/Elements/Beams/beambaseelement.h"
#include "sm/CrossSections/layeredcrosssection.h"
#include "dofmanager.h"

///@name Input fields for Beam2d
//@{
#define _IFT_Beam2d_Name "beam2d"
#define _IFT_Beam2d_dofstocondense "dofstocondense"
//@}

namespace oofem {
class FEI2dLineLin;
class FEI2dLineHermite;

/**
 * This class implements a 2-dimensional beam element
 * with cubic lateral displacement, quadratic rotations,
 * and linear longitudinal displacements and geometry.
 * This is an exact displacement approximation for beam with no
 * nonnodal loading.
 *
 * This class is not derived from linear beam or truss element, because it does not support
 * any material nonlinearities (if should, stiffness must be integrated)
 */
class Beam2d : public BeamBaseElement, public LayeredCrossSectionInterface
{
protected:
    double kappa, pitch, length;
    /**
     * Ghost nodes are used to introduce additional DOFs at element.
     * These are needed as we actually do not want to condense selected DOFs, but rather
     * allocate an extra equation to these. This allows to get cooresponding DOFs directly from
     * the global system, avoiding the need to postprocess local displacements at element.
     */
    DofManager *ghostNodes [ 2 ];
    /// number of condensed DOFs
    int numberOfCondensedDofs;

    static FEI2dLineLin interp_geom;
    static FEI2dLineHermite interp_beam;

public:
    Beam2d(int n, Domain *aDomain);
    virtual ~Beam2d();

    void computeConsistentMassMatrix(FloatMatrix &answer, TimeStep *tStep, double &mass, const double *ipDensity = NULL) override;
    void computeInitialStressMatrix(FloatMatrix &answer, TimeStep *tStep) override;
    void computeLumpedInitialStressMatrix(FloatMatrix &answer, TimeStep *tStep) override;

    void computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep) override;
    int giveLocalCoordinateSystem(FloatMatrix &answer) override;
    void giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord = 0) override;
    void giveEndForcesVector(FloatArray &answer, TimeStep *tStep);

    int testElementExtension(ElementExtension ext) override { return ( ext == Element_EdgeLoadSupport ); }

    Interface *giveInterface(InterfaceType) override;

    FEInterpolation *giveInterpolation() const override;
    FEInterpolation *giveInterpolation(DofIDItem id) const override { return nullptr; }

    int computeNumberOfDofs() override { return 6; }
    int computeNumberOfGlobalDofs() override { return 6 + this->numberOfCondensedDofs; }
    void giveDofManDofIDMask(int inode, IntArray &) const override;
    int giveNumberOfInternalDofManagers() const override { return ( ghostNodes [ 0 ] != nullptr ) + ( ghostNodes [ 1 ] != nullptr ); }
    DofManager *giveInternalDofManager(int i) const override {
        if ( i == 1 ) {
            if ( ghostNodes [ 0 ] ) { return ghostNodes [ 0 ]; } else { return ghostNodes [ 1 ]; }
        } else if ( i == 2 ) { // i==2
            return ghostNodes [ 1 ];
        } else {
            OOFEM_ERROR("No such DOF available on Element %d", number);
        }
    }
    void giveInternalDofManDofIDMask(int i, IntArray &answer) const override {
        if ( i == 1 ) {
            if ( ghostNodes [ 0 ] ) {
                ghostNodes [ 0 ]->giveCompleteMasterDofIDArray(answer);
            } else {
                ghostNodes [ 1 ]->giveCompleteMasterDofIDArray(answer);
            }
        } else if ( i == 2 ) { // i==2
            ghostNodes [ 1 ]->giveCompleteMasterDofIDArray(answer);
        } else {
            OOFEM_ERROR("No such DOF available on Element %d", number);
        }
    }

    void giveBoundaryLocationArray(IntArray &locationArray, const IntArray &bNodes, const UnknownNumberingScheme &s, IntArray *dofIds = NULL) override {
      giveLocationArray (locationArray, s, dofIds);
    }

    void giveBoundaryLocationArray(IntArray &locationArray, const IntArray &bNodes, const IntArray &dofIDMask, const UnknownNumberingScheme &s, IntArray *dofIds = NULL) override {
      giveLocationArray (locationArray, dofIDMask, s, dofIds);
    }

    double computeVolumeAround(GaussPoint *gp) override;
    void printOutputAt(FILE *file, TimeStep *tStep) override;

    Element_Geometry_Type giveGeometryType() const override {return EGT_line_1;}

    const char *giveClassName() const override { return "Beam2d"; }
    const char *giveInputRecordName() const override { return _IFT_Beam2d_Name; }
    void initializeFrom(InputRecord &ir) override;

#ifdef __OOFEG
    void drawRawGeometry(oofegGraphicContext &gc, TimeStep *tStep) override;
    void drawDeformedGeometry(oofegGraphicContext & gc, TimeStep * tStep, UnknownType) override;
#endif

    void computeStrainVectorInLayer(FloatArray &answer, const FloatArray &masterGpStrain,
                                    GaussPoint *masterGp, GaussPoint *slaveGp, TimeStep *tStep) override;

    int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep) override;

protected:
    void computeBoundaryEdgeLoadVector(FloatArray &answer, BoundaryLoad *load, int edge, CharType type, ValueModeType mode, TimeStep *tStep, bool global=true) override;
    void computeBmatrixAt(GaussPoint *, FloatMatrix &, int = 1, int = ALL_STRAINS) override;
    void computeNmatrixAt(const FloatArray &iLocCoord, FloatMatrix &) override;
    bool computeGtoLRotationMatrix(FloatMatrix &answer) override;

    void computeConstitutiveMatrixAt(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) override;
    void computeStressVector(FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep) override;

    void computeBodyLoadVectorAt(FloatArray &answer, Load *load, TimeStep *tStep, ValueModeType mode) override;

    double giveKappaCoeff(TimeStep *tStep);
    double computeLength() override;
    double givePitch();
    void computeGaussPoints() override;
    MaterialMode giveMaterialMode() override { return _2dBeam; }
    int giveNumberOfIPForMassMtrxIntegration() override { return 4; }

    bool hasDofs2Condense() { return ( ghostNodes [ 0 ] || ghostNodes [ 1 ] ); }
};
} // end namespace oofem
#endif // beam2d_h
