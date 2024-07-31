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

#ifndef beam3d_h
#define beam3d_h

#include "sm/Elements/Beams/beambaseelement.h"
#include "sm/CrossSections/fiberedcs.h"
#include "sm/Materials/winklermodel.h"
#include "dofmanager.h"
#include "vtkxmlexportmodule.h"

///@name Input fields for Beam3d
//@{
#define _IFT_Beam3d_Name "beam3d"
#define _IFT_Beam3d_dofstocondense "dofstocondense"
#define _IFT_Beam3d_refnode "refnode"
#define _IFT_Beam3d_refangle "refangle"
#define _IFT_Beam3d_yaxis "yaxis"
#define _IFT_Beam3d_zaxis "zaxis"
#define _IFT_Beam3d_subsoilmat "subsoilmat"
//@}

#define Beam3d_nSubBeams 10

namespace oofem {

class FEI3dLineLin;

/**
 * This class implements a 2-dimensional beam element
 * with cubic lateral displacement interpolation (rotations are quadratic)
 * and longitudial displacements are linear.
 * This is an exact displacement approximation for beam with no
 * nonnodal loading.
 *
 * This class is not derived from liBeam3d or truss element, because it does not support
 * any material nonlinearities (if should, stiffness must be integrated)
 * 
 * @author Giovanni
 * @author Mikael Öhman
 * @author (several other authors)
 */
 class Beam3d : public BeamBaseElement, public FiberedCrossSectionInterface, public Beam3dSubsoilMaterialInterface,  public VTKXMLExportModuleElementInterface
{
protected:
    /// Geometry interpolator only.
    static FEI3dLineLin interp;

    double kappay, kappaz, length;
    int referenceNode;
    FloatArray yaxis, zaxis;
    double referenceAngle = 0;
    //IntArray *dofsToCondense;
    /*
     * Ghost nodes are used to introduce additional DOFs at element.
     * These are needed as we actually do not want to condense selected DOFs, but rather
     * allocate an extra equation to these. This allows to get cooresponding DOFs directly from
     * the global system, avoiding the need to postprocess local displacements at element.
     */
    DofManager *ghostNodes [ 2 ];
    /// number of condensed DOFs
    int numberOfCondensedDofs;

    /// Subsoil material
    int subsoilMat;

public:
    Beam3d(int n, Domain *d);
    virtual ~Beam3d();

    FEInterpolation *giveInterpolation() const override;

    void computeConsistentMassMatrix(FloatMatrix &answer, TimeStep *tStep, double &mass, const double *ipDensity = NULL) override;
    void computeInitialStressMatrix(FloatMatrix &answer, TimeStep *tStep) override;
    void computeLumpedInitialStressMatrix(FloatMatrix &answer, TimeStep *tStep) override;
    void computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep) override;
    int giveLocalCoordinateSystem(FloatMatrix &answer) override;
    void giveInternalForcesVector(FloatArray &answer, TimeStep *, int useUpdatedGpRecord = 0) override;
    void giveEndForcesVector(FloatArray &answer, TimeStep *tStep);

    int computeGlobalCoordinates(FloatArray &answer, const FloatArray &lcoords) override;


    int testElementExtension(ElementExtension ext) override {
        return ( ( ext == Element_EdgeLoadSupport ) ? 1 : 0 );
    }
    //int hasLayeredSupport () {return 1;}

    int computeNumberOfDofs() override { return 12; }
    int computeNumberOfGlobalDofs() override { return 12 + this->numberOfCondensedDofs; }
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

    int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep) override;

    void printOutputAt(FILE *file, TimeStep *tStep) override;

    //
    // fibered cross section support functions
    //
    void FiberedCrossSectionInterface_computeStrainVectorInFiber(FloatArray &answer, const FloatArray &masterGpStrain,
                                                                 GaussPoint *slaveGp, TimeStep *tStep) override;

    Interface *giveInterface(InterfaceType it) override;

    // definition & identification
    const char *giveClassName() const override { return "Beam3d"; }
    const char *giveInputRecordName() const override { return _IFT_Beam3d_Name; }
    void initializeFrom(InputRecord &ir) override;
    ///@todo Introduce interpolator and remove these two:
    integrationDomain giveIntegrationDomain() const override { return _Line; }
    //Element_Geometry_Type giveGeometryType() const override { return EGT_line_1; }
    Element_Geometry_Type giveGeometryType() const override { return EGT_Composite; }
    void updateLocalNumbering(EntityRenumberingFunctor &f) override;

    /*
    /// Subsoil support implemented directly enabling the postprocessing of end-forces

    virtual void B3SSI_getNMatrix (FloatMatrix &answer, GaussPoint *gp) {
      this->computeNmatrixAt(gp->giveNaturalCoordinates(), answer);
    }
    virtual void B3SSI_getGtoLRotationMatrix (FloatMatrix &answer) {
      this->computeGtoLRotationMatrix(answer);
    }
    virtual double B3SSI_computeVolumeAround (GaussPoint* gp) {
      return this->computeVolumeAround(gp);
    }
    */
    FloatMatrixF<6,6> B3SSMI_getUnknownsGtoLRotationMatrix() const override;

    void giveCompositeExportData(std::vector< ExportRegion > &vtkPieces, IntArray &primaryVarsToExport, IntArray &internalVarsToExport, IntArray cellVarsToExport, TimeStep *tStep ) override;

#ifdef __OOFEG
    void drawRawGeometry(oofegGraphicContext &gc, TimeStep *tStep) override;
    void drawDeformedGeometry(oofegGraphicContext & gc, TimeStep * tStep, UnknownType) override;
#endif

protected:
    void computeBoundaryEdgeLoadVector(FloatArray &answer, BoundaryLoad *load, int edge, CharType type, ValueModeType mode, TimeStep *tStep, bool global=true) override;
    int computeLoadGToLRotationMtrx(FloatMatrix &answer) override;
    void computeBmatrixAt(GaussPoint *, FloatMatrix &, int = 1, int = ALL_STRAINS) override;
    void computeNmatrixAt(const FloatArray &iLocCoord, FloatMatrix &) override;
    bool computeGtoLRotationMatrix(FloatMatrix &answer) override;
    void computeBodyLoadVectorAt(FloatArray &answer, Load *load, TimeStep *tStep, ValueModeType mode) override;

    double giveKappayCoeff(TimeStep *tStep);
    double giveKappazCoeff(TimeStep *tStep);
    void computeKappaCoeffs(TimeStep *tStep);
    double computeLength() override;
    void computeGaussPoints() override;
    void computeConstitutiveMatrixAt(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) override;
    void computeStressVector(FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep) override;

    MaterialMode giveMaterialMode() override { return _3dBeam; }
    int giveNumberOfIPForMassMtrxIntegration() override { return 4; }

    bool hasDofs2Condense() { return ( ghostNodes [ 0 ] || ghostNodes [ 1 ] ); }


    void computeSubSoilNMatrixAt(GaussPoint *gp, FloatMatrix &answer);
    void computeSubSoilStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep);

    void giveInternalForcesVectorAtPoint(FloatArray &answer, TimeStep *tStep, FloatArray &coords);
    void computeInternalForcesFromBoundaryEdgeLoadVectorAtPoint(FloatArray &answer, BoundaryLoad *load, int edge, CharType type, ValueModeType mode,
                                                                TimeStep *tStep, FloatArray &pointCoords, double ds, bool global);
    void computeInternalForcesFromBodyLoadVectorAtPoint(FloatArray &answer, Load *forLoad, TimeStep *tStep, ValueModeType mode, FloatArray &pointCoords, double ds);
    virtual int computeGlobalCoordinates(FloatArray &answer, const FloatArray &lcoords, const FloatArray &point);

};
} // end namespace oofem
#endif // beam3d_h
