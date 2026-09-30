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
 *               Copyright (C) 1993 - 2026   Borek Patzak
 *
 *
 *
 *       Czech Technical University, Faculty of Civil Engineering,
 *   Department of Structural Mechanics, 166 29 Prague, Czech Republic
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

#ifndef lattice3d_h
#define lattice3d_h

#include "latticestructuralelement.h"

///@name Input fields for Lattice3d
//@{
#define _IFT_Lattice3d_Name "lattice3d"
#define _IFT_Lattice3d_mlength "mlength"
#define _IFT_Lattice3d_polycoords "polycoords"
#define _IFT_Lattice3d_couplingflag "couplingflag"
#define _IFT_Lattice3d_couplingnumber "couplingnumber"
#define _IFT_Lattice3d_pressures "pressures"
#define _IFT_Lattice3d_shellnormal "shellnormal"
#define _IFT_Lattice3d_zaxis "zaxis"
#define _IFT_Lattice3d_s "s"
//@}

namespace oofem {
/**
 * This class implements a 3-dimensional lattice element
 */

class Lattice3d : public LatticeStructuralElement
{
protected:
    double minLength = 1.e-20;
    double kappa, length;
  double I1, I2, Ip, J;
     double shearArea1, shearArea2;
    FloatArray polygonCoords;
    int numberOfPolygonVertices;
    FloatMatrix localCoordinateSystem;
    double eccS, eccT, area;
    FloatArray midPoint, centroid, globalCentroid, normal;
    int geometryFlag;
    int couplingFlag = 0;
    IntArray couplingNumbers;
    FloatArray pressures;

    /// Shell surface normal in global coordinates; presence flags the element as a shell.
    FloatArray shellNormal;

    /// Reference axis orienting the transverse frame of a property-defined (shape 3) cross-section.
    FloatArray zaxis;

    /// Gauss-point position along the element: parameter (1+s)/2 from node A (s=0 gives midpoint).
    double s = 0.;

    double shellH = 0.;
    double shellB = 0.;

public:
    Lattice3d(int n, Domain *);
    virtual ~Lattice3d();


    int giveLocalCoordinateSystem(FloatMatrix &answer) override;

    int computeGlobalCoordinates(Coordinates &answer, const FloatArray &lcoords) override;

    double giveLength() override;

   bool giveRectangularSectionDimensions(double &by, double &bz, GaussPoint *gp) const override;

    double giveNormalStress() override;

    double giveArea(GaussPoint *gp) override;

    int computeNumberOfDofs() override { return 12; }

    void giveDofManDofIDMask(int inode, IntArray &) const override;

    double computeVolumeAround(GaussPoint *) override;

    int giveNumberOfCrossSectionNodes() override { return numberOfPolygonVertices; }

    int giveCrackFlag() override;

    double giveCrackWidth() override;

    void givePlasticStrain(FloatArray &plas) override;
    void giveOldPlasticStrain(FloatArray &plas) override;

    void givePressures(FloatArray &pres) override { pres = pressures; }

    int giveCouplingFlag() override { return couplingFlag; }

    void giveCouplingNumbers(IntArray &numbers) override { numbers = this->couplingNumbers; }

    void giveCrossSectionCoordinates(FloatArray &coords) override { coords = polygonCoords; }

    /// Per-layer local offsets (y, z) from centroid and tributary areas.
    virtual void computeLayerPositions(FloatArray &yOffset, FloatArray &zOffset, FloatArray &areas);

    /// True if shellnormal was supplied (element is tagged as a shell).
    bool isShellElement() const { return shellNormal.giveSize() == 3; }

    /// Shell normal in world frame (size 3 for shell elements, 0 otherwise).
    const FloatArray &giveShellNormal() const { return shellNormal; }

    /// Plate thickness (dimension along shell normal). Zero for non-shell elements.
    double giveShellH() const { return shellH; }

    /// True if this is a shell element with a LatticeCrossSection (layer IPs active).
    bool isHybridShell();

    virtual void computeGeometryProperties();

    virtual void computeCrossSectionProperties();

    const char *giveInputRecordName() const override { return _IFT_Lattice3d_Name; }
    const char *giveClassName() const override { return "Lattice3d"; }
    void initializeFrom(const std::shared_ptr<InputRecord> &ir, int priority) override;
    void postInitialize() override;

    static ParamKey IPK_Lattice3d_mlength;
    static ParamKey IPK_Lattice3d_polycoords;
    static ParamKey IPK_Lattice3d_couplingflag;
    static ParamKey IPK_Lattice3d_couplingnumber;
    static ParamKey IPK_Lattice3d_pressures;
    static ParamKey IPK_Lattice3d_shellnormal;
    static ParamKey IPK_Lattice3d_zaxis;
    static ParamKey IPK_Lattice3d_s;

    void giveGpCoordinates(FloatArray &coords, GaussPoint *gp) override;

  double giveJ(GaussPoint *gp) override;
  double giveI1(GaussPoint *gp) override;
  double giveI2(GaussPoint *gp) override;
    double giveShearArea1(GaussPoint *gp) override;
    double giveShearArea2(GaussPoint *gp) override;

    double giveTributaryWidth(GaussPoint *gp) override;

  void giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord) override;


    Element_Geometry_Type giveGeometryType() const override { return EGT_line_1; }

    void saveContext(DataStream &stream, ContextMode mode) override;

    void restoreContext(DataStream &stream, ContextMode mode) override;

#ifdef __OOFEG
    void drawYourself(oofegGraphicContext &context, TimeStep *tStep) override;
    void drawRawGeometry(oofegGraphicContext &, TimeStep *tStep) override;
    void drawRawCrossSections(oofegGraphicContext &, TimeStep *tStep);
    void drawDeformedGeometry(oofegGraphicContext &, TimeStep *tStep, UnknownType) override;
#endif

protected:
    void computeBmatrixAt(GaussPoint *, FloatMatrix &, int = 1, int = ALL_STRAINS) override;
    bool computeGtoLRotationMatrix(FloatMatrix &) override;
    void computeLumpedMassMatrix(FloatMatrix &answer, TimeStep *tStep) override;
    void computeMassMatrix(FloatMatrix &answer, TimeStep *tStep) override
    { this->computeLumpedMassMatrix(answer, tStep); }    
    void computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep) override;
    void computeConstitutiveMatrixAt(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) override;
    void computeStressVector(FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep) override;

    void computeGaussPoints() override;
    integrationDomain giveIntegrationDomain() const override { return _Line; }
};
} // end namespace oofem
#endif
