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
 *               Copyright (C) 1993 - 2019   Borek Patzak
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
//@}

namespace oofem {
/**
 * This class implements a 3-dimensional lattice element
 */

class Lattice3d : public LatticeStructuralElement
{
protected:
    double minLength;
    double kappa, length;
    double I1, I2, Ip;
    FloatArray polygonCoords;
    int numberOfPolygonVertices;
    FloatMatrix localCoordinateSystem;
    double eccS, eccT, area;
    FloatArray midPoint, centroid, globalCentroid, normal;
    int geometryFlag;
    int couplingFlag;
    IntArray couplingNumbers;
    FloatArray pressures;

public:
    Lattice3d(int n, Domain *);
    virtual ~Lattice3d();


    int giveLocalCoordinateSystem(FloatMatrix &answer) override;

    int computeGlobalCoordinates(FloatArray &answer, const FloatArray &lcoords) override;

    double giveLength() override;

    double giveNormalStress() override;

    double giveArea() override;

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

    /**
     * This function gives the cross-section coordinates.
     */
    void giveCrossSectionCoordinates(FloatArray &coords) override { coords = polygonCoords; }

    virtual void giveGPCoordinates(FloatArray &coords);

    virtual void computeGeometryProperties();

    virtual void computeCrossSectionProperties();

    const char *giveInputRecordName() const override { return _IFT_Lattice3d_Name; }
    const char *giveClassName() const override { return "Lattice3d"; }
    void initializeFrom(InputRecord &ir) override;



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
    void computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep) override;
    void computeConstitutiveMatrixAt(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) override;
    void computeStressVector(FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep) override;

    /**
     * This computes the geometrical properties of the element once.
     */
    void computePropertiesOfCrossSection();

    void computeGaussPoints() override;
    integrationDomain giveIntegrationDomain() const override { return _Line; }
};
} // end namespace oofem
#endif
