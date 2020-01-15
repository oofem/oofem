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

#ifndef lattice2d_h
#define lattice2d_h

#include "sm/Elements/LatticeElements/latticestructuralelement.h"

///@name Input fields for Lattice2d
//@{
#define _IFT_Lattice2d_Name "lattice2d"
#define _IFT_Lattice2d_thick "thick"
#define _IFT_Lattice2d_width "width"
#define _IFT_Lattice2d_gpcoords "gpcoords"
#define _IFT_Lattice2d_couplingflag "couplingflag"
#define _IFT_Lattice2d_couplingnumber "couplingnumber"
//@}

namespace oofem {
/**
 * This class implements a 2-dimensional lattice element
 */
class Lattice2d : public LatticeStructuralElement
{
protected:
    double kappa, pitch, length;

    double width, thickness;
    FloatArray gpCoords;
    int couplingFlag;
    IntArray couplingNumbers;

public:
    Lattice2d(int n, Domain *d);
    virtual ~Lattice2d();

    int giveLocalCoordinateSystem(FloatMatrix &answer) override;

    /**
     * This function is different from the standard computeGlobalCorrdinates
     * function as it returns the global coordinates of the gausspoint
     * independent to the value of the lcoords.
     */
    int computeGlobalCoordinates(FloatArray &answer, const FloatArray &lcoords) override;

    double giveLength() override;

    double giveNormalStress() override;

    int hasBeenUpdated() override;

    double giveArea() override { return this->width * this->thickness; }

    int computeNumberOfDofs() override { return 6; }
    void giveDofManDofIDMask(int inode, IntArray &) const override;
    double computeVolumeAround(GaussPoint *gp) override;

    int giveCrackFlag() override;

    double giveCrackWidth() override;
    //double giveOldCrackWidth() override;

    double giveDissipation() override;
    double giveDeltaDissipation() override;

    int giveCouplingFlag() override { return couplingFlag; }

    void giveCouplingNumbers(IntArray &numbers) override { numbers = this->couplingNumbers; }
    //
    // definition & identification
    //
    const char *giveInputRecordName() const override { return _IFT_Lattice2d_Name; }
    const char *giveClassName() const override { return "Lattice2d"; }
    void initializeFrom(InputRecord &ir) override;
    Element_Geometry_Type giveGeometryType() const override { return EGT_line_1; }

#ifdef __OOFEG
    void drawYourself(oofegGraphicContext &gc, TimeStep *tStep) override;
    void drawRawGeometry(oofegGraphicContext &gc, TimeStep *tStep) override;
    void drawDeformedGeometry(oofegGraphicContext & gc, TimeStep *tStep, UnknownType) override;
    void drawSpecial(oofegGraphicContext &gc, TimeStep *tStep) override;
    void drawRawCrossSections(oofegGraphicContext &gc, TimeStep *tStep);
    void giveCrossSectionCoordinates(FloatArray &coords);
#endif

protected:
    void saveContext(DataStream &stream, ContextMode mode) override;
    void restoreContext(DataStream &stream, ContextMode mode) override;

    void computeBmatrixAt(GaussPoint *, FloatMatrix &, int = 1, int = ALL_STRAINS) override;
    bool computeGtoLRotationMatrix(FloatMatrix &) override;
    void computeStressVector(FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep) override;
    void computeConstitutiveMatrixAt(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) override;
    void computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep) override;

    int giveNumberOfCrossSectionNodes() override { return 2; }
    double givePitch();
    void computeGaussPoints() override;
    integrationDomain giveIntegrationDomain() const override { return _Line; }
    void giveGpCoordinates(FloatArray &coords) override;
};
} // end namespace oofem
#endif
