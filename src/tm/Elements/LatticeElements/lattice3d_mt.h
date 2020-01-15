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

#ifndef lattice3d_mt_h
#define lattice3d_mt_h

#include "tm/Elements/LatticeElements/latticetransportelement.h"
#include "spatiallocalizer.h"

///@name Input fields for Lattice3d_mt
//@{
#define _IFT_Lattice3d_mt_Name "latticemt3d"
#define _IFT_Lattice3DMT_polycoords "polycoords"
#define _IFT_Lattice3DMT_crackwidths "crackwidths"
#define _IFT_Lattice3DMT_couplingflag "couplingflag"
#define _IFT_Lattice3DMT_couplingnumber "couplingnumber"
#define _IFT_Lattice3DMT_dim "dim"
#define _IFT_Lattice3DMT_area "area"
#define _IFT_Lattice3DMT_ranarea "ranarea"
#define _IFT_Lattice3DMT_mlength "mlength"
//@}


namespace oofem {
/**
 * This class implements a 3-dimensional lattice mass transport element
 */

class Lattice3d_mt : public LatticeTransportElement
{
protected:

    double minLength = 0.;
    double length = 0.;
    double I1 = 0., I2 = 0., Ip = 0.;
    FloatArray polygonCoords;
    int numberOfPolygonVertices;
    FloatMatrix localCoordinateSystem;
    double eccS = 0., eccT = 0., area = 0.;
    FloatArray midPoint, centroid, globalCentroid;
    int geometryFlag = 0;
    FloatArray normal;

    int couplingFlag = 0;
    IntArray couplingNumbers;
    FloatArray crackWidths;
    FloatArray crackLengths;

    double dimension = 0.;

public:
    Lattice3d_mt(int, Domain *, ElementMode em = HeatTransferEM);

    void computeInternalSourceRhsVectorAt(FloatArray &answer, TimeStep *, ValueModeType mode) override;
    double computeVolumeAround(GaussPoint *) override;

    int computeGlobalCoordinates(FloatArray &answer, const FloatArray &lcoords) override;

    bool computeLocalCoordinates(FloatArray &answer, const FloatArray &gcoords) override;

    void computeConductivityMatrix(FloatMatrix &answer, MatResponseMode rmode, TimeStep *tStep) override;

    void computeGeometryProperties();

    void computeCrossSectionProperties();

    void computeSpecialCrossSectionProperties();

    void computeCapacityMatrix(FloatMatrix &answer, TimeStep *tStep) override;

    const char *giveInputRecordName() const override { return _IFT_Lattice3d_mt_Name; }
    const char *giveClassName() const override { return "Lattice3d_mt"; }

    Element_Geometry_Type giveGeometryType() const override { return EGT_line_1; }

    int computeNumberOfDofs() override { return 2; }
    void giveDofManDofIDMask(int inode, IntArray &) const override;
    void initializeFrom(InputRecord &ir) override;
    void updateInternalState(TimeStep *tStep) override;

    double giveLength() override;
    double giveArea() override;
    int giveCouplingFlag() override { return this->couplingFlag; }
    void giveCouplingNumbers(IntArray &numbers) override { numbers = this->couplingNumbers; }
    void giveCrackWidths(FloatArray &widths) override { widths = crackWidths; }
    void giveCrackLengths(FloatArray &lengths) override;

    void giveCrossSectionCoordinates(FloatArray &coords) override { coords = polygonCoords; }
    int giveNumberOfCrossSectionNodes() override { return numberOfPolygonVertices; }

    void computeFlow(FloatArray &answer, GaussPoint *gp, TimeStep *tStep) override;

#ifdef __OOFEG

    void drawYourself(oofegGraphicContext &gc, TimeStep *tStep) override;

    void  drawRawGeometry(oofegGraphicContext &, TimeStep *tStep) override;

    void drawRawCrossSections(oofegGraphicContext &gc, TimeStep *tStep) override;

#endif

protected:
    void  computeGaussPoints() override;

    void  computeGradientMatrixAt(FloatMatrix &answer, const FloatArray &lcoords) override;

    void computeBmatrixAt(FloatMatrix &answer, const FloatArray &lcoords) override { this->computeGradientMatrixAt(answer, lcoords); }

    void  computeNmatrixAt(FloatMatrix &n, const FloatArray &) override;

    double computeEdgeVolumeAround(GaussPoint *gp, int iEdge) override { return 0; }

    void computeNSubMatrixAt(FloatMatrix &n, const FloatArray &);

    int giveApproxOrder(int unknownIndx) override { return 1; }
};
} // end namespace oofem
#endif
