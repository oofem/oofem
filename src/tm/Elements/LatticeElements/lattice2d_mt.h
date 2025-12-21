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
 *               Copyright (C) 1993 - 2025   Borek Patzak
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

#ifndef lattice2d_mt_h
#define lattice2d_mt_h

#include "tm/Elements/LatticeElements/latticetransportelement.h"
#include "spatiallocalizer.h"

///@name Input fields for Lattice2d_mt
//@{
#define _IFT_Lattice2d_mt_Name "latticemt2d"
#define _IFT_Lattice2DMT_dim "dim"
#define _IFT_Lattice2DMT_thick "thick"
#define _IFT_Lattice2DMT_width "width"
#define _IFT_Lattice2DMT_gpcoords "gpcoords"
#define _IFT_Lattice2DMT_crackwidth "crackwidth"
#define _IFT_Lattice2DMT_couplingflag "couplingflag"
#define _IFT_Lattice2DMT_couplingnumber "couplingnumber"
//@}


namespace oofem {
class ParamKey;
/**
 * This class implements a 2-dimensional lattice mass transport element
 */

class Lattice2d_mt : public LatticeTransportElement
{
protected:
    double area = -1.;
    double length = 0.;

    int couplingFlag = 0;
    IntArray couplingNumbers;
    FloatArray crackWidths;
    FloatArray crackLengths;

    double dimension = 0., width = 0., thickness = 0.;
    FloatArray gpCoords;

    double crackWidth = 0.;

    static ParamKey IPK_Lattice2d_mt_dim;
    static ParamKey IPK_Lattice2d_mt_thickness;
    static ParamKey IPK_Lattice2d_mt_width;
    static ParamKey IPK_Lattice2d_mt_gpcoords;
    static ParamKey IPK_Lattice2d_mt_crackwidth;
    static ParamKey IPK_Lattice2d_mt_couplingflag;
    static ParamKey IPK_Lattice2d_mt_couplingnumber;

public:
    Lattice2d_mt(int, Domain *, ElementMode em = HeatTransferEM);

    /** Computes the contribution to balance equation(s) due to internal sources */
    void computeInternalSourceRhsVectorAt(FloatArray &answer, TimeStep *, ValueModeType mode) override;
    double computeVolumeAround(GaussPoint *) override;

    int giveCouplingFlag() override { return this->couplingFlag; }

    void giveCouplingNumbers(IntArray &numbers) override { numbers = this->couplingNumbers; }

    void giveCrackLengths(FloatArray &lengths) override { lengths = this->crackLengths; }

    void giveCrackWidths(FloatArray &widths) override { widths = crackWidths; }


    int computeGlobalCoordinates(FloatArray &answer, const FloatArray &lcoords) override;

    bool computeLocalCoordinates(FloatArray &answer, const FloatArray &gcoords) override;

    void computeConductivityMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep) override;

    void computeCapacityMatrix(FloatMatrix &answer, TimeStep *tStep) override;

    const char *giveInputRecordName() const override { return _IFT_Lattice2d_mt_Name; }
    const char *giveClassName() const override { return "Lattice2d_mtElement"; }

    Element_Geometry_Type giveGeometryType() const override { return EGT_line_1; }

    double giveWidth() { return width; }
    int computeNumberOfDofs() override { return 2; }
    void giveDofManDofIDMask(int inode, IntArray &) const override;
    void initializeFrom(InputRecord &ir, int priority) override;
    void postInitialize() override;

    void updateInternalState(TimeStep *tStep) override;

#ifdef __OOFEG
    // Graphics output
    void drawYourself(oofegGraphicContext &gc, TimeStep *tStep) override;
    void drawRawGeometry(oofegGraphicContext &gc, TimeStep *tStep) override;
    void drawRawCrossSections(oofegGraphicContext &gc, TimeStep *tStep);
    void giveCrossSectionCoordinates(FloatArray &coords) override;
#endif

protected:
    void computeGaussPoints() override;

    void computeBmatrixAt(FloatMatrix &answer, const FloatArray &lcoords) override { this->computeGradientMatrixAt(answer, lcoords); }
    void computeGradientMatrixAt(FloatMatrix &answer, const FloatArray &lcoords) override;
    void computeNmatrixAt(FloatMatrix &n, const FloatArray &) override;

    double givePressure() override;

    double giveOldPressure() override;

    double giveMass() override;

    /* computes the submatrix of interpolation matrix cooresponding to single unknown.
     */
    void computeNSubMatrixAt(FloatMatrix &n, const FloatArray &);

    double giveLength() override;

    double giveArea() override { return width * thickness; }

    void  giveGpCoordinates(FloatArray &coords) override;

    double computeEdgeVolumeAround(GaussPoint *gp, int iEdge) override { return 0; }

    int giveApproxOrder(int unknownIndx) override { return 1; }
};
} // end namespace oofem
#endif
