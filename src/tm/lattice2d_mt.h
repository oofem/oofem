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

#ifndef lattice2d_mt_h
#define lattice2d_mt_h

#include "latticetransportelement.h"
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
/**
 * This class implements a 2-dimensional lattice mass transport element
 */

class Lattice2d_mt : public LatticeTransportElement
{
protected:
    double area;
    double length;

    int couplingFlag;
    IntArray couplingNumbers;
    FloatArray crackWidths;
    FloatArray crackLengths;

    double dimension, width, thickness;
    FloatArray gpCoords;

    double crackWidth;

public:
    // constructor
    Lattice2d_mt(int, Domain *, ElementMode em = HeatTransferEM);
    virtual ~Lattice2d_mt();

    /** Computes the contribution to balance equation(s) due to internal sources */
    virtual void computeInternalSourceRhsVectorAt(FloatArray &answer, TimeStep *, ValueModeType mode);
    virtual double computeVolumeAround(GaussPoint *);

    virtual int giveCouplingFlag(){ return this->couplingFlag;}

    virtual void giveCouplingNumbers(IntArray &numbers){ numbers = this->couplingNumbers;}
    
    virtual void giveCrackLengths(FloatArray &lengths){ lengths = this->crackLengths;}

    virtual void giveCrackWidths(FloatArray &widths){ widths = crackWidths;}


    virtual int computeGlobalCoordinates(FloatArray &answer, const FloatArray &lcoords);

    virtual bool computeLocalCoordinates(FloatArray &answer, const FloatArray &gcoords);

    virtual void computeConductivityMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep);

    virtual void computeCapacityMatrix(FloatMatrix &answer, TimeStep *tStep);

    virtual const char *giveInputRecordName() const { return _IFT_Lattice2d_mt_Name; }
    virtual const char *giveClassName() const { return "Lattice2d_mtElement"; }

    virtual Element_Geometry_Type giveGeometryType() const { return EGT_line_1; }

    virtual double giveWidth() { return width; }
    virtual int computeNumberOfDofs() { return 2; }
    virtual void giveDofManDofIDMask(int inode, IntArray &) const;
    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual void updateInternalState(TimeStep *tStep);

#ifdef __OOFEG
    // Graphics output
    virtual void drawYourself(oofegGraphicContext &gc, TimeStep *tStep);
    virtual void drawRawGeometry(oofegGraphicContext &gc, TimeStep *tStep);
    virtual void drawRawCrossSections(oofegGraphicContext &gc, TimeStep *tStep);
    virtual void giveCrossSectionCoordinates(FloatArray &coords);
#endif

protected:
    virtual void computeGaussPoints();

    virtual void computeBmatrixAt(FloatMatrix &answer, const FloatArray &lcoords) { this->computeGradientMatrixAt(answer, lcoords); }
    virtual void  computeGradientMatrixAt(FloatMatrix &answer, const FloatArray &lcoords);
    virtual void  computeNmatrixAt(FloatMatrix &n, const FloatArray &);

    virtual double givePressure();

    virtual double giveOldPressure();

    virtual double giveMass();


    /* computes the submatrix of interpolation matrix cooresponding to single unknown.
     */
    virtual void  computeNSubMatrixAt(FloatMatrix &n, const FloatArray &);

    virtual double giveLength();

    virtual double giveArea() { return width * thickness; }


    virtual void  giveGpCoordinates(FloatArray &coords);

    virtual double computeEdgeVolumeAround(GaussPoint *gp, int iEdge) { return 0; }

    virtual int giveApproxOrder(int unknownIndx) { return 1; }
};
} // end namespace oofem
#endif
