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

#ifndef latticelink2d_h
#define latticelink2d_h

#include "latticestructuralelement.h"

///@name Input fields for LatticeLink2d
//@{
#define _IFT_LatticeLink2d_Name "latticelink2d"
#define _IFT_LatticeLink2d_length "length"
#define _IFT_LatticeLink2d_diameter "diameter"
#define _IFT_LatticeLink2d_dirvector "dirvector"
//@}

namespace oofem {
/**
 * This class implements a 2-dimensional link element, the 2D counterpart of LatticeLink3d.
 * It connects two nodes with u, v and the rotation about z (e.g. a frame node and a node
 * embedded in a continuum). The jump of node 2 relative to node 1 (with a rigid arm from
 * node 1 to node 2) is expressed in a local system with the first axis along dirvector:
 * slip along the bar, lateral jump across it and relative rotation. The material is
 * evaluated in _2dLattice mode, so the same lattice bond materials (e.g. latticeslip)
 * can be used in 2D and 3D.
 *
 * @author: Peter Grassl
 */

class LatticeLink2d : public LatticeStructuralElement
{
protected:
    double bondLength = 0.;

    FloatMatrix localCoordinateSystem;
    double bondDiameter = 0.;
    FloatArray directionVector;
    int geometryFlag = 0;
    FloatArray rigid;
    FloatArray globalCentroid;

public:
    LatticeLink2d(int n, Domain *);
    virtual ~LatticeLink2d();

    double computeVolumeAround(GaussPoint *aGaussPoint) override;

    double giveLength() override;

    double giveArea(GaussPoint *gp) override
    { return this->computeVolumeAround(gp) / this->giveLength(); }

    int giveLocalCoordinateSystem(FloatMatrix &answer) override;

    int computeGlobalCoordinates(Coordinates &answer, const FloatArray &lcoords) override;

    double giveBondLength() override;

    double giveBondDiameter() override;

    int computeNumberOfDofs() override { return 6; }

    void giveDofManDofIDMask(int inode, IntArray &) const override;

    void giveGpCoordinates(FloatArray &coords, GaussPoint *gp) override;

    virtual void computeGeometryProperties();

    void giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord = 0) override;

    const char *giveInputRecordName() const override { return _IFT_LatticeLink2d_Name; }
    const char *giveClassName()  const override { return "LatticeLink2d"; }
    void initializeFrom(const std::shared_ptr<InputRecord> &ir, int priority) override;
    void postInitialize() override;

    static ParamKey IPK_LatticeLink2d_length;
    static ParamKey IPK_LatticeLink2d_diameter;
    static ParamKey IPK_LatticeLink2d_dirvector;

    Element_Geometry_Type giveGeometryType() const override { return EGT_line_1; }

    void saveContext(DataStream &stream, ContextMode mode) override;

    void restoreContext(DataStream &stream, ContextMode mode) override;

#ifdef __OOFEG
    void drawYourself(oofegGraphicContext &context, TimeStep *tStep) override;
    void drawRawGeometry(oofegGraphicContext &, TimeStep *tStep) override;
    void drawDeformedGeometry(oofegGraphicContext &, TimeStep *tStep, UnknownType) override;
#endif


protected:
    void computeBmatrixAt(GaussPoint *, FloatMatrix &, int = 1, int = ALL_STRAINS) override;
    bool computeGtoLRotationMatrix(FloatMatrix &) override;
    void computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep) override;
    void computeConstitutiveMatrixAt(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) override;
    void computeStressVector(FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep) override;

    void computeGaussPoints() override;
    integrationDomain giveIntegrationDomain() const override { return _Line; }
};
} // end namespace oofem
#endif
