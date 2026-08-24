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

#ifndef lattice3dnl_h
#define lattice3dnl_h
#include "lattice3d.h"

///@name Input fields for lattice3dnl
//@{
#define _IFT_Lattice3dNL_Name "lattice3dnl"
#define _IFT_Lattice3dNL_mlength "mlength"
#define _IFT_Lattice3dNL_polycoords "polycoords"
#define _IFT_Lattice3dNL_couplingflag "couplingflag"
#define _IFT_Lattice3dNL_couplingnumber "couplingnumber"
#define _IFT_Lattice3dNL_pressures "pressures"
#define _IFT_Lattice3dNL_thickness "thickness"
//@}

namespace oofem {
/**
 * This class implements a geometric nonlinear 3-dimensional lattice element based on rigid body spring models.
Related reference (for the special case of a beam element with symmetric cross-section): A 3D frame element for large rotations based on the rigid-body-spring concept for analysing the failure of structures. International Journal of Solids and Structures. https://doi.org/10.1016/j.ijsolstr.2025.113812
 Author: Peter Grassl and Gumaa Abdelrhim, 2026
 */

class Lattice3dNL : public Lattice3d
{
protected:

public:
    Lattice3dNL(int n, Domain *);
    virtual ~Lattice3dNL();
  const char *giveInputRecordName() const override { return _IFT_Lattice3dNL_Name; }
    const char *giveClassName() const override { return "lattice3dnl"; }
  void giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord) override;

protected:
  bool computeGtoLRotationMatrix(FloatMatrix &) override;
    bool computeGtoLStrainRotationMatrix(FloatMatrix &answer);
    void computeNLBmatrixAt(GaussPoint *gp, FloatMatrix &, TimeStep *tStep);
    void updateYourself(TimeStep *tStep) override;
    void initForNewStep() override;
    void computeCurrentGtoLStrainRotationMatrix(FloatMatrix &GtoLCurrent, const FloatArray &u, const FloatArray &coordA, const FloatArray &coordB, const FloatArray &coordGP);
      void computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep) override;
    virtual void  computeStrainVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep) override;
};
} // end namespace oofem
#endif
