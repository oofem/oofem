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

#ifndef pdelta_h
#define pdelta_h

#include "sm/EngineeringModels/linearstatic.h"
#include "geneigvalsolvertype.h"
#include "sparsegeneigenvalsystemnm.h"
#include "sparselinsystemnm.h"
#include "sparsemtrx.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "nummet.h"

///@name Input fields for LinearStability
//@{
#define _IFT_Pdelta_Name "pdelta"
#define _IFT_Pdelta_rtolv "rtolv"
#define _IFT_Pdelta_stype "stype"
#define _IFT_Pdelta_maxiter "maxiter"
#define _IFT_Pdelta_lumpedInitialStressMatrix "lumped"
//@}

namespace oofem {

/**
 * This class implements p-delta analysis, where the effects of normal force on deformed configuration
 * is taken into account by means of initial stress matrix.
 *
 * Solution of this problem is base on equation in the form of: @f$ (K+K_\sigma(r))\cdot r=f @f$.
 * This is a nonlinear problem solved using simple iteration method.
 */
class Pdelta : public LinearStatic
{
private:
    std :: unique_ptr< SparseMtrx > initialStressMatrix;
    double rtolv;
    int maxiter;
    bool lumpedInitialStressMatrix; 

public:
    Pdelta(int i, EngngModel *master=nullptr);
    virtual ~Pdelta() { }

    void solveYourselfAt(TimeStep *tStep) override;
    void initializeFrom(InputRecord &ir) override;

    // identification
    const char *giveInputRecordName() const override { return _IFT_Pdelta_Name; }
    const char *giveClassName() const override { return "Pdelta"; }
};
} // end namespace oofem
#endif // pdelta_h
