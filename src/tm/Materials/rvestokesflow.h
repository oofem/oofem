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

#ifndef rvestokesflow_h
#define rvestokesflow_h

#ifdef __FM_MODULE

#include "floatarrayf.h"
#include "floatmatrixf.h"
#include "transportmaterial.h"
#include "stokesflowvelocityhomogenization.h"

#include <memory>

#define _IFT_RVEStokesFlow_Name "rvestokesflow"
#define _IFT_RVEStokesFlow_fileName "file"
#define _IFT_RVEStokesFlow_bctype "bctype"

namespace oofem {
/**
 * Material status class for the RVEStokesFlow class.
 *
 * @author Carl Sandström
 * @author Mikael Öhman
 */
class RVEStokesFlowMaterialStatus : public TransportMaterialStatus
{
protected:
    FloatMatrixF<3,3> temp_TangentMatrix, tangentMatrix;
    std :: unique_ptr< StokesFlowVelocityHomogenization > rve;

public:
    RVEStokesFlowMaterialStatus(int n, int rank, GaussPoint * g, const std :: string &inputfile);

    void setTimeStep(TimeStep *tStep);

    void initTempStatus() override;

    void updateYourself(TimeStep *tStep) override;

    void saveContext(DataStream &stream, ContextMode mode) override;
    void restoreContext(DataStream &stream, ContextMode mode) override;

    const FloatMatrixF<3,3> &giveTangentMatrix() const { return tangentMatrix; }
    const FloatMatrixF<3,3> &giveTempTangentMatrix() const { return temp_TangentMatrix; }
    void letTempTangentMatrixBe(const FloatMatrixF<3,3> &K) { temp_TangentMatrix = K; }

    StokesFlowVelocityHomogenization *giveRVE() { return rve.get(); }

    bool oldTangent;

    const char *giveClassName() const override { return "RVEStokesFlowMaterialStatus"; }
};


/**
 * Material class using an external .in file as a description of the substructure of a transport problem (in this case seepage). The external .in
 * file must be of the stokesflowvelocityhomogenization type.
 *
 * The response is computed by applying a pressure gradient over the domain and after solving, computing the mean velocity (which is the seepage
 * velocity).
 *
 * @author Carl Sandström
 * @author Mikael Öhman
 */
class RVEStokesFlow : public TransportMaterial
{
private:
    static int n;
    std :: string rveFilename;
    std :: string rveLogFilename;

public:
    RVEStokesFlow(int n, Domain * d);

    void initializeFrom(InputRecord &ir) override;

    FloatArrayF<3> computeFlux3D(const FloatArrayF<3> &grad, double field, GaussPoint *gp, TimeStep *tStep) const override;
    FloatMatrixF<3,3> computeTangent3D(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const override;

    double giveCharacteristicValue(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const override { return 0.0; }

    MaterialStatus *CreateStatus(GaussPoint *gp) const override;

    int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep) override;

    const char *giveClassName() const override { return "RVEStokesFlow"; }
    const char *giveInputRecordName() const override { return _IFT_RVEStokesFlow_Name; }
};
}

#endif // ifdef __FM_MODULE
#endif // rvestokesflow_h
