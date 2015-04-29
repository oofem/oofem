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

#include "floatarray.h"
#include "floatmatrix.h"
#include "transportmaterial.h"
#include "stokesflowvelocityhomogenization.h"

#include <memory>

#define _IFT_RVEStokesFlow_Name "rvestokesflow"
#define _IFT_RVEStokesFlow_fileName "file"
#define _IFT_RVEStokesFlow_bctype "bctype"
#define _IFT_RVEStokesFlow_supressoutput "supressoutput"

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
    FloatMatrix temp_TangentMatrix, tangentMatrix;
    std :: unique_ptr< StokesFlowVelocityHomogenization > rve;

public:
    RVEStokesFlowMaterialStatus(int n, Domain * d, GaussPoint * g, const std :: string &inputfile);

    virtual ~RVEStokesFlowMaterialStatus();

    void setTimeStep(TimeStep *tStep);

    virtual void initTempStatus();

    virtual void updateYourself(TimeStep *tStep);

    virtual contextIOResultType saveContext(DataStream &stream, ContextMode mode, void *obj = NULL);
    virtual contextIOResultType restoreContext(DataStream &stream, ContextMode mode, void *obj = NULL);

    const FloatMatrix &giveTangentMatrix() { return tangentMatrix; }
    const FloatMatrix &giveTempTangentMatrix() { return temp_TangentMatrix; }
    void letTempTangentMatrixBe(const FloatMatrix &K) { temp_TangentMatrix = K; }

    StokesFlowVelocityHomogenization *giveRVE() { return rve.get(); }

    bool oldTangent;

    virtual const char *giveClassName() const { return "RVEStokesFlowMaterialStatus"; }
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

    int SupressRVEoutput;

    void suppressStdout();
    void enableStdout();

public:
    RVEStokesFlow(int n, Domain * d);

    virtual ~RVEStokesFlow() { }

    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual void giveFluxVector(FloatArray &answer, GaussPoint *gp, const FloatArray &grad, const FloatArray &field, TimeStep *tStep);
    virtual void giveCharacteristicMatrix(FloatMatrix &answer, MatResponseMode, GaussPoint *gp, TimeStep *tStep);
    virtual double giveCharacteristicValue(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) { return 0.0; }

    virtual MaterialStatus *CreateStatus(GaussPoint *gp) const;

    virtual int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep);

    virtual const char *giveClassName() const { return "RVEStokesFlow"; }
    virtual const char *giveInputRecordName() const { return _IFT_RVEStokesFlow_Name; }
};
}

#endif // rvestokesflow_h
