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

#ifndef simpletransportcrosssection_h
#define simpletransportcrosssection_h

#include "transportcrosssection.h"

#define _IFT_SimpleTransportCrossSection_Name "simpletransportcs"
#define _IFT_SimpleTransportCrossSection_material "mat"
#define _IFT_SimpleTransportCrossSection_thickness "thickness"
#define _IFT_SimpleTransportCrossSection_area "area"



namespace oofem {
class TransportMaterial;

/**
 * Transort cross-section. It's functionality is to be a wrapper around the material behavior.
 * @todo There will eventually be a layered version of this cross-section, so it must capture all the values.
 *
 * @author Mikael Öhman
 */
class OOFEM_EXPORT SimpleTransportCrossSection : public TransportCrossSection
{
protected:
    int matNumber = 0;
    double thickness = 0.;

public:
    /**
     * Constructor. Creates cross section with number n belonging to domain d.
     * @param n Cross section number.
     * @param d Domain for cross section.
     */
    SimpleTransportCrossSection(int n, Domain * d);

    void initializeFrom(InputRecord &ir) override;
    void giveInputRecord(DynamicInputRecord &input) override;

    /// Temporary function that hands out the material. Must be removed for future layered support, but input files will still look the same.
    TransportMaterial *giveMaterial() const override;
    Material *giveMaterial(IntegrationPoint *ip) const override;

    bool isCharacteristicMtrxSymmetric(MatResponseMode rMode) const override;
    int giveIPValue(FloatArray &answer, GaussPoint *ip, InternalStateType type, TimeStep *tStep) override;

    int checkConsistency() override;

    int packUnknowns(DataStream &buff, TimeStep *tStep, GaussPoint *gp) override;
    int unpackAndUpdateUnknowns(DataStream &buff, TimeStep *tStep, GaussPoint *gp) override;
    int estimatePackSize(DataStream &buff, GaussPoint *gp) override;

    const char *giveClassName() const override { return "SimpleTransportCrossSection"; }
    const char *giveInputRecordName() const override { return _IFT_SimpleTransportCrossSection_Name; }
};
} // end namespace oofem
#endif // simpletransportcrosssection_h
