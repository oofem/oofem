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



#pragma once
#include "contactsegment.h"
#include "node.h"
#include "inputrecord.h"
#include "Elements/structuralelement.h"
#include "set.h"
#include "functioncontactsegment.h"
#include "classfactory.h"

#define _IFT_CircleContactSegment_Name "circlecontactsegment"
 //#define _IFT_FunctionContactSegment_function "function"
#define _IFT_CircleContactSegment_centerpoint "centerpoint"
#define _IFT_CircleContactSegment_radius "radius"

namespace oofem {
    class CircleContactSegment : public FunctionContactSegment
    {
    public:
        CircleContactSegment(int n, Domain *aDomain) : FunctionContactSegment(n, aDomain) { ; }
        ~CircleContactSegment() {};

        IRResultType initializeFrom(InputRecord * ir) override;

        const char *giveClassName() const override { return "Circlecontactsegment"; }
        const char *giveInputRecordName() const override { return _IFT_CircleContactSegment_Name; }

    private:
        double radius;
        FloatArray centerPoint;

    protected:

        void computeContactPoint(FloatArray& answer, FloatArray& normal, const FloatArray& nodeCoords) override;

    };

}