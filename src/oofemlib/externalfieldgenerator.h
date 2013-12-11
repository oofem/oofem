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

#ifndef externalfieldgenerator_h
#define externalfieldgenerator_h

#include "randomfieldgenerator.h"
#include "gausspoint.h"

#define NAME_MAX_LENGTH 200

///@name Input fields for ExternalFieldGenerator
//@{
#define _IFT_ExternalFieldGenerator_Name "externalfieldgenerator"
#define _IFT_ExternalFieldGenerator_name "name"
//@}

namespace oofem {
/**
 * This class implements a randomfieldgenerator which reads
 * an externally generated field interpolates to determine
 * at Gausspoints.
 */
class OOFEM_EXPORT ExternalFieldGenerator : public RandomFieldGenerator
{
protected:
    Domain *domain;
    FloatArray field;
    IntArray numberReal;
public:
    /// Constructor. Creates empty RandomFieldGenerator
    ExternalFieldGenerator(int num, Domain *d);
    /// Destructor
    virtual ~ExternalFieldGenerator();

    virtual void generateRandomValue(double &value, FloatArray *position);

    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual const char *giveClassName() const { return "ExternalFieldGenerator"; }
    virtual const char *giveInputRecordName() const { return _IFT_ExternalFieldGenerator_Name; }
};
} // end namespace oofem

#endif
