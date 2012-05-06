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
 *               Copyright (C) 1993 - 2012   Borek Patzak
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

#ifdef __SM_MODULE
 #include "scalarerrorindicator.h"
 #include "zzerrorestimator.h"
 #include "combinedzzsiee.h"
 #include "huertaerrorestimator.h"
 #include "simpleinterfacemat.h"
#endif

#ifdef __SM_MODULE
REGISTER_CLASS(ScalarErrorIndicator, EET_SEI)
REGISTER_CLASS(ZZErrorEstimator, EET_ZZEE)
REGISTER_CLASS(CombinedZZSIErrorEstimator, EET_CZZSI)
REGISTER_CLASS(HuertaErrorEstimator, EET_HEE)
#endif

