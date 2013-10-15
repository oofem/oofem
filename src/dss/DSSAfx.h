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

/**
 * DSSAfx.h : include file for standard system include files,
 * or project specific include files that are used frequently, but
 * are changed infrequently
 *
 * @author Richard Vondracek
 */

#ifndef _DSSAFX_H__
#define _DSSAFX_H__

#if _MSC_VER > 1000
 #pragma once
#endif // _MSC_VER > 1000

#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <ctime>
#include <cstring>

#ifndef min
 #define min(a, b) ( ( a ) < ( b ) ? ( a ) : ( b ) )
#endif

#define DSS_NAMESPASE_BEGIN
#define DSS_NAMESPASE_END
//#define DSS_NAMESPASE_BEGIN namespace DSS {
//#define DSS_NAMESPASE_END }

#endif // _DSSAFX_H__

