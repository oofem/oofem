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

#ifndef contextioerr_h
#define contextioerr_h

#include "oofemcfg.h"
#include "contextioresulttype.h"

namespace oofem {
/**
 * Context IO exception class
 * @todo Document more.
 */
class OOFEM_EXPORT ContextIOERR
{
    contextIOResultType error;
    const char *msg, *file;
    int line;

public:

    ContextIOERR(contextIOResultType e, const char *file, int line);
    ContextIOERR(contextIOResultType e, const char *msg, const char *file, int line);
    ~ContextIOERR();

    void print();
};

#define THROW_CIOERR(e) throw ContextIOERR(e, __FILE__, __LINE__);
#define THROW_CIOERRM(e, m) throw ContextIOERR(e, m, __FILE__, __LINE__);
} // end namespace oofem
#endif // contextioerr_h
