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

#ifndef communicatormode_h
#define communicatormode_h

namespace oofem {
/**
 * The communicator mode determines the communication.
 */
enum CommunicatorMode {
    /**
     * The mode can be static, meaning that each node can assemble its communication maps
     * independently (or by independent communication). This implies that the size of
     * communication buffers is known in advance. Also if no data are planned to send to remote node, there
     * is no communication with this node (both sender and receiver know that there will be no data to send).
     */
    CommMode_Static,
    /**
     * (Dynamic) In this case the communication pattern and the amount of data sent between nodes is
     * not known in advance. This requires to use dynamic (packeted) buffering.
     */
    CommMode_Dynamic,
};
} // end namespace oofem
#endif // communicatormode_h
