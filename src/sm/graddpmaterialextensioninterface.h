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

#ifndef graddpmaterialextensioninterface_h
#define graddpmaterialextensioninterface_h

#include "interface.h"
#include "matresponseform.h"
#include "matresponsemode.h"

namespace oofem {

class FloatMatrix;
class GaussPoint;
class TimeStep;

/**
 * Material interface for gradient material models.
 */
class GradDpMaterialExtensionInterface : public Interface
{
public:
    /// Left upper block
    virtual void givePDGradMatrix_uu(FloatMatrix &answer, MatResponseForm form, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) = 0;
    /// Left lower block
    virtual void givePDGradMatrix_ku(FloatMatrix &answer, MatResponseForm form, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) = 0;
    /// Right upper block
    virtual void givePDGradMatrix_uk(FloatMatrix &answer, MatResponseForm form, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) = 0;
    /// Right lower block
    virtual void givePDGradMatrix_kk(FloatMatrix &answer, MatResponseForm form, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) = 0;
    /// Stress-based averaging
    virtual void givePDGradMatrix_LD(FloatMatrix &answer, MatResponseForm form, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) = 0;
};

}
#endif