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
 #include "nlinearstatic.h"
 #include "nlineardynamic.h"
 #include "eigenvaluedynamic.h"
 #include "deidynamic.h"
 #include "nldeidynamic.h"
 #include "diidynamic.h"
 #include "incrementallinearstatic.h"
 #include "linearstability.h"
 #include "linearstatic.h"
 #include "adaptnlinearstatic.h"
 #include "adaptlinearstatic.h"
#endif //__SM_MODULE

#ifdef __TM_MODULE
 #include "stationarytransportproblem.h"
 #include "nonstationarytransportproblem.h"
 #include "nltransienttransportproblem.h"
 #include "staggeredproblem.h"
 #include "darcyflow.h"
#endif //__TM_MODULE

#ifdef __FM_MODULE
 #include "cbs.h"
 #include "stokesflow.h"
 #include "stokesflowstresshomogenization.h"
 #include "stokesflowvelocityhomogenization.h"
 #include "supg.h"
#endif // __FM_Module

#ifdef __SM_MODULE
REGISTER_CLASS(LinearStatic, "linearstatic", LinearStaticClass)
REGISTER_CLASS(EigenValueDynamic, "eigenvaluedynamic", EigenValueDynamicClass)
REGISTER_CLASS(NonLinearStatic, "nonlinearstatic", NonLinearStaticClass)
REGISTER_CLASS(NonLinearDynamic, "nonlineardynamic", NonLinearDynamicClass)
REGISTER_CLASS(NlDEIDynamic, "nldeidynamic", NlDEIDynamicClass)
REGISTER_CLASS(DEIDynamic, "deidynamic", DEIDynamicClass)
REGISTER_CLASS(DIIDynamic, "diidynamic", DIIDynamicClass)
REGISTER_CLASS(IncrementalLinearStatic, "incrlinearstatic", IncrementalLinearStaticClass)
REGISTER_CLASS(LinearStability, "linearstability", LinearStabilityClass)
REGISTER_CLASS(AdaptiveNonLinearStatic, "adaptnlinearstatic", AdaptiveNonLinearStaticClass)
REGISTER_CLASS(AdaptiveLinearStatic, "adaptlinearstatic", AdaptiveLinearStaticClass)
#endif //__SM_MODULE

#ifdef __TM_MODULE
REGISTER_CLASS(StationaryTransportProblem, "stationaryproblem", StationaryTransportProblemClass)
REGISTER_CLASS(NonStationaryTransportProblem, "nonstationaryproblem", NonStationaryTransportProblemClass)
REGISTER_CLASS(NLTransientTransportProblem, "nltransienttransportproblem", NLTransientTransportProblemClass)
REGISTER_CLASS(StaggeredProblem, "staggeredproblem", StaggeredProblemClass)
REGISTER_CLASS(DarcyFlow, "darcyflow", DarcyFlowClass)
#endif //__TM_MODULE

#ifdef __FM_MODULE
REGISTER_CLASS(CBS, "cbs", CBSClass)
REGISTER_CLASS(SUPG, "supg", SUPGClass)
REGISTER_CLASS(StokesFlow, "stokesflow", StokesFlowClass)
REGISTER_CLASS(StokesFlowStressHomogenization, "stokesflowstresshomogenization", StokesFlowStressHomogenizationClass)
REGISTER_CLASS(StokesFlowVelocityHomogenization, "stokesflowvelocityhomogenization", StokesFlowVelocityHomogenizationClass)
#endif //__FM_MODULE

