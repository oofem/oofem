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

#include "solutionstatusexportmodule.h"
#include "oofemcfg.h"
#include "engngm.h"
#include "timestep.h"
#include "classfactory.h"
#include "oofemtxtinputrecord.h"

#include <sstream>      // std::istringstream
#include <string>       // std::string
#include <cmath>

namespace oofem {
REGISTER_ExportModule(SolutionStatusExportModule)


SolutionStatusExportModule :: SolutionStatusExportModule(int n, EngngModel *e, FILE *out) : ExportModule(n, e), outputFile(out)
{
}

void
SolutionStatusExportModule :: initializeFrom(InputRecord &ir)
{
    ExportModule :: initializeFrom(ir);

    std::string formatStr("m:s:a:nite:t:dt:st:cr"); //default format
    IR_GIVE_OPTIONAL_FIELD(ir, formatStr, _IFT_SolutionStatusExportModule_format);
    // parse format string
    // use separator to read parts of the line
    std::istringstream input(formatStr);
    std::string token;
    while(std::getline(input, token, ':')) {
      recs.push_back(token);
    }
    this->checkRecs();
}

void
SolutionStatusExportModule :: doOutput(TimeStep *tStep, bool forcedOutput)
{
#if 0
    if ( !( testTimeStepOutput(tStep) || forcedOutput ) ) {
        return;
    }
#endif

    for (auto rec: this->recs) {
      if (rec == "m") {  // metstep number
        fprintf(outputFile, "%10d ", tStep->giveMetaStepNumber());
      } else if (rec == "s") { // solution step number
        fprintf(outputFile, "%10d ", tStep->giveNumber());
      } else if (rec == "a") { // attempt number
        fprintf(outputFile, "%10d ", tStep->numberOfAttempts);
      } else if (rec == "nite") { // number of iterations
        fprintf(outputFile, "%10d ", tStep->numberOfIterations); 
      } else if (rec == "t") {
        fprintf(outputFile, "%10.3e ", tStep->giveTargetTime());
      } else if (rec == "dt") {
        fprintf(outputFile, "%10.3e ", tStep->giveTimeIncrement());
      } else if (rec == "cr") {
        ConvergedReason r = tStep->convergedReason;
        std::string status;
        if (r == CR_CONVERGED) {
          status="Converged";
        } else if (r==CR_DIVERGED_ITS) {
          status="Diverged_I";
        } else if (r==CR_DIVERGED_TOL) {
          status="Diverged_T";
        } else if (r==CR_FAILED) {
          status="Failed";
        } else {
          status = "Unknown";
        }
        fprintf(outputFile, "%10s ", status.c_str());
      } else if (rec == "st") {
        int hours = static_cast<int>(tStep->solutionTime) / 3600;
        int minutes = (static_cast<int>(tStep->solutionTime) % 3600) / 60;
        double seconds = fmod(tStep->solutionTime, 60.0); // fractional seconds
        fprintf(outputFile, "%02d:%02d:%04.1f ", hours, minutes, seconds);
      } else if (rec == "-") {
        fprintf(outputFile, "%10s ", "-");
      }
    }
    fprintf(outputFile, "\n");
}

void
SolutionStatusExportModule::initialize()
{
  if (this->outputFile == nullptr) {
      char filename [100];
      sprintf( filename, "%s.m%d", this->emodel->giveOutputBaseFileName().c_str(), this->number);
      this->outputFile = fopen(filename, "w");
  }

  fprintf(outputFile,"%s (%s, %s)\nGitHash: %s\n", PRG_VERSION, HOST_TYPE, MODULE_LIST, OOFEM_GIT_HASH);
  time_t currtime = time(NULL);
  fprintf(outputFile,"Job: %s\n", this->emodel->giveReferenceFileName().c_str());
  fprintf(outputFile,"%s\n\n", ctime(& currtime));
  this->printRecsHeader();

}

void SolutionStatusExportModule::terminate()
{
  fprintf (outputFile, "\nAnalysis finished\n"); // would be nice to include completion status

  // get overall time consumed
  int rsec = 0, rmin = 0, rhrs = 0;
  int usec = 0, umin = 0, uhrs = 0;
  this->emodel->giveAnalysisTime(rhrs, rmin, rsec, uhrs, umin, usec);

  fprintf(outputFile, "Real time consumed: %03dh:%02dm:%02ds\n", rhrs, rmin, rsec);
  fprintf(outputFile, "User time consumed: %03dh:%02dm:%02ds\n\n", uhrs, umin, usec);
  
  int nw, ne;
  oofem_logger.getNumberOfWarningsAndErrors(nw, ne);
  fprintf (outputFile, "Total %d error(s) and %d warnings reported\n", ne, nw);

  fclose (this->outputFile);   
}


void
SolutionStatusExportModule :: checkRecs()
{
  std::string notrecognized;
  for (auto rec: this->recs) {
    if (!(rec == "m" || rec == "s" || rec =="a" || rec=="nite" || rec == "t" ||
          rec =="dt" || rec == "cr" || rec == "st" || rec=="-")) {
      notrecognized += " " + rec;
    }
  }
  if (notrecognized.size()) {
    OOFEM_WARNING ("SolutionStatusExportModule: invalid tokens detected: %s", notrecognized.c_str());
  }
}


void
SolutionStatusExportModule :: printRecsHeader()
{
  for (auto rec: this->recs) {
    std::string label;

    if (rec == "m") {  // metstep number
      label = "MetaStep";
    } else if (rec == "s") { // solution step number
      label = "Step";
    } else if (rec == "a") { // attempt number
      label = "Attempt";
    } else if (rec == "nite") { // number of iterations
      label = "Step.nite";
    } else if (rec == "t") {
      label = "Step.t";
    } else if (rec == "dt") {
      label = "Step.dt";
    } else if (rec == "cr") {
      label = "Status";
    } else if (rec == "st") {
      label = "Solver.t";
    } else if (rec == "-") {
      label = "-";
    }
    fprintf (outputFile, "%10s ", label.c_str());
  }
  fprintf(outputFile, "\n");
  std::string dashes (this->recs.size() * 11, '-');
  fprintf(outputFile, "%s\n", dashes.c_str());
}

} // end namespace oofem
