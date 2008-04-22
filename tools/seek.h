#include "tokenizer.h"

#define SEEK_TOL  1.e-4

enum elemRec { er_undefined, er_strain, er_stress, er_status };
struct stateType {
    double currStep;
    double currEigval;
    double currLoaLevel;
    int currNodeSide, currElement, currDof, currGP, currReactionRecFound;
    int currQuasiReactionRecFound;
    int currLoaLevelRecFound;
    elemRec currElementRec;
};

void resetStateExceptCurrStep();
void getCurrState(stateType *cs);
int seekSolutionStep(Tokenizer *t, double stepVal);
int seekNodeRecord(Tokenizer *t, int node);
int seekDofRecord(Tokenizer *t, int dof, char u, double &val);
int readDofUnknown(Tokenizer *t, int dof, char u, double &val);
int seekElementRecord(Tokenizer *t, int elem);
int seekGPRecord(Tokenizer *t, int gp);
int seekStressStrainGPRecord(Tokenizer *t, int ifstress, int ifsrain,
                             int stressstraincomp, double &val);
int seekBeamRecord(Tokenizer *t, int ifforce, int ifdispl,
                   int stressstraincomp, double &val);
int seekFirstReactionsRecord(Tokenizer *t);
int seekReactionRecord(Tokenizer *t, int node, int dof, double &val);
int seekEigvalSolutionStep(Tokenizer *t, double stepVal, double &eigVal);
int seekNLSolutionStep(Tokenizer *t, double stepVal, double &loadlevel, double &nite);
int seekAndParseMaterialStatusRecord(Tokenizer *t, char *keyword, int valIndex, double &val);
int seekFirstQuasiReactionRecord(Tokenizer *t);
int seekQuasiReactionRecord(Tokenizer *t, int node, int dof, double &val);
int seekStepUserTime(Tokenizer *t, double &utime);
