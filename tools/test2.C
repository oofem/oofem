#include "seek.h"

void main() {
    Tokenizer t(stdin);
    int i;
    double value;

    int res = seekSolutionStep(& t, 1.0);
    res = seekNodeRecord(& t, 2);
    res = seekDofRecord(& t, 2, 'd', value);
    printf("\n%lf ", value);
    res = seekElementRecord(& t, 3);
    res = seekGPRecord(& t, 1);
    res = seekStressStrainGPRecord(& t,  0, 1, 1, value);
    printf("%lf ", value);
    res = seekElementRecord(& t, 3);
    res = seekGPRecord(& t, 1);
    res = seekStressStrainGPRecord(& t,  0, 1, 2, value);
    printf("%lf ", value);
    res = seekElementRecord(& t, 3);
    res = seekGPRecord(& t, 1);
    res = seekStressStrainGPRecord(& t,  1, 0, 1, value);
    printf("%lf ", value);
    res = seekFirstReactionsRecord(& t);
    res = seekReactionRecord(& t, 2, 1, value);
    printf("%lf\n", value);
}

