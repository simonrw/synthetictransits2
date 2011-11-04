#ifndef MODEL_H

#define MODEL_H

#include <string>

struct Model
{
    int id;
    std::string name;
    int submodel_id;
    double period;
    double epoch;
    double a;
    double i;
    double rs;
    double rp;
    double mstar;
    double c1, c2, c3, c4;
    double teff;
};


#endif /* end of include guard: MODEL_H */
