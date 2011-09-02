#ifndef GENERATEMODEL_H

#define GENERATEMODEL_H


Lightcurve GenerateSyntheticFromParams(const vector<double> &Time, double period, double midpoint, const vector<double> &coeffs, 
        double semi, double rPlan, double rStar, double inclination, double dr, double noise);


#endif /* end of include guard: GENERATEMODEL_H */
