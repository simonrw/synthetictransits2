#include "GenerateModel.h"
#include <cmath>
#include <iostream>
#include <vector>
#include <algorithm>
#include "Model.h"
#include "constants.h"

using namespace std;

namespace
{
    inline double calcOmega(const vector<double> &coeffs)
    {
        double returnval;
        for (int n=0; n<=4; ++n)
        {
            returnval += coeffs.at(n) / (n + 4.);
        }
        return returnval;
    }

    template <typename T>
    inline double square(T val) { return val * val; }

    double I(double r, double c1, double c2, double c3, double c4)
    {
        double rprimed = 1. - square(r);

        double I = 1.;
        I -= c1 * (1. - pow(rprimed, 1./4.));
        I -= c2 * (1. - pow(rprimed, 2./4.));
        I -= c3 * (1. - pow(rprimed, 3./4.));
        I -= c4 * (1. - pow(rprimed, 4./4.));
        return I;
    }

    double I(double r, const std::vector<double> &coeffs)
    {
        double rprimed = 1. - square(r);


        double I = 1.;
        for (int i=1; i<=4; ++i)
        {
            double coeff = coeffs.at(i);

            I -= coeff * (1. - pow(rprimed, static_cast<double>(i)/4.));
        }


        return I;
    }

    double IntegratedI(double dr, double c1, double c2, double c3, double c4, double rlow, double rhigh)
    {
        double sum = 0;
        for (double r=rlow; r<=rhigh; r+=dr)
        {
            sum += I(r, c1, c2, c3, c4) * dr * 2. * r;
        }
        return sum;

    }


    double IntegratedI(double dr, const std::vector<double> &coeffs, double rlow, double rhigh)
    {
        double sum = 0;
        for (double r=rlow; r<=rhigh; r+=dr)
        {
            sum += I(r, coeffs) * dr * 2. * r;
        }
        return sum;

    }
}

vector<double> GenerateSynthetic(const vector<double> &jd, const Model &m)
{
    /* All values are in normalised units so have to scale always */
    
    vector<double> Flux;
    double normalisedDistance = m.a * AU / (m.rs * rSun);
//    cout << "Normalisation constant: " << normalisedDistance << endl;

    const double dr = 0.001;

    /* Package the coefficients into a vector */
    vector<double> coeffs(5);
    coeffs[1] = m.c1;
    coeffs[2] = m.c2;
    coeffs[3] = m.c3;
    coeffs[4] = m.c4;
    coeffs[0] = 1. - coeffs[1] - coeffs[2] - coeffs[3] - coeffs[4];


    double omega = calcOmega(coeffs);
    double angFreq = 2. * M_PI / (m.period * secondsInDay);
//    cout << "Angular frequency: " << angFreq << " rad per sec" << endl;

    ///* get the cosine of the inclination */
    double cosi = cos(m.i * radiansInDegree);

    const double p = (m.rp * rJup) / (m.rs * rSun);

    for (unsigned int i=0; i<jd.size(); ++i)
    {
        /* Seconds since epoch */
        double t = (jd.at(i) - m.epoch) * secondsInDay;

        /* get the normalised device coordinates */
        double firstTerm = square(sin(angFreq * t));
        double secondTerm = square(cosi * cos(angFreq * t));




        double z = normalisedDistance * sqrt(firstTerm + secondTerm);

        double intpart;
        double phase = fabs(modf(t  / (m.period * secondsInDay) , &intpart));
        phase = phase > 0.5 ? phase - 1.0 : phase;

        double F = 0;

        /* Hack to make sure the secondary eclipse is not created */
        if ((phase > -0.25) && (phase < 0.25))
        {


            if (z <= 1 - p)
            {
                F = 0.;
                //F = 1. - square(p);
                double norm = 1. / (4. * z * p);
                double integral = IntegratedI(dr, coeffs, z-p, z+p);
                integral *= norm;
                F = 1. - (square(p) * integral / 4. / omega);
            }
            else if (z > 1 + p)
            {
                F = 1.;
            }
            else
            {
                double startPoint = z - p;
                double a = square(startPoint);
                double norm = 1./(1 - a);

                /* Integrate the I*(z) function from startPoint to 1 */
                double integral = IntegratedI(dr, coeffs, startPoint, 1.);
                integral *= norm;

                double insideSqrt = square(p) - square(z - 1.);
                double sqrtVal = sqrt(insideSqrt);
                sqrtVal *= (z - 1.);

                double insideAcos = (z - 1.) / p;
                double firstTerm = square(p) * acos(insideAcos);

                F = 1. - (integral * (firstTerm - sqrtVal) / (4. * M_PI * omega));
            }

        }
        else
        {
            F = 1.;
        }


        ///* append the data to the vectors */
        //lc.jd[i] = t / secondsInDay + midpoint;
        //lc.flux[i] = F;
        Flux.push_back(F);

    }
    ////outfile.close();

    //lc.radius = rPlan / rJup;
    //return lc;
    return Flux;
}

