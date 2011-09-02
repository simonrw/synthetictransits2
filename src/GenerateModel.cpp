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
}

vector<double> GenerateSynthetic(const vector<double> &jd, const Model &m)
{
    /* All values are in normalised units so have to scale always */
    
    double normalisedDistance = m.a * AU / (m.rs * rSun);
    cout << "Normalisation constant: " << normalisedDistance << endl;

    /* Package the coefficients into a vector */
    vector<double> coeffs(5);
    coeffs[1] = m.c1;
    coeffs[2] = m.c2;
    coeffs[3] = m.c3;
    coeffs[4] = m.c4;
    coeffs[0] = 1. - coeffs[1] - coeffs[2] - coeffs[3] - coeffs[4];


    double omega = calcOmega(coeffs);
    double angFreq = 2. * M_PI / (m.period * secondsInDay);
    cout << "Angular frequency: " << angFreq << " rad per sec" << endl;

    //[> get the cosine of the inclination <]
    //double cosi = cos(inclination);

    //const double p = rPlan / rStar;

    //[> set up the random number generator <]
    //Normaldev_BM randGenerator(0., 1., time(NULL));



    //[> Output lightcurve <]
    //Lightcurve lc(Time.size());
    //lc.period = period;
    //lc.epoch = midpoint;

    //[> parallel process this part <]
//#pragma omp parallel for
    //for (unsigned int i=0; i<Time.size(); ++i)
    //{
        //double t = Time.at(i);

        //[> get the normalised device coordinates <]
        //double firstTerm = square(sin(angFreq * t));
        //double secondTerm = square(cosi * cos(angFreq * t));




        //double z = normalisedDistance * sqrt(firstTerm + secondTerm);

        //double intpart;
        //double phase = fabs(modf(t  / period , &intpart));
        //phase = phase > 0.5 ? phase - 1.0 : phase;

        //double F = 0;

        //[> Hack to make sure the secondary eclipse is not created <]
        //if ((phase > -0.25) && (phase < 0.25))
        //{


            //if (z <= 1 - p)
            //{
                //F = 0.;
                ////F = 1. - square(p);
                //double norm = 1. / (4. * z * p);
                //double integral = IntegratedI(dr, coeffs, z-p, z+p);
                //integral *= norm;
                //F = 1. - (square(p) * integral / 4. / omega);
            //}
            //else if (z > 1 + p)
            //{
                //F = 1.;
            //}
            //else
            //{
                //double startPoint = z - p;
                //double a = square(startPoint);
                //double norm = 1./(1 - a);

                //[> Integrate the I*(z) function from startPoint to 1 <]
                //double integral = IntegratedI(dr, coeffs, startPoint, 1.);
                //integral *= norm;

                //double insideSqrt = square(p) - square(z - 1.);
                //double sqrtVal = sqrt(insideSqrt);
                //sqrtVal *= (z - 1.);

                //double insideAcos = (z - 1.) / p;
                //double firstTerm = square(p) * acos(insideAcos);

                //F = 1. - (integral * (firstTerm - sqrtVal) / (4. * M_PI * omega));
            //}

        //}
        //else
        //{
            //F = 1.;
        //}

        //[> add the noise <]
        //F += noise * randGenerator.dev();

        //[> append the data to the vectors <]
        //lc.jd[i] = t / secondsInDay + midpoint;
        //lc.flux[i] = F;

    //}
    ////outfile.close();

    //lc.radius = rPlan / rJup;
    //return lc;
}

