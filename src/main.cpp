#include <iostream>
#include <cstring>
#include <cassert>
#include <stdexcept>
#include <sqlitepp/sqlitepp.hpp>
#include <fitsio.h>
#include <tclap/CmdLine.h>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <cstdlib>
#include <fstream>

/* Local includes */
#include "timer.h"
#include "GenerateModel.h"
#include "Model.h"
#include "FitsObject.h"
#include "constants.h"
#include "ObjectSkipDefs.h"



using namespace std;
using namespace sqlitepp;

typedef vector<string> stringlist;

enum
{
    add=101,
    sub
};

struct ConfigContainer
{
    bool isWASP;
    string DatabaseFilename;
    string SourceFilename;
    string OutputFilename;
};

template <typename T>
T square(T val) { return val * val; }

const double jd_ref = 2453005.5;

double wd2jd(double wd)
{
    return (wd / secondsInDay) + jd_ref;
}

double jd2wd(double jd)
{
    return (jd - jd_ref) * secondsInDay;
}


class ArithMeth
{
    
    int m_type;
    
public:
    ArithMeth(const string &type)
    {
        if (type == "+")
        {
            this->m_type = add;
        }
        else if (type == "-")
        {
            this->m_type = sub;
        }
        else
        {
            throw runtime_error("Unknown method specified (must be + or -)");
        }
    }
    
    int type() const { return this->m_type; }
};

double WidthFromParams(const Model &m)
{
    /* Returns the width of the full transit based on some lc parameters
     *
     * \frac{P}{\pi} \asin{\sqrt{(\frac{R_P + R_S}{a})^2 - \cos^2 i}}
     */
    const double Norm = (m.period * secondsInDay) / M_PI;
    if (m.a == 0)
    {
        /* Something hasn't updated the separation */
        return 0;
    }
    
    const double FirstTerm = square(((m.rp * rJup) + (m.rs * rSun)) / (m.a * AU));
    const double SecondTerm = square(cos(m.i * radiansInDegree));
    const double InsideSqrt = FirstTerm - SecondTerm;
    
    /* Check that InsideSqrt is not <= 0 */
    if (InsideSqrt <= 0.0)
    {
        return 0.0;
    }
    
    const double SqrtInsideSqrt = sqrt(InsideSqrt);
    
    if ((SqrtInsideSqrt < -1.0) || (SqrtInsideSqrt > 1.0))
    {
        return 0;
    }
    
    
    
    return Norm * asin(SqrtInsideSqrt);
    
}

template <typename T>
T WeightedMedian(const vector<T> &data, const double siglevel)
{
    vector<T> buffer;
    const size_t N = data.size();
    
    /* First calculate the mean */
    double av = 0;
    int ValidPoints = 0;
    for (size_t i=0; i<N; ++i)
    {
        if (!isnan(data.at(i)))
        {
            av += data.at(i);
            ValidPoints++;
        }
    }
    av /= (double)ValidPoints;
    
    /* Put in a check for if the av is 0 */
    if (av == 0)
    {
        throw runtime_error("Cannot calculate average - all points are 0");
    }

    if (isnan(av))
    {
        throw runtime_error("Cannot calculate average - average is NaN");
    }

    /* now the sigma */
    double sd = 0;
    for (size_t i=0; i<N; ++i)
    {
        if (!isnan(data.at(i)))
        {
            sd += (data.at(i) - av)*(data.at(i) - av);
        }
    }
    sd /= (double)ValidPoints;
    sd = sqrt(sd);

    const double upperlim = av + siglevel * sd;
    const double lowerlim = av - siglevel * sd;

    for (int i=0; i<N; ++i)
    {
        if ((data.at(i) < upperlim) && (data.at(i) > lowerlim) && !isnan(data.at(i)))
        {
            buffer.push_back(data.at(i));
        }
    }

//    cout << N - buffer.size() << " elements rejected" << endl;

    /* Sort the array */
    sort(buffer.begin(), buffer.end());

    /* Return the middle one */
    return buffer.at(buffer.size() / 2);


}



struct FalseColumnNumbers
{
    int skipdet, period, width, depth, epoch,
    rp, rs, a, i;
};


long indexOf(const stringlist &stringlist, const string &comp)
{
    for (size_t i=0; i<stringlist.size(); ++i)
    {
        if (stringlist.at(i) == comp)
        {
            return i;
        }
        else
        {
            string nameWithoutWhitespace(comp);
            nameWithoutWhitespace.erase(remove_if(nameWithoutWhitespace.begin(), nameWithoutWhitespace.end(), ::isspace), nameWithoutWhitespace.end());

            if (stringlist.at(i) == nameWithoutWhitespace)
            {
                return i;
            }
        }
    }


    /* If the loop gets here the object is not found */
    throw runtime_error("Cannot find object");
}

void AlterLightcurveData(Fits &f, const int startindex, const int length, const Model &m, const ArithMeth &arithtype, const ConfigContainer &Config)
{
    f.moveHDU("HJD");

    vector<double> jd(length);
    
    /* Fetch the jd data */
    fits_read_img(*f.fptr(), TDOUBLE, startindex, length, 0, &jd[0], 0, &f.status());

    /* if the Config.isWASP parameter is set then convert the array to jd */
    if (Config.isWASP)
    {
        for (int i=0; i<length; ++i)
        {
            jd[i] = wd2jd(jd[i]);
        }
    }
    
    /* Now get the addition model */
    vector<double> ModelFlux = GenerateSynthetic(jd, m);
    
    f.moveHDU("FLUX");
    vector<double> OriginalFlux(length);
    fits_read_img(*f.fptr(), TDOUBLE, startindex, length, 0, &OriginalFlux[0], 0, &f.status());
    
    /* Normalise to the weighted median flux level */
    double WeightedMed = WeightedMedian(OriginalFlux, 2.5);
    
    vector<double> TransitAdded(length);
    double result = 0;
    for (int i=0; i<length; ++i)
    {
        double fluxval = OriginalFlux.at(i);
        fluxval /= WeightedMed;
        
        double modelval = ModelFlux.at(i);

        /* Perform the arithmetic */
        
        if (arithtype.type() == add)
        {
            result = WeightedMed * (fluxval + modelval - 1.0);
        }
        else if (arithtype.type() == sub)
        {
            result = WeightedMed * (fluxval - modelval + 1.0);
        }
        TransitAdded[i] = result;
    }
    
    
    
    
    fits_write_img(*f.fptr(), TDOUBLE, startindex, length, &TransitAdded[0], &f.status());
    f.check();


}

stringlist &split(const string &s, char delim, stringlist &elems)
{
    stringstream ss(s);
    string item;
    while(getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}


stringlist split(const string &s, char delim)
{
    stringlist elems;
    return split(s, delim, elems);
}

string AlterObjectName(const string &OriginalName)
{

    /* Split the string at the J character */
    stringlist parts = split(OriginalName, 'J');

    /* check the length for validity */
    if (parts.size() != 2)
    {
        throw runtime_error("Unknown object name encountered");
    }

    stringstream NewName;
    NewName << "1SYNTH J" << parts[1];



    return NewName.str();
}

template <typename T>
void OverPrint(const T &val)
{
    cout << "\r" << val;
    cout.flush();

}

int main(int argc, char *argv[])
{
    try
    {
        TCLAP::CmdLine cmd("");
        TCLAP::ValueArg<string> infile_arg("i", "infile", "Input fits file", true, "", "Fits file", cmd);
        TCLAP::ValueArg<string> candidates_arg("c", "candidates", "Input candidates file", true, "", "SQLite3 database", cmd);
        TCLAP::ValueArg<string> output_arg("o", "output", "Output file", false, "output.fits", "Fits file", cmd);
        cmd.parse(argc, argv);

        Timer ts;
        ts.start("all");
        ts.start("copy");

        /* Set up the config variable */
        ConfigContainer Config;

        Fits infile(infile_arg.getValue());

        /* Make sure the project has the 'project' header key,
         * and if it's wasp make sure all jds are converted 
         * to wasp dates 
         */
        infile.moveHDU(1);

        char Project_cstr[FLEN_VALUE];
        fits_read_key(*infile.fptr(), TSTRING, "PROJECT", &Project_cstr, 0, &infile.status());

        /* check the key exists */
        if (infile.status() == KEY_NO_EXIST)
        {
            throw runtime_error("Project key not found, please set this");
        }
        
        /* Create string object for easy comparisons */
        string Project(Project_cstr);

        cout << Project << " project found" << endl;

        Config.isWASP = (Project == "WASP") ? true : false;
        Config.DatabaseFilename = candidates_arg.getValue();
        Config.SourceFilename = infile_arg.getValue();
        Config.OutputFilename = output_arg.getValue();

        if (Config.isWASP)
        {
            cout << "--- Converting times to WASP data" << endl;
        }


        NewFits outfile("!" + output_arg.getValue());

        /* Add the transinj key */
        outfile.moveHDU(1);
        bool transinj_val = true;
        fits_write_key(*outfile.fptr(), TLOGICAL, "TRANSINJ", &transinj_val, "Contains false transits", &outfile.status());
        outfile.check();

        /* Start by getting file information from the input */
        int nhdus = 0;
        fits_get_num_hdus(*infile.fptr(), &nhdus, &infile.status());
        infile.check();

        cout << nhdus << " hdus found" << endl;


        /* Open the sqlite3 database here */
        session conn(candidates_arg.getValue());
        statement st(conn);

        /* get the required number of new objects */
        int nextra = 0;
        st << "select count(*) from addmodels", into(nextra);
        if (!st.exec())
        {
            throw runtime_error("No input models found");
        }

        cout << "Inserting " << nextra << " extra models" << endl;


        for (int hdu=2; hdu<=nhdus; ++hdu)
        {
            int status = 0;
            infile.moveHDU(hdu);

            cout << "HDU: " << infile.hduname();

            /* Copy the data and header across */
            fits_copy_hdu(*infile.fptr(), *outfile.fptr(), 0, &status);

            /* Get the hdu type */
            int hdutype = 0;
            fits_get_hdu_type(*outfile.fptr(), &hdutype, &outfile.status());
            outfile.check();

            /* If an image hdu is found then just resize it */
            if (hdutype == IMAGE_HDU)
            {
                /* Get the current dimensions */
                long naxes[2];
                fits_get_img_size(*outfile.fptr(), 2, naxes, &outfile.status());
                outfile.check();

                cout << " - " << naxes[0] << "x" << naxes[1] << " pix";

                int bitpix = 0;
                fits_get_img_type(*outfile.fptr(), &bitpix, &outfile.status());
                outfile.check();

                long newnaxes[] = {naxes[0], naxes[1] + nextra};

                fits_resize_img(*outfile.fptr(), bitpix, 2, newnaxes, &outfile.status());
                outfile.check();

                cout << " -> " << newnaxes[0] << "x" << newnaxes[1] << " pix";
            }
            else
            {
                /* If the catalogue extension is found then add extra rows */
                const string hduname = outfile.hduname();
                long nrows = outfile.nrows();

                cout << " - " << nrows << " rows";

                if (hduname == "CATALOGUE")
                {
                    /* Need to add extra columns */
                    int ncols = 0;
                    fits_get_num_cols(*outfile.fptr(), &ncols, &outfile.status());
                    outfile.check();

                    cout << ", " << ncols << " columms";

                    /* Append extra rows */

                    fits_insert_rows(*outfile.fptr(), nrows, nextra, &outfile.status());
                    outfile.check();

                    cout << " -> " << nrows + nextra << " rows";

                    /* Need 9 extra columns */

                    char *ColumnNames[] = {"SKIPDET", "FAKE_PERIOD", "FAKE_WIDTH", "FAKE_DEPTH", "FAKE_EPOCH", "FAKE_RP", "FAKE_RS", "FAKE_A", "FAKE_I"};
                    char *ColumnFormats[] = {"1I", "1D", "1D", "1D", "1D", "1D", "1D", "1D", "1D"};

                    size_t nNewCols = sizeof(ColumnNames) / sizeof(char*);
                    assert((sizeof(ColumnFormats) / sizeof(char*)) == nNewCols);

                    fits_insert_cols(*outfile.fptr(), ncols + 1, nNewCols, ColumnNames, ColumnFormats, &outfile.status());
                    outfile.check();

                    fits_get_num_cols(*outfile.fptr(), &ncols, &outfile.status());
                    outfile.check();

                    cout << ", " << ncols << " columms";


                }
            }




            Fits::check(status);
            cout << endl;

        }


        ts.stop("copy");

        /* File copy finished */
        
        /* Prefetch the column numbers for later use */
        outfile.moveHDU("CATALOGUE");
        FalseColumnNumbers fcn;
        fcn.skipdet = outfile.columnNumber("SKIPDET");
        fcn.period = outfile.columnNumber("FAKE_PERIOD");
        fcn.width = outfile.columnNumber("FAKE_WIDTH");
        fcn.depth = outfile.columnNumber("FAKE_DEPTH");
        fcn.epoch = outfile.columnNumber("FAKE_EPOCH");
        fcn.rp = outfile.columnNumber("FAKE_RP");
        fcn.rs = outfile.columnNumber("FAKE_RS");
        fcn.a = outfile.columnNumber("FAKE_A");
        fcn.i = outfile.columnNumber("FAKE_I");
        outfile.check();




        ts.start("model.iterate");
        Model Current;
        const int NullSubIndex = -1;


        /* Get a list of the objects in the file */
        stringlist ObjectNames;
        infile.moveHDU("CATALOGUE");
        int obj_id_colno = infile.columnNumber("OBJ_ID");

        /* Read the data in as strings */
        int dispwidth;
        fits_get_col_display_width(*infile.fptr(), obj_id_colno, &dispwidth, &infile.status());

        long nrows;
        fits_get_num_rows(*infile.fptr(), &nrows, &infile.status());

        vector<char*> cstrnames(nrows);
        for (int i=0; i<nrows; ++i) cstrnames[i] = new char[dispwidth+1];
        fits_read_col_str(*infile.fptr(), obj_id_colno, 1, 1, nrows, 0, &cstrnames[0], 0, &infile.status());

        for (int i=0; i<nrows; ++i)
        {
            string CurrentName = cstrnames[i];
            /* Remove all whitespace */
            CurrentName.erase(remove_if(CurrentName.begin(), CurrentName.end(), ::isspace), CurrentName.end());
            ObjectNames.push_back(CurrentName);
            delete[] cstrnames[i];
        }

        /* Now iterate through every row adding a new lightcurve, and subtracting if necassary  */
        st << "select id, name, submodel_id, period, epoch, a, i, rs, rp, mstar, c1, c2, c3, c4, teff "
        " from addmodels", into(Current.id), into(Current.name), into(Current.submodel_id), 
        into(Current.period), into(Current.epoch), into(Current.a), into(Current.i), into(Current.rs),
        into(Current.rp), into(Current.mstar), into(Current.c1), into(Current.c2), into(Current.c3), 
        into(Current.c4), into(Current.teff);

        long counter = 0;
        cout << "Generating models" << endl;
        while (st.exec())
        {
            /* Location to write the data to */
            const long OutputIndex = nrows + counter;
            
            /* Need to append 1 for the catalogue information as the catalogue
             is 1 indexed */
            const long CatalogueIndex = OutputIndex + 1;
            
            /* Get the index of the original lightcurve */
            long SourceIndex = indexOf(ObjectNames, Current.name);

            
            /* Copy the original data to the new location */
            outfile.moveHDU("HJD");
            long naxes[2];
            fits_get_img_size(*outfile.fptr(), 2, naxes, &outfile.status());
            outfile.check();


            

            /* And copy the other data parts two */
            vector<double> buffer(naxes[0]);
            
            /* Load the jd data separately */
            vector<double> jd(naxes[0]);
            outfile.moveHDU("HJD");
            fits_read_img(*outfile.fptr(), TDOUBLE, (SourceIndex*naxes[0])+1, naxes[0], 0, &jd[0], 0, &outfile.status());
            fits_write_img(*outfile.fptr(), TDOUBLE, OutputIndex*naxes[0], naxes[0], &jd[0], &outfile.status());
            outfile.moveHDU("FLUX");
            fits_read_img(*outfile.fptr(), TDOUBLE, (SourceIndex*naxes[0])+1, naxes[0], 0, &buffer[0], 0, &outfile.status());
            fits_write_img(*outfile.fptr(), TDOUBLE, OutputIndex*naxes[0], naxes[0], &buffer[0], &outfile.status());
            outfile.moveHDU("FLUXERR");
            fits_read_img(*outfile.fptr(), TDOUBLE, (SourceIndex*naxes[0])+1, naxes[0], 0, &buffer[0], 0, &outfile.status());
            fits_write_img(*outfile.fptr(), TDOUBLE, OutputIndex*naxes[0], naxes[0], &buffer[0], &outfile.status());
            outfile.moveHDU("CCDX");
            fits_read_img(*outfile.fptr(), TDOUBLE, (SourceIndex*naxes[0])+1, naxes[0], 0, &buffer[0], 0, &outfile.status());
            fits_write_img(*outfile.fptr(), TDOUBLE, OutputIndex*naxes[0], naxes[0], &buffer[0], &outfile.status());
            outfile.moveHDU("CCDY");
            fits_read_img(*outfile.fptr(), TDOUBLE, (SourceIndex*naxes[0])+1, naxes[0], 0, &buffer[0], 0, &outfile.status());
            fits_write_img(*outfile.fptr(), TDOUBLE, OutputIndex*naxes[0], naxes[0], &buffer[0], &outfile.status());
            outfile.moveHDU("SKYBKG");
            fits_read_img(*outfile.fptr(), TDOUBLE, (SourceIndex*naxes[0])+1, naxes[0], 0, &buffer[0], 0, &outfile.status());
            fits_write_img(*outfile.fptr(), TDOUBLE, OutputIndex*naxes[0], naxes[0], &buffer[0], &outfile.status());
            outfile.moveHDU("QUALITY");
            fits_read_img(*outfile.fptr(), TDOUBLE, (SourceIndex*naxes[0])+1, naxes[0], 0, &buffer[0], 0, &outfile.status());
            fits_write_img(*outfile.fptr(), TDOUBLE, OutputIndex*naxes[0], naxes[0], &buffer[0], &outfile.status());
            outfile.check();
            
            /* Now get all of the data from Index: OutputIndex*naxes[0] */

            if (Current.submodel_id != NullSubIndex)
            {
                statement subst(conn);
                Model SubModel;
                
                subst << "select id, name, period, epoch, a, i, rs, rp, mstar, c1, c2, c3, c4, teff "
                " from submodels where id = " << Current.submodel_id, into(SubModel.id), into(SubModel.name),  
                into(SubModel.period), into(SubModel.epoch), into(SubModel.a), into(SubModel.i), into(SubModel.rs),
                into(SubModel.rp), into(SubModel.mstar), into(SubModel.c1), into(SubModel.c2), into(SubModel.c3), 
                into(SubModel.c4), into(SubModel.teff);
                
                if (!subst.exec())
                {
                    throw runtime_error("Cannot find subtraction object");
                }
                
                /* SubModel now contains the subtraction model */
                AlterLightcurveData(outfile, OutputIndex*naxes[0], naxes[0], SubModel, ArithMeth("-"), Config);
            }
            
            /* Add a transit model to the data */
            AlterLightcurveData(outfile, OutputIndex*naxes[0], naxes[0], Current, ArithMeth("+"), Config);

            
            /* And update the catalogue false transits information */
            outfile.moveHDU("CATALOGUE");

            /* Need to do some conversion but have to create a temp variable for this */
            double tmp = Current.period * secondsInDay;
            fits_write_col(*outfile.fptr(), TDOUBLE, fcn.period, CatalogueIndex, 1, 1, &tmp, &outfile.status());
            
            tmp = Current.epoch;
            fits_write_col(*outfile.fptr(), TDOUBLE, fcn.epoch, CatalogueIndex, 1, 1, &tmp, &outfile.status());
            
            tmp = Current.rp * rJup;
            fits_write_col(*outfile.fptr(), TDOUBLE, fcn.rp, CatalogueIndex, 1, 1, &tmp, &outfile.status());
            
            tmp = Current.rs * rSun;
            fits_write_col(*outfile.fptr(), TDOUBLE, fcn.rs, CatalogueIndex, 1, 1, &tmp, &outfile.status());
            
            tmp = Current.a * AU;
            fits_write_col(*outfile.fptr(), TDOUBLE, fcn.a, CatalogueIndex, 1, 1, &tmp, &outfile.status());
            
            tmp = Current.i * radiansInDegree;
            fits_write_col(*outfile.fptr(), TDOUBLE, fcn.i, CatalogueIndex, 1, 1, &tmp, &outfile.status());
            
            double TransitDepth = square((Current.rp * rJup) / (Current.rs * rSun));
            double TransitWidth = WidthFromParams(Current);
            
            /* These next two require calculation */
            fits_write_col(*outfile.fptr(), TDOUBLE, fcn.width, CatalogueIndex, 1, 1, &TransitWidth, &outfile.status());
            fits_write_col(*outfile.fptr(), TDOUBLE, fcn.depth, CatalogueIndex, 1, 1, &TransitDepth, &outfile.status());
            
            /* Now the skipdet flag */
            int SkipdetFlag = AlterDetrending::skipboth;
            fits_write_col(*outfile.fptr(), TINT, fcn.skipdet, CatalogueIndex, 1, 1, &SkipdetFlag, &outfile.status());

            /* Also write the original lightcurve's skipdet flag, but check what it is first*/
            int OldSkipdetFlag;
            fits_read_col(*outfile.fptr(), TINT, fcn.skipdet, SourceIndex+1, 1, 1, NULL, &OldSkipdetFlag, NULL, &outfile.status());

            /* Only update if it's set to AlterDetrending::include */
            if (OldSkipdetFlag == AlterDetrending::include)
            {
                /* New value of the flag */
                SkipdetFlag = AlterDetrending::skiptfa;
                fits_write_col(*outfile.fptr(), TINT, fcn.skipdet, SourceIndex + 1, 1, 1, &SkipdetFlag, &outfile.status());
            }



            /* Write the new name to the obj_id column */
            string NewName = AlterObjectName(Current.name);
            char *cstr = new char[26];
            strcpy(cstr, NewName.c_str());
            fits_write_col_str(*outfile.fptr(), obj_id_colno, CatalogueIndex, 1, 1, &cstr, &outfile.status());
            delete[] cstr;
            
            /* Now validate */
            outfile.check();



            
            stringstream ss;
            ss << counter + 1 << "/" << nextra;
            OverPrint(ss.str());
            
            ++counter;
        }

        cout << endl;



        ts.stop("model.iterate");

        ts.stop("all");
        return 0;
    }
    catch (TCLAP::ArgException &e)
    {
        cerr << "error: " << e.error() << " for arg " << e.argId() << endl;
    }
    catch (std::runtime_error &e)
    {
        cerr << "Runtime error: " << e.what() << endl;
    }
    catch (std::exception &e)
    {
        cerr << "std::exception: " << e.what() << endl;
    }

    return EXIT_FAILURE;
}
