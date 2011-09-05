#include <iostream>
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
    
    double FirstTerm = ((m.rp * rJup) + (m.rs * rSun)) / (m.a * AU);
    
    /* Square it */
    FirstTerm *= FirstTerm;
    
    const double InsideSqrt = FirstTerm - cos(m.i * radiansInDegree);
    
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
    for (size_t i=0; i<N; ++i)
    {
        av += data.at(i);
    }
    av /= (double)N;

    /* now the sigma */
    double sd = 0;
    for (size_t i=0; i<N; ++i)
    {
        sd += (data.at(i) - av)*(data.at(i) - av);
    }
    sd /= (double)N;
    sd = sqrt(sd);

    const double upperlim = av + siglevel * sd;
    const double lowerlim = av - siglevel * sd;

    for (int i=0; i<N; ++i)
    {
        if ((data.at(i) < upperlim) && (data.at(i) > lowerlim))
        {
            buffer.push_back(data.at(i));
        }
    }

    cout << N - buffer.size() << " elements rejected" << endl;

    /* Sort the array */
    sort(buffer.begin(), buffer.end());

    /* Return the middle one */
    return buffer.at(buffer.size() / 2);


}


template <typename T>
T square(T val) { return val * val; }

struct FalseColumnNumbers
{
    int skipdet, period, width, depth, epoch,
    rp, rs, a, i;
};


long indexOf(const vector<string> &stringlist, const string &comp)
{
    for (size_t i=0; i<stringlist.size(); ++i)
    {
        if (stringlist.at(i) == comp)
        {
            return i;
        }
    }

    /* If the loop gets here the object is not found */
    throw runtime_error("Cannot find object");
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

        Fits infile(infile_arg.getValue());
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
                long nrows = 0;
                fits_get_num_rows(*outfile.fptr(), &nrows, &outfile.status());
                outfile.check();

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
        fits_get_colnum(*outfile.fptr(), CASEINSEN, "SKIPDET", &fcn.skipdet, &outfile.status());
        fits_get_colnum(*outfile.fptr(), CASEINSEN, "FAKE_PERIOD", &fcn.period, &outfile.status());
        fits_get_colnum(*outfile.fptr(), CASEINSEN, "FAKE_WIDTH", &fcn.width, &outfile.status());
        fits_get_colnum(*outfile.fptr(), CASEINSEN, "FAKE_DEPTH", &fcn.depth, &outfile.status());
        fits_get_colnum(*outfile.fptr(), CASEINSEN, "FAKE_EPOCH", &fcn.epoch, &outfile.status());
        fits_get_colnum(*outfile.fptr(), CASEINSEN, "FAKE_RP", &fcn.rp, &outfile.status());
        fits_get_colnum(*outfile.fptr(), CASEINSEN, "FAKE_RS", &fcn.rs, &outfile.status());
        fits_get_colnum(*outfile.fptr(), CASEINSEN, "FAKE_A", &fcn.a, &outfile.status());
        fits_get_colnum(*outfile.fptr(), CASEINSEN, "FAKE_I", &fcn.i, &outfile.status());
        outfile.check();




        ts.start("model.iterate");
        Model Current;
        const int NullSubIndex = -1;


        /* Get a list of the objects in the file */
        vector<string> ObjectNames;
        infile.moveHDU("CATALOGUE");

        int obj_id_colno = -1;
        fits_get_colnum(*infile.fptr(), CASEINSEN, "OBJ_ID", &obj_id_colno, &infile.status());
        infile.check();

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


        int counter = 0;
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
            }

            /* Now get the addition model */
            vector<double> ModelFlux = GenerateSynthetic(jd, Current);

            outfile.moveHDU("FLUX");
            vector<double> OriginalFlux(naxes[0]);
            fits_read_img(*outfile.fptr(), TDOUBLE, OutputIndex*naxes[0], naxes[0], 0, &OriginalFlux[0], 0, &outfile.status());

            /* Normalise to the weighted median flux level */
            double WeightedMed = WeightedMedian(OriginalFlux, 2.5);
            
            vector<double> TransitAdded(naxes[0]);
            for (int i=0; i<naxes[0]; ++i)
            {
                double fluxval = OriginalFlux.at(i);
                fluxval /= WeightedMed;
                
                double modelval = ModelFlux.at(i);
                
                /* Add the flux */
                double result = WeightedMed * (fluxval + modelval - 1.0);
                TransitAdded[i] = result;
            }


            
            
            fits_write_img(*outfile.fptr(), TDOUBLE, OutputIndex * naxes[0], naxes[0], &TransitAdded[0], &outfile.status());
            outfile.check();
            
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
            
            /* Now validate */
            outfile.check();



            
            
            ++counter;
        }



        ts.stop("model.iterate");

        ts.stop("all");
        return 0;
    }
    catch (TCLAP::ArgException &e)
    {
        cerr << "error: " << e.error() << " for arg " << e.argId() << endl;
    }
    catch (std::exception &e)
    {
        cerr << e.what() << endl;
    }

    return EXIT_FAILURE;
}
