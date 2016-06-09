#include <iostream>
#include <cstring>
#include <cassert>
#include <stdexcept>
#include <algorithm>
#include <fitsio.h>
#include <tclap/CmdLine.h>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <sqlite3.h>

/* Local includes */
#include "timer.h"
#include "GenerateModel.h"
#include "Model.h"
#include "FitsObject.h"
#include "constants.h"
#include "ObjectSkipDefs.h"



using namespace std;

typedef vector<string> stringlist;

enum {
    add = 101,
    sub
};

struct ConfigContainer {
    bool isWASP;
    // Equivalent to isWASP...
    bool isNGTS;
    string DatabaseFilename;
    string SourceFilename;
    string OutputFilename;

    bool isWASPLike() const {
        return isWASP || isNGTS;
    }
};

template <typename T>
T square(T val) {
    return val * val;
}

const double jd_ref = 2456658.500000;

double wd2jd(double wd) {
    return (wd / secondsInDay) + jd_ref;
}

double jd2wd(double jd) {
    return (jd - jd_ref) * secondsInDay;
}


class ArithMeth {

        int m_type;

    public:
        ArithMeth(const string &type) {
            if (type == "+") {
                this->m_type = add;

            } else if (type == "-") {
                this->m_type = sub;

            } else {
                throw runtime_error("Unknown method specified (must be + or -)");
            }
        }

        int type() const {
            return this->m_type;
        }
};

double WidthFromParams(const Model &m) {
    /* Returns the width of the full transit based on some lc parameters
     *
     * \frac{P}{\pi} \asin{\sqrt{(\frac{R_P + R_S}{a})^2 - \cos^2 i}}
     */
    const double Norm = (m.period * secondsInDay) / M_PI;

    if (m.a == 0) {
        /* Something hasn't updated the separation */
        return 0;
    }

    const double FirstTerm = square(((m.rp * rJup) + (m.rs * rSun)) / (m.a * AU));
    const double SecondTerm = square(cos(m.i * radiansInDegree));
    const double InsideSqrt = FirstTerm - SecondTerm;

    /* Check that InsideSqrt is not <= 0 */
    if (InsideSqrt <= 0.0) {
        return 0.0;
    }

    const double SqrtInsideSqrt = sqrt(InsideSqrt);

    if ((SqrtInsideSqrt < -1.0) || (SqrtInsideSqrt > 1.0)) {
        return 0;
    }



    return Norm * asin(SqrtInsideSqrt);

}

template <typename T>
T WeightedMedian(const vector<T> &data, const double siglevel, long &npts) {

    vector<T> buffer;
    const size_t N = data.size();

    /* First calculate the mean */
    double av = 0;
    int ValidPoints = 0;

    for (size_t i = 0; i < N; ++i) {
        if (!isnan(data.at(i))) {

            av += data.at(i);
            ValidPoints++;

        }
    }

    av /= (double)ValidPoints;
    npts = ValidPoints;

    /* Put in a check for if the av is 0 */
    if (av == 0) {
        throw runtime_error("Cannot calculate average - all points are 0");
    }

    if (isnan(av)) {
        throw runtime_error("Cannot calculate average - average is NaN");
    }

    /* now the sigma */
    double sd = 0;

    for (size_t i = 0; i < N; ++i) {
        if (!isnan(data.at(i))) {
            sd += (data.at(i) - av) * (data.at(i) - av);
        }
    }

    sd /= (double)ValidPoints;
    sd = sqrt(sd);

    const double upperlim = av + siglevel * sd;
    const double lowerlim = av - siglevel * sd;

    for (int i = 0; i < N; ++i) {
        if ((data.at(i) < upperlim) && (data.at(i) > lowerlim) && !isnan(data.at(i))) {
            buffer.push_back(data.at(i));
        }
    }

//    cout << N - buffer.size() << " elements rejected" << endl;

    /* Sort the array */
    sort(buffer.begin(), buffer.end());

    /* Return the middle one */
    return buffer.at(buffer.size() / 2);


}



struct FalseColumnNumbers {
    int real_obj_id, skipdet, period, width, depth, epoch,
        rp, rs, a, i;
};

string zeroPad(const string &s, int width) {
    stringstream ss;
    ss << setw(width) << setfill('0') << s;
    return ss.str();
}


long indexOf(const stringlist &stringlist, const string &comp) {

    for (size_t i = 0; i < stringlist.size(); ++i) {

        if (stringlist.at(i) == comp) {
            return i;

        } else {
            string nameWithoutWhitespace(comp);
            nameWithoutWhitespace.erase(remove_if(nameWithoutWhitespace.begin(), nameWithoutWhitespace.end(), ::isspace), nameWithoutWhitespace.end());

            if (stringlist.at(i) == nameWithoutWhitespace) {
                return i;
            } else {
                // Zero pad
                string nameZeroPadded = zeroPad(nameWithoutWhitespace, 6);
                if (zeroPad(stringlist.at(i), 6) == zeroPad(nameWithoutWhitespace, 6)) {
                    return i;
                }
            }
        }
    }


    /* If the loop gets here the object is not found */
    throw runtime_error("Cannot find object");
}

pair<double, long> AlterLightcurveData(Fits &f, const long startindex, const int length, const Model &m, const ArithMeth &arithtype, const ConfigContainer &Config) {
    /* returns a pair of the mean flux and the number
     of valid points in the lightcurve */
    f.moveHDU("HJD");
    pair<double, long> OutputData;

    vector<double> jd(length);

    /* Fetch the jd data */
    fits_read_img(*f.fptr(), TDOUBLE, startindex, length, 0, &jd[0], 0, &f.status());

    // if the Config.isWASPLike parameter is set then convert the array to jd
    if (Config.isWASPLike()) {
        for (int i = 0; i < length; ++i) {
            jd[i] = wd2jd(jd[i]);
        }
    }

    /* Now get the addition model */
    vector<double> ModelFlux = GenerateSynthetic(jd, m);

    f.moveHDU("FLUX");
    vector<double> OriginalFlux(length);
    fits_read_img(*f.fptr(), TDOUBLE, startindex, length, 0, &OriginalFlux[0], 0, &f.status());

    /* Normalise to the weighted median flux level */
    double WeightedMed = WeightedMedian(OriginalFlux, 2.5, OutputData.second);
    OutputData.first = WeightedMed;

    vector<double> TransitAdded(length);
    double result = 0;

    for (int i = 0; i < length; ++i) {
        double fluxval = OriginalFlux.at(i);
        fluxval /= WeightedMed;

        double modelval = ModelFlux.at(i);

        /* Perform the arithmetic */

        if (arithtype.type() == add) {
            result = WeightedMed * (fluxval + modelval - 1.0);

        } else if (arithtype.type() == sub) {
            result = WeightedMed * (fluxval - modelval + 1.0);
        }

        TransitAdded[i] = result;
    }




    fits_write_img(*f.fptr(), TDOUBLE, startindex, length, &TransitAdded[0], &f.status());
    f.check();

    return OutputData;


}

stringlist &split(const string &s, char delim, stringlist &elems) {
    stringstream ss(s);
    string item;

    while (getline(ss, item, delim)) {
        elems.push_back(item);
    }

    return elems;
}


stringlist split(const string &s, char delim) {
    stringlist elems;
    return split(s, delim, elems);
}

/** Takes an object name and creates a new one

This function alters the object name to be unique in the
catalogue list. */
string AlterObjectName(const string &OriginalName) {
    /* The name MUST be unique as Orion uses a dictionary to see
    if the object already exists. Therefore a new naming scheme must
    be created.  */

    /* Counter variable to ensure uniqueness */
    static unsigned long counter;

    /* Fake naming scheme:
     *
     * NGTS objects are indexed by integers. We have 26 characters to play with
     * and need to have a unique name per object.
     */

    char buf[6];
    snprintf(buf, 6, "F%05d", counter);
    string ResultingString(buf);


#if 0

    /* Split the string at the J character */
    stringlist parts = split(OriginalName, 'J');

    /* check the length for validity */
    if (parts.size() != 2) {
        throw runtime_error("Unknown object name encountered");
    }

    stringstream NewName;
    NewName.flags(ios::left);
    NewName << setw(7) << counter << "J" << parts[1];
    string ResultingString = NewName.str();

    /* Allow up to 10 million objects to be created */
    if (ResultingString.size() > 26) {
        throw runtime_error("Too many objects in file (>10000000)");
    }

#endif

    /* Increment the counter */
    counter++;



    return ResultingString;
}

template <typename T>
void OverPrint(const T &val) {
    cout << "\r" << val;
    cout.flush();

}

/* Function to copy a row from a table across to a new position in the same table */
void CopyTableRow(Fits &infile, const long origindex, const long newindex) {
    /* Current hdu must be on a binary table */
    const long nrows = infile.nrows();

    /* Get the number of columns */
    int ncols;
    fits_get_num_cols(*infile.fptr(), &ncols, &infile.status());
    infile.check();

    for (int i = 1; i <= ncols; ++i) {
        /* Check for data type */
        int typecode;

        /* Get the repeat count for each cell */
        long repeat;

        /* Get the data width (unused) */
        long width;
        fits_get_coltype(*infile.fptr(), i, &typecode, &repeat, &width, &infile.status());
        infile.check();

        switch (typecode) {
            case TDOUBLE: {
                double tmp;
                fits_read_col(*infile.fptr(), TDOUBLE, i, origindex, 1, 1, NULL, &tmp, NULL, &infile.status());
                fits_write_col(*infile.fptr(), TDOUBLE, i, newindex, 1, 1, &tmp, &infile.status());
                infile.check();
                break;
            }

            case TSTRING: {
                char buf[repeat + 1], *bufptr = (char *)buf;

                fits_read_col_str(*infile.fptr(), i, origindex, 1, 1, "", &bufptr, NULL, &infile.status());
                fits_write_col_str(*infile.fptr(), i, newindex, 1, 1, &bufptr, &infile.status());
                infile.check();
                break;
            }

            case TLONG: {
                long tmp;
                fits_read_col(*infile.fptr(), TLONG, i, origindex, 1, 1, NULL, &tmp, NULL, &infile.status());
                fits_write_col(*infile.fptr(), TLONG, i, newindex, 1, 1, &tmp, &infile.status());
                infile.check();
                break;
            }

            case TINT: {
                int tmp;
                fits_read_col(*infile.fptr(), TINT, i, origindex, 1, 1, NULL, &tmp, NULL, &infile.status());
                fits_write_col(*infile.fptr(), TINT, i, newindex, 1, 1, &tmp, &infile.status());
                infile.check();
                break;
            }

            case TFLOAT: {
                float tmp;
                fits_read_col(*infile.fptr(), TFLOAT, i, origindex, 1, 1, NULL, &tmp, NULL, &infile.status());
                fits_write_col(*infile.fptr(), TFLOAT, i, newindex, 1, 1, &tmp, &infile.status());
                infile.check();
                break;
            }

            case TSHORT: {
                short int tmp;
                fits_read_col(*infile.fptr(), TSHORT, i, origindex, 1, 1, NULL, &tmp, NULL, &infile.status());
                fits_write_col(*infile.fptr(), TSHORT, i, newindex, 1, 1, &tmp, &infile.status());
                infile.check();
                break;
            }

            case TBYTE: {
                unsigned char tmp;
                fits_read_col(*infile.fptr(), TBYTE, i, origindex, 1, 1, NULL, &tmp, NULL, &infile.status());
                fits_write_col(*infile.fptr(), TBYTE, i, newindex, 1, 1, &tmp, &infile.status());
                infile.check();
                break;
            }

            default:
                stringstream ss;
                ss << "Unknown column type found: " << typecode;
                throw runtime_error(ss.str());
        }

    }


}

class FetchesParameters {
    public:

        FetchesParameters(const string &filename)
            : filename(filename) {
            sqlite3_open(filename.c_str(), &ppDb);
        }

        ~FetchesParameters() {
            if (ppDb != NULL) {
                sqlite3_close(ppDb);
            }
        }

        long nmodels() const {
            const string sql = "SELECT count(*) FROM addmodels;";
            sqlite3_stmt *stmt;
            sqlite3_prepare_v2(ppDb, sql.c_str(), sql.size(), &stmt, NULL);
            int step_state = 0;
            int nrows = 0;
            while (step_state != SQLITE_DONE) {
                switch (step_state) {
                    case SQLITE_ROW:
                        nrows = sqlite3_column_int(stmt, 0);
                        break;
                    case SQLITE_ERROR:
                        fprintf(stderr, "Error\n");
                        exit(1);
                }
                step_state = sqlite3_step(stmt);
            }

            sqlite3_finalize(stmt);

            return nrows;
        }

        Model model_from_statement(sqlite3_stmt *stmt) {
            Model model;
            int col = 0;
            model.id = sqlite3_column_int(stmt, col++);

            const unsigned char *model_name_cstr = sqlite3_column_text(stmt, col++);
            string model_name(reinterpret_cast<const char*>(model_name_cstr));
            model.name = model_name;

            model.submodel_id = sqlite3_column_int(stmt, col++);
            model.period = sqlite3_column_double(stmt, col++);
            model.epoch = sqlite3_column_double(stmt, col++);
            model.a =  sqlite3_column_double(stmt, col++);
            model.i = sqlite3_column_double(stmt, col++);
            model.rs = sqlite3_column_double(stmt, col++);
            model.rp = sqlite3_column_double(stmt, col++);
            model.mstar = sqlite3_column_double(stmt, col++);
            model.c1 = sqlite3_column_double(stmt, col++);
            model.c2 = sqlite3_column_double(stmt, col++);
            model.c3 = sqlite3_column_double(stmt, col++);
            model.c4 = sqlite3_column_double(stmt, col++);
            model.teff = sqlite3_column_double(stmt, col++);

            return model;
        }

        // TODO: why are some names NULL?
        vector<Model> fetch_models() {
            const string sql = "select id, name, submodel_id, period, epoch, a, "
                "i, rs, rp, mstar, c1, c2, c3, c4, teff "
                "from addmodels "
                "where name is not null "
                "order by name asc"
                ;
            sqlite3_stmt *stmt;
            sqlite3_prepare_v2(ppDb, sql.c_str(), sql.size(), &stmt, NULL);
            vector<Model> models;
            int step_state = 0;
            while (step_state != SQLITE_DONE) {
                switch (step_state) {
                    case SQLITE_ROW:
                        models.push_back(model_from_statement(stmt));
                        break;
                    case SQLITE_ERROR:
                        fprintf(stderr, "Error\n");
                        exit(1);
                }
                step_state = sqlite3_step(stmt);
            }

            sqlite3_finalize(stmt);

            return models;
        }

        struct StopIteration : public runtime_error {};

    private:
        string filename;
        sqlite3 *ppDb;

};

stringlist extract_object_names(ReadOnlyFits &infile, long &nrows) {
    stringlist ObjectNames;
    infile.moveHDU("CATALOGUE");
    int obj_id_colno = infile.columnNumber("OBJ_ID");

    /* Read the data in as strings */
    int dispwidth;
    fits_get_col_display_width(*infile.fptr(), obj_id_colno, &dispwidth, &infile.status());

    nrows = 0;
    fits_get_num_rows(*infile.fptr(), &nrows, &infile.status());

    vector<char *> cstrnames(nrows);

    for (int i = 0; i < nrows; ++i) {
        cstrnames[i] = new char[dispwidth + 1];
    }

    fits_read_col_str(*infile.fptr(), obj_id_colno, 1, 1, nrows, 0, &cstrnames[0], 0, &infile.status());

    for (int i = 0; i < nrows; ++i) {
        string CurrentName = cstrnames[i];
        /* Remove all whitespace */
        CurrentName.erase(remove_if(CurrentName.begin(), CurrentName.end(), ::isspace), CurrentName.end());
        ObjectNames.push_back(CurrentName);
        delete[] cstrnames[i];
    }

    return ObjectNames;
}

/* Copy from the documentation as ICC does not support this */
template< class InputIt, class UnaryPredicate >
bool any_of(InputIt first, InputIt last, UnaryPredicate p) {
    return std::find_if(first, last, p) != last;
}


/* A lightcurve is valid if any flux points are positive
 */
bool valid_lightcurve(const vector<double> &flux) {
    return any_of(flux.begin(), flux.end(), [](const double f) {
            return f > 0.;
    });
}

vector<Model> compute_valid_extra_models(const vector<Model> &models, ReadOnlyFits &infile) {
    vector<Model> out;

    long nrows;
    stringlist object_names = extract_object_names(infile, nrows);

    infile.moveHDU("FLUX");
    long naxes[2];
    fits_get_img_size(*infile.fptr(), 2, naxes, &infile.status());
    infile.check();

    vector<double> buffer(naxes[0]);
    cout << "Computing the list of valid lightcurves" << endl;
    for (auto itr=models.begin(); itr!=models.end(); itr++) {
        auto object_name = itr->name;
        /* Get the index of the original lightcurve */
        long SourceIndex = indexOf(object_names, object_name);

        fits_read_img(*infile.fptr(), TDOUBLE, (SourceIndex * naxes[0]) + 1, naxes[0], 0, &buffer[0], 0, &infile.status());
        infile.check();

        if (valid_lightcurve(buffer)) {
            out.push_back(*itr);
        }
    }

    return out;
}

int main(int argc, char *argv[]) {
    try {
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

        ReadOnlyFits infile(infile_arg.getValue());

        /* Make sure the project has the 'project' header key,
         * and if it's wasp make sure all jds are converted
         * to wasp dates
         */
        infile.moveHDU(1);

        char Project_cstr[FLEN_VALUE];
        fits_read_key(*infile.fptr(), TSTRING, "PROJECT", &Project_cstr, 0, &infile.status());

        /* check the key exists */
        if (infile.status() == KEY_NO_EXIST) {
            cout << "Project unspecified, defaulting to WASP" << endl;
            Config.isWASP = true;
            infile.status() = 0;

        } else {
            /* Create string object for easy comparisons */
            string Project(Project_cstr);

            cout << Project << " project found" << endl;

            Config.isWASP = (Project == "WASP") ? true : false;
            Config.isNGTS = (Project == "NGTS") ? true : false;
        }

        Config.DatabaseFilename = candidates_arg.getValue();
        Config.SourceFilename = infile_arg.getValue();
        Config.OutputFilename = output_arg.getValue();

        if (Config.isWASPLike()) {
            cout << "--- Converting times to WASP data" << endl;
        }


        NewFits outfile("!" + output_arg.getValue());

        /* Copy the primary hdu across */
        int status = 0;
        fits_copy_hdu(*infile.fptr(), *outfile.fptr(), 0, &status);
        Fits::check(status);

        /* Add the transinj key */
        bool transinj_val = true;
        fits_write_key(*outfile.fptr(), TLOGICAL, "TRANSINJ", &transinj_val, "Contains false transits", &outfile.status());
        outfile.check();


        /* Start by getting file information from the input */
        int nhdus = 0;
        fits_get_num_hdus(*infile.fptr(), &nhdus, &infile.status());
        infile.check();

        cout << nhdus << " hdus found" << endl;


        /* Open the sqlite3 database here */
        FetchesParameters param_fetcher(candidates_arg.getValue());

        /* get the new models to insert */
        vector<Model> models = param_fetcher.fetch_models();
        int nextra = models.size();

        if (nextra == 0) {
            throw runtime_error("No input models found");
        }

        /* Some lightcurves have no data, so compute the new number of extra
         * objects */
        vector<Model> valid_models = compute_valid_extra_models(models, infile);
        int valid_nextra = valid_models.size();

        cout << "Inserting " << valid_nextra << " valid extra models (found " << nextra << " in total)" << endl;

        for (int hdu = 2; hdu <= nhdus; ++hdu) {
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
            if (hdutype == IMAGE_HDU) {
                /* Get the current dimensions */
                long naxes[2];
                fits_get_img_size(*outfile.fptr(), 2, naxes, &outfile.status());
                outfile.check();

                cout << " - " << naxes[0] << "x" << naxes[1] << " pix";

                int bitpix = 0;
                fits_get_img_type(*outfile.fptr(), &bitpix, &outfile.status());
                outfile.check();

                long newnaxes[] = {naxes[0], naxes[1] + valid_nextra};

                fits_resize_img(*outfile.fptr(), bitpix, 2, newnaxes, &outfile.status());
                outfile.check();

                cout << " -> " << newnaxes[0] << "x" << newnaxes[1] << " pix";

            } else {
                /* If the catalogue extension is found then add extra rows */
                const string hduname = outfile.hduname();
                long nrows = outfile.nrows();

                cout << " - " << nrows << " rows";

                if (hduname == "CATALOGUE") {
                    /* Need to add extra columns */
                    int ncols = 0;
                    fits_get_num_cols(*outfile.fptr(), &ncols, &outfile.status());
                    outfile.check();

                    cout << ", " << ncols << " columms";

                    /* Append extra rows */

                    fits_insert_rows(*outfile.fptr(), nrows, valid_nextra, &outfile.status());
                    outfile.check();

                    cout << " -> " << nrows + valid_nextra << " rows";

                    /* Need 9 extra columns */

                    char *ColumnNames[] = {"REAL_OBJ_ID", "SKIPDET", "FAKE_PERIOD", "FAKE_WIDTH", "FAKE_DEPTH", "FAKE_EPOCH", "FAKE_RP", "FAKE_RS", "FAKE_A", "FAKE_I"};
                    char *ColumnFormats[] = {"6A", "1I", "1D", "1D", "1D", "1D", "1D", "1D", "1D", "1D"};

                    size_t nNewCols = sizeof(ColumnNames) / sizeof(char *);
                    assert((sizeof(ColumnFormats) / sizeof(char *)) == nNewCols);

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

        /* Prefetch the column numbers for later use */
        outfile.moveHDU("CATALOGUE");
        FalseColumnNumbers fcn;
        fcn.real_obj_id = outfile.columnNumber("REAL_OBJ_ID");
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
        const int NullSubIndex = -1;


        /* Get a list of the objects in the file */
        long nrows;
        stringlist ObjectNames = extract_object_names(infile, nrows);

        long counter = 0;

        for (vector<Model>::const_iterator itr = valid_models.begin();
                itr != valid_models.end();
                itr++) {
            const Model Current = *itr;

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
            fits_read_img(*outfile.fptr(), TDOUBLE, (SourceIndex * naxes[0]) + 1, naxes[0], 0, &jd[0], 0, &outfile.status());
            fits_write_img(*outfile.fptr(), TDOUBLE, OutputIndex * naxes[0], naxes[0], &jd[0], &outfile.status());
            outfile.moveHDU("FLUX");
            fits_read_img(*outfile.fptr(), TDOUBLE, (SourceIndex * naxes[0]) + 1, naxes[0], 0, &buffer[0], 0, &outfile.status());
            fits_write_img(*outfile.fptr(), TDOUBLE, OutputIndex * naxes[0], naxes[0], &buffer[0], &outfile.status());
            outfile.moveHDU("FLUXERR");
            fits_read_img(*outfile.fptr(), TDOUBLE, (SourceIndex * naxes[0]) + 1, naxes[0], 0, &buffer[0], 0, &outfile.status());
            fits_write_img(*outfile.fptr(), TDOUBLE, OutputIndex * naxes[0], naxes[0], &buffer[0], &outfile.status());
            outfile.moveHDU("CCDX");
            fits_read_img(*outfile.fptr(), TDOUBLE, (SourceIndex * naxes[0]) + 1, naxes[0], 0, &buffer[0], 0, &outfile.status());
            fits_write_img(*outfile.fptr(), TDOUBLE, OutputIndex * naxes[0], naxes[0], &buffer[0], &outfile.status());
            outfile.moveHDU("CCDY");
            fits_read_img(*outfile.fptr(), TDOUBLE, (SourceIndex * naxes[0]) + 1, naxes[0], 0, &buffer[0], 0, &outfile.status());
            fits_write_img(*outfile.fptr(), TDOUBLE, OutputIndex * naxes[0], naxes[0], &buffer[0], &outfile.status());
            outfile.moveHDU("SKYBKG");
            fits_read_img(*outfile.fptr(), TDOUBLE, (SourceIndex * naxes[0]) + 1, naxes[0], 0, &buffer[0], 0, &outfile.status());
            fits_write_img(*outfile.fptr(), TDOUBLE, OutputIndex * naxes[0], naxes[0], &buffer[0], &outfile.status());
            outfile.moveHDU("QUALITY");
            fits_read_img(*outfile.fptr(), TDOUBLE, (SourceIndex * naxes[0]) + 1, naxes[0], 0, &buffer[0], 0, &outfile.status());
            fits_write_img(*outfile.fptr(), TDOUBLE, OutputIndex * naxes[0], naxes[0], &buffer[0], &outfile.status());
            outfile.check();

            /* Add a transit model to the data */
            pair<double, long> LightcurveInfo = AlterLightcurveData(outfile, OutputIndex * naxes[0], naxes[0], Current, ArithMeth("+"), Config);


            /* And update the catalogue false transits information */
            outfile.moveHDU("CATALOGUE");

            /* Copy the original object data across */
            CopyTableRow(outfile, SourceIndex + 1, CatalogueIndex);

            /* Need to do some conversion but have to create a temp variable for this */
            double tmp = Current.period * secondsInDay;
            fits_write_col(*outfile.fptr(), TDOUBLE, fcn.period, CatalogueIndex, 1, 1, &tmp, &outfile.status());

            tmp = jd2wd(Current.epoch);
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
            fits_read_col(*outfile.fptr(), TINT, fcn.skipdet, SourceIndex + 1, 1, 1, NULL, &OldSkipdetFlag, NULL, &outfile.status());

            /* Only update if it's set to AlterDetrending::include */
            if (OldSkipdetFlag == AlterDetrending::include) {
                /* New value of the flag */
                SkipdetFlag = AlterDetrending::skiptfa;
                fits_write_col(*outfile.fptr(), TINT, fcn.skipdet, SourceIndex + 1, 1, 1, &SkipdetFlag, &outfile.status());
            }



            /* Write the new name to the obj_id column */
            string NewName = AlterObjectName(Current.name);
            char *cstr = new char[6];
            strncpy(cstr, NewName.c_str(), 6);
            fits_write_col_str(*outfile.fptr(), fcn.real_obj_id, CatalogueIndex, 1, 1, &cstr, &outfile.status());
            delete[] cstr;

            /* Update the flux mean and npts columns */
            int npts_col = outfile.columnNumber("NPTS");
            int flux_mean_col = outfile.columnNumber("FLUX_MEAN");
            fits_write_col(*outfile.fptr(), TLONG, npts_col, CatalogueIndex, 1, 1, &LightcurveInfo.second, &outfile.status());
            fits_write_col(*outfile.fptr(), TDOUBLE, flux_mean_col, CatalogueIndex, 1, 1, &LightcurveInfo.first, &outfile.status());

            /* Now validate */
            outfile.check();




            stringstream ss;
            ss << counter + 1 << "/" << valid_nextra;
            OverPrint(ss.str());

            ++counter;
        }

        cout << endl;



        ts.stop("model.iterate");

        ts.stop("all");
        return 0;

    } catch (TCLAP::ArgException &e) {
        cerr << "error: " << e.error() << " for arg " << e.argId() << endl;

    } catch (std::runtime_error &e) {
        cerr << "Runtime error: " << e.what() << endl;

    } catch (std::exception &e) {
        cerr << "std::exception: " << e.what() << endl;
    }

    return EXIT_FAILURE;
}
