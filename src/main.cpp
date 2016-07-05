#include <iostream>
#include <cstring>
#include <cassert>
#include <stdexcept>
#include <map>
#include <algorithm>
#include <fitsio.h>
#include <tclap/CmdLine.h>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <cstdlib>
#include <memory>
#include <fstream>
#include <sqlite3.h>
#include <batman.h>

/* Local includes */
#include "timer.h"
#include "GenerateModel.h"
#include "Model.h"
#include "FitsObject.h"
#include "constants.h"
#include "ObjectSkipDefs.h"
#include "fetches_parameters.h"



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
    /* Returns the width of the full transit in seconds based on some lc parameters
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

// Return the 6 character zero padded name
string sanitise_object_name(const string &name) {
    return zeroPad(name, 6);
}

vector<double> generate_model(const vector<double> &hjd, const Model &model);

pair<double, long> AlterLightcurveData(Fits &f, const string &models_filename, const long startindex, const int length, const Model &m, const ArithMeth &arithtype, const ConfigContainer &Config) {
    /* returns a pair of the mean flux and the number
     of valid points in the lightcurve */
    f.moveHDU("HJD");
    pair<double, long> OutputData;

    vector<double> jd(length);

    /* Fetch the jd data */
    fits_read_img(*f.fptr(), TDOUBLE, startindex, length, 0, &jd[0], 0, &f.status());

    /* Convert to seconds if required */
    if (!Config.isWASPLike()) {
        for (int i = 0; i < length; ++i) {
            jd[i] = jd2wd(jd[i]);
        }
    }

    /* Now get the addition model */
    vector<double> ModelFlux = generate_model(jd, m);

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
string GenerateNewObjectName(long counter) {
    /* The name MUST be unique as Orion uses a dictionary to see
    if the object already exists. Therefore a new naming scheme must
    be created.  */

    /* Fake naming scheme:
     *
     * NGTS objects are indexed by integers. We have 26 characters to play with
     * and need to have a unique name per object.
     */

    char buf[7];
    snprintf(buf, 7, "F%05lu", counter);
    string ResultingString(buf);

    /* Increment the counter */
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

map<string, unsigned int> extract_object_names(ReadOnlyFits &infile, long &nrows) {
    map<string, unsigned int> ObjectNames;
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
        string CurrentName = sanitise_object_name(cstrnames[i]);
        /* Remove all whitespace */
        CurrentName.erase(remove_if(CurrentName.begin(), CurrentName.end(), ::isspace), CurrentName.end());
        ObjectNames.insert(make_pair(CurrentName, i));
        delete[] cstrnames[i];
    }

    return ObjectNames;
}

/* Copy from the documentation as ICC does not support this */
template< class InputIt, class UnaryPredicate >
bool is_any_of(InputIt first, InputIt last, UnaryPredicate p) {
    return std::find_if(first, last, p) != last;
}


/* A lightcurve is valid if 70% of the flux points are positive
 */
bool valid_lightcurve(const vector<double> &hjd, const vector<double> &flux,
        double width_seconds, double period_days, double epoch_days) {

    int lc_size = flux.size();
    // First check the flux values
    int nvalid_points = std::count_if(flux.begin(), flux.end(), [](const double f) {
                return f > 0.0;
            });

    if (!(float(nvalid_points) / float(lc_size)) >= 0.7) {
        return false;
    }

    // If this is ok, make sure there is at least some data during the transit
    double epoch_seconds = jd2wd(epoch_days);
    double period_seconds = period_days * secondsInDay;
    double width_in_phase = width_seconds / period_seconds;

    int npoints_in_transit = 0;
    for (int i=0; i<hjd.size(); i++) {
        double phase = fmod((hjd[i] - epoch_seconds) / period_seconds, 1);

        if ((phase >= -width_in_phase / 2.) && (phase <= width_in_phase / 2.)) {
            npoints_in_transit++;
        }
    }

    if (npoints_in_transit < 10) {
        return false;
    }

    return true;
}

vector<Model> compute_valid_extra_models(const vector<Model> &models, ReadOnlyFits &infile) {
    vector<Model> out;

    long nrows;
    map<string, unsigned int> object_names = extract_object_names(infile, nrows);

    infile.moveHDU("FLUX");
    long naxes[2];
    fits_get_img_size(*infile.fptr(), 2, naxes, &infile.status());
    infile.check();

    vector<double> hjd_buffer(naxes[0]), flux_buffer(naxes[0]);
    double period = 0., epoch = 0.;
    for (auto itr=models.begin(); itr!=models.end(); itr++) {
        auto object_name = sanitise_object_name(itr->name);
        /* Get the index of the original lightcurve */
        unsigned int SourceIndex = object_names[object_name];

        infile.moveHDU("HJD");
        fits_read_img(*infile.fptr(), TDOUBLE, (SourceIndex * naxes[0]) + 1, naxes[0],
                0, &hjd_buffer[0], 0, &infile.status());
        infile.check();

        infile.moveHDU("FLUX");
        fits_read_img(*infile.fptr(), TDOUBLE, (SourceIndex * naxes[0]) + 1, naxes[0], 0, &flux_buffer[0], 0, &infile.status());
        infile.check();

        double transit_width = WidthFromParams(*itr);

        if (valid_lightcurve(hjd_buffer, flux_buffer, transit_width, itr->period, itr->epoch)) {
            out.push_back(*itr);
        }
    }

    return out;
}

class RunCommand {
    public:
        RunCommand(const char *cmd) : f(NULL) {
            f = popen(cmd, "r");
            if (!f) {
                throw runtime_error("Error spawning subprocess");
            }
        }

        ~RunCommand() {
            int status_code = pclose(f);
            if (status_code != 0) {
                throw runtime_error("Error running subcommand");
            }
        }

    private:
        FILE *f;

};

// Run shell command
void exec(const char *cmd) {
    RunCommand rc(cmd);
}

void exec(const vector<string> &cmd) {
    stringstream ss;
    for (auto word: cmd) {
        ss << word << " ";
    }
    exec(ss.str().c_str());
}

vector<double> generate_model(const vector<double> &hjd, const Model &model) {
    NonlinearLimbDarkeningParameters ldc;
    ldc.c1 = model.c1;
    ldc.c2 = model.c2;
    ldc.c3 = model.c3;
    ldc.c4 = model.c4;

    // Stellar radius in metres
    double stellar_radius = model.rs * rSun;

    Params params;
    params.t0 = jd2wd(model.epoch);
    params.per = model.period * 86400.;
    params.rp = (model.rp * rJup) / stellar_radius;
    params.a = (model.a * AU) / stellar_radius;
    params.inc = model.i;
    params.ecc = 0.;
    params.w = 90.;
    params.ldc = ldc;

    double *pflux = light_curve(&params, &hjd[0], hjd.size());
    vector<double> flux(pflux, pflux + hjd.size());
    free(pflux);

#if 0
    const string hjd_filename = "hjd.txt";
    const string flux_filename = "flux.txt";

    FILE *outfile = fopen(hjd_filename.c_str(), "w");
    if (!outfile) {
        throw std::runtime_error("Cannot open hjd file for writing");
    }

    for (auto value: hjd) {
        fprintf(outfile, "%.10lf\n", value);
    }

    fclose(outfile);

    // Generate the model
    vector<string> cmd;
    cmd.push_back("/usr/local/python/bin/python");
    cmd.push_back("generate_model.py");
    cmd.push_back(hjd_filename);
    cmd.push_back("--model-name");
    cmd.push_back(model.name);
    cmd.push_back("--candidates");
    cmd.push_back(models_filename);
    cmd.push_back("-o");
    cmd.push_back(flux_filename);
    exec(cmd);

    string line;
    ifstream infile(flux_filename.c_str());

    if (!infile.is_open()) {
        throw std::runtime_error("Cannot read from flux file");
    }

    while (getline(infile, line)) {
        istringstream buffer(line);
        double value;
        buffer >> value;
        flux.push_back(value);
    }
#endif

    return flux;
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
        FetchesParameters param_fetcher(Config.DatabaseFilename);

        /* get the new models to insert */
        cout << "Fetching models from database file" << endl;
        vector<Model> models = param_fetcher.fetch_models();
        int nextra = models.size();

        if (nextra == 0) {
            throw runtime_error("No input models found");
        }

        /* Some lightcurves have no data, so compute the new number of extra
         * objects */
        cout << "Computing the list of valid lightcurves" << endl;
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

                    char *ColumnNames[] = {"REAL_OBJ_ID", "SKIPDET", "TRANS_PERIOD", "TRANS_WIDTH", "TRANS_DEPTH", "TRANS_EPOCH", "TRANS_RP", "TRANS_RS", "TRANS_A", "TRANS_I"};
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
		int obj_id_column_number = outfile.columnNumber("OBJ_ID");

        FalseColumnNumbers fcn;
        fcn.real_obj_id = outfile.columnNumber("REAL_OBJ_ID");
        fcn.skipdet = outfile.columnNumber("SKIPDET");
        fcn.period = outfile.columnNumber("TRANS_PERIOD");
        fcn.width = outfile.columnNumber("TRANS_WIDTH");
        fcn.depth = outfile.columnNumber("TRANS_DEPTH");
        fcn.epoch = outfile.columnNumber("TRANS_EPOCH");
        fcn.rp = outfile.columnNumber("TRANS_RP");
        fcn.rs = outfile.columnNumber("TRANS_RS");
        fcn.a = outfile.columnNumber("TRANS_A");
        fcn.i = outfile.columnNumber("TRANS_I");
        outfile.check();




        ts.start("model.iterate");
        const int NullSubIndex = -1;


        /* Get a list of the objects in the file */
        long nrows;
        map<string, unsigned int> ObjectNames = extract_object_names(infile, nrows);

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
            long SourceIndex = ObjectNames[sanitise_object_name(Current.name)];


            /* Copy the original data to the new location */
            outfile.moveHDU("HJD");
            long naxes[2];
            fits_get_img_size(*outfile.fptr(), 2, naxes, &outfile.status());
            outfile.check();




            /* And copy the other data parts two */
            vector<double> buffer(naxes[0]);

            /* Load the jd data separately */
            vector<int> jd(naxes[0]);
            outfile.moveHDU("HJD");
            fits_read_img(*outfile.fptr(), TINT, (SourceIndex * naxes[0]) + 1, naxes[0], 0, &jd[0], 0, &outfile.status());
            fits_write_img(*outfile.fptr(), TINT, (OutputIndex * naxes[0]) + 1, naxes[0], &jd[0], &outfile.status());
            outfile.moveHDU("FLUX");
            fits_read_img(*outfile.fptr(), TDOUBLE, (SourceIndex * naxes[0]) + 1, naxes[0], 0, &buffer[0], 0, &outfile.status());
            fits_write_img(*outfile.fptr(), TDOUBLE, (OutputIndex * naxes[0]) + 1, naxes[0], &buffer[0], &outfile.status());
            outfile.moveHDU("FLUX_ERR");
            fits_read_img(*outfile.fptr(), TDOUBLE, (SourceIndex * naxes[0]) + 1, naxes[0], 0, &buffer[0], 0, &outfile.status());
            fits_write_img(*outfile.fptr(), TDOUBLE, (OutputIndex * naxes[0]) + 1, naxes[0], &buffer[0], &outfile.status());
            outfile.moveHDU("CCDX");
            fits_read_img(*outfile.fptr(), TDOUBLE, (SourceIndex * naxes[0]) + 1, naxes[0], 0, &buffer[0], 0, &outfile.status());
            fits_write_img(*outfile.fptr(), TDOUBLE, (OutputIndex * naxes[0]) + 1, naxes[0], &buffer[0], &outfile.status());
            outfile.moveHDU("CCDY");
            fits_read_img(*outfile.fptr(), TDOUBLE, (SourceIndex * naxes[0]) + 1, naxes[0], 0, &buffer[0], 0, &outfile.status());
            fits_write_img(*outfile.fptr(), TDOUBLE, (OutputIndex * naxes[0]) + 1, naxes[0], &buffer[0], &outfile.status());
            outfile.moveHDU("SKYBKG");
            fits_read_img(*outfile.fptr(), TDOUBLE, (SourceIndex * naxes[0]) + 1, naxes[0], 0, &buffer[0], 0, &outfile.status());
            fits_write_img(*outfile.fptr(), TDOUBLE, (OutputIndex * naxes[0]) + 1, naxes[0], &buffer[0], &outfile.status());
            outfile.moveHDU("FLAGS");
            fits_read_img(*outfile.fptr(), TDOUBLE, (SourceIndex * naxes[0]) + 1, naxes[0], 0, &buffer[0], 0, &outfile.status());
            fits_write_img(*outfile.fptr(), TDOUBLE, (OutputIndex * naxes[0]) + 1, naxes[0], &buffer[0], &outfile.status());
            outfile.check();

            /* Add a transit model to the data */
            pair<double, long> LightcurveInfo = AlterLightcurveData(outfile, Config.DatabaseFilename, OutputIndex * naxes[0] + 1, naxes[0], Current, ArithMeth("+"), Config);


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
			char *current_name_cstr = const_cast<char*>(Current.name.c_str());
            string NewName = GenerateNewObjectName(counter);
            char *cstr = new char[7];
            strncpy(cstr, NewName.c_str(), 7);
            /* Write the real object id to REAL_OBJ_ID column, and the fake object
             * id to the OBJ_ID column */
            fits_write_col_str(*outfile.fptr(), fcn.real_obj_id, CatalogueIndex, 1, 1, &current_name_cstr, &outfile.status());
            outfile.check();
			fits_write_col_str(*outfile.fptr(), obj_id_column_number, CatalogueIndex, 1, 1,
							   &cstr, &outfile.status());
            outfile.check();
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
