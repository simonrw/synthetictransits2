#include <iostream>
#include <batman.h>
#include <cmath>
#include <algorithm>
#include <sstream>
#include <iomanip>
#include <map>
#include "constants.h"
#include "fetches_parameters.h"
#include "GenerateModel.h"
#include "FitsObject.h"
#include "GenerateModel.h"

using namespace std;

void analyse_model(const Model &model, const vector<double> &hjd);
vector<Model> compute_valid_extra_models(const vector<Model> &models, ReadOnlyFits &infile);
map<string, unsigned int> extract_object_names(ReadOnlyFits &infile, long &nrows);
string zeroPad(const string &s, int width);
string sanitise_object_name(const string &name);
bool valid_lightcurve(const vector<double> &flux);

vector<double> read_hjd(ReadOnlyFits &infile);
vector<double> compute_batman_model(const Model&, const vector<double>&);
vector<double> compute_srw_model(const Model&, const vector<double>&);
vector<double> generate_model(const vector<double> &hjd, const Model &model);

const double jd_ref = 2456658.500000;
double jd2wd(double jd) {
    return (jd - jd_ref) * secondsInDay;
}

int main() {
    cout << "Fetching parameters" << endl;
    FetchesParameters param_fetcher("/Volumes/External/synthetic-transits/MODELS_NG0522-2518_802_2016_TEST16.db");
    auto models = param_fetcher.fetch_models();

    cout << "Rejecting invalid models" << endl;
    ReadOnlyFits infile("/Volumes/External/synthetic-transits/NG0522-2518.fits");
    vector<double> hjd = read_hjd(infile);
    vector<Model> valid_models = compute_valid_extra_models(models, infile);

    cout << "Iterating" << endl;
    int included_models = 0;
    for (auto model: valid_models) {
        if (included_models > 10) {
            break;
        }

        auto transit_depth = pow((model.rp * rJup) / (model.rs * rSun), 2);
        if ((transit_depth > 0.02) && (model.period < 7)) {
            analyse_model(model, hjd);
            included_models += 1;
        }
    }
}

vector<double> read_hjd(ReadOnlyFits &infile) {
    infile.moveHDU("IMAGELIST");
    infile.check();

    int ncols = 0;
    fits_get_num_cols(*infile.fptr(), &ncols, &infile.status());
    infile.check();

    int colno = infile.columnNumber("TMID");
    vector<double> out(ncols);
    fits_read_col(*infile.fptr(), TDOUBLE, colno, 1, 1, ncols, NULL, &out[0], NULL, &infile.status());
    infile.check();

    return out;
}

vector<double> compute_batman_model(const Model &model, const vector<double> &hjd) {
    return generate_model(hjd, model);
}

vector<double> compute_srw_model(const Model &model, const vector<double> &hjd) {
    return GenerateSynthetic(hjd, model);
}


void analyse_model(const Model &model, const vector<double> &hjd) {
    cout << "Including model " << model.name << endl;

    vector<double> batman_model = compute_batman_model(model, hjd);
    vector<double> srw_model = compute_srw_model(model, hjd);

    stringstream ss;
    ss << "lightcurves/" << model.name << ".txt";
    auto filename = ss.str();

    FILE *outfile = fopen(filename.c_str(), "w");
    for (int i=0; i<hjd.size(); i++) {
        fprintf(outfile, "%lf %.5lf %.5lf\n", hjd[i], batman_model[i], srw_model[i]);
    }
    fclose(outfile);

}

vector<Model> compute_valid_extra_models(const vector<Model> &models, ReadOnlyFits &infile) {
    vector<Model> out;

    long nrows;
    map<string, unsigned int> object_names = extract_object_names(infile, nrows);

    infile.moveHDU("FLUX");
    long naxes[2];
    fits_get_img_size(*infile.fptr(), 2, naxes, &infile.status());
    infile.check();

    vector<double> buffer(naxes[0]);
    cout << "Computing the list of valid lightcurves" << endl;
    for (auto itr=models.begin(); itr!=models.end(); itr++) {
        auto object_name = sanitise_object_name(itr->name);
        /* Get the index of the original lightcurve */
        unsigned int SourceIndex = object_names[object_name];

        fits_read_img(*infile.fptr(), TDOUBLE, (SourceIndex * naxes[0]) + 1, naxes[0], 0, &buffer[0], 0, &infile.status());
        infile.check();

        if (valid_lightcurve(buffer)) {
            out.push_back(*itr);
        }
    }

    return out;
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

string zeroPad(const string &s, int width) {
    stringstream ss;
    ss << setw(width) << setfill('0') << s;
    return ss.str();
}

// Return the 6 character zero padded name
string sanitise_object_name(const string &name) {
    return zeroPad(name, 6);
}

bool valid_lightcurve(const vector<double> &flux) {
    int nvalid_points = std::count_if(flux.begin(), flux.end(), [](const double f) {
                return f > 0.0;
            });
    return (float(nvalid_points) / float(flux.size())) >= 0.7;
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
