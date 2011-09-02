#include <iostream>
#include <cassert>
#include <stdexcept>
#include <sqlitepp/sqlitepp.hpp>
#include <fitsio.h>
#include <tclap/CmdLine.h>
#include <sstream>

/* Local includes */
#include "timer.h"
#include "GenerateModel.h"
#include "Model.h"



using namespace std;
using namespace sqlitepp;


class Fits
{
    public:
        Fits(const string &filename)
            : m_status(0), m_filename(filename)
        {
            fits_open_file(&this->m_fptr, this->m_filename.c_str(), READWRITE, &this->m_status);
            this->check();
        }


        virtual ~Fits()
        {
            fits_close_file(this->m_fptr, &this->m_status);
            this->check();
        }

        void check()
        {
            if (this->m_status)
            {
                char buf[FLEN_STATUS];
                fits_get_errstatus(this->m_status, buf);
                throw runtime_error(buf);
            }
        }

        void moveHDU(const string &hduname)
        {
            fits_movnam_hdu(this->m_fptr, ANY_HDU, const_cast<char*>(hduname.c_str()), 0, &this->m_status);
            this->check();
        }

        void moveHDU(int hdunum)
        {
            int hdutype;
            fits_movabs_hdu(this->m_fptr, hdunum, &hdutype, &this->m_status);
            this->check();
        }

        fitsfile **fptr() { return &this->m_fptr; }
        int &status() { return this->m_status; }

        static void check(int status)
        {
            if (status)
            {
                char buf[FLEN_STATUS];
                fits_get_errstatus(status, buf);
                throw runtime_error(buf);
            }
        }

        const string hduname()
        {
            char buf[FLEN_VALUE];
            fits_read_key(this->m_fptr, TSTRING, "EXTNAME", buf, NULL, &this->m_status);
            this->check();
            return string(buf);
        }



    protected:
        fitsfile *m_fptr;
        int m_status;
        string m_filename;

        /* Default constructor - does nothing */
        Fits() {}


};

class NewFits : public Fits
{
    public:
    NewFits(const string &filename)
    {
        this->m_filename = filename;
        this->m_status = 0;

        fits_create_file(&this->m_fptr, filename.c_str(), &this->m_status);
        this->check();

        /* Ensure the basic keywords are there */
        long naxes[] = {0, 0};
        fits_create_img(this->m_fptr, BYTE_IMG, 0, naxes, &this->m_status);
        this->check();
    }
};

template <typename T>
void CopyImageData(Fits &infile, Fits &outfile, float MemLimit)
{
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
                    char *ColumnFormats[] = {"1I", "1D", "1D", "1D", "1J", "1D", "1D", "1D", "1D"};

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

        ts.start("model.iterate");
        Model Current;
        const int NullSubIndex = -1;

        /* Now iterate through every row adding a new lightcurve, and subtracting if necassary  */
        st << "select id, name, submodel_id, period, epoch, a, i, rs, rp, mstar, c1, c2, c3, c4, teff "
            " from addmodels", into(Current.id), into(Current.name), into(Current.submodel_id), 
            into(Current.period), into(Current.epoch), into(Current.a), into(Current.i), into(Current.rs),
            into(Current.rp), into(Current.mstar), into(Current.c1), into(Current.c2), into(Current.c3), 
            into(Current.c4), into(Current.teff);

        /* Get a list of the objects in the file */
        list<string> ObjectNames;
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

        while (st.exec())
        {
            if (Current.submodel_id != NullSubIndex)
            {
                //cout << "Submodel required" << endl;
            }

            /* Generate the add model */
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
