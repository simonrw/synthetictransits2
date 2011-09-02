#include <iostream>
#include <stdexcept>
#include <sqlitepp/sqlitepp.hpp>
#include <fitsio.h>
#include <tclap/CmdLine.h>
#include <sstream>
#include "timer.h"



using namespace std;

struct Model
{
    int id;
    string name;
    int submodel_id;
    double period;
    double epoch;
    double a;
    double i;
    double rs;
    double rp;
    double ms;
    double c1, c2, c3, c4;
    double teff;
};

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

        /* Start by getting file information from the input */
        int nhdus = 0;
        fits_get_num_hdus(*infile.fptr(), &nhdus, &infile.status());
        infile.check();

        cout << nhdus << " hdus found" << endl;


        /* Open the sqlite3 database here */
        sqlitepp::session conn(candidates_arg.getValue());
        sqlitepp::statement st(conn);

        /* get the required number of new objects */
        int nextra = 0;
        st << "select count(*) from addmodels", sqlitepp::into(nextra);
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
                    /* Append extra rows */

                    fits_insert_rows(*outfile.fptr(), nrows, nextra, &outfile.status());
                    outfile.check();

                    cout << " -> " << nrows + nextra << " rows";


                }
            }




            Fits::check(status);
            cout << endl;

        }


        ts.stop("copy");

        /* File copy finished */

        /* Now iterate through every row adding a new lightcurve, and subtracting if necassary  */



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
