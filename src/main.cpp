#include <iostream>
#include <stdexcept>
#include <sqlitepp/sqlitepp.hpp>
#include <fitsio.h>
#include <tclap/CmdLine.h>
#include <sstream>
#include "timer.h"



using namespace std;

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
        Timer ts;
        ts.start("all");

        Fits infile("../data.fits");
        NewFits outfile("!output.fits");
        infile.moveHDU("FLUX");
        ts.stop("all");
        return 0;
    }
    catch (runtime_error &e)
    {
        cerr << e.what() << endl;
    }

    return EXIT_FAILURE;
}
