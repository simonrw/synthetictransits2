#include <iostream>
#include <stdexcept>
#include <sqlitepp/sqlitepp.hpp>
#include <fitsio.h>
#include <tclap/CmdLine.h>
#include "timer.h"



using namespace std;

class Fits
{
    public:
        Fits(const string &filename)
            : m_status(0), m_filename(filename)
        {
            fits_open_file(&this->m_fptr, this->m_filename.c_str(), READWRITE, &this->m_status);
            check();
        }


        virtual ~Fits()
        {
            fits_close_file(this->m_fptr, &this->m_status);
            check();
        }

        void check()
        {
            if (this->m_status)
            {
                throw runtime_error("Error with fitsio");
            }
        }

        void moveHDU(const string &hduname)
        {
            fits_movnam_hdu(this->m_fptr, ANY_HDU, const_cast<char*>(hduname.c_str()), 0, &this->m_status);
            check();
        }

        fitsfile **fptr() { return &this->m_fptr; }
        int &status() { return this->m_status; }

    protected:
        fitsfile *m_fptr;
        int m_status;
        string m_filename;


};

template <typename T>
void CopyImageData(Fits &infile, Fits &outfile, float MemLimit)
{
}



int main(int argc, char *argv[])
{
    Timer ts;
    ts.start("all");
    ts.stop("all");
    return 0;
}
