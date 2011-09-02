#include <iostream>
#include <stdexcept>
#include <sqlitepp/sqlitepp.hpp>
#include <fitsio.h>
#include <tclap/CmdLine.h>


using namespace std;

struct Fits
{
    fitsfile *fptr;
    int status;
    Fits(const string &filename)
        : status(0)
    {
        fits_open_file(&fptr, filename.c_str(), READWRITE, &status);
        check();
    }


    virtual ~Fits()
    {
        fits_close_file(fptr, &status);
        check();
    }

    void check()
    {
        if (status)
        {
            throw runtime_error("Error with fitsio");
        }
    }
};


int main(int argc, char *argv[])
{
    cout << "Hello world!" << endl;
    return 0;
}
