#ifndef FITSOBJECT_H

#define FITSOBJECT_H

#include <string>
#include <fitsio.h>

class Fits
{
    public:
    Fits(const std::string &filename);
    virtual ~Fits();
    void check();
    void moveHDU(const std::string &hduname);
    void moveHDU(int hdunum);
    fitsfile **fptr();
    int &status();
    static void check(int status);
    const std::string hduname();


    protected:
    fitsfile *m_fptr;
    int m_status;
    std::string m_filename;


    /* Default constructor - does nothing */
    Fits();
};


class NewFits : public Fits
{
    public:
    NewFits(const std::string &filename);
};



#endif /* end of include guard: FITSOBJECT_H */
