#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
sys.path.insert(0, '/home/astro/phrfbf/build/lib/python2.6/site-packages')
sys.path.insert(0, '/home/astro/phrebf/Python/JG')
import pyfits
from pylab import *
import argparse
from matplotlib.backends.backend_pdf import PdfPages
from jg.ctx import wd2jd

wasp12 = {'p': 1.0914222,
        'e': 2454508.97605
        }



rJup = 71492E3
rSun = 6.995E8
AU = 1.496E11
mSun = 1.9891E30
secondsInMinute = 60.
secondsInHour = 60. * secondsInMinute
secondsInDay = 24. * secondsInHour
radiansInDegree = 2. * 3.14 / 360.
degreesInRadian = 360. / 2. / 3.14
GravConst = 6.67E-11
tSun = 5778.


def PJWMethod(m):
    ''' 
    Returns the transit width based on Pete's lecture notes
    '''
    Norm = m['period'] / pi
    FirstTerm = ((m['rp'] + m['rs']) / m['a'])**2
    SecondTerm = (cos(m['i'])**2)

    InsideSqrt = FirstTerm - SecondTerm
    Width = Norm * arcsin(sqrt(InsideSqrt))
    return Width

def JWMethod(m):
    '''
    Returns the transit width based on Joshua Winn's method
    '''
    Norm = m['period'] / pi
    Norm2 = m['rs'] / m['a']
    k = m['rp'] / m['rs']
    b = m['a'] * cos(m['i']) / m['rs']
    InsideSqrt = ((1+k)**2 - b**2) / sin(m['i'])
    Width = Norm * arcsin(Norm2 * sqrt(InsideSqrt))
    return Width


def main(args):
    f  = pyfits.open(args.file)

    Catalogue = f['catalogue'].data
    #Flux = f['flux'].data
    #HJD = f['hjd'].data

    ## get the objects where the width is not 0
    Widths = Catalogue.field("FAKE_WIDTH")
    Epochs = Catalogue.field("FAKE_EPOCH")
    Periods = Catalogue.field("FAKE_PERIOD")
    Depths = Catalogue.field("FAKE_DEPTH")
    RPlanets = Catalogue.field('FAKE_RP')
    RStars = Catalogue.field('FAKE_RS')
    Inclinations = Catalogue.field("FAKE_I")
    Separations = Catalogue.field("FAKE_A")
    Names = Catalogue.field('OBJ_ID')
    #Index, = where(Widths!=0)

    Index = []
    for i in range(Names.size):
        if "SYNTH" in Names[i]:
            Index.append(i)



    Index = array(Index)

    pp = PdfPages("output.pdf")
    pp2 = PdfPages('wasp12phase.pdf')

    ##ion()
    Reversed = Index[::-1]

    nObjects = 20
    skipStep = 100
    for i in Reversed[:nObjects*skipStep:skipStep]:

        CurrentModel = {'epoch': Epochs[i],
                'period': Periods[i],
                'rp': RPlanets[i],
                'rs': RStars[i],
                'i': Inclinations[i],
                'a': Separations[i],
                }

        print i, CurrentModel

        cla()
        Time = wd2jd(f['hjd'].section[i])
        Lightcurve = f['flux'].section[i]

        # remove nans
        GoodIndex = (Lightcurve==Lightcurve) & (Time==Time)
        Lightcurve = Lightcurve[GoodIndex]
        Time = Time[GoodIndex]

        Phase = ((Time - CurrentModel['epoch']) / (CurrentModel['period'] / secondsInDay)) % 1.0
        Phase[Phase>0.5] -= 1.0

        Lightcurve /= median(Lightcurve)

        plot(Phase, Lightcurve, 'r,')
        title("Depth: %f" % (Depths[i],))

        axhline(1. - Depths[i])
        axvline(-Widths[i]/2. / CurrentModel['period'], color='k')
        axvline(Widths[i]/2. / CurrentModel['period'], color='k', label='Original')


        xlim(-0.3, 0.3)
        ylim(0.5, 1.5)

        pp.savefig()

        cla()
        Phase = ((Time - wasp12['e']) / wasp12['p']) % 1.0
        Phase[Phase>0.5] -= 1.0

        plot(Phase, Lightcurve, 'r,')
        title("WASP-12b phase")
        xlim(-0.3, 0.3)
        pp2.savefig()

    pp.close()
    pp2.close()

    #print Epochs[Index]


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("file")
    args = parser.parse_args()

    main(args)

