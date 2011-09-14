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
    Widths = Catalogue.field("FAKE_WIDTH")
    Index, = where(Widths!=0)

    pp = PdfPages("output.pdf")
    pp2 = PdfPages('wasp12phase.pdf')

    ##ion()
    Reversed = Index[::-1]
    for i in Reversed[:20]:
        print i
        cla()
        Time = wd2jd(f['hjd'].section[i])
        Epoch = Epochs[i]
        Period = Periods[i]

        Phase = ((Time - Epoch) / (Period / 86400.)) % 1.0
        Phase[Phase>0.5] -= 1.0
        Lightcurve = f['flux'].section[i]

        Lightcurve /= median(Lightcurve)

        plot(Phase, Lightcurve, 'r,')
        title("Depth: %f, width: %f" % (Depths[i], Widths[i] / Period))

        axhline(1. - Depths[i])
        axvline(-Widths[i] / Period)
        axvline(Widths[i] / Period)

        xlim(-0.3, 0.3)

        pp.savefig()

    pp.close()
    pp2.close()

    #print Epochs[Index]


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("file")
    args = parser.parse_args()

    main(args)

