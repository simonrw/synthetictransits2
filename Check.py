#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pyfits
from pylab import *
import argparse

def main(args):
    f  = pyfits.open(args.file)

    Catalogue = f['catalogue'].data
    Flux = f['flux'].data
    HJD = f['hjd'].data

    # get the objects where the width is not 0
    Widths = Catalogue.field("FAKE_WIDTH")
    Periods = Catalogue.field("FAKE_PERIOD")
    Index, = where(Widths!=0)

    ion()
    for i in Index[::-1]:
        cla()
        Time = HJD[i]
        Lightcurve = Flux[i]
        Period = (Periods[i] / 86400.)
        plot((Time/Period) % 1.0, Lightcurve, 'r,')
        draw()

        raw_input("Press enter to continue")



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("file")
    args = parser.parse_args()

    main(args)

