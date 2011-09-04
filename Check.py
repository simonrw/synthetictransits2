#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pyfits
from pylab import *
import argparse
from matplotlib.backends.backend_pdf import PdfPages


def main(args):
    f  = pyfits.open(args.file)

    Catalogue = f['catalogue'].data
    Flux = f['flux'].data
    HJD = f['hjd'].data

    # get the objects where the width is not 0
    Widths = Catalogue.field("FAKE_WIDTH")
    Epochs = Catalogue.field("FAKE_EPOCH")
    Periods = Catalogue.field("FAKE_PERIOD")
    Depths = Catalogue.field("FAKE_DEPTH")
    Widths = Catalogue.field("FAKE_WIDTH")
    Index, = where(Widths!=0)

    pp = PdfPages("output.pdf")

    #ion()
    Reversed = Index[::-1]
    for i in Reversed[:20]:
        print i
        cla()
        Time = HJD[i]
        Epoch = Epochs[i]
        Period = Periods[i]

        Phase = ((Time - Epoch) / (Period / 86400.)) % 1.0
        Phase[Phase>0.5] -= 1.0
        Lightcurve = Flux[i]

        plot(Phase, Lightcurve, 'r,')
        title("Depth: %f, width: %f" % (Depths[i], Widths[i] / Period))

        xlim(-0.3, 0.3)

        pp.savefig()

    pp.close()



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("file")
    args = parser.parse_args()

    main(args)

