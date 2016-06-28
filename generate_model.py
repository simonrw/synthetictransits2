#!/usr/local/python/bin/python
# -*- coding: utf-8 -*-

from __future__ import absolute_import, division, print_function
import numpy as np
import batman
from astropy import units as u, constants as c
import argparse
import sqlite3

def dict_factory(cursor, row):
    d = {}
    for idx, col in enumerate(cursor.description):
        d[col[0]] = row[idx]
    return d


def jd2wd(jd):
    jd_ref = 2456658.500000
    return (jd - jd_ref) * 86400.


def build_params(model):
    stellar_radius = model['rs'] * u.R_sun
    params = batman.TransitParams()
    params.t0 = jd2wd(model['epoch'])
    params.per = model['period'] * 86400.
    params.rp = ((model['rp'] * c.R_jup).to(u.m) / stellar_radius.to(u.m)).to(
        u.dimensionless_unscaled).value
    params.a = (model['a'] * u.AU / stellar_radius).to(u.dimensionless_unscaled).value
    params.inc = model['i']
    params.ecc = 0.
    params.w = 90.
    params.limb_dark = 'nonlinear'
    params.u = [
        model[key] for key in ['c1', 'c2', 'c3', 'c4']
    ]

    return params


def main(args):
    hjd = jd2wd(np.loadtxt(args.hjd_filename))

    name = args.model_name
    with sqlite3.connect(args.candidates) as con:
        con.row_factory = dict_factory
        cur = con.cursor()
        cur.execute('''select * from addmodels where name = ?''', (name, ))
        rows = cur.fetchall()

    assert len(rows) == 1
    model_params = rows[0]

    params = build_params(model_params)
    model = batman.TransitModel(params, hjd)
    flux = model.light_curve(params)
    np.savetxt(args.output, flux)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('hjd_filename')
    parser.add_argument('-n', '--model-name', required=True)
    parser.add_argument('-c', '--candidates', required=True)
    parser.add_argument('-o', '--output', required=True)
    main(parser.parse_args())
