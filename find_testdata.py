#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import absolute_import, division, print_function
import os
import socket

hostname = socket.gethostname()
if 'ngts10' in hostname.lower():
    print('/local/srw/synthetic-testdata')
elif 'ngtshead' in hostname.lower():
    print(os.path.realpath(
        os.path.join(os.path.dirname(__file__), '..', 'testdata')))
else:
    # Laptop
    if os.path.isdir('/Volumes/Extra'):
        print('/Volumes/Extra/NGTS/synthetic-transits/testdata')
    elif os.path.isdir('/Volumes/External'):
        raise NotImplementedError
    else:
        raise OSError('Cannot find test data directory')
