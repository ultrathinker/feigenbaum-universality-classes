#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# extract.py -- tiny helper: print lines of a log matching any keyword. [oid:deg6B] NOT RIGOROUS
import sys

data = open(sys.argv[1], 'rb').read().decode('utf-8', errors='replace')
keys = sys.argv[2:]
for line in data.splitlines():
    if any(k in line for k in keys):
        print(line.strip())
