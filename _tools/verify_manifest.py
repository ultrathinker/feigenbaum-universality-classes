# -*- coding: utf-8 -*-
"""Verify MANIFEST.sha256 of the package: python _tools/verify_manifest.py  (run from package root or anywhere)."""
import hashlib, os, sys
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
bad = 0; n = 0
for line in open(os.path.join(ROOT, 'MANIFEST.sha256'), encoding='utf-8'):
    line = line.rstrip('\n')
    if not line.strip(): continue
    h, rel = line.split('  ', 1)
    p = os.path.join(ROOT, rel)
    n += 1
    if not os.path.exists(p):
        print('MISSING', rel); bad += 1; continue
    if hashlib.sha256(open(p, 'rb').read()).hexdigest() != h:
        print('MISMATCH', rel); bad += 1
print('checked %d files, %d problems' % (n, bad))
sys.exit(1 if bad else 0)
