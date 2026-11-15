# helper to record timing appendix
import time
import os
import sys

def append(label, t):
    with open('timing_steps.txt', 'a', encoding='utf-8') as f:
        f.write('%s\t%.1f\t%s\n' % (time.strftime('%H:%M:%S'), t, label))
    print('timing: %s %.1fs' % (label, t))

if __name__ == '__main__':
    label = sys.argv[1]
    t = float(sys.argv[2])
    append(label, t)
