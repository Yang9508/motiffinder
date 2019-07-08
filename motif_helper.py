#!/usr/bin/env python
# encoding: utf-8

import random
from random import uniform
from math import log
import numpy as np

THRESHOLD = 0.001

def f(a, b, c):
    d = 1 - a - b -c
    if c < 0 or c > 1 or d < 0 or d > 1:
        return -999
    res = 0
    if a != 0 : res += a*log(4*a,2)
    if b != 0: res += b*log(4*b,2)
    if c != 0: res += c*log(4*c, 2)
    if d != 0: res += d*log(4*d, 2)
    return res

def solve(a, b, ICPC):
    left, right = 1e-5, (1-a-b)/2
    fmax = f(a, b, left)
    fmin = f(a, b, right)
    if ICPC > fmax or ICPC < fmin:
        #print("no solution:", a, b, fmin, fmax)
        return -1,-1
    while True:
        mid = (left + right) / 2
        fmid = f(a, b, mid)
        if abs(fmid-ICPC) < THRESHOLD or abs(left-right) < 1e-3:
            break
        if ICPC > fmid:
            left = mid
        else:
            right = mid
    return mid, fmid

# Only call this method from other files
def get_motif_column(ICPC):
    if ICPC == 2:
        motif = [0] * 4
        motif[random.randint(0,3)] = 1
        return motif
    if ICPC == 0:
        return [0.25]*4
    while True:
        a = uniform(0, 1)
        b = uniform(0, 1)
        c, fc = solve(a, b, ICPC)
        if c >= 0 and abs(fc-ICPC)<THRESHOLD:
            return [a,b,c,1-a-b-c]
