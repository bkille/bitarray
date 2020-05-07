#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 31 01:08:05 2020

@author: bkille
"""

from bitarray import bitarray, bits2bytes
import time
from random import randint
import os
from tqdm import tqdm

def randombitarray(start=1):
    for n in list(range(start, 25)) + [randint(1000, 2000)]:
        a = bitarray(endian='little')
        a.frombytes(os.urandom(bits2bytes(n)))
        del a[n:]
        yield a
        
        

n = 2**31
N = 2
a = bitarray("1"*(n//2) + "0"*(n//2))
# a.frombytes(os.urandom(bits2bytes(n)))
t0 = time.time()
for i in range(N):
    a.reverse()
total_time = time.time() - t0
print(f"Reversed in {total_time/N}")
t0 = time.time()
for i in range(N):
    a.index(0)
total_time = time.time() - t0
print(f"Indexed in {total_time/ N}")

t0 = time.time()
for i in range(N):
    a.count(0)
total_time = time.time() - t0
print(f"Counted in {total_time/N}")

t0 = time.time()
for i in range(N):
    a.search(bitarray("11111111111111111111111111111111111111111111111111111111111100000000000"))
total_time = time.time() - t0
print(f"Searched in {total_time/N}")


t0 = time.time()
for i in range(N):
    b = a & a
total_time = time.time() - t0
print(f"Bitwise OP in {total_time/N}")


