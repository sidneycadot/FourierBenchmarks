#! /usr/bin/env python3

import numpy as np

size = 16

k = np.arange(size)

z = np.cos(1 + k ** 2) + 1j * np.sin(1 + k **3)

z = np.fft.fft(z)

for (k, zk) in enumerate(z):
    print(k, zk)
