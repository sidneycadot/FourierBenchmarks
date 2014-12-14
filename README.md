
* Functionality offered

  - Number of transforms
  - Forward, inverse
  - Scaling
  - Real vs complex
  - Basic data type (float, double, long double)
  - Dimensionality of transform
  - Memory layout
  - In-place or not
  - Preserve input or not

# Quality

## Precision of calculations

## Speed of calculations

  - For different sizes, memory layouts, etc.
  - Cache
  - Threading
  - Optimization (planning) level (FFTW)
  - Different HW plaforms? (Intel/AMD, AMD/nVidia)

## Quality of implementation

  - Memory leaks?


# Implementations

## Overview

|------------|-----------------------------------------|----------------|--------|---------|--------------------------------------------|---------|
| name       | long name                               | CPU/GPU        | vendor | license | homepage                                   | version |
|------------|-----------------------------------------|----------------|--------|---------|--------------------------------------------|---------|
| MKL        | Intel Math Kernel Library               | CPU            | Intel  |         | https://software.intel.com/en-us/intel-mkl |         |
| IPP        | Intel Integrated Performance Primitives | CPU            | Intel  |         | https://software.intel.com/en-us/intel-ipp |         |
| FFTW       | Fastest Fourier Transfer in the West    | CPU            |        |         |                                            |         |
| FFTS       | Fastest Fourier Transfer in the South   | CPU            |        |         | https://github.com/anthonix/ffts           |         |
| KissFFT    | Kiss FFT                                | CPU            |        | BSD     | http://sourceforge.net/projects/kissfft/   |         |
| cuFFT      | nVidia CUDA FFT                         | GPU (nVidia)   | nVidia |         | https://developer.nvidia.com/cuFFT         |         |
|------------|-----------------------------------------|----------------|--------|---------|--------------------------------------------|---------|

## MKL

### API

## IPP

### API

## FFTW

### API

## FFTS

### API

## KissFFT

### API

## cuFFT

### API
