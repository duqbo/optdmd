# optdmd 

A MATLAB package for computing the optimized DMD

## About

Information about the optimized DMD and this implementation
can be found in [this preprint](https://arxiv.org/abs/1704.02343v1).
The source files contain documentation; for 
an interactive look at this documentation, open 
the index.html file from the docs folder in
a web browser.

## Set-up

For the most efficient version of this code, you
must compile a mex binary. It can be difficult to 
avoid platform-dependent issues for mex binaries,
so please bear with us. 

The first thing to try is to Open MatLab and run 
setup.m to install. If the test runs and returns
a small error, then the install was likely successful.

If setup.m fails and you have experience with
mex files, consider editing the compile flags in
src/buildqrmex

Please create an issue using the issues tab on
the [optdmd GitHub page](https://github.com/duqbo/optdmd)
if you are unable to get these files compiled on your
system.

If everything has failed, you can still run this code
but with slower, less memory efficient routines. 
To do this, copy varpro2.m.slow into varpro2.m in the
src folder (the original fast version is saved in 
varpro2.m.fast).

## How to use

If you'd like to see how to use the OPTDMD wrappers
(for computing the optimized dynamic mode decomposition)
the best place to start is to check out simple_example.m
(in "examples" folder)

## Updates

Feel free to submit bug-fix/feature requests through
the issues tab on GitHub.

## Citing this software

To cite this software, follow the Zenodo link here:
[![DOI](https://zenodo.org/badge/86738651.svg)](https://zenodo.org/badge/latestdoi/86738651)


## License 

The files in the "src" and "examples" directories are available under the MIT license unless noted otherwise (see license* files in src directory).

The MIT License (MIT)

Copyright (c) 2017 Travis Askham

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.