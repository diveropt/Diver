Diver
=====

A fast parameter sampler and optimiser based on differential evolution.


About
--

Diver is written in Fortran2008 and uses MPI.  We originally wrote it with applications in particle physics and astronomy in mind, but it has since found broad applicability in other fields.

The code and its options are documented in the project's [GitHub wiki](https://github.com/diveropt/Diver/wiki), which has been updated from the description in the original release paper (included in the distribution as ScannerBit.pdf). Any papers that use results or insights obtained with Diver should cite this paper:
  1. Martinez, McKay, Farmer, Scott, Roebber, Putze & Conrad 2017, European Physical Journal C 77 (2017) 761, [arXiv:1705.07959](https://arxiv.org/abs/1705.07959).
One can also find detailed performance comparisons of Diver with other samplers and optimisers in both the above paper and in:
  2. DarkMachines High Dimensional Sampling Group, JHEP 05 (2021) 108, [arXiv:2101.04525](https://arxiv.org/abs/2101.04525).
In particular, this latter paper demonstrates that Diver significantly outperforms Scipy's implementation of differential evolution.

Compilation
--

The Diver build system is not really complex enough to require autotools or cmake.  Just change [the makefile](makefile) by hand to suit your system, or call it from another makefile or the command line with overrides.

To build Diver as a static library, and build all examples, do
```
  make
```

To instead build Diver as a shared library, do
```
  make libdiver.so
```

To build only the static library, do
```
  make libdiver.a
```

Testing
--

To build only the Fortran, C and C++ examples, do
```
  make example_f
  make example_c
  make example_cpp
```

The executables will appear in the Diver root directory. (You should run them.)


Licensing
--

The actual license is below.  In simple terms: Diver is free for academic use (no need to ask us), potentially free for other non-profit use (but you need to ask us explicitly), and not free for commercial use (you'll likely need to pay for a license in that case).

Contact: diver.optimisation@gmail.com


License
--

Copyright 2013-2023 Elinore Roebber and Pat Scott

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Use of the software must occur exclusively in an academic context.

2. Publications that include results or insights obtained using the software must cite the literature listed above.

3. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

4. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

5. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

Use of the software outside academic contexts requires contacting the copyright holders to make alternative licensing arrangements.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
