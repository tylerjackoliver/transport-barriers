# Transport Barriers

[![Contributors][contributors-shield]][contributors-url]
[![Issues](https://img.shields.io/github/issues/tylerjackoliver/transport-barriers)](https://github.com/tylerjackoliver/transport-barriers/issues)
![License](https://img.shields.io/github/license/tylerjackoliver/transport-barriers)


<!-- PROJECT LOGO -->
<br />
<p align="center">
  <h3 align="center">Transport Barriers: Find LCS in generic dynamical systems</h3>

  <p align="center">
    This program identifies reduced strainlines and stretchlines in generic (mostly astrodynamic) systems, as per <a href="https://doi.org/10.1016/j.physd.2014.01.007">Blazevski & Haller, 2014</a>.
    <br />
    <a href="https://github.com/tylerjackoliver/transport-barriers/issues">Report Bug</a>
  </p>
</p>



<!-- TABLE OF CONTENTS -->
## Table of Contents

* [About the Project](#about-the-project)
  * [Built With](#built-with)
* [Getting Started](#getting-started)
  * [Prerequisites](#prerequisites)
  * [Installation](#building)
* [Contact](#contact)



<!-- ABOUT THE PROJECT -->
## About The Project

As part of ongoing research into the use of Lagrangian Coherent Structures (LCSs) in Astrodynamics applications, there is a requirement for a generic tool for finding LCSs in dynamical systems. This program will take a three-dimensional arbitrary dynamical system, defined on some dense grid, and determine the reduced strainlines and stretchlines for the flow.

### Built With

* [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page) for obtaining eigenvalues and eigendirections of the right Cauchy-Green strain tensor.
* [BOOST](https://www.boost.org/) for data containers and numerical integration.
* [OpenMP](https://www.openmp.org/) for shared memory parallelism.
* [Intel MKL](https://software.intel.com/content/www/us/en/develop/tools/math-kernel-library.html) for accelerating Eigen linear algebra operations.
* [Intel MPI](https://software.intel.com/content/www/us/en/develop/tools/mpi-library.html), although any MPI implementation that supports the Intel C++ compiler will do.
* [CMake](https://cmake.org/) for building the software.
* [SPLINTER](https://github.com/bgrimstad/splinter), for Helicity field interpolation.

<!-- GETTING STARTED -->
## Getting Started

### Prerequisites

The dependencies given above are required prior to build. The program is currently built and tested with the Intel C++ compiler for better compiler-directed vectorisation and optimization support. GNU compilers have been tested to work, but are not officially supported.

Any compiler used must support the C++14 standard.
### Building

0. Install any required dependencies

1. Clone the repo
```sh
git clone https://github.com/tylerjackoliver/transport-barriers.git
```
2. Run the CMake wrapper
```sh
sh build.sh
```

### Running
TBW.

<!-- CONTACT -->
## Contact

Jack Tyler - [@tylerjackoliver](https://twitter.com/tylerjackoliver) - jack.tyler@soton.ac.uk

Project Link: [https://github.com/tylerjackoliver/ACROBAT](https://github.com/tylerjackoliver/transport-barriers)

<!-- MARKDOWN LINKS & IMAGES -->
<!-- https://www.markdownguide.org/basic-syntax/#reference-style-links -->
[contributors-shield]: https://img.shields.io/github/contributors/tylerjackoliver/transport-barriers
[contributors-url]: https://img.shields.io/github/contributors/tylerjackoliver/transport-barriers
[issues-shield]: https://github.com/tylerjackoliver/transport-barriers/issues
[issues-url]: https://img.shields.io/github/issues/tylerjackoliver/transport-barriers
[license-shield]: ""
[license-url]: https://img.shields.io/github/license/tylerjackoliver/transport-barriers
[product-screenshot]: images/screenshot.png
