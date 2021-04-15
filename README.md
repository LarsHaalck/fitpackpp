fitpackpp: C++ B-Spline Smoothing using FITPACK 
=================================================

Fitpackpp supplies some C++-bindings of Fortran subroutines from [FITPACK](http://www.netlib.org/dierckx) as used in [SciPy](https://www.scipy.org).

The class BSplineCurve implements the functions [splrep](https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.splrep.html) (via the Constructor of `BSplineCurve`) and [splev](https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.splev.html#scipy.interpolate.splev) (via `operator()`) as described on their respective pages in the SciPy-Documentation.

Usage
=================================================
```
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ../
make
make install
```

Then use the library in your CMake-project with
```
find_package(fitpackpp)
target_link_libraries(tgt fitpackpp::fitpackpp)
```

Todos:
=================================================
* [x] implement splrep
* [x] implement slpev
* [ ] implement splreps `task = 1`
* [ ] implement splder
* [ ] implement spalde
* [ ] implement splantider
* [ ] implement sproot
* [ ] implement subroutines for n-d functions and surfaces

