* Cosmic ray particle generator written in C++
* Uses PARMA Fortran library ([https://phits.jaea.go.jp/expacs/]).
* It generates single cosmic rays (type, energy, momentum direction).
* Input settings are: wanted type, altitude, and energy range. Can be changed in `src/src/Settings.hh` .
* Tested on Linux only (Ubuntu 18.04)
* To build: open a terminal in the folder `build/` and run  `cmake ../` and then `make`
* To run : execute `./CR_GE` 

