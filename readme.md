Required libs:
- CGAL v5.0.4 (apt-get install libcgal-dev)
- CLI11 v2.2.0 (download https://github.com/CLIUtils/CLI11/releases/download/v2.2.0/CLI11.hpp into /usr/include)
- Eigen v3.3.7 (download https://gitlab.com/libeigen/eigen/-/releases/3.3.7 then copy "Eigen" folder into /usr/include/Eigen)

rlk/obj is statically linked (sources copied to this repository)

compiler:
- gcc 9.4.0

OS: 
- ubuntu 20.04

`/meshtool $ cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_STANDARD=20 CMakeLists.txt`

`/meshtool $ make`


everything else should be included in CGAL. If you want to use this project academically or for any other reason, and you need help with it, tweet or DM me on twitter @themightyzek. 
