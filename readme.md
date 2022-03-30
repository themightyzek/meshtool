install CGAL (apt-get install libcgal-dev)
install CLI11 (download https://github.com/CLIUtils/CLI11/releases/download/v2.2.0/CLI11.hpp into /usr/include)

/meshtool $ cgal_create_CMakeLists
/meshtool $ cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_STANDARD=20 CMakeLists.txt
/meshtool $ make
