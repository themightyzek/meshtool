#include <iostream>
#include <vector>
#include <string>
#include "CLI11.hpp"

using namespace std;

int main(int argc, char** argv)
{
    cout << "You have entered " << argc
         << " arguments:" << "\n";
  
    for (int i = 0; i < argc; ++i)
        cout << argv[i] << "\n";

    CLI::App app{"App description"};
    string filename = "default";
    app.add_option("-i,--input", filename, "The input file");

    CLI11_PARSE(app, argc, argv);

    cout << "filename is: " << filename << endl;

    return 0;
}