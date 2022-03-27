#include <iostream>
#include <vector>
#include <string>
#include "CLI11.hpp"
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/convex_hull_2.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point_2;

using namespace std;

int main(int argc, char **argv)
{
    cout << "You have entered " << argc
         << " arguments:"
         << "\n";

    for (int i = 0; i < argc; ++i)
        cout << argv[i] << "\n";

    CLI::App app{"App description"};
    string filename = "default";
    app.add_option("-i,--input", filename, "The input file");

    CLI11_PARSE(app, argc, argv);

    cout << "filename is: " << filename << endl;

    Point_2 points[5] = {Point_2(0, 0), Point_2(10, 0), Point_2(10, 10), Point_2(6, 5), Point_2(4, 1)};
    Point_2 result[5];
    Point_2 *ptr = CGAL::convex_hull_2(points, points + 5, result);
    std::cout << ptr - result << " points on the convex hull:" << std::endl;
    for (int i = 0; i < ptr - result; i++)
    {
        std::cout << result[i] << std::endl;
    }

    return 0;
}