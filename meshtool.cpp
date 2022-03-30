#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <utility>
#include <map>

#include "CLI11.hpp"

#include "lodepng/png.h"

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh/Surface_mesh.h>
#include <CGAL/property_map.h>
#include <CGAL/IO/read_ply_points.h>
#include <CGAL/Surface_mesh_parameterization/parameterize.h>
#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/Triangulation_2.h>
#include <CGAL/boost/graph/iterator.h>

typedef CGAL::Simple_cartesian<double> K;
typedef K::FT FT;
typedef K::Point_3 Point;
typedef K::Point_2 Point_2;
typedef K::Vector_3 Vector;
typedef std::array<unsigned char, 3> Color;

// Point with normal, color and intensity
typedef std::tuple<Point, Vector, Color, int> PNCI;
typedef CGAL::Nth_of_tuple_property_map<0, PNCI> Point_map;
typedef CGAL::Nth_of_tuple_property_map<1, PNCI> Normal_map;
typedef CGAL::Nth_of_tuple_property_map<2, PNCI> Color_map;
typedef CGAL::Nth_of_tuple_property_map<3, PNCI> Intensity_map;

typedef CGAL::Surface_mesh<Point> SurfaceMesh;
typedef CGAL::Triangulation_2<K> Triangulation;

typedef boost::graph_traits<SurfaceMesh>::vertex_descriptor vertex_descriptor;
typedef boost::graph_traits<SurfaceMesh>::halfedge_descriptor halfedge_descriptor;
typedef boost::graph_traits<SurfaceMesh>::face_descriptor face_descriptor;

using namespace std;
namespace SMP = CGAL::Surface_mesh_parameterization;

int main(int argc, char **argv)
{
    CLI::App app{"App description"};
    string meshFilename = "default";
    string cloudFilename = "default";
    app.add_option("-m,--in-mesh", meshFilename, "The input mesh in .ply format");
    app.add_option("-c,--in-cloud", cloudFilename, "The input cloud in .ply format");
    CLI11_PARSE(app, argc, argv);

    // --------------CGAL----------------

    // read point cloud
    vector<PNCI> points;
    ifstream in_c(cloudFilename);
    if (!in_c ||
        !CGAL::read_ply_points_with_properties(
            in_c,
            back_inserter(points),
            CGAL::make_ply_point_reader(Point_map()),
            make_tuple(Color_map(),
                       CGAL::Construct_array(),
                       CGAL::PLY_property<unsigned char>("red"),
                       CGAL::PLY_property<unsigned char>("green"),
                       CGAL::PLY_property<unsigned char>("blue")),
            CGAL::make_ply_normal_reader(Normal_map())))
    {
        cerr << "Error: unable to read from file " << cloudFilename << endl;
        return EXIT_FAILURE;
    }

    // read mesh
    SurfaceMesh mesh;
    ifstream in_m(meshFilename);
    if (!in_m ||
        !CGAL::read_ply(
            in_m,
            mesh))
    {
        cerr << "Error: unable to read from file " << meshFilename << endl;
        return EXIT_FAILURE;
    }
    else
    {
        cout << "Success: loaded mesh " << meshFilename << endl;
        cout << mesh.number_of_vertices() << " verts, " << mesh.number_of_faces() << " tris" << endl;
        cout << "Vertex list:" << endl;
        for(vertex_descriptor v : mesh.vertices())
        {
            cout << v << ": " << mesh.point(v) << endl;
        }

        cout << "Face list:" << endl;
        for(face_descriptor face_i : mesh.faces())
        {
            cout << face_i << ": ";
            for(vertex_descriptor v : mesh.vertices_around_face(mesh.halfedge(face_i)))
                cout << v << " ";
            cout << endl;
        }
    }

    // UV unwrap

    halfedge_descriptor bhd = CGAL::Polygon_mesh_processing::longest_border(mesh).first;

    typedef SurfaceMesh::Property_map<vertex_descriptor, Point_2> UV_pmap;
    UV_pmap uv_map = mesh.add_property_map<vertex_descriptor, Point_2>("h:uv").first;
    SMP::parameterize(mesh, bhd, uv_map);

    cout << "Success: created UV map for mesh" << endl << "UV Coords: " << endl;
    for(vertex_descriptor v : mesh.vertices())
    {
        cout << v << ": " << uv_map[v] << endl;
    }

    // create triangulation from UV map
    Triangulation Tr;
    map<vertex_descriptor, Triangulation::Vertex_handle> mesh2uv_tr;

    for(vertex_descriptor v : mesh.vertices())
    {
        mesh2uv_tr[v] = Tr.insert(uv_map[v]);
    }

    cout << "Success: Created Triangulation from UV Map" << endl;
    cout << Tr.number_of_vertices() << " vertices, " << Tr.number_of_faces() << " faces" << endl;
    cout << "Point list:" << endl; 

    for(Point_2 p : Tr.points())
        cout << p << endl;

    vector<Triangulation::Face_handle> face_handles;
    cout << "Face handles:" << endl;
    for(Triangulation::Face_handle f : Tr.finite_face_handles())
    {
        Triangulation::Face face = *f;
        cout << face << ": " << *face.vertex(0) << " | " << *face.vertex(1) << " | " << *face.vertex(2) << endl;
        face_handles.push_back(f);
    }

    // Triangulation has correct vertices, but its faces are not correct
    for(Triangulation::Face_handle f : face_handles)
    {
        Tr.delete_face(f);
    }

    cout << "Deleted Triangulation faces" << endl; 

    for(face_descriptor face_d : mesh.faces())
    {
        vector<vertex_descriptor> verts;
        for(vertex_descriptor vert : mesh.vertices_around_face(mesh.halfedge(face_d)))
            verts.push_back(vert);
    
        if(verts.size() == 3)
        {
            Triangulation::Face_handle newface = Tr.create_face(mesh2uv_tr[verts[0]], mesh2uv_tr[verts[1]], mesh2uv_tr[verts[2]]);
            cout << "created new face " << *newface << " with vertices " << verts[0] << verts[1] << verts[2] << endl;
        }
    }

    // create a texture

    const int x_size = 256;
    const int y_size = 256;

    PNG texture(256, 256);

    return 0;

    // iterate pixels and sample
    for (int y = 0; y < y_size; y++)
        for (int x = 0; x < x_size; x++)
        {
            for (SurfaceMesh::Face_index face_i : mesh.faces())
            {
                cout << "mesh face index: " << face_i << endl;
                Triangulation T;
                auto h = mesh.halfedge(face_i);
                for (vertex_descriptor v : mesh.vertices_around_face(h))
                {
                    cout << "found mesh vertex: " << mesh.point(v) << endl;
                    T.insert(uv_map[v]);
                    cout << "adding UV vertex: " << uv_map[v] << endl;
                }
                Triangulation::Face_handle huh = T.locate(Point_2(x, y));
                bool ok = false;
                for (auto f : T.all_face_handles())
                {
                    if (f == huh)
                        ok = true;
                }
            }
        }

    

    // print points
    // for (size_t i = 0; i < points.size(); ++i)
    // {
    //     const Point &p = get<0>(points[i]);
    //     const Vector &n = get<1>(points[i]);
    //     const Color &c = get<2>(points[i]);
    //     cout << "Point (" << p << ") with normal (" << n
    //          << ") and color (" << int(c[0]) << " " << int(c[1]) << " " << int(c[2])
    //          << ")" << endl;
    // }

    return 0;
}