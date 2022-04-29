#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <map>
#include <utility>

// arg parsing
#include "CLI11.hpp"

// output
#include "png.h"
extern "C"
{
#include "obj/obj.h"
}

// CGAL
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/boost/iterator/counting_iterator.hpp>
#include <CGAL/boost/graph/iterator.h>

// point cloud processing
#include <CGAL/IO/read_ply_points.h>
#include <CGAL/compute_average_spacing.h>
#include <CGAL/remove_outliers.h>
#include <CGAL/grid_simplify_point_set.h>
#include <CGAL/jet_smooth_point_set.h>
#include <CGAL/pca_estimate_normals.h>
#include <CGAL/mst_orient_normals.h>

// poisson recon
#include <CGAL/Polyhedron_3.h>
// #include <CGAL/Surface_mesh_default_triangulation_3.h>
// #include <CGAL/make_surface_mesh.h>
// #include <CGAL/Implicit_surface_3.h>
//#include <CGAL/IO/facets_in_complex_2_to_triangle_mesh.h>
//#include <CGAL/Poisson_reconstruction_function.h>
#include <CGAL/poisson_surface_reconstruction.h>

// UV unwrapping
#include <CGAL/Surface_mesh/Surface_mesh.h>
#include <CGAL/property_map.h>
#include <CGAL/Surface_mesh_parameterization/parameterize.h>
#include <CGAL/Surface_mesh_parameterization/Mean_value_coordinates_parameterizer_3.h>
#include <CGAL/Surface_mesh_parameterization/Square_border_parameterizer_3.h>
#include <CGAL/Surface_mesh_parameterization/IO/File_off.h>

// point cloud sampling
#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/Search_traits_adapter.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

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

// Point with normal vector
typedef std::pair<Point, Vector> PN;

typedef CGAL::Surface_mesh<Point> SurfaceMesh;

typedef boost::graph_traits<SurfaceMesh>::vertex_descriptor vertex_descriptor;
typedef boost::graph_traits<SurfaceMesh>::halfedge_descriptor halfedge_descriptor;
typedef boost::graph_traits<SurfaceMesh>::face_descriptor face_descriptor;

class My_point_property_map
{
    const std::vector<Point> &points;

public:
    typedef Point value_type;
    typedef const value_type &reference;
    typedef std::size_t key_type;
    typedef boost::lvalue_property_map_tag category;
    My_point_property_map(const std::vector<Point> &pts) : points(pts) {}
    reference operator[](key_type k) const { return points[k]; }
    friend reference get(const My_point_property_map &ppmap, key_type i)
    {
        return ppmap[i];
    }
};

typedef CGAL::Random_points_in_cube_3<Point> Random_points_iterator;
typedef CGAL::Search_traits_3<K> Traits_base;
typedef CGAL::Search_traits_adapter<std::size_t, My_point_property_map, Traits_base> Traits;
typedef CGAL::Orthogonal_k_neighbor_search<Traits> Neighbor_search;
typedef Neighbor_search::Tree Tree;
typedef Tree::Splitter Splitter;
typedef Neighbor_search::Distance Distance;

namespace SMP = CGAL::Surface_mesh_parameterization;
typedef SMP::Square_border_arc_length_parameterizer_3<SurfaceMesh> Border_parameterizer;
typedef SMP::Mean_value_coordinates_parameterizer_3<SurfaceMesh, Border_parameterizer> Parameterizer;

// simple surface recon

// advanced surface recon
typedef CGAL::First_of_pair_property_map<PN> PN_point_map;
typedef CGAL::Second_of_pair_property_map<PN> PN_normal_map;
typedef CGAL::Polyhedron_3<K> Polyhedron;
// typedef CGAL::Poisson_reconstruction_function<K> Poisson_recon_function;
// typedef CGAL::Surface_mesh_default_triangulation_3 SM_triangulation;
// typedef CGAL::Surface_mesh_complex_2_in_triangulation_3<SM_triangulation> C2tri3;
// typedef CGAL::Implicit_surface_3<K, Poisson_recon_function> Surface_3;

using namespace std;

int main(int argc, char **argv)
{
#pragma region argument parsing

    CLI::App app{"App description"};
    string meshFilename = "default";
    string cloudFilename = "default";
    app.add_option("-m,--in-mesh", meshFilename, "The input mesh in .ply format");
    app.add_option("-c,--in-cloud", cloudFilename, "The input cloud in .ply format");
    CLI11_PARSE(app, argc, argv);

#pragma endregion

#pragma region read point cloud

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
    else
    {
        cout << "Success: cloud loaded. " << points.size() << " points" << endl;
    }

    list<PN> simple_points;
    for (auto &&p : points)
    {
        simple_points.push_back(PN(get<0>(p), Vector(CGAL::NULL_VECTOR)));
    }

#pragma endregion

#pragma region outlier removal & simplification

    // calculate average spacing
    const unsigned int average_spacing_neighbors = 6;
    double average_spacing = CGAL::compute_average_spacing<CGAL::Sequential_tag>(
        simple_points,
        average_spacing_neighbors,
        CGAL::parameters::point_map(CGAL::First_of_pair_property_map<PN>()));
    cout << "Average spacing: " << average_spacing << endl;

    // Point with distance above 2*average_spacing are considered outliers
    // remove outliers
    std::list<PN>::iterator first_to_remove = CGAL::remove_outliers(
        simple_points,
        average_spacing_neighbors,
        CGAL::parameters::point_map(CGAL::First_of_pair_property_map<PN>()));
    cout << (100. * std::distance(first_to_remove, simple_points.end()) / (double)(simple_points.size()))
         << "% of the points are considered outliers when using a distance threshold of "
         << 2. * average_spacing << endl;
    simple_points.erase(first_to_remove, simple_points.end());

    // grid simplification
    const double cell_size = 0.01;
    first_to_remove = CGAL::grid_simplify_point_set(simple_points, cell_size, CGAL::parameters::point_map(CGAL::First_of_pair_property_map<PN>()));
    cout << (100. * std::distance(first_to_remove, simple_points.end()) / (double)(simple_points.size()))
         << "% of the points are culled by grid simplfication" << endl;
    simple_points.erase(first_to_remove, simple_points.end());

    cout << "point array size after SOR + grid simplification: " << simple_points.size() << endl;


#pragma endregion

#pragma region smoothing

    const unsigned int smoothing_neighbors = 12;
    CGAL::jet_smooth_point_set<CGAL::Sequential_tag>(
        simple_points,
        smoothing_neighbors,
        CGAL::parameters::point_map(CGAL::First_of_pair_property_map<PN>()));

#pragma endregion

#pragma region normals estimation

    unsigned int normal_est_neighbors = 12;
    CGAL::pca_estimate_normals<CGAL::Sequential_tag>(
        simple_points,
        normal_est_neighbors,
        CGAL::parameters::point_map(CGAL::First_of_pair_property_map<PN>()).normal_map(CGAL::Second_of_pair_property_map<PN>()));

    CGAL::mst_orient_normals(
        simple_points,
        normal_est_neighbors,
        CGAL::parameters::point_map(CGAL::First_of_pair_property_map<PN>()).normal_map(CGAL::Second_of_pair_property_map<PN>()));

#pragma endregion

#pragma region spacial search tree init

    vector<Point> rawpoints;
    rawpoints.reserve(points.size());
    for (size_t i = 0; i < points.size(); ++i)
    // for (size_t i = 0; i < 1000; ++i) // debug only - dont copy full cloud for less wait time
    {
        rawpoints.push_back(get<0>(points[i]));
    }
    My_point_property_map ppmap(rawpoints);

    Tree tree(
        boost::counting_iterator<size_t>(0),
        boost::counting_iterator<size_t>(rawpoints.size()),
        Splitter(),
        Traits(ppmap));
    Distance tr_dist(ppmap);

    cout << "Success: spatial search tree created. Initializing...";
    Neighbor_search init_search(tree, Point(0, 0, 0), 1, 0, true, tr_dist);
    for (Neighbor_search::iterator it = init_search.begin(); it != init_search.end(); it++)
    {
        cout << " d(q, nearest neighbor)=  "
             << tr_dist.inverse_of_transformed_distance(it->second) << " "
             << rawpoints[it->first] << " " << it->first << endl;
    }

    cout << " done." << endl;

#pragma endregion

#pragma region read mesh

    SurfaceMesh mesh;
    // ifstream in_m(meshFilename);
    // if (!in_m ||
    //     !CGAL::read_ply(
    //         in_m,
    //         mesh))
    // {
    //     cerr << "Error: unable to read from file " << meshFilename << endl;
    //     return EXIT_FAILURE;
    // }
    // else
    // {
    //     cout << "Success: loaded mesh " << meshFilename << endl;
    //     cout << mesh.number_of_vertices() << " verts, " << mesh.number_of_faces() << " tris" << endl;
    //     /*
    //     cout << "Vertex list:" << endl;
    //     for (vertex_descriptor v : mesh.vertices())
    //     {
    //         cout << v << ": " << mesh.point(v) << endl;
    //     }

    //     cout << "Face list:" << endl;
    //     for (face_descriptor face_i : mesh.faces())
    //     {
    //         cout << face_i << ": ";
    //         for (vertex_descriptor v : mesh.vertices_around_face(mesh.halfedge(face_i)))
    //             cout << v << " ";
    //         cout << endl;
    //     }
    //     */
    // }

#pragma endregion

#pragma region surface recon

    Polyhedron output_mesh;
    double average_spacing_2 = CGAL::compute_average_spacing<CGAL::Sequential_tag>(
        simple_points,
        24,
        CGAL::parameters::point_map(CGAL::First_of_pair_property_map<PN>()));

    CGAL::poisson_surface_reconstruction_delaunay
      (simple_points.begin(), simple_points.end(),
       CGAL::First_of_pair_property_map<PN>(),
       CGAL::Second_of_pair_property_map<PN>(),
       output_mesh, average_spacing, 5.0, 30.0, 5);
    
    if(
        true
    )
    {
        CGAL::copy_face_graph(output_mesh, mesh);
        cout << "surface reconstruction successful." << endl;
        cout << "number of vertices: " << mesh.num_vertices() << endl;
    }
    else
    {
        cerr << "Surface reconstruction failed. Exiting." << endl;
        return EXIT_FAILURE;
    }

#pragma endregion

#pragma region UV unwrap

    halfedge_descriptor bhd = CGAL::Polygon_mesh_processing::longest_border(mesh).first;

    typedef SurfaceMesh::Property_map<vertex_descriptor, Point_2> UV_pmap;
    UV_pmap uv_map = mesh.add_property_map<vertex_descriptor, Point_2>("h:uv").first;
    SMP::Error_code UV_err = SMP::parameterize(mesh,
                                               Parameterizer(),
                                               bhd,
                                               uv_map);

    if (UV_err != SMP::OK)
    {
        cerr << "Error: " << SMP::get_error_message(UV_err) << endl;
        return EXIT_FAILURE;
    }
    else
    {
        cout << "Success: created UV map for mesh" << endl;
        // cout << "UV Coords: " << endl;
        // for (vertex_descriptor v : mesh.vertices())
        // {
        //     cout << v << ": " << uv_map[v] << endl;
        // }
    }

#pragma endregion

#pragma region sample texture

    // create a texture
    const int x_size = 2048;
    const int y_size = 2048;
    PNG texture(x_size, y_size);
    cout << "Success: created " << x_size << "x" << y_size << " PNG texture" << endl
         << "Sampling cloud...";

    // iterate pixels and sample
    for (int y = 0; y < y_size; y++)
        for (int x = 0; x < x_size; x++)
        {
            Point_2 uv_point((float)x / float(x_size), (float)y / (float)y_size);
            bool found = false;
            Point_2 v0, v1, v2;
            vector<vertex_descriptor> vertices;
            for (face_descriptor face_d : mesh.faces())
            {
                vertices.clear();
                for (auto v : mesh.vertices_around_face(mesh.halfedge(face_d)))
                    vertices.push_back(v);

                v0 = uv_map[vertices[0]];
                v1 = uv_map[vertices[1]];
                v2 = uv_map[vertices[2]];

                K::Triangle_2 uv_face(v0, v1, v2);
                if (CGAL::do_intersect(uv_point, uv_face))
                {
                    found = true;
                    break;
                }
            }
            if (found)
            {
                double denominator = ((v1.y() - v2.y()) * (v0.x() - v2.x()) + (v2.x() - v1.x()) * (v0.y() - v2.y()));
                double a = ((v1.y() - v2.y()) * (uv_point.x() - v2.x()) + (v2.x() - v1.x()) * (uv_point.y() - v2.y())) / denominator;
                double b = ((v2.y() - v0.y()) * (uv_point.x() - v2.x()) + (v0.x() - v2.x()) * (uv_point.y() - v2.y())) / denominator;
                double c = 1.f - a - b;

                Point sample_location = CGAL::barycenter(mesh.point(vertices[0]),
                                                         a,
                                                         mesh.point(vertices[1]),
                                                         b,
                                                         mesh.point(vertices[2]));

                Neighbor_search search(tree, sample_location, 1, 0, true, tr_dist);
                for (Neighbor_search::iterator it = search.begin(); it != search.end(); ++it)
                {
                    Color c = get<2>(points[it->first]);
                    //                   vvvvvvvvvvvv Fill the image from bottom to top, since that is the way UV coords are oriented
                    texture.set_pixel(x, (y_size - 1) - y, c[0], c[1], c[2]);
                }
            }
        }

#pragma endregion

#pragma region save png texture

    const char *tex_filename = "out/tex.png";
    texture.save(tex_filename);
    cout << "Success: texture sampled and saved as " << tex_filename << endl;

#pragma endregion

#pragma region save obj model

    obj *o = obj_create(nullptr);

    for (vertex_descriptor v : mesh.vertices())
    {
        Point p = mesh.point(v);
        int vi = obj_add_vert(o);
        float position[3];
        position[0] = p.x();
        position[1] = p.y();
        position[2] = p.z();
        obj_set_vert_v(o, vi, position);
        float uv[2];
        uv[0] = uv_map[v].x();
        uv[1] = uv_map[v].y();
        obj_set_vert_t(o, vi, uv);
    }

    int o_surface_index = obj_add_surf(o);

    for (face_descriptor f : mesh.faces())
    {
        int i = 0;
        int fis[3];
        for (vertex_descriptor v : mesh.vertices_around_face(mesh.halfedge(f)))
        {
            fis[i++] = v.idx();
        }
        int fi = obj_add_poly(o, o_surface_index);
        obj_set_poly(o, o_surface_index, fi, fis);
    }

    int mi = obj_add_mtrl(o);
    obj_set_mtrl_name(o, mi, "texture");

    const char *obj_filename = "out/objtest.obj";
    const char *mtl_filename = "out/trash.mtl";

    obj_write(o, obj_filename, mtl_filename, 4);

    cout << "Success: object saved as " << obj_filename << endl;

#pragma endregion

}