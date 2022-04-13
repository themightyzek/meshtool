#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <map>

#include "CLI11.hpp"

#include "png.h"

extern "C" {
    #include "obj/obj.h"
}

#include "assimp/Importer.hpp"
#include "assimp/Exporter.hpp"
#include "assimp/scene.h"
#include "assimp/postprocess.h"
#include "assimp/vector3.h"
#include "assimp/material.h"
#include "assimp/mesh.h"
#include "assimp/types.h"
#include "assimp/DefaultLogger.hpp"

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh/Surface_mesh.h>
#include <CGAL/property_map.h>
#include <CGAL/IO/read_ply_points.h>
#include <CGAL/Surface_mesh_parameterization/parameterize.h>
#include <CGAL/Surface_mesh_parameterization/Mean_value_coordinates_parameterizer_3.h>
#include <CGAL/Surface_mesh_parameterization/Square_border_parameterizer_3.h>
#include <CGAL/Surface_mesh_parameterization/IO/File_off.h>

#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/Search_traits_adapter.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/boost/iterator/counting_iterator.hpp>
#include <CGAL/boost/graph/iterator.h>

#include <utility>

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

using namespace std;

int main(int argc, char **argv)
{
    // Assimp::Importer importer;
    // const aiScene* testScene = importer.ReadFile("test_in/vodka.obj", aiProcess_ValidateDataStructure);
    // cout << importer.GetErrorString() << endl;

    // aiTexture* testTex;
    // if(testScene->mMaterials[0]->GetTexture(aiTextureType_, 0, "vodkadöner.png") != AI_SUCCESS)
    // {

    // }


    // --------------------------------------------
    CLI::App app{"App description"};
    string meshFilename = "default";
    string cloudFilename = "default";
    app.add_option("-m,--in-mesh", meshFilename, "The input mesh in .ply format");
    app.add_option("-c,--in-cloud", cloudFilename, "The input cloud in .ply format");
    CLI11_PARSE(app, argc, argv);

    // --------------CGAL----------------5

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
    else
    {
        cout << "Success: cloud loaded. " << points.size() << " points" << endl;
    }

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
    // Tree tree(rawpoints.begin(), rawpoints.end());

    cout << "Success: spatial search tree created. Initializing...";
    Neighbor_search init_search(tree, Point(0, 0, 0), 1, 0, true, tr_dist);
    for (Neighbor_search::iterator it = init_search.begin(); it != init_search.end(); it++)
    {
        cout << " d(q, nearest neighbor)=  "
             << tr_dist.inverse_of_transformed_distance(it->second) << " "
             << rawpoints[it->first] << " " << it->first << endl;
    }

    cout << " done." << endl;

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
        for (vertex_descriptor v : mesh.vertices())
        {
            cout << v << ": " << mesh.point(v) << endl;
        }

        cout << "Face list:" << endl;
        for (face_descriptor face_i : mesh.faces())
        {
            cout << face_i << ": ";
            for (vertex_descriptor v : mesh.vertices_around_face(mesh.halfedge(face_i)))
                cout << v << " ";
            cout << endl;
        }
    }

    // UV unwrap

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
        cout << "Success: created UV map for mesh" << endl
             << "UV Coords: " << endl;
        for (vertex_descriptor v : mesh.vertices())
        {
            cout << v << ": " << uv_map[v] << endl;
        }
    }

    // create a texture

    const int x_size = 1024;
    const int y_size = 1024;

    PNG texture(x_size, y_size);

    // iterate pixels and sample
    for (int y = 0; y < y_size; y++)
        for (int x = 0; x < x_size; x++)
        {
            Point_2 uv_point((float)x / float(x_size), (float)y / (float)y_size);
            // cout << "sampling x" << x << " y" << y << " = " << uv_point << endl;
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
                    // cout << "found face containing P: " << face_d << endl;
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

                // cout << "calculated 3d sample location for P: " << sample_location << endl;
                Neighbor_search search(tree, sample_location, 1, 0, true, tr_dist);
                for (Neighbor_search::iterator it = search.begin(); it != search.end(); ++it)
                {
                    // cout << get<0>(points[it->first]) << " found with index " << it->first
                    //      << " . squared distance: " << it->second << endl;
                    Color c = get<2>(points[it->first]);
                    // cout << "color at target: R" << int(c[0]) << " G" << int(c[1]) << " B" << int(c[2]) << endl;
                    //                   vvvvvvvvvvvv Fill the image from bottom to top, since that is the way UV coords are oriented 
                    texture.set_pixel(x, (y_size-1)-y, c[0], c[1], c[2]);
                }
            }
        }

    texture.save("out/tex.png");

    
    
    obj* o = obj_create(nullptr);

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

    obj_write(o, "out/objtest.obj", "out/objtest.mtl", 4);



    /* ASSIMP STUFF
    
    Assimp::DefaultLogger::create("", Assimp::Logger::VERBOSE);

    vector<aiVector3D> assVertices;
    vector<aiVector3D> assUVs;
    for (vertex_descriptor v : mesh.vertices())
    {
        Point p = mesh.point(v);
        assVertices.push_back(
            aiVector3D(p.x(),
                       p.y(),
                       p.z()));
        assUVs.push_back(
            aiVector3D(
                uv_map[v].x(),
                uv_map[v].y(),
                0.0));
    }
    vector<aiFace> assFaces;
    for (face_descriptor f : mesh.faces())
    {
        aiFace face;
        unsigned i = 0;
        face.mIndices = new unsigned[3];
        for (vertex_descriptor v : mesh.vertices_around_face(mesh.halfedge(f)))
        {
            face.mIndices[i++] = (unsigned)v.idx();
        }
        face.mNumIndices = i;
        assFaces.push_back(face);
    }

    aiMesh *assMesh = new aiMesh();
    assMesh->mNumVertices = mesh.number_of_vertices();
    assMesh->mVertices = &assVertices.front();
    assMesh->mNumFaces = mesh.number_of_faces();
    assMesh->mFaces = &assFaces.front();
    assMesh->mPrimitiveTypes = aiPrimitiveType_TRIANGLE; // workaround, issue #3778
    assMesh->mTextureCoords[0] = &assUVs.front();
    assMesh->mNumUVComponents[0] = 2;

    aiTexture* assTexture = new aiTexture();
    assTexture->mFilename = "tex.png";

    aiMaterial* assMat = new aiMaterial();
    assMat->AddProperty("texture", 8, AI_MATKEY_NAME);
    aiString test;
    if(AI_SUCCESS == assMat->Get(AI_MATKEY_NAME, test))
    {
        cout << "mat name is " << test.C_Str() << endl;
    }

    aiScene scene;
    scene.mNumMeshes = 1;
    scene.mMeshes = new aiMesh *[1] { assMesh };
    scene.mNumTextures = 1;
    scene.mTextures = new aiTexture* [1] {assTexture};
    scene.mNumMaterials = 1;
    scene.mMaterials = new aiMaterial *[1]
    { new aiMaterial() };
    scene.mRootNode = new aiNode();
    scene.mRootNode->mNumMeshes = 1;
    scene.mRootNode->mMeshes = new unsigned[1]{0};
    scene.mMetaData = new aiMetadata(); // workaround, issue #3781

    Assimp::Exporter exporter;

    exporter.Export(&scene, "obj", "out/assimptest.obj", 0U, (const Assimp::ExportProperties *)nullptr);
    

    ofstream out("out/uv.off");
    SMP::IO::output_uvmap_to_off(mesh, bhd, uv_map, out);

    cout << "done, sadly, a segmentation fault will occur on exit" << endl;


    // ~aiScene() causes SIGSEGV, i dont know why
    // probably something going out of scope that ~aiScene then attempts to delete

    */

    return 0;
}