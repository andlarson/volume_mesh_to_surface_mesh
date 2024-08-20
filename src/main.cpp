#include <iostream>
#include <fstream>
#include <string>
#include <unordered_map>
#include <tuple>

#include "CLI/CLI.hpp"

#include "gmsh.h"

using namespace std;

typedef tuple<double, double, double> Point3D_t;

// A surface mesh is a collection of vertices and a collection of faces. A single
//     face is composed of 3 indices into the collection of faces.
typedef tuple<vector<Point3D_t>, vector<tuple<uint64_t, uint64_t, uint64_t>>> SurfaceMesh_t;


/*
    Checks that a string has .stl suffix. If it does have .stl suffix, returns
        empty string. Otherwise, it returns an error message.
*/
string has_stl_suffix(const string& str)
{
    if (!str.ends_with(".stl"))
        return "Desired output file doesn't end with .stl suffix!"; 
    else
        return "";
}


/*
    Writes the surface mesh data structure to the desired .stl file path. If the
        .stl file already exists, it is overwritten.
*/
void write_to_stl(const string& stl_file_path,
                  const SurfaceMesh_t& surface_mesh)
{

}


/*
    Checks that the currently open Gmsh model contains only a single entity
        of type Discrete Volume, composed of a mesh of tetrahedrons. 
*/
void check_gmsh_model()
{
    // Get all the elementary entities in the model.
    vector<pair<int, int>> entities;
    gmsh::model::getEntities(entities);
    
    // Make sure the model contains what we expect. 
    // I can't find documentation for some return values... 
    assert(entities.size() == 1);
    int dim {entities[0].first};
    int tag {entities[0].second};
    assert(dim == 3);
    string type;
    gmsh::model::getType(dim, tag, type);
    assert(type == "Discrete Volume");
    vector<int> element_types;
    gmsh::model::mesh::getElementTypes(element_types, dim, tag);
    assert(element_types.size() == 1);

    // DEBUG
    cout << "The single element type in the model has number: " << element_types[0] << endl;
}


/*
    Returns the surface mesh derived from the volumetric mesh contained in the
        first entity of the currently open Gmsh model. Assumes that the volumetric
        mesh of the currently open Gmsh model contains only triangular faces. 
    The surface mesh is represented in the form of a collection of non-repeating 
        nodes and a collection of triangular faces. The triangular faces are specified
        as indices into the collection of nodes.
*/
SurfaceMesh_t extract_surface_mesh()
{
    // Extract basic information about the model.
    vector<pair<int, int>> entities;
    gmsh::model::getEntities(entities);
    int dim {entities[0].first};
    int tag {entities[0].second};
    vector<int> element_types;
    gmsh::model::mesh::getElementTypes(element_types, dim, tag);

    // Create all the faces in this model. By default, the data structure for the
    //     model doesn't store faces.
    gmsh::model::mesh::createFaces(entities);
    
    // Create a mapping: For each face, the mapping stores the number of distinct
    //     elements that the face appears in. The faces on the surface of the
    //     volumetric mesh are those faces that belong to only a single volumetric
    //     mesh element.
    vector<size_t> face_nodes, face_tags;
    gmsh::model::mesh::getAllFaces(3, face_tags, face_nodes);
    unordered_map<size_t, uint64_t> face_to_appearence_cnt {face_tags.begin(), face_tags.end()};

    // Get all the elements in the model.
    vector<vector<size_t>> per_elem_node_tags, element_tags;
    gmsh::model::mesh::getElements(element_types, element_tags, per_elem_node_tags, dim, tag);
    
    // For each element of every type, figure out which nodes are in the element
    //     and use the nodes to determine which faces are in the element. Every
    //     time a face is encountered, increment its count.
    vector<size_t> node_tags;
    int element_type;
    vector<int> face_orientations;
    for (int type : element_types)
        for (size_t element_tag : element_tags[type])
        {
            gmsh::model::mesh::getElement(element_tag, element_type, node_tags, dim, tag); 

            // For the collection of nodes that make up an element, get all the
            //     faces that belong to the element.
            gmsh::model::mesh::getFaces(3, node_tags, face_tags, face_orientations);

            for (auto tag : face_tags)
                face_to_appearence_cnt.at(tag) += 1;             
        }
    
    // Iterate through all the faces again, and for each face that belongs to
    //     only a single volumetric element, add the face and the nodes belonging
    //     to the face to the surface mesh.

    // Maps a node to a unique identifier assigned to the node. Used to make
    //     sure that there are no duplicate nodes and also that each node has
    //     a unique identifier associated with it.
    uint64_t id {0};
    unordered_map<Point3D_t, uint64_t> nodes_to_ids;

    SurfaceMesh_t ret;
    for (

}


/*
    Given a Gmsh .msh file containing a volumetric mesh composed of tetrahedrons, 
        this function extracts the associated surface mesh and writes the result
        to the specified .stl file.
*/
void volume_mesh_to_surface_mesh(const string& msh_file_path,
                                 const string& stl_file_path)
{
    cout << "Extracting the surface mesh from the file: " << msh_file_path << endl;
    
    gmsh::open(msh_file_path);

    check_gmsh_model();
    
    auto surface_mesh {extract_surface_mesh()}; 

    write_to_stl(stl_file_path, surface_mesh);
}


int main(int argc, char *argv[])
{
    CLI::App app{"Converts a volumetric mesh in Gmsh's .msh format (version 4) "
                 "to a surface mesh in .stl format."};
    argv = app.ensure_utf8(argv);
    
    string msh_file_path;
    app.add_option("-i,--input", msh_file_path, "Path to the .msh file "
                   "containing the volumetric mesh from which the surface mesh "
                   "should be extracted. The .msh file should be in Gmsh's "
                   "version 4 format. No guarantees about the resulting surface "
                   "mesh are made: it may or may not contain a single connected "
                   "component, etc.")->required()->check(CLI::ExistingFile);

    string desired_stl_file_path;
    app.add_option("-o,--output", desired_stl_file_path, "Path to desired .stl "
                   "file. If this file already exists, it will be overwritten. "
                   "Must have .stl suffix.")->required()->check(has_stl_suffix);
    
    CLI11_PARSE(app, argc, argv);
    
    gmsh::initialize();
    volume_mesh_to_surface_mesh(msh_file_path, desired_stl_file_path);
    gmsh::finalize();    

    return EXIT_SUCCESS;
}
