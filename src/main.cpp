#include <iostream>
#include <fstream>
#include <string>
#include <unordered_map>
#include <array>

#include "CLI/CLI.hpp"

#include "gmsh.h"

using namespace std;

typedef array<double, 3>  Point3D_t;
typedef array<Point3D_t, 3> TriangularFace3D_t;

// A surface mesh is a collection of faces. Each face is composed of exactly
//     3 vertices. This representation permits vertices to be repeated and
//     doesn't explicitly include normals.
typedef vector<TriangularFace3D_t> SurfaceMesh_t;


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
    Surface normals are not written to the .stl file. 
*/
void write_to_stl(const string& stl_file_path,
                  const SurfaceMesh_t& surface_mesh)
{
    ofstream f {stl_file_path};
    
    const string SOLID_NAME {"from_volumetric_gmsh"};
    f << "solid " << SOLID_NAME << endl;

    for (const TriangularFace3D_t& face : surface_mesh)
    {
        // Don't produce normals in the .stl file.
        f << "facet normal 0, 0, 0" << endl;
        f << "outer loop" << endl;

        for (const Point3D_t& point : face)
        {
            f << "vertex";
            for (const double coord : point)
            {
                f << " " << coord;
            }
            f << endl;
        }

        f << "end loop" << endl;
        f << "end facet" << endl;
    }

    f << "endsolid " << SOLID_NAME;
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
    // I can't find named constants in API.... 
    assert(entities.size() == 1);
    int dim {entities[0].first};
    int tag {entities[0].second};
    assert(dim == 3);
    string type;
    gmsh::model::getType(dim, tag, type);
    const string MESH_ENTITY_TYPE {"Discrete volume"};
    assert(type == MESH_ENTITY_TYPE);
    vector<int> element_types;
    gmsh::model::mesh::getElementTypes(element_types, dim, tag);
    assert(element_types.size() == 1);
    const size_t TETRAHEDRON_ELEMENT_TYPE {4};
    assert(element_types[0] == TETRAHEDRON_ELEMENT_TYPE);
}


/*
    Returns the surface mesh derived from the volumetric mesh contained in the
        first entity of the currently open Gmsh model. Assumes that the volumetric
        mesh of the currently open Gmsh model contains only triangular faces. 
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
    //     elements that the face appears in. 
    vector<size_t> face_nodes, face_tags;
    const uint64_t NODES_PER_TRIANGLE {3};
    gmsh::model::mesh::getAllFaces(NODES_PER_TRIANGLE, face_tags, face_nodes);
    unordered_map<size_t, uint64_t> face_to_appearence_cnt;
    for (const size_t& face_tag : face_tags)
        assert(face_to_appearence_cnt.insert({face_tag, 0}).second);

    // Get all the elements in the model.
    vector<vector<size_t>> per_elem_type_node_tags, per_elem_type_element_tags;
    gmsh::model::mesh::getElements(element_types, per_elem_type_element_tags, per_elem_type_node_tags, dim, tag);
    
    // For each element of every type, figure out which nodes are in the element
    //     and use the nodes to determine which faces are in the element. Every
    //     time a face is encountered, increment its count.
    vector<size_t> node_tags;
    int element_type;
    vector<int> face_orientations;
    for (uint64_t type_idx {0}; type_idx < element_types.size(); ++type_idx)
    {
        const int type {element_types[type_idx]};
        for (const size_t element_tag : per_elem_type_element_tags[type_idx])
        {
            gmsh::model::mesh::getElement(element_tag, element_type, node_tags, dim, tag); 

            // For the collection of nodes that make up an element, get all the
            //     faces that are part of the element.
            gmsh::model::mesh::getFaces(NODES_PER_TRIANGLE, node_tags, face_tags, face_orientations);

            for (const size_t tag : face_tags)
                face_to_appearence_cnt.at(tag) += 1;             
        }
    }
    
    // For each face that belongs to only a single volumetric element, we would
    //     like to identify the nodes on the face and add the nodes to the surface
    //     mesh.
    // This isn't so easy with Gmsh's API. The API does not offer a function that
    //     does: face tag -> node tags.
    // The API does offer:
    //     getNodes(): entity tag -> nodes in entity
    //     getNodesbyElementType(): element type -> nodes in elements of the specified 
    //         type
    //     getNode(): node tag -> coordinates of the node
    //     getElements(): entity tag -> for each element type in entity, element 
    //         tags of the individual elements and node tags of all the nodes that 
    //         belong to all the elements of a particular type
    //     getElement(): element tag -> element type and node tags of nodes that 
    //         makeup element
    //     getFaces(): collection of node tags -> face tags corresponding to the
    //         nodes
    //     getAllFaces(): face type -> face tags and node tags of all the faces of
    //         the specified type
    //     getElementFaceNodes(): face type, element type -> node tags of all the
    //         nodes that belong to the specified face types on the specified
    //         element type
    // The example "x7" provtages another hint:
    //     The getElementsByType() function returns element tags that map nicely
    //         to edge tags and face tags. Precisely, element tag 0 corresponds to 
    //         face tag 0 -> face tag 3, element tag 1 corresponds to face tag 4 ->
    //         face tag 7, etc. This makes it easy to map from a face tag -> element 
    //         tag and vice versa.
    //     I checked to see if there might be a relationship between the face tag
    //         and the node tag, based on the "x7" hint, but it doesn't seem like
    //         there is.
    // Given Gmsh's API constraints, we there use the sequence of operations:
    //     (1) Identify the collection of surface elements. A surface element is
    //             an element with one or more surface faces. Do this by using the
    //             getElement() function to consider all the nodes of each element
    //             individually and then by using getFaces() to identify the faces
    //             on each element. If at least one of the faces of the element is
    //             a surface face, the element is a surface element.
    //     (2) For each surface element, use getElement() to identify the nodes
    //             that compose the surface element. For each group of three nodes,
    //             use getFaces() to check if the corresponding face is a surface
    //             face. If it is, save off the nodes in the surface mesh data
    //             structure.
    
    // Identify the surface elements using the surface faces.
    vector<size_t> surface_elements;
    for (uint64_t type_idx {0}; type_idx < element_types.size(); ++type_idx)
    {
        const int type {element_types[type_idx]};
        for (const size_t element_tag : per_elem_type_element_tags[type_idx])
        {
            gmsh::model::mesh::getElement(element_tag, element_type, node_tags, dim, tag); 

            // For the collection of nodes that make up an element, get all the
            //     faces that are part of the element.
            gmsh::model::mesh::getFaces(NODES_PER_TRIANGLE, node_tags, face_tags, face_orientations);

            // If one of the faces is a surface face, then this is a surface
            //     element.
            for (const size_t face_tag : face_tags)
            {
                if (face_to_appearence_cnt.at(face_tag) == 1)
                {
                    surface_elements.push_back(element_tag);

                    // Even if a single element has multiple surface faces, it
                    //     should only be added to the surface elements collection
                    //     once.
                    break;
                }
            }
        }
    }
    
    // Identify and store the nodes on the surface faces using the surface
    //     elements.
    SurfaceMesh_t ret;
    const vector<vector<size_t>> groups_of_three_idxs {{0, 1, 2}, {0, 1, 3}, {0, 2, 3}, {1, 2, 3}};
    for (const size_t element_tag : surface_elements)
    {
        gmsh::model::mesh::getElement(element_tag, element_type, node_tags, dim, tag);

        // Consider each group of three nodes that make up the element.
        for (const auto& idxs : groups_of_three_idxs)
        {
            vector<size_t> potential_surface_node_tags;
            for (const auto& idx : idxs)
                potential_surface_node_tags.push_back(node_tags[idx]); 
            
            // Get the face the three nodes are a part of.
            gmsh::model::mesh::getFaces(NODES_PER_TRIANGLE, potential_surface_node_tags, face_tags, face_orientations);
            assert(face_tags.size() == 1);

            // If these nodes correspond to a surface face, hurrah!
            if (face_to_appearence_cnt.count(face_tags[0]) == 1)
            {
                vector<double> coord;
                vector<double> parametricCoord;
                TriangularFace3D_t face;
                
                // Populate the face with three nodes.
                for (size_t idx {0}; idx < NODES_PER_TRIANGLE; ++idx)
                {
                    gmsh::model::mesh::getNode(potential_surface_node_tags[idx], coord, parametricCoord, dim, tag);
                    Point3D_t surface_node {coord.at(0), coord.at(1), coord.at(2)};
                    face[idx] = surface_node;
                }

                // Add the face to the surface mesh.
                ret.push_back(face);
            }
        }
    }

    return ret;
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
    CLI::App app{"Converts a volumetric mesh in Gmsh's .msh format to a surface "
                 "mesh in .stl format."};
    argv = app.ensure_utf8(argv);
    
    string msh_file_path;
    app.add_option("-i,--input", msh_file_path, "Path to the .msh file "
                   "containing the volumetric mesh from which the surface mesh "
                   "should be extracted. Supports all versions of Gmsh's .msh file "
                   "format. No guarantees about the resulting surface "
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
