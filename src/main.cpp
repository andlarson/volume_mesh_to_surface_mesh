#include <iostream>
#include <fstream>
#include <string>
#include <unordered_map>
#include <array>

#include "CLI/CLI.hpp"

#include "gmsh.h"

using namespace std;

typedef array<double, 3>  Point3D_t;
typedef array<double, 3>  Vector3D_t;
typedef array<Point3D_t, 3> TriangularFace3D_t;

// A surface mesh is a collection of faces with associated normals. Each face 
//     has an associated normal.
typedef vector<pair<TriangularFace3D_t, Vector3D_t>> SurfaceTriangleMesh_t;


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
                  const SurfaceTriangleMesh_t& triangle_surface_mesh)
{
    ofstream f {stl_file_path};
    
    const string SOLID_NAME {"from_volumetric_gmsh"};
    f << "solid " << SOLID_NAME << endl;

    for (const auto& triangle : triangle_surface_mesh)
    {
        // TODO: Add normals!
        // Don't produce normals in the .stl file.
        f << "facet normal 0, 0, 0" << endl;
        f << "outer loop" << endl;
        for (const Point3D_t& point : triangle.first)
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
    Computes the cross product between two vectors.
*/
Vector3D_t cross_product(const Vector3D_t& v1, const Vector3D_t& v2)
{
    // The components of the two vectors to be crossed. 
    double v1_x {v1.at(0)};
    double v1_y {v1.at(1)};
    double v1_z {v1.at(2)};
    double v2_x {v2.at(0)};
    double v2_y {v2.at(1)};
    double v2_z {v2.at(2)};

    Vector3D_t ret;
    ret.at(0) = v1_y * v2_z - v1_z * v2_y;
    ret.at(1) = v1_z * v2_x - v1_x * v2_z;
    ret.at(1) = v1_x * v2_y - v1_y * v2_x;

    return ret;
}


/*
    Computes the vector pointing from one point to another. Returns the vector
        from the first point to the second point.
*/
Vector3D_t vector_between_points(const Point3D_t& p1, const Point3D_t& p2)
{
    return {p2.at(0) - p1.at(0), p2.at(1) - p1.at(1), p2.at(2) - p1.at(2)};
}


/*
    Computes the normal vector to a triangular face.
    The orientation of the returned normal with respect to the face depends 
        upon the order of the vertices stored in the face. This function respects
        the right hand rule.
    For example, if the triangular face comes from a tetrahedron and the vertices
        that compose the face are in counter clockwise order (when looking
        at the face on the tetrahedron), the resulting normal will point outside
        of the tetrahedron.
*/
Vector3D_t normal_to_tri_face(const TriangularFace3D_t& face)
{
    Vector3D_t v1 {vector_between_points(face.at(0), face.at(1))};
    Vector3D_t v2 {vector_between_points(face.at(0), face.at(2))};
    return cross_product(v1, v2);
}


/*
    Returns the surface mesh derived from the volumetric mesh contained in the
        first entity of the currently open Gmsh model. Assumes that the volumetric
        mesh of the currently open Gmsh model contains only triangular faces. 
*/
SurfaceTriangleMesh_t extract_surface_mesh()
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
    
    vector<size_t> face_nodes, face_tags;
    const uint64_t NODES_PER_TRIANGLE {3};
    gmsh::model::mesh::getAllFaces(NODES_PER_TRIANGLE, face_tags, face_nodes);

    // Create a mapping: For each face, the mapping will store the number of distinct
    //     elements that the face appears in. 
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
    const vector<vector<size_t>> groups_of_three_idxs {{0, 1, 2}, {0, 1, 3}, {0, 2, 3}, {1, 2, 3}};
    for (uint64_t type_idx {0}; type_idx < element_types.size(); ++type_idx)
    {
        const int type {element_types[type_idx]};
        for (const size_t element_tag : per_elem_type_element_tags[type_idx])
        {
            gmsh::model::mesh::getElement(element_tag, element_type, node_tags, dim, tag); 

            // For the collection of nodes that make up an element, get all the
            //     faces that are part of the element.
            vector<size_t> per_face_node_tags;
            for (const vector<size_t>& idxs : groups_of_three_idxs)
                for (const size_t idx : idxs)
                    per_face_node_tags.push_back(node_tags[idx]);
            gmsh::model::mesh::getFaces(NODES_PER_TRIANGLE, per_face_node_tags, face_tags, face_orientations);

            for (const size_t tag : face_tags)
                face_to_appearence_cnt.at(tag) += 1;             
        }
    }

    /* 
        For each face that belongs to only a single volumetric element, we would
            like to identify the nodes on the face and add the nodes to the surface
            mesh.
        This isn't so easy with Gmsh's API. The API does not offer a function that
            does: face tag -> node tags.
        The API does offer:
            getNodes(): entity tag -> nodes in entity
            getNodesbyElementType(): element type -> nodes in elements of the specified 
                type
            getNode(): node tag -> coordinates of the node
            getElements(): entity tag -> for each element type in entity, element 
                tags of the individual elements and node tags of all the nodes that 
                belong to all the elements of a particular type
            getElement(): element tag -> element type and node tags of nodes that 
                makeup element
            getFaces(): collection of node tags -> face tags corresponding to the
                nodes
            getAllFaces(): face type -> face tags and node tags of all the faces of
                the specified type
            getElementFaceNodes(): face type, element type -> node tags of all the
                nodes that belong to the specified face types on the specified
                element type
        The example "x7" provtages another hint:
            The getElementsByType() function returns element tags that map nicely
                to edge tags and face tags. Precisely, element tag 0 corresponds to 
                face tag 0 -> face tag 3, element tag 1 corresponds to face tag 4 ->
                face tag 7, etc. This makes it easy to map from a face tag -> element 
                tag and vice versa.
            I checked to see if there might be a relationship between the face tag
                and the node tag, based on the "x7" hint, but it doesn't seem like
                there is.
        Given Gmsh's API constraints, we there use the sequence of operations:
            (1) Identify the collection of surface elements. A surface element is
                    an element with one or more surface faces. Do this by using the
                    getElement() function to consider all the nodes of each element
                    individually and then by using getFaces() to identify the faces
                    on each element. If at least one of the faces of the element is
                    a surface face, the element is a surface element.
            (2) For each surface element, use getElement() to identify the nodes
                    that compose the surface element. For each group of three nodes,
                    use getFaces() to check if the corresponding face is a surface
                    face. If it is, save off the nodes in the surface mesh data
                    structure.
    */
    
    // Identify the surface elements using the surface faces.
    vector<size_t> surface_elements;
    for (uint64_t type_idx {0}; type_idx < element_types.size(); ++type_idx)
    {
        const int type {element_types[type_idx]};
        for (const size_t element_tag : per_elem_type_element_tags[type_idx])
        {
            // Get the nodes associated with the element.
            gmsh::model::mesh::getElement(element_tag, element_type, node_tags, dim, tag); 
            
            // For the collection of nodes that make up an element, get all the
            //     faces that are part of the element.
            vector<size_t> per_face_node_tags;
            for (const vector<size_t>& idxs : groups_of_three_idxs)
                for (const size_t idx : idxs)
                    per_face_node_tags.push_back(node_tags[idx]);
            gmsh::model::mesh::getFaces(NODES_PER_TRIANGLE, per_face_node_tags, face_tags, face_orientations);

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
    SurfaceTriangleMesh_t ret;
    
    /*
        The surface mesh representation requires consistent face normals. However, Gmsh's
            API doesn't directly offer a way to get consistent face normals.
            It does offer getNormal(), but that function requires a parametric
            representation of the coordinate at which to get the normal. It might
            be possible to compute the parametric representation of a point on
            each surface, and then pass that parametric representation to the
            getNormal() function, but the following approach turns out to be easier.
        It turns out the getElement() function offered by Gmsh's API always returns
            the node tags of the nodes that compose the element in a special order.
            In the case of a tetrahedron, the index permutations shown below
            guarantee that the nodes will be considered in the counter-clockwise 
            direction when looking into the tetrahedron.
                (0, 1, 3)
                (2, 3, 1)
                (0, 3, 2)
                (0, 2, 1)
        Knowing that the nodes on a face are always listed in counter-clockwise order,
            when looking into the element from the outside, makes it possible to
            compute the normal and ensure that the normal points outside. This
            makes it possible to compute consistent face normals.
    */
    const vector<vector<size_t>> gmsh_cw_tet_node_perms {{0, 1, 3}, {2, 3, 1}, {0, 3, 2}, {0, 2, 1}};

    // DEBUG: Consider one surface element only.
    // for (size_t i {1}; i < 2; i++)
    for (const size_t element_tag : surface_elements)
    {
        // DEBUG: Consider one surfaace element only.
        // const size_t element_tag {surface_elements[i]};

        gmsh::model::mesh::getElement(element_tag, element_type, node_tags, dim, tag);
        
        /*
        // DEBUG: Consider one surface element only.
        cout << "For this element, the order of the node tags is: ";
        for (const size_t node_tag : node_tags)
            cout << node_tag << " "; 
        cout << endl;
        vector<double> coord;
        vector<double> parametric_coord;
        cout << "For this element, the coordinates of the nodes are: " << endl;
        for (const size_t node_tag : node_tags)
        {
            gmsh::model::mesh::getNode(node_tag, coord, parametric_coord, dim, tag);
            cout << "Node tag " << node_tag << " : " << "(" << coord.at(0) << ", " << coord.at(1) << ", " << coord.at(2) << ")" << endl;
        }
        */

        // Consider each group of three nodes that make up the element.
        for (const vector<size_t>& idxs : gmsh_cw_tet_node_perms)
        {
            vector<size_t> potential_surface_node_tags;
            for (const size_t idx : idxs)
                potential_surface_node_tags.push_back(node_tags[idx]); 
            
            // Get the face the three nodes are a part of.
            gmsh::model::mesh::getFaces(NODES_PER_TRIANGLE, potential_surface_node_tags, face_tags, face_orientations);
            assert(face_tags.size() == 1);

            // DEBUG: Consider one surface element only.
            // if (face_to_appearence_cnt.at(face_tags[0]) != 0)

            // If these nodes correspond to a surface face, hurrah!
            if (face_to_appearence_cnt.at(face_tags[0]) == 1)
            {
                vector<double> coord;
                vector<double> parametric_coord;
                TriangularFace3D_t face;

                // DEBUG: Retrieving normal with parametric coords.
                // array<vector<double>, NODES_PER_TRIANGLE> per_node_parametric_coord;
                
                // Populate the face with three nodes.
                for (size_t idx {0}; idx < NODES_PER_TRIANGLE; ++idx)
                {
                    gmsh::model::mesh::getNode(potential_surface_node_tags[idx], coord, parametric_coord, dim, tag);
                    Point3D_t surface_node {coord.at(0), coord.at(1), coord.at(2)};
                    face[idx] = surface_node;

                    // DEBUG: Retrieving normal with parametric coords.
                    // per_node_parametric_coord.at(idx) = parametric_coord;
                }

                // DEBUG: Retrieving normal with parametric coords.
                // Compute the normal to this face using the parametric
                //     coordinates of the nodes that compose this face.
                //     Since the getNormal() only accepts parametric coordinates,
                //     we can compute the normal at the center of the surface
                //     face by computing the average parametric coordinate of
                //     the nodes that compose the face. This assumes that the
                //     parameterization is linear.
                /*
                vector<double> normals;
                const size_t PARAMETRIC_COORD_CNT {2};
                array<double, PARAMETRIC_COORD_CNT> average_parametric_coord;
                for (const vector<double>& pc : per_node_parametric_coord)
                {
                    assert(pc.size() == PARAMETRIC_COORD_CNT);                     
                    average_parametric_coord.at(0) += pc.at(0); 
                    average_parametric_coord.at(1) += pc.at(1); 
                }
                average_parametric_coord.at(0) /= per_node_parametric_coord.size();
                average_parametric_coord.at(1) /= per_node_parametric_coord.size();
                gmsh::model::getNormal(tag, {average_parametric_coord.begin(), average_parametric_coord.end()}, normals);
                assert(normals.size() == 1);
                Vector3D_t normal {normals.at(0), normals.at(1), normals.at(2)};

                ret.push_back({face, normal});
                */

                ret.push_back({face, normal_to_tri_face(face)});
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
