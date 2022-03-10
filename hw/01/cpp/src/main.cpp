#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <unordered_map>

#include "Gmap.h"

void readObj(std::string &file_in, std::vector<Vertex> &vertices, std::vector<std::vector<int>> &face_indices) {
    /*
        readObj reads .obj file formats, and separates the vertices (x,y,z) from the faces('s indices) into
        two vectors. These two vectors can then be used for further processing.

        Input:  - (filepath to) .obj
                - two empty vectors for the vertices and face indices of types Vertex and int, respectively

        Output: - void function, stores output into passed vectors.
     */
    std::ifstream stream_in;
    stream_in.open(file_in); // open file

    if (stream_in.is_open()) {
        std::string line;

        while (getline(stream_in, line)) {
            std::istringstream iss(line);
            std::string word;
            iss >> word;

            // Extract all the vertices in the .obj, and store in vector vertices
            if (word == "v") {
                std::vector<float> coordinates;
                while (iss >> word) coordinates.push_back(std::stof(word));
                if (coordinates.size() == 3) vertices.emplace_back(coordinates[0], coordinates[1], coordinates[2]);
                else vertices.emplace_back();
            }

            // Extract all the faces('s indices) in the .obj, and store in vector face_indices
            if (word == "f") {
                std::vector<int> face;
                while (iss >> word) face.emplace_back(std::stoi(word));
                face_indices.push_back(face);
            }
        }
    }
    stream_in.close();
}

void extractCells(std::vector<Vertex> &vertices, std::vector<std::vector<int>> &face_indices,
                  std::unordered_map<std::string, Vertex> &vertexMap, std::unordered_map<std::string, Edge> &edgeMap,
                  std::vector<Face> &faceVec) {
    /*
        extractCells is a function that loops the raw data of the .obj file per 2-cell (face), and for each 2-cell
        extract the corresponding 1- and 0-cells. Each n-cell is stored into an appropriate n-container.

        Input:  - raw input of .obj file (vertices and face_indices)
        Output: - void function, stores output into appropriate containers (either vector or unordered_map)
     */

    for (int i = 0; i < face_indices.size(); i++) {

        // 2-cells
        Face face_cur; // container for single 2-cell

        // 0-cells
        for (int j = 0; j < face_indices[i].size(); j++) {

            // Initialise variables
            Vertex vertex_cur;
            std::string xyz;

            // Construct Vertex from current visiting point
            vertex_cur = Vertex(vertices[face_indices[i][j]-1]);

            // For storing the cells, we use unordered_map, this will prevent multiple addition of same cells
            // We convert the Point coordinates to a string, as this will ease the use of the unordered_map
            xyz = vertex_cur.xyz_tostring(vertex_cur.point.x,vertex_cur.point.y,vertex_cur.point.z);
            vertexMap.insert({xyz, vertex_cur});

            // 1-cells
            // container for single 1-cell
            Edge edge_cur;
            std::string edgeS;

            // 1. For finding all the connecting edges between the face vertices, it is needed to connect the last vertex
            // with first vertex. This is what the first conditional statement is for.
            if (face_indices[i][j] == face_indices[i].back()) {

                // 2. To get the same origin-vertex->end-vertex combination for each edge, we need to order them.
                // In this case, it is chosen to make the combination in ascending order.
                if ((face_indices[i][j]-1) > (face_indices[i][0]-1)) {
                    edge_cur = Edge(face_indices[i][0]-1, face_indices[i][j]-1);
                    edgeS = edge_cur.edge_tostring(edge_cur.origin_v, edge_cur.end_v);
                    edgeMap.insert({edgeS, edge_cur});
                }
                else {
                    edge_cur = Edge(face_indices[i][j]-1, face_indices[i][0]-1);
                    edgeS = edge_cur.edge_tostring(edge_cur.origin_v, edge_cur.end_v);
                    edgeMap.insert({edgeS, edge_cur});
                }
            }

            // 1. (see comments above)
            else {
                // 2. (see comments above)
                if ((face_indices[i][j]-1) > (face_indices[i][j+1]-1)) {
                    edge_cur = Edge(face_indices[i][j+1]-1, face_indices[i][j]-1);
                    edgeS = edge_cur.edge_tostring(edge_cur.origin_v, edge_cur.end_v);
                    edgeMap.insert({edgeS, edge_cur});
                }
                else {
                    edge_cur = Edge(face_indices[i][j]-1, face_indices[i][j+1]-1);
                    edgeS = edge_cur.edge_tostring(edge_cur.origin_v, edge_cur.end_v);
                    edgeMap.insert({edgeS, edge_cur});
                }
            }

            // For each index, push vertex into face vertices
            face_cur.face_vertices.push_back(face_indices[i][j]-1);

            // When it is the last iteration of the j loop, push the whole Face
            // into the vector of faces
            if (j == (face_indices[i].size() - 1)) {
                faceVec.push_back(face_cur);
            }
        }
    }
}

void generateCombinatorial(std::vector<Vertex> &vertices, std::vector<std::vector<int>> &face_indices,
                    std::unordered_map<std::string, Edge> &edgeMap, std::vector<Face> &faceVec, Volume &volume,
                    std::vector<Dart*> &darts) {
    /*
        initCombiStruct initialises the combinatorial structure that describes the primitives and the relationships
        between them. The loop mechanism follows the same logic as the one seen in extractCells, with the difference
        that darts are created, instead of cells. For each 1-cell, two darts are created. Each dart is assigned its
        corresponding vertex, edge, face and volume. Finally, involutions alpha_i are performed on all darts.
     */

    // For each 2-cell
    for (int i = 0; i < faceVec.size(); i++) {

        //For each 2-cell's vertices
        for (int j = 0; j < faceVec[i].face_vertices.size(); j++) {

            // Conditional to check if last 0-cell of 2-cell
            if (face_indices[i][j]-1 == faceVec[i].face_vertices.back()) {
                std::string origin_vS, end_vS, edgeS;

                // Fix ordering of edge key
                if ((face_indices[i][j]-1) > (faceVec[i].face_vertices.front())) {
                    origin_vS = std::to_string(faceVec[i].face_vertices.front());
                    end_vS = std::to_string(face_indices[i][j]-1);
                    edgeS = origin_vS + end_vS;
                    //construct first dart
                    Dart* dart_a = new Dart();
                    // assign vertex
                    dart_a->v = &vertices[face_indices[i][j]-1];
                    // assign edge
                    std::unordered_map<std::string, Edge>::iterator e_dart_a = edgeMap.find(edgeS);
                    dart_a->e = &e_dart_a->second;
                    // assign face
                    dart_a->f = &faceVec[i];
                    // assign volume
                    dart_a->vo = &volume;
                    darts.emplace_back(dart_a);
                    //construct second dart
                    Dart* dart_b = new Dart();
                    // assign vertex
                    dart_b->v = &vertices[faceVec[i].face_vertices.front()];
                    // assign edge
                    std::unordered_map<std::string, Edge>::iterator e_dart_b = edgeMap.find(edgeS);
                    dart_b->e = &e_dart_b->second;
                    // assign face
                    dart_b->f = &faceVec[i];
                    // assign volume
                    dart_b->vo = &volume;
                    darts.emplace_back(dart_b);
                }
                else {
                    origin_vS = std::to_string(face_indices[i][j]-1);
                    end_vS = std::to_string(faceVec[i].face_vertices.front());
                    edgeS = origin_vS + end_vS;
                    // construct first dart
                    Dart* dart_a = new Dart();
                    // assign vertex
                    dart_a->v = &vertices[face_indices[i][j]-1];
                    // assign edge
                    std::unordered_map<std::string, Edge>::iterator e_dart_a = edgeMap.find(edgeS);
                    dart_a->e = &e_dart_a->second;
                    // assign face
                    dart_a->f = &faceVec[i];
                    // assign volume
                    dart_a->vo = &volume;
                    darts.emplace_back(dart_a);
                    //construct second dart
                    Dart* dart_b = new Dart();
                    // assign vertex
                    dart_b->v = &vertices[faceVec[i].face_vertices.front()];
                    // assign edge
                    std::unordered_map<std::string, Edge>::iterator e_dart_b = edgeMap.find(edgeS);
                    dart_b->e = &e_dart_b->second;
                    // assign face
                    dart_b->f = &faceVec[i];
                    // assign volume
                    dart_b->vo = &volume;
                    darts.emplace_back(dart_b);
                }
            }
            else {
                std::string origin_vS, end_vS, edgeS;
                if ((face_indices[i][j]-1) > (face_indices[i][j+1]-1)) {
                    origin_vS = std::to_string(face_indices[i][j+1]-1);
                    end_vS = std::to_string(face_indices[i][j]-1);
                    edgeS = origin_vS + end_vS;
                    // construct first dart
                    Dart* dart_a = new Dart();
                    // assign vertex
                    dart_a->v = &vertices[face_indices[i][j]-1];
                    // assign edge
                    std::unordered_map<std::string, Edge>::iterator e_dart_a = edgeMap.find(edgeS);
                    dart_a->e = &e_dart_a->second;
                    // assign face
                    dart_a->f = &faceVec[i];
                    // assign volume
                    dart_a->vo = &volume;
                    darts.emplace_back(dart_a);
                    //construct second dart
                    Dart* dart_b = new Dart();
                    // assign vertex
                    dart_b->v = &vertices[face_indices[i][j+1]-1];
                    // assign edge
                    std::unordered_map<std::string, Edge>::iterator e_dart_b = edgeMap.find(edgeS);
                    dart_b->e = &e_dart_b->second;
                    // assign face
                    dart_b->f = &faceVec[i];
                    // assign volume
                    dart_b->vo = &volume;
                    darts.emplace_back(dart_b);
                }

                else {
                    origin_vS = std::to_string(face_indices[i][j]-1);
                    end_vS = std::to_string(face_indices[i][j+1]-1);
                    edgeS = origin_vS + end_vS;
                    //construct first dart
                    Dart* dart_a = new Dart();
                    // assign vertex
                    dart_a->v = &vertices[face_indices[i][j]-1];
                    // assign edge
                    std::unordered_map<std::string, Edge>::iterator e_dart_a = edgeMap.find(edgeS);
                    dart_a->e = &e_dart_a->second;
                    // assign face
                    dart_a->f = &faceVec[i];
                    // assign volume
                    dart_a->vo = &volume;
                    darts.emplace_back(dart_a);
                    //construct second dart
                    Dart* dart_b = new Dart();
                    // assign vertex
                    dart_b->v = &vertices[face_indices[i][j+1]-1];
                    // assign edge
                    std::unordered_map<std::string, Edge>::iterator e_dart_b = edgeMap.find(edgeS);
                    dart_b->e = &e_dart_b->second;
                    // assign face
                    dart_b->f = &faceVec[i];
                    // assign volume
                    dart_b->vo = &volume;
                    darts.emplace_back(dart_b);
                }
            }
        }
    }

    // Perform involutions for all darts
    for (auto i: darts) {
        for (auto j: darts) {
            i->invol_a0(i,j);
            i->invol_a1(i,j);
            i->invol_a2(i,j);
        }
    }
}

void generateEmbedding(std::unordered_map<std::string, Vertex> &vertexMap, std::unordered_map<std::string, Edge> &edgeMap,
                       std::vector<Face> &faceVec, Volume &volume, std::vector<Dart*> &darts) {
    /*
        The generateEmbedding function creates the embedding part of the Generalised map (Gmap). It consists of links
        to the cell structures of each dimension (ie Vertex, Edge, Face, Volume) and the darts they belong to.

        Input:  - the containers the cell structures and darts are stored in, passed by reference.
        Output: - (void) The input containers will be processed and contain the embedding information within
     */

    // Generate Vertex Embedding
    std::unordered_map<std::string, Vertex>:: iterator itrV;
    for (itrV = vertexMap.begin(); itrV != vertexMap.end(); itrV++) {
        for (auto j:darts) {
            if (itrV->second.xyz_tostring(itrV->second.point.x,
                                          itrV->second.point.y,
                                          itrV->second.point.z) == j->v->xyz_tostring(j->v->point.x,
                                                                                         j->v->point.y,
                                                                                         j->v->point.z)) {
                itrV->second.dart = j;
            }
        }
    }

    // Generate Edge Embedding
    std::unordered_map<std::string, Edge>:: iterator itrE;
    for (itrE = edgeMap.begin(); itrE != edgeMap.end(); itrE++) {
        for (auto j:darts) {
            if (itrE->second.edgeS == j->e->edgeS) {
                itrE->second.dart = j;
            }
        }
    }

    // Generate Face Embedding
    std::vector<Face>::iterator itrF;
    for (itrF = faceVec.begin(); itrF != faceVec.end(); itrF++) {
        for (auto j : darts) {
            if (itrF->face_vertices == j->f->face_vertices) {
                itrF->dart = j;
            }
        }
    }

    // Volume embedding, as we can assume there is only one volume we simply assign one arbitrary darts (ie the first)
    volume.dart = darts[0];
}


void calcBarycenters(std::unordered_map<std::string,Edge>&edgeMap, std::vector<Face>&faceVec, std::vector<Vertex>&vertices){

    /*
        This function calculates the barycenters of the 1- and 2-cells.

        Input:  - vector of faces, map of edges, vector of vertices
        Output: - the updated vectors and map with the calculated barycenters
    */

    std::unordered_map<std::string, Edge>::iterator itrE;
    for (itrE = edgeMap.begin(); itrE != edgeMap.end(); itrE++) {
        itrE->second.barycenter(vertices);
    }

    std::vector<Face>::iterator itrF;
    for (itrF = faceVec.begin(); itrF != faceVec.end(); itrF++) {
        itrF->barycenter(vertices);
    }
}

void storePointsOrdered (std::vector<Face> &faceVec, std::unordered_map<std::string, Point> &outObjVerticesUnique,
                         std::vector<Point> &outObjVertices, std::vector<std::string> &keys, std::vector<Vertex> &vertices,
                         std::unordered_map<std::string, Edge> &edgeMap) {
    /*
        In order to perform triangulation, we need to merge the barycenter points with the original points into a vector
        This function uses an unordered_map to prevent addition of duplicate points, and uses a vector of string keys to
        perform lookups

        Input:  - Containers holding all n-cells structs, raw input and additional containers as unordered_map for
        duplicate prevention and a vector of strings for lookup
        Output: - (void)
     */
    for (auto i : faceVec) {

        if (outObjVerticesUnique.find(i.barCF.xyz_tostring()) == outObjVerticesUnique.end()) {
            outObjVerticesUnique.insert({i.barCF.xyz_tostring(), i.barCF});
            outObjVertices.emplace_back(i.barCF);
            keys.emplace_back(i.barCF.xyz_tostring());
        }
        for (int j = 0; j < i.face_vertices.size(); j++) {
            if (i.face_vertices[j] == i.face_vertices.back()) {
                if (i.face_vertices[j] > i.face_vertices.front()) {
                    // vertices of triangle
                    Point a, b, c, d;
                    a = i.barCF;
                    b = vertices[i.face_vertices[j]].point;
                    // add vertex b if not added yet
                    if (outObjVerticesUnique.find(b.xyz_tostring()) == outObjVerticesUnique.end()) {
                        outObjVerticesUnique.insert({b.xyz_tostring(), b});
                        outObjVertices.emplace_back(b);
                        keys.emplace_back(b.xyz_tostring());
                    }
                    std::string origin_v, end_v, edgeKey;
                    origin_v = std::to_string(i.face_vertices.front());
                    end_v  = std::to_string(i.face_vertices[j]);
                    edgeKey = origin_v + end_v;   // make key of edge
                    std::unordered_map<std::string, Edge>::iterator edgeItr = edgeMap.find(edgeKey);
                    c = edgeItr->second.barCE;
                    // add vertex c if not added yet
                    if (outObjVerticesUnique.find(c.xyz_tostring()) == outObjVerticesUnique.end()) {
                        outObjVerticesUnique.insert({c.xyz_tostring(), c});
                        outObjVertices.emplace_back(c);
                        keys.emplace_back(c.xyz_tostring());
                    }
                    d = vertices[i.face_vertices.front()].point;
                    // add vertex d if not added yet
                    if (outObjVerticesUnique.find(d.xyz_tostring()) == outObjVerticesUnique.end()) {
                        outObjVerticesUnique.insert({d.xyz_tostring(), d});
                        outObjVertices.emplace_back(d);
                        keys.emplace_back(d.xyz_tostring());
                    }
                }
                else {
                    // vertices of triangle
                    Point a, b, c, d;
                    a = i.barCF;
                    b = vertices[i.face_vertices[j]].point;
                    // add vertex b if not added yet
                    if (outObjVerticesUnique.find(b.xyz_tostring()) == outObjVerticesUnique.end()) {
                        outObjVerticesUnique.insert({b.xyz_tostring(), b});
                        outObjVertices.emplace_back(b);
                        keys.emplace_back(b.xyz_tostring());
                    }
                    std::string origin_v, end_v, edgeKey;
                    origin_v = std::to_string(i.face_vertices[j]);
                    end_v = std::to_string(i.face_vertices.front());
                    edgeKey = origin_v + end_v; // make key of edge
                    std::unordered_map<std::string, Edge>::iterator edgeItr = edgeMap.find(edgeKey);
                    c = edgeItr->second.barCE;
                    if (outObjVerticesUnique.find(c.xyz_tostring()) == outObjVerticesUnique.end()) {
                        outObjVerticesUnique.insert({c.xyz_tostring(), c});
                        outObjVertices.emplace_back(c);
                        keys.emplace_back(c.xyz_tostring());
                    }
                    // add vertex d if not added yet
                    d = vertices[i.face_vertices.front()].point;
                    if (outObjVerticesUnique.find(d.xyz_tostring()) == outObjVerticesUnique.end()) {
                        outObjVerticesUnique.insert({d.xyz_tostring(), d});
                        outObjVertices.emplace_back(d);
                        keys.emplace_back(d.xyz_tostring());
                    }
                }
            }
            else {
                if (i.face_vertices[j] > i.face_vertices[j+1]) {
                    // vertices of two triangle
                    Point a, b, c, d;
                    a = i.barCF;
                    b = vertices[i.face_vertices[j]].point;
                    // add vertex b if not added yet
                    if (outObjVerticesUnique.find(b.xyz_tostring()) == outObjVerticesUnique.end()) {
                        outObjVerticesUnique.insert({b.xyz_tostring(), b});
                        outObjVertices.emplace_back(b);
                        keys.emplace_back(b.xyz_tostring());
                    }
                    std::string origin_v, end_v, edgeKey;
                    origin_v = std::to_string(i.face_vertices[j+1]);
                    end_v = std::to_string(i.face_vertices[j]);
                    edgeKey = origin_v + end_v;   // make key of edge
                    std::unordered_map<std::string, Edge>::iterator edgeItr = edgeMap.find(edgeKey);
                    c = edgeItr->second.barCE;
                    // add vertex c if not added yet
                    if (outObjVerticesUnique.find(c.xyz_tostring()) == outObjVerticesUnique.end()) {
                        outObjVerticesUnique.insert({c.xyz_tostring(), c});
                        outObjVertices.emplace_back(c);
                        keys.emplace_back(c.xyz_tostring());
                    }
                    d = vertices[i.face_vertices[j+1]].point;
                    // add vertex d if not added yet
                    if (outObjVerticesUnique.find(d.xyz_tostring()) == outObjVerticesUnique.end()) {
                        outObjVerticesUnique.insert({d.xyz_tostring(), d});
                        outObjVertices.emplace_back(d);
                        keys.emplace_back(d.xyz_tostring());
                    }
                }
                else {
                    // vertices of two triangle
                    Point a, b, c, d;
                    a = i.barCF;
                    b = vertices[i.face_vertices[j]].point;
                    // add vertex b if not added yet
                    if (outObjVerticesUnique.find(b.xyz_tostring()) == outObjVerticesUnique.end()) {
                        outObjVerticesUnique.insert({b.xyz_tostring(), b});
                        outObjVertices.emplace_back(b);
                        keys.emplace_back(b.xyz_tostring());
                    }
                    std::string origin_v, end_v, edgeKey;
                    origin_v = std::to_string(i.face_vertices[j]);
                    end_v = std::to_string(i.face_vertices[j+1]);
                    edgeKey = origin_v + end_v;  // make key of edge
                    std::unordered_map<std::string, Edge>::iterator edgeItr = edgeMap.find(edgeKey);
                    c = edgeItr->second.barCE;
                    // add vertex c if not added yet
                    if (outObjVerticesUnique.find(c.xyz_tostring()) == outObjVerticesUnique.end()) {
                        outObjVerticesUnique.insert({c.xyz_tostring(), c});
                        outObjVertices.emplace_back(c);
                        keys.emplace_back(c.xyz_tostring());
                    }

                    d = vertices[i.face_vertices[j+1]].point;
                    // add vertex d if not added yet
                    if (outObjVerticesUnique.find(d.xyz_tostring()) == outObjVerticesUnique.end()) {
                        outObjVerticesUnique.insert({d.xyz_tostring(), d});
                        outObjVertices.emplace_back(d);
                        keys.emplace_back(d.xyz_tostring());
                    }
                }
            }
        }
    }
}

void connectTriangles(std::vector<Face> &faceVec, std::vector<std::string> &keys,
                      std::vector<Vertex> &vertices, std::unordered_map<std::string, Edge> &edgeMap,
                      std::vector<std::vector<int>> & outObjTriangles) {
    /*
        This function connects all the triangles with their corresponding vertices, and stores them in CCW order
        in a 1-based index vector suitable for writing to .obj

        Input:  - containers with raw input, edge and face structs, output vector for indexed triangles, and vector containing keys
        Output: - void
     */

    for (auto i : faceVec) {

        auto itra = std::find(keys.begin(), keys.end(), i.barCF.xyz_tostring());

        for (int j = 0; j < i.face_vertices.size(); j++) {
            if (i.face_vertices[j] == i.face_vertices.back()) {
                if (i.face_vertices[j] > i.face_vertices.front()) {

                    // vertices of triangle
                    Point a, b, c, d;
                    a = i.barCF;
                    b = vertices[i.face_vertices[j]].point;

                    std::string origin_v, end_v, edgeKey;
                    origin_v = std::to_string(i.face_vertices.front());
                    end_v  = std::to_string(i.face_vertices[j]);
                    edgeKey = origin_v + end_v;
                    std::unordered_map<std::string, Edge>::iterator edgeItr = edgeMap.find(edgeKey);
                    c = edgeItr->second.barCE;

                    d = vertices[i.face_vertices.front()].point;

                    // Store triangle's face indices and push into outObjTriangles
                    // abc
                    auto itrb = std::find(keys.begin(), keys.end(), b.xyz_tostring());
                    auto itrc = std::find(keys.begin(), keys.end(), c.xyz_tostring());
                    auto itrd = std::find(keys.begin(), keys.end(), d.xyz_tostring());

                    int ia, ib, ic, id;
                    std::vector<int> tri1, tri2;
                    ia = itra - keys.begin() + 1;
                    ib = itrb - keys.begin() + 1;
                    ic = itrc - keys.begin() + 1;
                    id = itrd - keys.begin() + 1;

                    tri1.emplace_back(ia); tri1.emplace_back(ib); tri1.emplace_back(ic);
                    outObjTriangles.emplace_back(tri1);
                    //acb
                    tri2.emplace_back(ia); tri2.emplace_back(ic); tri2.emplace_back(id);
                    outObjTriangles.emplace_back(tri2);


                }
                else {

                    // vertices of triangle
                    Point a, b, c, d;
                    a = i.barCF;
                    b = vertices[i.face_vertices[j]].point;

                    std::string origin_v, end_v, edgeKey;
                    origin_v = std::to_string(i.face_vertices[j]);
                    end_v = std::to_string(i.face_vertices.front());
                    edgeKey = origin_v + end_v;
                    std::unordered_map<std::string, Edge>::iterator edgeItr = edgeMap.find(edgeKey);
                    c = edgeItr->second.barCE;

                    d = vertices[i.face_vertices.front()].point;

                    // Store triangle's face indices and push into outObjTriangles
                    // abc
                    auto itrb = std::find(keys.begin(), keys.end(), b.xyz_tostring());
                    auto itrc = std::find(keys.begin(), keys.end(), c.xyz_tostring());
                    auto itrd = std::find(keys.begin(), keys.end(), d.xyz_tostring());

                    int ia, ib, ic, id;
                    std::vector<int> tri1, tri2;
                    ia = itra - keys.begin() + 1;
                    ib = itrb - keys.begin() + 1;
                    ic = itrc - keys.begin() + 1;
                    id = itrd - keys.begin() + 1;

                    tri1.emplace_back(ia); tri1.emplace_back(ib); tri1.emplace_back(ic);
                    outObjTriangles.emplace_back(tri1);
                    //acb
                    tri2.emplace_back(ia); tri2.emplace_back(ic); tri2.emplace_back(id);
                    outObjTriangles.emplace_back(tri2);
                }


            }
            else {
                if (i.face_vertices[j] > i.face_vertices[j+1]) {

                    // vertices of two triangle
                    Point a, b, c, d;
                    a = i.barCF;
                    b = vertices[i.face_vertices[j]].point;

                    std::string origin_v, end_v, edgeKey;
                    origin_v = std::to_string(i.face_vertices[j+1]);
                    end_v = std::to_string(i.face_vertices[j]);
                    edgeKey = origin_v + end_v;
                    std::unordered_map<std::string, Edge>::iterator edgeItr = edgeMap.find(edgeKey);
                    c = edgeItr->second.barCE;

                    d = vertices[i.face_vertices[j+1]].point;

                    // Store triangle's face indices and push into outObjTriangles
                    // abc
                    auto itrb = std::find(keys.begin(), keys.end(), b.xyz_tostring());
                    auto itrc = std::find(keys.begin(), keys.end(), c.xyz_tostring());
                    auto itrd = std::find(keys.begin(), keys.end(), d.xyz_tostring());

                    int ia, ib, ic, id;
                    std::vector<int> tri1, tri2;
                    ia = itra - keys.begin() + 1;
                    ib = itrb - keys.begin() + 1;
                    ic = itrc - keys.begin() + 1;
                    id = itrd - keys.begin() + 1;

                    tri1.emplace_back(ia); tri1.emplace_back(ib); tri1.emplace_back(ic);
                    outObjTriangles.emplace_back(tri1);
                    //acb
                    tri2.emplace_back(ia); tri2.emplace_back(ic); tri2.emplace_back(id);
                    outObjTriangles.emplace_back(tri2);
                }
                else {

                    // vertices of two triangle
                    Point a, b, c, d;
                    a = i.barCF;
                    b = vertices[i.face_vertices[j]].point;

                    // make key of edge
                    std::string origin_v, end_v, edgeKey;
                    origin_v = std::to_string(i.face_vertices[j]);
                    end_v = std::to_string(i.face_vertices[j+1]);
                    edgeKey = origin_v + end_v;
                    std::unordered_map<std::string, Edge>::iterator edgeItr = edgeMap.find(edgeKey);
                    c = edgeItr->second.barCE;

                    d = vertices[i.face_vertices[j+1]].point;

                    // Store triangle's face indices and push into outObjTriangles
                    // abc
                    auto itrb = std::find(keys.begin(), keys.end(), b.xyz_tostring());
                    auto itrc = std::find(keys.begin(), keys.end(), c.xyz_tostring());
                    auto itrd = std::find(keys.begin(), keys.end(), d.xyz_tostring());

                    int ia, ib, ic, id;
                    std::vector<int> tri1, tri2;
                    ia = itra - keys.begin() + 1;
                    ib = itrb - keys.begin() + 1;
                    ic = itrc - keys.begin() + 1;
                    id = itrd - keys.begin() + 1;

                    tri1.emplace_back(ia); tri1.emplace_back(ib); tri1.emplace_back(ic);
                    outObjTriangles.emplace_back(tri1);
                    //acb
                    tri2.emplace_back(ia); tri2.emplace_back(ic); tri2.emplace_back(id);
                    outObjTriangles.emplace_back(tri2);
                }
            }
        }
    }
}

void outDartsCSV(std::string darts_file, std::vector<Dart*> &darts){
  /*
      Output of darts into csv file containing information of the darts, id is pointer based


      Input:    - darts_file => "darts.csv"
      Output:   - darts.csv
   */

    std::ofstream out_darts( darts_file );
    out_darts << "id;a0;a1;a2;a3;v;e;f" << std::endl;
    for (auto i : darts) {
        out_darts << i<< ";";
        out_darts << i->a0 << ";";
        out_darts << i->a1 << ";";
        out_darts << i->a2 << ";";
        out_darts << i->v << ";";
        out_darts << i->e << ";";
        out_darts << i->f;
        out_darts<<"\n";
    }
    out_darts.close();
}

void outVerticesCSV(std::string vertex_file, std::unordered_map<std::string, Vertex> &vertexMap){
    /*
        Output of vertices into the vertices csv, id is 0-index based
        Info: id, relating dart, coordinates (x, y, z)

        Input:  - vertex_file and file location of vertex csv
        Output: - vertex.csv
     */

    // one CSV file for the 0-cells (ending on vertices.csv) with at least the columns ID, dart, x, y, and z,
    std::ofstream out_vert (vertex_file);
    out_vert<<"id;dart;x;y;z\n";
    int vertind=0;
    for (auto i : vertexMap) {
        out_vert << vertind << ";";
        out_vert << i.second.dart << ";";
        out_vert <<i.second.point.x  << ";";
        out_vert<<i.second.point.y << ";";
        out_vert<<i.second.point.z;
        out_vert<<"\n";
        vertind++;
    }
    out_vert.close();
}

void outEdgesCSV(std::string edges_file, std::unordered_map<std::string, Edge> &edgeMap){
    /*
        Output of edges into the edges csv, id is 0-index based
        Info: id, relating dart

        Input:  - edge_file and file location of edges csv
        Output: - edges.csv
     */
    std::ofstream out_edges (edges_file);
    out_edges<<"id;dart"<<std::endl;
    int edgeIdx = 0;
    for (auto i : edgeMap) {
        out_edges << edgeIdx << ";";
        out_edges << i.second.dart;
        out_edges<<"\n";
        edgeIdx++;
    }
    out_edges.close();
}

void outFacesCSV(std::string faces_file, std::vector<Face> &faceVec){
    /*
        Output of faces into the faces csv, id is 0-index based
        Info: id, relating dart

        Input:  - faces_file and file location of faces csv
        Output: - faces.csv
     */
    std::ofstream out_face (faces_file);
    out_face<<"id;dart"<<"\n";
    int idFace = 0;
    for (auto i : faceVec) {
        out_face << idFace << ";";
        out_face << i.dart;
        out_face<<"\n";
        idFace++;
    }
    out_face.close();
}

void outVolumeCSV(std::string volume_file, Volume &volume){
    /*
        Output of volume into the volume csv, id is 0-index based
        Info: id, relating dart

        Input:  - volume_file and file location of volume csv
        Output: - volume.csv
     */
    std::ofstream out_vol (volume_file);
    out_vol<<"id;dart\n";
    out_vol << 0 << ";";
    out_vol << volume.dart;
    out_vol.close();
}

void objWriter(std::string file_in, std::string file_out_obj, std::vector<Point>&objPoints,std::vector<std::vector<int>>&outObjTriangles) {
    /*
        objWriter writes Vertex' Points and the triangle face index mapping to .obj
        Note: index mapping should be 1-index based and in CCW order

        Input:  - file name, vector of object Points, vector of triangle pairs
        Output: - triangulated.obj
     */
    std::ofstream triObj(file_out_obj);
    for (auto i: objPoints) {
        triObj << "v ";
        triObj << std::fixed << std::setprecision(6) << i.x << " ";
        triObj << std::fixed << std::setprecision(6) << i.y << " ";
        triObj << std::fixed << std::setprecision(6) << i.z << "\n";
    }

    for (int i = 0; i < outObjTriangles.size(); i++) {
        triObj << "f ";
        for (int j = 0; j < outObjTriangles[i].size(); j++) {
            triObj << outObjTriangles[i][j] << " ";

        }
        triObj << '\n';
    }

    triObj.close();
}

int main(int argc, const char * argv[]) {

    std::string file_in =        "/Users/danieldobson/Library/CloudStorage/OneDrive-Personal/GEOMATICS/GEO1004/assignments/geo1004.2122_hw01/hw/01/data/torus.obj";
    std::string cube_test =      "/Users/danieldobson/Library/CloudStorage/OneDrive-Personal/GEOMATICS/GEO1004/assignments/geo1004.2122_hw01/hw/01/data/cube.obj";
    std::string file_out_obj =   "/Users/danieldobson/Library/CloudStorage/OneDrive-Personal/GEOMATICS/GEO1004/assignments/geo1004.2122_hw01/hw/01/data/torus_triangulated.obj";
    std::string file_out_test =   "/Users/danieldobson/Library/CloudStorage/OneDrive-Personal/GEOMATICS/GEO1004/assignments/geo1004.2122_hw01/hw/01/data/cube_triangulated.obj";
    std::string file_out_csv_d = "/Users/danieldobson/Library/CloudStorage/OneDrive-Personal/GEOMATICS/GEO1004/assignments/geo1004.2122_hw01/hw/01/data/torus_darts.csv";
    std::string file_out_csv_0 = "/Users/danieldobson/Library/CloudStorage/OneDrive-Personal/GEOMATICS/GEO1004/assignments/geo1004.2122_hw01/hw/01/data/torus_vertices.csv";
    std::string file_out_csv_1 = "/Users/danieldobson/Library/CloudStorage/OneDrive-Personal/GEOMATICS/GEO1004/assignments/geo1004.2122_hw01/hw/01/data/torus_edges.csv";
    std::string file_out_csv_2 = "/Users/danieldobson/Library/CloudStorage/OneDrive-Personal/GEOMATICS/GEO1004/assignments/geo1004.2122_hw01/hw/01/data/torus_faces.csv";
    std::string file_out_csv_3 = "/Users/danieldobson/Library/CloudStorage/OneDrive-Personal/GEOMATICS/GEO1004/assignments/geo1004.2122_hw01/hw/01/data/torus_volume.csv";

    // ## Read OBJ file ##
    // The vertices and faces are read and stored into vectors.

    std::vector<Vertex> vertices;
    std::vector<std::vector<int>> face_indices;

    readObj(file_in, vertices, face_indices);

    // ## Construct generalised map using the structures from Gmap.h ##
    // The Gmap will be created in two parts
    // 1. A combinatorial structure
    // 2. An embedding structure

    // 1. Combinatorial structure
    // Initialisation of containers for darts and n-cells.
    std::vector<Dart*> darts;
    std::unordered_map<std::string, Vertex> vertexMap;
    std::unordered_map<std::string, Edge> edgeMap;
    std::vector<Face> faceVec;
    Volume volume;

    // Start building the combinatorial structure by extracting all cells from the object
    extractCells(vertices, face_indices, vertexMap, edgeMap, faceVec);

    // Continue generation of the combinatorial structure by initializing darts and
    // performing corresponding involutions alpha_i
    generateCombinatorial(vertices,face_indices,edgeMap,faceVec,volume,darts);

    // 2. Embedding structure
    // Consists of; Vertex, Edge and Face embedding
    generateEmbedding(vertexMap, edgeMap, faceVec, volume, darts);


    // ## Create triangles from the darts ##
    // The triangulation is achieved by calculating all the barycenters of the darts (edges in this case) and faces
    calcBarycenters(edgeMap, faceVec, vertices);

    // Next, we traverse all the faces and their corresponding vertices edge-wise
    // By doing so, we can connect the face barycenter with the two edge vertices
    // and the edge barycenter resulting in two triangles.

    std::vector<Point> outObjVertices; // Store vertices as Point for .obj
    std::unordered_map<std::string, Point> outObjVerticesUnique; // Store vertices as Point for .obj without duplicates
    std::vector<std::string> keys; // Store the vertex keys for lookup operations later
    std::vector<std::vector<int>> outObjTriangles; // Store triangle indices as 1-based index for .obj

    // Store all Points, including the barycenters, ordered
    storePointsOrdered(faceVec, outObjVerticesUnique, outObjVertices, keys, vertices, edgeMap);

    // Connect triangles, and store in 1-based index
    connectTriangles(faceVec, keys, vertices, edgeMap, outObjTriangles);

    // ## Output generalised map to CSV ##
    outDartsCSV(file_out_csv_d, darts);
    outVerticesCSV(file_out_csv_0, vertexMap);
    outEdgesCSV(file_out_csv_1, edgeMap);
    outFacesCSV(file_out_csv_2, faceVec);
    outVolumeCSV(file_out_csv_3, volume);

    // ## Write triangles to obj ##
    objWriter(file_in, file_out_obj, outObjVertices, outObjTriangles);

    return 0;
}