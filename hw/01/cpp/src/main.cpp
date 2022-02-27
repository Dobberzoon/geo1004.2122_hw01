// some standard libraries that are helpfull for reading/writing text files
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <unordered_map>

#include "Gmap.h"

void readObj(std::string &file_in, std::vector<Vertex> &vertices, std::vector<std::vector<int>> &face_indices) {
    std::ifstream stream_in;
    stream_in.open(file_in);

    if (stream_in.is_open()) {
        std::string line;
        while (getline(stream_in, line)) {
            std::istringstream iss(line);
            std::string word;
            iss >> word;
            if (word == "v") {
                std::vector<float> coordinates;
                while (iss >> word) coordinates.push_back(std::stof(word));
                if (coordinates.size() == 3) vertices.emplace_back(coordinates[0], coordinates[1], coordinates[2]);
                else vertices.emplace_back();
            }

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
                  std::vector<Face> &faceVec, Volume &volume) {

    // This loop traverses all faces (per indices), and in the double loop we traverse the vertices
    // that make up each face.

    for (int i = 0; i < face_indices.size(); i++) {

        // 2-cells
        Face face_cur;


        // the std::cout's are only for visualising the loop process
        //std::cout << "( ";


        // 0-cells
        for (int j = 0; j < face_indices[i].size(); j++) {
            //std::cout << face_indices[i][j] << " ";

            // Initialise variables
            Vertex vertex_cur;
            std::string xyz;

            // Construct Vertex from current visiting point
            //std::cout << "vertices[j-1].point.x: " << vertices[j-1].point.x << "\n";
            vertex_cur = Vertex(vertices[face_indices[i][j]-1]);

            // For storing the cells, we use unordered_map, this will prevent multiple addition of same cells
            xyz = vertex_cur.xyz_tostring(vertex_cur.point.x,vertex_cur.point.y,vertex_cur.point.z);
            vertexMap.insert({xyz, vertex_cur});

            // visualising first vertex repetition
            //if (face_indices[i][j] == face_indices[i].back()) {std::cout << face_indices[i][0];}

            // 1-cells
            Edge edge_cur;
            std::string edgeS;

            if (face_indices[i][j] == face_indices[i].back()) {
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

            else {
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

            // 2-cells
            face_cur.face_vertices.push_back(face_indices[i][j]-1);

            if (j == (face_indices[i].size() - 1)) {
                faceVec.push_back(face_cur);
            }
        }
        //std::cout << " ) \n";
    }
}

void initCombiStruct(std::vector<Vertex> &vertices, std::vector<std::vector<int>> &face_indices,
                    std::unordered_map<std::string, Edge> &edgeMap, std::vector<Face> &faceVec, Volume &volume,
                    std::vector<Dart*> &darts) {

    int countDart = 0;
    for (int i = 0; i < faceVec.size(); i++) {
        //std::cout << "f" << i << " \n"; //this code keeps track which face we are

        for (int j = 0; j < faceVec[i].face_vertices.size(); j++) {

            if (face_indices[i][j]-1 == faceVec[i].face_vertices.back()) {
                //std::cout << faceVec[i].face_vertices.front();
                std::string origin_vS, end_vS, edgeS;
                //origin_vS = std::to_string(face_indices[i][j]-1);
                //end_vS = std::to_string(faceVec[i].face_vertices.front());
                //edgeS = origin_vS + end_vS;

                if ((face_indices[i][j]-1) > (faceVec[i].face_vertices.front())) {
                    origin_vS = std::to_string(faceVec[i].face_vertices.front());
                    end_vS = std::to_string(face_indices[i][j]-1);
                    edgeS = origin_vS + end_vS;

                    // vertex track
                    countDart++;
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
                    std::cout << "dart " << countDart << ": v" << face_indices[i][j]-1 << ", e" << edgeS << ", f" << i << "\n";

                    countDart++;
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

                    // edge track

                    std::cout << "dart " << countDart << ": v" << faceVec[i].face_vertices.front() << ", e" << edgeS << ", f" << i << "\n";

                }
                else {
                    origin_vS = std::to_string(face_indices[i][j]-1);
                    end_vS = std::to_string(faceVec[i].face_vertices.front());
                    edgeS = origin_vS + end_vS;

                    // vertex track
                    countDart++;
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
                    std::cout << "dart " << countDart << ": v" << face_indices[i][j]-1 << ", e" << edgeS << ", f" << i << "\n";
                    countDart++;
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
                    std::cout << "dart " << countDart << ": v" << faceVec[i].face_vertices.front() << ", e" << edgeS << ", f" << i << "\n";
                    // edge track
                }
            }
            else {
                std::string origin_vS, end_vS, edgeS;

                if ((face_indices[i][j]-1) > (face_indices[i][j+1]-1)) {
                    origin_vS = std::to_string(face_indices[i][j+1]-1);
                    end_vS = std::to_string(face_indices[i][j]-1);
                    edgeS = origin_vS + end_vS;

                    // vertex track
                    countDart++;
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
                    std::cout << "dart " << countDart << ": v" << face_indices[i][j]-1 << ", e" << edgeS << ", f" << i << "\n";
                    countDart++;
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
                    std::cout << "dart " << countDart << ": v" << face_indices[i][j+1]-1 << ", e" << edgeS << ", f" << i << "\n";

                }
                else {
                    origin_vS = std::to_string(face_indices[i][j]-1);
                    end_vS = std::to_string(face_indices[i][j+1]-1);
                    edgeS = origin_vS + end_vS;

                    // vertex track
                    countDart++;
                    //construct first dart
                    Dart* dart_a = new Dart();
                    // assign vertex
                    dart_a->v = &vertices[face_indices[i][j]-1];
                    // assign edge
                    std::unordered_map<std::string, Edge>::iterator e_dart_a = edgeMap.find(edgeS);
                    //std::cout << "edge dart_a: " << e_dart_a->second.edgeS;
                    dart_a->e = &e_dart_a->second;
                    // assign face
                    dart_a->f = &faceVec[i];
                    // assign volume
                    dart_a->vo = &volume;
                    darts.emplace_back(dart_a);
                    //construct second dart
                    std::cout << "dart " << countDart << ": v" << face_indices[i][j]-1 << ", e" << edgeS << ", f" << i << "\n";
                    countDart++;
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
                    std::cout << "dart " << countDart << ": v" << face_indices[i][j+1]-1 << ", e" << edgeS << ", f" << i << "\n";
                }
            }
        }
        std::cout << "\n";
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

int main(int argc, const char * argv[]) {

    std::string file_in = "/Users/danieldobson/Library/CloudStorage/OneDrive-Personal/GEOMATICS/GEO1004/assignments/geo1004.2022/hw/01/data/torus.obj";
    std::string cube_test = "/Users/danieldobson/Library/CloudStorage/OneDrive-Personal/GEOMATICS/GEO1004/assignments/hw01/hw/01/data/cube.obj";
    std::string cube_test2 = "/Users/danieldobson/Library/CloudStorage/OneDrive-Personal/GEOMATICS/GEO1004/assignments/geo1004.2022/hw/01/data/cube2.obj";
    std::string file_out_obj = "/Users/danieldobson/Library/CloudStorage/OneDrive-Personal/GEOMATICS/GEO1004/assignments/geo1004.2022/hw/01/data/torus_triangulated.obj";
    std::string file_out_csv_d = "/Users/danieldobson/Library/CloudStorage/OneDrive-Personal/GEOMATICS/GEO1004/assignments/geo1004.2022/hw/01/data/torus_darts.csv";
    std::string file_out_csv_0 = "/Users/danieldobson/Library/CloudStorage/OneDrive-Personal/GEOMATICS/GEO1004/assignments/geo1004.2022/hw/01/data/torus_vertices.csv";
    std::string file_out_csv_1 = "/Users/danieldobson/Library/CloudStorage/OneDrive-Personal/GEOMATICS/GEO1004/assignments/geo1004.2022/hw/01/data/torus_edges.csv";
    std::string file_out_csv_2 = "/Users/danieldobson/Library/CloudStorage/OneDrive-Personal/GEOMATICS/GEO1004/assignments/geo1004.2022/hw/01/data/torus_faces.csv";
    std::string file_out_csv_3 = "/Users/danieldobson/Library/CloudStorage/OneDrive-Personal/GEOMATICS/GEO1004/assignments/geo1004.2022/hw/01/data/torus_volume.csv";

    // ## Read OBJ file ##
    // The vertices and faces are read and stored into vectors.

    std::vector<Vertex> vertices;
    std::vector<std::vector<int>> face_indices;

    // ## Construct generalised map using the structures from Gmap.h ##

    // where we store darts and cells
    std::vector<Dart*> darts;
    std::unordered_map<std::string, Vertex> vertexMap;
    std::unordered_map<std::string, Edge> edgeMap;
    std::vector<Face> faceVec;
    Volume volume;

    readObj(cube_test, vertices, face_indices);

    extractCells(vertices, face_indices, vertexMap, edgeMap, faceVec, volume);

    std::cout << "\nvertexMap.size() = " << vertexMap.size() << "\n";
    std::cout << "edgeMap.size() = " << edgeMap.size() << "\n";
    std::cout << "faceVec.size() = " << faceVec.size() << "\n";


    /*
     * code loop to visualise the edgeMap and its contents
    int countEdges = 0;
    std::unordered_map<std::string, Edge>:: iterator itrE;
    std::cout << "\nEdges : \n";
    for (itrE = edgeMap.begin(); itrE != edgeMap.end(); itrE++) {
        // itrE works as a pointer to pair<string, double>
        // type itrE->first stores the key part  and
        // itrE->second stores the value part
        std::cout << countEdges << " key: " << itrE->first << ", value: " << itrE->second.edgeS << std::endl;
        countEdges++;
    }
    */

    /*
     * code loop to visualise the contents of faceVec, more exact: to reveal how to access its vertices information
    std::cout << "\nfaces : \n";
    for (auto i:faceVec) {
        std::cout << "the face vertices: ";
        for (auto j:i.face_vertices) {
            std::cout << j << " ";
            if (j == i.face_vertices.back()) {std::cout << i.face_vertices.front();}
        }
        std::cout << "\n";
    }
    */

    std::cout << "\ndarts : \n";

    initCombiStruct(vertices,face_indices,edgeMap,faceVec,volume,darts);

    std::cout << "size of darts: " << darts.size() << "\n\n";

    // loop to visualise the contents of the combinatorial structure part of the gmap
    int countDarts = 0;
    for (auto i: darts) {
    countDarts++;
    //std::cout << "dart " << countDarts << ": \t" << i->v->point << "\t\t, e: \t" << i->e->edgeS << ", f: \t" << i->f << ", vo: \t" << i->vo << "\n";
    std::cout << "dart " << countDarts << "\t" << i << ": \t" << "a0: " << i->a0 << "\t\ta1: " << i->a1 << "\t\ta2: " << i->a2 << "\t\ta3: " << i->a3 << "\n";
    }

    int countVertices = 0;
    for (auto i: vertexMap) {
        countVertices++;
        std::cout << "Vertex " << countVertices << ":\tdart\t" << i.second.dart << "\n";
    }

    // ## Output generalised map to CSV ##



    // ## Create triangles from the darts ##

    // ## Write triangles to obj ##

    return 0;
}

