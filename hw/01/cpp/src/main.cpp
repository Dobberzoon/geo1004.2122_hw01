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
        std::cout << "( ";


        // 0-cells
        for (int j = 0; j < face_indices[i].size(); j++) {
            std::cout << face_indices[i][j] << " ";

            // Initialise variables
            Vertex vertex_cur;
            std::string xyz;

            // Construct Vertex from current visiting point
            //std::cout << "vertices[j-1].point.x: " << vertices[j-1].point.x << "\n";
            vertex_cur = Vertex(vertices[face_indices[i][j]-1]);

            // For storing the cells, we use unordered_map, this will prevent multiple addition of same cells

            xyz = vertex_cur.xyz_tostring(vertex_cur.point.x,vertex_cur.point.y,vertex_cur.point.z);
            vertexMap.insert({xyz, vertex_cur});
            //std::cout << "vertex_cur: " << vertex_cur.point << "\n";
            //std::cout << "vertex_cur STRING: " << xyz << "\n";
            if (face_indices[i][j] == face_indices[i].back()) {std::cout << face_indices[i][0];}

            // 1-cells
            Edge edge_cur;
            std::string edgeS;

            if (face_indices[i][j] == face_indices[i].back()) {
                //std::cout << "This should be the last: " << face_indices[i][0] - 1  << "\n";
                //std::cout << "1st vertex: " << face_indices[i][j]-1 << " , second vertex: " << face_indices[i][0] - 1 << "\n";
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

                //std::cout << "edge in string: " << edgeS << "\n";
            }

            else {
                //std::cout << "1st vertex: " << face_indices[i][j]-1 << " , second vertex: " << face_indices[i][j+1]-1 << "\n";
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

            face_cur.face_vertices.push_back(face_indices[i][j]-1);

            if (j == (face_indices[i].size() - 1)) {
                //std::cout << "end is reached! \n";
                //std::cout << "size face_cur.vertices: " << face_cur.face_vertices.size() << "\n";
                //std::cout << "face_cur.vertices: " << face_cur.face_vertices[0] << ", " << face_cur.face_vertices[1] << ", " << face_cur.face_vertices[2] << ", " << face_cur.face_vertices[3] << "\n";
                faceVec.push_back(face_cur);
            }
        }
        std::cout << " ) \n";
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
  std::vector<Dart> darts;
  std::unordered_map<std::string, Vertex> vertexMap;
  std::unordered_map<std::string, Edge> edgeMap;
  std::vector<Face> faceVec;
  Volume volume;

  readObj(cube_test, vertices, face_indices);

  extractCells(vertices, face_indices, vertexMap, edgeMap, faceVec, volume);

  std::cout << "vertexMap.size() = " << vertexMap.size() << "\n";
  std::cout << "edgeMap.size() = " << edgeMap.size() << "\n";
  std::cout << "faceVec.size() = " << faceVec.size() << "\n";



  //    iterating over all value of vertexMap
  std::unordered_map<std::string, Vertex>:: iterator itrV;
  std::unordered_map<std::string, Edge>:: iterator itrE;
  std::cout << "\nAll Elements : \n";

  std::cout << "\nVertices : \n";
  for (itrV = vertexMap.begin(); itrV != vertexMap.end(); itrV++) {
      // itrV works as a pointer to pair<string, double>
      // type itrV->first stores the key part  and
      // itrV->second stores the value part
      std::cout << "key: " << itrV->first << ", value: " << itrV->second.point << std::endl;
  }

    std::cout << "\nEdges : \n";
  for (itrE = edgeMap.begin(); itrE != edgeMap.end(); itrE++) {
      // itrE works as a pointer to pair<string, double>
      // type itrE->first stores the key part  and
      // itrE->second stores the value part
      std::cout << "key: " << itrE->first << ", value: " << itrE->second.edgeS << std::endl;
  }

  std::cout << "\nfaces : \n";
  for (auto i:faceVec) {
      std::cout << "the face vertices: ";
      for (auto j:i.face_vertices) {
          std::cout << j << " ";
      }
      std::cout << "\n";
  }

  std::cout << "\ndarts : \n";

  for (auto i:faceVec) {
      std::cout << "the face vertices: ";
      for (auto j:i.face_vertices) {
          std::cout << j << " ";
      }
      std::cout << "\n";
  }
  

  
  // ## Output generalised map to CSV ##

  // ## Create triangles from the darts ##

  // ## Write triangles to obj ##
  
  return 0;
}

