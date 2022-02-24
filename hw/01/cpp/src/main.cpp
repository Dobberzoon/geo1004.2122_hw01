// some standard libraries that are helpfull for reading/writing text files
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <unordered_map>

#include "Gmap.h"

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

    std::ifstream stream_in;
    stream_in.open(cube_test);
    std::vector<Vertex> vertices;
    std::vector<std::vector<int>> face_indices;
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

    std::vector<Dart> darts;
    std::unordered_map<std::string, Vertex> vertexMap;
    std::vector<Edge> edgeVec;
    std::vector<Face> faceVec;

    /*
    int count = 1;

    for (auto i : vertices) {
        std::cout << "point " << count << ": " << i.point << "\n";
        count++;
    }

    for (auto i : face_indices) {
        std::cout << "( ";
        for (auto j: i) {
            std::cout << j << " ";
            if (j == i.back()) {std::cout << i.front() << " ";}

        }
        std::cout << ") \n";

        for (auto j: i) {
            std::cout << j << " " << vertices[j-1].point << "\n";
        }
    }
    */

    // This loop traverses all faces (per indices), and in the double loop we traverse the vertices
    // that make up each face.

    for (int i = 0; i < face_indices.size(); i++) {
        // the std::cout's are only for visualising the loop process
        std::cout << "( ";

        for (int j = 0; j < face_indices[i].size(); j++) {
            std::cout << face_indices[i][j] << " ";

            // Initialise variables
            Vertex vertex_cur;
            std::string xyz;

            // Construct Vertex from current visiting point
            //std::cout << "vertices[j-1].point.x: " << vertices[j-1].point.x << "\n";
            vertex_cur = Vertex(vertices[face_indices[i][j]-1]);

            // For storing the cells, we use unordered_map, this will prevent multiple addition of same cells
            // 0-cells
            xyz = vertex_cur.xyz_tostring(vertex_cur.point.x,vertex_cur.point.y,vertex_cur.point.z);
            vertexMap.insert({xyz, vertex_cur});
            //std::cout << "vertex_cur: " << vertex_cur.point << "\n";
            //std::cout << "vertex_cur STRING: " << xyz << "\n";

            // 1-cells
            if (face_indices[i][j] == face_indices[i].back()) {std::cout << face_indices[i][0];}

        }
        std::cout << ") \n";

        /*
        for (int j = 0; j < face_indices[i].size(); j++) {
            std::cout << "point " << face_indices[i][j] << ": " << vertices[face_indices[i][j]-1].point << "\n";
        }
        std::cout << "\n";
        */
    }
    std::cout << "vertexMap.size() = " << vertexMap.size();


  
  // ## Construct generalised map using the structures from Gmap.h ##
  
  // ## Output generalised map to CSV ##

  // ## Create triangles from the darts ##

  // ## Write triangles to obj ##
  
  return 0;
}
