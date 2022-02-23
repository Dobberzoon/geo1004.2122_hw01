// some standard libraries that are helpfull for reading/writing text files
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>

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
          else vertices.push_back(Vertex());
        }

        if (word == "f") {
          std::vector<int> face;
          while (iss >> word) face.push_back(std::stoi(word));
          face_indices.push_back(face);
        }

      }
    }
    stream_in.close();

    std::vector<Dart> darts;
    std::vector<Vertex> vertexVec;
    std::vector<Edge> edgeVec;
    std::vector<Face> faceVec;

    int count = 0;

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


  
  // ## Construct generalised map using the structures from Gmap.h ##
  
  // ## Output generalised map to CSV ##

  // ## Create triangles from the darts ##

  // ## Write triangles to obj ##
  
  return 0;
}
