# Triangulating a polygonal mesh with generalised maps

Assignment for geo1004.2022

Authors: 		      Daniël Dobson, Irina Gheorghiu, Leo Kan
Student number:   5152739, 5627834, 5505801

This program writes a mesh (.obj) into a datastructure (generalised map), and outputs a triangulated mesh (.obj). This is achieved in four steps:

1. read a 3D polygonal mesh from an OBJ file,
2. store the polygonal mesh in a generalised map,
3. output all the darts and cells from the generalised map to .csv files, and finally
4. output a triangulation of the polygonal mesh to a new OBJ file.

## How to build from command line

Below is an example of how to build the project on Unix based machine, however to clone you need to ask permission first to access this private repo @Dobberzoon on github.com.

Make sure you edit the CMakelists.txt, such that it finds the libraries from their relative paths accordingly.

```
$ git clone https://github.com/Dobberzoon/geo1004.2122_hw01
$ cd geo1004.2122_hw01/hw/01/cpp
$ mkdir build
$ cd build
$ cmake ..
$ make
```

Or with zip file, extract and:

```
$ cd geo1004.2122_hw01/hw/01/cpp
$ mkdir build
$ cd build
$ cmake ..
$ make
```

Now, to start the program from the build folder in cmd line:

```
$ ./hw01
```


## Alternative build: simply use IDE of choice (e.g. CLion 2021.3)


## Add own path for input dataset
Change path for input_data location at line 810-817 to link your local path accordingly.

## Provided dataset

cube.obj -> for testing purposes
torus.obj 
