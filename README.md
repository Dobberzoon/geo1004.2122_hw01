# Geo1004.2022

Assignments for geo1004.2022

Authors: 		Daniël Dobson, Irina, Leo
Student number: 5152739

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

Runtime on Apple Mac Air M1:  

## Alternative build: simply use IDE of choice (e.g. CLion 2021.3)

Runtime on Apple Mac Air M1: 133ms

## How to use


## Add own path for input dataset
Change path for input_data location at line 313-321 to link your local path accordingly.

## Provided dataset

cube.obj -> for testing purposes
torus.obj 
