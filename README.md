Converts a tetrahedral volumetric mesh in Gmsh's .msh file format (version
4) to a surface mesh in .stl format.

This tool depends upon:
```
CLI11
GMsh
```

### WARNING

The GMsh developers don't integrate GMsh with CMake. As such, a custom CMake
Find module must be used to link against GMsh.

The build system does not offer installation. As such, `cmake --install .`
should not be used. 

### Building with CMake.
```
mkdir build/
cd build/
cmake -D CMAKE_MODULE_PATH=/path/to/your/gmsh/find/module/ ..
ccmake .
cmake --build .
```
