Converts a tetrahedral volumetric mesh in Gmsh's .msh file format (version
4) to a surface mesh in .stl format.

This tool depends upon the CLI11 and Gmsh libraries.

### WARNING

This tool has only been tested on macOS and currently includes some machine-specific 
behavior in the build system. As such, if you want to use this tool, you *must* modify
the CMakeLists.txt file.

### Building with CMake.

```
mkdir build/
cd build/
cmake ..
cmake --build .
```
