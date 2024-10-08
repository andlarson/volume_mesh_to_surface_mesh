Converts a tetrahedral volumetric mesh in Gmsh's .msh file format (version
4) to a surface mesh in .stl format.

This tool depends upon the CLI11 and Gmsh libraries.

### WARNING
This tool has only been tested on macOS.

Furthermore, it is necessary to write your own CMake Find module to use this
tool because Gmsh is not integrated with Gmsh.

The build system does not offer installation. As such, `cmake --install .`
should not be issued. The binary can be used in the created `build/` directory.

### Building with CMake.
```
mkdir build/
cd build/
cmake -D CMAKE_FIND_MODULE=/path/to/your/gmsh/find/module/ ..
ccmake .
cmake --build .
```
