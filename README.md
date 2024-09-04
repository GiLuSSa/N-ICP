# N-ICP

This is a rough implementation of the Optimal Step Nonrigid ICP Algorithms  ([Amberg et al., 2007](https://doi.org/10.1109/CVPR.2007.383165)) for Rhino 8.

This implementation relies on CPython recently added in Grasshopper and uses only standard libraries like numpy and scipy, as well as rhinoscriptsyntax and Rhino. To achieve this, the solver used in this implementation is scipy.sparse.linalg.lsmr that, to be used, required the row-major flattening of the 4n x 3 unknown matrix X to express it into a 12n x 1 vector. 

# Installation
No installation is needed, just download the .gh file and have fun.

# Usage
The algorithm requires:
- a template (S) give as a triangular mesh;
- a target (T) surface given as a mesh, a point cloud or a brep;
- the optional landmarks (L) given as two separate lists of points;
- a bunch of numerical parameters that will be explained later.

The point cloud can be referred using the new native Point Cloud block or constructed with the Point Cloud Attributes block. Nor Volvox nor Cockroach point clouds can be directly used and must be converted. 



