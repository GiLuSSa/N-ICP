# N-ICP

This is a rough implementation of the Optimal Step Nonrigid ICP Algorithms  ([Amberg et al., 2007](https://doi.org/10.1109/CVPR.2007.383165)) for Rhino 8.

This implementation relies on CPython recently added in Grasshopper and uses only standard libraries like numpy and scipy, as well as RhinoCommon API. To achieve this, the solver used in this implementation is scipy.sparse.linalg.lsmr that, to be used, required the row-major flattening of the 4n x 3 unknown matrix X to express it into a 12n x 1 vector. 

# Installation
No installation is needed, just download the .gh file and have fun.

# Usage
The algorithm requires:
- a template (S) give as a triangular mesh;
- a target (T) surface given as a mesh, a point cloud or a brep;
- the optional landmarks (L) given as two separate lists, one containing the indexes of verteces of the template and the other containing the corresponding points coordinates to aim at;
- a bunch of numerical parameters that will be explained later.

The point cloud can be referred using the new native Point Cloud block or constructed with the Point Cloud Attributes block. Nor Volvox nor Cockroach point clouds can be directly used and must be converted. 
Landmarks can be defined in several ways and in the grasshopper definition some alternatives are given.

The algorithm returns:
- the aligned template mesh;
- the transformation matrix of each point as a 4x4 rhino affine transformations 
- the displacement field.
- debug.

Some post processing is provided to better understand what's happened.

Due to the current limitations of Python in Grasshopper, and the unavailability of the console while a script is running, some debug information are printed to a file (StatusDebug.txt) to monitor the registration. I assure that this is really handy for long registrations.

# Various
I've also added the .py file that's inside the .gh file so that it can be used by non Rhino users. To implement the code I used RhinoCommon but open3d can be used to deal with the meshed and sklearn can be used to find the nearest neighbors between S(X) and T

