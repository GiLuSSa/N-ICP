# N-ICP

This is a rough implementation of the Optimal Step Nonrigid ICP Algorithms  ([Amberg et al., 2007](https://doi.org/10.1109/CVPR.2007.383165)) for Rhino 8.

This implementation relies on CPython added in Grasshopper in Rhino 8 and uses only standard libraries like numpy and scipy, as well as RhinoCommon API. To achieve this, the solver used in this implementation is scipy.sparse.linalg.lsmr that, to be used, required the row-major flattening of the 4n x 3 unknown matrix X to express it into a 12n x 1 vector. 
My code is not really pythonic, I'm sorry for that.

# Installation
No installation is needed, just download the .gh file and have fun.

# Usage
The algorithm requires:
- a template (S) given as a triangular mesh;
- a target (T) surface given as a mesh, a point cloud or a brep;
- the landmarks (L) given as two separate lists, one containing the indexes of verteces of the template and the other containing the corresponding points coordinates to aim at;
- a bunch of numerical parameters that will be explained later.

The point cloud can be referred using the new (Rhino 8) native Point Cloud block or constructed with the Point Cloud Attributes block. Nor Volvox nor Cockroach point clouds can be directly used and must be converted. 
Landmarks can be defined in several ways and in the grasshopper definitions some alternatives are given.

The algorithm returns:
- the aligned template mesh;
- the displacement field;
- debug (currently the last set of closest points).

Due to the current limitations of Python in Grasshopper, and the unavailability of the terminal while a script is running, some information are printed also to a file (StatusDebug.txt) to monitor the registration. I assure that this is really handy for long registrations.

# Various
I've also added the .py file that's inside the .gh file so that it can be viewed by non Rhino users. To implement the code I used RhinoCommon but open3d can be used to deal with the meshed and sklearn can be used to find the nearest neighbors between S(X) and T.

# Demo
Right now I only provided a small demo to test the code. It's not really the usecase for which the N-ICP was originally developet but was reason I (under the supervision of prof. Sinan Acikgoz) tried to use the algorithm and the reason why I also output the displacement field. 

The .3dm contains:
- a mesh representing the undeformed untradox surface of a barrel vault used as a source mesh;
- the deformed intradox surface of the same barrel vault obtained with FE analysis used as a target (T);
- the lines at the spring of the vaul in the two configuration that are used to evaluate the landmarks.

My purpose (and hope) was to use this algorithm to reconstruct the displacement field corresponding to the damage mechanism of the vault. The results are better than what it's possible to obtain with just a C2M algorithm but not satisfying, hence i tried to modify the N-ICP but eventually I decided to use another regularization approach wich ended in the development of the Pi-ICP I'll publish soon.
