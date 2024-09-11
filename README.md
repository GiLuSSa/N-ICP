# N-ICP

This is a rough implementation of the Optimal Step Nonrigid ICP Algorithms  ([Amberg et al., 2007](https://doi.org/10.1109/CVPR.2007.383165)) for Rhino 8.

This implementation relies on CPython added in Grasshopper in Rhino 8 and uses only standard libraries such as NumPy and SciPy, as well as the Rhino Common API. To achieve this, the solver used in this implementation is scipy.sparse.linalg.lsmr, which requires the row-major flattening of the 4n x 3 unknown matrix X to express it into a 12n x 1 vector. Matrix A and vector B were adjusted accordingly. 
My code is not really pythonic; I apologise for this.

# Installation
No installation is required, only downloading the .gh file and have fun.

# Usage
The algorithm requires:
- a template (S) given as a triangular mesh;
- a target (T) surface given as a mesh, a point cloud or a brep;
- the landmarks (L) given as two separate lists, one containing the indices of vertices of the template and the other containing the corresponding point coordinates to aim at;
- a bunch of numerical parameters:
  - *Gamma* - can be used to weight differently the rotational and the skew part of the deformation against the translation part;
  - *Threshhold* - the maximum distance between a vertex of S and T to include that vertex in the distance term;
  - *maxAlpha* - maximum (and first) value for alpha; minAlpha, ninimum (and last) value of alpha;
  - *iterations* - number of values of alpha to iterate in the main loop;
  - *maxSubIterations* - maximum number of iterations of the inner loop;
  - *epsilon* - is used to break the sub-iteration before the maximum number of iterations if *epsilon < ||X<sub>n</sub> - X<sub>n-1</sub>||*.

The point cloud can be referred using the new (Rhino 8) native Point Cloud block or constructed with the Point Cloud Attributes block. Nor Volvox nor Cockroach point clouds can be directly used and must be converted. 
Landmarks can be defined in several ways and in the grasshopper definitions some alternatives are given.

The algorithm returns:
- the aligned template mesh;
- the displacement field;
- debug (currently the last set of closest points).

Due to the current limitations of Python in Grasshopper, and the unavailability of the terminal while a script is running, some information is also printed to a file (StatusDebug.txt) to monitor the registration. I assure that this is really handy for long registrations.

# Various
I've also added the .py file that's inside the .gh file such that it can be viewed by non-Rhino users. To implement the code I used Rhino Common but open3d can be used to deal with the meshed and scikit-learn can be used to find the nearest neighbors
between S(X) and T.

# Demo
Right now, I only provided a small demo to test the code. It's not really the use case for which the N-ICP was originally developed, but was the use case I (under the supervision of prof. Sinan Acikgoz) tried to use the algorithm on. 

The .3dm contains:
- a mesh representing the undeformed intrados surface of a barrel vault, used as a source mesh;
- the deformed intrados surface of the same barrel vault obtained with FE analysis, used as a target (T);
- the lines at the spring of the vault in the two configurations used to evaluate the landmarks.

My purpose (and hope) was to use this algorithm to reconstruct the displacement field associated with the vault damage mechanism. The results are better than what is possible to obtain with a C2M algorithm, but not satisfying. I tried to modify the N-ICP, but eventually decided to use another regularisation approach which ended in the development of the Pi-ICP I'll publish soon.
