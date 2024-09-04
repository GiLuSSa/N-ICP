# N-ICP
N-ICP for Rhino 8

This is a rough implementation of the Optimal Step Nonrigid ICP Algorithms of Amberg et al. (2007, https://doi.org/10.1109/CVPR.2007.383165) for Rhino 8.
This implementation relies on CPython recently added in Grasshopper and uses only standard libraries like numpy and scipy, as well as rhinoscriptsyntax and Rhino. To achieve this, the solver used in this implementation is scipy.sparse.linalg.lsmr that, to be used, required the row-major flattening of the 4n x 3 unknown matrix X to express it into a 12n x 1 vector. 
