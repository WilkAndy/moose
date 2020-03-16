# GeneratedMeshGenerator

## Description

The `GeneratedMeshGenerator` object is the built-in mesh generation capable of creating points, lines, rectangles, and rectangular
prisms ("boxes"). The mesh automatically creates side sets that are logically named and numbered as follows:

- In 0D, the point belongs to a node set with ID = 0
- In 1D, left = 0, right = 1
- In 2D, bottom = 0, right = 1, top = 2, left = 3
- In 3D, back = 0, bottom = 1, right = 2, top = 3, left = 4, front = 5

The length, width, and height of the domain, as well as the number of elements in each direction can be specified
independently.

In the zero-dimensional case, one node is created at position `(xmin, ymin, zmin)`, and one EDGE2 element is created that connects this node with itself.  The node belongs to a node set with ID = 0, and a side set with ID = 0.  Both boundary sets are called "the_point".

## Further GeneratedMeshGenerator Documentation

!syntax parameters /MeshGenerators/GeneratedMeshGenerator

!syntax inputs /MeshGenerators/GeneratedMeshGenerator

!syntax children /MeshGenerators/GeneratedMeshGenerator
