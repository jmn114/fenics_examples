// GMSH allows you to define variables
l = 0.1;
lsmall = 0.01;
// Points are defined by 3 coordinates
// The 4-th entry describes the required
// element size around this point
Point(1) = {0, 0, 0, l};
Point(2) = {0, 1, 0, l};
Point(3) = {1, 1, 0, l};
Point(4) = {1, 0, 0, l};
Point(5) = {0.2, 1, 0, lsmall};
Point(6) = {0.8, 1, 0, lsmall};
Point(7) = {0.1, 1, 0, l};
Point(8) = {0.9, 1, 0, l};
Point(9) = {-0, 0.7, 0, l};
Point(10) = {1, 0.7, 0, l};
Line(1) = {2, 7};
Line(2) = {7, 5};
Line(3) = {5, 6};
Line(4) = {6, 8};
Line(5) = {8, 3};
Line(6) = {3, 10};
Line(7) = {10, 4};
Line(8) = {4, 1};
Line(9) = {1, 9};
Line(10) = {9, 2};
Line(11) = {9, 10};
Line Loop(12) = {1, 2, 3, 4, 5, 6, -11, 10};
Plane Surface(13) = {12};
Line Loop(14) = {11, 7, 8, 9};
Plane Surface(15) = {14};
// These Physical Lines will be used to define
// the boundaries; if you start counting from 1
// then all internal boundaries will be labeled 0
Physical Line(1) = {10, 1, 9, 8, 7, 6, 5};
Physical Line(2) = {3};
Physical Line(3) = {2, 4};
// Physical Surfaces define regions
// I usually start counting from 0
// You need at least 2 or dolfin-convert ignores these
Physical Surface(0) = {13};
Physical Surface(1) = {15};
