D = 0.1;
size = 5e-3;

Point (1) = {10*D, 0.05*D + D/2, 0, size};
Point (2) = {10*D + D/2, 0.05*D + D/2, 0, size};
Point (3) = {10*D, 0.05*D + D, 0, size};
Point (4) = {10*D - D/2, 0.05*D + D/2, 0, size};
Point (5) = {10*D, 0.05*D, 0, size/4};
Point (6) = {0, 0, 0, size};
Point (7) = {10*D - D/2, 0, 0, size/4};
Point (8) = {10*D + D/2, 0, 0, size/4};
Point (9) = {30*D, 0, 0, size};
Point (10) = {30*D, 4*D, 0, size*10};
Point (11) = {0, 4*D, 0, size*10};

Circle (1) = {2, 1, 3};
Circle (2) = {3, 1, 4};
Circle (3) = {4, 1, 5};
Circle (4) = {5, 1, 2};
Line (5) = {6, 11};
Line (6) = {11, 10};
Line (7) = {6, 7};
Line (8) = {7,8};
Line (9) = {8, 9};
Line (10) = {9, 10};

Physical Line (1) = {5}; //left
Physical Line (2) = {6}; //top
Physical Line (3) = {7,8,9}; //bottom
Physical Line (4) = {10}; //right
Physical Line (5) = {1, 2, 3, 4}; //circle

Line Loop (11) = {1, 2, 3, 4};
Line Loop (12) = {-5, -6, 7, 8, 9, 10};
Plane Surface (13) = {12, 11};
Physical Surface (14) = {13};

Field[1]=Attractor;
Field[1].EdgesList={3,4,8};

Field[2]=MathEval;
Field[2].F="0.005+0.5*F1";
Background Field=2;

View "comments" {T2(10, 15, 0){ StrCat("File created on ", Today) };};
