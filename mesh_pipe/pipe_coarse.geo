D = 0.1;
size = 1e-0;

Point (1) = {5*D, 0.05*D + D/2, 0, size/50};
Point (2) = {5*D + D/2, 0.05*D + D/2, 0, size/50};
Point (3) = {5*D, 0.05*D + D, 0, size/50};
Point (4) = {5*D - D/2, 0.05*D + D/2, 0, size/50};
Point (5) = {5*D, 0.05*D, 0, size/50};
Point (6) = {0, 0, 0, size/25};
Point (7) = {5*D - D, 0, 0, size/50};
Point (8) = {5*D + D, 0, 0, size/50};
Point (9) = {30*D, 0, 0, size/25};
Point (10) = {30*D, 4*D, 0, size/10};
Point (11) = {0, 4*D, 0, size/10};

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
//Physical Line (3) = {7}; //bottom_left
Physical Line (4) = {7, 8, 9}; //bottom
//Physical Line (5) = {9}; //bottom_right
Physical Line (6) = {10}; //right
Physical Line (7) = {1, 2, 3, 4}; //circle

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
