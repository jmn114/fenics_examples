edgeLength = 0.5;

Point(1) = {0.0, 0.0, 0.0, edgeLength};

Extrude {1, 0.0, 0.0} {
  Point{1}; Layers{1/edgeLength};
}

Extrude {0.0, 1, 0.0} {
  Line{1}; Layers{1/edgeLength};
}

Physical Line(1) = {1}; //bottom 
Physical Line(2) = {2}; //top
Physical Line(3) = {3}; //inlet
Physical Line(4) = {4}; //outlet
Physical Surface(5) = {5};