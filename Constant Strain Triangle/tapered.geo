//+
Point(1) = {0, -0, 0, 1.0};
//+
Point(2) = {0, 0.1, 0, 1.0};
//+
Point(3) = {0.1, 0.1, 0, 1.0};
//+
Point(4) = {0.1, -0, 0, 1.0};
//+
Point(5) = {0.2, 0.075, 0, 1.0};
//+
Point(6) = {0.2, 0.025, 0, 1.0};
//+
Point(7) = {0.25, 0.075, 0, 1.0};
//+
Point(8) = {0.25, 0.025, 0, 1.0};
//+
Point(9) = {0.1, 0.05, 0, 1.0};
//+
Point(10) = {0.12, 0.05, 0, 0.1};
//+
Point(11) = {0.1, 0.07, 0, 0.1};
//+
Point(12) = {0.08, 0.05, 0, 0.1};
//+
Point(13) = {0.1, 0.03, 0, 0.1};
//+
Line(1) = {1, 4};
//+
Line(2) = {4, 6};
//+
Line(3) = {6, 8};
//+
Line(4) = {8, 7};
//+
Line(5) = {7, 5};
//+
Line(6) = {5, 3};
//+
Line(7) = {3, 2};
//+
Line(8) = {2, 1};
//+
Circle(9) = {10, 9, 11};
//+
Circle(10) = {11, 9, 12};
//+
Circle(11) = {12, 9, 13};
//+
Circle(12) = {13, 9, 10};
//+
Curve Loop(1) = {10, 11, 12, 9};
//+
Curve Loop(2) = {7, 8, 1, 2, 3, 4, 5, 6};
//+
Plane Surface(1) = {1, 2};
//+
Physical Curve("FIXED", 13) = {8};
