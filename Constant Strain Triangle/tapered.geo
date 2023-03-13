//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {0.1, 0, 0, 1.0};
//+
Point(3) = {0.175, 0.025, 0, 1.0};
//+
Point(4) = {0.225, 0.025, 0, 1.0};
//+
Point(5) = {0.225, 0.075, 0, 1.0};
//+
Point(6) = {0.175, 0.075, 0, 1.0};
//+
Point(7) = {0.1, 0.1, 0, 1.0};
//+
Point(8) = {0, 0.1, 0, 1.0};
//+
Point(9) = {0.1, 0.05, 0, 1.0};
//+
Point(10) = {0.1125, 0.05, 0, 0.1};
//+
Point(11) = {0.1, 0.0625, 0, 0.1};
//+
Point(12) = {0.0875, 0.05, 0, 0.1};
//+
Point(13) = {0.1, 0.0375, 0, 0.1};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 5};
//+
Line(5) = {5, 6};
//+
Line(6) = {6, 7};
//+
Line(7) = {7, 8};
//+
Line(8) = {8, 1};
//+
Circle(9) = {10, 9, 11};
//+
Circle(10) = {11, 9, 12};
//+
Circle(11) = {12, 9, 13};
//+
Circle(12) = {13, 9, 10};
//+
Curve Loop(1) = {9, 10, 11, 12};
//+
Curve Loop(2) = {7, 8, 1, 2, 3, 4, 5, 6};
//+
Plane Surface(1) = {1, 2};
//+
Physical Curve("FIXED", 13) = {8};