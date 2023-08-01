//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {0.04, 0, 0, 1.0};
//+
Point(3) = {0.08, 0, 0, 1.0};
//+
Point(4) = {0.12, 0, 0, 1.0};
//+
Point(5) = {0.16, 0, 0, 1.0};
//+
Point(6) = {0.2, 0, 0, 1.0};
//+
Point(7) = {0.2, 0.04, 0, 1.0};
//+
Point(8) = {0.16, 0.04, 0, 1.0};
//+
Point(9) = {0.12, 0.04, 0, 1.0};
//+
Point(10) = {0.08, 0.04, 0, 1.0};
//+
Point(11) = {0.04, 0.04, 0, 1.0};
//+
Point(12) = {0, 0.04, 0, 1.0};
//+
Point(13) = {0.02, 0.02, 0, 1.0};
//+
Point(14) = {0.03, 0.02, 0, 0.1};
//+
Point(15) = {0.025, 0.028, 0, 0.1};
//+
Point(16) = {0.015, 0.028, 0, 0.1};
//+
Point(17) = {0.01, 0.02, 0, 0.1};
//+
Point(18) = {0.015, 0.012, 0, 0.1};
//+
Point(19) = {0.025, 0.012, 0, 0.1};
//+
Point(20) = {0.1, 0.02, 0, 1.0};
//+
Point(21) = {0.11, 0.02, 0, 0.1};
//+
Point(22) = {0.105, 0.028, 0, 0.1};
//+
Point(23) = {0.095, 0.028, 0, 0.1};
//+
Point(24) = {0.09, 0.02, 0, 0.1};
//+
Point(25) = {0.095, 0.012, 0, 0.1};
//+
Point(26) = {0.105, 0.012, 0, 0.1};
//+
Point(27) = {0.18, 0.02, 0, 1.0};
//+
Point(28) = {0.19, 0.02, 0, 0.1};
//+
Point(29) = {0.185, 0.028, 0, 0.1};
//+
Point(30) = {0.175, 0.028, 0, 0.1};
//+
Point(31) = {0.17, 0.02, 0, 0.1};
//+
Point(32) = {0.175, 0.012, 0, 0.1};
//+
Point(33) = {0.185, 0.012, 0, 0.1};

//+
Circle(1) = {1, 13, 2};
//+
Line(2) = {2, 3};
//+
Circle(3) = {3, 20, 4};
//+
Line(4) = {4, 5};
//+
Circle(5) = {5, 27, 6};
//+
Circle(6) = {6, 27, 7};
//+
Circle(7) = {7, 27, 8};
//+
Line(8) = {8, 9};
//+
Circle(9) = {9, 20, 10};
//+
Line(10) = {10, 11};
//+
Circle(11) = {11, 13, 12};
//+
Circle(12) = {12, 13, 1};

//+
Line(13) = {14, 15};
//+
Line(14) = {15, 16};
//+
Line(15) = {16, 17};
//+
Line(16) = {17, 18};
//+
Line(17) = {18, 19};
//+
Line(18) = {19, 14};

//+
Line(19) = {21, 22};
//+
Line(20) = {22, 23};
//+
Line(21) = {23, 24};
//+
Line(22) = {24, 25};
//+
Line(23) = {25, 26};
//+
Line(24) = {26, 21};

//+
Line(25) = {28, 29};
//+
Line(26) = {29, 30};
//+
Line(27) = {30, 31};
//+
Line(28) = {31, 32};
//+
Line(29) = {32, 33};
//+
Line(30) = {33, 28};
//+
Curve Loop(1) = {12, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
//+
Curve Loop(2) = {13, 14, 15, 16, 17, 18};
//+
Curve Loop(3) = {19, 20, 21, 22, 23, 24};
//+
Curve Loop(4) = {25, 26, 27, 28, 29, 30};
//+
Plane Surface(1) = {1, 2, 3, 4};
//+
Physical Curve("FIXED", 31) = {13, 14, 15, 16, 17, 18};
