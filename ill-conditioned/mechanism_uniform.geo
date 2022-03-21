// This code was created by pygmsh vunknown.
p0 = newp;
Point(p0) = {0.0, 0.0, 0.0, 1.0};
p1 = newp;
Point(p1) = {120, 0.0, 0.0, 1.0};
p2 = newp;
Point(p2) = {120, 50, 0.0, 1.0};
p3 = newp;
Point(p3) = {120.0, 60.0, 0.0, 1.0};
p4 = newp;
Point(p4) = {0.0, 60.0, 0.0, 0.08};
p5 = newp;
Point(p5) = {0.0, 50.0, 0.0, 0.08};
p6 = newp;
Point(p6) = {0.0, 2.0, 0.0, 1.0};
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
Line(7) = {7, 1};
//+
Curve Loop(1) = {4, 5, 6, 7, 1, 2, 3};
//+
Plane Surface(1) = {1};
//+
Transfinite Surface {1} = {5, 1, 2, 4};
//+
Transfinite Curve {4, 1} = 61 Using Progression 1;
//+
Transfinite Curve {2} = 26 Using Progression 1;
//+
Transfinite Curve {5, 3} = 6 Using Progression 1;
//+
Transfinite Curve {7} = 2 Using Progression 1;
//+
Transfinite Curve {6} = 25 Using Progression 1;
//+
Physical Curve(1) = {7};
//+
Physical Curve(2) = {5};
//+
Physical Curve(3) = {3};
//+
Physical Curve(4) = {4};
//+
Physical Surface(1) = {1};
