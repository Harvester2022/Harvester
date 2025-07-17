// Gmsh project created on Fri Jan 24 17:11:53 2025
SetFactory("OpenCASCADE");

fix = 8;   //Distancia fija a la plataforma
Lt  = 60;  //Longitud total Fe frame
w1  = 0.4; //Grosor frame

xb = 12; //Posicion lado izquierdo bobina. 
d  = 10; //longitud bobina. 
wc = 1;  //Grosor Bobina

wm2 = 5; //Lados iman2
x2  = 5; //posicion centro iman 2

wm1 = 5;     //Lados iman1 
x1  = 35;    //posicion centro iman 1

xGaFe = xb-2.5; //Posicion lado izquierdo lamina GaFe, se puede poner en funcion de xb, la de la bobina
dGaFe = 15;  //LongitudlaminaGaFE 
wGaFe = 0.5; // Grosor lamina GaFe

//+
Point(1) = {-fix , 0  , 0 , 1.0};
Point(2) = {0    , 0  , 0 ,1.0};
Point(3) = {xb   , 0  , 0 ,1.0};
Point(4) = {xb+d , 0  , 0 ,1.0};
Point(5) = {Lt-fix , 0  , 0 ,1.0};
Point(6) = {Lt-fix , w1  , 0 ,1.0};
Point(7) = {x2+wm2/2 , w1  , 0 ,w1/6};
Point(8) = {x2-wm2/2 , w1  , 0 ,w1/6};
Point(9) = {xGaFe   , w1 , 0 ,1.0};
Point(10) = {xGaFe+dGaFe , w1  , 0 ,1.0};
Point(11) = {x1+wm1/2 , w1  , 0 ,w1/6};
Point(12) = {x1-wm1/2 , w1  , 0 ,w1/6};
Point(13) = {-fix , w1  , 0 ,1.0};

Point(14) = {x1+wm1/2 , w1+wm1  , 0 ,1.0};
Point(15) = {x1-wm1/2 , w1+wm1  , 0 ,1.0};

Point(16) = {x2+wm2/2 , w1+wm2  , 0 ,1.0};
Point(17) = {x2-wm2/2 , w1+wm2  , 0 ,1.0};

Point(21) = {xb   , w1 + wGaFe  , 0 ,1.0};
Point(22) = {xb , w1 + wGaFe + wc , 0 ,1.0};
Point(23) = {xb+d   , w1 + wGaFe  , 0 ,1.0};
Point(24) = {xb+d , w1 + wGaFe + wc , 0 ,1.0};

Point(25) = {xb     , -wc  , 0 ,1.0};
Point(26) = {xb+d   , -wc  , 0 ,1.0};

Point(27) = {xGaFe    , w1 + wGaFe  , 0 ,1.0};
Point(28) = {xGaFe + dGaFe    , w1 + wGaFe  , 0 ,1.0};

//+
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 11};
Line(7) = {12, 10};
Line(8) = {10, 9};
Line(9) = {9, 7};
Line(10) = {7, 8};
Line(11) = {8, 13};
Line(12) = {13, 1};
Line(13) = {7, 16};
Line(14) = {16, 17};
Line(15) = {17, 8};
Line(16) = {11, 14};
Line(17) = {14, 15};
Line(18) = {15, 12};
Line(19) = {28, 23};
Line(20) = {23, 21};
Line(21) = {27, 9};
Line(22) = {3, 25};
Line(23) = {25, 26};
Line(24) = {26, 4};
Line(25) = {23, 24};
Line(26) = {24, 22};
Line(27) = {22, 21};
Line(29) = {12, 11};
Line(30) = {22, 21};
Line(31) = {12, 11};
Line(32) = {27, 21};
Line(33) = {28, 10};


//+
Circle(28) = {0, 0, 0, Lt*5, 0, 2*Pi};
//+
//+
Curve Loop(13) = {28};
Curve Loop(14) = {1, 2, 22, 23, 24, 4, 5, 6, 16, 17, 18, 7, -33, 19, 25, 26, 27, -32, 21, 9, 13, 14, 15, 11, 12};
Plane Surface(1) = {13, 14};


//+
Curve Loop(3) = {11, 12, 1, 2, 3, 4, 5, 6, -29, 7, 8, 9, 10};
Plane Surface(2) = {3};
//+
Curve Loop(4) = {15, -10, 13, 14};
Plane Surface(3) = {4};
//+
Curve Loop(12) = {18, 29, 16, 17};
Plane Surface(4) = {12};
//+
Curve Loop(9) = {20, -32, 21, -8, 33, 19};
Plane Surface(5) = {9};
//+
Curve Loop(7) = {20, -27, -26, -25};
Plane Surface(6) = {7};
//+
Curve Loop(8) = {3, -24, -23, -22};
Plane Surface(7) = {8};

//ElementSizes
Characteristic Length{ PointsOf{ Curve{28};} } = Lt;
Characteristic Length{ PointsOf{ Surface{2};Line{29};Line{10}; } } = w1/6;
Characteristic Length{ PointsOf{ Surface{5};Surface{6};Surface{7};} } = w1/4;
Characteristic Length{ PointsOf{ Surface{3};Surface{4};} } = wm1/8;


//+
Physical Curve("infty", 1) = {28};
Physical Curve("fix  ", 2) = {1};
Physical Surface("aire", 1) = {1};
Physical Surface("Frame", 2) = {2};
Physical Surface("mag1", 3) = {3};
Physical Surface("mag2", 4) = {4};
Physical Surface("GaFe", 5) = {5};
Physical Surface("Coil+", 6) = {6};
Physical Surface("Coil-", 7) = {7};

