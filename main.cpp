#include <iostream>
#include <eigen3/Eigen/Dense>
#include <cmath>
using Eigen::Matrix4d;
using Eigen::MatrixXd;
using Eigen::Matrix3d;
using Eigen::Vector3d;
using Eigen::VectorXd;    
#define LENGTH_OF_LEGPART 0.21
#define TRIANGLE_SIDE 0.4
/*T_base3 =
 
[cos(LegA(0))*cos(LegA(2))*sin(LegA(1)) - cos(LegA(0))*cos(LegA(1))*sin(LegA(2)), - cos(LegA(0))*cos(LegA(1))*cos(LegA(2)) - cos(LegA(0))*sin(LegA(1))*sin(LegA(2)),  sin(LegA(0)), LENGTH_OF_LEGPART*cos(LegA(0))*cos(LegA(1))*sin(LegA(2)) - LENGTH_OF_LEGPART*cos(LegA(0))*sin(LegA(1)) - (sqrt(3)*TRIANGLE_SIDE)/3 - LENGTH_OF_LEGPART*cos(LegA(0))*cos(LegA(2))*sin(LegA(1))]
[cos(LegA(2))*sin(LegA(0))*sin(LegA(1)) - cos(LegA(1))*sin(LegA(0))*sin(LegA(2)), - sin(LegA(0))*sin(LegA(1))*sin(LegA(2)) - cos(LegA(1))*cos(LegA(2))*sin(LegA(0)), -cos(LegA(0)),                 LENGTH_OF_LEGPART*cos(LegA(1))*sin(LegA(0))*sin(LegA(2)) - LENGTH_OF_LEGPART*sin(LegA(0))*sin(LegA(1)) - LENGTH_OF_LEGPART*cos(LegA(2))*sin(LegA(0))*sin(LegA(1))]
[                cos(LegA(1))*cos(LegA(2)) + sin(LegA(1))*sin(LegA(2)),                   cos(LegA(2))*sin(LegA(1)) - cos(LegA(1))*sin(LegA(2)),        0,                                   - LENGTH_OF_LEGPART - LENGTH_OF_LEGPART*cos(LegA(1)) - LENGTH_OF_LEGPART*cos(LegA(1))*cos(LegA(2)) - LENGTH_OF_LEGPART*sin(LegA(1))*sin(LegA(2))]
[                                                0,                                                   0,        0,                                                                                         1]
 
 
T_base7 =
 
[cos(LegB(2))*sin(LegB(1))*(cos(LegB(0))/2 - (sqrt(3)*sin(LegB(0)))/2) - cos(LegB(1))*sin(LegB(2))*(cos(LegB(0))/2 - (sqrt(3)*sin(LegB(0)))/2), - cos(LegB(1))*cos(LegB(2))*(cos(LegB(0))/2 - (sqrt(3)*sin(LegB(0)))/2) - sin(LegB(1))*sin(LegB(2))*(cos(LegB(0))/2 - (sqrt(3)*sin(LegB(0)))/2), sin(LegB(0))/2 + (sqrt(3)*cos(LegB(0)))/2, (sqrt(3)*TRIANGLE_SIDE)/6 - LENGTH_OF_LEGPART*sin(LegB(1))*(cos(LegB(0))/2 - (sqrt(3)*sin(LegB(0)))/2) + LENGTH_OF_LEGPART*cos(LegB(1))*sin(LegB(2))*(cos(LegB(0))/2 - (sqrt(3)*sin(LegB(0)))/2) - LENGTH_OF_LEGPART*cos(LegB(2))*sin(LegB(1))*(cos(LegB(0))/2 - (sqrt(3)*sin(LegB(0)))/2)]
[cos(LegB(2))*sin(LegB(1))*(sin(LegB(0))/2 + (sqrt(3)*cos(LegB(0)))/2) - cos(LegB(1))*sin(LegB(2))*(sin(LegB(0))/2 + (sqrt(3)*cos(LegB(0)))/2), - cos(LegB(1))*cos(LegB(2))*(sin(LegB(0))/2 + (sqrt(3)*cos(LegB(0)))/2) - sin(LegB(1))*sin(LegB(2))*(sin(LegB(0))/2 + (sqrt(3)*cos(LegB(0)))/2), (sqrt(3)*sin(LegB(0)))/2 - cos(LegB(0))/2,           TRIANGLE_SIDE/2 - LENGTH_OF_LEGPART*sin(LegB(1))*(sin(LegB(0))/2 + (sqrt(3)*cos(LegB(0)))/2) + LENGTH_OF_LEGPART*cos(LegB(1))*sin(LegB(2))*(sin(LegB(0))/2 + (sqrt(3)*cos(LegB(0)))/2) - LENGTH_OF_LEGPART*cos(LegB(2))*sin(LegB(1))*(sin(LegB(0))/2 + (sqrt(3)*cos(LegB(0)))/2)]
[                                                                    cos(LegB(1))*cos(LegB(2)) + sin(LegB(1))*sin(LegB(2)),                                                                       cos(LegB(2))*sin(LegB(1)) - cos(LegB(1))*sin(LegB(2)),                               0,                                                                                                                 - LENGTH_OF_LEGPART - LENGTH_OF_LEGPART*cos(LegB(1)) - LENGTH_OF_LEGPART*cos(LegB(1))*cos(LegB(2)) - LENGTH_OF_LEGPART*sin(LegB(1))*sin(LegB(2))]
[                                                                                                    0,                                                                                                       0,                               0,                                                                                                                                                                       1]
 
 
 
T_base12 =
 
[cos(LegC(1))*sin(LegC(2))*(cos(LegC(0))/2 + (sqrt(3)*sin(LegC(0)))/2) - cos(LegC(2))*sin(LegC(1))*(cos(LegC(0))/2 + (sqrt(3)*sin(LegC(0)))/2), cos(LegC(1))*cos(LegC(2))*(cos(LegC(0))/2 + (sqrt(3)*sin(LegC(0)))/2) + sin(LegC(1))*sin(LegC(2))*(cos(LegC(0))/2 + (sqrt(3)*sin(LegC(0)))/2), (sqrt(3)*cos(LegC(0)))/2 - sin(LegC(0))/2, LENGTH_OF_LEGPART*sin(LegC(1))*(cos(LegC(0))/2 + (sqrt(3)*sin(LegC(0)))/2) - (sqrt(3)*TRIANGLE_SIDE)/6 - LENGTH_OF_LEGPART*cos(LegC(1))*sin(LegC(2))*(cos(LegC(0))/2 + (sqrt(3)*sin(LegC(0)))/2) + LENGTH_OF_LEGPART*cos(LegC(2))*sin(LegC(1))*(cos(LegC(0))/2 + (sqrt(3)*sin(LegC(0)))/2)]
[cos(LegC(1))*sin(LegC(2))*(sin(LegC(0))/2 - (sqrt(3)*cos(LegC(0)))/2) - cos(LegC(2))*sin(LegC(1))*(sin(LegC(0))/2 - (sqrt(3)*cos(LegC(0)))/2), cos(LegC(1))*cos(LegC(2))*(sin(LegC(0))/2 - (sqrt(3)*cos(LegC(0)))/2) + sin(LegC(1))*sin(LegC(2))*(sin(LegC(0))/2 - (sqrt(3)*cos(LegC(0)))/2), cos(LegC(0))/2 + (sqrt(3)*sin(LegC(0)))/2,           TRIANGLE_SIDE/2 + LENGTH_OF_LEGPART*sin(LegC(1))*(sin(LegC(0))/2 - (sqrt(3)*cos(LegC(0)))/2) - LENGTH_OF_LEGPART*cos(LegC(1))*sin(LegC(2))*(sin(LegC(0))/2 - (sqrt(3)*cos(LegC(0)))/2) + LENGTH_OF_LEGPART*cos(LegC(2))*sin(LegC(1))*(sin(LegC(0))/2 - (sqrt(3)*cos(LegC(0)))/2)]
[                                                                    cos(LegC(1))*cos(LegC(2)) + sin(LegC(1))*sin(LegC(2)),                                                                     cos(LegC(2))*sin(LegC(1)) - cos(LegC(1))*sin(LegC(2)),                               0,                                                                                                                 - LENGTH_OF_LEGPART - LENGTH_OF_LEGPART*cos(LegC(1)) - LENGTH_OF_LEGPART*cos(LegC(1))*cos(LegC(2)) - LENGTH_OF_LEGPART*sin(LegC(1))*sin(LegC(2))]
[                                                                                                    0,                                                                                                     0,                               0,                                                                                                                                                                       1]*/
/*J_A =
 
[LENGTH_OF_LEGPART*sin(LegA(0))*sin(LegA(1)) - LENGTH_OF_LEGPART*cos(LegA(1))*sin(LegA(0))*sin(LegA(2)) + LENGTH_OF_LEGPART*cos(LegA(2))*sin(LegA(0))*sin(LegA(1)), - LENGTH_OF_LEGPART*cos(LegA(0))*cos(LegA(1)) - LENGTH_OF_LEGPART*cos(LegA(0))*cos(LegA(1))*cos(LegA(2)) - LENGTH_OF_LEGPART*cos(LegA(0))*sin(LegA(1))*sin(LegA(2)), LENGTH_OF_LEGPART*cos(LegA(0))*cos(LegA(1))*cos(LegA(2)) + LENGTH_OF_LEGPART*cos(LegA(0))*sin(LegA(1))*sin(LegA(2))]
[LENGTH_OF_LEGPART*cos(LegA(0))*cos(LegA(1))*sin(LegA(2)) - LENGTH_OF_LEGPART*cos(LegA(0))*sin(LegA(1)) - LENGTH_OF_LEGPART*cos(LegA(0))*cos(LegA(2))*sin(LegA(1)), - LENGTH_OF_LEGPART*cos(LegA(1))*sin(LegA(0)) - LENGTH_OF_LEGPART*cos(LegA(1))*cos(LegA(2))*sin(LegA(0)) - LENGTH_OF_LEGPART*sin(LegA(0))*sin(LegA(1))*sin(LegA(2)), LENGTH_OF_LEGPART*cos(LegA(1))*cos(LegA(2))*sin(LegA(0)) + LENGTH_OF_LEGPART*sin(LegA(0))*sin(LegA(1))*sin(LegA(2))]
[                                                                        0,                           LENGTH_OF_LEGPART*sin(LegA(1)) - LENGTH_OF_LEGPART*cos(LegA(1))*sin(LegA(2)) + LENGTH_OF_LEGPART*cos(LegA(2))*sin(LegA(1)),                 LENGTH_OF_LEGPART*cos(LegA(1))*sin(LegA(2)) - LENGTH_OF_LEGPART*cos(LegA(2))*sin(LegA(1))]
[                                                                        0,                                                              sin(LegA(0) + pi/2),                                        sin(LegA(0) + pi/2)]
[                                                                        0,                                                             -cos(LegA(0) + pi/2),                                       -cos(LegA(0) + pi/2)]
[                                                                        1,                                                                           0,                                                     0]
 
 
J_B =
 
[LENGTH_OF_LEGPART*sin(LegB(1))*(sin(LegB(0))/2 + (sqrt(3)*cos(LegB(0)))/2) - LENGTH_OF_LEGPART*cos(LegB(1))*sin(LegB(2))*(sin(LegB(0))/2 + (sqrt(3)*cos(LegB(0)))/2) + LENGTH_OF_LEGPART*cos(LegB(2))*sin(LegB(1))*(sin(LegB(0))/2 + (sqrt(3)*cos(LegB(0)))/2), - LENGTH_OF_LEGPART*cos(LegB(1))*(cos(LegB(0))/2 - (sqrt(3)*sin(LegB(0)))/2) - LENGTH_OF_LEGPART*cos(LegB(1))*cos(LegB(2))*(cos(LegB(0))/2 - (sqrt(3)*sin(LegB(0)))/2) - LENGTH_OF_LEGPART*sin(LegB(1))*sin(LegB(2))*(cos(LegB(0))/2 - (sqrt(3)*sin(LegB(0)))/2), LENGTH_OF_LEGPART*cos(LegB(1))*cos(LegB(2))*(cos(LegB(0))/2 - (sqrt(3)*sin(LegB(0)))/2) + LENGTH_OF_LEGPART*sin(LegB(1))*sin(LegB(2))*(cos(LegB(0))/2 - (sqrt(3)*sin(LegB(0)))/2)]
[LENGTH_OF_LEGPART*cos(LegB(1))*sin(LegB(2))*(cos(LegB(0))/2 - (sqrt(3)*sin(LegB(0)))/2) - LENGTH_OF_LEGPART*sin(LegB(1))*(cos(LegB(0))/2 - (sqrt(3)*sin(LegB(0)))/2) - LENGTH_OF_LEGPART*cos(LegB(2))*sin(LegB(1))*(cos(LegB(0))/2 - (sqrt(3)*sin(LegB(0)))/2), - LENGTH_OF_LEGPART*cos(LegB(1))*(sin(LegB(0))/2 + (sqrt(3)*cos(LegB(0)))/2) - LENGTH_OF_LEGPART*cos(LegB(1))*cos(LegB(2))*(sin(LegB(0))/2 + (sqrt(3)*cos(LegB(0)))/2) - LENGTH_OF_LEGPART*sin(LegB(1))*sin(LegB(2))*(sin(LegB(0))/2 + (sqrt(3)*cos(LegB(0)))/2), LENGTH_OF_LEGPART*cos(LegB(1))*cos(LegB(2))*(sin(LegB(0))/2 + (sqrt(3)*cos(LegB(0)))/2) + LENGTH_OF_LEGPART*sin(LegB(1))*sin(LegB(2))*(sin(LegB(0))/2 + (sqrt(3)*cos(LegB(0)))/2)]
[                                                                                                                                                      0,                                                                                                         LENGTH_OF_LEGPART*sin(LegB(1)) - LENGTH_OF_LEGPART*cos(LegB(1))*sin(LegB(2)) + LENGTH_OF_LEGPART*cos(LegB(2))*sin(LegB(1)),                                                                     LENGTH_OF_LEGPART*cos(LegB(1))*sin(LegB(2)) - LENGTH_OF_LEGPART*cos(LegB(2))*sin(LegB(1))]
[                                                                                                                                                cos(LegB(0)),                                                                                                                                                         0,                                                                                                   sin(LegB(0))]
[                                                                                                                                                sin(LegB(0)),                                                                                                                                                         0,                                                                                                  -cos(LegB(0))]
[                                                                                                                                                      1,                                                                                                                                                         0,                                                                                                         0]
 
 
J_C =
 
[LENGTH_OF_LEGPART*cos(LegC(1))*sin(LegC(2))*(sin(LegC(0))/2 - (sqrt(3)*cos(LegC(0)))/2) - LENGTH_OF_LEGPART*sin(LegC(1))*(sin(LegC(0))/2 - (sqrt(3)*cos(LegC(0)))/2) - LENGTH_OF_LEGPART*cos(LegC(2))*sin(LegC(1))*(sin(LegC(0))/2 - (sqrt(3)*cos(LegC(0)))/2), LENGTH_OF_LEGPART*cos(LegC(1))*(cos(LegC(0))/2 + (sqrt(3)*sin(LegC(0)))/2) + LENGTH_OF_LEGPART*cos(LegC(1))*cos(LegC(2))*(cos(LegC(0))/2 + (sqrt(3)*sin(LegC(0)))/2) + LENGTH_OF_LEGPART*sin(LegC(1))*sin(LegC(2))*(cos(LegC(0))/2 + (sqrt(3)*sin(LegC(0)))/2), - LENGTH_OF_LEGPART*cos(LegC(1))*cos(LegC(2))*(cos(LegC(0))/2 + (sqrt(3)*sin(LegC(0)))/2) - LENGTH_OF_LEGPART*sin(LegC(1))*sin(LegC(2))*(cos(LegC(0))/2 + (sqrt(3)*sin(LegC(0)))/2)]
[LENGTH_OF_LEGPART*sin(LegC(1))*(cos(LegC(0))/2 + (sqrt(3)*sin(LegC(0)))/2) - LENGTH_OF_LEGPART*cos(LegC(1))*sin(LegC(2))*(cos(LegC(0))/2 + (sqrt(3)*sin(LegC(0)))/2) + LENGTH_OF_LEGPART*cos(LegC(2))*sin(LegC(1))*(cos(LegC(0))/2 + (sqrt(3)*sin(LegC(0)))/2), LENGTH_OF_LEGPART*cos(LegC(1))*(sin(LegC(0))/2 - (sqrt(3)*cos(LegC(0)))/2) + LENGTH_OF_LEGPART*cos(LegC(1))*cos(LegC(2))*(sin(LegC(0))/2 - (sqrt(3)*cos(LegC(0)))/2) + LENGTH_OF_LEGPART*sin(LegC(1))*sin(LegC(2))*(sin(LegC(0))/2 - (sqrt(3)*cos(LegC(0)))/2), - LENGTH_OF_LEGPART*cos(LegC(1))*cos(LegC(2))*(sin(LegC(0))/2 - (sqrt(3)*cos(LegC(0)))/2) - LENGTH_OF_LEGPART*sin(LegC(1))*sin(LegC(2))*(sin(LegC(0))/2 - (sqrt(3)*cos(LegC(0)))/2)]
[                                                                                                                                                      0,                                                                                                       LENGTH_OF_LEGPART*sin(LegC(1)) - LENGTH_OF_LEGPART*cos(LegC(1))*sin(LegC(2)) + LENGTH_OF_LEGPART*cos(LegC(2))*sin(LegC(1)),                                                                       LENGTH_OF_LEGPART*cos(LegC(1))*sin(LegC(2)) - LENGTH_OF_LEGPART*cos(LegC(2))*sin(LegC(1))]
[                                                                                                                                                cos(LegC(0)),                                                                                                                                                       0,                                                                                                     sin(LegC(0))]
[                                                                                                                                                sin(LegC(0)),                                                                                                                                                       0,                                                                                                    -cos(LegC(0))]
[                                                                                                                                                      1,                                                                                                                                                       0,                                                                                                           0]                                                                                                                                                    1,                                                                                                                                                         0,                                                                                                         0]
*/
/*Пересчет T_base7 =
 
[cos(LegB(1))*sin(LegB(2))*(cos(LegB(0))/2 + (sqrt(3)*sin(LegB(0)))/2) - cos(LegB(2))*sin(LegB(1))*(cos(LegB(0))/2 + (sqrt(3)*sin(LegB(0)))/2), cos(LegB(1))*cos(LegB(2))*(cos(LegB(0))/2 + (sqrt(3)*sin(LegB(0)))/2) + sin(LegB(1))*sin(LegB(2))*(cos(LegB(0))/2 + (sqrt(3)*sin(LegB(0)))/2), (sqrt(3)*cos(LegB(0)))/2 - sin(LegB(0))/2, LENGTH_OF_LEGPART*sin(LegB(1))*(cos(LegB(0))/2 + (sqrt(3)*sin(LegB(0)))/2) - (sqrt(3)*TRIANGLE_SIDE)/6 - LENGTH_OF_LEGPART*cos(LegB(1))*sin(LegB(2))*(cos(LegB(0))/2 + (sqrt(3)*sin(LegB(0)))/2) + LENGTH_OF_LEGPART*cos(LegB(2))*sin(LegB(1))*(cos(LegB(0))/2 + (sqrt(3)*sin(LegB(0)))/2)]
[cos(LegB(1))*sin(LegB(2))*(sin(LegB(0))/2 - (sqrt(3)*cos(LegB(0)))/2) - cos(LegB(2))*sin(LegB(1))*(sin(LegB(0))/2 - (sqrt(3)*cos(LegB(0)))/2), cos(LegB(1))*cos(LegB(2))*(sin(LegB(0))/2 - (sqrt(3)*cos(LegB(0)))/2) + sin(LegB(1))*sin(LegB(2))*(sin(LegB(0))/2 - (sqrt(3)*cos(LegB(0)))/2), cos(LegB(0))/2 + (sqrt(3)*sin(LegB(0)))/2,           TRIANGLE_SIDE/2 + LENGTH_OF_LEGPART*sin(LegB(1))*(sin(LegB(0))/2 - (sqrt(3)*cos(LegB(0)))/2) - LENGTH_OF_LEGPART*cos(LegB(1))*sin(LegB(2))*(sin(LegB(0))/2 - (sqrt(3)*cos(LegB(0)))/2) + LENGTH_OF_LEGPART*cos(LegB(2))*sin(LegB(1))*(sin(LegB(0))/2 - (sqrt(3)*cos(LegB(0)))/2)]
[                                                                    cos(LegB(1))*cos(LegB(2)) + sin(LegB(1))*sin(LegB(2)),                                                                     cos(LegB(2))*sin(LegB(1)) - cos(LegB(1))*sin(LegB(2)),                               0,                                                                                                                 - LENGTH_OF_LEGPART - LENGTH_OF_LEGPART*cos(LegB(1)) - LENGTH_OF_LEGPART*cos(LegB(1))*cos(LegB(2)) - LENGTH_OF_LEGPART*sin(LegB(1))*sin(LegB(2))]
[                                                                                                    0,                                                                                                     0,                               0,                                                                                                                                                                       1]
 
 
Jv =
 
[LENGTH_OF_LEGPART*cos(LegB(1))*sin(LegB(2))*(sin(LegB(0))/2 - (sqrt(3)*cos(LegB(0)))/2) - LENGTH_OF_LEGPART*sin(LegB(1))*(sin(LegB(0))/2 - (sqrt(3)*cos(LegB(0)))/2) - LENGTH_OF_LEGPART*cos(LegB(2))*sin(LegB(1))*(sin(LegB(0))/2 - (sqrt(3)*cos(LegB(0)))/2), LENGTH_OF_LEGPART*cos(LegB(1))*(cos(LegB(0))/2 + (sqrt(3)*sin(LegB(0)))/2) + LENGTH_OF_LEGPART*cos(LegB(1))*cos(LegB(2))*(cos(LegB(0))/2 + (sqrt(3)*sin(LegB(0)))/2) + LENGTH_OF_LEGPART*sin(LegB(1))*sin(LegB(2))*(cos(LegB(0))/2 + (sqrt(3)*sin(LegB(0)))/2), - LENGTH_OF_LEGPART*cos(LegB(1))*cos(LegB(2))*(cos(LegB(0))/2 + (sqrt(3)*sin(LegB(0)))/2) - LENGTH_OF_LEGPART*sin(LegB(1))*sin(LegB(2))*(cos(LegB(0))/2 + (sqrt(3)*sin(LegB(0)))/2)]
[LENGTH_OF_LEGPART*sin(LegB(1))*(cos(LegB(0))/2 + (sqrt(3)*sin(LegB(0)))/2) - LENGTH_OF_LEGPART*cos(LegB(1))*sin(LegB(2))*(cos(LegB(0))/2 + (sqrt(3)*sin(LegB(0)))/2) + LENGTH_OF_LEGPART*cos(LegB(2))*sin(LegB(1))*(cos(LegB(0))/2 + (sqrt(3)*sin(LegB(0)))/2), LENGTH_OF_LEGPART*cos(LegB(1))*(sin(LegB(0))/2 - (sqrt(3)*cos(LegB(0)))/2) + LENGTH_OF_LEGPART*cos(LegB(1))*cos(LegB(2))*(sin(LegB(0))/2 - (sqrt(3)*cos(LegB(0)))/2) + LENGTH_OF_LEGPART*sin(LegB(1))*sin(LegB(2))*(sin(LegB(0))/2 - (sqrt(3)*cos(LegB(0)))/2), - LENGTH_OF_LEGPART*cos(LegB(1))*cos(LegB(2))*(sin(LegB(0))/2 - (sqrt(3)*cos(LegB(0)))/2) - LENGTH_OF_LEGPART*sin(LegB(1))*sin(LegB(2))*(sin(LegB(0))/2 - (sqrt(3)*cos(LegB(0)))/2)]
[                                                                                                                                                      0,                                                                                                       LENGTH_OF_LEGPART*sin(LegB(1)) - LENGTH_OF_LEGPART*cos(LegB(1))*sin(LegB(2)) + LENGTH_OF_LEGPART*cos(LegB(2))*sin(LegB(1)),                                                                       LENGTH_OF_LEGPART*cos(LegB(1))*sin(LegB(2)) - LENGTH_OF_LEGPART*cos(LegB(2))*sin(LegB(1))]
 
 
T_base12 =
 
[cos(LegC(1))*sin(LegC(2))*(cos(LegC(0))/2 - (sqrt(3)*sin(LegC(0)))/2) - cos(LegC(2))*sin(LegC(1))*(cos(LegC(0))/2 - (sqrt(3)*sin(LegC(0)))/2), cos(LegC(1))*cos(LegC(2))*(cos(LegC(0))/2 - (sqrt(3)*sin(LegC(0)))/2) + sin(LegC(1))*sin(LegC(2))*(cos(LegC(0))/2 - (sqrt(3)*sin(LegC(0)))/2), - sin(LegC(0))/2 - (sqrt(3)*cos(LegC(0)))/2, LENGTH_OF_LEGPART*sin(LegC(1))*(cos(LegC(0))/2 - (sqrt(3)*sin(LegC(0)))/2) - (sqrt(3)*TRIANGLE_SIDE)/6 - LENGTH_OF_LEGPART*cos(LegC(1))*sin(LegC(2))*(cos(LegC(0))/2 - (sqrt(3)*sin(LegC(0)))/2) + LENGTH_OF_LEGPART*cos(LegC(2))*sin(LegC(1))*(cos(LegC(0))/2 - (sqrt(3)*sin(LegC(0)))/2)]
[cos(LegC(1))*sin(LegC(2))*(sin(LegC(0))/2 + (sqrt(3)*cos(LegC(0)))/2) - cos(LegC(2))*sin(LegC(1))*(sin(LegC(0))/2 + (sqrt(3)*cos(LegC(0)))/2), cos(LegC(1))*cos(LegC(2))*(sin(LegC(0))/2 + (sqrt(3)*cos(LegC(0)))/2) + sin(LegC(1))*sin(LegC(2))*(sin(LegC(0))/2 + (sqrt(3)*cos(LegC(0)))/2),   cos(LegC(0))/2 - (sqrt(3)*sin(LegC(0)))/2,           LENGTH_OF_LEGPART*sin(LegC(1))*(sin(LegC(0))/2 + (sqrt(3)*cos(LegC(0)))/2) - TRIANGLE_SIDE/2 - LENGTH_OF_LEGPART*cos(LegC(1))*sin(LegC(2))*(sin(LegC(0))/2 + (sqrt(3)*cos(LegC(0)))/2) + LENGTH_OF_LEGPART*cos(LegC(2))*sin(LegC(1))*(sin(LegC(0))/2 + (sqrt(3)*cos(LegC(0)))/2)]
[                                                                    cos(LegC(1))*cos(LegC(2)) + sin(LegC(1))*sin(LegC(2)),                                                                     cos(LegC(2))*sin(LegC(1)) - cos(LegC(1))*sin(LegC(2)),                                 0,                                                                                                                 - LENGTH_OF_LEGPART - LENGTH_OF_LEGPART*cos(LegC(1)) - LENGTH_OF_LEGPART*cos(LegC(1))*cos(LegC(2)) - LENGTH_OF_LEGPART*sin(LegC(1))*sin(LegC(2))]
[                                                                                                    0,                                                                                                     0,                                 0,                                                                                                                                                                       1]
 
 
Jv =
 
[LENGTH_OF_LEGPART*cos(LegC(1))*sin(LegC(2))*(sin(LegC(0))/2 + (sqrt(3)*cos(LegC(0)))/2) - LENGTH_OF_LEGPART*sin(LegC(1))*(sin(LegC(0))/2 + (sqrt(3)*cos(LegC(0)))/2) - LENGTH_OF_LEGPART*cos(LegC(2))*sin(LegC(1))*(sin(LegC(0))/2 + (sqrt(3)*cos(LegC(0)))/2), LENGTH_OF_LEGPART*cos(LegC(1))*(cos(LegC(0))/2 - (sqrt(3)*sin(LegC(0)))/2) + LENGTH_OF_LEGPART*cos(LegC(1))*cos(LegC(2))*(cos(LegC(0))/2 - (sqrt(3)*sin(LegC(0)))/2) + LENGTH_OF_LEGPART*sin(LegC(1))*sin(LegC(2))*(cos(LegC(0))/2 - (sqrt(3)*sin(LegC(0)))/2), - LENGTH_OF_LEGPART*cos(LegC(1))*cos(LegC(2))*(cos(LegC(0))/2 - (sqrt(3)*sin(LegC(0)))/2) - LENGTH_OF_LEGPART*sin(LegC(1))*sin(LegC(2))*(cos(LegC(0))/2 - (sqrt(3)*sin(LegC(0)))/2)]
[LENGTH_OF_LEGPART*sin(LegC(1))*(cos(LegC(0))/2 - (sqrt(3)*sin(LegC(0)))/2) - LENGTH_OF_LEGPART*cos(LegC(1))*sin(LegC(2))*(cos(LegC(0))/2 - (sqrt(3)*sin(LegC(0)))/2) + LENGTH_OF_LEGPART*cos(LegC(2))*sin(LegC(1))*(cos(LegC(0))/2 - (sqrt(3)*sin(LegC(0)))/2), LENGTH_OF_LEGPART*cos(LegC(1))*(sin(LegC(0))/2 + (sqrt(3)*cos(LegC(0)))/2) + LENGTH_OF_LEGPART*cos(LegC(1))*cos(LegC(2))*(sin(LegC(0))/2 + (sqrt(3)*cos(LegC(0)))/2) + LENGTH_OF_LEGPART*sin(LegC(1))*sin(LegC(2))*(sin(LegC(0))/2 + (sqrt(3)*cos(LegC(0)))/2), - LENGTH_OF_LEGPART*cos(LegC(1))*cos(LegC(2))*(sin(LegC(0))/2 + (sqrt(3)*cos(LegC(0)))/2) - LENGTH_OF_LEGPART*sin(LegC(1))*sin(LegC(2))*(sin(LegC(0))/2 + (sqrt(3)*cos(LegC(0)))/2)]
[                                                                                                                                                      0,                                                                                                       LENGTH_OF_LEGPART*sin(LegC(1)) - LENGTH_OF_LEGPART*cos(LegC(1))*sin(LegC(2)) + LENGTH_OF_LEGPART*cos(LegC(2))*sin(LegC(1)),                                                                       LENGTH_OF_LEGPART*cos(LegC(1))*sin(LegC(2)) - LENGTH_OF_LEGPART*cos(LegC(2))*sin(LegC(1))]*/
MatrixXd ForwardKinematicsRelatedToRobotCOM(Vector3d LegA, Vector3d LegB, Vector3d LegC) {
    MatrixXd res(3,3);//coordinates of end effectors in system related with robot
    res(0,0) = LENGTH_OF_LEGPART * cos(LegA(0))*cos(LegA(1))*sin(LegA(2)) - LENGTH_OF_LEGPART * cos(LegA(0))*sin(LegA(1)) - (sqrt(3))/3 * TRIANGLE_SIDE - LENGTH_OF_LEGPART * cos(LegA(0)) * cos(LegA(2)) * sin(LegA(1));
    res(1,0) = LENGTH_OF_LEGPART * cos(LegA(1))*sin(LegA(0))*sin(LegA(2)) - LENGTH_OF_LEGPART * sin(LegA(0))*sin(LegA(1)) - LENGTH_OF_LEGPART*cos(LegA(2))*sin(LegA(0))*sin(LegA(1));
    res(2,0) = -LENGTH_OF_LEGPART - LENGTH_OF_LEGPART*cos(LegA(1)) - LENGTH_OF_LEGPART*cos(LegA(1))*cos(LegA(2)) - LENGTH_OF_LEGPART*sin(LegA(1))*sin(LegA(2));
    res(0,1) = LENGTH_OF_LEGPART*sin(LegB(1))*(cos(LegB(0))/2 + (sqrt(3)*sin(LegB(0)))/2) - (sqrt(3)*TRIANGLE_SIDE)/6 - LENGTH_OF_LEGPART*cos(LegB(1))*sin(LegB(2))*(cos(LegB(0))/2 + (sqrt(3)*sin(LegB(0)))/2) + LENGTH_OF_LEGPART*cos(LegB(2))*sin(LegB(1))*(cos(LegB(0))/2 + (sqrt(3)*sin(LegB(0)))/2);
    res(1,1) = TRIANGLE_SIDE/2 + LENGTH_OF_LEGPART*sin(LegB(1))*(sin(LegB(0))/2 - (sqrt(3)*cos(LegB(0)))/2) - LENGTH_OF_LEGPART*cos(LegB(1))*sin(LegB(2))*(sin(LegB(0))/2 - (sqrt(3)*cos(LegB(0)))/2) + LENGTH_OF_LEGPART*cos(LegB(2))*sin(LegB(1))*(sin(LegB(0))/2 - (sqrt(3)*cos(LegB(0)))/2);
    res(2,1) = -LENGTH_OF_LEGPART - LENGTH_OF_LEGPART*cos(LegB(1)) - LENGTH_OF_LEGPART*cos(LegB(1))*cos(LegB(2)) - LENGTH_OF_LEGPART*sin(LegB(1))*sin(LegB(2));
    res(0,2) = LENGTH_OF_LEGPART*sin(LegC(1))*(cos(LegC(0))/2 - (sqrt(3)*sin(LegC(0)))/2) - (sqrt(3)*TRIANGLE_SIDE)/6 - LENGTH_OF_LEGPART*cos(LegC(1))*sin(LegC(2))*(cos(LegC(0))/2 - (sqrt(3)*sin(LegC(0)))/2) + LENGTH_OF_LEGPART*cos(LegC(2))*sin(LegC(1))*(cos(LegC(0))/2 - (sqrt(3)*sin(LegC(0)))/2);
    res(1,2) = LENGTH_OF_LEGPART*sin(LegC(1))*(sin(LegC(0))/2 + (sqrt(3)*cos(LegC(0)))/2) - TRIANGLE_SIDE/2 - LENGTH_OF_LEGPART*cos(LegC(1))*sin(LegC(2))*(sin(LegC(0))/2 + (sqrt(3)*cos(LegC(0)))/2) + LENGTH_OF_LEGPART*cos(LegC(2))*sin(LegC(1))*(sin(LegC(0))/2 + (sqrt(3)*cos(LegC(0)))/2);
    res(2,2) = - LENGTH_OF_LEGPART - LENGTH_OF_LEGPART*cos(LegC(1)) - LENGTH_OF_LEGPART*cos(LegC(1))*cos(LegC(2)) - LENGTH_OF_LEGPART*sin(LegC(1))*sin(LegC(2));
    return res;
}
Vector3d NumericalInverseKinematicsOnTheNextStepLegA(Vector3d LegALinearSpeed, Vector3d LegA, double step) {
    MatrixXd Jv(3,3);
    Jv(0,0) = LENGTH_OF_LEGPART*sin(LegA(0))*sin(LegA(1)) - LENGTH_OF_LEGPART*cos(LegA(1))*sin(LegA(0))*sin(LegA(2)) + LENGTH_OF_LEGPART*cos(LegA(2))*sin(LegA(0))*sin(LegA(1));
    Jv(0,1) = - LENGTH_OF_LEGPART*cos(LegA(0))*cos(LegA(1)) - LENGTH_OF_LEGPART*cos(LegA(0))*cos(LegA(1))*cos(LegA(2)) - LENGTH_OF_LEGPART*cos(LegA(0))*sin(LegA(1))*sin(LegA(2));
    Jv(0,2) = LENGTH_OF_LEGPART*cos(LegA(0))*cos(LegA(1))*cos(LegA(2)) + LENGTH_OF_LEGPART*cos(LegA(0))*sin(LegA(1))*sin(LegA(2));
    Jv(1,0) = LENGTH_OF_LEGPART*cos(LegA(0))*cos(LegA(1))*sin(LegA(2)) - LENGTH_OF_LEGPART*cos(LegA(0))*sin(LegA(1)) - LENGTH_OF_LEGPART*cos(LegA(0))*cos(LegA(2))*sin(LegA(1));
    Jv(1,1) = - LENGTH_OF_LEGPART*cos(LegA(1))*sin(LegA(0)) - LENGTH_OF_LEGPART*cos(LegA(1))*cos(LegA(2))*sin(LegA(0)) - LENGTH_OF_LEGPART*sin(LegA(0))*sin(LegA(1))*sin(LegA(2));
    Jv(1,2) = LENGTH_OF_LEGPART*cos(LegA(1))*cos(LegA(2))*sin(LegA(0)) + LENGTH_OF_LEGPART*sin(LegA(0))*sin(LegA(1))*sin(LegA(2));
    Jv(2,0) = 0;
    Jv(2,1) = LENGTH_OF_LEGPART*sin(LegA(1)) - LENGTH_OF_LEGPART*cos(LegA(1))*sin(LegA(2)) + LENGTH_OF_LEGPART*cos(LegA(2))*sin(LegA(1));
    Jv(2,2) = LENGTH_OF_LEGPART*cos(LegA(1))*sin(LegA(2)) - LENGTH_OF_LEGPART*cos(LegA(2))*sin(LegA(1));
    Vector3d LegANext = LegA;
    if (Jv.determinant() != 0) {
        Matrix3d Jv_inv = Jv.inverse();
        LegANext+=LegALinearSpeed.transpose()*Jv_inv*step;
    }
    return LegANext;
}
Vector3d NumericalInverseKinematicsOnTheNextStepLegB(Vector3d LegBLinearSpeed, Vector3d LegB, double step) {
    MatrixXd Jv(3,3);
    Jv(0,0) = LENGTH_OF_LEGPART*cos(LegB(1))*sin(LegB(2))*(sin(LegB(0))/2 - (sqrt(3)*cos(LegB(0)))/2) - LENGTH_OF_LEGPART*sin(LegB(1))*(sin(LegB(0))/2 - (sqrt(3)*cos(LegB(0)))/2) - LENGTH_OF_LEGPART*cos(LegB(2))*sin(LegB(1))*(sin(LegB(0))/2 - (sqrt(3)*cos(LegB(0)))/2);
    Jv(0,1) = LENGTH_OF_LEGPART*cos(LegB(1))*(cos(LegB(0))/2 + (sqrt(3)*sin(LegB(0)))/2) + LENGTH_OF_LEGPART*cos(LegB(1))*cos(LegB(2))*(cos(LegB(0))/2 + (sqrt(3)*sin(LegB(0)))/2) + LENGTH_OF_LEGPART*sin(LegB(1))*sin(LegB(2))*(cos(LegB(0))/2 + (sqrt(3)*sin(LegB(0)))/2);
    Jv(0,2) = - LENGTH_OF_LEGPART*cos(LegB(1))*cos(LegB(2))*(cos(LegB(0))/2 + (sqrt(3)*sin(LegB(0)))/2) - LENGTH_OF_LEGPART*sin(LegB(1))*sin(LegB(2))*(cos(LegB(0))/2 + (sqrt(3)*sin(LegB(0)))/2);
    Jv(1,0) = LENGTH_OF_LEGPART*sin(LegB(1))*(cos(LegB(0))/2 + (sqrt(3)*sin(LegB(0)))/2) - LENGTH_OF_LEGPART*cos(LegB(1))*sin(LegB(2))*(cos(LegB(0))/2 + (sqrt(3)*sin(LegB(0)))/2) + LENGTH_OF_LEGPART*cos(LegB(2))*sin(LegB(1))*(cos(LegB(0))/2 + (sqrt(3)*sin(LegB(0)))/2);
    Jv(1,1) = LENGTH_OF_LEGPART*cos(LegB(1))*(sin(LegB(0))/2 - (sqrt(3)*cos(LegB(0)))/2) + LENGTH_OF_LEGPART*cos(LegB(1))*cos(LegB(2))*(sin(LegB(0))/2 - (sqrt(3)*cos(LegB(0)))/2) + LENGTH_OF_LEGPART*sin(LegB(1))*sin(LegB(2))*(sin(LegB(0))/2 - (sqrt(3)*cos(LegB(0)))/2);
    Jv(1,2) = - LENGTH_OF_LEGPART*cos(LegB(1))*cos(LegB(2))*(sin(LegB(0))/2 - (sqrt(3)*cos(LegB(0)))/2) - LENGTH_OF_LEGPART*sin(LegB(1))*sin(LegB(2))*(sin(LegB(0))/2 - (sqrt(3)*cos(LegB(0)))/2);
    Jv(2,0) = 0;
    Jv(2,1) = LENGTH_OF_LEGPART*sin(LegB(1)) - LENGTH_OF_LEGPART*cos(LegB(1))*sin(LegB(2)) + LENGTH_OF_LEGPART*cos(LegB(2))*sin(LegB(1));
    Jv(2,2) = LENGTH_OF_LEGPART*cos(LegB(1))*sin(LegB(2)) - LENGTH_OF_LEGPART*cos(LegB(2))*sin(LegB(1));
    Vector3d LegBNext = LegB;
    if (Jv.determinant() != 0) {
        Matrix3d Jv_inv = Jv.inverse();
        LegBNext+=LegBLinearSpeed.transpose()*Jv_inv*step;
    }
    return LegBNext;
}
Vector3d NumericalInverseKinematicsOnTheNextStepLegC(Vector3d LegCLinearSpeed, Vector3d LegC, double step) {
    MatrixXd Jv(3,3);
    Jv(0,0) = LENGTH_OF_LEGPART*cos(LegC(1))*sin(LegC(2))*(sin(LegC(0))/2 + (sqrt(3)*cos(LegC(0)))/2) - LENGTH_OF_LEGPART*sin(LegC(1))*(sin(LegC(0))/2 + (sqrt(3)*cos(LegC(0)))/2) - LENGTH_OF_LEGPART*cos(LegC(2))*sin(LegC(1))*(sin(LegC(0))/2 + (sqrt(3)*cos(LegC(0)))/2);
    Jv(0,1) = LENGTH_OF_LEGPART*cos(LegC(1))*(cos(LegC(0))/2 - (sqrt(3)*sin(LegC(0)))/2) + LENGTH_OF_LEGPART*cos(LegC(1))*cos(LegC(2))*(cos(LegC(0))/2 - (sqrt(3)*sin(LegC(0)))/2) + LENGTH_OF_LEGPART*sin(LegC(1))*sin(LegC(2))*(cos(LegC(0))/2 - (sqrt(3)*sin(LegC(0)))/2);
    Jv(0,2) = - LENGTH_OF_LEGPART*cos(LegC(1))*cos(LegC(2))*(cos(LegC(0))/2 - (sqrt(3)*sin(LegC(0)))/2) - LENGTH_OF_LEGPART*sin(LegC(1))*sin(LegC(2))*(cos(LegC(0))/2 - (sqrt(3)*sin(LegC(0)))/2);
    Jv(1,0) = LENGTH_OF_LEGPART*sin(LegC(1))*(cos(LegC(0))/2 - (sqrt(3)*sin(LegC(0)))/2) - LENGTH_OF_LEGPART*cos(LegC(1))*sin(LegC(2))*(cos(LegC(0))/2 - (sqrt(3)*sin(LegC(0)))/2) + LENGTH_OF_LEGPART*cos(LegC(2))*sin(LegC(1))*(cos(LegC(0))/2 - (sqrt(3)*sin(LegC(0)))/2);
    Jv(1,1) = LENGTH_OF_LEGPART*cos(LegC(1))*(sin(LegC(0))/2 + (sqrt(3)*cos(LegC(0)))/2) + LENGTH_OF_LEGPART*cos(LegC(1))*cos(LegC(2))*(sin(LegC(0))/2 + (sqrt(3)*cos(LegC(0)))/2) + LENGTH_OF_LEGPART*sin(LegC(1))*sin(LegC(2))*(sin(LegC(0))/2 + (sqrt(3)*cos(LegC(0)))/2);
    Jv(1,2) = - LENGTH_OF_LEGPART*cos(LegC(1))*cos(LegC(2))*(sin(LegC(0))/2 + (sqrt(3)*cos(LegC(0)))/2) - LENGTH_OF_LEGPART*sin(LegC(1))*sin(LegC(2))*(sin(LegC(0))/2 + (sqrt(3)*cos(LegC(0)))/2);
    Jv(2,0) = 0;
    Jv(2,1) = LENGTH_OF_LEGPART*sin(LegC(1)) - LENGTH_OF_LEGPART*cos(LegC(1))*sin(LegC(2)) + LENGTH_OF_LEGPART*cos(LegC(2))*sin(LegC(1));
    Jv(2,2) = LENGTH_OF_LEGPART*cos(LegC(1))*sin(LegC(2)) - LENGTH_OF_LEGPART*cos(LegC(2))*sin(LegC(1));
    Vector3d LegCNext = LegC;
    if (Jv.determinant() != 0) {
        Matrix3d Jv_inv = Jv.inverse();
        LegCNext+=LegCLinearSpeed.transpose()*Jv_inv*step;
    }
    return LegCNext;
}
int main()
{
    /*using namespace std;
    for (int q1 = -90; q1 < 91; q1++) {
        for (int q2 = -90; q2 < 91; q2++) {
            for(int q3 = -90; q3 < 91; q3++) {
                for (int LegB(0) = -90; LegB(0) < 91; LegB(0)++) {
                    for (int LegB(1) = -90; LegB(1) < 91; LegB(1)++) {
                        for(int LegB(2) = -90; LegB(2) < 91; LegB(2)++) {
                            for (int LegC(0) = -90; LegC(0) < 91; LegC(0)++) {
                                for (int LegC(1) = -90; LegC(1) < 91; LegC(1)++) {
                                    for(int LegC(2) = -90; LegC(2) < 91; LegC(2)++) {
                                        std::cout << ForwardKinematicsRelatedToRobotCOM({q1*M_PI/180, q2*M_PI/180, q3*M_PI/180}, {LegB(0)*M_PI/180, LegB(1)*M_PI/180, LegB(2)*M_PI/180}, {LegC(0)*M_PI/180, LegC(1)*M_PI/180, LegC(2)*M_PI/180}) << std::endl << std::endl;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }*/

    return 0;
}