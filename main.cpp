#include <iostream>
#include <eigen3/Eigen/Dense>
#include <cmath>
using Eigen::Matrix4d;
using Eigen::MatrixXd;
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
 
[cos(LegC(2))*sin(LegC(1))*(cos(LegC(0))/2 - (sqrt(3)*sin(LegC(0)))/2) - cos(LegC(1))*sin(LegC(2))*(cos(LegC(0))/2 - (sqrt(3)*sin(LegC(0)))/2), - cos(LegC(1))*cos(LegC(2))*(cos(LegC(0))/2 - (sqrt(3)*sin(LegC(0)))/2) - sin(LegC(1))*sin(LegC(2))*(cos(LegC(0))/2 - (sqrt(3)*sin(LegC(0)))/2), sin(LegC(0))/2 + (sqrt(3)*cos(LegC(0)))/2, (sqrt(3)*TRIANGLE_SIDE)/6 - LENGTH_OF_LEGPART*sin(LegC(1))*(cos(LegC(0))/2 - (sqrt(3)*sin(LegC(0)))/2) + LENGTH_OF_LEGPART*cos(LegC(1))*sin(LegC(2))*(cos(LegC(0))/2 - (sqrt(3)*sin(LegC(0)))/2) - LENGTH_OF_LEGPART*cos(LegC(2))*sin(LegC(1))*(cos(LegC(0))/2 - (sqrt(3)*sin(LegC(0)))/2)]
[cos(LegC(2))*sin(LegC(1))*(sin(LegC(0))/2 + (sqrt(3)*cos(LegC(0)))/2) - cos(LegC(1))*sin(LegC(2))*(sin(LegC(0))/2 + (sqrt(3)*cos(LegC(0)))/2), - cos(LegC(1))*cos(LegC(2))*(sin(LegC(0))/2 + (sqrt(3)*cos(LegC(0)))/2) - sin(LegC(1))*sin(LegC(2))*(sin(LegC(0))/2 + (sqrt(3)*cos(LegC(0)))/2), (sqrt(3)*sin(LegC(0)))/2 - cos(LegC(0))/2,           TRIANGLE_SIDE/2 - LENGTH_OF_LEGPART*sin(LegC(1))*(sin(LegC(0))/2 + (sqrt(3)*cos(LegC(0)))/2) + LENGTH_OF_LEGPART*cos(LegC(1))*sin(LegC(2))*(sin(LegC(0))/2 + (sqrt(3)*cos(LegC(0)))/2) - LENGTH_OF_LEGPART*cos(LegC(2))*sin(LegC(1))*(sin(LegC(0))/2 + (sqrt(3)*cos(LegC(0)))/2)]
[                                                                    cos(LegC(1))*cos(LegC(2)) + sin(LegC(1))*sin(LegC(2)),                                                                       cos(LegC(2))*sin(LegC(1)) - cos(LegC(1))*sin(LegC(2)),                               0,                                                                                                                 - LENGTH_OF_LEGPART - LENGTH_OF_LEGPART*cos(LegC(1)) - LENGTH_OF_LEGPART*cos(LegC(1))*cos(LegC(2)) - LENGTH_OF_LEGPART*sin(LegC(1))*sin(LegC(2))]
[                                                                                                    0,                                                                                                       0,                               0,                                                                                                                                                                       1]
*/
MatrixXd ForwardKinematicsRelatedToRobotCOM(Vector3d LegA, Vector3d LegB, Vector3d LegC) {
    MatrixXd res(3,3);//coordinates of end effector
    res(0,0) = LENGTH_OF_LEGPART * cos(LegA(0))*cos(LegA(1))*sin(LegA(2)) - LENGTH_OF_LEGPART * cos(LegA(0))*sin(LegA(1)) - (sqrt(3))/3 * TRIANGLE_SIDE - LENGTH_OF_LEGPART * cos(LegA(0)) * cos(LegA(2)) * sin(LegA(1));
    res(1,0) = LENGTH_OF_LEGPART * cos(LegA(1))*sin(LegA(0))*sin(LegA(2)) - LENGTH_OF_LEGPART * sin(LegA(0))*sin(LegA(1)) - LENGTH_OF_LEGPART*cos(LegA(2))*sin(LegA(0))*sin(LegA(1));
    res(2,0) = -LENGTH_OF_LEGPART - LENGTH_OF_LEGPART*cos(LegA(1)) - LENGTH_OF_LEGPART*cos(LegA(1))*cos(LegA(2)) - LENGTH_OF_LEGPART*sin(LegA(1))*sin(LegA(2));
    res(0,1) = (sqrt(3)*TRIANGLE_SIDE)/6 - LENGTH_OF_LEGPART*sin(LegB(1))*(cos(LegB(0))/2 - (sqrt(3)*sin(LegB(0)))/2) + LENGTH_OF_LEGPART*cos(LegB(1))*sin(LegB(2))*(cos(LegB(0))/2 - (sqrt(3)*sin(LegB(0)))/2) - LENGTH_OF_LEGPART*cos(LegB(2))*sin(LegB(1))*(cos(LegB(0))/2 - (sqrt(3)*sin(LegB(0)))/2);
    res(1,1) = TRIANGLE_SIDE/2 - LENGTH_OF_LEGPART*sin(LegB(1))*(sin(LegB(0))/2 + (sqrt(3)*cos(LegB(0)))/2) + LENGTH_OF_LEGPART*cos(LegB(1))*sin(LegB(2))*(sin(LegB(0))/2 + (sqrt(3)*cos(LegB(0)))/2) - LENGTH_OF_LEGPART*cos(LegB(2))*sin(LegB(1))*(sin(LegB(0))/2 + (sqrt(3)*cos(LegB(0)))/2);
    res(2,1) = -LENGTH_OF_LEGPART - LENGTH_OF_LEGPART*cos(LegB(1)) - LENGTH_OF_LEGPART*cos(LegB(1))*cos(LegB(2)) - LENGTH_OF_LEGPART*sin(LegB(1))*sin(LegB(2));
    res(0,2) = (sqrt(3)*TRIANGLE_SIDE)/6 - LENGTH_OF_LEGPART*sin(LegC(1))*(cos(LegC(0))/2 - (sqrt(3)*sin(LegC(0)))/2) + LENGTH_OF_LEGPART*cos(LegC(1))*sin(LegC(2))*(cos(LegC(0))/2 - (sqrt(3)*sin(LegC(0)))/2) - LENGTH_OF_LEGPART*cos(LegC(2))*sin(LegC(1))*(cos(LegC(0))/2 - (sqrt(3)*sin(LegC(0)))/2);
    res(1,2) = TRIANGLE_SIDE/2 - LENGTH_OF_LEGPART*sin(LegC(1))*(sin(LegC(0))/2 + (sqrt(3)*cos(LegC(0)))/2) + LENGTH_OF_LEGPART*cos(LegC(1))*sin(LegC(2))*(sin(LegC(0))/2 + (sqrt(3)*cos(LegC(0)))/2) - LENGTH_OF_LEGPART*cos(LegC(2))*sin(LegC(1))*(sin(LegC(0))/2 + (sqrt(3)*cos(LegC(0)))/2);
    res(2,2) = -LENGTH_OF_LEGPART - LENGTH_OF_LEGPART*cos(LegC(1)) - LENGTH_OF_LEGPART*cos(LegC(1))*cos(LegC(2)) - LENGTH_OF_LEGPART*sin(LegC(1))*sin(LegC(2));
    return res;
}

int main()
{
    using namespace std;
    for (int q1 = -90; q1 < 91; q1++) {
        for (int q2 = -90; q2 < 91; q2++) {
            for(int q3 = -90; q3 < 91; q3++) {
                for (int q4 = -90; q4 < 91; q4++) {
                    for (int q5 = -90; q5 < 91; q5++) {
                        for(int q6 = -90; q6 < 91; q6++) {
                            for (int q7 = -90; q7 < 91; q7++) {
                                for (int q8 = -90; q8 < 91; q8++) {
                                    for(int q9 = -90; q9 < 91; q9++) {
                                        std::cout << ForwardKinematicsRelatedToRobotCOM({q1*M_PI/2, q2*M_PI/2, q3*M_PI/2}, {q4*M_PI/2, q5*M_PI/2, q6*M_PI/2}, {q7*M_PI/2, q8*M_PI/2, q9*M_PI/2}) << std::endl << std::endl;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    return 0;
}