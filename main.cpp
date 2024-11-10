#include <iostream>
#include <eigen3/Eigen/Dense>
#include <cmath>
#include <gnuplot-iostream.h>
#include <vector>
#include <thread>
#include <chrono>
#include <mutex>
#include "graywolvesalgo.hpp"
using Eigen::Matrix4d;
using Eigen::MatrixXd;
using Eigen::Matrix3d;
using Eigen::Vector3d;
using Eigen::VectorXd;   
using std::vector;
#define LENGTH_OF_LEGPART 0.21
#define TRIANGLE_SIDE 0.4
#define TIME_ON_STEP 1
#define TIME_STEP_IN_MAIN TIME_ON_STEP/(2*NUMBER_OF_ITERATIONS)
#define NUMBER_OF_ITERATIONS 1000
#define HEIGHT_OF_STEP 0.28
#define LENGTH_OF_STEP 0.065
#define NUMBER_OF_COORDS 3
#define NUMBER_OF_THREADS 10
vector<double> wolves_algo(vector<double> func_coeffs, FUNC_PTR func, int number_of_coeffs) {
    vector<double> max(3);
    srand(time(NULL));
    max[0] = 0;
    max[1] = 3*M_PI/4;
    max[2] = 0;
    vector<double> min(3);
    min[0] = 0;
    min[1] = 0;
    min[2] = -M_PI/2;
    vector<vector<double>> matrix_of_wolves = create_matrix(NUMBER_OF_WOLVES, number_of_coeffs, max, min);
    matrix_sort(&matrix_of_wolves, NUMBER_OF_WOLVES, func_coeffs, func, 0);
    for (int t = 0; t < NUMBER_OF_ITER; t++) {
        vector<vector<double>> three_best_wolves{matrix_of_wolves[0], matrix_of_wolves[1], matrix_of_wolves[2]};
        vector<double> a (number_of_coeffs, (double)(2.0*(2.0-(double)t*ITER_STEP)));
        for (unsigned int i = SAVING_NUMBER_OF_WOLVES; i < NUMBER_OF_WOLVES;i++) {
            vector<vector<double>> A = {vector_div(vector_prod(a,get_random_vector(number_of_coeffs),number_of_coeffs),a, number_of_coeffs), vector_div(vector_prod(a,get_random_vector(number_of_coeffs),number_of_coeffs),a, number_of_coeffs), vector_div(vector_prod(a,get_random_vector(number_of_coeffs),number_of_coeffs),a, number_of_coeffs), vector_div(vector_prod(a,get_random_vector(number_of_coeffs),number_of_coeffs),a, number_of_coeffs)};
            vector<vector<double>> C = {vect_prod_num(get_random_vector(number_of_coeffs), 2, number_of_coeffs), vect_prod_num(get_random_vector(number_of_coeffs), 2, number_of_coeffs), vect_prod_num(get_random_vector(number_of_coeffs), 2, number_of_coeffs)};
            vector<vector<double>> D = {vector_div(vector_prod(C[0], three_best_wolves[0], number_of_coeffs),matrix_of_wolves[i], number_of_coeffs), vector_div(vector_prod(C[1], three_best_wolves[1], number_of_coeffs),matrix_of_wolves[i], number_of_coeffs), vector_div(vector_prod(C[2], three_best_wolves[2], number_of_coeffs),matrix_of_wolves[i], number_of_coeffs)};
            vector<vector<double>> X_next = {vector_div(three_best_wolves[0], vector_prod(A[0], D[0], number_of_coeffs), number_of_coeffs), vector_div(three_best_wolves[1], vector_prod(A[1], D[1], number_of_coeffs), number_of_coeffs), vector_div(three_best_wolves[2], vector_prod(A[2], D[2], number_of_coeffs), number_of_coeffs)};
            matrix_of_wolves[i] = vect_prod_num(vector_sum(vector_sum(X_next[0], X_next[1], number_of_coeffs), X_next[2], number_of_coeffs), 1.0/3.0, number_of_coeffs);
            for (int j = 0; j < number_of_coeffs; j++) {
                if (matrix_of_wolves[i][j] > max[j]) matrix_of_wolves[i][j] = max[j];
                if (matrix_of_wolves[i][j] < min[j]) matrix_of_wolves[i][j] = min[j];
            }
        }
        matrix_sort(&matrix_of_wolves, NUMBER_OF_WOLVES, func_coeffs, func, 0);
    }
    vector<double> val(0);
    for (int i = 0; i < number_of_coeffs; i++) {
        val.push_back(matrix_of_wolves[0][i]);
    }
    printf(" error : %lf coords %lf %lf %lf ", (*func)(func_coeffs, val, 0), matrix_of_wolves[0][0], matrix_of_wolves[0][1], matrix_of_wolves[0][2]);
    return val;
}
MatrixXd ForwardKinematicsMatrixOfRobotRelatedToFrameOnEarthLegA(Vector3d LegA){//world frame is moving with robot
    MatrixXd val(4,4);
    val(0,0) = cos(LegA[0])*sin(LegA[1]+LegA[2]);
    val(0,1) = cos(LegA[0])*cos(LegA[1]+LegA[2]);
    val(0,2) = sin(LegA[0]);
    val(0,3) = TRIANGLE_SIDE/sqrt(3) + LENGTH_OF_LEGPART*cos(LegA[0])*(sin(LegA[1])+sin(LegA[1]+LegA[2]));
    val(1,0) = sin(LegA[0])*sin(LegA[1]+LegA[2]);
    val(1,1) = sin(LegA[0])*cos(LegA[1]+LegA[2]);
    val(1,2) = -cos(LegA[0]);
    val(1,3) = LENGTH_OF_LEGPART*sin(LegA[0])*(sin(LegA[1])+sin(LegA[1]+LegA[2]));
    val(2,0) = -cos(LegA[1]+LegA[2]);
    val(2,1) = sin(LegA[1]+LegA[2]);
    val(2,2) = 0;
    val(2,3) = -LENGTH_OF_LEGPART*(-2+cos(LegA[1]+LegA[2])+cos(LegA[1]));
    val(3,0) = 0;
    val(3,1) = 0;
    val(3,2) = 0;
    val(3,3) = 1;
    return val;
}
MatrixXd ForwardKinematicsMatrixOfRobotRelatedToFrameOnEarthLegB(Vector3d LegB){//world frame is moving with robot
    MatrixXd val(4,4);
    val(0,0) = -cos(LegB[0])*sin(LegB[1]+LegB[2])/2+sqrt(3)*sin(LegB[0]);
    val(0,1) = -cos(LegB[0])*cos(LegB[1]+LegB[2])/2+sqrt(3)*sin(LegB[0]);
    val(0,2) = -sin(LegB[0])/2 + sqrt(3)*cos(LegB[0])/2;
    val(0,3) = -TRIANGLE_SIDE*sqrt(3)/6 + sqrt(3)*sin(LegB[0]) - LENGTH_OF_LEGPART*cos(LegB[0])*(sin(LegB[1])+sin(LegB[1]+LegB[2]))/2;
    val(1,0) = -sin(LegB[0])*sin(LegB[1]+LegB[2])/2 - sqrt(3)*cos(LegB[0]);
    val(1,1) = -sin(LegB[0])*cos(LegB[1]+LegB[2])/2 - sqrt(3)*cos(LegB[0]);
    val(1,2) = cos(LegB[0])/2 + sqrt(3) * sin(LegB[0])/2;
    val(1,3) = TRIANGLE_SIDE/2 - sqrt(3)*cos(LegB[0]) - LENGTH_OF_LEGPART*sin(LegB[0])*(sin(LegB[1])+sin(LegB[1]+LegB[2]))/2;
    val(2,0) = -cos(LegB[1]+LegB[2]);
    val(2,1) = sin(LegB[1]+LegB[2]);
    val(2,2) = 0;
    val(2,3) = -LENGTH_OF_LEGPART*(-2+cos(LegB[1]+LegB[2])+cos(LegB[1]));
    val(3,0) = 0;
    val(3,1) = 0;
    val(3,2) = 0;
    val(3,3) = 1;
    return val;
}
MatrixXd ForwardKinematicsMatrixOfRobotRelatedToFrameOnEarthLegC(Vector3d LegC){//world frame is moving with robot
    MatrixXd val(4,4);
    val(0,0) = -cos(LegC[0])*sin(LegC[1]+LegC[2])/2 - sqrt(3)*sin(LegC[0]);
    val(0,1) = -cos(LegC[0])*cos(LegC[1]+LegC[2])/2 - sqrt(3)*sin(LegC[0]);
    val(0,2) = -sin(LegC[0])/2 - sqrt(3)*cos(LegC[0])/2;
    val(0,3) = -TRIANGLE_SIDE*sqrt(3)/6 - sqrt(3)*sin(LegC[0]) + LENGTH_OF_LEGPART*cos(LegC[0])*(-sin(LegC[1])+sin(LegC[1]+LegC[2]))/2;
    val(1,0) = -sin(LegC[0])*sin(LegC[1]+LegC[2])/2 + sqrt(3)*cos(LegC[0]);
    val(1,1) = -sin(LegC[0])*sin(LegC[1]+LegC[2])/2 + sqrt(3)*cos(LegC[0]);
    val(1,2) = cos(LegC[0])/2 - sqrt(3) * sin(LegC[0])/2;
    val(1,3) = -TRIANGLE_SIDE/2 + sqrt(3)*cos(LegC[0]) + LENGTH_OF_LEGPART*sin(LegC[0])*(sin(LegC[1] - sin(LegC[1]+LegC[2])))/2;
    val(2,0) = -cos(LegC[1]+LegC[2]);
    val(2,1) = sin(LegC[1]+LegC[2]);
    val(2,2) = 0;
    val(2,3) = -LENGTH_OF_LEGPART*(-2+cos(LegC[1]+LegC[2])+cos(LegC[1]));
    val(3,0) = 0;
    val(3,1) = 0;
    val(3,2) = 0;
    val(3,3) = 1;
    return val;
}
Vector3d NumericalInverseKinematicsOnTheNextStepLegA(VectorXd LegASpeed, Vector3d LegA, double step) {//it doesn't work(GG)
    MatrixXd J(6,3);
    J(0,0) = -LENGTH_OF_LEGPART * sin(LegA(0)) * (sin(LegA(1)) + sin(LegA(1)+LegA(2)));
    J(0,1) = LENGTH_OF_LEGPART * cos(LegA(0)) * (cos(LegA(1)) + cos(LegA(1) + LegA(2)));
    J(0,2) = LENGTH_OF_LEGPART * cos(LegA(0)) * cos(LegA(1) + LegA(2));
    J(1,0) = LENGTH_OF_LEGPART * cos(LegA(0)) * (sin(LegA(1)) + sin(LegA(1)+LegA(2)));
    J(1,1) = LENGTH_OF_LEGPART * sin(LegA(0)) * (cos(LegA(1)) + cos(LegA(1) + LegA(2)));
    J(1,2) = LENGTH_OF_LEGPART * sin(LegA(0)) * cos(LegA(1)+LegA(2));
    J(2,0) = 0;
    J(2,1) = LENGTH_OF_LEGPART * (sin(LegA(1)) + sin(LegA(1)+LegA(2)));
    J(2,2) = LENGTH_OF_LEGPART * sin(LegA(1) + LegA(2));
    J(3,0) = -sin(LegA(0)) * sin(LegA(1)+LegA(2));
    J(3,1) = cos(LegA(0)) * cos(LegA(1)+LegA(2));
    J(3,2) = cos(LegA(0)) * cos(LegA(1)+LegA(2));
    J(4,0) = cos(LegA(0)) * sin(LegA(1) + LegA(2));
    J(4,1) = sin(LegA(0)) * cos(LegA(1)+LegA(2));
    J(4,2) = sin(LegA(0)) * cos(LegA(1)+LegA(2));
    J(5,0) = 0;
    J(5,1) = sin(LegA(1) + LegA(2));
    J(5,2) = sin(LegA(1) + LegA(2));
    Vector3d LegANext = LegA;
    MatrixXd Jv = (J.transpose()*J);
    std::cout << J << std::endl << std::endl;
    std::cout << Jv << std::endl;
    if (Jv.determinant() != 0) {
        MatrixXd Jv_inv = Jv.inverse();
        LegANext+=LegASpeed.transpose()*Jv_inv*step;
    }
    return LegANext;
}
vector<double> getCoords(vector<double> LegA) {
    vector<double> val(NUMBER_OF_COORDS);
    Vector3d tmp;
    for (int i = 0; i < NUMBER_OF_COORDS; i++) {
        tmp(i) = LegA[i];
    }
    for (int i = 0; i < NUMBER_OF_COORDS; i++) {
        val[i] = ForwardKinematicsMatrixOfRobotRelatedToFrameOnEarthLegA(tmp)(i, 3);
    }
    return val;
}
double getCoordsDiff(vector<double> coords_etalon, vector<double> LegA, double tf) {
    double val = 0;
    for (int i = 0; i < NUMBER_OF_COORDS; i++) {
        val += (getCoords(LegA)[i] - coords_etalon[i])*(getCoords(LegA)[i] - coords_etalon[i]);
    }
    return sqrt(val);
}
vector<double> InverseKinematicsOnTheNextStepLegA(vector<double> needed_coords) {
    return wolves_algo(needed_coords, &getCoordsDiff, NUMBER_OF_COORDS);
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
    Jv(0,2) = -LENGTH_OF_LEGPART*cos(LegC(1))*cos(LegC(2))*(cos(LegC(0))/2 - (sqrt(3)*sin(LegC(0)))/2) - LENGTH_OF_LEGPART*sin(LegC(1))*sin(LegC(2))*(cos(LegC(0))/2 - (sqrt(3)*sin(LegC(0)))/2);
    Jv(1,0) = LENGTH_OF_LEGPART*sin(LegC(1))*(cos(LegC(0))/2 - (sqrt(3)*sin(LegC(0)))/2) - LENGTH_OF_LEGPART*cos(LegC(1))*sin(LegC(2))*(cos(LegC(0))/2 - (sqrt(3)*sin(LegC(0)))/2) + LENGTH_OF_LEGPART*cos(LegC(2))*sin(LegC(1))*(cos(LegC(0))/2 - (sqrt(3)*sin(LegC(0)))/2);
    Jv(1,1) = LENGTH_OF_LEGPART*cos(LegC(1))*(sin(LegC(0))/2 + (sqrt(3)*cos(LegC(0)))/2) + LENGTH_OF_LEGPART*cos(LegC(1))*cos(LegC(2))*(sin(LegC(0))/2 + (sqrt(3)*cos(LegC(0)))/2) + LENGTH_OF_LEGPART*sin(LegC(1))*sin(LegC(2))*(sin(LegC(0))/2 + (sqrt(3)*cos(LegC(0)))/2);
    Jv(1,2) = -LENGTH_OF_LEGPART*cos(LegC(1))*cos(LegC(2))*(sin(LegC(0))/2 + (sqrt(3)*cos(LegC(0)))/2) - LENGTH_OF_LEGPART*sin(LegC(1))*sin(LegC(2))*(sin(LegC(0))/2 + (sqrt(3)*cos(LegC(0)))/2);
    Jv(2,0) = 0;
    Jv(2,1) = LENGTH_OF_LEGPART*sin(LegC(1)) - LENGTH_OF_LEGPART*cos(LegC(1))*sin(LegC(2)) + LENGTH_OF_LEGPART*cos(LegC(2))*sin(LegC(1));
    Jv(2,2) = LENGTH_OF_LEGPART*cos(LegC(1))*sin(LegC(2)) - LENGTH_OF_LEGPART*cos(LegC(2))*sin(LegC(1));
    Vector3d LegCNext = LegC;
    LegCNext+=Jv*LegCLinearSpeed*step;
    return LegCNext;
}
void visualisation_for_q1(vector<double> data) {
    Gnuplot gp("gnuplot -persist");
    gp << "set xlabel 'Время'\n";
    gp << "set ylabel 'q1'\n";
    gp << "set xrange [0:2]\n";//NUMBER_OF_ITERATIONS
    gp << "plot '-' with lines title 'Alpha Wolf Path'\n";
    gp.send1d(data);
    //gp << "plot 0\n";
}
void visualisation_for_q2(vector<double> data) {
    Gnuplot gp("gnuplot -persist");
    gp << "set xlabel 'Время'\n";
    gp << "set ylabel 'q2'\n";
    gp << "set xrange [0:2]\n";//NUMBER_OF_ITERATIONS
    gp << "plot '-' with lines title 'Alpha Wolf Path'\n";
    gp.send1d(data);
    //gp << "plot 0\n";
}
void visualisation_for_q3(vector<double> data) {
    Gnuplot gp("gnuplot -persist");
    gp << "set xlabel 'Время'\n";
    gp << "set ylabel 'q1'\n";
    gp << "set xrange [0:2]\n";//NUMBER_OF_ITERATIONS
    gp << "plot '-' with lines title 'Alpha Wolf Path'\n";
    gp.send1d(data);
    //gp << "plot 0\n";
}
class Leg {
    private:
    double x, y, z, q1, q2, q3;
    double speedx, speedy, speedz, speednx, speedny, speednz;
    public:
    Leg(double x, double y, double z, double speedx, double speedy, double speedz, double speednx, double speedny, double speednz, double q1, double q2, double q3) {
        this->x = x;
        this->y = y;
        this->z = z;
        this->speedx = speedx;
        this->speedy = speedy;
        this->speedz = speedz;
        this->speednx = speednx;
        this->speedny = speedny;
        this->speednz = speednz;
        this->q1 = q1;
        this->q2 = q2;
        this->q3 = q3;
    }
    void setCoords(double x, double y, double z) {
        this->x = x;
        this->y = y;
        this->z = z;
    }
    void setSpeed(double speedx, double speedy, double speedz, double speednx, double speedny, double speednz) {
        this->speedx = speedx;
        this->speedy = speedy;
        this->speedz = speedz;
        this->speednx = speednx;
        this->speedny = speedny;
        this->speednz = speednz;
    }
    void setQ(double q1, double q2, double q3) {
        this->q1 = q1;
        this->q2 = q2;
        this->q3 = q3;
    }
    void setQ(Vector3d vec) {
        this->q1 = vec[0];
        this->q2 = vec[1];
        this->q3 = vec[2];
    }
    Vector3d getQ() {
        Vector3d result;
        result[0] = this->q1;
        result[1] = this->q2;
        result[2] = this->q3;
        return result;
    }
    Vector3d getSpeed() {
        Vector3d result;
        result[0] = this->speedx;
        result[1] = this->speedy;
        result[2] = this->speedz;
        return result;  
    }
};
std::mutex mtx;
void calculateLegCoords(int start, int end, vector<double>& x, vector<double>& y, vector<double>& z, vector<double>& q1, vector<double>& q2, vector<double>& q3, FILE* fptr) {
    for (int i = start; i < end; i++) {
        vector<double> coords = InverseKinematicsOnTheNextStepLegA({x[i], y[i], z[i]});
        q1[i] = coords[0];
        q2[i] = coords[1];
        q3[i] = coords[2];
        {
            std::lock_guard<std::mutex> lock(mtx);
            fprintf(fptr, "%d %lf %lf %lf\n", i, q1[i], q2[i], q3[i]);
        }
        
        Vector3d tmp(3);
        tmp(0) = q1[i];
        tmp(1) = q2[i];
        tmp(2) = q3[i];
        printf("xyz: %lf %lf %lf\n", ForwardKinematicsMatrixOfRobotRelatedToFrameOnEarthLegA(tmp)(0,3), ForwardKinematicsMatrixOfRobotRelatedToFrameOnEarthLegA(tmp)(1,3), ForwardKinematicsMatrixOfRobotRelatedToFrameOnEarthLegA(tmp)(2,3));
    }
}
int main()
{
    vector<double> x(NUMBER_OF_ITERATIONS);
    vector<double> y(NUMBER_OF_ITERATIONS);
    vector<double> z(NUMBER_OF_ITERATIONS);
    vector<double> q1(NUMBER_OF_ITERATIONS);
    vector<double> q2(NUMBER_OF_ITERATIONS);
    vector<double> q3(NUMBER_OF_ITERATIONS);
    vector<std::thread> threads(0);
    FILE* fptr = fopen("coords.dat", "w+");
    q1[0] = 0;
    q2[0] = 0;
    q3[0] = 0;
    for (int i = 0; i < NUMBER_OF_ITERATIONS; i++) {
        x[i] = TRIANGLE_SIDE/sqrt(3) + LENGTH_OF_STEP * (1 + cos(M_PI*i/(2*NUMBER_OF_ITERATIONS*TIME_ON_STEP)));
        y[i] = 0;
        z[i] = HEIGHT_OF_STEP*sin(M_PI*i/(2*NUMBER_OF_ITERATIONS*TIME_ON_STEP));
    }
    for (int i = 0; i < NUMBER_OF_THREADS - 1; i++) {
        int start = i * NUMBER_OF_ITERATIONS/NUMBER_OF_THREADS + 1;
        int end = start + NUMBER_OF_ITERATIONS/NUMBER_OF_THREADS;
        threads.emplace_back(calculateLegCoords, start, end, std::ref(x), std::ref(y), std::ref(z), std::ref(q1), std::ref(q2), std::ref(q3), fptr);
    }
    for (auto& t : threads) {
        t.join();
    }
    visualisation_for_q1(q1);
    visualisation_for_q2(q2);
    visualisation_for_q3(q3);
    fclose(fptr);
    // for (int i = 0; i < NUMBER_OF_ITERATIONS; i++) std::cout << q2[i] << std::endl;
    // std::cout << std::endl << std::endl;
    // for (int i = 0; i < NUMBER_OF_ITERATIONS; i++) std::cout << q3[i] << std::endl;
    return 0;
}
//Переписать уравнение движения ноги. Реализовать в main алгоритм с помощью функции InverseKinematicsForLeg