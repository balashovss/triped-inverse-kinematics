#include <gnuplot-iostream.h>
#include <iostream>
#include <vector>
#include <random>
using namespace std;
using FUNC_PTR = double(*)(vector<double>, vector<double>, double);
using FUNC_PTR_FOR_SQR = char(*)(vector<double>);
#define NUMBER_OF_ITER 300// количество итераций
#define ITER_STEP 2.0/(NUMBER_OF_ITER-1) // шаг вектора
#define NUMBER_OF_WOLVES 120 // количество волков
#define NUMBER_OF_COEFS 3 // число координат
#define SAVING_NUMBER_OF_WOLVES 3
#define TIME_STEP 0.1
#define PREDICTION_TIME MAX*TIME_STEP
#define EPS 0.005
double get_random_num(double, double);
double get_random_num_in_range_zero_to_one();
double func_for_testing(vector<double> a, vector<double> q, double tf);
vector<vector<double>> create_matrix(unsigned int height, unsigned int width, vector<double> max, vector<double> min);
void matrix_sort(vector<vector<double>>* matrix, unsigned int height, vector<double> a, FUNC_PTR func, double tf);
double control_functional(vector<double> x_tf, vector<double> x_f, double tf);
double vector_abs(vector<double> vect);
vector<double> vector_div(vector<double> vect1, vector<double> vect2, unsigned int width);
vector<double> vector_sum(vector<double> vect1, vector<double> vect2, unsigned int width);
vector<double> vector_prod(vector<double> vect1, vector<double> vect2, unsigned int width);
vector<double> vect_prod_num(vector<double> vect, double num, unsigned int width);
vector<double> get_random_vector(unsigned int width);
void print_matrix(vector<vector<double>> matrix, unsigned int height, unsigned int width, vector<double> func_coeffs, FUNC_PTR func, double tf);
void visualisation(vector<double> data);
vector<vector<double>> create_coord_matrix(unsigned int height, unsigned int width);
void visualisation_of_robot_coords(double time);
struct robot_state {
    double x[NUMBER_OF_COEFS+1];
    double speed[NUMBER_OF_COEFS+1];
    double control[NUMBER_OF_COEFS];
    robot_state(double x[NUMBER_OF_COEFS+1], double speed[NUMBER_OF_COEFS+1], double control[NUMBER_OF_COEFS]) {
        this->x[0] = x[0];
        this->x[1] = x[1];
        this->x[2] = x[2];
        this->speed[0] = speed[0];
        this->speed[1] = speed[1];
        this->speed[2] = speed[2];
        this->control[0] = control[0];
        this->control[1] = control[1];
    }
    robot_state() {
        this->x[0] = 0;
        this->x[1] = 0;
        this->x[2] = 0;
        this->speed[0] = 0;
        this->speed[1] = 0;
        this->speed[2] = 0;
        this->control[0] = 0;
        this->control[1] = 0;
    }
    void moveRobot(double time_step){
        this->speed[2] = 0.5 * (this->control[0] - this->control[1]);
        this->x[2] += this->speed[2]*time_step;
        this->speed[0] = 0.5 * (this->control[0] + this->control[1]) * cos(this->x[2]);
        this->speed[1] = 0.5 * (this->control[0] + this->control[1]) * sin(this->x[2]);
        this->x[0] += this->speed[0]*time_step;
        this->x[1] += this->speed[1]*time_step;
    }
};
void robots_sort(robot_state* Robots, int height, vector<double> final_point, double tf, FUNC_PTR func);
void robots_sort(robot_state* Robots, int height, vector<double> final_point, double tf, double start_time, double final_time);
char sqrt_heviside(vector<double> vec);
double integrate(robot_state Robot, double time_step, double start_time, double final_time, FUNC_PTR_FOR_SQR func);
double control_functional(vector<double> x_tf, vector<double> x_f, double tf, robot_state Robot, double start_time, double final_time);
