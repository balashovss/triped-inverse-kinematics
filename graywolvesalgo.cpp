#include <gnuplot-iostream.h>
#include <iostream>
#include <vector>
#include <random>
#include "graywolvesalgo.hpp"
//множество допустимых решений [MIN; MAX] // количество не двигающихся в текущей итерации волков(почему-то без этого не работает)
double get_random_num(double max, double min) {return (double)(rand())/RAND_MAX*(max-min) + min;}
double get_random_num_in_range_zero_to_one() {return (double)(rand())/RAND_MAX;}
double func_for_testing(vector<double> a, vector<double> q,double tf) {
    double func = 0;
    for (unsigned int i = 0; i < NUMBER_OF_COEFS; i++) {
        func+= (double)((a[i]-q[i])*(a[i]-q[i]));
    }
    return func+tf;
}
double control_functional(vector<double> x_tf, vector<double> x_f, double tf) {
    double func = 0;
    for (unsigned int i = 0; i < NUMBER_OF_COEFS + 1; i++) {
        func+= (double)((x_tf[i]-x_f[i])*(x_tf[i]-x_f[i]));
    }
    double first_sqrt = 0;
    double second_sqrt = 0;
    for (unsigned int i = 0; i < NUMBER_OF_COEFS; i++) {
        first_sqrt += (x_tf[i] - 2.5) * (x_tf[i] - 2.5);
        second_sqrt += (x_tf[i] - 7.5) * (x_tf[i] - 7.5);
    }
    return sqrt(func) + (sqrt(first_sqrt) < 2.5) * 10000 + (sqrt(second_sqrt) < 2.5) * 10000;
}
double control_functional(vector<double> x_tf, vector<double> x_f, double tf, robot_state Robot, double start_time, double final_time) {
    double func = 0;
    for (unsigned int i = 0; i < NUMBER_OF_COEFS + 1; i++) {
        func+= (double)((x_tf[i]-x_f[i])*(x_tf[i]-x_f[i]));
    }
    return sqrt(func) + integrate(Robot, TIME_STEP, start_time, final_time, &sqrt_heviside);
}
double integrate(robot_state Robot, double time_step, double start_time, double final_time, FUNC_PTR_FOR_SQR func) {
    robot_state TempRobot = Robot;
    double x1_start = Robot.x[0];
    double x2_start = Robot.x[1];
    double val = 0;
    double temp_time = start_time;
    while (temp_time <= final_time) {
        TempRobot.moveRobot(time_step);
        val += func({TempRobot.x[0], TempRobot.x[1]});
        temp_time += time_step;
    }
    return val;
}
vector<vector<double>> create_matrix(unsigned int height, unsigned int width, vector<double> max, vector<double> min) {
    vector<vector<double>> matrix_of_wolves(height);
    for (unsigned int i = 0; i < height; i++) {
        for (unsigned int j = 0; j < width; j++) {
            matrix_of_wolves[i].push_back(get_random_num(max[j], min[j]));
        }
    };
    return matrix_of_wolves;
}
vector<vector<double>> create_coord_matrix(unsigned int height, unsigned int width){
    vector<vector<double>> matrix_of_wolves(height);
    for (unsigned int i = 0; i < height; i++) {
        for (unsigned int j = 0; j < width - 1; j++) {
            matrix_of_wolves[i].push_back(10);
        }
        matrix_of_wolves[i][width - 1] = 0;
    };
    return matrix_of_wolves;
}
void matrix_sort(vector<vector<double>>* matrix, unsigned int height, vector<double> a, FUNC_PTR func, double tf) {
    for (unsigned int i = 0; i < height - 1; i++) {
        for (unsigned int j = 0; j < height - 1 - i; j++) {
            double temp1 = (*func)(a, (*matrix)[j], tf);
            double temp2 = (*func)(a, (*matrix)[j+1], tf);
            if (temp1 > temp2) {
                (*matrix)[j].swap((*matrix)[j+1]);
            }
        }
    }
}
double vector_abs(vector<double> vect) {
    double sum = 0;
    for (unsigned int i = 0; i < 10; i++) {
        sum += vect[i]*vect[i];
    }
    return sqrt(sum);
}
vector<double> vector_div(vector<double> vect1, vector<double> vect2, unsigned int width) {
    vector<double> res_vect(width);
    for (unsigned int i = 0; i < width; i++) {
        res_vect[i] = vect1[i] - vect2[i];
    }
    return res_vect;
}
vector<double> vector_sum(vector<double> vect1, vector<double> vect2, unsigned int width) {
    vector<double> res_vect(width);
    for (unsigned int i = 0; i < width; i++) {
        res_vect[i] = vect1[i] + vect2[i];
    }
    return res_vect;
}
vector<double> vector_prod(vector<double> vect1, vector<double> vect2, unsigned int width) {
    vector<double> res_vect(width);
    for (unsigned int i = 0; i < width; i++) {
        res_vect[i] = vect1[i] * vect2[i];
    }
    return res_vect;
}
vector<double> vect_prod_num(vector<double> vect, double num, unsigned int width) {
    vector<double> res_vect(width);
    for (unsigned int i = 0; i < width; i++) {
        res_vect[i] = vect[i]*num;
    }
    return res_vect;
}
vector<double> get_random_vector(unsigned int width) {
    vector<double> r1(width);
    for (unsigned int i = 0; i < width; i++) {
        r1[i] = get_random_num_in_range_zero_to_one();
    }
    return r1;
}
void print_matrix(vector<vector<double>> matrix, unsigned int height, unsigned int width, vector<double> func_coeffs, FUNC_PTR func,double tf) {
    for (unsigned int i = 0; i < height; i++) {
        printf("Wolf%u:", i+1);
        for (unsigned int j = 0; j < width; j++) {
            printf("%lf ", matrix[i][j]);
        }
        printf("F = %lf\n", (*func)(func_coeffs, matrix[i], tf));
    }
}
void visualisation(vector<double> data) {
    Gnuplot gp("gnuplot -persist");
    gp << "set xlabel 'Количество итераций'\n";
    gp << "set ylabel 'Значения функций для альфа волка'\n";
    gp << "set xrange [1:500]\n";
    gp << "plot '-' with lines title 'Alpha Wolf Path'\n";
    gp.send1d(data);
    //gp << "plot 0\n";
}
void visualisation_of_robot_coords(double time) {
    Gnuplot gp("gnuplot -persist");
    gp << "set xlabel 'x1'\n";
    gp << "set ylabel 'x2'\n";
    gp << "set autoscale\n";
    gp << "plot 'data.dat' with lines title 'Robot Path', 'walls.dat' with lines, 'walls2.dat' with lines\n";
    //gp << "plot 'walls.dat' with lines\n";
    //gp << "plot 'walls2.dat' with lines\n";
    printf("%.3lf\n", time);
    //gp << "plot 0\n";
}
void robots_sort(robot_state* Robots, int height, vector<double> final_point, double tf, FUNC_PTR func) {
    for (unsigned int i = 0; i < height - 1; i++) {
        for (unsigned int j = 0; j < height - 1 - i; j++) {
            double temp1 = (*func)(final_point, {Robots[j].x[0], Robots[j].x[1], Robots[j].x[2]}, tf);
            double temp2 = (*func)(final_point, {Robots[j+1].x[0], Robots[j+1].x[1], Robots[j+1].x[2]}, tf);
            if (temp1 > temp2) {
                robot_state TempRobot = Robots[j];
                Robots[j] = Robots[j+1];
                Robots[j+1] = TempRobot;
            }
        }
    }
}
void robots_sort(robot_state* Robots, int height, vector<double> final_point, double tf, double start_time, double final_time) {
    for (unsigned int i = 0; i < height - 1; i++) {
        for (unsigned int j = 0; j < height - 1 - i; j++) {
            double temp1 = control_functional(final_point, {Robots[j].x[0], Robots[j].x[1], Robots[j].x[2]}, tf, Robots[j], start_time, final_time);
            double temp2 = control_functional(final_point, {Robots[j+1].x[0], Robots[j+1].x[1], Robots[j+1].x[2]}, tf, Robots[j+1], start_time, final_time);
            if (temp1 > temp2) {
                robot_state TempRobot = Robots[j];
                Robots[j] = Robots[j+1];
                Robots[j+1] = TempRobot;
            }
        }
    }
}
char sqrt_heviside(vector<double> vec) {
    return (sqrt((vec[0] - 2.5) * (vec[0] - 2.5) + (vec[1] - 2.5) * (vec[1] - 2.5)) < 2.5 + sqrt((vec[0] - 7.5) * (vec[0] - 7.5) + (vec[1] - 7.5) * (vec[1] - 7.5)) < 2.5);
}
