#include <math.h>
#include <stdio.h>
#include <time.h>

#define NUM_PROCS 4 //2 4 8
#define NUM_INTERVALS 100 // 25 50 75 125 150 175
#define TAU 0.2 // Maximum allowed error ratio between new_trap_1 + new_trap_2 and old_trap

static double total_time_ver_1;
static double time_ver_1_per_proc[NUM_PROCS];
static double total_time_ver_2;


double inv_function(double x) {
    return 1.0/x;
}

void print_times() {
    printf("\nNUM_PROCS = %d", NUM_PROCS);
    printf("Total time for version 1: %f\n", total_time_ver_1);
    printf("Time taken to finish:\n");
    for (int proc_index = 0; proc_index < NUM_PROCS; ++proc_index) {
        printf("\tproc #%d: %f\n", proc_index, time_ver_1_per_proc[proc_index]);
    }

    printf("\nNUM_INTERVALS = %d", NUM_INTERVALS);
    printf("Total time for version 2: %f\n", total_time_ver_2);
}

inline double trap_area (double left, double f_left, double right, double f_right) {
    return 0.5 * (right - left) * (f_left + f_right);
}

inline double error_ratio (double approx_left_half, double approx_right_half, double prev_approx) {
    return fabs( ((approx_left_half + approx_right_half) / prev_approx) - 1.0 );
}

double adap_quadrature(double (*function_ptr)(double),
                       double left,
                       double f_left,
                       double right,
                       double f_right,
                       double prev_approx) {
    double center = 0.5 * (left + right);
    double f_center = (*function_ptr)(center);
    double approx_left_half = trap_area(left, f_left, center, f_center);
    double approx_right_half = trap_area(center, f_center, right, f_right);
    if (error_ratio(approx_left_half, approx_right_half, prev_approx) > TAU) {
        double left_area = adap_quadrature(function_ptr,
                                           left,
                                           f_left,
                                           center,
                                           f_center,
                                           approx_left_half);
        double right_area = adap_quadrature(function_ptr,
                                            center,
                                            f_center,
                                            right,
                                            f_right,
                                            approx_right_half);
        return left_area + right_area;
    } else {
        return approx_left_half + approx_right_half;
    }
}

double version_1(double (*function_ptr)(double), double left, double right) {
    // Unimplemented!
    double f_left = (*function_ptr)(left);
    double f_right = (*function_ptr)(right);
    double orig_approx = trap_area(left,
                                   f_left,
                                   right,
                                   f_right);
    return orig_approx;
}

double version_2(double (*function_ptr)(double), double left, double right) {
    // Unimplemented!
    double f_left = (*function_ptr)(left);
    double f_right = (*function_ptr)(right);
    double orig_approx = trap_area(left,
                                   f_left,
                                   right,
                                   f_right);
    return orig_approx;
}

int main() {
    // Define here the function we want to integrate
    double (*function_ptr)(double) = &inv_function;
    double left = 0.1;
    double right = 10.0;

    double version_1_area = version_1(function_ptr, left, right);
    double version_2_area = version_2(function_ptr, left, right);

    printf("Area from version 1: %f\n", version_1_area);
    printf("Area from version 2: %f\n", version_2_area);
    //print_times();
}
