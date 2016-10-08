#include <math.h>
#include <stdio.h>
#include <time.h>

#include <mpi.h>


#define NUM_WORKER_PROCS 4 //2 4 8
#define NUM_INTERVALS 100 // 25 50 75 125 150 175
#define TAU 0.2 // Maximum allowed error ratio between new_trap_1 + new_trap_2 and old_trap

static double total_time_ver_1;
static double time_ver_1_per_proc[NUM_WORKER_PROCS];
static double total_time_ver_2;


double inv_function(double x) {
    return 1.0/x;
}

double const_function(double x) {
    return 1.0;
}

void print_times() {
    printf("\nNUM_WORKER_PROCS = %d", NUM_WORKER_PROCS);
    printf("Total time for version 1: %f\n", total_time_ver_1);
    printf("Time taken to finish:\n");
    for (int proc_index = 0; proc_index < NUM_WORKER_PROCS; ++proc_index) {
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

void version_1(double (*function_ptr)(double), double left, double right,
               double* returned_total_area) {
    // Initialize the MPI environment
    MPI_Init(NULL, NULL);

    // Get the rank of the process
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    double total_area;

    if (rank == 0) {
        // Master process
        clock_t start_time = clock();

        double local_area = 0.0;
        double interval_step_size = (left + right) / (double)NUM_WORKER_PROCS;

        // Send a task to each process
        for (int proc_number = 0; proc_number < NUM_WORKER_PROCS; ++proc_number) {
            double send_data[2];
            send_data[0] = left + interval_step_size * proc_number;
            send_data[1] = left + interval_step_size * (proc_number + 1);
            MPI_Send((void*) &send_data,
                     2,
                     MPI_DOUBLE,
                     proc_number,
                     0,
                     MPI_COMM_WORLD);
        }

        MPI_Reduce((void*) &local_area,
           (void*) &total_area,
           1,
           MPI_DOUBLE,
           MPI_SUM,
           0,
           MPI_COMM_WORLD);

        *returned_total_area = total_area;

        clock_t finish_time = clock();
        total_time_ver_1 = (finish_time - start_time) / CLOCKS_PER_SEC;

    } else {
        // Worker process
        clock_t proc_start_time = clock();

        double received_data[2];
        MPI_Recv((void*) &received_data,
                 2,
                 MPI_DOUBLE,
                 0,
                 MPI_ANY_TAG,
                 MPI_COMM_WORLD,
                 MPI_STATUS_IGNORE);
        double curr_left = received_data[0];
        double curr_right = received_data[1];

        // DEBUG:
        printf("Process %d received interval (%f, %f) from process 0\n",
               rank,
               curr_left,
               curr_right);

        // Initial area calculation
        double f_curr_left = (*function_ptr)(curr_left);
        double f_curr_right = (*function_ptr)(curr_right);
        double orig_approx = trap_area(curr_left,
                                       f_curr_left,
                                       curr_right,
                                       f_curr_right);
        double local_area = adap_quadrature(function_ptr,
                                            curr_left,
                                            f_curr_left,
                                            curr_right,
                                            f_curr_right,
                                            orig_approx);

        clock_t proc_finish_time = clock();
        time_ver_1_per_proc[rank - 1] = (proc_finish_time - proc_start_time) / CLOCKS_PER_SEC;

        MPI_Reduce((void*) &local_area,
                   (void*) &total_area,
                   1,
                   MPI_DOUBLE,
                   MPI_SUM,
                   0,
                   MPI_COMM_WORLD);
    }

    MPI_Finalize();
}

void version_2(double (*function_ptr)(double), double left, double right,
               double* returned_total_area) {
    clock_t start_time = clock();
    // Initialize the MPI environment
    MPI_Init(NULL, NULL);

    // Get the rank of the process
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank == 0) {
        // Master process
        clock_t start_time = clock();

        double total_area = 0.0;
        double interval_step_size = (left + right) / (double)NUM_INTERVALS;

        // At the beginning each worker process will send 0.0 to the master process.
        int intervals_done = -NUM_WORKER_PROCS;
        int tasks_sent = 0;
        int flag, source;
        MPI_Status status;
        double received_data;
        double send_data[3]; // {bool_continue, left, right}

        while (intervals_done < NUM_INTERVALS) {
            MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &flag, &status);
            if (flag) {
                source = status.MPI_SOURCE;
                MPI_Recv((void*) &received_data,
                         1,
                         MPI_DOUBLE,
                         source,
                         0,
                         MPI_COMM_WORLD,
                         &status);
                ++intervals_done;
                total_area += received_data;

                // DEBUG
                printf("Worker %d sent %f\n", source, received_data);


                if (tasks_sent < NUM_INTERVALS) {
                    send_data[0] = 1.0;
                    send_data[1] = left + interval_step_size * tasks_sent;
                    send_data[2] = left + interval_step_size * (tasks_sent + 1);
                    MPI_Send((void*) &send_data,
                             3,
                             MPI_DOUBLE,
                             source,
                             1,
                             MPI_COMM_WORLD);
                    ++tasks_sent;
                } else {
                    // End worker processes
                    send_data[0] = 0.0;
                    MPI_Send((void*) &send_data,
                             3,
                             MPI_DOUBLE,
                             source,
                             1,
                             MPI_COMM_WORLD);
                }
            }
        }
        *returned_total_area = total_area;

        clock_t finish_time = clock();
        total_time_ver_2 = (finish_time - start_time) / CLOCKS_PER_SEC;

    } else {
        // Worker process
        clock_t proc_start_time = clock();

        double local_area = 0.0;
        double received_data[3] = {1.0, 0.0, 0.0};
        double curr_left, curr_right;

        while(received_data[0]) {
            curr_left = received_data[0];
            curr_right = received_data[1];

            if (curr_left != curr_right) {
                // We added the `if` above just to make sure we don't try to calculate the function
                // at the origin in case it is not well-defined.

                // DEBUG:
                printf("Process %d received interval (%f, %f) from process 0\n",
                       rank,
                       curr_left,
                       curr_right);

                // Initial area calculation
                double f_curr_left = (*function_ptr)(curr_left);
                double f_curr_right = (*function_ptr)(curr_right);
                double orig_approx = trap_area(curr_left,
                                               f_curr_left,
                                               curr_right,
                                               f_curr_right);
                local_area = adap_quadrature(function_ptr,
                                             curr_left,
                                             f_curr_left,
                                             curr_right,
                                             f_curr_right,
                                             orig_approx);
            }
            MPI_Send((void*) &local_area,
                     1,
                     MPI_DOUBLE,
                     0,
                     0,
                     MPI_COMM_WORLD);
            MPI_Recv((void*) &received_data,
                     3,
                     MPI_DOUBLE,
                     0,
                     1,
                     MPI_COMM_WORLD,
                     MPI_STATUS_IGNORE);
        }

        clock_t proc_finish_time = clock();
        time_ver_1_per_proc[rank - 1] = (proc_finish_time - proc_start_time) / CLOCKS_PER_SEC;
    }

    MPI_Finalize();
}

int main(int argc, char** argv) {
    // Define here the function we want to integrate
    // double (*function_ptr)(double) = &inv_function;
    double (*function_ptr)(double) = &const_function;
    double left = 0.1;
    double right = 10.0;

    double version_1_area = 0.0;
    version_1(function_ptr, left, right, &version_1_area);
    printf("Area from version 1: %f\n", version_1_area);

    double version_2_area = 0.0;
    version_2(function_ptr, left, right, &version_2_area);
    printf("Area from version 2: %f\n", version_2_area);
    print_times();
}
