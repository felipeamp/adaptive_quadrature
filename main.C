#include <math.h>
#include <stdio.h>
#include <time.h>

#include <mpi.h>

// Define colors for printf use
#define RED     "\x1b[31m"
#define GREEN   "\x1b[32m"
#define YELLOW  "\x1b[33m"
#define RESET "\x1b[0m"

#define DEBUG 0

///////////////////////////////////////////////////////////////////////////////////////////////////
////                                                                                           ////
////                                 CONFIGURATION VARIABLES                                   ////
////                                                                                           ////
///////////////////////////////////////////////////////////////////////////////////////////////////

#define NUM_WORKER_PROCS 4 //2 4 8
#define NUM_INTERVALS 100 // 25 50 75 125 150 175
#define TAU 0.001 // Maximum allowed error ratio between new_trap_1 + new_trap_2 and old_trap

double inv_function(double x) {
    return 1.0/x;
}

double const_function(double x) {
    return 1.0;
}

static double (*function_ptr)(double) = &inv_function;
//static double (*function_ptr)(double) = &const_function;

static double left = 0.001;
static double right = 10.0;



///////////////////////////////////////////////////////////////////////////////////////////////////
////                                                                                           ////
////                                   AUXILIARY FUNCTIONS                                     ////
////                                                                                           ////
///////////////////////////////////////////////////////////////////////////////////////////////////

static double area_found_ver_1 = 0.0;
static double area_found_ver_2 = 0.0;

static double total_time_ver_1;
static double total_time_per_proc_ver_1[NUM_WORKER_PROCS];
static double total_time_ver_2;
static double total_time_per_proc_ver_2[NUM_WORKER_PROCS];

void print_times() {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank == 0) {
        printf("\n============================================================================\n");
        printf("\nNUM_WORKER_PROCS = %d\n", NUM_WORKER_PROCS);
        printf("NUM_INTERVALS = %d\n\n", NUM_INTERVALS);

        printf("Version 1\n");
        printf("\tTotal area found: %f\n", area_found_ver_1);
        printf("\tTime taken:\n");
        printf("\t\tMaster process:    %fs\n", total_time_ver_1);
        for (int worker_proc_index = 0;
             worker_proc_index < NUM_WORKER_PROCS;
             ++worker_proc_index) {
            printf("\t\tWorker process #%d: %fs\n",
                   worker_proc_index,
                   total_time_per_proc_ver_1[worker_proc_index]);
        }

        printf("Version 2\n");
        printf("\tTotal area found: %f\n", area_found_ver_2);
        printf("\tTime taken:\n");
        printf("\t\tMaster process:    %fs\n", total_time_ver_2);
        for (int worker_proc_index = 0;
             worker_proc_index < NUM_WORKER_PROCS;
             ++worker_proc_index) {
            printf("\t\tWorker process #%d: %fs\n",
                   worker_proc_index,
                   total_time_per_proc_ver_2[worker_proc_index]);
}
        printf("\n============================================================================\n");
    }
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
                       double prev_approx,
                       int rank) {
    double center = 0.5 * (left + right);
    double f_center = (*function_ptr)(center);
    double approx_left_half = trap_area(left, f_left, center, f_center);
    double approx_right_half = trap_area(center, f_center, right, f_right);
    if ((prev_approx == 0.0 && (approx_left_half != 0.0 || approx_right_half != 0.0))
        || error_ratio(approx_left_half, approx_right_half, prev_approx) > TAU) {
        #if DEBUG
            printf(RED "\tWorker proc %d needed an extra level of recursion to calculate the "
                   "integral on the interval [%f, %f].\n" RESET, rank, left, right);
        #endif
        double left_area = adap_quadrature(function_ptr,
                                           left,
                                           f_left,
                                           center,
                                           f_center,
                                           approx_left_half,
                                           rank);
        double right_area = adap_quadrature(function_ptr,
                                            center,
                                            f_center,
                                            right,
                                            f_right,
                                            approx_right_half,
                                            rank);
        return left_area + right_area;
    } else {
        return approx_left_half + approx_right_half;
    }
}



///////////////////////////////////////////////////////////////////////////////////////////////////
////                                                                                           ////
////                                        VERSION 1                                          ////
////                                                                                           ////
///////////////////////////////////////////////////////////////////////////////////////////////////

void master_process_ver_1(double left, double right) {
    clock_t start_time = clock();

    double local_area = 0.0; // this variable exists only to be passed to MPI_Reduced
    double interval_step_size = (right - left) / (double)NUM_WORKER_PROCS;

    // Send a task to each process
    for (int worker_proc_number = 0;
         worker_proc_number < NUM_WORKER_PROCS;
         ++worker_proc_number) {
        double send_data[2];
        send_data[0] = left + interval_step_size * worker_proc_number;
        send_data[1] = left + interval_step_size * (worker_proc_number + 1);
        MPI_Send((void*) &send_data,
                 2,
                 MPI_DOUBLE,
                 worker_proc_number + 1,
                 0,
                 MPI_COMM_WORLD);

        #if DEBUG
            printf(GREEN "Sent interval [%f, %f] to worker_proc_number %d.\n" RESET,
                   send_data[0],
                   send_data[1],
                   worker_proc_number);
        #endif
    }

    MPI_Reduce((void*) &local_area,
       (void*) &area_found_ver_1,
       1,
       MPI_DOUBLE,
       MPI_SUM,
       0,
       MPI_COMM_WORLD);

    clock_t finish_time = clock();
    total_time_ver_1 = (double)(finish_time - start_time) / CLOCKS_PER_SEC;

    for (int worker_proc_number = 0; worker_proc_number < NUM_WORKER_PROCS; ++worker_proc_number) {
        MPI_Recv((void*) &total_time_per_proc_ver_1[worker_proc_number],
                 1,
                 MPI_DOUBLE,
                 1 + worker_proc_number,
                 MPI_ANY_TAG,
                 MPI_COMM_WORLD,
                 MPI_STATUS_IGNORE);
    }
}

void worker_process_ver_1(double (*function_ptr)(double), int rank) {
    double received_data[2];
    MPI_Recv((void*) &received_data,
             2,
             MPI_DOUBLE,
             0,
             MPI_ANY_TAG,
             MPI_COMM_WORLD,
             MPI_STATUS_IGNORE);

    clock_t proc_start_time = clock();

    double curr_left = received_data[0];
    double curr_right = received_data[1];

    #if DEBUG
        printf(YELLOW "Process %d received interval [%f, %f] from process 0.\n" RESET,
               rank,
               curr_left,
               curr_right);
    #endif

    // Initial area approximation
    double f_curr_left = (*function_ptr)(curr_left);
    double f_curr_right = (*function_ptr)(curr_right);
    double orig_approx = trap_area(curr_left,
                                   f_curr_left,
                                   curr_right,
                                   f_curr_right);
    // Area calculation
    double local_area = adap_quadrature(function_ptr,
                                        curr_left,
                                        f_curr_left,
                                        curr_right,
                                        f_curr_right,
                                        orig_approx,
                                        rank);

    clock_t proc_finish_time = clock();
    double time_taken = (double)(proc_finish_time - proc_start_time) / (double)CLOCKS_PER_SEC;

    MPI_Reduce((void*) &local_area,
               NULL,
               1,
               MPI_DOUBLE,
               MPI_SUM,
               0,
               MPI_COMM_WORLD);

    MPI_Send((void*) &time_taken,
             1,
             MPI_DOUBLE,
             0,
             0,
             MPI_COMM_WORLD);
}

void version_1(double (*function_ptr)(double), double left, double right) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) {
        master_process_ver_1(left, right);
    } else {
        worker_process_ver_1(function_ptr, rank);
    }
}



///////////////////////////////////////////////////////////////////////////////////////////////////
////                                                                                           ////
////                                        VERSION 2                                          ////
////                                                                                           ////
///////////////////////////////////////////////////////////////////////////////////////////////////

void master_process_ver_2(double left, double right) {
    clock_t start_time = clock();

    double interval_step_size = (right - left) / (double)NUM_INTERVALS;

    // At the beginning each worker process will send 0.0 to the master process. Thus, there will
    // be NUM_WORKER_PROCS spurious communications at the beginning.
    int intervals_done = -NUM_WORKER_PROCS;
    int tasks_sent = 0;
    int flag, source;
    MPI_Status status;
    double received_data;
    double send_data[3]; // {bool_continue_working, left, right}

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
            area_found_ver_2 += received_data;

            #if DEBUG
                printf(GREEN "Worker %d sent local_area = %f.\n" RESET, source, received_data);
            #endif

            if (tasks_sent < NUM_INTERVALS) {
                // There is still some task to be done.
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
                // All tasks have been send to a process. Hence, end this worker process
                send_data[0] = 0.0; // bool_continue_working
                MPI_Send((void*) &send_data,
                         3,
                         MPI_DOUBLE,
                         source,
                         1,
                         MPI_COMM_WORLD);
                MPI_Recv((void*) &total_time_per_proc_ver_2[source - 1],
                         1,
                         MPI_DOUBLE,
                         source,
                         MPI_ANY_TAG,
                         MPI_COMM_WORLD,
                         MPI_STATUS_IGNORE);
            }
        }
    }

    clock_t finish_time = clock();
    total_time_ver_2 = (double)(finish_time - start_time) / CLOCKS_PER_SEC;
}

void worker_process_ver_2(double (*function_ptr)(double), int rank) {
    clock_t proc_start_time = clock();

    double local_area = 0.0;
    double received_data[3] = {1.0, 0.0, 0.0};
    double curr_left, curr_right;

    while(received_data[0]) {
        curr_left = received_data[1];
        curr_right = received_data[2];

        if (curr_left != curr_right) {
            // We added the `if` above just to make sure we don't try to calculate the function
            // at the origin in case it is not well-defined. That could happen only at the
            // beginning of a run.

            #if DEBUG
                printf(YELLOW "Process %d received interval [%f, %f] from process 0.\n" RESET,
                       rank,
                       curr_left,
                       curr_right);
            #endif

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
                                         orig_approx,
                                         rank);
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
    double total_time = (double)(proc_finish_time - proc_start_time) / CLOCKS_PER_SEC;
    MPI_Send((void*) &total_time,
         1,
         MPI_DOUBLE,
         0,
         0,
         MPI_COMM_WORLD);
}

void version_2(double (*function_ptr)(double), double left, double right) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) {
        master_process_ver_2(left, right);
    } else {
        worker_process_ver_2(function_ptr, rank);
    }
}



///////////////////////////////////////////////////////////////////////////////////////////////////
////                                                                                           ////
////                                           MAIN                                            ////
////                                                                                           ////
///////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv) {
    setbuf(stdout, NULL);
    MPI_Init(NULL, NULL);
    version_1(function_ptr, left, right);
    version_2(function_ptr, left, right);
    print_times();
    MPI_Finalize();
}
