# adaptive_quadrature
Solution to the 2nd project of the "Concurrent and Parallel Programming" course.


Problem Specification
---------------------

Use *mpich* to implement the Trapezoidal Rule, where we approximate the integral of a function in a given interval by the sum of trapeziums formed by points xi, f(xi), xi+1 and f(xi+1).

The function must be reasonably complex to think about parallelizing its integral calculation. The program must easily accept different functions and intervals.

The program must have two versions:

1. The master process partitions the total interval in p intervals, where p is the number of executed processes. Each process receives one interval and calculates its area using the adaptive quadrature. By the end of calculation, the processes execute a reduction operation to obtain the sum of the result of each one of them. Try to find functions that cause an unbalanced load between processes, with very different number of iterations needed to convergence on each interval. Without this feature, the second version will always be much worse than this one.

2. The master process will initially create a list of tasks, containing the intervals extrema (for instance, 100 tasks). Each thread executes a task until it finishes and asks the master process for a new task, sending in this request the result of the previous task. That happens until there are no more tasks. The main process waits until all process hand over their final results.

Execute both versions with different number of processes (2, 4 and 8) and different numbers of tasks in the initial pool, doing a few time measurements for each combination.


How to use and Configure
------------------------

To change the number of worker processes, number of intervals (for version 2) and maximum error ratio allowed, just change the global variables NUM_WORKER_PROCS, NUM_INTERVALS and TAU at the beginning of the code. To change the functio or interval of integration, change the first lines of the main function.

To compile and execute:

`mpicc main.C -o main.out`

`mpiexec -np <NUM_WORKER_PROCS + 1> main.out`
