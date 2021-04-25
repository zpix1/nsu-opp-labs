#include <mpi.h>

#include <iostream>
#include <thread>
#include <vector>
#include <mutex>
#include <cmath>
#include <cassert>

int p_rank;
int p_count;

std::mutex list_mutex;
std::mutex no_tasks_left;
std::mutex end_iteration;

std::vector<int> task_list;

bool all_done = false;

#define ITER_COUNT 10
#define TASKS_PER_ITER 144
// TASKS_PER_ITER % p_count == 0
#define ROOT 0
#define MESSAGE_END_ITERATION -1
#define MESSAGE_WHAT_TO_DO -2

#define DEBUG(var) \
            do { std::cout << p_rank << " has " << #var << ": " << var << std::endl; } while (0)

int get_random_task() {
    return rand() % 1000000;
}

void task_distributor() {
    for (int iter = 0; iter < ITER_COUNT; iter++) {
        DEBUG(iter);
        for (int i = 0; i < TASKS_PER_ITER; i++) {
            int message;
            MPI_Status status;
            MPI_Recv(&message, 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
            if (message > 0) {
                DEBUG(message);
                DEBUG(status.MPI_SOURCE);
            } else
            if (message == MESSAGE_WHAT_TO_DO) {
                int send_message = get_random_task() * (status.MPI_SOURCE + 1);
                MPI_Send(&send_message, 1, MPI_INT, status.MPI_SOURCE, 0, MPI_COMM_WORLD);
            } else {
                // no other message types allowed
                assert(false);
            }
        }
        int end_iter_sent_count = 0;
        while (end_iter_sent_count < p_count) {
            int message;
            MPI_Status status;
            MPI_Recv(&message, 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);

            if (message == MESSAGE_WHAT_TO_DO) {
                int send_message = MESSAGE_END_ITERATION;
                MPI_Send(&send_message, 1, MPI_INT, status.MPI_SOURCE, 0, MPI_COMM_WORLD);
            } else {
                assert(false);
            }
            end_iter_sent_count++;
        }
    }
}

void receiver() {
    for (int iter = 0; iter < ITER_COUNT; iter++) {
        while (1) {
            no_tasks_left.lock();
            if (task_list.size() > 0) {
                no_tasks_left.unlock();
                continue;
            }
            list_mutex.lock();

            // Request work
            int send_message = MESSAGE_WHAT_TO_DO;
            MPI_Send(&send_message, 1, MPI_INT, ROOT, 0, MPI_COMM_WORLD);

            // Listen for answer
            int message;
            MPI_Recv(&message, 1, MPI_INT, ROOT, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            // It's a task
            if (message > 0) {
                task_list.push_back(message);
            } else if (message == MESSAGE_END_ITERATION) {
                list_mutex.unlock();
                no_tasks_left.unlock();
                break;
            }

            list_mutex.unlock();
            no_tasks_left.unlock();
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
    all_done = true;
}


void do_task(int task) {
    double ans = 0.0;
    for (int i = 0; i < task; i++) {
        ans += sin(i);
    }
}

void worker() {
    int tasks_done = 0;
    while (!all_done) {
        no_tasks_left.lock();
        while (task_list.size() != 0) {
            list_mutex.lock();
            int task = task_list.back();
            task_list.pop_back();
            list_mutex.unlock();
            tasks_done++;
            do_task(task);
        }
        no_tasks_left.unlock();
    }
    std::cout << p_rank << " has done " << tasks_done << std::endl;
}

int main(int argc, char** argv) {
    int provided = 0;

    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
    MPI_Comm_rank(MPI_COMM_WORLD, &p_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &p_count);

    // checks
    assert(TASKS_PER_ITER % p_count == 0);

    
    if (provided != MPI_THREAD_MULTIPLE) {
        printf("provided is not MPI_THREAD_MULTIPLE\n");
        exit(1);
    }
    
    
    // start distributor on root process
    std::thread task_distributor_thread;
    
    if (p_rank == ROOT) {
        task_distributor_thread = std::thread(task_distributor);
    }

    std::thread receiver_thread(receiver);
    std::thread worker_thread(worker);
    
    if (p_rank == 0) {
        task_distributor_thread.join();
    }

    receiver_thread.join();
    worker_thread.join();

    MPI_Finalize();
    
    return 0;
}