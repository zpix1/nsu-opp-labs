#include <mpi.h>
#include <iostream>
#include <thread>
#include <mutex>

int p_rank;

std::mutex list_mutex;
std::mutex no_tasks_left;

std::vector<int> task_list;

void task_distributor() {
    while (1) {
        // MPI_Recv(...)
    }
}

void receiver() {
    for (int iteration = 0; iteration < 100; iteration++) {
        list_mutex.lock();
        for (int i = 0; i < 10; i++) {
            task_list.push_back(rand());
        }
        list_mutex.unlock();

        while (1) {
            no_tasks_left.lock();
            list_mutex.lock();
            // MPI_Send()
            // MPI_Recv()
            list_mutex.unlock();
            no_tasks_left.unlock();
        }
    }
}

void worker() {
    while (true) {
        no_tasks_left.lock();
        while (task_list.size() != 0) {
            list_mutex.lock();
            int task = task_list.pop();
            list_mutex.unlock();
            do_task(task);
        }
        no_tasks_left.unlock();
    }
}

void do_task(int task) {
    // 
}

int main(int argc, char** argv) {
    int provided = 0;

    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
    MPI_Comm_rank(MPI_COMM_WORLD, &p_rank);
    
    if (provided != MPI_THREAD_MULTIPLE) {
        printf("provided is not MPI_THREAD_MULTIPLE\n");
        exit(1);
    }
    
    // first process distributes tasks
    if (p_rank == 0) {
        // std::thread task_distributor_thread(task_distributor);
        exit(1);
    } else {

    }

    MPI_Finalize();
    
    return 0;
}