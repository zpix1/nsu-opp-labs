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
    while (1) {
        list_mutex.lock();
        // MPI_Send()
        // MPI_Recv()
        list_mutex.unlock();
    }
}

void worker() {
    while (true) {
        // TODO: add list size mutex
        while (task_list.size() != 0) {
            list_mutex.lock();
            int task = task_list.pop();
            list_mutex.unlock();
            do_task(task);
        }
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