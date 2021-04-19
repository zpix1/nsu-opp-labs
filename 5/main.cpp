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
    while(1) {
        no_tasks_left.lock();
    }
}

void worker() {

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