#include <bits/stdint-uintn.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdint.h>
#include <mpi.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

typedef int arr_type;

int arr_type_cmp(void* a, void* b) {

    if(*(arr_type*)a < *(arr_type*)b) {
        return 1;
    } else if(*(arr_type*)a > *(arr_type*)b) {
        return -1;
    } else {
        return 0;
    }
}

void swap(void* a, void* b, unsigned long element_size) {

    if(a==b) return;

    uint8_t tmp[element_size];
    memcpy(tmp, a, element_size);
    memcpy(a, b, element_size);
    memcpy(b, tmp, element_size);
}

void* get_element(void* arr, unsigned long element_size, unsigned long idx) {

    return (uint8_t*)arr+(element_size * idx);
}

// cmp_func: 1 if a < b, 0 if a == b, -1 if a > b
void quicksort(void* arr, unsigned long element_size, unsigned long num_elements, arr_type (*cmp_func)(void* a, void* b)) {

    if( num_elements < 2 ) return;

    unsigned long lower_index = 0;
    for(unsigned long i = 1; i < num_elements; i++) {
        void* current_element = get_element(arr, element_size, i);
        switch( cmp_func(current_element, arr) ) {
            case 1: 
                swap(current_element, get_element(arr, element_size, ++lower_index), element_size);
                break;
            case 0:
            case -1:
                break;
        }
    }

    swap(get_element(arr, element_size, lower_index), arr, element_size);

    quicksort(arr, element_size, lower_index, cmp_func);
    quicksort(get_element(arr, element_size, lower_index + 1), element_size, num_elements - lower_index - 1, cmp_func);
}

unsigned long shallow_quicksort(void* arr, unsigned long element_size, unsigned long num_elements, arr_type (*cmp_func)(void* a, void* b), void* pivot) {

    if(!num_elements) return 0;
    // if the first element is equal or greater than the pivot, return 0
    if(cmp_func(pivot, arr) != -1) return 0;

    // make a binary search to find the index of the first element that the pivot is smaller
    unsigned long first = 0;
    unsigned long last = num_elements - 1;

    while(first != last) {

        unsigned long median_idx = (first + last) >> 1;
        switch( cmp_func(get_element(arr, element_size, median_idx), pivot) ) {

            // if median is smaller than the pivot
            case 1:
                // if this is true, there is only the 'first' and 'last' elements left
                if(first == median_idx) {
                    first = last;
                } else {
                    first = median_idx;
                }
                break;
            // if median is equal or greater than the pivot
            case 0:
            case -1:
                last = median_idx;
                break;
        }
    }
    return first + (cmp_func(get_element(arr, element_size, first), pivot) == 1);
}

void* scatter_array(void* arr, unsigned long element_size, unsigned long num_elements, unsigned long* local_num_elements, int comm_size, int rank) {

    void* local_arr = calloc(element_size, num_elements);
    *local_num_elements = num_elements / comm_size;
    unsigned long total_bytes = (*local_num_elements) * element_size;
    MPI_Scatter(arr, total_bytes, MPI_UINT8_T, local_arr, total_bytes, MPI_UINT8_T, 0, MPI_COMM_WORLD);

    // get any rest left (if the size of the array is not divisible by the number of processors)
    unsigned long rest = num_elements % comm_size;
    if(rest && (!rank || rank == comm_size - 1)) {
        total_bytes = element_size * rest;
        if(!rank) {
            MPI_Send(get_element(arr, element_size, num_elements - rest), total_bytes, MPI_UINT8_T, comm_size - 1, 0, MPI_COMM_WORLD);
        } else if(rank == comm_size - 1) {
            MPI_Recv(get_element(local_arr, element_size, *local_num_elements), total_bytes, MPI_UINT8_T, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            *local_num_elements += rest;
        }
    }
    return local_arr;
}

void gather_array(void* arr, unsigned long element_size, unsigned long local_num_elements, void* local_arr, int comm_size, int rank) {

    // get all local_arrs back into arr
    if(!rank) {

        memcpy(arr, local_arr, element_size * local_num_elements);
        // get pointer to where we will be writing next
        void* next_partition_write = get_element(arr, element_size, local_num_elements);
        for(unsigned long i = 1; i < (unsigned long)comm_size; i++) {

            // get partition size
            MPI_Recv(&local_num_elements, 1, MPI_UNSIGNED_LONG, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            if(local_num_elements) {
                // get partition
                MPI_Recv(next_partition_write, local_num_elements * element_size, MPI_UINT8_T, i, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                next_partition_write = ((uint8_t*)next_partition_write) + (local_num_elements * element_size);
            }
        }
    } else {
        MPI_Send(&local_num_elements, 1, MPI_UNSIGNED_LONG, 0, 0, MPI_COMM_WORLD);
        if(local_num_elements)
            MPI_Send(local_arr, local_num_elements * element_size, MPI_UINT8_T, 0, 1, MPI_COMM_WORLD);
    }
}

void broadcast_pivot(void* arr, unsigned long element_size, unsigned long num_elements, int rank, MPI_Comm comm, void* pivot) {
    if(!rank) {
        if(num_elements) {
            memcpy(pivot, get_element(arr, element_size, num_elements / 2), element_size);
        } else {
            memset(pivot, 0x00, element_size);
        }
    }
    MPI_Bcast(pivot, element_size, MPI_UINT8_T, 0, comm);
}

unsigned long merge_arrays(void* result, void* a, void* b, unsigned long a_size, unsigned long b_size, unsigned long element_size, int (*cmp_func)(void*,void*)) {

    unsigned long a_idx = 0;
    unsigned long b_idx = 0;
    unsigned long result_size = a_size + b_size;
    for(unsigned long result_idx = 0; result_idx < result_size; result_idx++) {

        int comp = 0;
        if( a_idx >= a_size ) {
            comp = -1;
        } else if( b_idx >= b_size ){
            comp = 1;
        } else {
            comp = cmp_func( get_element(a, element_size, a_idx), get_element(b, element_size, b_idx) );
        }
        switch( comp ) {

            // b < a
            case -1:
                memcpy( get_element(result, element_size, result_idx), get_element(b, element_size, b_idx++), element_size );
                break;
            // a < b
            case 0:
            case 1:
                memcpy( get_element(result, element_size, result_idx), get_element(a, element_size, a_idx++), element_size );
                break;
        }
    }
    return result_size;
}

void swap_with_pair(void* local_arr, unsigned long element_size, unsigned long* local_num_elements, unsigned long low_partition, unsigned long num_elements, int (*cmp_func)(void*,void*), int rank, int pair, MPI_Comm new_comm) {

    uint8_t recv_buf[ element_size * num_elements ];
    uint8_t merge_buf[ element_size * num_elements ];

    unsigned long high_partition = (*local_num_elements) - low_partition;
    void* original_partition_arr = local_arr;
    // we can get the number of elements we received from this status struct
    MPI_Status status;

    if(rank > pair) { // high processes
        // send lower partition
        MPI_Send(local_arr, low_partition * element_size, MPI_UINT8_T, pair, 0, new_comm);
        // receive pair's high partition (number of bytes is written to 'status')
        MPI_Recv(recv_buf, num_elements * element_size, MPI_UINT8_T, pair, 0, new_comm, &status);
        // we change the value of 'original_partition_arr' here so it will point to the partition we will use later
        original_partition_arr = get_element(local_arr, element_size, low_partition);
        // update 'local_num_elements' for the size of the partition from the original array
        *local_num_elements = high_partition;
    } else { // low processes
        // receive pair's high partition (number of bytes is written to 'status')
        MPI_Recv(recv_buf, num_elements * element_size, MPI_UINT8_T, pair, 0, new_comm, &status);
        // send higher partition
        MPI_Send(get_element(local_arr, element_size, low_partition), high_partition * element_size, MPI_UINT8_T, pair, 0, new_comm);
        // update 'local_num_elements' for the size of the partition from the original array
        *local_num_elements = low_partition;
    }

    // get number of elements in the recv_buf (we need to divide by the element_size, since what we get is the number of bytes from 'status')
    int recv_buf_elements;
    MPI_Get_count(&status, MPI_UINT8_T, &recv_buf_elements);
    recv_buf_elements /= element_size;

    // merge recv_buf and local_arr, putting result in 'merge_bugf' and returning its size
    *local_num_elements = merge_arrays(merge_buf, original_partition_arr, recv_buf, *local_num_elements, recv_buf_elements, element_size, cmp_func);
    // copy 'merge_buf' to 'local_arr'
    memcpy(local_arr, merge_buf, (*local_num_elements) * element_size);
}

void _hyperquicksort(void* local_arr, unsigned long element_size, unsigned long* local_num_elements, unsigned long num_elements, int (*cmp_func)(void*,void*), int comm_size, int world_rank, unsigned long iter) {

    int color = pow(2, iter)*world_rank/comm_size;

    // get a new communicator and rank, dividing the process into groups (16, 8, 4, 2, ...)
    MPI_Comm new_comm;
    MPI_Comm_split(MPI_COMM_WORLD, color, world_rank, &new_comm);
    int rank, sz;
    MPI_Comm_rank(new_comm, &rank);
    MPI_Comm_size(new_comm, &sz);

    // broadcast pivot to all process in our group
    uint8_t pivot[element_size];
    broadcast_pivot(local_arr, element_size, *local_num_elements, rank, new_comm, pivot);

    // since our local arr is always ordered, we can just do a binary search to find the size of the low partition
    unsigned long low_partition = shallow_quicksort(local_arr, element_size, *local_num_elements, cmp_func, pivot);
    unsigned long sz_halfed = sz >> 1;

    // this calculation will return the right pair process (the modulo is used to deal with high id processes)
    int pair_process = (rank+sz_halfed)%sz;
    // a merge is done here, leaving the final local_arr still sorted
    swap_with_pair(local_arr, element_size, local_num_elements, low_partition, num_elements, cmp_func, rank, pair_process, new_comm);

    MPI_Comm_free(&new_comm);
}

void hyperquicksort(void* arr, unsigned long element_size, unsigned long num_elements, arr_type (*cmp_func)(void* a, void* b), int comm_size, int world_rank) {

    if(num_elements <= (unsigned long)comm_size) {
        if(!world_rank) quicksort(arr, element_size, num_elements, arr_type_cmp);
        return;
    }

    // scatter the original arr into a local_arr for each process
    unsigned long local_num_elements;
    void* local_arr = scatter_array(arr, element_size, num_elements, &local_num_elements, comm_size, world_rank);

    // sort local array
    quicksort(local_arr, element_size, local_num_elements, cmp_func);

    unsigned long logp = round(log(comm_size)/log(2));
    for(unsigned long i = 0; i < logp; i++) {

        _hyperquicksort(local_arr, element_size, &local_num_elements, num_elements, cmp_func, comm_size, world_rank, i);

#ifdef DEBUG
        MPI_Barrier(MPI_COMM_WORLD);
        gather_array(arr, element_size, local_num_elements, local_arr, comm_size, world_rank);
        printf("iter: %lu, rank %d, number of elements: %lu\n",i,  world_rank, local_num_elements);
        MPI_Barrier(MPI_COMM_WORLD);
        if(!world_rank) {
            for(unsigned long i = 0; i < num_elements; i++) {
                printf("%d ", *(int*)get_element(arr, element_size, i));
            }
            printf("\n");
        }
        MPI_Barrier(MPI_COMM_WORLD);
#endif
    }

    gather_array(arr, element_size, local_num_elements, local_arr, comm_size, world_rank);
    free(local_arr);
}

double quicksort_test(void* original_arr, void* arr, unsigned long element_size, unsigned long num_elements, arr_type (*cmp_func)(void* a, void* b)) {

    double quicksort_time = MPI_Wtime();
    quicksort(arr, element_size, num_elements, cmp_func);
    quicksort_time = MPI_Wtime() - quicksort_time;

#ifdef DEBUG
    // check that the array is indeed sorted
    for(unsigned long i = 1; i < num_elements; i++) {

        if(cmp_func(get_element(arr, element_size, i-1), get_element(arr, element_size, i)) == -1) {

            printf("ERROR in quicksort in index %lu and %lu!\n", i-1, i);
            printf("original arr:\n");
            for(unsigned long j = 0; j < num_elements; j++) {
                printf("%d ", *(arr_type*)get_element(original_arr, element_size, j));
            }
            printf("\n");
            printf("arr:\n");
            for(unsigned long j = 0; j < num_elements; j++) {
                printf("[%lu]=%d ", j, *(arr_type*)get_element(arr, element_size, j));
            }
            printf("\n");
            exit(1);
        }
    }
#endif
    return quicksort_time;
}

void hyperquicksort_test(void* original_arr, void* arr, unsigned long element_size, unsigned long num_elements, arr_type (*cmp_func)(void* a, void* b), double* time, int comm_size, int rank) {
    MPI_Barrier(MPI_COMM_WORLD);
    if(!rank) {
        *time = MPI_Wtime();
    }
    hyperquicksort(arr, element_size, num_elements, cmp_func, comm_size, rank);
    if(!rank) {
        *time = MPI_Wtime() - *time;
    }

#ifdef DEBUG
    if(!rank) {
        // check that the array is indeed sorted
        for(unsigned long i = 1; i < num_elements; i++) {

            if(cmp_func(get_element(arr, element_size, i-1), get_element(arr, element_size, i)) == -1) {

                printf("ERROR in hyperquicksort between %lu and %lu!\n", i-1, i);
                printf("original arr:\n");
                for(unsigned long j = 0; j < num_elements; j++) {
                    printf("%d ", *(arr_type*)get_element(original_arr, element_size, j));
                }
                printf("\n");
                printf("arr:\n");
                for(unsigned long j = 0; j < num_elements; j++) {
                    printf("[%lu]=%d ", j, *(arr_type*)get_element(arr, element_size, j));
                }
                printf("\n");
                exit(1);
            }
        }
    }
#endif
}

void test() {

    unsigned long min_arr_size, max_arr_size, n_tests_per_size, step_size;
    arr_type* arr = NULL;
    arr_type* quicksort_copy = NULL;
    arr_type* hyperquicksort_copy = NULL;

    double quicksort_time;
    double hyperquicksort_time;
    double tmp_time;

    FILE* results_file = NULL;

    int comm_size, rank;
    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if(!rank) {
        srand(time(0));
        printf("mininum array size:\n");
        scanf("%lu", &min_arr_size);
        printf("maximum array size:\n");
        scanf("%lu", &max_arr_size);
        printf("number of tests per size:\n");
        scanf("%lu", &n_tests_per_size);
        printf("step size:\n");
        scanf("%lu", &step_size);
        char filename[50];
        sprintf(filename, "results_mpi_%d.plt", comm_size);
        results_file = fopen(filename, "w");
    }
    MPI_Bcast(&min_arr_size, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
    MPI_Bcast(&max_arr_size, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
    MPI_Bcast(&n_tests_per_size, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
    MPI_Bcast(&step_size, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);

    for(unsigned long size = min_arr_size; size <= max_arr_size; size += step_size) {
        if(!rank) {
            arr = calloc(sizeof(arr_type), size);
            quicksort_copy = calloc(sizeof(arr_type), size);
            hyperquicksort_copy = calloc(sizeof(arr_type), size);

            for(unsigned long i = 0; i < size; i++) {
                arr[i] = rand() - RAND_MAX/2;
            }

            quicksort_time = 0;
            hyperquicksort_time = 0;
            for(unsigned long i = 0; i < n_tests_per_size; i++) {
                memcpy(quicksort_copy, arr, size * sizeof(arr_type));
                quicksort_time += quicksort_test(arr, quicksort_copy, sizeof(arr_type), size, arr_type_cmp);
            }
        }

        for(unsigned long i = 0; i < n_tests_per_size; i++) {

            if(!rank) {
                memcpy(hyperquicksort_copy, arr, size * sizeof(arr_type));
            }
            hyperquicksort_test(arr, hyperquicksort_copy, sizeof(arr_type), size, arr_type_cmp, &tmp_time, comm_size, rank);
            if(!rank) {
                hyperquicksort_time += tmp_time;
            }
        }

        if(!rank) {

            free(arr);
            free(quicksort_copy);
            free(hyperquicksort_copy);

            hyperquicksort_time /= n_tests_per_size;
            quicksort_time /= n_tests_per_size;
            fprintf(results_file, "%lu %.12F %.12F\n", size, quicksort_time, hyperquicksort_time);
            fflush(results_file);
        }
    }
    MPI_Finalize();
}

int main() {

    test();
    /*int arr[] = {730547560, -226810937, 607950954, 640895092, 884005970, -649503488, -353856437, 576018669, -477225174, 115899598};
    MPI_Init(NULL, NULL);
    int comm_size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    hyperquicksort(arr, sizeof(arr[0]), sizeof(arr)/sizeof(arr[0]), int_cmp, comm_size, rank);

    if(!rank) {

        for(unsigned long i = 0; i < sizeof(arr)/sizeof(arr[0]); i++) {

            printf("%d ", arr[i]);
        }
        printf("\n");
    }

    MPI_Finalize();*/
    /*int arr[] = {23, 30, 62, 82};
    quicksort(arr, sizeof(arr[0]), sizeof(arr)/sizeof(arr[0]), arr_type_cmp);
    int pivot = 68;
    printf("%lu\n", shallow_quicksort(arr, sizeof(arr[0]), sizeof(arr)/sizeof(arr[0]), arr_type_cmp, &pivot));*/
}
