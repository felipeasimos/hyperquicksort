#include <bits/stdint-uintn.h>
#include <omp.h>
#include <stdint.h>
#include <string.h>
#include <time.h>
#include <stdio.h>
#include <math.h>

#define debug fprintf(stderr, "---------------------\n%s:%d\n---------------------\n", __FILE__, __LINE__)

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

double quicksort_test(void* original_arr, void* arr, unsigned long element_size, unsigned long num_elements, arr_type (*cmp_func)(void* a, void* b)) {

    double quicksort_time = omp_get_wtime();
    quicksort(arr, element_size, num_elements, cmp_func);
    quicksort_time = omp_get_wtime() - quicksort_time;

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

void* scatter_array(void* arr, unsigned long element_size, unsigned long num_elements, unsigned long* local_num_elements, int num_threads, int id) {

    void* local_arr = calloc(element_size, num_elements);
    *local_num_elements = num_elements / num_threads;

    memcpy(local_arr, get_element(arr, element_size, (*local_num_elements) * id), element_size * (*local_num_elements));

    // get any rest left (if the size of the array is not divisible by the number of threads)
    unsigned long rest = num_elements % num_threads;
    if(rest && id == num_threads-1) {
        memcpy(get_element(local_arr, element_size, *local_num_elements), get_element(arr, element_size, (*local_num_elements) * num_threads), element_size * rest);
        *local_num_elements += rest;
    }


    return local_arr;
}

void gather_array_seq(void* arr, unsigned long element_size, int num_threads, unsigned long local_num_elements[num_threads], void* local_arrs[num_threads]) {

    unsigned long next_idx = 0;
    for(int i = 0; i < num_threads; i++) {

        memcpy(get_element(arr, element_size, next_idx), local_arrs[i], local_num_elements[i] * element_size);
        next_idx += local_num_elements[i];
    }
}

void broadcast_pivot(unsigned long element_size, void* local_arr, unsigned long local_num_elements, uint8_t* pivot) {
    if(local_num_elements) {
        memcpy(pivot, get_element(local_arr, element_size, local_num_elements / 2), element_size);
    } else {
        memset(pivot, 0x00, element_size);
    }
}

unsigned long shallow_quicksort(void* arr, unsigned long element_size, unsigned long num_elements, arr_type (*cmp_func)(void* a, void* b), void* pivot) {

    if(!num_elements) return 0;
    // if the first element is equal or greater than the pivot, return 0

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

void swap_with_pair(unsigned long num_elements, unsigned long element_size, int num_threads, void* local_arrs[num_threads], unsigned long local_num_elements[num_threads], unsigned long low_partitions[num_threads], int (*cmp_func)(void*,void*), int id, int pair_thread) {

    uint8_t recv_buf[ element_size * num_elements ];
    unsigned long n_recv_buf_elements;
    uint8_t merge_buf[ element_size * num_elements ];

    unsigned long high_partition = local_num_elements[id] - low_partitions[id];
    unsigned long pair_high_partition = local_num_elements[pair_thread] - low_partitions[pair_thread];
    void* original_partition_arr = local_arrs[id];

    #pragma omp barrier
    if(id > pair_thread) { // high processes
        // 1. get pair's high partition
        memcpy(recv_buf, get_element(local_arrs[pair_thread], element_size, low_partitions[pair_thread]), element_size * pair_high_partition);
        n_recv_buf_elements = pair_high_partition;
        // 2. we change the value of 'original_partition_arr' here so it will point to the partition we will use later
        original_partition_arr = get_element(local_arrs[id], element_size, low_partitions[id]);
        // 3. update 'local_num_elements' for the size of the partition from the original array
        local_num_elements[id] = high_partition;
    } else { // low processes
        // 1. receive pair's low partition
        memcpy(recv_buf, local_arrs[pair_thread], element_size * low_partitions[pair_thread]);
        n_recv_buf_elements = low_partitions[pair_thread];
        // 2; update 'local_num_elements' for the size of the partition from the original array
        local_num_elements[id] = low_partitions[id];
    }

    #pragma omp barrier
    // merge recv_buf and local_arr, putting result in 'merge_buf' and returning its size
    local_num_elements[id] = merge_arrays(merge_buf, original_partition_arr, recv_buf, local_num_elements[id], n_recv_buf_elements, element_size, cmp_func);
    // copy 'merge_buf' to 'local_arr'
    memcpy(local_arrs[id], merge_buf, local_num_elements[id] * element_size);
}

void _hyperquicksort(unsigned long num_elements, unsigned long element_size, int num_threads, void* local_arrs[num_threads], unsigned long local_num_elements[num_threads], int (*cmp_func)(void*,void*), unsigned long iter) {

    int num_threads_per_group = num_threads>>iter;
    int num_threads_per_group_halfed = num_threads_per_group>>1;
    int num_thread_groups = pow(2, iter);
    uint8_t pivot[num_thread_groups][element_size];
    unsigned long low_partitions[num_threads];

    #pragma omp parallel
    {
        int id = omp_get_thread_num();
        int group_id = id / num_threads_per_group;
        int group_root_thread_id = group_id * num_threads_per_group;

        // each group's root thread should write to the pivot
        if(id == group_root_thread_id) {
            broadcast_pivot(element_size, local_arrs[id], local_num_elements[id], (uint8_t*)&pivot[group_id]);
        }
        #pragma omp barrier
        // since our local arr is always ordered, we can just do a binary search to find the size of the low partition
        low_partitions[id] = shallow_quicksort(local_arrs[id], element_size, local_num_elements[id], cmp_func, pivot[group_id]);

        #pragma omp barrier
        // this calculation will return the right pair process (the modulo is used to deal with high id processes)
        int pair_thread = group_root_thread_id+((id+num_threads_per_group_halfed)%num_threads_per_group);

        // a merge is done here, leaving the final local_arr still sorted
        swap_with_pair(num_elements, element_size, num_threads, local_arrs, local_num_elements, low_partitions, cmp_func, id, pair_thread);
    }
}

void hyperquicksort(void* arr, unsigned long element_size, unsigned long num_elements, arr_type (*cmp_func)(void* a, void* b)) {

    int num_threads;
    #pragma omp parallel
    {
        num_threads = omp_get_num_threads();
    }
    if(num_elements <= (unsigned long)num_threads) {
        quicksort(arr, element_size, num_elements, arr_type_cmp);
        return;
    }

    void* local_arrs[num_threads];
    unsigned long local_num_elements[num_threads];
    unsigned long logp = round(log(num_threads)/log(2));
 
    #pragma omp parallel
    {
        // actual hyperquicksort algorithm
        int id = omp_get_thread_num();
        // scatter the original arr into a local_arr for each process
        local_arrs[id] = scatter_array(arr, element_size, num_elements, &local_num_elements[id], num_threads, id);

        // sort local array
        quicksort(local_arrs[id], element_size, local_num_elements[id], cmp_func);
    }

    for(unsigned long i = 0; i < logp; i++) {

        _hyperquicksort(num_elements, element_size, num_threads, local_arrs, local_num_elements, cmp_func, i);
#ifdef DEBUG
        gather_array_seq(arr, element_size, num_threads, local_num_elements, local_arrs);
        for(int j = 0; j < num_threads; j++) {
            printf("iter: %lu, id: %d, number of elements: %lu\n", i,  j, local_num_elements[j]);
        }

        for(unsigned long j = 0; j < num_elements; j++) {
            printf("%d ", *(int*)get_element(arr, element_size, j));
        }
        printf("\n");
#endif
    }

    gather_array_seq(arr, element_size, num_threads, local_num_elements, local_arrs);
    #pragma omp parallel
    {
        int id = omp_get_thread_num();
        free(local_arrs[id]);
    }
}

void hyperquicksort_test(void* original_arr, void* arr, unsigned long element_size, unsigned long num_elements, arr_type (*cmp_func)(void* a, void* b), double* time) {
    *time = omp_get_wtime();
    hyperquicksort(arr, element_size, num_elements, cmp_func);
    *time = omp_get_wtime() - *time;

#ifdef DEBUG
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
    /*
    // check that the array is equal to the original one
    for(unsigned long i = 1; i < num_elements; i++) {

        if(memcmp(get_element(arr, element_size, i), get_element(original_arr, element_size, i), element_size) ) {

            printf("ERROR in hyperquicksort at index %lu!\n", i);
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
    */
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
    unsigned long num_threads;

    #pragma omp parallel
    {
        if(omp_get_thread_num() == 0) {
            num_threads = omp_get_num_threads();
        }
    }

    sprintf(filename, "results_openmp_%d.plt", (int)num_threads);
    results_file = fopen(filename, "w");

    for(unsigned long size = min_arr_size; size <= max_arr_size; size += step_size) {
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

        for(unsigned long i = 0; i < n_tests_per_size; i++) {

            memcpy(hyperquicksort_copy, arr, size * sizeof(arr_type));
            hyperquicksort_test(arr, hyperquicksort_copy, sizeof(arr_type), size, arr_type_cmp, &tmp_time);
            hyperquicksort_time += tmp_time;
        }


        free(arr);
        free(quicksort_copy);
        free(hyperquicksort_copy);

        hyperquicksort_time /= n_tests_per_size;
        quicksort_time /= n_tests_per_size;
        fprintf(results_file, "%lu %.12F %.12F\n", size, quicksort_time, hyperquicksort_time);
        fflush(results_file);
    }
}

int main() {

    test();
    return 0;
}
