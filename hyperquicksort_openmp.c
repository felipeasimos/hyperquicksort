#include <bits/stdint-uintn.h>
#include <omp.h>
#include <stdint.h>
#include <string.h>

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

int main() {

    return 0;
}
