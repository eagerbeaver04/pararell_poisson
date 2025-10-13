#ifndef VECTOR_HEADER_H
#define VECTOR_HEADER_H

#include <stdio.h>

typedef struct
{
    int* data;
    size_t size;
    size_t capacity;
} Vector;

Vector* create_vector(void);

int reserve(Vector* v, size_t new_capacity);
int push_back(Vector* v, int value);
int pop_back(Vector* v);
int at(Vector* v, size_t index, int* value);
int resize(Vector* v, size_t new_size);
void clear(Vector* v);
void clear_vector(Vector* v);
int dot_product(Vector* v1, Vector* v2, int* result);
int multiply_by_scalar(Vector* v, int scalar);
Vector* multiply_by_scalar_copy(Vector* v, int scalar);

#endif
