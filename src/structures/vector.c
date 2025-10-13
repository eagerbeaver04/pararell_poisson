#include "includes/vector.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

Vector* create_vector(void)
{
    Vector* v = (Vector*)malloc(sizeof(Vector));
    if(!v)
        return NULL;
    v->data = NULL;
    v->size = 0;
    v->capacity = 0;
    return v;
}

int reserve(Vector* v, size_t new_capacity)
{
    if(new_capacity <= v->capacity)
        return 0;

    int* new_data = (int*)realloc(v->data, new_capacity * sizeof(int));
    if(!new_data)
        return -1;

    v->data = new_data;
    v->capacity = new_capacity;
    return 0;
}

int push_back(Vector* v, int value)
{
    if(v->size >= v->capacity)
    {
        size_t new_cap = (v->capacity == 0) ? 4 : v->capacity * 2;
        if(reserve(v, new_cap) != 0)
            return -1;
    }

    v->data[v->size++] = value;
    return 0;
}

int pop_back(Vector* v)
{
    if(v->size == 0)
        return -1;
    v->size--;
    return 0;
}

int at(Vector* v, size_t index, int* value)
{
    if(index >= v->size)
        return -1;
    *value = v->data[index];
    return 0;
}

int resize(Vector* v, size_t new_size)
{
    if(new_size > v->capacity)
    {
        if(reserve(v, new_size) != 0)
            return -1;
    }
    v->size = new_size;
    return 0;
}

void clear(Vector* v) { v->size = 0; }

void clear_vector(Vector* v)
{
    free(v->data);
    free(v);
}

int dot_product(Vector* v1, Vector* v2, int* result)
{
    if(v1->size != v2->size)
    {
        return -1;
    }

    if(v1->size == 0)
    {
        *result = 0;
        return 0;
    }

    int sum = 0;
    for(size_t i = 0; i < v1->size; ++i)
    {
        sum += v1->data[i] * v2->data[i];
    }

    *result = sum;
    return 0;
}

int multiply_by_scalar(Vector* v, int scalar)
{
    if(v->size == 0)
        return 0;

    for(size_t i = 0; i < v->size; ++i)
    {
        v->data[i] *= scalar;
    }
    return 0;
}

Vector* multiply_by_scalar_copy(Vector* v, int scalar)
{
    Vector* result = create_vector();
    if(!result)
        return NULL;

    if(reserve(result, v->size) != 0)
    {
        clear_vector(result);
        return NULL;
    }

    for(size_t i = 0; i < v->size; ++i)
    {
        push_back(result, v->data[i] * scalar);
    }

    return result;
}

double norm(Vector* v)
{
    double val = 0;
    for(size_t i = 0; i < v->size; ++i)
    {
        val += pow(v->data[i], 2);
    }
    return sqrt(val);
}
