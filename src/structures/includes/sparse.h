#ifndef SPARSE_HEADER_H
#define SPARSE_HEADER_H

#include "list.h"

struct Sparse;

typedef struct Sparse Sparse;

struct Sparse
{
    size_t rows;
    size_t columns;
    List** row_lists;
};

void add_node_to_sparse(Sparse* sparse, Node* node, size_t row);
void print_sparse(Sparse*);
Sparse* create_sparse(size_t rows, size_t columns);
void clear_sparse(Sparse* sparse);

#endif // SPARSE_HEADER_H
