#include <stdio.h>
#include <stdlib.h>

#include "includes/list.h"
#include "includes/sparse.h"

void print_sparse(Sparse* sparse)
{
    for(size_t i = 0; i < sparse->rows; ++i)
    {
        print_list_with_separator(sparse->row_lists[i], "\n");
    }
}

void add_node_to_sparse(Sparse* sparse, Node* node, size_t row)
{
    if(sparse->rows < row)
    {
        clear_node(node);
        return;
    }
    add_node_to_list(sparse->row_lists[row], node);
}

void clear_sparse(Sparse* sparse)
{
    if(sparse != NULL)
    {
        for(size_t i = 0; i < sparse->rows; ++i)
        {
            clear_list(sparse->row_lists[i]);
        }
    }
}

Sparse* create_sparse(size_t rows, size_t columns)
{
    if(rows == 0 || columns == 0)
    {
        fprintf(stderr, "Error: rows and columns must be positive integers\n");
        return NULL;
    }

    struct Sparse* sparse = (struct Sparse*)malloc(sizeof(struct Sparse));
    if(sparse == NULL)
    {
        fprintf(stderr, "Error: Memory allocation failed for sparse matrix\n");
        return NULL;
    }

    sparse->rows = rows;
    sparse->columns = columns;

    sparse->row_lists = (List**)calloc(rows, sizeof(List*));
    if(sparse->row_lists == NULL)
    {
        fprintf(stderr, "Error: Memory allocation failed for row lists\n");
        clear_sparse(sparse);
        return NULL;
    }

    for(size_t i = 0; i < rows; i++)
    {
        sparse->row_lists[i] = NULL;
    }

    return sparse;
}
