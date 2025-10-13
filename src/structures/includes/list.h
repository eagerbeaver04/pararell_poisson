#ifndef LIST_HEADER_H
#define LIST_HEADER_H

#include <stdio.h>

struct Node;

typedef struct Node Node;

struct Node
{
    size_t column;
    double value;
    Node* next;
};

void print_node(Node* node);
void print_node_with_separator(Node* node, char* separator);
Node* create_node(size_t column);
Node* create_node_with_value(size_t columnn, double value);
void clear_node(Node* node);

struct List;

typedef struct List List;

struct List
{
    Node* head;
    Node* tail;
};

void print_list(List* list);
void print_list_with_separator(List* list, char* separator);
List* create_list(void);
void clear_list(List* list);

void add_node_to_list(List* list, Node* node);

#endif // LIST_HEADER_H
