#include "includes/list.h"
#include <stdio.h>
#include <stdlib.h>

void print_node(Node* node) { printf("%.2lf", node->value); }

void print_node_with_separator(Node* node, char* separator)
{
    if(separator)
        printf("%.2lf %s", node->value, separator);
    else
        print_node(node);
}

Node* create_node(size_t column)
{
    Node* node = malloc(sizeof(Node));
    if(node == NULL)
    {
        fprintf(stderr, "Error: Memory allocation failed for node\n");
        return NULL;
    }
    node->column = column;
    node->next = NULL;

    return node;
}

Node* create_node_with_value(size_t column, double value)
{
    Node* node = create_node(column);
    node->value = value;

    return node;
}

void clear_node(Node* node) { free(node); }

void print_list(List* list)
{
    if(list)
    {
        Node* node = list->head;
        while(node)
            print_node_with_separator(node, " ");
    }
}

void print_list_with_separator(List* list, char* separator)
{
    print_list(list);
    if(separator)
        printf("%s", separator);
}

List* create_list(void)
{
    List* list = malloc(sizeof(List));
    if(list == NULL)
    {
        fprintf(stderr, "Error: Memory allocation failed for list\n");
        return NULL;
    }
    return list;
}

void clear_list(List* list)
{
    Node* node = list->head;
    Node* cur_node = NULL;
    while(node)
    {
        cur_node = node->next;
        clear_node(node);
        node = cur_node;
    }
    free(list);
}

void add_node_to_already_existed_node(Node* exsisted_node, Node* new_node)
{
    fprintf(stderr,
            "Error: trying to emplace node on already existing position\n");
    exsisted_node->value += new_node->value; // TBD
    clear_node(new_node);
}

void add_node_to_list(
    List* list,
    Node* node) // node may become invalid after usage in this function
{
    size_t column = node->column;
    Node* list_node = list->head;

    if(!list_node)
    {
        list->head = node;
        list->tail = node;
        return;
    }

    if(list_node->column > column)
    {
        list->head = node;
        node->next = list_node;
        return;
    }
    else if(list_node->column == column)
    {
        add_node_to_already_existed_node(list_node, node);
        return;
    }

    Node* after_list_node = list_node->next;

    while(after_list_node && after_list_node->column <= column)
    {
        if(after_list_node->column == column)
        {
            add_node_to_already_existed_node(after_list_node, node);
            return;
        }
        list_node = after_list_node;
        after_list_node = after_list_node->next;
    }

    // insert between list_node and after_list_node
    list_node->next = node;
    node->next = after_list_node;

    // end of list reached
    if(!after_list_node)
    {
        list->tail = node;
    }
}
