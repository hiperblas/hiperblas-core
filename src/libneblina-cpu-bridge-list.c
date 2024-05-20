#include "neblina_list.h"
#include <stdio.h>
#include <stdlib.h>

list_t *list_new()
{
    return (list_t *)NULL;
}

list_t *list_append(list_t *list, object_t o)
{
    list_t *newnode = (list_t *)malloc(sizeof(list_t));
    newnode->obj = o;
    newnode->next = NULL;

    if (list == NULL)
    {
        return newnode;
    }
    else
    {
        list_t *ptr = list;
        while (ptr->next != NULL)
            ptr = ptr->next;
        ptr->next = newnode;
        return list;
    }
}

int list_len(list_t *list)
{
    list_t *ptr = list;
    int ret = 0;
    while (ptr != NULL)
    {
        ret++;
        ptr = ptr->next;
    }
    return ret;
}

object_t list_get(list_t *list, int i)
{
    list_t *ptr = list;
    int ret = 1;
    while (ptr != NULL)
    {
        if (ret == i)
            return ptr->obj;
        ret++;
        ptr = ptr->next;
    }
    fprintf(stderr, "Runtime error: list get() index out of range\n");
    exit(1);
}
