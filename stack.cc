/* $Id: stack.cc,v 1.1 2008/06/10 02:18:26 garrett Exp $ */

/* tiny integer stack manager MP (Michael Pique) */
/* note that "top" points to next free location, not last used location */

#include <stdlib.h>
#include <stdio.h>
#include "stack.h"

stack stack_create(int maxsize) 
{
    stack s;
    s = (integer_stack_t *) malloc( sizeof(integer_stack_t) );
    if (s == NULL) return s;
    s->base = (int *) calloc( maxsize, sizeof(int) );
    s->size = maxsize;
    s->top = 0;
    s->trace = NULL;
    return s;
}

int stack_pop(stack s)
{
    (s->top)--;
    if (s->top < 0) {
        fprintf(stderr, "Stack empty pop\n");

        if (s->trace != NULL) fprintf(s->trace, "Stack empty pop\n");
        s->top = 0;
        return 0;
        }
    if (s->trace != NULL) fprintf(s->trace, "Stack(%d/%d) pop %d\n",
        s->top, s->size, s->base[s->top]);
    return s->base[s->top];
}

void  stack_push(stack s, int i)
{
    if (s->trace != NULL) fprintf(s->trace, "Stack(%d/%d) push %d\n",
        s->top, s->size,i);
    if (s->top > s->size) {
        fprintf(stderr, "Stack full (%d) push\n", s->size);

        if (s->trace != NULL) fprintf(s->trace, "Stack full push\n");
        s->top = s->size - 1;
    }
    s->base[s->top] = i;
    s->top++;
    return;
}

int stack_depth(stack s)
{
    if (s->trace != NULL) fprintf(s->trace, "Stack depth %d/%d\n",
        s->top, s->size);
    return s->top;
}

int stack_size(stack s)
{
    return s->size;
}

void stack_trace(stack s, FILE *f)
{
    if (f==NULL && s->trace != NULL) fprintf(s->trace, "Stack trace off\n");
    s->trace = f;
    if (f != NULL) fprintf(s->trace, "Stack (%d/%d) trace on\n", 
        s->top, s->size);
    return;
}

int stack_test(void)
{
    stack s;
    s = stack_create(10);
    stack_push(s, 3);
    stack_push(s, 1);
    stack_push(s, 4);
    stack_push(s, 1);
    stack_push(s, 5);
    printf("Depth of stack = %d\n", stack_depth(s));
    printf("Pop stack gives %d\n", stack_pop(s));
    printf("Pop stack gives %d\n", stack_pop(s));
    printf("Pop stack gives %d\n", stack_pop(s));
    printf("Pop stack gives %d\n", stack_pop(s));
    printf("Pop stack gives %d\n", stack_pop(s));
    printf("Pop stack gives %d\n", stack_pop(s));
    return 0;
}
