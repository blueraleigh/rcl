#include <R.h>
#include "node.h"
#include "sb.h"


struct node *node_alloc()
{
    struct node *x = Calloc(1, struct node);
    x->label = NULL;
    x->degree = 0;
    x->index = 0;
    x->brlen = 0.0;
    x->parent = NULL;
    x->children = NULL;
    x->data = NULL;
    x->data_free = NULL;
    return x;
}


void node_free(struct node *x)
{
    if (x == NULL)
        return;
    for (int i = 0; i < x->degree; ++i)
        node_free(*(x->children + i));
    Free(x->label);
    Free(x->children);
    if (x->data_free != NULL)
        x->data_free(x->data);
    Free(x);
}


void node_add_data(struct node *x, void *data, void (*data_free)(void*))
{
    if (x->data && x->data_free) {
        x->data_free(x->data);
        x->data = NULL;
        x->data_free = NULL;
    }
    x->data = data;
    x->data_free = data_free;
}


void node_add_child(struct node *parent, struct node *child)
{
    parent->degree += 1;
    struct node **children = Realloc(parent->children, parent->degree, struct node*);
    children[parent->degree - 1] = child;
    parent->children = children;
    child->parent = parent;
}


static void node_read_label(struct node *x, char **newick)
{
    char *buffer = NULL;

    int read = 1;

    while (read) {
        switch (*(*newick)) {
            case ':':
            case ',':
            case '(':
            case ')':
            case '[':
            case ']':
            case ';':
                read = 0;
                break;
            default:
                sb_push(buffer, (*(*newick)++));
                break;
        }
    }
    if (sb_count(buffer) > 0) {
        sb_push(buffer, '\0');
        char *label = Calloc(sb_count(buffer), char);
        memset(label, '\0', sizeof(char) * (sb_count(buffer)));
        strcpy(label, buffer);
        x->label = label;
        sb_free(buffer);
    }
}


static void node_read_brlen(struct node *x, char **newick)
{
    if (*(*newick) == ':') {
        ++(*newick);

        char *buffer = NULL;

        int read = 1;

        while (read) {
            switch (*(*newick)) {
                case ':':
                case ',':
                case '(':
                case ')':
                case '[':
                case ']':
                case ';':
                    read = 0;
                    break;
                default:
                    sb_push(buffer, (*(*newick)++));
                    break;
            }
        }
        sb_push(buffer, '\0');
        x->brlen = atof(buffer);
        sb_free(buffer);
    }
}


static struct node *node_read_node(char **newick)
{
    struct node *parent = node_alloc();

    if (*(*newick) == '(') {

        ++(*newick);

        struct node *child = node_read_node(newick);
        node_add_child(parent, child);

        if (*(*newick) == ',') {
            ++(*newick);
            struct node *child = node_read_node(newick);
            node_add_child(parent, child);

        } else {
            error("invalid tree");
        }

        while (*(*newick) == ',') {
            ++(*newick);
            struct node *child = node_read_node(newick);
            node_add_child(parent, child);
        }

        ++(*newick);

    }

    node_read_label(parent, newick);
    node_read_brlen(parent, newick);

    return parent;
}


struct node *node_read_newick(char **newick)
{
    return node_read_node(newick);
}


static void node_write_newick(struct node *x, char **buffer)
{
    char brlen[50], *buf = *buffer;
    memset(brlen, '\0', sizeof(char) * 50);
    if (x->degree == 0) {
        sprintf(brlen, "%f", x->brlen);
        sb_add(buf, strlen(brlen)+1);
        sb_add(buf, strlen(x->label)+1);
        sb_add(buf, 2);
        strcat(buf, x->label);
        strcat(buf, ":");
        strcat(buf, brlen);
    } else {
        sb_add(buf, 2);
        strcat(buf, "(");
        *buffer = buf;
        for (int i = 0; i < x->degree; ++i) {
            node_write_newick(*(x->children+i), buffer);
            buf = *buffer;
            sb_add(buf, 2);
            if (i < (x->degree - 1))
                strcat(buf, ",");
            *buffer = buf;
        }
        sb_add(buf, 3);
        strcat(buf, "):");
        sprintf(brlen, "%f", x->brlen);
        sb_add(buf, strlen(brlen)+1);
        strcat(buf, brlen);
    }
    *buffer = buf;
}


char *node_print_newick(struct node *x)
{
    char *buffer = NULL;
    sb_push(buffer, '\0');
    node_write_newick(x, &buffer);
    sb_add(buffer, 2);
    strcat(buffer, ";");
    char *newick = Calloc(sb_count(buffer)+1, char);
    strcpy(newick, buffer);
    sb_free(buffer);
    return newick;
}
