#ifndef RCL_NODE_H_
#define RCL_NODE_H_


struct node;
struct node {
    char *label;
    int degree;
    int index;
    double brlen;
    struct node *parent;
    struct node **children;
    void *data;
    void (*data_free)(void*);   // the destructor for the data pointer
};


struct node *node_alloc();
void node_free(struct node*);
void node_add_data(struct node*, void*, void (*data_free)(void*));
void node_add_child(struct node*, struct node*);
struct node *node_read_newick(char**);
char *node_print_newick(struct node*);

#endif
