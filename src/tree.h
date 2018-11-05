#ifndef RCL_TREE_H_
#define RCL_TREE_H_

#include "node.h"

#define PREORDER 0
#define POSTORDER 1
#define INTERNAL_NODES_ONLY 1
#define ALL_NODES 0

struct tree;
struct tree_traversal;

struct tree {
    /* number of terminal nodes */
    int ntip;

    /* number of nodes (internal and terminal) */
    int nnode;

    /*
    ** lastvisit_all[nodeix] returns the node index of the final node
    ** visited by the node with index `nodeix` when traversing the
    ** node's subtree in pre-order sequence.
    */
    int *lastvisit_all;


    /*
    ** lastvisit_int[nodeix - tree->ntip] returns the node index of the final internal node
    ** visited by the internal node with index `nodeix` when traversing the
    ** node's subtree in pre-order sequence.
    */
    int *lastvisit_int;

    /*
    ** seq_pos[nodeix] returns the index position in the
    ** seq_all array of the node with index `nodeix`.
    */
    int *seq_pos;

    struct node *root;

    /*
    ** Array of all nodes. node[nodeix] returns the node with index
    ** `nodeix`. Terminal nodes are number 1...ntip-1 and internal
    ** nodes are numbered ntip...nnode-1. Node indices follow the
    ** sequence of nodes as visited in a pre-order traversal. The
    ** root node will always have index ntip.
    */
    struct node **node;

    /*
    ** Array of all nodes in pre-order sequence traversal. Its reverse
    ** is a post-order traversal.
    */
    struct node **seq_all;

    /*
    ** Array of internal nodes in pre-order sequence traversal. Its reverse
    ** is a post-order traversal.
    */
    struct node **seq_int;
};


struct tree_traversal {
    int order;           /* Flag to indicate if the traversal is pre-order or post-order */
    int visit;           /* Flag to indicate if traversal visits all nodes or just internal nodes */
    struct tree *tree;
    struct node *begin;
    struct node *node;
    struct node *next;
};


void tree_free(struct tree*);
struct tree *tree_read(const char*);
int tree_isbinary(struct tree*);
struct tree *tree_build_tree(struct node*);
struct node *tree_mrca(struct node*, struct node*);
struct node *tree_step(struct tree_traversal*);
struct node *tree_jump(struct node*, struct tree_traversal*);
struct tree_traversal tree_traverse(int, int, struct node*, struct tree*);
void tree_reset(struct tree_traversal*);

#endif
