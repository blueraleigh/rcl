#include <R.h>
#include "tree.h"
#include "sb.h"


static struct tree *tree_alloc()
{
    struct tree *tree = Calloc(1, struct tree);
    tree->ntip = 0;
    tree->nnode = 0;
    tree->lastvisit_all = NULL;
    tree->lastvisit_int = NULL;
    tree->seq_pos = NULL;
    tree->root = NULL;
    tree->node = NULL;
    tree->seq_all = NULL;
    tree->seq_int = NULL;
    return tree;
}


void tree_free(struct tree *tree)
{
    if (tree == NULL)
        return;
    node_free(tree->root);
    Free(tree->lastvisit_all);
    Free(tree->lastvisit_int);
    Free(tree->seq_pos);
    Free(tree->node);
    Free(tree->seq_all);
    Free(tree->seq_int);
    Free(tree);
}


// tips are numbered 0 -> ntip-1; internal nodes, ntip -> nnode-1
// consecutively by pre-order traversal order.
static void tree_index_nodes(struct node *node, struct tree *tree, int *ix1, int *ix2)
{
    if (node->degree > 0)
        node->index = tree->ntip + (*ix2)++;
    else
        node->index = (*ix1)++;
    for (int i = 0; i < node->degree; ++i)
        tree_index_nodes(node->children[i], tree, ix1, ix2);
}


static void tree_preorder_nodes(struct node *node, struct tree *tree, int *ix1, int *ix2)
{
    tree->seq_pos[node->index] = *ix1;
    tree->node[node->index] = node;
    if (node->degree > 0)
        tree->seq_int[(*ix2)++] = node;
    tree->seq_all[(*ix1)++] = node;
    for (int i = 0; i < node->degree; ++i)
        tree_preorder_nodes(node->children[i], tree, ix1, ix2);
    tree->lastvisit_all[node->index] = (tree->seq_all[*ix1-1])->index;
    if (node->degree > 0)
        tree->lastvisit_int[node->index - tree->ntip] = (tree->seq_int[*ix2-1])->index;
}


static void tree_count_nodes(struct node *node, int *nnode, int *ntip)
{
    ++(*nnode);
    if (node->degree == 0)
        ++(*ntip);
    for (int i = 0; i < node->degree; ++i)
        tree_count_nodes(node->children[i], nnode, ntip);
}


static struct node *tree_next_node(struct tree_traversal *t)
{
    struct tree *tree = t->tree;
    struct node *begin = t->begin;
    struct node *node = t->node;

    int seq_pos, stop;

    if (begin->index < tree->ntip)  // this is a terminal node
        return NULL;

    if (node == NULL)
        return NULL;

    if (t->order == PREORDER) {
        if (t->visit == INTERNAL_NODES_ONLY) {
            seq_pos = node->index - tree->ntip;
            stop = tree->lastvisit_int[begin->index - tree->ntip];

            if (node->index == stop)
                return NULL;

            return tree->seq_int[++seq_pos];

        } else {
            seq_pos = tree->seq_pos[node->index];
            stop = tree->lastvisit_all[begin->index];

            if (node->index == stop)
                return NULL;

            return tree->seq_all[++seq_pos];
        }

    } else if (t->order == POSTORDER) {
        stop = begin->index;

        if (node->index == stop)
            return NULL;

        if (t->visit == INTERNAL_NODES_ONLY) {
            seq_pos = node->index - tree->ntip;
            return tree->seq_int[--seq_pos];

        } else {
            seq_pos = tree->seq_pos[node->index];
            return tree->seq_all[--seq_pos];
        }
    }

    return NULL;
}


struct tree *tree_build_tree(struct node *root)
{
    struct tree *tree = tree_alloc();
    tree->root = root;
    tree_count_nodes(tree->root, &tree->nnode, &tree->ntip);
    tree->lastvisit_all = Calloc(tree->nnode, int);
    tree->lastvisit_int = Calloc(tree->nnode - tree->ntip, int);
    tree->seq_pos = Calloc(tree->nnode, int);
    tree->node = Calloc(tree->nnode, struct node*);
    tree->seq_all = Calloc(tree->nnode, struct node*);
    tree->seq_int = Calloc(tree->nnode - tree->ntip, struct node*);
    int j = 0, k = 0;
    tree_index_nodes(tree->root, tree, &j, &k);
    j = 0, k = 0;
    tree_preorder_nodes(tree->root, tree, &j, &k);
    return tree;
}


struct tree *tree_read(const char *newick)
{
    char *newick2read = (char *)newick;
    struct tree *tree = tree_build_tree(node_read_newick(&newick2read));
    return tree;
}


int tree_isbinary(struct tree *tree)
{
    return (tree->ntip == ((tree->nnode + 1) / 2)) ? 1 : 0;
}


struct node *tree_mrca(struct node *a, struct node *b)
{
    if (a == NULL || b == NULL)
        return NULL;
    if (a->index == b->index)
        return a;
    struct node **path = NULL, *mrca = NULL, *p = NULL;
    sb_push(path, a);
    p = a->parent;
    while (p != NULL) {
        sb_push(path, p);
        p = p->parent;
    }
    p = b;
    while (p != NULL) {
        for (int i = 0; i < sb_count(path); ++i) {
            if (path[i]->index == p->index) {
                mrca = path[i];
                sb_free(path);
                return mrca;
            }
        }
        p = p->parent;
    }
    sb_free(path);
    return mrca;
}


/*
** Initialize the structure for performing pre-order and post-order tree traversal.
**
** Given an initial node, this function returns a structure that can be used to
** perform pre-order and post-order traversal of the node's subtree in combination with
** the tree_step function that follows this function. The canonical way to perform
** traversals is like so,
**
**     struct tree_traversal;
**     struct node *node = tree->root;
**     t = tree_traverse(PREORDER, INTERNAL_NODES_ONLY, node, tree);
**     // t = tree_traverse(PREORDER, ALL_NODES, node, tree);
**     // t = tree_traverse(POSTORDER, INTERNAL_NODES_ONLY, node, tree);
**     // t = tree_traverse(POSTORDER, ALL_NODES, node, tree);
**     node = tree_step(&t);
**     while (node != NULL) {
**          // perform some operations
**          node = tree_step(&t);
**     }
*/
struct tree_traversal tree_traverse(int order, int visit, struct node *node, struct tree *tree)
{
    struct tree_traversal t;
    t.order = order;
    if (visit == INTERNAL_NODES_ONLY || visit == ALL_NODES)
        t.visit = visit;
    else
        error("Unrecognized visitation type");
    t.tree = tree;
    if (order == PREORDER) {
        t.begin = node;
        t.node = node;
        t.next = node;
    } else if (order == POSTORDER) {
        if (visit == ALL_NODES) {
            t.begin = node;
            t.node = tree->node[tree->lastvisit_all[node->index]];
            t.next = t.node;
        } else {
            t.begin = node;
            t.node = tree->node[tree->lastvisit_int[node->index - tree->ntip]];
            t.next = t.node;
        }
    } else {
        error("Unrecognized tree traversal sequence");
    }
    return t;
}


struct node *tree_step(struct tree_traversal *t)
{
    t->node = t->next;
    t->next = tree_next_node(t);
    return t->node;
}


/*
** Skip traversal of *node's subtree and go straight to the next
** node that would be visited after completing a full traversal of
** *node's subtree.
*/
struct node *tree_jump(struct node *node, struct tree_traversal *t)
{
    int stop;
    int node_stop;
    int seq_pos;

    if (t->order == PREORDER) {
        if (t->visit == INTERNAL_NODES_ONLY) {
            stop = t->tree->lastvisit_int[t->begin->index - t->tree->ntip];
            node_stop = t->tree->lastvisit_int[node->index - t->tree->ntip];
            t->node = t->tree->node[node_stop];
            if (node_stop == stop) {
                t->next = NULL;
            } else {
                seq_pos = node->index - t->tree->ntip;
                t->next = t->tree->seq_int[++seq_pos];
            }
        } else {
            stop = t->tree->lastvisit_all[t->begin->index];
            node_stop = t->tree->lastvisit_all[node->index];
            t->node = t->tree->node[node_stop];
            if (node_stop == stop) {
                t->next = NULL;
            } else {
                seq_pos = t->tree->seq_pos[node_stop];
                t->next = t->tree->seq_all[++seq_pos];
            }
        }
    } else {
        error("Cannot perform postorder traversal jump.");
    }

    return tree_step(t);
}


void tree_reset(struct tree_traversal *t)
{
    if (t->order == PREORDER) {
        t->node = t->begin;
        t->next = t->node;
    } else if (t->order == POSTORDER) {
        if (t->visit == ALL_NODES) {
            t->node = t->tree->node[t->tree->lastvisit_all[t->begin->index]];
            t->next = t->node;
        } else {
            t->node = t->tree->node[t->tree->lastvisit_int[t->begin->index - t->tree->ntip]];
            t->next = t->node;
        }
    }
}
