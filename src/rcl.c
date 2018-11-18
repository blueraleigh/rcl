#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>

#include "rcl.h"

#define CALLDEF(name, n)  {#name, (DL_FUNC) &name, n}


static const R_CallMethodDef CallEntries[] = {
    CALLDEF(rcl_read_newick, 1),
    CALLDEF(rcl_build_tree, 1),
    CALLDEF(rcl_write_newick, 1),
    CALLDEF(rcl_subtree, 2),
    CALLDEF(rcl_tiplabels, 1),
    CALLDEF(rcl_node_brlens, 1),
    CALLDEF(rcl_node_ages, 1),

    CALLDEF(rcl_plot_ctree, 3),
    CALLDEF(rcl_plot_ptree, 2),

    CALLDEF(rcl_ancestors, 2),
    CALLDEF(rcl_children, 2),
    CALLDEF(rcl_mrca, 3),
    CALLDEF(rcl_subgraph, 2),
    CALLDEF(rcl_subgraph_brlen, 2),
    CALLDEF(rcl_subgraph_newick, 4),

    CALLDEF(rcl_tree_traverse, 4),
    CALLDEF(rcl_tree_step, 1),
    CALLDEF(rcl_tree_reset, 1),
    CALLDEF(rcl_tree_jump, 2),

    CALLDEF(rcl_mkp_loglk, 3),

    CALLDEF(rcl_mkp_asr, 4),

    CALLDEF(rcl_mkp_cache_build, 3),
    CALLDEF(rcl_mkp_smap_clade, 6),
    CALLDEF(rcl_mkp_smap_branchset, 6),

    CALLDEF(rcl_fitch_pscore, 3),
    CALLDEF(rcl_fitch_mpr, 4),
    CALLDEF(rcl_fitch_history, 5),
    CALLDEF(rcl_fitch_count, 3),
    CALLDEF(rcl_sankoff_pscore, 4),
    CALLDEF(rcl_sankoff_mpr, 6),
    CALLDEF(rcl_sankoff_history, 7),
    CALLDEF(rcl_sankoff_count, 5),
    CALLDEF(rcl_sankoff_cost, 5),

    CALLDEF(rcl_dmm_mcmc, 10),

    CALLDEF(rcl_bmm_mcmc, 10),

    {NULL, NULL, 0}
};


void attribute_visible R_init_rcl(DllInfo *info)
{
    R_registerRoutines(info, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(info, FALSE);
    R_forceSymbols(info, TRUE);
}
