#ifndef RCL_H_
#define RCL_H_

#include <R.h>
#include <Rinternals.h>


/* Declarations for all native routine entry points */

/* treeio.c */
SEXP rcl_read_newick(SEXP);
SEXP rcl_build_tree(SEXP);
SEXP rcl_write_newick(SEXP);
SEXP rcl_subtree(SEXP, SEXP);
SEXP rcl_tiplabels(SEXP);
SEXP rcl_node_brlens(SEXP);
SEXP rcl_node_ages(SEXP);


/* treeplot.c */
SEXP rcl_plot_ctree(SEXP, SEXP, SEXP);
SEXP rcl_plot_ptree(SEXP, SEXP);


/* treeutil.c */
SEXP rcl_ancestors(SEXP, SEXP);
SEXP rcl_children(SEXP, SEXP);
SEXP rcl_mrca(SEXP, SEXP, SEXP);
SEXP rcl_subgraph(SEXP, SEXP);
SEXP rcl_subgraph_brlen(SEXP, SEXP);
SEXP rcl_subgraph_newick(SEXP, SEXP, SEXP, SEXP);
SEXP rcl_tree_traverse(SEXP, SEXP, SEXP, SEXP);
SEXP rcl_tree_step(SEXP);
SEXP rcl_tree_reset(SEXP);
SEXP rcl_tree_jump(SEXP, SEXP);


/* mkp.c */
SEXP rcl_mkp_loglk(SEXP, SEXP, SEXP);

/* mkp_asr.c */
SEXP rcl_mkp_asr(SEXP, SEXP, SEXP, SEXP);


/* mkp_smap.c */
SEXP rcl_mkp_cache_build(SEXP, SEXP, SEXP);
SEXP rcl_mkp_smap_clade(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP rcl_mkp_smap_branchset(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);


/* parsimony.c */
SEXP rcl_fitch_pscore(SEXP, SEXP, SEXP);
SEXP rcl_fitch_mpr(SEXP, SEXP, SEXP, SEXP);
SEXP rcl_fitch_history(SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP rcl_fitch_count(SEXP, SEXP, SEXP);
SEXP rcl_sankoff_pscore(SEXP, SEXP, SEXP, SEXP);
SEXP rcl_sankoff_mpr(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP rcl_sankoff_history(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP rcl_sankoff_count(SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP rcl_sankoff_cost(SEXP, SEXP, SEXP, SEXP, SEXP);


/* dmm_mcmc.c */
SEXP rcl_dmm_mcmc(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);


#endif
