
#--------------------------------------------------------------------------------------
#
#    Functions for the eQTL sharing pipeline
#
#--------------------------------------------------------------------------------------
#
# Functions: get_parent_cell_type()     L1 parent of a subcluster (GABA-1 -> GABA)
#
# Packages:  base R only.
#
# Used in:   devspec_extract_universe.R   get_parent_cell_type
#            devspec_eqtl_effect_sizes.R  get_parent_cell_type
#
# Notes:     Functions shared with the single-cell replication pipeline (Jang
#            artifact filter, allele harmonisation, peak enrichment) are in
#            repl_sc_functions.R.
#
#--------------------------------------------------------------------------------------
 
# L1 parent of a cell type: L1 types return themselves; subclusters drop their
# last "-<suffix>" (GABA-1 -> GABA, Glu-UL-0 -> Glu-UL)
get_parent_cell_type <- function(cell_type, cell_types_L1) {
  if (cell_type %in% cell_types_L1) return(cell_type)
  sub("-[^-]+$", "", cell_type)
}
 
#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
 
