# A script to illustrate working with recon_1 and reporter interactions
import nampy
from nampy.io import networkio
from nampy.multipartiteanalysis import reporternodes
from nampy.statistics import networkstatistics


# Let's look at the impact of IFNb on metabolism in
# human airway epithelial cells.  The sample data
# for IFNb stimulation and control has already been 
# background corrected by RMA and is a subset of the 
# datasets from:
# Shapira, S. D. et al. (2009). A physical and regulatory map of host-influenza interactions 
# reveals pathways in H1N1 infection. Cell, 139(7), 1255–67.
data_dir = nampy.__path__[0] + '/data/'
the_rma_values = networkio.read_table_file_to_dict(data_dir + 'GSE19392_IFNb18hr.txt', force_to_float = True)

# Prepare p-values
# for each of the probes
# We don't need to do BH correction
# up front since we use this as input to
# reporter nodes, which does its own
# normalization for the null distribution
probe_stats_dict = {}
for the_probe_id in the_rma_values.keys():
    list_2 = [log2(the_rma_values[the_probe_id]['GSM528691_IFNb18hr_1']), log2(the_rma_values[the_probe_id]['GSM528692_IFNb18hr_2'])]
    list_1 = [log2(the_rma_values[the_probe_id]['GSM528749_media18hr_1']), log2(the_rma_values[the_probe_id]['GSM528750_media18hr_2'])]
    the_result = networkstatistics.t_uneqvar(list_1, list_2)
    probe_stats_dict[the_probe_id] = the_result

# Now aggregate the probe significance
# at the level of the Entrez gene ID.
# This is an abbreviated form of HG_U133A
# annotation file with probe : entrez gene mappings
# downloaded from Affymetrix (http://www.affymetrix.com)
# on Feb 4, 2014
data_dir = nampy.__path__[0] + '/data/'
the_probe_key_gene_value_dict = networkio.read_table_file_to_dict(data_dir + 'HT_HG_U133A_na34_annot_abbrev_140204.txt', force_to_float = False, comment_char = '#')
for the_probe_id in the_probe_key_gene_value_dict.keys():
    the_line = the_probe_key_gene_value_dict[the_probe_id]['Entrez Gene']
    the_line = the_line.split(' /// ')
    the_probe_key_gene_value_dict[the_probe_id] = the_line

# Combine p-values for 1 gene : many probes using Stouffer's method
aggregated_pvalue_dict = networkstatistics.aggregate_probe_statistics_to_gene(probe_stats_dict, the_probe_key_gene_value_dict, aggregation_type = 'stouffer')

# Let's work with recon 1 from:
# Duarte, N. C., Becker, S. a, Jamshidi, N., Thiele, I., Mo, M. L., Vo, T. D., … Palsson, B. Ø. (2007). 
# Global reconstruction of the human metabolic network based on genomic and bibliomic data. 
# PNAS, 104(6), 1777–82. doi:10.1073/pnas.0610772104
# Note if we want to be more sophisticated we
# could build a cell-line specific model to work with.
the_network = networkio.load_pickled_network(data_dir + 'recon_1_nampy_v01')
the_network.convert_to_multipartite()

# We don't have the RefSeq  mRNA mappings. 
# Assign p-values to each gene-level entrez id, 
# which can then be mapped back to transcripts.
the_entrez_ids = list(set([x.id.split('.')[0] for x in the_network.nodetypes.get_by_id("gene").nodes]))

# Pick out the model genes that successfully mapped
the_model_gene_pval_dict = {}
for the_gene_id in the_entrez_ids:
    # This indicates spontaneous reactions in some COBRA
    # models, we can ignore it if we like
    if the_gene_id != 's0001':
        if the_gene_id in aggregated_pvalue_dict['mapped'].keys():
            the_model_gene_pval_dict[the_gene_id] = aggregated_pvalue_dict['mapped'][the_gene_id]

# Check how many mapped
len(the_model_gene_pval_dict.keys())
len(the_entrez_ids)

# Now we need to map back to the model transcripts
the_transcript_pval_dict = {}
for the_transcript in the_network.nodetypes.get_by_id("gene").nodes:
    the_gene_id = the_transcript.id.split(".")[0]
    if the_gene_id in the_model_gene_pval_dict.keys():
        the_transcript_pval_dict[the_transcript.id] = the_model_gene_pval_dict[the_gene_id]
        
# Now we map to model reactions and then run reporter metabolites
hyperedge_score_dict = reporternodes.evaluate_reaction_pvalues(the_network, the_transcript_pval_dict)
the_reporter_dict = reporternodes.calculate_reporter_scores(the_network, hyperedge_score_dict['p'], 'metabolite', 'reaction', number_of_randomizations = 10000)


# We can prepare summar calculations to
# help with network visualization. 
# These files can be imported to Cytoscape
# to visualize the results.
from math import log10
node_property_dict = {}

# Uncorrected p's for the nodes
from copy import deepcopy
node_property_dict['uncorrected_p'] = deepcopy(hyperedge_score_dict['p'])
node_property_dict['dir'] = deepcopy(hyperedge_score_dict['dir'])
transcript_uncorrected_p = {}
for the_transcript_id in the_transcript_pval_dict.keys():
    node_property_dict['uncorrected_p'][the_transcript_id] = deepcopy(the_transcript_pval_dict[the_transcript_id]['p'])
    node_property_dict['dir'][the_transcript_id] = deepcopy(the_transcript_pval_dict[the_transcript_id]['dir'])
    transcript_uncorrected_p[the_transcript_id] = deepcopy(the_transcript_pval_dict[the_transcript_id]['p'])

# Corrected p's.  Note reporters are already corrected by their null distribution.
# We will split transcript and reaction p's for bh adjustment since we
# have a distinct fdr correction for each nodetype and for the edges
# We also will want to make a judicious selection of cutoff p since
# we only have 2 replicates to establish statistics.
from nampy.statistics import networkstatistics
node_property_dict['corrected_p_reporter'] = deepcopy(the_reporter_dict['p_values'])
node_property_dict['corrected_p_gene'] = (networkstatistics.mtcorrect(transcript_uncorrected_p, method = '"BH"'))
node_property_dict['corrected_p_reaction'] = (networkstatistics.mtcorrect(hyperedge_score_dict['p'], method = '"BH"'))
the_reactions = node_property_dict['corrected_p_reaction'].keys()
the_metabolites = node_property_dict['corrected_p_reporter'].keys()

# -log10 p's are nice to help with visualization
node_property_dict['dir_x_-log10(corrected_p_reaction)'] = {}
node_property_dict['dir_x_-log10(corrected_p_gene)'] = {}
nodes_to_calculate = set([x for x in (node_property_dict['corrected_p_gene'].keys() + node_property_dict['corrected_p_reaction'].keys())]).intersection(set([x for x in node_property_dict['dir'].keys()]))
for the_node_id in nodes_to_calculate:
    the_dir = node_property_dict['dir'][the_node_id]
    if the_dir == '-':
        the_dir = -1.
    else:
        the_dir = 1.
    if the_node_id in node_property_dict['corrected_p_gene']:
        the_p = node_property_dict['corrected_p_gene'][the_node_id]
        node_property_dict['dir_x_-log10(corrected_p_gene)'][the_node_id] = the_dir * -1. * log10(the_p)
    if the_node_id in node_property_dict['corrected_p_reaction']:
        the_p = node_property_dict['corrected_p_reaction'][the_node_id]
        node_property_dict['dir_x_-log10(corrected_p_reaction)'][the_node_id] = the_dir * -1. * log10(the_p)

node_property_dict['mlog10(corrected_p_reporter)'] = {}
node_property_dict['mlog10(corrected_p_reaction)'] = {}
node_property_dict['mlog10(corrected_p_gene)'] = {}
nodes_to_calculate = set(node_property_dict['corrected_p_reporter'].keys() + node_property_dict['corrected_p_gene'].keys() + node_property_dict['corrected_p_reaction'].keys())

for the_node_id in nodes_to_calculate:
    if the_node_id in node_property_dict['corrected_p_gene']:
        the_p = node_property_dict['corrected_p_gene'][the_node_id]
        node_property_dict['mlog10(corrected_p_gene)'][the_node_id] = -1. * log10(the_p)
    if the_node_id in node_property_dict['corrected_p_reaction']:
        the_p = node_property_dict['corrected_p_reaction'][the_node_id]
        node_property_dict['mlog10(corrected_p_reaction)'][the_node_id] = -1. * log10(the_p)
    if the_node_id in node_property_dict['corrected_p_reporter']:
        the_p = node_property_dict['corrected_p_reporter'][the_node_id]
        node_property_dict['mlog10(corrected_p_reporter)'][the_node_id] = -1. * log10(the_p)
    

# Also add the -log10(p) * dir to the edge
# property dict so these can assist 
# visualization in Cytoscape
edge_property_dict = {}
edge_property_dict['dir_x_-log10(corrected_p_edge)'] = {}
for the_edge in the_network.edges:
    the_node_list = the_edge.get_node_pair()
    # Only use edges that connect metabolites to reactions
    if 'reaction' in [x.get_nodetype() for x in the_node_list]:
        if 'metabolite' in [x.get_nodetype() for x in the_node_list]:
            the_reaction_id = [x.id for x in the_node_list if x.get_nodetype() == 'reaction'][0]
            if the_reaction_id in node_property_dict['corrected_p_reaction'].keys():
                edge_property_dict['dir_x_-log10(corrected_p_edge)'][the_edge.id] = (node_property_dict['dir_x_-log10(corrected_p_reaction)'][the_reaction_id])

# Don't write the genes to file as nodes,
# this isn't needed for visualization
from nampy.io import cytoscapeio
cytoscapeio.write_network_textfile(the_network, properties_dict = edge_property_dict, exclude_nodetypes = ['gene'])
cytoscapeio.write_node_attributes_to_textfile(the_network, properties_dict = node_property_dict, exclude_nodetypes = ['gene'])
# A .cys file made from this analysis is
# available for viewing.

# Optionally, check the histograms for 
# the reporters to pick out cutoffs
# You'll need matplotlib for this.
hist(node_property_dict['corrected_p_reporter'].values(),100)
# A ctuoff around 0.05 appears OK
hist(node_property_dict['corrected_p_reaction'].values(),100)
# We won't cut out reactions based on p but we will visualize up/down in cytoscape
# Note N is only 2 replicates for the examples here.
