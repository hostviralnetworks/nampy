# A script to illustrate working with recon_1 and reporter interactions
import nampy
from nampy.io import networkio
from nampy.cobra import reporternodes

# Let's work with recon 1 from:
# Duarte, N. C., Becker, S. a, Jamshidi, N., Thiele, I., Mo, M. L., Vo, T. D., … Palsson, B. Ø. (2007). 
# Global reconstruction of the human metabolic network based on genomic and bibliomic data. 
# PNAS, 104(6), 1777–82. doi:10.1073/pnas.0610772104
data_dir = nampy.__path__[0] + '/data/'
the_network = networkio.load_pickled_network(data_dir + 'recon_1_nampy_v01')
the_network.convert_to_multipartite()

# Say we don't have transcript-level data. Besides, we don't have 100% of the RefSeq
# mRNA mappings.  Let's assign p-values to each gene-level entrez id, which can then be mapped
# back to transcripts.
the_entrez_ids = list(set([x.id.split('.')[0] for x in the_network.nodetypes.get_by_id("gene").nodes]))

# In the absence of real data, we can just randomly generate
# the p-value to illustrate the method
the_gene_pval_dict = {}
for the_gene_id in the_entrez_ids:
    # This indicates spontaneous reactions in some COBRA
    # models, we can ignore it if we like
    if the_gene_id != 's0001':
        the_gene_pval_dict[the_gene_id] = rand()

# Now we need to map back to the model transcripts
# This is valid since the different isoforms are largely
# ascribed the same function
the_transcript_pval_dict = {}
for the_transcript in the_network.nodetypes.get_by_id("gene").nodes:
    the_gene_id = the_transcript.id.split(".")[0]
    the_transcript_pval_dict[the_transcript.id] = the_gene_pval_dict[the_gene_id]

hyperedge_score_dict = reporternodes.evaluate_reaction_pvalues(the_network, the_transcript_pval_dict)

the_reporter_dict = reporternodes.calculate_reporter_scores(the_network, hyperedge_score_dict, 'metabolite', 'reaction', number_of_randomizations = 10000)

# Now write these to file.  These files can be imported to Cytoscape
# to visualize the results.
from nampy.io import cytoscapeio
cytoscapeio.write_network_textfile(the_network)
node_property_dict = {}
node_property_dict['p_values'] = the_reporter_dict['p_values']
node_property_dict['p_values'].update(hyperedge_score_dict)
node_property_dict['p_values'].update(the_transcript_pval_dict)
cytoscapeio.write_node_attributes_to_textfile(the_network, properties_dict = node_property_dict)
