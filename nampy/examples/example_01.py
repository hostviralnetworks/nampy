# A script to illustrate nampy and run network propagation
import nampy
from nampy.networkio import networkio
from nampy.annotation import idmapping
from nampy.manipulation import manipulation
from nampy.monopartiteanalysis import prince

# Let's work with a version of HumanNet
# Lee, I., Blom, U. M., Wang, P. I., Shim, J. E., & Marcotte, E. M. (2011). 
# Prioritizing candidate disease genes by network-based boosting of genome-wide association data. 
# Genome research, 21(7), 1109–21. doi:10.1101/gr.118992.110
data_dir = nampy.__path__[0] + '/data/'

network_file = data_dir + "HumanNet_v1_join_networkonly.txt"
humannet = networkio.create_network_model_from_textfile('humannet', network_file, verbose = True)

# Add ids courtesy of bioservices
# note you may need to "easy_install bioservices" first
# also, you may want to add email = 'me@myemail.me'
# so ncbi can contact you if there are issues with the query
# Note there may some errors returned while querying but 
# each query is generally successful within the three tries.
humannet = idmapping.get_more_node_ids(humannet, node_id_type = "Entrez Gene (GeneID)", mapping_types = ['Entrez Gene (GeneID)', "UniProtKB ACC", 'UniProtKB ID'], verbose = True)

# Note we may miss a few this way, make sure we at least get all of the Entrez Gene IDs
# e.g. the database might not requrn the id's we query with
counter = 0
# Just have one nodetype
for the_node in humannet.nodetypes[0].nodes:
    if len(the_node.notes["Entrez Gene (GeneID)"]) == 0:
        the_node.notes["Entrez Gene (GeneID)"].append(the_node.id)
        counter +=1
print(counter)

# APMS host targets from
# Jäger, S., Cimermancic, P., Gulbahce, N., Johnson, J. R., McGovern, K. E., Clarke, S. C.,
#  … Krogan, N. J. (2012). Global landscape of HIV-human protein complexes. 
# Nature, 481(7381), 365–70. doi:10.1038/nature10719
hiv_apms_file = data_dir + "published_hiv_apms_factors.txt"
apms_source = networkio.create_source_dict_from_textfile(hiv_apms_file)

apms_source_dict = idmapping.get_more_source_dict_ids(apms_source, "UniProtKB ACC/ID", mapping_types = ["Entrez Gene (GeneID)", "UniProtKB ACC", "UniProtKB ID"])
counter = 0
for the_key in apms_source_dict.keys():
    if len(apms_source_dict[the_key]["UniProtKB ACC"]) == 0:
        apms_source_dict[the_key]["UniProtKB ACC"].append(the_key)
        counter +=1
print(counter)

humannet, id_matching_dict = manipulation.add_source(humannet, apms_source_dict, match_key_type = 'Entrez Gene (GeneID)')

# Just 100 permutations done here for demonstration purposes, to establish more reliable statistics you will need more
the_result, the_permutations = prince.prince(humannet, n_permutations = 100, verbose = True)

# The propagation scores are stored in the_result.  Here is how to save/load
# prince.save_prince_result(humannet, the_result, the_permutations, filename)
# the_result, the_permutations = prince.load_prince_result(humannet, filename)

from nampy.statistics import networkstatistics
ttp_dict, side_dict = networkstatistics.get_pvalue_from_scores(the_result, the_permutations, verbose = True)

otp_dict = networkstatistics.ttp_to_otp(ttp_dict, side_dict)

# This is the dict of p-values from the propagations...
otp_dict_corrected = networkstatistics.mtcorrect(otp_dict, method = '"BH"')
 
