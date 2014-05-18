# A script to illustrate nampy and run network propagation
# We follow the data and analysis published in:
# Emig-Agius, D., Olivieri, K., Pache, L., Shih, C., Pustovalova, O., 
# Bessarabova, M., … Ideker, T. (2014). An integrated map of HIV-human 
# protein complexes that facilitate viral infection. PLOS ONE, 9(5), e96687.

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
humannet = idmapping.get_more_node_ids(humannet, node_id_type = "Entrez Gene (GeneID)", mapping_types = ['Entrez Gene (GeneID)', "UniProtKB ACC", 'UniProtKB ID', 'Symbol'], verbose = True)

# APMS host targets from
# Jäger, S., Cimermancic, P., Gulbahce, N., Johnson, J. R., McGovern, K. E., Clarke, S. C.,
#  … Krogan, N. J. (2012). Global landscape of HIV-human protein complexes. 
# Nature, 481(7381), 365–70. doi:10.1038/nature10719
hiv_apms_file = data_dir + "published_hiv_apms_factors.txt"
apms_source = networkio.create_source_dict_from_textfile(hiv_apms_file)

apms_source_dict = idmapping.get_more_source_dict_ids(apms_source, "UniProtKB ACC/ID", mapping_types = ["Entrez Gene (GeneID)", "UniProtKB ACC", "UniProtKB ID"])
# Make sure note to lose information in the UniProt query.
# get_more_source_dict_ids didn't check this since
# we don't specify "UniProtKB ACC" or "UniProtKB ID"
# for UniProt queries
counter = 0
for the_key in apms_source_dict.keys():
    if len(apms_source_dict[the_key]["UniProtKB ACC"]) == 0:
        apms_source_dict[the_key]["UniProtKB ACC"].append(the_key)
        counter +=1

# This step automatically paired 2,597 when it was run while testing this script.
# The number may change based on ID database updates.
humannet, id_matching_dict = manipulation.add_source(humannet, apms_source_dict, match_key_type = 'Entrez Gene (GeneID)')

# Note that you can check the status of ID's that failed a 1:1 pairing in id_matching_dict['unmatched_ids'].  For example,
# you can look at id_matching_dict['unmatched_ids']['network_one_to_source_dict_many'] and evaluate how you
# would like these multiple protein products from a single gene to influence the propagation on the gene-based
# network.  We will neglect these for this example.

# Just 100 permutations done here for demonstration purposes.
# To establish more reliable statistics you will need more permutations.
# FYI, on an iMac Corei7/8GB test system it took 12-18 hr to run 10,000 permutations.
the_result, the_permutations = prince.prince(humannet, n_permutations = 100, verbose = True)

# The propagation scores are stored in the_result.  Here is how to save/load
# filename = 'test_apms_propagations'
# prince.save_prince_result(humannet, the_result, the_permutations, filename)
# networkio.pickle_network(humannet, 'humannet')
# humannet = networkio.load_pickled_network('humannet')
# the_result, the_permutations = prince.load_prince_result(humannet, filename)

from nampy.statistics import networkstatistics
ttp_dict_apms, side_dict_apms = networkstatistics.get_pvalue_from_scores(the_result, the_permutations, verbose = True)

otp_dict_apms = networkstatistics.ttp_to_otp(ttp_dict_apms, side_dict_apms)

# This is the dict of p-values from the propagations...
otp_dict_corrected_apms = networkstatistics.mtcorrect(otp_dict_apms, method = 'BH')
 
# We are done with the APMS propagation.  Clear the permutations, which can utilize a lot of memory.
del the_permutations


###
# Now, repeat the analysis with the RNAi data
hiv_rnai_file = data_dir + 'published_hiv_host_restriction_factors.txt'
rnai_source = networkio.create_source_dict_from_textfile(hiv_rnai_file)

# Query for more ID's, though this isn't strictly needed since Entrez ID's
# were used for the nodes in HumanNet
rnai_source_dict = idmapping.get_more_source_dict_ids(rnai_source, "Entrez Gene (GeneID)", mapping_types = ["Entrez Gene (GeneID)", "UniProtKB ACC", "UniProtKB ID"])
# We can check the recovery of the Entrez ID's from the query.
counter = 0
for the_key in rnai_source_dict.keys():
    if len(rnai_source_dict[the_key]["Entrez Gene (GeneID)"]) == 0:
        rnai_source_dict[the_key]["Entrez Gene (GeneID)"].append(the_key)
        counter +=1

# This step automatically paired 676 when it was run while testing this script.
humannet, id_matching_dict = manipulation.add_source(humannet, rnai_source_dict, match_key_type = 'Entrez Gene (GeneID)', reset_source = True)

the_result, the_permutations = prince.prince(humannet, n_permutations = 100, verbose = True)

# filename = 'test_rnai_propagations'
# prince.save_prince_result(humannet, the_result, the_permutations, filename)
# humannet = networkio.load_pickled_network('humannet')
# the_result, the_permutations = prince.load_prince_result(humannet, filename)

from nampy.statistics import networkstatistics
ttp_dict_rnai, side_dict_rnai = networkstatistics.get_pvalue_from_scores(the_result, the_permutations, verbose = True)

otp_dict_rnai = networkstatistics.ttp_to_otp(ttp_dict_rnai, side_dict_rnai)

# This is the dict of p-values from the propagations...
otp_dict_corrected_rnai = networkstatistics.mtcorrect(otp_dict_rnai, method = 'BH')
 
# We are done with the APMS propagation.  Clear the permutations, which can utilize a lot of memory.
del the_permutations

##
# Now take a look at the size of our high confidence set.
# If you increase the number of permutations for the APMS and RNAi
# datasets, you can compare this to the set of 554 that was found
# in Emig-Agius 2014

p_dict_sig_apms = {x: otp_dict_corrected_apms[x] for x in otp_dict_corrected_apms.keys() if otp_dict_corrected_apms[x] < 0.0001}
len(p_dict_sig_apms.keys())

p_dict_sig_rnai = {x: otp_dict_corrected_rnai[x] for x in otp_dict_corrected_rnai.keys() if otp_dict_corrected_rnai[x] < 0.0001}
len(p_dict_sig_rnai.keys())

high_confidence_set = set(p_dict_sig_apms.keys()).intersection(set(p_dict_sig_rnai.keys()))
len(high_confidence_set)

# Another example that integrates network results with CORUM complexes will be presented
# in example_03
