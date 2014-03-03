from ..core.shared_functions import test_kwarg
from ..core.parameters import available_mapping_target, available_mapping_source


default_mapping_target_list = ["UniProtKB ACC", 
    "UniProtKB ID",
    "EMBL/GenBank/DDBJ",
    "EMBL/GenBank/DDBJ CDS",
    "UniGene",
    "Entrez Gene (GeneID)",
    "GI number*", 
    "RefSeq Protein",
    "RefSeq Nucleotide",
    "PDB",
    "Ensembl",
    "Ensembl Protein",
    "Ensembl Transcript",
    "Ensembl Genomes",
    "Ensembl Genomes Protein",
    "Ensembl Genomes Transcript",
    "KEGG",
    "PATRIC",
    "GeneCards",
    "BioCyc",
    "DrugBank",
    "GenomeRNAi",
    "NextBio"]


def get_more_node_ids(the_network, **kwargs):
    """ Script to add more identifiers to model notes
    based on the node.id

    Arguments:
     the_network: a Network object, modified in place

    kwargs:
     node_id_type: current type of ids used for the nodes.
      Currently can be Entrez Gene (GeneID) or any of the 
      options in the BioServices UniProt mappings
     mapping_types: a list of target mapping id types to include
     verbose:

    Returns:
     the_network 

    TODO: determine the best source db/module
    for pairings from bioservices
                  
    """
    continue_flag = True

    try:
        from bioservices import UniProt
        u = UniProt(verbose=False)
    except:
        print("No bioservices module installed or cannot connect, exiting...")
        print("e.g. if you are using pip, did you 'pip install bioservices'?")
        continue_flag = False

    the_node_locations = the_network.get_node_locations()
    if len(the_node_locations) == 0:
        print 'The network has no nodes, exiting...'
        continue_flag = False

    if 'node_id_type' in kwargs:
        node_id_type = kwargs['node_id_type'] 
    else:
        node_id_type = "Entrez Gene (GeneID)"

    if 'mapping_types' in kwargs:
        mapping_types = kwargs['mapping_types'] 
    else:
        mapping_types = default_mapping_target_list

    verbose = test_kwarg('verbose', kwargs, [True, False])

    # Maximum number of items to
    # query at a time 
    # Note there is a length limit in bioservices 1.2.1
    # for the web-based query string.
    # Trial-and-error suggests the most
    # id's that can be queried are
    # between 100 and 1000
    max_query_length = 500

    if continue_flag:
        query_string = ''
        model_node_ids = []
        for the_nodetype in the_network.nodetypes:
            model_node_ids += [x.id for x in the_nodetype.nodes]

        the_node_id_list_list = [[]]
        i = 0
        j = 0
        for the_node_id in model_node_ids:
            if (j + 1) % max_query_length == 0:
                the_node_id_list_list.append([])
                i += 1
                the_node_id_list_list[i] = []
                j = 0
            the_node_id_list_list[i].append(the_node_id)
            j += 1

        query_string_list = []
        for i, the_node_id_list in enumerate(the_node_id_list_list):
            query_string = ''
            for the_node_id in the_node_id_list:
                if len(query_string) > 0:
                    query_string = query_string + ' ' + the_node_id
                else:
                    query_string = the_node_id
            query_string_list.append(query_string)

        
        for the_target_type in mapping_types:
            the_result = {}
            for the_query_string in query_string_list:
                the_result.update(u.mapping(fr = available_mapping_source[node_id_type], to = available_mapping_target[the_target_type], query = the_query_string))
            if verbose:
                print("**Finished mapping for %s to %s.**" % (node_id_type, the_target_type))
            for the_nodetype in the_network.nodetypes:
                for the_node in the_nodetype.nodes:
                    if (the_node.id in the_result.keys()):
                        if len(the_result[the_node.id]) > 0:
                            the_node.notes[the_target_type] = the_result[the_node.id]
                        else:
                            the_node.notes[the_target_type] = []
                    else:
                        the_node.notes[the_target_type] = []

    return the_network



def get_more_source_dict_ids(source_dict, primary_key, **kwargs):
    """ Script to add more ids to source dict nodes
    to facilitate pairing to a network

    Arguments:
     source_dict: id_key: value

     primary_key: current type of ids used for the nodes.
      Currently can be 'Entrez Gene (GeneID)' or any of the options 
      in the BioServices UniProt mappings.

    kwargs:
     mapping_types: a list of mapping types to include
     verbose

    Returns:
     source_dict, also modified in place

    
    """

    continue_flag = True

    file_key = primary_key
    if primary_key not in available_mapping_source.keys():
        continue_flag = False
        print "Error, you must specify a valid primary_key descriptor to match to in the available database, exiting..."

    if 'mapping_types' in kwargs:
        mapping_types = kwargs['mapping_types'] 
    else:
        mapping_types = default_mapping_target_list

    try:
        from bioservices import UniProt
        u = UniProt(verbose=False)
    except ImportError:
        print("No BioServices module installed or cannot connect, exiting...")
        print("e.g. if you are using pip, did you 'pip install bioservices'?")
        continue_flag = False

    verbose = test_kwarg('verbose', kwargs, [False, True])

    if 'node_id_type' in kwargs:
        node_id_type = kwargs['node_id_type'] 
    else:
        node_id_type = "Entrez Gene (GeneID)"

    verbose = test_kwarg('verbose', kwargs, [False, True])

    # Maximum number of items to
    # query at a time 
    # Note there is a length limit in bioservices 1.2.1
    # for the web-based query string.
    # Trial-and-error suggests the most
    # id's that can be queried are
    # between 100 and 1000
    max_query_length = 500

    if continue_flag:

        the_query_id_list_list = [[]]
        i = 0
        j = 0
        for the_query_id in source_dict.keys():
            if (j + 1) % max_query_length == 0:
                the_query_id_list_list.append([])
                i += 1
                the_query_id_list_list[i] = []
                j = 0
            the_query_id_list_list[i].append(the_query_id)
            j += 1

        the_query_string_list = []
        for i, the_query_id_list in enumerate(the_query_id_list_list):
            query_string = ''
            for the_query_id in the_query_id_list:
                if len(query_string) > 0:
                    query_string = query_string + ' ' + the_query_id
                else:
                    query_string = the_query_id
            the_query_string_list.append(query_string)

        for the_key in source_dict.keys():
            if type(source_dict[the_key]) != dict:
                the_value = source_dict[the_key]
                source_dict[the_key] = {}
                source_dict[the_key]['value'] = the_value
                
        for the_target_type in mapping_types:
            the_result = {}
            for the_query_string in the_query_string_list:
                the_result.update(u.mapping(fr = available_mapping_source[file_key], to = available_mapping_target[the_target_type], query = the_query_string))
            if verbose:
                print("** Finished mapping for %s to %s. **" % (file_key, the_target_type))
            for the_query_id in source_dict.keys():
                if the_query_id in the_result.keys():
                    if len(the_result[the_query_id]) > 0:
                        source_dict[the_query_id][the_target_type] = the_result[the_query_id]
                    else:
                        source_dict[the_query_id][the_target_type] = []
                else:
                    source_dict[the_query_id][the_target_type] = []

    return source_dict



def get_entrez_annotation(id_list, **kwargs):
 
    """Annotates Entrez Gene IDs using Bio.Entrez, in particular epost (to
    submit the data to NCBI) and esummary to retrieve the information. 
    Returns a list of dictionaries with the annotations.

    credit: this function based on code found on biopython.org on Dec 14 2013

    Arguments:
     id_list: a list of Entrez id's each as a string
     
    kwargs:
     e-mail: for ncbi

    Returns: 
     annotations: a dict of Entrez annotations 

     
    """
    import sys
    continue_flag = True
    try:
        from Bio import Entrez
    except:
        print "A functional Biopython is needed for this function."
        print("e.g. if you are using pip, did you 'pip install biopython'?")
        continue_flag = False
        
    if 'email' in kwargs:
        email = kwargs['email']
    else:
        email = ''

    verbose = test_kwarg('verbose', kwargs, [False, True])

    # limit to how many id's can be queried at once 
    max_query = 10000
    
    if continue_flag:
        query_all = False
        query_counter = 0
        annotations = {}
        while not query_all:
            current_query = id_list[(query_counter * max_query) : min(((query_counter + 1) * max_query), len(id_list))]
            request = Entrez.epost("gene",id=",".join(current_query))
            try:
                result = Entrez.read(request)

                webEnv = result["WebEnv"]
                queryKey = result["QueryKey"]
                data = Entrez.esummary(db="gene", webenv=webEnv, query_key =
                                   queryKey)
                annotation_list = Entrez.read(data)
            
                for the_entry in annotation_list:
                    if 'Id' in the_entry.keys():
                        annotations[the_entry['Id']] = the_entry
                        annotations[the_entry['Id']].pop('Id')
                query_counter += 1
                if (query_counter * max_query) > len(id_list):
                    query_all = True
                
            except RuntimeError as e:
                print "An error occurred while retrieving some of the annotations."
                print "The error returned was %s" % e
                
            if verbose:
                print "Retrieved %d annotations for %d genes" % (len(annotations), len(current_query))
 
        return annotations
    else:
        return {}


def get_go_terms(uniprot_acc_id, **kwargs):
    """ Query for GO terms given a 'UniProt ACC'.

    Arguments: 
     uniprot_acc_id: UniProt accession identifier

     kwargs:
      get_descriptions: [True (default), False]
       whether to get include descriptions for the
       GO terms.

    Output:
     the_go_terms: a list if get_descriptions == False
      or a dict if get_descriptions == True

    """
    continue_flag = True
    get_descriptions = test_kwarg('get_descriptions', kwargs, [True, False])

    try:
        from bioservices import QuickGO
        go = QuickGO(verbose=False)
    except ImportError:
        print("No BioServices module installed or cannot connect, exiting...")
        print("e.g. if you are using pip, did you 'pip install bioservices'?")
        continue_flag = False

    if not get_descriptions:
        the_go_terms = []
    else:
        the_go_terms = {}
    
    if continue_flag:
        if not get_descriptions:
            the_go_query = go.Annotation(protein=uniprot_acc_id, format="tsv",
            source="UniProt", col="goID")
            the_go_query = the_go_query.split('\n')
            # The first element is a header
            the_go_query.pop(0)
            # The last element is empty
            the_go_query.pop(len(the_go_terms) - 1)
            the_go_terms = list(set(the_go_query))
            the_go_terms.sort()
        else:
            the_go_query = go.Annotation(protein=uniprot_acc_id, format="tsv",
            source="UniProt", col="goID,goName")
            the_go_query = the_go_query.split('\n')
            # The first element is a header
            the_go_query.pop(0)
            # The last element is empty
            the_go_query.pop(len(the_go_query) - 1)           
            for the_line in the_go_query:
                the_line = the_line.split("\t")
                the_go_id = the_line[0]
                the_description = the_line[1]
                the_go_terms[the_go_id] = the_description
    return the_go_terms



def get_go_descriptions(the_go_id, **kwargs):
    """ Query for GO descriptions given a GO identifier.

    Arguments: 
     go_id: A valid GO identifier of the form: 'GO:XXXXXXX'

     kwargs:
      just a pass-through

    Output:
     the_go_string

    """
    continue_flag = True

    try:
        from bioservices import QuickGO
        go = QuickGO(verbose=False)
    except ImportError:
        print("No BioServices module installed or cannot connect, exiting...")
        print("e.g. if you are using pip, did you 'pip install bioservices'?")
        continue_flag = False

    the_go_string = ''
    if continue_flag:
        the_go_bs = go.Term(the_go_id)

        for the_header in the_go_bs.getchildren():
            if the_header.tag == 'term':
                for the_entry in the_header.getchildren():
                    # Here are some known headers we should get back but we only
                    # really care about one
                    # the_entry.tag == 'id'
                    if the_entry.tag == 'name':
                        the_go_string += the_entry.text
                    # This might be useful
                    # at some point
                    # the_entry.tag == 'namespace'
                    # the_entry.tag == 'def'
                    # the_entry.tag == 'synonym'
                    # We don't care about cross-references
                    # the_entry.tag == 'xref'
                    # print the_entry.text
                    # the_entry.tag == 'is_a':
                    # the_entry.text

    return the_go_string
