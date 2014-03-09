from ..core.shared_functions import test_kwarg
from ..core.parameters import available_mapping_target, available_mapping_source


default_mapping_target_list = ["UniProtKB ACC", 
    "UniProtKB ID",
    "Entrez Gene (GeneID)",
    "Symbol"]


def get_more_node_ids(the_network, **kwargs):
    """ Script to add more identifiers to model notes
    based on the node.id

    Arguments:
     the_network: a Network object, modified in place

    kwargs:
     node_id_type: current type of ids used for the nodes.
      Currently can be 'Entrez Gene (GeneID)' or any of the 
      options in the BioServices UniProt mappings
     mapping_types: a list of target mapping id types to include.
      Options can be viewed in core.parameters.py
      Note "Symbol" is an additional option for the
      officieal gene nomenclature symbol.
     email: optional, for NCBI queries.
     verbose: [True (default), False]
     

    Returns:
     the_network 

    TODO: determine the best source db/module
    for pairings from bioservices
                  
    """
    continue_flag = True
    valid_mapping_targets = available_mapping_target.keys() + ['Symbol']
    verbose = test_kwarg('verbose', kwargs, [True, False])

    try:
        from bioservices import UniProt
        # Don't want verbosity at this low of a level
        u = UniProt(verbose = False)
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
        if node_id_type == 'Symbol':
            print "'Symbol' is a special case, not yet able to query with this option, exiting..."
            continue_flag = False            
    else:
        print "No node id type specified, attempting to use 'Entrez Gene (GeneID)'"
        node_id_type = 'Entrez Gene (GeneID)'

    if 'mapping_types' in kwargs:
        mapping_types = [x for x in kwargs['mapping_types'] if x in valid_mapping_targets]
        if len(mapping_types) == 0:
            print('No valid mapping_types selected, exiting...')
            continue_flag = False
        elif 'Symbol' in mapping_types:
            if (('Entrez Gene (GeneID)' not in mapping_types) & (node_id_type != 'Entrez Gene (GeneID)')):
                print "'Symbol' mapping type needs 'Entrez Gene (GeneID)', exiting..."
                continue_flag = False
    else:
        mapping_types = default_mapping_target_list

    if 'email' in kwargs:
        email = kwargs['email'] 
    else:
        email = ''

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
        model_nodes = []
        for the_nodetype in the_network.nodetypes:
            model_nodes += [x for x in the_nodetype.nodes]

        the_node_id_list_list = [[]]
        i = 0
        j = 0
        for the_node in model_nodes:
            if (j + 1) % max_query_length == 0:
                the_node_id_list_list.append([])
                i += 1
                the_node_id_list_list[i] = []
                j = 0
            the_node_id_list_list[i].append(the_node.id)
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
            if the_target_type != 'Symbol':
                the_result = {}
                for the_query_string in query_string_list:
                    the_result.update(u.mapping(fr = available_mapping_source[node_id_type], to = available_mapping_target[the_target_type], query = the_query_string))
                if verbose:
                    print("**Finished mapping for %s to %s.**" % (node_id_type, the_target_type))
                for the_node in model_nodes:
                    if (the_node.id in the_result.keys()):
                        if len(the_result[the_node.id]) > 0:
                            the_node.notes[the_target_type] = the_result[the_node.id]
                        else:
                            the_node.notes[the_target_type] = []
                    else:
                        the_node.notes[the_target_type] = []

        # To avoid a loss of information, we should also make 
        # sure queried IDs are returned in the appropriate 
        # field in case they weren't available in the database.
        if node_id_type in mapping_types:
            # Not yet supported anyway, but can leave this here.
            if node_id_type != 'Symbol':
                for the_node in model_nodes:
                    if the_node.id not in the_node.notes[node_id_type]:
                        the_node.notes[node_id_type].append(the_node.id)
                    
        if "Symbol" in mapping_types:
            if ((node_id_type == "Entrez Gene (GeneID)") | ("Entrez Gene (GeneID)" in mapping_types)):
                the_entrez_to_query = []
                query_dict = {}
                for the_node in model_nodes:
                    query_dict[the_node.id] = {}
                    query_dict[the_node.id]["Entrez Gene (GeneID)"] = []
                    if node_id_type == "Entrez Gene (GeneID)":
                        query_dict[the_node.id]["Entrez Gene (GeneID)"].append(the_node.id)
                    if "Entrez Gene (GeneID)" in mapping_types:
                        the_entrez_list = the_node.notes["Entrez Gene (GeneID)"]
                        if len(the_entrez_list) > 0:
                            for the_entrez_id in the_entrez_list:
                                if the_entrez_id not in query_dict[the_node.id]["Entrez Gene (GeneID)"]:
                                    query_dict[the_node.id]["Entrez Gene (GeneID)"].append(the_entrez_id)
                    the_entrez_to_query += query_dict[the_node.id]["Entrez Gene (GeneID)"]
                the_entrez_to_query = list(set(the_entrez_to_query))
                the_symbol_dict = get_entrez_annotation(the_entrez_to_query, email = email, verbose = verbose)
                for the_node in model_nodes:
                    the_node.notes["Symbol"] = []
                    for the_entrez_id in query_dict[the_node.id]["Entrez Gene (GeneID)"]:
                        the_symbol_id = the_symbol_dict[the_entrez_id]['NomenclatureSymbol']
                        if len(the_symbol_id) > 0:
                            the_node.notes["Symbol"].append(the_symbol_id)
                print("**Finished mapping for %s to %s.**" % (node_id_type, "Symbol"))
            elif verbose:
                print "'Entrez Gene (GeneID)' mappings are needed first in order to query symbols, skipping..."

    return the_network



def get_more_source_dict_ids(source_dict, primary_key_type, **kwargs):
    """ Script to add more ids to source dict nodes
    to facilitate pairing to a network

    Arguments:
     source_dict: id_key: value

     primary_key: current type of ids used for the top level dict key.
      Currently can be 'Entrez Gene (GeneID)' or any of the options 
      in the BioServices UniProt mappings.

    kwargs:
     mapping_types: a list of mapping types to include.
      See core.parameters for the full list.  Note
      'Symbol' is a special case for querying that depends
       on Entrez ID availability.
     verbose: [False (default), True]
     email: optional, for NCBI if querying for 'Symbol'

    Returns:
     source_dict, also modified in place

    
    """

    continue_flag = True
    verbose = test_kwarg('verbose', kwargs, [False, True])
    valid_mapping_targets = available_mapping_target.keys() + ['Symbol']

    if primary_key_type not in available_mapping_source.keys():
        if primary_key_type == 'Symbol':
            print "'Symbol' is a special case, not yet able to query with this as a primary key."
            print "Error, you must specify a valid primary_key_type descriptor to match to in the available database, exiting..."
        continue_flag = False   

    if 'mapping_types' in kwargs:
        mapping_types = [x for x in kwargs['mapping_types'] if x in valid_mapping_targets]
        if len(mapping_types) == 0:
            print('No valid mapping_types selected, exiting...')
            continue_flag = False
        elif 'Symbol' in mapping_types:
            if (('Entrez Gene (GeneID)' not in mapping_types) & (primary_key_type != 'Entrez Gene (GeneID)')):
                print "'Symbol' mapping type needs 'Entrez Gene (GeneID)', exiting..."
                continue_flag = False
    else:
        mapping_types = default_mapping_target_list

    if 'email' in kwargs:
        email = kwargs['email'] 
    else:
        email = ''

    try:
        from bioservices import UniProt
        # Don't want verbosity at this low of a level
        u = UniProt(verbose = False)
    except ImportError:
        print("No BioServices module installed or cannot connect, exiting...")
        print("e.g. if you are using pip, did you 'pip install bioservices'?")
        continue_flag = False

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
            if the_target_type != 'Symbol':
                the_result = {}
                for the_query_string in the_query_string_list:
                    the_result.update(u.mapping(fr = available_mapping_source[primary_key_type], to = available_mapping_target[the_target_type], query = the_query_string))
                if verbose:
                    print("** Finished mapping for %s to %s. **" % (primary_key_type, the_target_type))
                for the_query_id in source_dict.keys():
                    if the_query_id in the_result.keys():
                        if len(the_result[the_query_id]) > 0:
                            source_dict[the_query_id][the_target_type] = the_result[the_query_id]
                        else:
                            source_dict[the_query_id][the_target_type] = []
                    else:
                        source_dict[the_query_id][the_target_type] = []

        # To avoid a loss of information, we should also make 
        # sure queried IDs are returned in the appropriate 
        # field in case they weren't available in the database.
        if primary_key_type in mapping_types:
            # Not yet supported but we can check to avoid breaking this
            if primary_key_type != 'Symbol':
                for the_source_dict_id in source_dict.keys():
                    if the_source_dict_id not in source_dict[the_source_dict_id][primary_key_type]:
                        source_dict[the_source_dict_id][primary_key_type].append(the_source_dict_id)
                    
        if "Symbol" in mapping_types:
            if ((primary_key_type == "Entrez Gene (GeneID)") | ("Entrez Gene (GeneID)" in mapping_types)):
                the_entrez_to_query = []
                # Make query_dict in case "Entrez Gene (GeneID)" was 
                # a primary_key_type but not in mapping_types
                query_dict = {}
                for the_source_dict_id in source_dict.keys():
                    query_dict[the_source_dict_id] = {}
                    query_dict[the_source_dict_id]["Entrez Gene (GeneID)"] = []
                    if primary_key_type == "Entrez Gene (GeneID)":
                        query_dict[the_source_dict_id]["Entrez Gene (GeneID)"].append(the_source_dict_id)
                    if "Entrez Gene (GeneID)" in mapping_types:
                        the_entrez_list = source_dict[the_source_dict_id]["Entrez Gene (GeneID)"]
                        if len(the_entrez_list) > 0:
                            for the_entrez_id in the_entrez_list:
                                if the_entrez_id not in query_dict[the_source_dict_id]["Entrez Gene (GeneID)"]:
                                    query_dict[the_source_dict_id]["Entrez Gene (GeneID)"].append(the_entrez_id)
                    the_entrez_to_query += query_dict[the_source_dict_id]["Entrez Gene (GeneID)"]
                the_entrez_to_query = list(set(the_entrez_to_query))
                the_symbol_dict = get_entrez_annotation(the_entrez_to_query, email = email, verbose = verbose)
                for the_source_dict_id in source_dict.keys():
                    source_dict[the_source_dict_id]["Symbol"] = []
                    for the_entrez_id in query_dict[the_source_dict_id]["Entrez Gene (GeneID)"]:
                        the_symbol_id = the_symbol_dict[the_entrez_id]['NomenclatureSymbol']
                        if len(the_symbol_id) > 0:
                            source_dict[the_source_dict_id]["Symbol"].append(the_symbol_id)
                print("**Finished mapping for %s to %s.**" % (primary_key_type, "Symbol"))
                        

    return source_dict



def get_entrez_annotation(id_list, **kwargs):
 
    """Annotates Entrez Gene IDs using Bio.Entrez, in particular epost (to
    submit the data to NCBI) and esummary to retrieve the information. 
    Returns a list of dictionaries with the annotations.

    credit: this function based on code found on biopython.org on Dec 14 2013

    Arguments:
     id_list: a list of Entrez id's each as a string
     
    kwargs:
     email: for ncbi

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

                webenv = result["WebEnv"]
                query_key = result["QueryKey"]
                data = Entrez.esummary(db="gene", webenv=webenv, query_key =
                                   query_key)
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


def get_go_ids(uniprot_acc_id, **kwargs):
    """ Query for GO ids given a 'UniProt ACC'.

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

    try:
        from bioservices import QuickGO
        # We really don't want verbosity at this low level.
        go = QuickGO(verbose = False)
    except ImportError:
        print("No BioServices module installed or cannot connect, exiting...")
        print("e.g. if you are using pip, did you 'pip install bioservices'?")
        continue_flag = False
    
    continue_flag = True
    get_descriptions = test_kwarg('get_descriptions', kwargs, [True, False])

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



def get_go_terms(the_go_id, **kwargs):
    """ Query for GO terms given a GO identifier.

    Arguments: 
     go_id: A valid GO identifier of the form: 'GO:XXXXXXX'

     kwargs:
      just a pass through

    Output:
     the_go_string

    """
    verbose = test_kwarg('verbose', kwargs, [True, False])
    continue_flag = True

    try:
        from bioservices import QuickGO
        # We really don't want verbosity at this low level.
        go = QuickGO(verbose = False)
    except ImportError:
        print("No BioServices module installed or cannot connect, exiting...")
        print("e.g. if you are using pip, did you 'pip install bioservices'?")
        continue_flag = False

    the_go_string = ''
    if continue_flag:
        the_go_bs = go.Term(the_go_id)

        # bs for beautiful_soup
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


def get_node_go_ids(the_network, **kwargs):
    """ Add Gene Ontology controlled vocabulary IDs and 
    terms to network nodes based on 'UniProtKB ACC'.

    Arguments: 
     the_network.  Note the_network must include 
      a notes field titled 'UniProtKB ACC' in order
      to get go terms for the node.

    kwargs:
     verbose: [True (default), False]
     include_term: [True (default), False]

    Returns:
     the_network: modified in place.  Note that
     all 'GO ID' fields and 'GO Term' fields will be
     overwritten for all network nodes.
     

    """
    continue_flag = True
    verbose = test_kwarg('verbose', kwargs, [True, False])
    include_term = test_kwarg('include_term', kwargs, [True, False])
    # For argument passing
    kwargs['include_term'] = True

    the_node_locations = the_network.get_node_locations()
    if len(the_node_locations) == 0:
        print 'The network has no nodes, exiting...'
        continue_flag = False

    the_uniprot_id_key_node_list_value_dict = {}
    for the_nodetype in the_network.nodetypes:
        for the_node in the_nodetype.nodes:
            if 'UniProtKB ACC' in the_node.notes.keys():
                if len(the_node.notes['UniProtKB ACC']) > 0:
                    for the_uniprot_id in the_node.notes['UniProtKB ACC']:
                        if len(the_uniprot_id) > 0:
                            # used_uniprot_ids.add(the_uniprot_id)
                            if the_uniprot_id not in the_uniprot_id_key_node_list_value_dict.keys():
                                the_uniprot_id_key_node_list_value_dict[the_uniprot_id] = []
                            the_uniprot_id_key_node_list_value_dict[the_uniprot_id].append(the_node)

    uniprot_id_list = the_uniprot_id_key_node_list_value_dict.keys()
    the_n_uniprot_ids = len(uniprot_id_list)
    if the_n_uniprot_ids > 50:
        if verbose:
            print "Warning, querying for GO ID's is currently done individually.  It is recommended to use get_node_go_ids() with smaller networks.  There are currently %i terms to query." %(the_n_uniprot_ids)

    uniprot_ids_without_go = set([])
    the_uniprot_key_golist_value_dict = {}
    the_go_description_dict = {}
    for the_index, the_uniprot_id in enumerate(uniprot_id_list):
        try_flag = True
        try_counter = 1
        while try_flag:
            try:
                current_go_dict = get_go_ids(the_uniprot_id, **kwargs)
                if len(current_go_dict.keys()) > 0:
                    the_uniprot_key_golist_value_dict[the_uniprot_id] = current_go_dict.keys()
                    the_go_description_dict.update(current_go_dict)
                else:
                    uniprot_ids_without_go.add(the_uniprot_id)
                if verbose:
                    if (the_index + 1) % 100 == 0:
                        print 'Completed %i of %i.' %(the_index + 1, the_n_uniprot_ids)
                try_flag = False
            except:
                try_counter += 1
            if try_counter > 10:
                try_flag = False
                if verbose:
                    print 'Unable to retrieve GO ID for %s.' %(the_uniprot_id)
                uniprot_ids_without_go.add(the_uniprot_id)

    uniprot_id_list = list(set(uniprot_id_list) - uniprot_ids_without_go)
    nodes_with_go_id = set([])
    for the_uniprot_id in uniprot_id_list:
        the_go_list = the_uniprot_key_golist_value_dict[the_uniprot_id]
        #for the_go_id in the_go_list:
        #    the_go_id_key_term_value_dict = the_go_description_dict[the_go_id]
        for the_node in the_uniprot_id_key_node_list_value_dict[the_uniprot_id]:
            nodes_with_go_id.add(the_node)
            # Overwrite existing IDs
            the_node.notes['GO ID'] = set([])
            for the_go_id in the_go_list:
                the_node.notes['GO ID'].add(the_go_id)

    for the_node in nodes_with_go_id:
        the_node.notes['GO ID'] = list(the_node.notes['GO ID'])
        the_node.notes['GO ID'].sort()
        the_node.notes['GO Term'] = []
        for the_go_id in the_node.notes['GO ID']:
            the_node.notes['GO Term'].append(the_go_description_dict[the_go_id])
        if not include_term:
            the_node.notes.pop('GO Term')

    # Assign empty lists to remaining nodes
    for the_nodetype in the_network.nodetypes:
        for the_node in the_nodetype.nodes:
            if the_node not in nodes_with_go_id:
                the_node.notes['GO ID'] = []
                the_node.notes['GO Term'] = []
                if not include_term:
                    the_node.notes.pop('GO Term')

    return the_network
