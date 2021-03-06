from .. core.Edge import Edge
from .. core.Node import Node
from .. core.NodeType import NodeType
from .. core.Network import Network
from ..core.shared_functions import test_kwarg

from ..core.parameters import saved_node_attribute_list
from ..core.parameters import saved_edge_attribute_list
from ..core.parameters import saved_network_attribute_list


def is_float_try(the_string):
    """ Simple test to check if a string can be represented as a 
    floating point.

    """
    try:
        float(the_string)
        return True
    except ValueError:
        return False

    
def string_to_list(the_string, force_to_float = True, sep_char =  ", "):
    """ Convert strings read from file to Python lists
    """
    the_string = the_string.lstrip("[")
    the_string = the_string.rstrip("]")
    the_list = the_string.split(sep_char)
    if the_list[0].startswith("'") & the_list[0].endswith("'"):
        the_list = [x.lstrip("'").rstrip("'") for x in the_list]
    elif the_list[0].startswith("\'") & the_list[0].endswith("\'"):
        the_list = [x.lstrip("\'").rstrip("\'") for x in the_list]    
    for index, the_value in enumerate(the_list):
        
        if is_float_try(the_value) & (force_to_float == True):
            the_list[index] = float(the_value)
        else:
            the_list[index] = the_value
    return the_list    


def pickle_network(the_network, the_filename, **kwargs):
    """ Break apart and save the network

    Arguments:
     the_network: object to save
     filename: effectively a prefix, don't include .pickle, .npy, etc..., these will be added

    kwargs:
     path: directory to save to

    
    """
    # collections does not appear to be well preserved
    # across platforms, so switched to a
    # list to improve pickle portability
    import cPickle
    # don't need to copy since the object
    # will be written to file.
    # from copy import deepcopy
    from numpy import save
    from numpy import load
    import os

    if 'path' in kwargs:
        path = kwargs['path']
    else:
        path = ""

    # We avoid dependence on
    # package classes.  This will ensure models can be
    # fully ported even with base class updates.
    the_node_list_dict = [] 
    the_edge_list_dict = []
    if len(the_network.nodetypes) == 1:
        the_nodes = the_network.nodetypes[0].nodes
        for the_node in the_nodes:
            the_dict = {}
            the_dict[the_node.id] = {}
            for the_attribute in saved_node_attribute_list:
                the_dict[the_node.id][the_attribute] = getattr(the_node, the_attribute)
            the_node_list_dict.append(the_dict)
        for the_edge in the_network.edges:
            the_dict = {}
            the_dict[the_edge.id] = {}
            for the_attribute in saved_edge_attribute_list:
                the_dict[the_edge.id][the_attribute] = getattr(the_edge, the_attribute)
            the_dict[the_edge.id]['nodes'] = [the_edge._nodes[0].id, the_edge._nodes[1].id]
            the_edge_list_dict.append(the_dict)        
        the_network_dict = {}
        the_network_dict['nodes'] = the_node_list_dict
        the_network_dict['edges'] = the_edge_list_dict
        for the_attribute in saved_network_attribute_list:
            if the_attribute in dir(the_network):
                the_network_dict[the_attribute] = getattr(the_network, the_attribute)
        the_network_dict['id'] = the_network.id
        if not the_filename.endswith(".pickle"):
            the_filename += ".pickle"
        fp = open(path + the_filename, "wb")
        cPickle.dump(the_network_dict, fp)
        fp.close()
    else:
        print "Save error. Convert to monopartite network to save."


def load_pickled_network(the_filename, **kwargs):
    """ Load network object from file

    Arguments:
     the_filename: effectively a prefix, don't include .pickle, .npy, etc..., these will be added
    kwargs:
     path: directory to load from

    
    """
    import cPickle
    import os

    if 'path' in kwargs:
        path = kwargs['path']
    else:
        path = ""

    verbose = test_kwarg('verbose', kwargs, [False, True])

    if not the_filename.endswith(".pickle"):
        the_filename += ".pickle"
    
    if verbose:
        print 'Loading network file %s ...' %(path + the_filename)
        
    fp = open(path + the_filename, "rb")
    the_network_dict = cPickle.load(fp)
    fp.close()

    if verbose:
        print '... file loaded. Creating network nodes ...'    

    the_network = Network(the_network_dict['id'])
    for the_attribute in saved_network_attribute_list:
        if the_attribute in dir(the_network):
            setattr(the_network, the_attribute, the_network_dict[the_attribute])

    the_node_list_dict = the_network_dict['nodes']
    the_edge_list_dict = the_network_dict['edges']

    the_node_ids = [x.keys()[0] for x in the_node_list_dict]
    the_nodetype = the_network.nodetypes[0]
    the_nodetype.add_nodes(the_node_ids)
    for the_index, the_node in enumerate(the_nodetype.nodes):
        for the_attribute in saved_node_attribute_list:
            if the_attribute in dir(the_node):
                setattr(the_node, the_attribute, the_node_list_dict[the_index][the_node.id][the_attribute])

    if verbose:
        print '... nodes created. Linking nodes ...'  

    the_node_pair_list = []
    for the_index, the_edge_dict in enumerate(the_edge_list_dict):
        # We allow for updating the edge id here with the new default separation character
        the_edge_dict = the_edge_dict[the_edge_dict.keys()[0]]
        the_node_pair = [the_network.nodetypes[0].nodes.get_by_id(the_edge_dict['nodes'][0]), the_network.nodetypes[0].nodes.get_by_id(the_edge_dict['nodes'][1])]
        the_node_pair_list.append(the_node_pair)

    the_edge_list = the_network.connect_node_pair_list(the_node_pair_list)
                
    for the_index, the_edge in enumerate(the_edge_list):
        the_edge_dict = the_edge_list_dict[the_index]
        the_edge_dict = the_edge_dict[the_edge_dict.keys()[0]]
        for the_attribute in saved_edge_attribute_list:
            if the_attribute in dir(the_edge):
                setattr(the_edge, the_attribute, the_edge_dict[the_attribute])
        
    the_network.update()
    
    return the_network


def create_network_model_from_textfile(network_id, network_file, **kwargs):
    """ Script to build a basic network model
    from a network model. Will reserve annotation,
    initialization for helper functions to call later

    Arguments:
     network_id: string, name of the model
     file_name: name of the network file.  2 columns tab delimited, with node1 node2

    kwargs:
     none, just a pass-through

    Returns:
     network_model


    """

    verbose = test_kwarg('verbose', kwargs, [False, True])

    if verbose:
        print "Reading the network file..."
        
    fp = open(network_file, 'rU')
    the_list = fp.readlines()
    the_list = [x.rstrip("\n") for x in the_list]
    fp.close()

    if verbose:
        print "     ... read the file."  
        print "Creating the nodes..."
    
    the_nodes_1 = []
    the_nodes_2 = []

    for the_entry in the_list:
        the_entry = the_entry.split("\t")
        the_nodes_1.append(the_entry[0])
        the_nodes_2.append(the_entry[1])
    the_nodes = list(set(the_nodes_1 + the_nodes_2))
    the_nodes.sort()

    the_network = Network(network_id)
    the_nodetype = the_network.nodetypes[0]
    the_nodetype.add_nodes(the_nodes)
    # we will define nodes as the default 'monopartite' nodetype
    for the_node in the_nodetype.nodes:
        the_node.set_nodetype(the_nodetype.id)
        
    if verbose:
        print "     ... the nodes are created."
        print "Linking the nodes, this may take a while..."

    the_nodes_1 = [the_network.nodetypes[0].nodes.get_by_id(the_id) for the_id in the_nodes_1]
    the_nodes_2 = [the_network.nodetypes[0].nodes.get_by_id(the_id) for the_id in the_nodes_2]
    # Get rid of duplicate entries
    the_node_pairs = [[the_nodes_1[i], the_nodes_2[i]] for i in range(0, len(the_nodes_1))]
    # Get rid of pairs where node pairs with itself
    the_node_pairs = [x for x in the_node_pairs if (x[0] != x[1])]
        
    the_edge_list = the_network.connect_node_pair_list(the_node_pairs, **kwargs)
    
    the_network.update()

    return the_network


def create_source_dict_from_textfile(source_file, **kwargs):
    """ Script to create a dict of source values.  This
    will attempt to pull more ids with bioservices if more ids are
    specified.

    Arguments:
     source_file: name of the file that contains the sources.
      1 or 2 columns tab delimited, with node /t value_1
      node: id will be used for pairing to netowrk model nodes
      value_1: optional, assumed to be 1.

    kwargs: 
     none, just a pass through

    Returns:
     source_dict

                  
    """
    continue_flag = True
    
    try:
        fp = open(source_file, 'rU')
        the_list = fp.readlines()
        the_list = [x.rstrip("\n") for x in the_list]
        limit_dict = {}
    except:
        print('Cannot load the source file.  Exiting...')
        continue_flag = False

    source_dict = {}
    if continue_flag:
        source_dict = {}
        for the_entry in the_list:
            the_entry = the_entry.split("\t")
            file_node_id = the_entry[0]
            source_dict[file_node_id] = {}
            if len(the_entry) == 1:
                the_value = 1.
            else:
                the_value = float(the_entry[1])
            source_dict[file_node_id]['value'] = the_value

    return source_dict


def read_table_file_to_dict(filename, **kwargs):
    """simple function to extract a dict of dicts as
    {top_key:{key: value}} from an 
    Unix / tab-delineated file.  This is useful for getting 
    node / edge attributes stored in text format.  The first row
    is assumed to be header information.

    Arguments:
     filename

    kwargs:
     top_key: if left blank, the first entry in the first line is taken to
      be the highest level key for the 
     subfield_key_list: a list of fields to include.  Leave blank if all
                        fields are to be included.
     comment_char: Lines starting with this character/string
                  will be ignored
     sep_char: Separation character
     interpret_lists: [True (default), False]
     columns_on_top: [False (default), True]
      if True, the columns will be used
      as the top level keys.

    Returns:
     organized_return_dict
      
    
    """

    if 'top_key' in kwargs:
        top_key = kwargs['top_key']
    else:
        top_key = ''

    if 'subfield_key_list' in kwargs:
        subfield_key_list = kwargs['subfield_key_list']
    else:
        subfield_key_list = []

    if 'comment_char' in kwargs:
        comment_char = kwargs['comment_char']
    else:
        comment_char = ""

    if 'sep_char' in kwargs:
        sep_char = kwargs['sep_char']
    else:
        sep_char = "\t"

    columns_on_top = test_kwarg('columns_on_top',kwargs,[False, True])
        
    force_to_float = test_kwarg('force_to_float', kwargs, [True, False])
    interpret_lists = test_kwarg('interpret_lists', kwargs, [True, False])
    
    fp = open(filename, 'rU')
    the_list = fp.readlines()
    
    the_list = [x.rstrip("\n") for x in the_list]
    
    if (len(comment_char) > 0):
        found_first = False
        while not(found_first):
            if the_list[0].startswith(comment_char):
                the_list.pop(0)
            else:
                found_first = True
                
    if top_key == '':
        firstline = the_list[0].split(sep_char)
        key_index = 0
        top_key = firstline[0]
    else:
        firstline = the_list[0].split(sep_char)
        key_index = firstline.index(top_key)     

    if (len(subfield_key_list) > 0):
        valueindices = [firstline.index(subfield_key) for subfield_key in subfield_key_list]
    else:
        valueindices = [x for x in range(0, len(firstline)) if (x != key_index)]
        subfield_key_list = [firstline[x] for x in valueindices]

    the_list.pop(0)
    
    return_dict={}
    for the_line in the_list:
        if ((len(comment_char) == 0) | (not(the_line.startswith(comment_char)))):
            the_line = the_line.split(sep_char)
            return_dict[the_line[key_index]] = {}
            for index, subfield_key in enumerate(subfield_key_list):
                cur_value = the_line[valueindices[index]]
                if is_float_try(cur_value) & (force_to_float == True):
                    cur_value = float(cur_value)
                else:
                    if ((cur_value.startswith('"')) & (cur_value.endswith(('"')))):
                        # The double-quotes are sometimes inserted before
                        # and after long string/lists
                        # and aren't needed
                        cur_value = cur_value.lstrip('"').rstrip('"')
                    if interpret_lists:
                        if type(cur_value) != list:
                            if ((cur_value.startswith('[')) & (cur_value.endswith((']')))):
                                cur_value = string_to_list(cur_value, force_to_float = force_to_float)
                return_dict[the_line[key_index]][subfield_key] = cur_value

    if columns_on_top:
        organized_return_dict = {}
        for the_column_header in subfield_key_list:
            organized_return_dict[the_column_header] = {}
            for the_row_name in return_dict.keys():
                organized_return_dict[the_column_header][the_row_name] = return_dict[the_row_name][the_column_header]
            
    else:
        organized_return_dict = return_dict
            
    return organized_return_dict


def write_dict_to_textfile(the_filename, the_output_dict, **kwargs):
    """ Write a dict of dicts to a tab-delimited textfile table.

    Arguments:
     the_filename: file to write to
     the_output_dict: the dict to write

    kwargs:
     top_key: string to describe the top key / id in the table.
     subfields_on_top: [False (default), True]
      whether to write to flip the order
      so the subfields become the column headers in
      the table and the top levels keys become the
      row names.

    Returns:
     nothing, just writes to file
     

    """
    from copy import deepcopy

    subfields_on_top = test_kwarg('subfields_on_top',kwargs,[False, True])

    if 'top_key' in kwargs:
        top_key = kwargs['top_key']
    else:
        top_key = 'name'

    if not subfields_on_top:
        all_top_ids = the_output_dict.keys()
        all_subfield_ids = set([])
        for the_key in all_top_ids:
            all_subfield_ids.update(set(the_output_dict[the_key].keys()))
        all_subfield_ids = list(all_subfield_ids)
        all_subfield_ids.sort()
        all_top_ids.sort()
        the_original_dict = deepcopy(the_output_dict)
        the_output_dict = {}
        for the_new_top_key in all_subfield_ids:
            the_output_dict[the_new_top_key] = {}
            for the_new_subfield_key in all_top_ids:
                the_output_dict[the_new_top_key][the_new_subfield_key] = the_original_dict[the_new_subfield_key][the_new_top_key]

    # First we convert to a table to align
    # the columns
    # First scan through and get all the keys                
    all_entries = the_output_dict.keys()
    all_entries.sort()    

    col_names = set([])
    for cur_entry in the_output_dict.keys():
        col_names.update(set(the_output_dict[cur_entry].keys()))
    col_names = [x for x in col_names]
    col_names.sort()

    the_table = []
    the_table.append([])
    the_table[0].append(top_key)
    the_table[0].extend(col_names)
    for row_index, entry in enumerate(all_entries):
        the_table.append([])
        target_row=[]
        target_row.append(entry)
        for target_col_index, cur_col in enumerate(col_names):
            try:
                cur_value = the_output_dict[entry][cur_col]
            except:
                cur_value=''
            target_row.append(cur_value)
        the_table[row_index + 1] = target_row
    
    # Now write the table
    fp = file(the_filename, 'w')
    for row in the_table:
        for colindex, entry in enumerate(row):
            fp.writelines("%s" % entry)
            if colindex == (len(row)-1):
                fp.writelines("\n")
            else:
                fp.writelines("\t")        
    fp.close()


