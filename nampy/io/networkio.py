from .. core.Edge import Edge
from .. core.Node import Node
from .. core.NodeType import NodeType
from .. core.Network import Network
from ..core.shared_functions import test_kwarg

saved_node_attribute_list = ['notes', 'annotation', 'source', '_nodetype']
saved_edge_attribute_list = ['weight', 'notes', 'annotation']
saved_network_attribute_list = ['notes', 'annotation']

def is_float_try(str):
    """ Simple test to check if a string can be represented as a 
    floating point.

    """
    try:
        float(str)
        return True
    except ValueError:
        return False


def pickle_network(the_network, the_filename, path = ""):
    """ Break apart and save the network

    Arguments:
     the_network: object to save
     filename: effectively a prefix, don't include .pickle, .npy, etc..., these will be added
     path: directory to save to

    
    """
    from collections import OrderedDict
    import cPickle
    # don't need to copy since the object
    # will be written to file.
    # from copy import deepcopy
    from numpy import save
    from numpy import load
    import os

    # We use an ordered dict here to avoid dependence on
    # package classes.  This will ensure models can be
    # fully ported even with base class updates.
    the_node_ordered_dict = OrderedDict()
    the_edge_ordered_dict = OrderedDict()
    if len(the_network.nodetypes) == 1:
        the_nodes = the_network.nodetypes[0].nodes
        for the_node in the_nodes:
            the_dict = {}
            for the_attribute in saved_node_attribute_list:
                the_dict[the_attribute] = getattr(the_node, the_attribute)
            the_node_ordered_dict[the_node.id] = the_dict
        for the_edge in the_network.edges:
            the_dict = {}
            for the_attribute in saved_edge_attribute_list:
                the_dict[the_attribute] = getattr(the_edge, the_attribute)
            the_dict['nodes'] = [the_edge._nodes[0].id, the_edge._nodes[1].id]
            the_edge_ordered_dict[the_edge.id] = the_dict        
        the_network_dict = {}
        the_network_dict['nodes'] = the_node_ordered_dict
        the_network_dict['edges'] = the_edge_ordered_dict
        for the_attribute in saved_network_attribute_list:
            if the_attribute in dir(the_network):
                the_network_dict[the_attribute] = getattr(the_network, the_attribute)
        the_network_dict['id'] = the_network.id
        fp = open(path + the_filename + ".pickle", "wb")
        cPickle.dump(the_network_dict, fp)
        fp.close()
    else:
        print "Save error. Convert to monopartite network to save."


def load_pickled_network(filename, **kwargs):
    """ Load network object from file

    Arguments:
     filename: effectively a prefix, don't include .pickle, .npy, etc..., these will be added
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

    if verbose:
        print 'Loading network file %s ...' %(path + filename + ".pickle")
        
    fp = open(path + filename + ".pickle", "rb")
    the_network_dict = cPickle.load(fp)
    fp.close()

    if verbose:
        print '... file loaded. Creating network nodes ...'    

    the_network = Network(the_network_dict['id'])
    for the_attribute in saved_network_attribute_list:
        if the_attribute in dir(the_network):
            setattr(the_network, the_attribute, the_network_dict[the_attribute])

    the_node_ordered_dict = the_network_dict['nodes']
    the_edge_ordered_dict = the_network_dict['edges']

    the_node_ids = the_node_ordered_dict.keys()
    the_nodetype = the_network.nodetypes[0]
    the_nodetype.add_nodes(the_node_ids)
    for the_node in the_nodetype.nodes:
        for the_attribute in saved_node_attribute_list:
            if the_attribute in dir(the_node):
                setattr(the_node, the_attribute, the_node_ordered_dict[the_node.id][the_attribute])

    if verbose:
        print '... nodes created. Linking nodes ...'  

    the_node_pair_list = []
    for i, the_edge_id in enumerate(the_edge_ordered_dict.keys()):
        the_edge_dict = the_edge_ordered_dict[the_edge_id]
        the_node_pair = [the_network.nodetypes[0].nodes.get_by_id(the_edge_dict['nodes'][0]), the_network.nodetypes[0].nodes.get_by_id(the_edge_dict['nodes'][1])]
        the_node_pair_list.append(the_node_pair)

    the_edge_list = the_network.connect_node_pair_set(the_node_pair_list, **kwargs)
                
    for the_edge in the_edge_list:
        for the_attribute in saved_edge_attribute_list:
            if the_attribute in dir(the_edge):
                setattr(the_edge, the_attribute, the_edge_ordered_dict[the_edge.id][the_attribute])

    # if verbose:
    #         if i % 10000 == 0:
    #             print "Completed linking %s of %s node pairs" %(str(i), str(len(the_edge_ordered_dict.keys())))
        
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

    TODO: might want to compartmentalize this in the future
     if this data becomes available.

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
        
    the_edge_list = the_network.connect_node_pair_set(the_node_pairs, **kwargs)
    
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
    """simple function to extract key:value from an Unix /
    tab-delineated file.  This is useful for getting node / edge
    attributes stored in text format.

    Arguments:
     filename

    kwargs:
     top_key: if left blank, it is assumed the table lacks headers (e.g. the
              first line is just data)
     subfield_key_list: a list of fields to include.  Leave blank if all
                        fields are to be included.
     commentchar: Lines starting with this character/string
                  will be ignored
    
    """

    if 'top_key' in kwargs:
        top_key = kwargs['top_key']
    else:
        top_key = ''

    if 'subfield_key_list' in kwargs:
        subfield_key_list = kwargs['subfield_key_list']
    else:
        subfield_key_list = []

    if 'commentchar' in kwargs:
        commentchar = kwargs['commentchar']
    else:
        commentchar = ""

    force_to_float = test_kwarg('force_to_float', kwargs, [True, False])
    
    fp = open(filename, 'rU')
    the_list = fp.readlines()
    
    the_list = [x.rstrip("\n") for x in the_list]
    
    if (len(commentchar) > 0):
        found_first = False
        while not(found_first):
            if the_list[0].startswith(commentchar):
                the_list.pop(0)
            else:
                found_first = True
                
    if top_key == '':
        firstline = the_list[0].split("\t")
        firstline = [str(i) for i, x in enumerate(firstline)]
        keyindex = 0
    else:
        firstline = the_list[0].split("\t")
        keyindex = firstline.index(top_key)        

    if (len(subfield_key_list) > 0):
        valueindices = [firstline.index(subfield_key) for subfield_key in subfield_key_list]
    else:
        valueindices = [x for x in range(0, len(firstline)) if (x != keyindex)]
        subfield_key_list = [firstline[x] for x in valueindices]

    if top_key != '':
        the_list.pop(0)
    
    return_dict={}
    for the_line in the_list:
        if ((len(commentchar) == 0) | (not(the_line.startswith(commentchar)))):
            the_line = the_line.split("\t")
            return_dict[the_line[keyindex]] = {}
            for index, subfield_key in enumerate(subfield_key_list):
                cur_value = the_line[valueindices[index]]
                if force_to_float:
                    if is_float_try(cur_value):
                        cur_value = float(cur_value)
                return_dict[the_line[keyindex]][subfield_key] = cur_value
            
    return return_dict


def write_dict_to_textfile(the_filename, the_output_dict, top_key):
    """ Write a dict to a tab-delimited textfile.

    Arguments:
     the_filename: file to write to
     the_output_dict: the dict to write
     top_key: string to describe the top key / id in the table.

    """
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
