from .. core.Edge import Edge
from .. core.Node import Node
from .. core.NodeType import NodeType
from .. core.Network import Network

saved_node_attribute_list = ['notes', 'annotation', 'source', '_nodetype']
saved_edge_attribute_list = ['weight', 'notes', 'annotation']
saved_network_attribute_list = ['notes', 'annotation']

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


def load_pickled_network(filename, path = ""):
    """ Load network object from file

    Arguments:
     filename: effectively a prefix, don't include .pickle, .npy, etc..., these will be added
     dir: directory to load from

    
    """
    import cPickle
    import os
    
    fp = open(path + filename + ".pickle", "rb")
    the_network_dict = cPickle.load(fp)
    fp.close()

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
    for the_edge_id in the_edge_ordered_dict.keys():
        the_edge_dict = the_edge_ordered_dict[the_edge_id]
        the_edge = the_network.connect_node_pair([the_network.nodetypes[0].nodes.get_by_id(the_edge_dict['nodes'][0]), the_network.nodetypes[0].nodes.get_by_id(the_edge_dict['nodes'][1])], the_weight = the_edge_dict['weight'])
        the_edge.id = the_edge_id
        for the_attribute in saved_edge_attribute_list:
            if the_attribute in dir(the_edge):
                setattr(the_edge, the_attribute, the_edge_ordered_dict[the_edge.id][the_attribute])
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

    if 'verbose' in kwargs:
        verbose = kwargs['verbose']
    else:
        verbose = False

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

    if verbose:
        print "     ... the nodes are created."
        print "Linking the nodes, this may take a while..."

    for i, node_1 in enumerate(the_nodes_1):
        the_network.connect_node_pair([the_network.nodetypes[0].nodes.get_by_id(the_nodes_1[i]), the_network.nodetypes[0].nodes.get_by_id(the_nodes_2[i])])
        if verbose:
            if i % 10000 == 0:
                print "Completed linking %s of %s node pairs" %(str(i), str(len(the_nodes_1)))

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
