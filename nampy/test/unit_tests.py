from __future__ import with_statement
import sys
from warnings import warn
from unittest import TestCase, TestLoader, TextTestRunner
from tempfile import gettempdir
from os import unlink
from os.path import join
from copy import deepcopy
from  nampy import __version__ as nampy_version
try:
    from unittest import skipIf
except:
    try:
        from unittest2 import skipIf
    except:
        skipIf = None

# deal with absolute imports by adding the appropriate directory to the path
if __name__ == "__main__":
    sys.path.insert(0, "../..")
    from nampy.test import data_directory, create_test_model, nampy_directory
    # import any functions here?   
    from nampy.core.Network import Network  
    from nampy.core.Node import Node 
    from nampy.core.Object import Object
    from nampy.core.DictList import DictList    
    from nampy.core.parameters import edge_id_separator 
    from nampy.test import humannet_filename as test_network_text_filename    
    from nampy.test import recon_1_filename as test_reporter_input_network_filename
    sys.path.pop(0)
    assert 0
else:
    from ..core.Network import Network  
    from ..core.Node import Node
    from ..core.Object import Object
    from ..core.DictList import DictList    
    from ..core.parameters import edge_id_separator         
    from . import data_directory, create_test_model, nampy_directory
    from . import humannet_filename as test_network_text_filename
    from . import recon_1_filename as test_reporter_input_network_filename      
    sys.path.insert(0, nampy_directory)
    # import any functions here?
    sys.path.pop(0)

# libraries which may or may not be installed
libraries = ["rpy2", "Bio", "bioservices"]
# these are not required bu they are useful

installed_library_list = []
for library in libraries:
    try:
        exec("import %s" % library)
        installed_library_list.append(library) 
    except ImportError:
        exec("%s = None" % library)


class NampyTestCase(TestCase):
    def setUp(self):
        # Make a small network model that can be used 
        # for subsequent tests
        self.model = Network('toy_network')
        the_node_id_list = [str(x) for x in range(0,10)]
        the_node_number_list = [Node(the_number) for the_number in the_node_id_list]
        for the_node in the_node_number_list:
            the_node.set_nodetype('metabolite')
        the_node_id_list = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j']
        the_node_letter_list = [Node(the_letter) for the_letter in the_node_id_list]
        for the_node in the_node_letter_list:
            the_node.set_nodetype('reaction')
        self.model.nodetypes[0].add_nodes(the_node_number_list + the_node_letter_list)
        the_node_id_pair_list = [['0', 'a'],['a', '9'],['1', 'b'],['b', '8'],['2', 'c'],['c', '7'],['3', 'd'],['d', '6'],['4', 'e'],['e', '5'],['5', 'f'],['f', '4'],['6', 'g'],['g', '3'],['7', 'h'],['h', '2'],['8', 'i'],['i', '1'],['9','j'],['j','0']]
        for the_index, the_pair in enumerate(the_node_id_pair_list):
            the_node_id_pair_list[the_index] = [self.model.nodetypes[0].nodes.get_by_id(the_id) for the_id in the_pair]
        self.model.connect_node_pair_list(the_node_id_pair_list)
        self.model.update()
        

class TestDictList(NampyTestCase):
    """ Test core elements ported from COBRApy"""
    def test_object(self):
        self.object = Object("test_object_id")
        self.dictlist = DictList()
        self.dictlist.append(self.object)
        self.assertEqual(self.dictlist[0].id, "test_object_id")

class TestNampyCore(NampyTestCase):
    """ Test core NAMpy functions"""
    def test_add_node(self):     
        the_old_node_id_list = [x.id for x in self.model.nodetypes[0].nodes]
        the_new_node = Node("pi")
        self.model.nodetypes[0].add_nodes([the_new_node])
        the_old_node_id_list.append("pi")
        self.assertEqual(the_old_node_id_list, [x.id for x in self.model.nodetypes[0].nodes])

    def test_remove_node(self):  
        the_old_node_id_list = [x.id for x in self.model.nodetypes[0].nodes]
        the_old_node_ids_to_remove = [x for x in the_old_node_id_list[0:2]]
        self.model.nodetypes[0].remove_nodes(the_old_node_ids_to_remove)
        the_old_node_id_list.pop(0)
        the_old_node_id_list.pop(0)
        self.assertEqual(the_old_node_id_list, [x.id for x in self.model.nodetypes[0].nodes])

    def test_add_edge(self):
        the_new_node = Node("pi")
        old_edge_length = len(self.model.edges)
        self.model.nodetypes[0].add_nodes([the_new_node])
        the_first_node = self.model.nodetypes[0].nodes[0]
        the_edge = self.model.connect_node_pair([the_first_node, the_new_node], the_weight = 1)
        the_expected_edge_id = the_first_node.id + edge_id_separator + the_new_node.id
        self.assertEqual(the_edge.id, the_expected_edge_id)
        self.assertEqual((old_edge_length + 1), len(self.model.edges))

    def test_remove_edge(self):
        the_edge_ids = [x.id for x in self.model.edges]
        self.model.remove_edges([the_edge_ids[0]])
        the_edge_ids.pop(0)
        self.assertEqual(the_edge_ids, [x.id for x in self.model.edges])        
    
    def test_copy(self):
        self.copied_model = self.model.copy()
        self.assertEqual([x.id for x in self.model.nodetypes[0].nodes], [x.id for x in self.copied_model.nodetypes[0].nodes])
        self.assertEqual([x.id for x in self.model.edges], [x.id for x in self.copied_model.edges])
        self.assertEqual(self.model.matrix, self.copied_model.matrix)
        the_new_node = Node("pi")
        self.copied_model.nodetypes[0].add_nodes([the_new_node])
        self.copied_model.connect_node_pair([the_new_node, self.copied_model.nodetypes[0].nodes[0]])
        self.assertNotEqual(len(self.model.nodetypes[0].nodes), len(self.copied_model.nodetypes[0].nodes))
        self.assertNotEqual(len(self.model.edges), len(self.copied_model.edges))

    def test_node_copy(self):
        the_node = self.model.nodetypes[0].nodes[0]
        the_copied_node = the_node.copy()
        self.assertEqual(the_copied_node.id, the_node.id)
        the_copied_node.notes['test_note'] = ['pi','foo']
        the_copied_node.source = 0.3456
        self.assertNotEqual(the_node.notes, the_copied_node.notes)
        self.assertNotEqual(the_node.source, the_copied_node.source)        

    def test_convert_to_multipartite(self):
        the_nodetype_strings = set([])
        for the_nodetype in self.model.nodetypes:
            for the_node in the_nodetype.nodes:
                the_nodetype_strings.add(the_node.get_nodetype())
        self.model.convert_to_multipartite()
        the_nodetype_strings_after_conversion = set([x.id for x in self.model.nodetypes])
        self.assertEqual(the_nodetype_strings_after_conversion, the_nodetype_strings)

    def test_node_locations(self):
        the_location_dict = self.model.get_node_locations()
        the_node_ids = set([x.id for x in self.model.nodetypes[0].nodes])
        the_nodetype_ids = set([x.id for x in self.model.nodetypes])
        self.assertEqual(the_nodetype_ids, set(the_location_dict.values()))
        self.assertEqual(the_node_ids, set(the_location_dict.keys()))
        self.model.convert_to_multipartite()
        the_location_dict = self.model.get_node_locations()
        the_nodetype_ids = set([x.id for x in self.model.nodetypes])
        self.assertEqual(the_nodetype_ids, set(the_location_dict.values()))
        self.assertEqual(the_node_ids, set(the_location_dict.keys()))

    def test_merge_networks(self):
        # Test network manipulation functions
        from nampy.manipulation import manipulation
        the_network = self.model
        the_edges = the_network.edges[0:16]
        the_node_pairs = [x.get_node_pair() for x in the_edges]
        the_node_id_set_both = set([])
        the_node_id_set_1 = set([])
        the_node_id_set_2 = set([])
        the_node_pairs_1 = the_node_pairs[0:10]
        the_node_pairs_2 = the_node_pairs[8:16]
        the_node_ids_1 = []
        all_node_ids = []
        for the_node_pair in the_node_pairs_1:
            the_node_ids_1 += [the_node_pair[0].id, the_node_pair[1].id]
        the_node_ids_2 = []
        for the_node_pair in the_node_pairs_2:            
            the_node_ids_2 += [the_node_pair[0].id, the_node_pair[1].id]
        the_node_ids_1 = list(set(the_node_ids_1))
        the_node_ids_2 = list(set(the_node_ids_2))
        all_node_ids = list(set(the_node_ids_1 + the_node_ids_2))
        the_subnetwork_1 = manipulation.make_subnetwork(the_network, the_node_ids_1, 'subnetwork_1')
        the_subnetwork_2 = manipulation.make_subnetwork(the_network, the_node_ids_2, 'subnetwork_2')
        the_merged_network = manipulation.merge_networks_by_node(the_subnetwork_1, the_subnetwork_2, 'merged_network')
        the_node_ids = [x.id for x in the_merged_network.nodetypes[0].nodes]
        self.assertEqual(set(the_node_ids), set(all_node_ids))        
              
    
        

class TestNampyMore(NampyTestCase):
    """extended test for additional nampy functions
    used in the examples"""

    def test_read_pickled_and_run_reporter(self):
        from nampy.multipartiteanalysis import reporterfeatures
        from nampy.networkio import networkio
        the_network = networkio.load_pickled_network(test_reporter_input_network_filename)
        the_network.convert_to_multipartite()
        p_values_dict = {}
        for the_reaction in the_network.nodetypes.get_by_id('reaction').nodes:
            p_values_dict[the_reaction.id] = 1E-3
        the_reporter_dict = reporterfeatures.calculate_reporter_scores(the_network, p_values_dict, 'metabolite', 'reaction', number_of_randomizations = 100, verbose = False)
        self.assertEqual(len(the_reporter_dict['p_values']), len(the_network.nodetypes.get_by_id('metabolite').nodes))

    def test_textfile_read_and_propagate(self):
        from nampy.manipulation import manipulation
        from nampy.monopartiteanalysis import prince
        the_network = create_test_model(test_network_file = test_network_text_filename)
        the_node_pairs = [x.get_node_pair() for x in the_network.edges]
        the_node_pairs_1 = the_node_pairs[0:16]
        the_node_ids_1 = []
        for the_node_pair in the_node_pairs_1:
            the_node_ids_1 += [the_node_pair[0].id, the_node_pair[1].id]
        the_node_ids_1 = list(set(the_node_ids_1))
        the_subnetwork_1 = manipulation.make_subnetwork(the_network, the_node_ids_1, 'subnetwork_1')
        for the_node in the_subnetwork_1.nodetypes[0].nodes:
            the_node.source = 1
        the_result, the_permutations = prince.prince(the_subnetwork_1, n_permutations = 10, verbose = False)    
        self.assertEqual(len(the_permutations[the_node_ids_1[0]]), 10)



# make a test suite to run all of the tests
loader = TestLoader()
suite = loader.loadTestsFromModule(sys.modules[__name__])

def test_all():
    TextTestRunner(verbosity=2).run(suite)

if __name__ == "__main__":
    test_all()
