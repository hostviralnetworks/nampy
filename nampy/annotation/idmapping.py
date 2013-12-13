# There is a web page that gives the list of correct `database identifiers
#  <http://www.uniprot.org/faq/28>

# ID mapping options from bioservoces.uniprot
available_mapping_target = {"UniProtKB ACC": "ACC",
               "UniProtKB ID": "ID",
               "UniParc": "UPARC",
               "UniRef50": "NF50",
               "UniRef90": "NF90",
               "UniRef100": "NF100",
               "EMBL/GenBank/DDBJ": "EMBL_ID",
               "EMBL/GenBank/DDBJ CDS": "EMBL",
               "PIR": "PIR",
               "UniGene": "UNIGENE_ID",
               "Entrez Gene (GeneID)": "P_ENTREZGENEID",
               "GI number*":"P_GI", 
               "IPI": "P_IPI",
               "RefSeq Protein": "P_REFSEQ_AC",
               "RefSeq Nucleotide": "REFSEQ_NT_ID",
               "PDB": "PDB_ID",
               "DisProt": "DISPROT_ID",
               "HSSP": "HSSP_ID",
               "DIP": "DIP_ID",
               "MINT": "MINT_ID",
               "Allergome": "ALLERGOME_ID",
               "MEROPS": "MEROPS_ID",
               "mycoCLAP": "MYCOCLAP_ID",
               "PeroxiBase": "PEROXIBASE_ID",
               "PptaseDB": "PPTASEDB_ID",
               "REBASE": "REBASE_ID",
               "TCDB": "TCDB_ID",
               "PhosSite": "PHOSSITE_ID",
               "DMDM": "DMDM_ID",
               "Aarhus/Ghent-2DPAGE": "AARHUS_GHENT_2DPAGE_ID",
               "World-2DPAGE": "WORLD_2DPAGE_ID",
               "DNASU": "DNASU_ID",
               "Ensembl": "ENSEMBL_ID",
               "Ensembl Protein": "ENSEMBL_PRO_ID",
               "Ensembl Transcript": "ENSEMBL_TRS_ID",
               "Ensembl Genomes": "ENSEMBLGENOME_ID",
               "Ensembl Genomes Protein": "ENSEMBLGENOME_PRO_ID",
               "Ensembl Genomes Transcript": "ENSEMBLGENOME_TRS_ID",
               "GeneID": "P_ENTREZGENEID",
               "GenomeReviews": "GENOMEREVIEWS_ID",
               "KEGG": "KEGG_ID",
               "PATRIC": "PATRIC_ID",
               "UCSC": "UCSC_ID",
               "VectorBase": "VECTORBASE_ID",
               "AGD": "AGD_ID",
               "ArachnoServer": "ARACHNOSERVER_ID",
               "CGD": "CGD",
               "ConoServer": "CONOSERVER_ID",
               "CYGD": "CYGD_ID",
               "dictyBase": "DICTYBASE_ID",
               "EchoBASE": "ECHOBASE_ID",
               "EcoGene": "ECOGENE_ID",
               "euHCVdb": "EUHCVDB_ID",
               "EuPathDB": "EUPATHDB_ID",
               "FlyBase": "FLYBASE_ID",
               "GeneCards": "GENECARDS_ID",
               "GeneFarm": "GENEFARM_ID",
               "GenoList": "GENOLIST_ID",
               "H-InvDB": "H_INVDB_ID",
               "HGNC": "HGNC_ID",
               "HPA": "HPA_ID",
               "LegioList": "LEGIOLIST_ID",
               "Leproma": "LEPROMA_ID",
               "MaizeGDB": "MAIZEGDB_ID",
               "MIM": "MIM_ID",
               "MGI": "MGI_ID",
               "neXtProt": "NEXTPROT_ID",
               "Orphanet": "ORPHANET_ID",
               "PharmGKB": "PHARMGKB_ID",
               "PomBase": "POMBASE_ID",
               "PseudoCAP": "PSEUDOCAP_ID",
               "RGD": "RGD_ID",
               "SGD": "SGD_ID",
               "TAIR": "TAIR_ID",
               "TubercuList": "TUBERCULIST_ID",
               "WormBase": "WORMBASE_ID",
               "WormBase Transcript": "WORMBASE_TRS_ID",
               "WormBase Protein": "WORMBASE_PRO_ID",
               "Xenbase": "XENBASE_ID",
               "ZFIN": "ZFIN_ID",
               "eggNOG": "EGGNOG_ID",
               "GeneTree": "GENETREE_ID",
               "HOGENOM": "HOGENOM_ID",
               "HOVERGEN": "HOVERGEN_ID",
               "KO": "KO_ID",
               "OMA": "OMA_ID",
               "OrthoDB": "ORTHODB_ID",
               "ProtClustDB": "PROTCLUSTDB_ID",
               "BioCyc": "BIOCYC_ID",
               "Reactome": "REACTOME_ID",
               "UniPathWay": "UNIPATHWAY_ID",
               "CleanEx": "CLEANEX_ID",
               "GermOnline": "GERMONLINE_ID",
               "ChEMBL": "CHEMBL_ID",
               "ChiTaRS": "CHITARS_ID",
               "DrugBank": "DRUGBANK_ID",
               "GenomeRNAi": "GENOMERNAI_ID",
               "NextBio": "NEXTBIO_ID"}

available_mapping_source = {"UniProtKB ACC/ID": "ACC+ID",
               "UniParc": "UPARC",
               "UniRef50": "NF50",
               "UniRef90": "NF90",
               "UniRef100": "NF100",
               "EMBL/GenBank/DDBJ": "EMBL_ID",
               "EMBL/GenBank/DDBJ CDS": "EMBL",
               "PIR": "PIR",
               "UniGene": "UNIGENE_ID",
               "Entrez Gene (GeneID)": "P_ENTREZGENEID",
               "GI number*":"P_GI", 
               "IPI": "P_IPI",
               "RefSeq Protein": "P_REFSEQ_AC",
               "RefSeq Nucleotide": "REFSEQ_NT_ID",
               "PDB": "PDB_ID",
               "DisProt": "DISPROT_ID",
               "HSSP": "HSSP_ID",
               "DIP": "DIP_ID",
               "MINT": "MINT_ID",
               "Allergome": "ALLERGOME_ID",
               "MEROPS": "MEROPS_ID",
               "mycoCLAP": "MYCOCLAP_ID",
               "PeroxiBase": "PEROXIBASE_ID",
               "PptaseDB": "PPTASEDB_ID",
               "REBASE": "REBASE_ID",
               "TCDB": "TCDB_ID",
               "PhosSite": "PHOSSITE_ID",
               "DMDM": "DMDM_ID",
               "Aarhus/Ghent-2DPAGE": "AARHUS_GHENT_2DPAGE_ID",
               "World-2DPAGE": "WORLD_2DPAGE_ID",
               "DNASU": "DNASU_ID",
               "Ensembl": "ENSEMBL_ID",
               "Ensembl Protein": "ENSEMBL_PRO_ID",
               "Ensembl Transcript": "ENSEMBL_TRS_ID",
               "Ensembl Genomes": "ENSEMBLGENOME_ID",
               "Ensembl Genomes Protein": "ENSEMBLGENOME_PRO_ID",
               "Ensembl Genomes Transcript": "ENSEMBLGENOME_TRS_ID",
               "GeneID": "P_ENTREZGENEID",
               "GenomeReviews": "GENOMEREVIEWS_ID",
               "KEGG": "KEGG_ID",
               "PATRIC": "PATRIC_ID",
               "UCSC": "UCSC_ID",
               "VectorBase": "VECTORBASE_ID",
               "AGD": "AGD_ID",
               "ArachnoServer": "ARACHNOSERVER_ID",
               "CGD": "CGD",
               "ConoServer": "CONOSERVER_ID",
               "CYGD": "CYGD_ID",
               "dictyBase": "DICTYBASE_ID",
               "EchoBASE": "ECHOBASE_ID",
               "EcoGene": "ECOGENE_ID",
               "euHCVdb": "EUHCVDB_ID",
               "EuPathDB": "EUPATHDB_ID",
               "FlyBase": "FLYBASE_ID",
               "GeneCards": "GENECARDS_ID",
               "GeneFarm": "GENEFARM_ID",
               "GenoList": "GENOLIST_ID",
               "H-InvDB": "H_INVDB_ID",
               "HGNC": "HGNC_ID",
               "HPA": "HPA_ID",
               "LegioList": "LEGIOLIST_ID",
               "Leproma": "LEPROMA_ID",
               "MaizeGDB": "MAIZEGDB_ID",
               "MIM": "MIM_ID",
               "MGI": "MGI_ID",
               "neXtProt": "NEXTPROT_ID",
               "Orphanet": "ORPHANET_ID",
               "PharmGKB": "PHARMGKB_ID",
               "PomBase": "POMBASE_ID",
               "PseudoCAP": "PSEUDOCAP_ID",
               "RGD": "RGD_ID",
               "SGD": "SGD_ID",
               "TAIR": "TAIR_ID",
               "TubercuList": "TUBERCULIST_ID",
               "WormBase": "WORMBASE_ID",
               "WormBase Transcript": "WORMBASE_TRS_ID",
               "WormBase Protein": "WORMBASE_PRO_ID",
               "Xenbase": "XENBASE_ID",
               "ZFIN": "ZFIN_ID",
               "eggNOG": "EGGNOG_ID",
               "GeneTree": "GENETREE_ID",
               "HOGENOM": "HOGENOM_ID",
               "HOVERGEN": "HOVERGEN_ID",
               "KO": "KO_ID",
               "OMA": "OMA_ID",
               "OrthoDB": "ORTHODB_ID",
               "ProtClustDB": "PROTCLUSTDB_ID",
               "BioCyc": "BIOCYC_ID",
               "Reactome": "REACTOME_ID",
               "UniPathWay": "UNIPATHWAY_ID",
               "CleanEx": "CLEANEX_ID",
               "GermOnline": "GERMONLINE_ID",
               "ChEMBL": "CHEMBL_ID",
               "ChiTaRS": "CHITARS_ID",
               "DrugBank": "DRUGBANK_ID",
               "GenomeRNAi": "GENOMERNAI_ID",
               "NextBio": "NEXTBIO_ID"}

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
    """ Script to add more ids to model notes

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
        print("No bioservices module or cannot connect, exiting...")
        continue_flag = False

    if 'node_id_type' in kwargs:
        node_id_type = kwargs['node_id_type'] 
    else:
        node_id_type = "Entrez Gene (GeneID)"

    if 'mapping_types' in kwargs:
        mapping_types = kwargs['mapping_types'] 
    else:
        mapping_types = default_mapping_target_list

    if 'verbose' in kwargs:
        verbose = kwargs['verbose'] 
    else:
        verbose = False

    if continue_flag:
        query_string = ''
        model_node_ids = []
        for the_nodetype in the_network.nodetypes:
            model_node_ids += [x.id for x in the_nodetype.nodes]

        # Note there is a length limit in bioservices 1.1.3
        # on read lengths back from the server from queries.
        # You may want to adjust the modulo term,
        # to ~500 query items (returned list must 
        # have < 2000 items) or adjust this read(2000)
        # to read() in bioservices.uniprot
        the_node_id_list_list = [[]]
        i = 0
        for the_node_id in model_node_ids:
            if (i + 1) % 1000000 == 0:
                the_node_id_list_list.append([])
                i += 1
                node_id_list_list[i] = []
            the_node_id_list_list[i].append(the_node_id)

        query_string_list = []
        for i, the_node_id_list in enumerate(the_node_id_list_list):
            query_string = ''
            for the_node_id in the_node_id_list:
                if len(query_string) > 0:
                    query_string = query_string + ' ' + the_node_id
                else:
                    query_string = the_node_id
            query_string_list.append(query_string)
        
        #for the_node_id in model_node_ids:
        #    if len(query_string) > 0:
        #        query_string = query_string + ' ' + the_node_id
        #    else:
        #        query_string = the_node_id
        target_id_list = []
        for the_target_type in mapping_types:
            the_result = {}
            for the_query_string in query_string_list:
                the_result.update(u.mapping(fr = available_mapping_source[node_id_type], to = available_mapping_target[the_target_type], query = the_query_string))
            if verbose:
                print("Got mapping for %s to %s." % (node_id_type, the_target_type))
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



def get_more_source_dict_ids(source_dict, **kwargs):
    """ Script to add more ids to source dict nodes
    to facilitate pairing to a network

    Arguments:
     source_dict: id_key: value

    kwargs:
     primary_key: current type of ids used for the nodes.
      Currently can be 'Entrez Gene (GeneID)' or any of the options 
      in the BioServices UniProt mappings.
     mapping_types: a list of mapping types to include

    Returns:
     source_dict, also modified in place

    
    """

    continue_flag = True
    
    if 'primary_key' in kwargs:
        file_key = kwargs['primary_key'] 
    else:
        file_key = 'default'

    if 'mapping_types' in kwargs:
        mapping_types = kwargs['mapping_types'] 
    else:
        mapping_types = default_mapping_target_list

    try:
        from bioservices import UniProt
        u = UniProt(verbose=False)
    except:
        print("No bioservices module or cannot connect, exiting...")
        continue_flag = False

    if 'node_id_type' in kwargs:
        node_id_type = kwargs['node_id_type'] 
    else:
        node_id_type = "Entrez Gene (GeneID)"

    if continue_flag:
        if file_key != 'default':
            mapping_dict = {}
            for the_type in mapping_types:
                mapping_dict[the_type] = available_mapping_target[the_type]
            mapping_dict[file_key] = available_mapping_source[file_key]
            if len(mapping_types) > 0:
                from bioservices import UniProt
                u = UniProt(verbose=False)
                query_string = ''
                query_ids = [x for x in source_dict.keys()]
                for the_query_id in query_ids:
                    if len(query_string) > 0:
                        query_string = query_string + ' ' + the_query_id
                    else:
                        query_string = the_query_id
                for the_key in source_dict.keys():
                    if type(source_dict[the_key]) != dict:
                        the_value = source_dict[the_key]
                        source_dict[the_key] = {}
                        source_dict[the_key]['value'] = the_value
                for the_target_type in mapping_types:
                    the_result = u.mapping(fr = available_mapping_source[file_key], to = available_mapping_target[the_target_type], query = query_string)
                    print("Got mapping for %s to %s." % (file_key, the_target_type))
                    for the_query_id in query_ids:
                        if (the_query_id in the_result.keys()):
                            if len(the_result[the_query_id]) > 0:
                                source_dict[the_query_id][the_target_type] = the_result[the_query_id]
                            else:
                                source_dict[the_query_id][the_target_type] = []
                        else:
                            source_dict[the_query_id][the_target_type] = []

    return source_dict
