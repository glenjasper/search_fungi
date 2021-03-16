#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import re
import sys
import time
import traceback
import pandas as pd
from tqdm import tqdm
from Bio import Entrez
from pprint import pprint

class Parse:

    def __init__(self):
        self.VERSION = 1.0
        self.ROOT = os.path.dirname(os.path.realpath(__file__))

        # Log
        self.LOG_NAME = "log_%s_%s.log" % (os.path.splitext(os.path.basename(__file__))[0], time.strftime('%Y%m%d'))
        self.LOG_FILE = None # os.path.join(self.ROOT, self.LOG_NAME)
        self.NODES_FILE = 'network_nodes.csv'
        self.EDGES_FILE = 'network_edges.csv'

        self.OUTPUT_PATH = os.path.join(self.ROOT, 'output_parse')
        self.TAXONOMIC_RANK_FILE = 'taxonomic_rank.csv'
        self.TAXONOMIC_RANK_FILE_FILLED = 'taxonomic_rank_filled.csv'
        self.TAXONOMIC_SP_FILE = 'taxonomic_sp.txt'
        self.TAXONOMIC_OTHER_NAMES_FILE = 'taxonomic_other_names.txt'
        self.TAXONOMIC_INCONSISTENCIES = 'taxonomic_inconsistencies.txt'

        self.STATUS_OK = 'Ok'
        self.STATUS_REVIEWED = 'Reviewed'
        self.COLUMN_STATUS = 'status'

        # Entrez
        self.EMAIL = 'youremail@domain.com'
        self.DATABASE = 'taxonomy'
        self.FORMAT = 'xml'

        self.KEY_COUNT = 'count'

        # Tags
        self.KEY_WARNINGLIST = 'WarningList'
        self.KEY_OUTPUTMESSAGE = 'OutputMessage'

        self.KEY_ID_LIST = 'IdList'
        self.KEY_LINEAGE_EX = 'LineageEx'

        self.KEY_TAX_ID = 'TaxId'
        self.KEY_RANK = 'Rank'
        self.KEY_SCIENTIFIC_NAME = 'ScientificName'

        self.RANK_TAXON_ID = 'taxon_oid'
        self.RANK_DOMAIN = 'domain'             # [1]
        self.RANK_SUPERKINGDOM = 'superkingdom' # [1]
        self.RANK_KINGDOM = 'kingdom'           # [2]
        self.RANK_PHYLUM = 'phylum'             # [3]
        self.RANK_CLASS = 'class'               # [4]
        self.RANK_ORDER = 'order'               # [5]
        self.RANK_FAMILY = 'family'             # [6]
        self.RANK_GENUS = 'genus'               # [7]
        self.RANK_SPECIES = 'species'           # [8]

        self.SUPERKINGDOM_EUKARYOTA = 'Eukaryota'
        self.KINGDOM_FUNGI = 'Fungi'

        self.VALUE_UNCLASSIFIED = 'unclassified'

        self.UNASSIGNED_LABEL = 'Unassigned'
        self.UNASSIGNED_INDEX = 1

        # Fonts
        self.RED = '\033[31m'
        self.BIRED = '\033[1;91m'
        self.BIGREEN = '\033[1;92m'
        self.END = '\033[0m'

    def show_print(self, message, logs = None, showdate = True, font = None):
        msg_print = message
        msg_write = message

        if font is not None:
            msg_print = "%s%s%s" % (font, msg_print, self.END)

        if showdate is True:
            _time = time.strftime('%Y-%m-%d %H:%M:%S')
            msg_print = "%s %s" % (_time, msg_print)
            msg_write = "%s %s" % (_time, message)

        print(msg_print)
        if logs is not None:
            for log in logs:
                if log is not None:
                    with open(log, 'a') as f:
                        f.write("%s\n" % msg_write)
                        f.close()

    def start_time(self):
        return time.time()

    def finish_time(self, start, message = None):
        finish = time.time()
        runtime = time.strftime("%H:%M:%S", time.gmtime(finish - start))
        if message is None:
            return runtime
        else:
            return "%s: %s" % (message, runtime)

    def check_path(self, path):
        if len(path) > 0 and os.path.exists(path):
            return True
        else:
            return False

    def create_directory(self, path):
        output = True
        try:
            if len(path) > 0 and not os.path.exists(path):
                os.makedirs(path)
        except Exception as e:
            output = False
        return output

    def get_num_lines(self, input_file):
        num_lines = sum(1 for line in open(input_file))
        return num_lines

    def read_exported_file_img(self, file):

        def is_unclassified(term):
            if term == self.VALUE_UNCLASSIFIED:
                term = ''
            return term

        self.show_print("Raw file: %s" % file, [self.LOG_FILE])
        self.show_print("", [self.LOG_FILE])

        df = pd.read_csv(filepath_or_buffer = file, sep = '\t', header = 0, index_col = False)
        df = df.where(pd.notnull(df), '')

        biosamples = {}
        for _, row in df.iterrows():
            _species = row[self.RANK_SPECIES.capitalize()]

            if _species not in biosamples:
                current = {}
                current[self.RANK_TAXON_ID] = row[self.RANK_TAXON_ID]
                current[self.RANK_SUPERKINGDOM] = is_unclassified(row[self.RANK_DOMAIN.capitalize()])
                current[self.RANK_KINGDOM] = ''
                current[self.RANK_PHYLUM] = is_unclassified(row[self.RANK_PHYLUM.capitalize()])
                current[self.RANK_CLASS] = is_unclassified(row[self.RANK_CLASS.capitalize()])
                current[self.RANK_ORDER] = is_unclassified(row[self.RANK_ORDER.capitalize()])
                current[self.RANK_FAMILY] = is_unclassified(row[self.RANK_FAMILY.capitalize()])
                current[self.RANK_GENUS] = is_unclassified(row[self.RANK_GENUS.capitalize()])
                current[self.RANK_SPECIES] = _species
                current[self.KEY_COUNT] = 1
            else:
                current = biosamples[_species].copy()
                current.update({self.KEY_COUNT: current[self.KEY_COUNT] + 1})

            biosamples.update({_species: current})

        return biosamples

    def separate_types_classifications(self, dict_biosamples):
        species = {}
        species_sp = {}
        others = {}
        for _, data in dict_biosamples.items():
            _count = data[self.KEY_COUNT]
            organism = data[self.RANK_SPECIES]
            _organism = organism.split()

            if re.search('^[a-z]', organism) or \
               len(_organism) == 1 or \
               (len(_organism) >= 2 and (re.search('[A-Z]{2}', _organism[0]) or re.search('[A-Z]{2}', _organism[1]) or _organism[1].lower() == 'x')):
                if organism not in others:
                    others.update({organism: _count})
                else:
                    others.update({organism: others[organism] + 1})
            else:
                '''
                    Cadophora sp. DSE1049
                    Neotyphodium sp. FaTG-4
                    Ostreococcus sp. 'lucimarinus'

                    unclassified
                '''
                if re.search('sp[.]', _organism[1]):
                    if organism not in species_sp:
                        species_sp.update({organism: _count})
                    else:
                        species_sp.update({organism: species_sp[organism] + 1})
                    continue
                else:
                    _organism = ' '.join(_organism[0:2])

                    if _organism not in species:
                        # Only the first id is registered
                        species.update({_organism: data})
                    else:
                        current = species[_organism].copy()
                        current[self.KEY_COUNT] = current[self.KEY_COUNT] + data[self.KEY_COUNT]
                        species.update({_organism: current})

        # pprint(species)
        # pprint(species_sp)
        # pprint(others)
        # exit()

        return species, species_sp, others

    def create_taxonomic_rank_file(self, dict_species, dict_species_sp, dict_others):
        collect_phylums = {}
        dict_species_status = dict_species.copy()
        for term, data in dict_species.items():
            current = data.copy()
            current[self.COLUMN_STATUS] = None
            dict_species_status.update({term: current})

            _phylum = data[self.RANK_PHYLUM]
            if _phylum:
                collect_phylums.update({_phylum: {self.RANK_KINGDOM: None}})

        # pprint(dict_species_status)
        # pprint(collect_phylums)
        # exit()

        detail_csv = {}
        if os.path.exists(self.TAXONOMIC_RANK_FILE):
            df = pd.read_csv(filepath_or_buffer = self.TAXONOMIC_RANK_FILE, sep = '\t', header = 0, index_col = False)
            df = df.where(pd.notnull(df), '')

            for index, row in df.iterrows():
                _term = row[self.RANK_SPECIES]
                _phylum = row[self.RANK_PHYLUM]
                _kingdom = row[self.RANK_KINGDOM]

                current = dict_species_status[_term]
                current[self.COLUMN_STATUS] = self.STATUS_REVIEWED
                current[self.RANK_KINGDOM] = _kingdom
                dict_species_status.update({_term: current})

                _row = '\t'.join(str(value) for value in row.to_numpy())
                detail_csv.update({index + 1: _row})

                if _phylum in collect_phylums:
                    if _kingdom:
                        collect_phylums[_phylum].update({self.RANK_KINGDOM: _kingdom})

        # pprint(dict_species_status)
        # pprint(collect_phylums)
        # exit()

        self.show_print("Obtaining taxonomic information of %s organisms from Entrez..." % len(dict_species_status), [self.LOG_FILE])

        with open(self.TAXONOMIC_RANK_FILE, 'w') as fw:
            fw.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (self.RANK_TAXON_ID,
                                                                       self.KEY_COUNT,
                                                                       self.RANK_SUPERKINGDOM,
                                                                       self.RANK_KINGDOM,
                                                                       self.RANK_PHYLUM,
                                                                       self.RANK_CLASS,
                                                                       self.RANK_ORDER,
                                                                       self.RANK_FAMILY,
                                                                       self.RANK_GENUS,
                                                                       self.RANK_SPECIES,
                                                                       self.COLUMN_STATUS))

            with tqdm(total = len(dict_species_status)) as pbar:
                if detail_csv:
                    for _, line in detail_csv.items():
                        fw.write('%s\n' % line)
                    pbar.update(len(detail_csv))

                for term, detail in dict_species_status.items():
                    _status = detail[self.COLUMN_STATUS]
                    _phylum = detail[self.RANK_PHYLUM]

                    # print('>>>sp: %s | status: %s | phylum: %s' % (term, _status, _phylum))
                    if _status is None:
                        if _phylum:
                            _kingdom = collect_phylums[_phylum][self.RANK_KINGDOM]
                            if not _kingdom:
                                taxonomic_rank = self.get_taxonomic_rank(_phylum)
                                # pprint(taxonomic_rank)

                                _kingdom = ''
                                if self.RANK_KINGDOM in taxonomic_rank:
                                    _kingdom = taxonomic_rank[self.RANK_KINGDOM]
                                    collect_phylums[_phylum].update({self.RANK_KINGDOM: _kingdom})
                        else:
                            '''
                            {'OutputMessage': 'Ok',
                             'class': 'Cryptophyceae',
                             'family': 'Geminigeraceae',
                             'genus': 'Guillardia',
                             'order': 'Pyrenomonadales',
                             'superkingdom': 'Eukaryota'}
                            '''
                            taxonomic_rank = self.get_taxonomic_rank(term)
                            # pprint(taxonomic_rank)
                            
                            _kingdom = ''
                            if self.RANK_KINGDOM in taxonomic_rank:
                                _kingdom = taxonomic_rank[self.RANK_KINGDOM]

                        fw.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (detail[self.RANK_TAXON_ID],
                                                                                   detail[self.KEY_COUNT],
                                                                                   detail[self.RANK_SUPERKINGDOM],
                                                                                   _kingdom, # detail[self.RANK_KINGDOM],
                                                                                   detail[self.RANK_PHYLUM],
                                                                                   detail[self.RANK_CLASS],
                                                                                   detail[self.RANK_ORDER],
                                                                                   detail[self.RANK_FAMILY],
                                                                                   detail[self.RANK_GENUS],
                                                                                   term,
                                                                                   self.STATUS_OK))
                        pbar.update(1)
        fw.close()

        with open(self.TAXONOMIC_SP_FILE, 'w') as fw:
            fw.write('Term\tQuantity\n')
            for term, quant in sorted(dict_species_sp.items()):
                fw.write('%s\t%s\n' % (term, quant))
        fw.close()

        with open(self.TAXONOMIC_OTHER_NAMES_FILE, 'w') as fw:
            fw.write('Term\tQuantity\n')
            for term, quant in sorted(dict_others.items()):
                fw.write('%s\t%s\n' % (term, quant))
        fw.close()
 
        self.show_print("  File with taxonomic information: %s" % self.TAXONOMIC_RANK_FILE, [self.LOG_FILE])
        self.show_print("  File with sp. names: %s" % self.TAXONOMIC_SP_FILE, [self.LOG_FILE])
        self.show_print("  File with other names: %s" % self.TAXONOMIC_OTHER_NAMES_FILE, [self.LOG_FILE])
        self.show_print("", [self.LOG_FILE])

    def get_taxonomic_rank(self, term):
        # print('>>>', term)
        taxonomic_rank = {}
        Entrez.email = self.EMAIL
        handle = Entrez.esearch(term = term, db = self.DATABASE, retmode = self.FORMAT)
        record = Entrez.read(handle)
        # pprint(record)
        # exit()

        valid = True
        id_list = 0
        warning_message = self.STATUS_OK
        if self.KEY_ID_LIST in record:
            id_list = record[self.KEY_ID_LIST]
            if not id_list:
                valid = False
                warning_message = 'ID not found'

        if self.KEY_WARNINGLIST in record:
            valid = False
            warning_message = record[self.KEY_WARNINGLIST][self.KEY_OUTPUTMESSAGE][0]

        if valid:
            # print('  >>>', id_list)
            handle = Entrez.efetch(id = id_list, db = self.DATABASE, retmode = self.FORMAT)
            record = Entrez.read(handle)
            # pprint(record)
            # exit()

            tax_id = record[0][self.KEY_TAX_ID]
            scientific_name = record[0][self.KEY_SCIENTIFIC_NAME]

            for data in record[0][self.KEY_LINEAGE_EX]:
                _rank = data[self.KEY_RANK]
                _scientific_name = data[self.KEY_SCIENTIFIC_NAME]

                if _rank == self.RANK_SUPERKINGDOM:
                    taxonomic_rank.update({self.RANK_SUPERKINGDOM: _scientific_name})
                if _rank == self.RANK_KINGDOM:
                    taxonomic_rank.update({self.RANK_KINGDOM: _scientific_name})
                if _rank == self.RANK_PHYLUM:
                    taxonomic_rank.update({self.RANK_PHYLUM: _scientific_name})
                if _rank == self.RANK_CLASS:
                    taxonomic_rank.update({self.RANK_CLASS: _scientific_name})
                if _rank == self.RANK_ORDER:
                    taxonomic_rank.update({self.RANK_ORDER: _scientific_name})
                if _rank == self.RANK_FAMILY:
                    taxonomic_rank.update({self.RANK_FAMILY: _scientific_name})
                if _rank == self.RANK_GENUS:
                    taxonomic_rank.update({self.RANK_GENUS: _scientific_name})

                '''
                [DictElement({'TaxId': '4890',
                              'ScientificName': 'Ascomycota',
                              'OtherNames': DictElement({'Includes': [],
                                                         'Inpart': [],
                                                         'Anamorph': [],
                                                         'Name': [],
                                                         'GenbankAnamorph': [],
                                                         'CommonName': ['sac fungi'],
                                                         'Acronym': [],
                                                         'EquivalentName': [],
                                                         'Misnomer': [],
                                                         'Teleomorph': [],
                                                         'GenbankSynonym': [],
                                                         'Misspelling': [],
                                                         'Synonym': [],
                                                         'GenbankCommonName': 'ascomycetes',
                                                         'BlastName': 'ascomycetes'}, attributes={}),
                              'ParentTaxId': '451864',
                              'Rank': 'phylum',
                              'Division': 'Plants and Fungi',
                              'GeneticCode': DictElement({'GCId': '1',
                                                          'GCName': 'Standard'}, attributes={}),
                              'MitoGeneticCode': DictElement({'MGCId': '4',
                                                              'MGCName': 'Mold Mitochondrial; Protozoan Mitochondrial; Coelenterate Mitochondrial; Mycoplasma; Spiroplasma'}, attributes={}),
                              'Lineage': 'cellular organisms; Eukaryota; Opisthokonta; Fungi; Dikarya',
                              'LineageEx': [DictElement({'TaxId': '131567',
                                                         'ScientificName': 'cellular organisms',
                                                         'Rank': 'no rank'}, attributes={}),
                                            DictElement({'TaxId': '2759',
                                                         'ScientificName': 'Eukaryota',
                                                         'Rank': 'superkingdom'}, attributes={}),
                                            DictElement({'TaxId': '33154',
                                                         'ScientificName': 'Opisthokonta',
                                                         'Rank': 'clade'}, attributes={}),
                                            DictElement({'TaxId': '4751',
                                                         'ScientificName': 'Fungi',
                                                         'Rank': 'kingdom'}, attributes={}), 
                                            DictElement({'TaxId': '451864',
                                                         'ScientificName': 'Dikarya', 
                                                         'Rank': 'subkingdom'}, attributes={})], 
                              'CreateDate': '1995/02/27 09:24:00', 
                              'UpdateDate': '2017/06/14 10:56:24', 
                              'PubDate': '1993/04/20 01:00:00'}, attributes={})]
                '''

        taxonomic_rank.update({self.KEY_OUTPUTMESSAGE: warning_message})
        # taxonomic_rank.update({self.RANK_SPECIES: term})

        return taxonomic_rank

    def get_dictionaries_index(self):
        self.show_print("Generating IDs for nodes...", [self.LOG_FILE])

        df = pd.read_csv(filepath_or_buffer = self.TAXONOMIC_RANK_FILE, sep = '\t', header = 0, index_col = False)
        df = df.where(pd.notnull(df), None)

        dict_superkingdom = {}
        dict_kingdom = {}
        dict_phylum = {}
        dict_class = {}
        dict_order = {}
        dict_family = {}
        dict_genus = {}
        dict_species = {}
        for _, row in df.iterrows():
            _species = row[self.RANK_SPECIES]
            _superkingdom = row[self.RANK_SUPERKINGDOM]
            _kingdom = row[self.RANK_KINGDOM]
            _phylum = row[self.RANK_PHYLUM]
            _class = row[self.RANK_CLASS]
            _order = row[self.RANK_ORDER]
            _family = row[self.RANK_FAMILY]
            _genus = row[self.RANK_GENUS]
            _status = row[self.COLUMN_STATUS]

            if _status == self.STATUS_OK and _superkingdom == self.SUPERKINGDOM_EUKARYOTA and _kingdom == self.KINGDOM_FUNGI:
                if _superkingdom is not None and _superkingdom not in dict_superkingdom:
                    dict_superkingdom.update({_superkingdom: None})
                if _kingdom is not None and _kingdom not in dict_kingdom:
                    dict_kingdom.update({_kingdom: None})
                if _phylum is not None and _phylum not in dict_phylum:
                    dict_phylum.update({_phylum: None})
                if _class is not None and _class not in dict_class:
                    dict_class.update({_class: None})
                if _order is not None and _order not in dict_order:
                    dict_order.update({_order: None})
                if _family is not None and _family not in dict_family:
                    dict_family.update({_family: None})
                if _genus is not None and _genus not in dict_genus:
                    dict_genus.update({_genus: None})
                if _species is not None and _species not in dict_species:
                    dict_species.update({_species: None})

        index = 1
        for key, _ in dict_superkingdom.items():
            dict_superkingdom[key] = index
            index += 1
        for key, _ in dict_kingdom.items():
            dict_kingdom[key] = index
            index += 1
        for key, _ in dict_phylum.items():
            dict_phylum[key] = index
            index += 1
        for key, _ in dict_class.items():
            dict_class[key] = index
            index += 1
        for key, _ in dict_order.items():
            dict_order[key] = index
            index += 1
        for key, _ in dict_family.items():
            dict_family[key] = index
            index += 1
        for key, _ in dict_genus.items():
            dict_genus[key] = index
            index += 1
        for key, _ in dict_species.items():
            dict_species[key] = index
            index += 1

        # pprint(dict_superkingdom)
        # pprint(dict_kingdom)
        # pprint(dict_phylum)
        # pprint(dict_class)
        # pprint(dict_order)
        # pprint(dict_family)
        # pprint(dict_genus)
        # pprint(dict_species)

        return dict_superkingdom, dict_kingdom, dict_phylum, dict_class, dict_order, dict_family, dict_genus, dict_species

    def create_ghepi_files(self, dict_superkingdom, dict_kingdom, dict_phylum, dict_class, dict_order, dict_family, dict_genus, dict_species):
        self.show_print("Generating Ghepi files...", [self.LOG_FILE])

        detail_label = {}
        for label, index in dict_superkingdom.items():
            detail_label.update({index: [label, 'color_1']})
        for label, index in dict_kingdom.items():
            detail_label.update({index: [label, 'color_2']})
        for label, index in dict_phylum.items():
            detail_label.update({index: [label, 'color_3']})
        for label, index in dict_class.items():
            detail_label.update({index: [label, 'color_4']})
        for label, index in dict_order.items():
            detail_label.update({index: [label, 'color_5']})
        for label, index in dict_family.items():
            detail_label.update({index: [label, 'color_6']})
        for label, index in dict_genus.items():
            detail_label.update({index: [label, 'color_7']})
        for label, index in dict_species.items():
            detail_label.update({index: [label, 'color_8']})
        # pprint(detail_label)

        with open(self.NODES_FILE, 'w') as fw:
            fw.write('Id,Label,Color\n')
            for index, data in detail_label.items():
                line = '%s,%s,%s\n' % (index, data[0], data[1])
                fw.write(line)
        fw.close()

        tuples_nodes = self.get_detail_tuples()

        with open(self.EDGES_FILE, 'w') as fw:
            fw.write('Id,Source,Target,Type\n')

            for index, _tuple in enumerate(tuples_nodes, start = 1):
                rank_1 = _tuple[0][0]
                label_1 = _tuple[0][1]
                rank_2 = _tuple[1][0]
                label_2 = _tuple[1][1]
                # print(rank_1, label_1, rank_2, label_2)

                if rank_1 == self.RANK_SUPERKINGDOM:
                    _source = dict_superkingdom[label_1]
                elif rank_1 == self.RANK_KINGDOM:
                    _source = dict_kingdom[label_1]
                elif rank_1 == self.RANK_PHYLUM:
                    _source = dict_phylum[label_1]
                elif rank_1 == self.RANK_CLASS:
                    _source = dict_class[label_1]
                elif rank_1 == self.RANK_ORDER:
                    _source = dict_order[label_1]
                elif rank_1 == self.RANK_FAMILY:
                    _source = dict_family[label_1]
                elif rank_1 == self.RANK_GENUS:
                    _source = dict_genus[label_1]
                elif rank_1 == self.RANK_SPECIES:
                    _source = dict_species[label_1]

                if rank_2 == self.RANK_SUPERKINGDOM:
                    _target = dict_superkingdom[label_2]
                elif rank_2 == self.RANK_KINGDOM:
                    _target = dict_kingdom[label_2]
                elif rank_2 == self.RANK_PHYLUM:
                    _target = dict_phylum[label_2]
                elif rank_2 == self.RANK_CLASS:
                    _target = dict_class[label_2]
                elif rank_2 == self.RANK_ORDER:
                    _target = dict_order[label_2]
                elif rank_2 == self.RANK_FAMILY:
                    _target = dict_family[label_2]
                elif rank_2 == self.RANK_GENUS:
                    _target = dict_genus[label_2]
                elif rank_2 == self.RANK_SPECIES:
                    _target = dict_species[label_2]

                line = '%s,%s,%s,Undirected\n' % (index, _source, _target)
                fw.write(line)
        fw.close()

        self.show_print("Ghepi files:", [self.LOG_FILE])
        self.show_print("  Nodes file: %s" % self.NODES_FILE, [self.LOG_FILE])
        self.show_print("  Edges file: %s" % self.EDGES_FILE, [self.LOG_FILE])
        self.show_print("", [self.LOG_FILE])

    def get_detail_tuples(self):
        df = pd.read_csv(filepath_or_buffer = self.TAXONOMIC_RANK_FILE, sep = '\t', header = 0, index_col = False)
        df = df.where(pd.notnull(df), None)
        # print(df)

        tuples_nodes = []
        for _, row in df.iterrows():
            _species = row[self.RANK_SPECIES]
            _superkingdom = row[self.RANK_SUPERKINGDOM]
            _kingdom = row[self.RANK_KINGDOM]
            _phylum = row[self.RANK_PHYLUM]
            _class = row[self.RANK_CLASS]
            _order = row[self.RANK_ORDER]
            _family = row[self.RANK_FAMILY]
            _genus = row[self.RANK_GENUS]
            _status = row[self.COLUMN_STATUS]

            if _status == self.STATUS_OK and _superkingdom == self.SUPERKINGDOM_EUKARYOTA and _kingdom == self.KINGDOM_FUNGI:
                taxonomy = [(self.RANK_SUPERKINGDOM, _superkingdom),
                            (self.RANK_KINGDOM, _kingdom),
                            (self.RANK_PHYLUM, _phylum),
                            (self.RANK_CLASS, _class),
                            (self.RANK_ORDER, _order),
                            (self.RANK_FAMILY, _family),
                            (self.RANK_GENUS, _genus),
                            (self.RANK_SPECIES, _species)]
         
                _begin = None
                _end = None
                for index, item in enumerate(taxonomy, start = 1):
                    if _begin is None and item[1] is not None:
                        _begin = item
                        continue

                    if _end is None and item[1] is not None:
                        _end = item

                        _tuple = (_begin, _end)
                        if _tuple not in tuples_nodes:
                            tuples_nodes.append(_tuple)

                        _begin = _end
                        _end = None

        return tuples_nodes

    def fill_dataframe(self):

        def update_df(rank_child, rank_parent):
            dict_tupla = {}
            for rowid, row in filtered_df.iterrows():
                _parent = row[rank_parent]
                _child = row[rank_child]

                if _child not in dict_tupla:
                    dict_tupla.update({_child: {'parent': [{rank_parent: _parent, 'rowid': rowid}]}})
                else:
                    current = dict_tupla[_child]['parent']
                    current.append({rank_parent: _parent, 'rowid': rowid})
                    dict_tupla.update({_child: {'parent': current}})
            # pprint(dict_tupla)

            if rank_child == self.RANK_GENUS:
                genus_unassigned = ''
                if genus_unassigned in dict_tupla:
                    dict_tupla.update({self.UNASSIGNED_LABEL: dict_tupla[genus_unassigned]})
                    del dict_tupla[genus_unassigned]
            # pprint(dict_tupla)

            for _, data in dict_tupla.items():
                parents = data['parent']

                update_index = False
                for item in parents:
                    rowid = item['rowid']
                    parent = item[rank_parent]

                    if parent == '':
                        item[rank_parent] = '%s%s' % (self.UNASSIGNED_LABEL, self.UNASSIGNED_INDEX)
                        update_index = True

                if update_index:
                    self.UNASSIGNED_INDEX += 1
            # pprint(dict_tupla)

            for key, data in dict_tupla.items():
                parents = data['parent']

                current_parents = {}
                for item in parents:
                    rowid = item['rowid']
                    parent = item[rank_parent]

                    # Update dataframe (Only genus)
                    if rank_child == self.RANK_GENUS:
                        if key.startswith(self.UNASSIGNED_LABEL):
                            filtered_df.at[rowid, rank_child] = key

                    # Update dataframe
                    if parent.startswith(self.UNASSIGNED_LABEL):
                        filtered_df.at[rowid, rank_parent] = parent

                    # Check parents
                    if key not in current_parents:
                        current_parents.update({key: [parent]})
                    else:
                        current = current_parents[key].copy()
                        if parent not in current:
                            current.append(parent)
                        current_parents.update({key: current})

                if len(current_parents[key]) > 1:
                    line = '{level_child}: {child}\n'.format(level_child=rank_child, child=key)
                    for parent in current_parents[key]:
                        line = line + '\t{level_parent}: {parent}\n'.format(level_parent=rank_parent, parent=parent)
                    handle_inconsistencies.write('%s\n' % line)
            # print(filtered_df)

        df = pd.read_csv(filepath_or_buffer = self.TAXONOMIC_RANK_FILE, sep = '\t', header = 0, index_col = False)
        df = df.where(pd.notnull(df), '')

        _cond_status = df[self.COLUMN_STATUS] == self.STATUS_OK
        _cond_superkingdom = df[self.RANK_SUPERKINGDOM] == self.SUPERKINGDOM_EUKARYOTA
        _cond_kingdom = df[self.RANK_KINGDOM] == self.KINGDOM_FUNGI

        _cond_phylum = df[self.RANK_PHYLUM] != ''
        _cond_class = df[self.RANK_CLASS] != ''
        _cond_order = df[self.RANK_ORDER] != ''
        _cond_family = df[self.RANK_FAMILY] != ''
        _cond_genus = df[self.RANK_GENUS] != ''

        _cond_null_phylum = df[self.RANK_PHYLUM] == ''
        _cond_null_class = df[self.RANK_CLASS] == ''
        _cond_null_order = df[self.RANK_ORDER] == ''
        _cond_null_family = df[self.RANK_FAMILY] == ''
        _cond_null_genus = df[self.RANK_GENUS] == ''

        _conditions = _cond_status & _cond_superkingdom & _cond_kingdom
        # _conditions = _cond_status & _cond_superkingdom & _cond_kingdom & _cond_phylum & _cond_class & _cond_order & _cond_family & _cond_genus
        # _conditions = _cond_status & _cond_superkingdom & _cond_kingdom & (_cond_null_phylum | _cond_null_class | _cond_null_order | _cond_null_family | _cond_null_genus)
        filtered_df = df[_conditions]
        # print(filtered_df)
        # print(filtered_df.shape)

        handle_inconsistencies = open(self.TAXONOMIC_INCONSISTENCIES, 'w')
        update_df(self.RANK_GENUS, self.RANK_FAMILY)
        update_df(self.RANK_FAMILY, self.RANK_ORDER)
        update_df(self.RANK_ORDER, self.RANK_CLASS)
        update_df(self.RANK_CLASS, self.RANK_PHYLUM)
        handle_inconsistencies.close()
        # print(filtered_df[filtered_df.columns[3:11]])

        if self.get_num_lines(self.TAXONOMIC_INCONSISTENCIES) > 0:
            self.show_print("File with inconsistencies: %s" % self.TAXONOMIC_INCONSISTENCIES, [self.LOG_FILE])
            self.show_print("", [self.LOG_FILE])
        else:
            os.remove(self.TAXONOMIC_INCONSISTENCIES)

        filtered_df.to_csv (self.TAXONOMIC_RANK_FILE_FILLED, index = False, header = True, sep = '\t')

    def create_json(self):

        def get_tuplas(rank_1, rank_2, dict_support):
            dict_tupla = {}
            for _, row in filtered_df.iterrows():
                _level_1 = row[rank_1]
                _level_2 = row[rank_2]

                if _level_1 not in dict_tupla:
                    dict_tupla.update({_level_1: {_level_2: dict_support[_level_2]}})
                else:
                    current = dict_tupla[_level_1]
                    current.update({_level_2: dict_support[_level_2]})
                    dict_tupla.update({_level_1: current})
            return dict_tupla

        # df = pd.read_csv(filepath_or_buffer = self.TAXONOMIC_RANK_FILE, sep = '\t', header = 0, index_col = False)
        df = pd.read_csv(filepath_or_buffer = self.TAXONOMIC_RANK_FILE_FILLED, sep = '\t', header = 0, index_col = False)
        df = df.where(pd.notnull(df), '')

        _cond_status = df[self.COLUMN_STATUS] == self.STATUS_OK
        _cond_superkingdom = df[self.RANK_SUPERKINGDOM] == self.SUPERKINGDOM_EUKARYOTA
        _cond_kingdom = df[self.RANK_KINGDOM] == self.KINGDOM_FUNGI
        _cond_phylum = df[self.RANK_PHYLUM] != ''
        _cond_class = df[self.RANK_CLASS] != ''
        _cond_order = df[self.RANK_ORDER] != ''
        _cond_family = df[self.RANK_FAMILY] != ''
        _cond_genus = df[self.RANK_GENUS] != ''

        # _conditions = _cond_status & _cond_superkingdom & _cond_kingdom
        _conditions = _cond_status & _cond_superkingdom & _cond_kingdom & _cond_phylum & _cond_class & _cond_order & _cond_family & _cond_genus
        filtered_df = df[_conditions]
        # print(filtered_df)
        # print(filtered_df.shape)

        dict_genus = {}
        for _, row in filtered_df.iterrows():
            _genus = row[self.RANK_GENUS]
            _species = row[self.RANK_SPECIES]

            if _genus not in dict_genus:
                dict_genus.update({_genus: {_species: {}}})
            else:
                current = dict_genus[_genus]
                current.update({_species: {}})
                dict_genus.update({_genus: current})
        # pprint(dict_genus)

        dict_family = get_tuplas(self.RANK_FAMILY, self.RANK_GENUS, dict_genus)
        # pprint(dict_family)

        dict_order = get_tuplas(self.RANK_ORDER, self.RANK_FAMILY, dict_family)
        # pprint(dict_order)

        dict_class = get_tuplas(self.RANK_CLASS, self.RANK_ORDER, dict_order)
        # pprint(dict_class)

        dict_phylum = get_tuplas(self.RANK_PHYLUM, self.RANK_CLASS, dict_class)
        # pprint(dict_phylum)

        dict_kingdom = get_tuplas(self.RANK_KINGDOM, self.RANK_PHYLUM, dict_phylum)
        # pprint(dict_kingdom)

        return dict_kingdom

    def create_json_d3(self, json, level):
        json_d3 = {}
        for key_kingdom, data_phylum in json.items():
            # print(key_kingdom)

            json_d3.update({'name': key_kingdom,
                            'size': len(data_phylum),
                            'children': []})

            _index_phylum = 0
            for key_phylum, data_class in data_phylum.items():
                # print(key_phylum)

                _children1 = {'name': key_phylum,
                              'size': len(data_class),
                              'children': []}
                # json_d3['children'].append(_children1)
                json_d3['children'].insert(_index_phylum, _children1)

                _index_class = 0
                for key_class, data_order in data_class.items():
                    # print(key_class)

                    _children2 = {'name': key_class,
                                  'size': len(data_order),
                                  'children': []}
                    json_d3['children'][_index_phylum]['children'].insert(_index_class, _children2)

                    _index_order = 0
                    for key_order, data_family in data_order.items():
                        # print(key_order)

                        _children3 = {'name': key_order,
                                      'size': len(data_family),
                                      'children': []}
                        json_d3['children'][_index_phylum]['children'][_index_class]['children'].insert(_index_order, _children3)

                        _index_family = 0
                        for key_family, data_genus in data_family.items():
                            # print(key_family)

                            _children4 = {'name': key_family,
                                          'size': len(data_genus),
                                          'children': []}
                            json_d3['children'][_index_phylum]['children'][_index_class]['children'][_index_order]['children'].insert(_index_family, _children4)

                            _index_genus = 0
                            for key_genus, data_species in data_genus.items():
                                # print(key_genus)

                                if level == 'genus':
                                    # Without species
                                    _children5 = {'name': key_genus,
                                                  'size': len(data_species)}
                                    json_d3['children'][_index_phylum]['children'][_index_class]['children'][_index_order]['children'][_index_family]['children'].insert(_index_genus, _children5)
                                elif level == 'species':
                                    # With species
                                    _children5 = {'name': key_genus,
                                                  'size': len(data_species),
                                                  'children': []}
                                    json_d3['children'][_index_phylum]['children'][_index_class]['children'][_index_order]['children'][_index_family]['children'].insert(_index_genus, _children5)

                                    _index_species = 0
                                    for key_species, _ in data_species.items():
                                        # print(key_species)

                                        _children6 = {'name': key_species,
                                                      'size': 1}
                                        json_d3['children'][_index_phylum]['children'][_index_class]['children'][_index_order]['children'][_index_family]['children'][_index_genus]['children'].insert(_index_species, _children6)

                                        _index_species += 1

                                _index_genus += 1

                            _index_family += 1

                        _index_order += 1

                    _index_class += 1

                _index_phylum += 1

        return json_d3

    def create_json_d3_file(self, level = 'species'):
        self.fill_dataframe()
        json = self.create_json()
        json_d3_raw = self.create_json_d3(json, level)

        if level == 'species':
            json_name = 'network-species.json'
        elif level == 'genus':
            json_name = 'network-genus.json'

        json_d3 = os.path.join(self.OUTPUT_PATH, json_name)
        with open(json_d3, 'w') as fw:
            fw.write(str(json_d3_raw).replace('\'', '"'))
        fw.close()

        self.show_print("File JSON: %s" % json_d3, [self.LOG_FILE])
        self.show_print("", [self.LOG_FILE])

def main(args):
    try:
        start = oparse.start_time()
        oparse.create_directory(oparse.OUTPUT_PATH)
        oparse.LOG_FILE = os.path.join(oparse.OUTPUT_PATH, oparse.LOG_NAME)
        oparse.EDGES_FILE = os.path.join(oparse.OUTPUT_PATH, oparse.EDGES_FILE)
        oparse.NODES_FILE = os.path.join(oparse.OUTPUT_PATH, oparse.NODES_FILE)
        oparse.TAXONOMIC_RANK_FILE = os.path.join(oparse.OUTPUT_PATH, oparse.TAXONOMIC_RANK_FILE)
        oparse.TAXONOMIC_RANK_FILE_FILLED = os.path.join(oparse.OUTPUT_PATH, oparse.TAXONOMIC_RANK_FILE_FILLED)
        oparse.TAXONOMIC_SP_FILE = os.path.join(oparse.OUTPUT_PATH, oparse.TAXONOMIC_SP_FILE)
        oparse.TAXONOMIC_OTHER_NAMES_FILE = os.path.join(oparse.OUTPUT_PATH, oparse.TAXONOMIC_OTHER_NAMES_FILE)
        oparse.TAXONOMIC_INCONSISTENCIES = os.path.join(oparse.OUTPUT_PATH, oparse.TAXONOMIC_INCONSISTENCIES)

        oparse.show_print("###########################################################", [oparse.LOG_FILE], font = oparse.BIGREEN)
        oparse.show_print("########################### RUN ###########################", [oparse.LOG_FILE], font = oparse.BIGREEN)
        oparse.show_print("###########################################################", [oparse.LOG_FILE], font = oparse.BIGREEN)

        ############################################
        # Crear el fill_file y el .json
        # Disponer el taxonomic_rank.csv curado
        oparse.create_json_d3_file()
        oparse.create_json_d3_file(level = 'genus')
        exit()
        ############################################

        # raw_path = 'C:\\Users\\Glen\\Dropbox\\Developmet\\network\\scripts\\img-m\\data'
        raw_path = '/home/g13nj45p3r/Dropbox/Developmet/network/scripts/img-m/data'

        # Minimus
        # raw_file = 'taxontable_eukarya-min.xls'

        # All
        raw_file = 'taxontable_eukarya.xls'

        raw_file = os.path.join(raw_path, raw_file)

        biosamples = oparse.read_exported_file_img(raw_file)
        # pprint(biosamples)
        dict_species, dict_species_sp, dict_others = oparse.separate_types_classifications(biosamples)
        # print(len(dict_species), len(dict_species_sp), len(dict_others))
        oparse.create_taxonomic_rank_file(dict_species, dict_species_sp, dict_others)

        dict_superkingdom, dict_kingdom, dict_phylum, dict_class, dict_order, dict_family, dict_genus, dict_species = oparse.get_dictionaries_index()
        oparse.create_ghepi_files(dict_superkingdom, dict_kingdom, dict_phylum, dict_class, dict_order, dict_family, dict_genus, dict_species)

        oparse.create_json_d3_file()
        oparse.create_json_d3_file(level = 'genus')

        oparse.show_print(oparse.finish_time(start, "Elapsed time"), [oparse.LOG_FILE])
        oparse.show_print("Done.", [oparse.LOG_FILE])
    except Exception as e:
        oparse.show_print("\n%s" % traceback.format_exc(), [oparse.LOG_FILE], font = oparse.RED)
        oparse.show_print(oparse.finish_time(start, "Elapsed time"), [oparse.LOG_FILE])
        oparse.show_print("Done.", [oparse.LOG_FILE])

if __name__ == '__main__':
    oparse = Parse()
    main(sys.argv)
