#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import re
import sys
import time
import requests
import traceback
import statistics
import pandas as pd
import xml.etree.ElementTree as ET
from tqdm import tqdm
from Bio import Entrez
from colorama import init
init()
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
        self.TAXONOMIC_NOM_INVAL_FILE = 'taxonomic_nom_inval.txt'
        self.TAXONOMIC_OTHER_NAMES_FILE = 'taxonomic_other_names.txt'
        self.TAXONOMIC_INCONSISTENCIES = 'taxonomic_inconsistencies.txt'

        self.STATUS_OK = 'Ok'
        self.COLUMN_STATUS = 'status'
        self.STATUS_REVIEWED = 'Reviewed'

        # Entrez
        self.EMAIL = 'youremail@domain.com'
        self.DATABASE = 'taxonomy'
        self.FORMAT = 'xml'

        self.KEY_BIOSAMPLE = 'biosample'
        self.KEY_COUNT = 'count'

        # Tags
        self.KEY_WARNINGLIST = 'WarningList'
        self.KEY_OUTPUTMESSAGE = 'OutputMessage'

        self.KEY_ID_LIST = 'IdList'
        self.KEY_LINEAGE_EX = 'LineageEx'

        self.KEY_TAX_ID = 'TaxId'
        self.KEY_RANK = 'Rank'
        self.KEY_SCIENTIFIC_NAME = 'ScientificName'

        self.RANK_SUPERKINGDOM = 'superkingdom' # [1] Domain
        self.RANK_KINGDOM = 'kingdom'           # [2]
        self.RANK_SUBKINGDOM = 'subkingdom'
        self.RANK_SUPERPHYLUM = 'superphylum'
        self.RANK_PHYLUM = 'phylum'             # [3]
        self.RANK_SUBPHYLUM = 'subphylum'
        self.RANK_SUPERCLASS = 'superclass'
        self.RANK_CLASS = 'class'               # [4]
        self.RANK_SUBCLASS = 'subclass'
        self.RANK_SUPERORDER = 'superorder'
        self.RANK_ORDER = 'order'               # [5]
        self.RANK_SUBORDER = 'suborder'
        self.RANK_SUPERFAMILY = 'superfamily'
        self.RANK_FAMILY = 'family'             # [6]
        self.RANK_SUBFAMILY = 'subfamily'
        self.RANK_SUPERGENUS = 'supergenus'
        self.RANK_GENUS = 'genus'               # [7]
        self.RANK_SUBGENUS = 'subgenus'
        self.RANK_SPECIES = 'species'           # [8]
        self.RANK_CLADE = 'clade'
        self.RANK_NO_RANK = 'no rank'

        self.RANK_SPECIES_DB = 'species (Taxonomy)'

        self.SUPERKINGDOM_EUKARYOTA = 'Eukaryota'
        self.KINGDOM_FUNGI = 'Fungi'

        self.UNASSIGNED_LABEL = 'Unassigned'
        self.UNASSIGNED_INDEX = 1

        # Specials
        # https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=119167
        self.SPECIAL_NAME_MIXED_EST = 'mixed EST library'
        # 'Nom. inval.' (= nomen invalidum, = invalid name) refers to a name not published in accordance with rules enumerated
        #                                                   in the International Code for algae, fungi and plants.
        # 'Nom. nud.' (= nomen nudum) refers to a name published without a diagnosis.
        self.SPECIAL_NAME_NOM_INVAL = '(nom. inval.)'

        # Index Fungorum
        self.if_flag = 'false'
        self.if_hits = '10'
        self.URL_SEARCHTEXT = 'http://www.indexfungorum.org/ixfwebservice/fungus.asmx/NameSearch?SearchText=<SPECIES_NAME>&AnywhereInText=%s&MaxNumber=<HITS>' % (self.if_flag)
        self.URL_SEARCHKEY = 'http://www.indexfungorum.org/ixfwebservice/fungus.asmx/NameByKey?NameKey=<NAME_KEY>'

        self.SPECIAL_NAME_NOM_ILLEGITIMATE = 'Nom. illegit.'
        self.SPECIAL_NAME_NOM_INVALID = 'Nom. inval.'

        self.IF_COLUMN_KEY = 'if_current_id'
        self.IF_COLUMN_BASIONYM_KEY = 'if_basionym_id'
        self.IF_COLUMN_COMMENT = 'comment'
        self.IF_COLUMN_RANK = 'rank'
        self.IF_COLUMN_BASIONYM = 'basionym'

        self.IF_TAG_CURRENT_NAME_KEY = 'CURRENT_x0020_NAME_x0020_RECORD_x0020_NUMBER'
        self.IF_TAG_BASIONYM_NAME_KEY = 'BASIONYM_x0020_RECORD_x0020_NUMBER'
        self.IF_TAG_NAME_KEY = 'RECORD_x0020_NUMBER'
        self.IF_TAG_RANK = 'INFRASPECIFIC_x0020_RANK'

        self.IF_TAG_COMMENT = 'NOMENCLATURAL_x0020_COMMENT'

        self.IF_TAG_KINGDOM = 'Kingdom_x0020_name'
        self.IF_TAG_PHYLUM = 'Phylum_x0020_name'
        self.IF_TAG_CLASS = 'Class_x0020_name'
        self.IF_TAG_ORDER = 'Order_x0020_name'
        self.IF_TAG_FAMILY = 'Family_x0020_name'
        self.IF_TAG_GENUS = 'Genus_x0020_name'
        self.IF_TAG_QUERY_NAME = 'NAME_x0020_OF_x0020_FUNGUS'

        self.IF_RANK_SPECIES = 'sp.'
        self.IF_RANK_GENUS = 'gen.'
        self.IF_RANK_FAMILY = 'fam.'
        self.IF_RANK_ORDER = 'ord.'
        self.IF_RANK_CLASS = 'class.'
        self.IF_RANK_PHYLUM = 'phyl.'
        self.IF_RANK_KINGDOM = 'regn.'

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
                    with open(log, 'a', encoding = 'utf-8') as f:
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

    def read_exported_file_ncbi(self, file):
        '''
        1: Naganishia_liquefaciens_N6
        Identifiers: BioSample: SAMD00235376; SRA: DRS144800
        Organism: Naganishia liquefaciens
        Attributes:
            /sample name="Naganishia_liquefaciens_N6"
            /collection date="1995-07-19/1995-07-25"
            /broad-scale environmental context="marine biome"
            /local-scale environmental context="trench"
            /environmental medium="sediment"
            /estimated size="19.5 Mbp"
            /geographic location="Pacific Ocean"
            /isolation and growth condition="missing"
            /latitude and longitude="40 N 144 E"
            /locus_tag_prefix="NliqN6"
            /number of replicons="missing"
            /ploidy="diploid"
            /project name="Genome assembly of a deep-sea yeast, Naganishia (Cryptococcus) liquefaciens strain N6"
            /propagation="asexual"
            /strain="N6"
        Description:
        Keywords: GSC:MIxS;MIGS:5.0
        Accession: SAMD00235376 ID: 15496769
        '''

        self.show_print("Raw file: %s" % file, [self.LOG_FILE])
        self.show_print("", [self.LOG_FILE])

        biosamples = {}
        current = {}
        index = None
        with open(file, 'r', encoding = "utf8") as fr:
            for line in fr:
                line = line.strip()

                # if re.search('^[0-9]{1,6}[:]', line):
                if re.search('^[0-9]+[:]', line):
                    if current:
                        biosamples.update({index: current})
                        current = {}

                    index = int(line.split(':')[0].strip())
                    sample_name = line.split(':')[1].strip()
                    current.update({'sample-name': sample_name})
                    continue

                if re.search('^Organism:', line):
                    current.update({'organism': line.split(':')[1].strip()})

                _search = re.search('BioSample[:][\s]+[A-Z]+[0-9]+[;]*', line)
                biosample_id = None
                if _search:
                    biosample_id = _search.group()
                    biosample_id = biosample_id.replace('BioSample:', '')
                    biosample_id = biosample_id.replace(';', '')
                    biosample_id = biosample_id.strip()
                    current.update({'biosample': biosample_id})
        fr.close()

        if current:
            biosamples.update({index: current})

        return biosamples

    def separate_types_classifications(self, dict_biosamples):

        def save_specie(this_organism, this_biosample, dict_species):
            if this_organism not in dict_species:
                # Only the first sample_id is registered
                dict_species.update({this_organism: {self.KEY_COUNT: 1, self.KEY_BIOSAMPLE: this_biosample}})
            else:
                current = dict_species[this_organism].copy()
                current[self.KEY_COUNT] = current[self.KEY_COUNT] + 1
                dict_species.update({this_organism: current})

        species = {}
        species_sp_remove = {}
        species_nom_inval = {}
        others = {}
        for _, data in dict_biosamples.items():
            '''
                Cercophora sp.
                Corticiaceae sp.
                Monographella sp.
                Monosporascus sp. 5C6A
                Fusarium sp. F.17.10

                Amanita aff. grandis
                Boletellus aff. emodensis fb171
                Aspergillus aff. floccosus IMV 01167
                Borofutus cf. dhakanus fb36
                Cercospora cf. flagellaris
                Microbotryum cf. violaceum FR02 5.2 A2D

                Saccharomycopsis fibuligera x Saccharomycopsis cf. fibuligera
                Metschnikowia matae var. maris
                Metschnikowia matae var. matae
                Chaetomium thermophilum var. dissitum
                Puccinia coronata var. avenae f. sp. avenae
                Clonostachys rosea f. catenulata
                Cronartium quercuum f. sp. banksianae
                Cronartium quercuum f. sp. fusiforme
                Alternaria alternata f. sp. lycopersici

                Stachybotrys elegans/Rhizoctonia solani mixed EST library
                Triticum aestivum/Phaeosphaeria nodorum mixed EST library
                Glycine max/Fusarium solani f. sp. glycines mixed EST library

                Pseudopestalotiopsis camelliae-sinensis
                Microbotryum tragopogonis-pratensis
                Raffaelea quercus-mongolicae

                # Tylopilus sp.
                # Tylopilus aff.
                # Tylopilus cf.

                # affinis: sp. aff. affin.
                # confer: cf. cfr.
            '''
            biosample = data['biosample']
            organism = data['organism']
            organism = organism.replace('[', '').replace(']', '')

            _organism = organism.split()

            # print('>>>', organism)
            if re.search('^Candidatus[\s]', organism):
                _organism = _organism[1:]
                organism = ' '.join(_organism)

            # print('>>>', organism)
            if re.search(self.SPECIAL_NAME_NOM_INVAL, organism):
                if organism not in species_nom_inval:
                    species_nom_inval.update({organism: 1})
                else:
                    species_nom_inval.update({organism: species_nom_inval[organism] + 1})
            elif re.search('^[a-z]', organism) or \
                 re.search('^[A-Z][a-z]+[\s][A-Z]+', organism) or \
                 (len(_organism) >= 2 and (re.search('[A-Z]{2}', _organism[0]) or \
                                           re.search('[A-Z]{2}', _organism[1]) or \
                                           _organism[1].lower() == 'x')) or \
                 len(_organism) == 1:
                if organism not in others:
                    others.update({organism: 1})
                else:
                    others.update({organism: others[organism] + 1})
            else:
                # print('   >>>', organism)
                if len(_organism) >= 2 and re.search('^[A-Z]', organism) and re.search('^sp[.]', _organism[1]):
                    if re.search('sp[.][\s][\'][a-z]+[\']', organism):
                        _organism = '%s %s' % (_organism[0], _organism[2].replace('\'', ''))
                        save_specie(_organism, biosample, species)
                    elif re.search('sp[.]($|[\s][A-Z]+|[\s][0-9]+|[\s][\']|[\s][(]|[\s][a-z][\w]*[-]?[\w]*[0-9]+)', organism):
                        _organism = '%s %s' % (_organism[0], _organism[1])
                        save_specie(_organism, biosample, species)
                    else:
                        _organism = '%s %s' % (_organism[0], _organism[1])
                        save_specie(_organism, biosample, species)
                elif len(_organism) >= 2 and re.search('^[A-Z]', organism) and re.search('^(gen|n)[.][\s]sp[.]', ' '.join(_organism[1:])):
                    _organism = '%s sp.' % _organism[0]
                    save_specie(_organism, biosample, species)
                elif re.search('aff[.]', _organism[1]) or re.search('cf[.]', _organism[1]):
                    _organism = '%s %s' % (_organism[0], _organism[2])
                    save_specie(_organism, biosample, species)
                elif re.search(self.SPECIAL_NAME_MIXED_EST, organism):
                    _mid = _organism[1].split('/')
                    _organism1 = '%s %s' % (_organism[0], _mid[0])
                    _organism2 = '%s %s' % (_mid[1], _organism[2])

                    save_specie(_organism1, biosample, species)
                    save_specie(_organism2, biosample, species)
                else:
                    _organism = ' '.join(_organism[0:2])
                    save_specie(_organism, biosample, species)

        species_sp_genus = {}
        for _species, data in species.copy().items():
            if re.search('[\s]sp[.]$', _species):
                _count = data[self.KEY_COUNT]
                species_sp_genus.update({_species: {self.COLUMN_STATUS: False, self.KEY_COUNT: _count}})

        for _species, data in species.copy().items():
            if not re.search('[\s]sp[.]$', _species):
                _species_sp = '%s sp.' % _species.split()[0]
                if _species_sp in species_sp_genus:
                    current = species_sp_genus[_species_sp].copy()
                    current.update({self.COLUMN_STATUS: True})
                    species_sp_genus.update({_species_sp: current})

        for _species, data in species_sp_genus.items():
            _status = data[self.COLUMN_STATUS]
            _count = data[self.KEY_COUNT]
            if _status:
                del species[_species]
                species_sp_remove.update({_species: _count})

        return species, species_sp_remove, species_nom_inval, others

    def create_taxonomic_rank_file(self, dict_species, dict_species_sp, dict_nom_inval, dict_others):
        '''
            {'superkingdom': 'Eukaryota', 
             'kingdom': 'Fungi', 
             'phylum': 'Ascomycota', 
             'class': 'Saccharomycetes', 
             'order': 'Saccharomycetales', 
             'family': 'Dipodascaceae', 
             'genus': 'Yarrowia', 
             'species': 'Yarrowia lipolytica'}

            # No items found.
            [Candida] auris
            Tylopilus aff.
        '''
        dict_species_status = dict_species.copy()
        for term, data in dict_species.items():
            current = data.copy()
            current[self.COLUMN_STATUS] = None
            dict_species_status.update({term: current})

        detail_csv = {}
        if os.path.exists(self.TAXONOMIC_RANK_FILE):
            df = pd.read_csv(filepath_or_buffer = self.TAXONOMIC_RANK_FILE, sep = '\t', header = 0, index_col = False, dtype = str)
            df = df.where(pd.notnull(df), '')

            for index, row in df.iterrows():
                if not row[self.COLUMN_STATUS].startswith('[ERROR]:'):
                    _term = row[self.RANK_SPECIES]

                    current = dict_species_status[_term]
                    current[self.COLUMN_STATUS] = self.STATUS_REVIEWED
                    dict_species_status.update({_term: current})

                    _row = '\t'.join(str(value) for value in row.to_numpy())
                    detail_csv.update({index + 1: _row})

        self.show_print("Obtaining taxonomic information of %s organisms from Entrez..." % len(dict_species_status), [self.LOG_FILE])

        with open(self.TAXONOMIC_RANK_FILE, 'w', encoding = 'utf-8') as fw:
            fw.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (self.KEY_BIOSAMPLE,
                                                                                                                               self.KEY_COUNT,
                                                                                                                               self.RANK_SUPERKINGDOM,
                                                                                                                               self.RANK_KINGDOM,
                                                                                                                               self.RANK_PHYLUM,
                                                                                                                               self.RANK_CLASS,
                                                                                                                               self.RANK_ORDER,
                                                                                                                               self.RANK_FAMILY,
                                                                                                                               self.RANK_GENUS,
                                                                                                                               self.RANK_SPECIES_DB,
                                                                                                                               self.RANK_SPECIES,
                                                                                                                               self.COLUMN_STATUS,
                                                                                                                               self.IF_COLUMN_KEY,
                                                                                                                               self.IF_COLUMN_BASIONYM_KEY,
                                                                                                                               self.IF_COLUMN_RANK,
                                                                                                                               '%s (if)' % self.RANK_KINGDOM,
                                                                                                                               '%s (if)' % self.RANK_PHYLUM,
                                                                                                                               '%s (if)' % self.RANK_CLASS,
                                                                                                                               '%s (if)' % self.RANK_ORDER,
                                                                                                                               '%s (if)' % self.RANK_FAMILY,
                                                                                                                               '%s (if)' % self.RANK_GENUS,
                                                                                                                               '%s (if)' % self.RANK_SPECIES,
                                                                                                                               '%s (if)' % self.IF_COLUMN_BASIONYM,
                                                                                                                               '%s (if)' % self.COLUMN_STATUS,
                                                                                                                               self.IF_COLUMN_COMMENT))
            '''
            fw.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (self.KEY_BIOSAMPLE,
                                                                                                                       self.KEY_COUNT,
                                                                                                                       self.RANK_SUPERKINGDOM,
                                                                                                                       self.RANK_KINGDOM,
                                                                                                                       self.RANK_SUBKINGDOM,
                                                                                                                       self.RANK_SUPERPHYLUM,
                                                                                                                       self.RANK_PHYLUM,
                                                                                                                       self.RANK_SUBPHYLUM,
                                                                                                                       self.RANK_SUPERCLASS,
                                                                                                                       self.RANK_CLASS,
                                                                                                                       self.RANK_SUBCLASS,
                                                                                                                       self.RANK_SUPERORDER,
                                                                                                                       self.RANK_ORDER,
                                                                                                                       self.RANK_SUBORDER,
                                                                                                                       self.RANK_SUPERFAMILY,
                                                                                                                       self.RANK_FAMILY,
                                                                                                                       self.RANK_SUBFAMILY,
                                                                                                                       self.RANK_SUPERGENUS,
                                                                                                                       self.RANK_GENUS,
                                                                                                                       self.RANK_SUBGENUS,
                                                                                                                       self.RANK_SPECIES_DB,
                                                                                                                       self.RANK_SPECIES,
                                                                                                                       self.COLUMN_STATUS))
            '''

            with tqdm(total = len(dict_species_status)) as pbar:
                if detail_csv:
                    for _, line in detail_csv.items():
                        fw.write('%s\n' % line)
                    pbar.update(len(detail_csv))

                for term, detail in dict_species_status.items():
                    status = detail[self.COLUMN_STATUS]

                    if status is None:
                        _term = term
                        if re.search('[\s]sp[.]$', _term):
                            _term = _term.split()[0]

                        taxonomic_rank = self.get_taxonomic_rank(_term)
                        # pprint(taxonomic_rank)

                        _superkingdom = ''
                        if self.RANK_SUPERKINGDOM in taxonomic_rank:
                            _superkingdom = taxonomic_rank[self.RANK_SUPERKINGDOM]

                        _kingdom = ''
                        if self.RANK_KINGDOM in taxonomic_rank:
                            _kingdom = taxonomic_rank[self.RANK_KINGDOM]

                        _phylum = ''
                        if self.RANK_PHYLUM in taxonomic_rank:
                            _phylum = taxonomic_rank[self.RANK_PHYLUM]

                        _class = ''
                        if self.RANK_CLASS in taxonomic_rank:
                            _class = taxonomic_rank[self.RANK_CLASS]

                        _order = ''
                        if self.RANK_ORDER in taxonomic_rank:
                            _order = taxonomic_rank[self.RANK_ORDER]

                        _family = ''
                        if self.RANK_FAMILY in taxonomic_rank:
                            _family = taxonomic_rank[self.RANK_FAMILY]

                        _genus = ''
                        if self.RANK_GENUS in taxonomic_rank:
                            _genus = taxonomic_rank[self.RANK_GENUS]

                        _species = ''
                        if self.RANK_SPECIES in taxonomic_rank:
                            _species = taxonomic_rank[self.RANK_SPECIES]

                        # Plus
                        '''
                        _subkingdom = ''
                        if self.RANK_SUBKINGDOM in taxonomic_rank:
                            _subkingdom = taxonomic_rank[self.RANK_SUBKINGDOM]

                        _superphylum = ''
                        if self.RANK_SUPERPHYLUM in taxonomic_rank:
                            _superphylum = taxonomic_rank[self.RANK_SUPERPHYLUM]

                        _subphylum = ''
                        if self.RANK_SUBPHYLUM in taxonomic_rank:
                            _subphylum = taxonomic_rank[self.RANK_SUBPHYLUM]

                        _superclass = ''
                        if self.RANK_SUPERCLASS in taxonomic_rank:
                            _superclass = taxonomic_rank[self.RANK_SUPERCLASS]

                        _subclass = ''
                        if self.RANK_SUBCLASS in taxonomic_rank:
                            _subclass = taxonomic_rank[self.RANK_SUBCLASS]

                        _superorder = ''
                        if self.RANK_SUPERORDER in taxonomic_rank:
                            _superorder = taxonomic_rank[self.RANK_SUPERORDER]

                        _suborder = ''
                        if self.RANK_SUBORDER in taxonomic_rank:
                            _suborder = taxonomic_rank[self.RANK_SUBORDER]

                        _superfamily = ''
                        if self.RANK_SUPERFAMILY in taxonomic_rank:
                            _superfamily = taxonomic_rank[self.RANK_SUPERFAMILY]

                        _subfamily = ''
                        if self.RANK_SUBFAMILY in taxonomic_rank:
                            _subfamily = taxonomic_rank[self.RANK_SUBFAMILY]

                        _supergenus = ''
                        if self.RANK_SUPERGENUS in taxonomic_rank:
                            _supergenus = taxonomic_rank[self.RANK_SUPERGENUS]

                        _subgenus = ''
                        if self.RANK_SUBGENUS in taxonomic_rank:
                            _subgenus = taxonomic_rank[self.RANK_SUBGENUS]
                        '''

                        # Index Fungorum
                        taxonomic_rank_if = self.get_taxonomic_rank_if(_term)

                        fw.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (detail[self.KEY_BIOSAMPLE],
                                                                                                                                           detail[self.KEY_COUNT],
                                                                                                                                           _superkingdom,
                                                                                                                                           _kingdom,
                                                                                                                                           _phylum,
                                                                                                                                           _class,
                                                                                                                                           _order,
                                                                                                                                           _family,
                                                                                                                                           _genus,
                                                                                                                                           _species,
                                                                                                                                           term,
                                                                                                                                           taxonomic_rank[self.KEY_OUTPUTMESSAGE],
                                                                                                                                           taxonomic_rank_if[self.IF_COLUMN_KEY],
                                                                                                                                           taxonomic_rank_if[self.IF_COLUMN_BASIONYM_KEY],
                                                                                                                                           taxonomic_rank_if[self.IF_COLUMN_RANK],
                                                                                                                                           taxonomic_rank_if[self.RANK_KINGDOM],
                                                                                                                                           taxonomic_rank_if[self.RANK_PHYLUM],
                                                                                                                                           taxonomic_rank_if[self.RANK_CLASS],
                                                                                                                                           taxonomic_rank_if[self.RANK_ORDER],
                                                                                                                                           taxonomic_rank_if[self.RANK_FAMILY],
                                                                                                                                           taxonomic_rank_if[self.RANK_GENUS],
                                                                                                                                           taxonomic_rank_if[self.RANK_SPECIES],
                                                                                                                                           taxonomic_rank_if[self.IF_COLUMN_BASIONYM],
                                                                                                                                           taxonomic_rank_if[self.COLUMN_STATUS],
                                                                                                                                           taxonomic_rank_if[self.IF_COLUMN_COMMENT]))
                        '''
                        fw.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (detail[self.KEY_BIOSAMPLE],
                                                                                                                                   detail[self.KEY_COUNT],
                                                                                                                                   _superkingdom,
                                                                                                                                   _kingdom,
                                                                                                                                   _subkingdom,
                                                                                                                                   _superphylum,
                                                                                                                                   _phylum,
                                                                                                                                   _subphylum,
                                                                                                                                   _superclass,
                                                                                                                                   _class,
                                                                                                                                   _subclass,
                                                                                                                                   _superorder,
                                                                                                                                   _order,
                                                                                                                                   _suborder,
                                                                                                                                   _superfamily,
                                                                                                                                   _family,
                                                                                                                                   _subfamily,
                                                                                                                                   _supergenus,
                                                                                                                                   _genus,
                                                                                                                                   _subgenus,
                                                                                                                                   _species,
                                                                                                                                   term,
                                                                                                                                   taxonomic_rank[self.KEY_OUTPUTMESSAGE]))
                        '''
                        pbar.update(1)
        fw.close()

        with open(self.TAXONOMIC_SP_FILE, 'w', encoding = 'utf-8') as fw:
            fw.write('Term\tQuantity\n')
            for term, quant in sorted(dict_species_sp.items()):
                fw.write('%s\t%s\n' % (term, quant))
        fw.close()

        with open(self.TAXONOMIC_NOM_INVAL_FILE, 'w', encoding = 'utf-8') as fw:
            fw.write('Term\tQuantity\n')
            for term, quant in sorted(dict_nom_inval.items()):
                fw.write('%s\t%s\n' % (term, quant))
        fw.close()

        with open(self.TAXONOMIC_OTHER_NAMES_FILE, 'w', encoding = 'utf-8') as fw:
            fw.write('Term\tQuantity\n')
            for term, quant in sorted(dict_others.items()):
                fw.write('%s\t%s\n' % (term, quant))
        fw.close()
 
        self.show_print("  File with taxonomic information: %s" % self.TAXONOMIC_RANK_FILE, [self.LOG_FILE])
        self.show_print("  File with sp. names: %s" % self.TAXONOMIC_SP_FILE, [self.LOG_FILE])
        self.show_print("  File with invalid names: %s" % self.TAXONOMIC_NOM_INVAL_FILE, [self.LOG_FILE])
        self.show_print("  File with other names: %s" % self.TAXONOMIC_OTHER_NAMES_FILE, [self.LOG_FILE])
        self.show_print("", [self.LOG_FILE])

    def get_taxonomic_rank(self, term):
        # print('>>>', term)
        taxonomic_rank = {}
        try:
            Entrez.email = self.EMAIL
            handle = Entrez.esearch(term = term, db = self.DATABASE, retmode = self.FORMAT)
            record = Entrez.read(handle)
            # pprint(record)

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

                tax_id = record[0][self.KEY_TAX_ID]
                scientific_name = record[0][self.KEY_SCIENTIFIC_NAME]
                rank = record[0][self.KEY_RANK]

                if rank == self.RANK_SUPERKINGDOM:
                    taxonomic_rank.update({self.RANK_SUPERKINGDOM: scientific_name})
                elif rank == self.RANK_KINGDOM:
                    taxonomic_rank.update({self.RANK_KINGDOM: scientific_name})
                elif rank == self.RANK_PHYLUM:
                    taxonomic_rank.update({self.RANK_PHYLUM: scientific_name})
                elif rank == self.RANK_CLASS:
                    taxonomic_rank.update({self.RANK_CLASS: scientific_name})
                elif rank == self.RANK_ORDER:
                    taxonomic_rank.update({self.RANK_ORDER: scientific_name})
                elif rank == self.RANK_FAMILY:
                    taxonomic_rank.update({self.RANK_FAMILY: scientific_name})
                elif rank == self.RANK_GENUS:
                    taxonomic_rank.update({self.RANK_GENUS: scientific_name})
                elif rank == self.RANK_SPECIES:
                    taxonomic_rank.update({self.RANK_SPECIES: scientific_name})

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
                    # Plus
                    if _rank == self.RANK_SUBKINGDOM:
                        taxonomic_rank.update({self.RANK_SUBKINGDOM: _scientific_name})
                    if _rank == self.RANK_SUPERPHYLUM:
                        taxonomic_rank.update({self.RANK_SUPERPHYLUM: _scientific_name})
                    if _rank == self.RANK_SUBPHYLUM:
                        taxonomic_rank.update({self.RANK_SUBPHYLUM: _scientific_name})
                    if _rank == self.RANK_SUPERCLASS:
                        taxonomic_rank.update({self.RANK_SUPERCLASS: _scientific_name})
                    if _rank == self.RANK_SUBCLASS:
                        taxonomic_rank.update({self.RANK_SUBCLASS: _scientific_name})
                    if _rank == self.RANK_SUPERORDER:
                        taxonomic_rank.update({self.RANK_SUPERORDER: _scientific_name})
                    if _rank == self.RANK_SUBORDER:
                        taxonomic_rank.update({self.RANK_SUBORDER: _scientific_name})
                    if _rank == self.RANK_SUPERFAMILY:
                        taxonomic_rank.update({self.RANK_SUPERFAMILY: _scientific_name})
                    if _rank == self.RANK_SUBFAMILY:
                        taxonomic_rank.update({self.RANK_SUBFAMILY: _scientific_name})
                    if _rank == self.RANK_SUPERGENUS:
                        taxonomic_rank.update({self.RANK_SUPERGENUS: _scientific_name})
                    if _rank == self.RANK_SUBGENUS:
                        taxonomic_rank.update({self.RANK_SUBGENUS: _scientific_name})

            taxonomic_rank.update({self.KEY_OUTPUTMESSAGE: warning_message})
        except Exception as e:
            warning_message = '[ERROR]: %s' % e
            taxonomic_rank.update({self.KEY_OUTPUTMESSAGE: warning_message})

        return taxonomic_rank

    def get_taxonomic_rank_if(self, term):

        def get_taxonomic_rank_if_details(name_id):
            dict_species = {}
            url_searchkey = self.URL_SEARCHKEY.replace('<NAME_KEY>', name_id)
            # print(url_searchkey)
            r = requests.get(url_searchkey)
            tree = ET.fromstring(r.content)

            for element in tree:
                for child in element:
                    _text = child.text
                    if _text:
                        _text = str(_text).strip()
                    dict_species[child.tag] = _text

            return dict_species

        # Get ID
        if len(term.split()) > 1:
            _term = '+'.join(term.split())
            self.if_hits = '10'
        else:
            _term = term
            self.if_hits = '50'

        url_searchtext = self.URL_SEARCHTEXT.replace('<SPECIES_NAME>', _term).replace('<HITS>', self.if_hits)
        r = requests.get(url_searchtext)
        tree = ET.fromstring(r.content)

        msg_status = self.STATUS_OK
        dict_rank_all = {}
        for index, element in enumerate(tree):
            xml_child = {}
            for child in element:
                xml_child[child.tag] = child.text
            dict_rank_all[index] = xml_child

        # print(url_searchtext)
        # pprint(len(dict_rank_all))
        # pprint(dict_rank_all)

        # Get Detail
        basionym_name_id = None
        basionym_name = ''
        if_comment = ''
        msg_status = msg_status if dict_rank_all else 'Not found'
        dict_rank = {}
        if dict_rank_all:
            current_name_id = None
            for _, item in dict_rank_all.items():
                # pprint(item)
                if self.IF_TAG_COMMENT in item:
                    if_comment = item[self.IF_TAG_COMMENT]
                    if self.SPECIAL_NAME_NOM_ILLEGITIMATE in if_comment:
                        continue
                    if self.SPECIAL_NAME_NOM_INVALID in if_comment:
                        pass

                if self.IF_TAG_CURRENT_NAME_KEY in item:
                    current_name_id = item[self.IF_TAG_CURRENT_NAME_KEY]
                else:
                    current_name_id = item[self.IF_TAG_NAME_KEY]

                basionym_name_id = item[self.IF_TAG_BASIONYM_NAME_KEY]
                break

            if current_name_id:
                dict_rank = get_taxonomic_rank_if_details(current_name_id)
            # pprint(dict_rank)

            msg_status = msg_status if dict_rank else 'Detail not found'
            if dict_rank:
                if current_name_id != basionym_name_id:
                    dict_species_basionym = get_taxonomic_rank_if_details(basionym_name_id)

                    if self.IF_TAG_QUERY_NAME in dict_species_basionym:
                        basionym_name = dict_species_basionym[self.IF_TAG_QUERY_NAME]
                else:
                    if self.IF_TAG_QUERY_NAME in dict_rank:
                        basionym_name = dict_rank[self.IF_TAG_QUERY_NAME]

        taxonomic_rank = {self.COLUMN_STATUS: msg_status,
                          self.RANK_KINGDOM: '',
                          self.RANK_PHYLUM: '',
                          self.RANK_CLASS: '',
                          self.RANK_ORDER: '',
                          self.RANK_FAMILY: '',
                          self.RANK_GENUS: '',
                          self.RANK_SPECIES: '',
                          self.IF_COLUMN_RANK: '',
                          self.IF_COLUMN_KEY: '',
                          self.IF_COLUMN_BASIONYM_KEY: '',
                          self.IF_COLUMN_BASIONYM: '',
                          self.IF_COLUMN_COMMENT: ''}
        if dict_rank:
            if self.IF_TAG_NAME_KEY in dict_rank:
                taxonomic_rank.update({self.IF_COLUMN_KEY: dict_rank[self.IF_TAG_NAME_KEY]})
            if self.IF_TAG_KINGDOM in dict_rank:
                taxonomic_rank.update({self.RANK_KINGDOM: dict_rank[self.IF_TAG_KINGDOM]})
            if self.IF_TAG_PHYLUM in dict_rank:
                taxonomic_rank.update({self.RANK_PHYLUM: dict_rank[self.IF_TAG_PHYLUM]})
            if self.IF_TAG_CLASS in dict_rank:
                taxonomic_rank.update({self.RANK_CLASS: dict_rank[self.IF_TAG_CLASS]})
            if self.IF_TAG_ORDER in dict_rank:
                taxonomic_rank.update({self.RANK_ORDER: dict_rank[self.IF_TAG_ORDER]})
            if self.IF_TAG_FAMILY in dict_rank:
                taxonomic_rank.update({self.RANK_FAMILY: dict_rank[self.IF_TAG_FAMILY]})
            if self.IF_TAG_GENUS in dict_rank:
                taxonomic_rank.update({self.RANK_GENUS: dict_rank[self.IF_TAG_GENUS]})

            current_rank = ''
            if self.IF_TAG_RANK in dict_rank:
                current_rank = dict_rank[self.IF_TAG_RANK]

            query_name = ''
            if self.IF_TAG_QUERY_NAME in dict_rank:
                query_name = dict_rank[self.IF_TAG_QUERY_NAME]

            if current_rank == self.IF_RANK_SPECIES:
                if taxonomic_rank[self.RANK_SPECIES] == '':
                    taxonomic_rank.update({self.RANK_SPECIES: query_name})
            elif current_rank == self.IF_RANK_GENUS:
                if taxonomic_rank[self.RANK_GENUS] == '':
                    taxonomic_rank.update({self.RANK_GENUS: query_name})
            elif current_rank == self.IF_RANK_FAMILY:
                if taxonomic_rank[self.RANK_FAMILY] == '':
                    taxonomic_rank.update({self.RANK_FAMILY: query_name})
            elif current_rank == self.IF_RANK_ORDER:
                if taxonomic_rank[self.RANK_ORDER] == '':
                    taxonomic_rank.update({self.RANK_ORDER: query_name})
            elif current_rank == self.IF_RANK_CLASS:
                if taxonomic_rank[self.RANK_CLASS] == '':
                    taxonomic_rank.update({self.RANK_CLASS: query_name})
            elif current_rank == self.IF_RANK_PHYLUM:
                if taxonomic_rank[self.RANK_PHYLUM] == '':
                    taxonomic_rank.update({self.RANK_PHYLUM: query_name})
            elif current_rank == self.IF_RANK_KINGDOM:
                if taxonomic_rank[self.RANK_KINGDOM] == '':
                    taxonomic_rank.update({self.RANK_KINGDOM: query_name})

            taxonomic_rank.update({self.IF_COLUMN_BASIONYM_KEY: basionym_name_id})
            taxonomic_rank.update({self.IF_COLUMN_BASIONYM: basionym_name})
            taxonomic_rank.update({self.IF_COLUMN_COMMENT: if_comment})
            taxonomic_rank.update({self.IF_COLUMN_RANK: current_rank})

        # pprint(taxonomic_rank)

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

        with open(self.NODES_FILE, 'w', encoding = 'utf-8') as fw:
            fw.write('Id,Label,Color\n')
            for index, data in detail_label.items():
                line = '%s,%s,%s\n' % (index, data[0], data[1])
                fw.write(line)
        fw.close()

        tuples_nodes = self.get_detail_tuples()

        with open(self.EDGES_FILE, 'w', encoding = 'utf-8') as fw:
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
            _superkingdom = row[self.RANK_SUPERKINGDOM]
            _kingdom = row[self.RANK_KINGDOM]
            _phylum = row[self.RANK_PHYLUM]
            _class = row[self.RANK_CLASS]
            _order = row[self.RANK_ORDER]
            _family = row[self.RANK_FAMILY]
            _genus = row[self.RANK_GENUS]
            _species = row[self.RANK_SPECIES]
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

    def create_json(self, level):

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

        def get_count_nodes(dictionary, level_count):
            for key, data in dictionary.items():
                _sub_count = 0
                if level_count != self.RANK_GENUS:
                    for sub_key, sub_data in data.items():
                        _sub_count += sub_data['count-nodes']
                data.update({'count-nodes': len(data) + _sub_count,
                             'level': level_count})
                dictionary.update({key: data})

        def get_count_edges(dictionary, level_count):
            for key, data in dictionary.copy().items():
                data_count = data.copy()
                _count = len(data_count)
                if level_count != self.RANK_KINGDOM:
                    _count += 1
                data_count.update({'count-edges': _count,
                                   'level': level_count})
                dictionary.update({key: data_count})

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

        # For Modularity
        dict_kingdom_count = {}
        if level == self.RANK_SPECIES:
            dict_genus_count = dict_genus.copy()
            for key, data in dict_genus.copy().items():
                data_count = data.copy()
                data_count.update({'count-nodes': len(data_count),
                                   'level': self.RANK_GENUS})
                dict_genus_count.update({key: data_count})

            dict_family_count = get_tuplas(self.RANK_FAMILY, self.RANK_GENUS, dict_genus_count)
            get_count_nodes(dict_family_count, self.RANK_FAMILY)
            dict_order_count = get_tuplas(self.RANK_ORDER, self.RANK_FAMILY, dict_family_count)
            get_count_nodes(dict_order_count, self.RANK_ORDER)
            dict_class_count = get_tuplas(self.RANK_CLASS, self.RANK_ORDER, dict_order_count)
            get_count_nodes(dict_class_count, self.RANK_CLASS)
            dict_phylum_count = get_tuplas(self.RANK_PHYLUM, self.RANK_CLASS, dict_class_count)
            get_count_nodes(dict_phylum_count, self.RANK_PHYLUM)
            dict_kingdom_count = get_tuplas(self.RANK_KINGDOM, self.RANK_PHYLUM, dict_phylum_count)
            get_count_nodes(dict_kingdom_count, self.RANK_KINGDOM)

        # For k (degree) network
        dict_kingdom_edges = {}
        if level == self.RANK_SPECIES:
            dict_genus_edges = dict_genus.copy()
            get_count_edges(dict_genus_edges, self.RANK_GENUS)
            dict_family_edges = get_tuplas(self.RANK_FAMILY, self.RANK_GENUS, dict_genus_edges)
            get_count_edges(dict_family_edges, self.RANK_FAMILY)
            dict_order_edges = get_tuplas(self.RANK_ORDER, self.RANK_FAMILY, dict_family_edges)
            get_count_edges(dict_order_edges, self.RANK_ORDER)
            dict_class_edges = get_tuplas(self.RANK_CLASS, self.RANK_ORDER, dict_order_edges)
            get_count_edges(dict_class_edges, self.RANK_CLASS)
            dict_phylum_edges = get_tuplas(self.RANK_PHYLUM, self.RANK_CLASS, dict_class_edges)
            get_count_edges(dict_phylum_edges, self.RANK_PHYLUM)
            dict_kingdom_edges = get_tuplas(self.RANK_KINGDOM, self.RANK_PHYLUM, dict_phylum_edges)
            get_count_edges(dict_kingdom_edges, self.RANK_KINGDOM)

        return dict_kingdom, dict_kingdom_count, dict_kingdom_edges

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

                                if level == self.RANK_GENUS:
                                    # Without species
                                    _children5 = {'name': key_genus,
                                                  'size': len(data_species)}
                                    json_d3['children'][_index_phylum]['children'][_index_class]['children'][_index_order]['children'][_index_family]['children'].insert(_index_genus, _children5)
                                elif level == self.RANK_SPECIES:
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
        json, json_count, json_count_edges = self.create_json(level)
        json_d3_raw = self.create_json_d3(json, level)

        if level == self.RANK_SPECIES:
            json_name = 'network-species.json'
        elif level == self.RANK_GENUS:
            json_name = 'network-genus.json'

        json_d3 = os.path.join(self.OUTPUT_PATH, json_name)
        with open(json_d3, 'w', encoding = 'utf-8') as fw:
            fw.write(str(json_d3_raw).replace('\'', '"'))
        fw.close()

        # For Modularity
        if json_count:
            # pprint(json_count)
            self.get_coherence_profile(json_count)

        # For k (degree) network
        if json_count_edges:
            # pprint(json_count_edges)
            self.get_degree_stats(json_count_edges)

        self.show_print("File JSON: %s" % json_d3, [self.LOG_FILE])
        self.show_print("", [self.LOG_FILE])

    def get_coherence_profile(self, dict_json):
        network_edges = {self.RANK_KINGDOM: [],
                         self.RANK_PHYLUM: [],
                         self.RANK_CLASS: [],
                         self.RANK_ORDER: [],
                         self.RANK_FAMILY: [],
                         self.RANK_GENUS: []}
        network_n = 1
        for _kingdom, data_kingdom in dict_json.items():
            _level = data_kingdom['level']
            _count = data_kingdom['count-nodes']
            network_n += _count
            network_edges[_level].append(_count)
            for _phylum, data_phylum in data_kingdom.items():
                if _phylum not in ['level', 'count-nodes']:
                    _level = data_phylum['level']
                    _count = data_phylum['count-nodes']
                    network_edges[_level].append(_count)
                    for _class, data_class in data_phylum.items():
                        if _class not in ['level', 'count-nodes']:
                            _level = data_class['level']
                            _count = data_class['count-nodes']
                            network_edges[_level].append(_count)
                            for _order, data_order in data_class.items():
                                if _order not in ['level', 'count-nodes']:
                                    _level = data_order['level']
                                    _count = data_order['count-nodes']
                                    network_edges[_level].append(_count)
                                    for _family, data_family in data_order.items():
                                        if _family not in ['level', 'count-nodes']:
                                            _level = data_family['level']
                                            _count = data_family['count-nodes']
                                            network_edges[_level].append(_count)
                                            for _genus, data_genus in data_family.items():
                                                if _genus not in ['level', 'count-nodes']:
                                                    _level = data_genus['level']
                                                    _count = data_genus['count-nodes']
                                                    network_edges[_level].append(_count)
        # print(network_edges)
        # print(network_edges[self.RANK_PHYLUM])
        # print(network_n)

        network_coherence = {self.RANK_KINGDOM: [],
                             self.RANK_PHYLUM: [],
                             self.RANK_CLASS: [],
                             self.RANK_ORDER: [],
                             self.RANK_FAMILY: [],
                             self.RANK_GENUS: []}
        for level, array in network_edges.items():
            for edges in array:
                # Modular Coherence for each module
                coherence = 2 * edges
                coherence = coherence/(network_n * (network_n - 1))
                coherence = coherence + ((edges + 1)/network_n)
                network_coherence[level].append(coherence)
        # print(network_coherence)
        # print(network_coherence[self.RANK_PHYLUM])

        network_x = {self.RANK_KINGDOM: {'mean': 0, 'sd': 0},
                     self.RANK_PHYLUM: {'mean': 0, 'sd': 0},
                     self.RANK_CLASS: {'mean': 0, 'sd': 0},
                     self.RANK_ORDER: {'mean': 0, 'sd': 0},
                     self.RANK_FAMILY: {'mean': 0, 'sd': 0},
                     self.RANK_GENUS: {'mean': 0, 'sd': 0}}
        for level, array in network_coherence.items():
            _mean = statistics.mean(array)
            try:
                _stdev = statistics.stdev(array)
            except Exception as e:
                _stdev = 0

            network_x[level].update({'mean': round(_mean, 4)})
            network_x[level].update({'sd': round(_stdev, 4)})
        # pprint(network_x)

        coherence_profile = os.path.join(self.OUTPUT_PATH, 'coherence_profile.txt')
        with open(coherence_profile, 'w', encoding = 'utf-8') as fw:
            line = 'row\tcategory\tmean\tsd\n'
            fw.write(line)
            # line = '%s\t%s\t%s\t%s\n' % (0, self.RANK_KINGDOM, network_x[self.RANK_KINGDOM]['mean'], network_x[self.RANK_KINGDOM]['sd'])
            # fw.write(line)
            line = '%s\t%s\t%s\t%s\n' % (1, self.RANK_PHYLUM, network_x[self.RANK_PHYLUM]['mean'], network_x[self.RANK_PHYLUM]['sd'])
            fw.write(line)
            line = '%s\t%s\t%s\t%s\n' % (2, self.RANK_CLASS, network_x[self.RANK_CLASS]['mean'], network_x[self.RANK_CLASS]['sd'])
            fw.write(line)
            line = '%s\t%s\t%s\t%s\n' % (3, self.RANK_ORDER, network_x[self.RANK_ORDER]['mean'], network_x[self.RANK_ORDER]['sd'])
            fw.write(line)
            line = '%s\t%s\t%s\t%s\n' % (4, self.RANK_FAMILY, network_x[self.RANK_FAMILY]['mean'], network_x[self.RANK_FAMILY]['sd'])
            fw.write(line)
            line = '%s\t%s\t%s\t%s\n' % (5, self.RANK_GENUS, network_x[self.RANK_GENUS]['mean'], network_x[self.RANK_GENUS]['sd'])
            fw.write(line)
        fw.close()

        self.show_print("Coherence Profile File: %s" % coherence_profile, [self.LOG_FILE])
        self.show_print("", [self.LOG_FILE])

    def get_degree_stats(self, dict_json):
        network_edges = {self.RANK_KINGDOM: [],
                         self.RANK_PHYLUM: [],
                         self.RANK_CLASS: [],
                         self.RANK_ORDER: [],
                         self.RANK_FAMILY: [],
                         self.RANK_GENUS: []}
        for _kingdom, data_kingdom in dict_json.items():
            _level = data_kingdom['level']
            _count = data_kingdom['count-edges']
            network_edges[_level].append(_count)
            for _phylum, data_phylum in data_kingdom.items():
                if _phylum not in ['level', 'count-edges']:
                    _level = data_phylum['level']
                    _count = data_phylum['count-edges']
                    network_edges[_level].append(_count)
                    for _class, data_class in data_phylum.items():
                        if _class not in ['level', 'count-edges']:
                            _level = data_class['level']
                            _count = data_class['count-edges']
                            network_edges[_level].append(_count)
                            for _order, data_order in data_class.items():
                                if _order not in ['level', 'count-edges']:
                                    _level = data_order['level']
                                    _count = data_order['count-edges']
                                    network_edges[_level].append(_count)
                                    for _family, data_family in data_order.items():
                                        if _family not in ['level', 'count-edges']:
                                            _level = data_family['level']
                                            _count = data_family['count-edges']
                                            network_edges[_level].append(_count)
                                            for _genus, data_genus in data_family.items():
                                                if _genus not in ['level', 'count-edges']:
                                                    _level = data_genus['level']
                                                    _count = data_genus['count-edges']
                                                    network_edges[_level].append(_count)
        # print(network_edges)
        # print(network_edges[self.RANK_PHYLUM])

        edges_degree = []
        gamma = 2
        for level, array in network_edges.items():
            for k in array:
                edges_degree.append(k)
        edges_degree.sort(reverse = True)
        # print(edges_degree)

        k_profile = os.path.join(self.OUTPUT_PATH, 'degree_profile.txt')
        with open(k_profile, 'w', encoding = 'utf-8') as fw:
            line = 'index\tvalue\n'
            fw.write(line)
            for index, value in enumerate(edges_degree, start = 1):
                line = '%s\t%s\n' % (index, value)
                fw.write(line)
        fw.close()

        self.show_print("Degree File: %s" % k_profile, [self.LOG_FILE])
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
        oparse.TAXONOMIC_NOM_INVAL_FILE = os.path.join(oparse.OUTPUT_PATH, oparse.TAXONOMIC_NOM_INVAL_FILE)
        oparse.TAXONOMIC_OTHER_NAMES_FILE = os.path.join(oparse.OUTPUT_PATH, oparse.TAXONOMIC_OTHER_NAMES_FILE)
        oparse.TAXONOMIC_INCONSISTENCIES = os.path.join(oparse.OUTPUT_PATH, oparse.TAXONOMIC_INCONSISTENCIES)

        oparse.show_print("###########################################################", [oparse.LOG_FILE], font = oparse.BIGREEN)
        oparse.show_print("########################### RUN ###########################", [oparse.LOG_FILE], font = oparse.BIGREEN)
        oparse.show_print("###########################################################", [oparse.LOG_FILE], font = oparse.BIGREEN)

        option_extract = False
        option_gephi = True
        option_json_d3 = True

        if option_extract:
            # raw_path = 'C:\\Users\\Glen\\Dropbox\\Developmet\\network\\scripts\\ncbi-biosamples\\data'
            # raw_path = '/home/g13nj45p3r/Dropbox/Developmet/network/scripts/ncbi-biosamples/data'
            # raw_path = '/home/g13nj45p3r/Descargas/ncbi'
            # raw_path = '/media/g13nj45p3r/Data/network'
            raw_path = 'D:\\network'

            # raw_file = 'biosample_result-min.txt'
            raw_file = 'biosample_result_20211113.txt'
            # raw_file = 'biosample_result.txt'

            raw_file = os.path.join(raw_path, raw_file)

            biosamples = oparse.read_exported_file_ncbi(raw_file)
            dict_species, dict_species_sp, dict_species_nom_inval, dict_others = oparse.separate_types_classifications(biosamples)
            oparse.create_taxonomic_rank_file(dict_species, dict_species_sp, dict_species_nom_inval, dict_others)

        if option_gephi:
            # Crear los archivos .csv
            # Disponer el taxonomic_rank.csv curado
            dict_superkingdom, dict_kingdom, dict_phylum, dict_class, dict_order, dict_family, dict_genus, dict_species = oparse.get_dictionaries_index()
            oparse.create_ghepi_files(dict_superkingdom, dict_kingdom, dict_phylum, dict_class, dict_order, dict_family, dict_genus, dict_species)

        if option_json_d3:
            # Crear el fill_file y el .json
            # Disponer el taxonomic_rank.csv curado
            oparse.create_json_d3_file()
            oparse.create_json_d3_file(level = oparse.RANK_GENUS)

        oparse.show_print(oparse.finish_time(start, "Elapsed time"), [oparse.LOG_FILE])
        oparse.show_print("Done.", [oparse.LOG_FILE])
    except Exception as e:
        oparse.show_print("\n%s" % traceback.format_exc(), [oparse.LOG_FILE], font = oparse.RED)
        oparse.show_print(oparse.finish_time(start, "Elapsed time"), [oparse.LOG_FILE])
        oparse.show_print("Done.", [oparse.LOG_FILE])

if __name__ == '__main__':
    oparse = Parse()
    main(sys.argv)
