#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import re
import sys
import json
import time
import requests
import traceback
import pandas as pd
from tqdm import tqdm
from Bio import Entrez
from urllib.parse import quote
from urllib.request import urlretrieve
from requests.exceptions import Timeout
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
        self.TAXONOMIC_RANK_FILE_RAW = 'taxonomic_rank_raw.csv'
        self.TAXONOMIC_RANK_FILE = 'taxonomic_rank.csv'
        self.TAXONOMIC_RANK_FILE_FILLED = 'taxonomic_rank_filled.csv'
        self.TAXONOMIC_RANK_FILE_VALIDATED = 'taxonomic_rank_validated.csv'
        self.TAXONOMIC_SP_FILE = 'taxonomic_sp.txt'
        self.TAXONOMIC_OTHER_NAMES_FILE = 'taxonomic_other_names.txt'
        self.TAXONOMIC_INCONSISTENCIES = 'taxonomic_inconsistencies.txt'

        self.BASE_SAMPLE_URL = 'https://api.mg-rast.org/sample/<MGS>?verbosity=full'
        self.BASE_METAGENOME_URL = 'https://api.mg-rast.org/metagenome/<MGM>?verbosity=full'
        self.BASE_GTDB_BACTERIA = 'https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/bac120_taxonomy.tsv'
        self.TIMEOUT = 120 # 2 min

        self.STATUS_OK = 'Ok'
        self.STATUS_NO_OK = 'No Ok'
        self.STATUS_ERROR = 'Error'
        self.STATUS_TIMEOUT = 'Timeout'
        self.STATUS_REVIEWED = 'Reviewed'
        self.COLUMN_STATUS = 'status'

        self.KEY_COUNT = 'count'

        # Exported file
        self.COLUMN_PROJECT_ID = 'project_id'
        self.COLUMN_SAMPLE_ID = 'sample_id'
        self.COLUMN_ORGANIZATION = 'PI_organization'
        self.COLUMN_COUNTRY = 'PI_organization_country'
        self.COLUMN_LOCATION = 'location'
        self.COLUMN_LONGITUDE = 'longitude'
        self.COLUMN_LATITUDE = 'latitude'
        self.COLUMN_COLLECTION_DATE = 'collection_date'
        self.COLUMN_CREATED_ON = 'created_on'
        self.COLUMN_ENV_PACKAGE_TYPE = 'env_package_type'
        self.COLUMN_FEATURE = 'feature'
        self.COLUMN_MATERIAL = 'material'
        self.COLUMN_INVESTIGATION_TYPE = 'investigation_type'
        self.COLUMN_SEQ_METHOD = 'seq_meth'
        self.COLUMN_SEQUENCE_TYPE = 'sequence_type'
        self.COLUMN_PUBLIC = 'public'

        self.COLUMN_METAGENOME_ID = 'metagenome_id'
        self.COLUMN_DOMAIN_GROUP = 'domain_group'
        self.COLUMN_PHYLUM_GROUP = 'phylum_group'
        self.COLUMN_CLASS_GROUP = 'class_group'
        self.COLUMN_ORDER_GROUP = 'order_group'
        self.COLUMN_FAMILY_GROUP = 'family_group'
        self.COLUMN_GENUS_GROUP = 'genus_group'
        self.COLUMN_SPECIES_GROUP = 'species_group'

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

        self.KEY_METAGENOMES = 'metagenomes'
        self.KEY_STATISTICS = 'statistics'
        self.KEY_TAXONOMY = 'taxonomy'

        self.RANK_DOMAIN = 'domain'
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

        self.RANK_SPECIES_DB = 'species (db)'

        self.SUPERKINGDOM_EUKARYOTA = 'Eukaryota'
        self.SUPERKINGDOM_BACTERIA = 'Bacteria'
        self.KINGDOM_FUNGI = 'Fungi'

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

    def read_exported_file_mgrast(self, file):
        self.show_print("Raw file: %s" % file, [self.LOG_FILE])
        self.show_print("", [self.LOG_FILE])

        df = pd.read_csv(filepath_or_buffer = file, sep = '\t', header = 1, index_col = False)
        df = df.where(pd.notnull(df), '')

        biosamples = {}
        for _, row in df.iterrows():
            _sample_id = row[self.COLUMN_SAMPLE_ID]

            if _sample_id.startswith('mgs'):
                _project_id = row[self.COLUMN_PROJECT_ID]
                _organization = row[self.COLUMN_ORGANIZATION]
                _country = row[self.COLUMN_COUNTRY]
                _location = row[self.COLUMN_LOCATION]
                _longitude = row[self.COLUMN_LONGITUDE]
                _latitude = row[self.COLUMN_LATITUDE]
                _collection_date = row[self.COLUMN_COLLECTION_DATE]
                _created_on = row[self.COLUMN_CREATED_ON]
                _env_package_type = row[self.COLUMN_ENV_PACKAGE_TYPE]
                _feature = row[self.COLUMN_FEATURE]
                _material = row[self.COLUMN_MATERIAL]
                _investigation_type = row[self.COLUMN_INVESTIGATION_TYPE]
                _seq_method = row[self.COLUMN_SEQ_METHOD]
                _sequence_type = row[self.COLUMN_SEQUENCE_TYPE]
                _public = row[self.COLUMN_PUBLIC]

                current = {}
                current[self.COLUMN_PROJECT_ID] = _project_id
                current[self.COLUMN_SAMPLE_ID] = _sample_id
                current[self.COLUMN_ORGANIZATION] = _organization
                current[self.COLUMN_COUNTRY] = _country
                current[self.COLUMN_LOCATION] = _location
                current[self.COLUMN_LONGITUDE] = _longitude
                current[self.COLUMN_LATITUDE] = _latitude
                current[self.COLUMN_COLLECTION_DATE] = _collection_date
                current[self.COLUMN_CREATED_ON] = _created_on
                current[self.COLUMN_ENV_PACKAGE_TYPE] = _env_package_type
                current[self.COLUMN_FEATURE] = _feature
                current[self.COLUMN_MATERIAL] = _material
                current[self.COLUMN_INVESTIGATION_TYPE] = _investigation_type
                current[self.COLUMN_SEQ_METHOD] = _seq_method
                current[self.COLUMN_SEQUENCE_TYPE] = _sequence_type
                current[self.COLUMN_PUBLIC] = _public
                biosamples.update({_sample_id: current})

        return biosamples

    def create_taxonomic_rank_file_raw(self, dict_biosamples):
        dict_biosamples_status = dict_biosamples.copy()
        for sample_id, sample_detail in dict_biosamples.items():
            current = sample_detail.copy()
            current[self.COLUMN_STATUS] = None
            dict_biosamples_status.update({sample_id: current})

        detail_csv = {}
        if os.path.exists(self.TAXONOMIC_RANK_FILE_RAW):
            with open(self.TAXONOMIC_RANK_FILE_RAW, 'r') as fr:
                index = 1
                for line in fr:
                    line = line.strip()
                    if not line.startswith(self.COLUMN_SAMPLE_ID):
                        _sample_id = line.split('\t')[0]

                        current = dict_biosamples_status[_sample_id]
                        current[self.COLUMN_STATUS] = self.STATUS_OK

                        dict_biosamples_status.update({_sample_id: current})
                        detail_csv.update({index: line})
                        index += 1
            fr.close()

        self.show_print("Obtaining a collection of species from %s samples..." % len(dict_biosamples_status), [self.LOG_FILE])

        with open(self.TAXONOMIC_RANK_FILE_RAW, 'w') as fw:
            fw.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (self.COLUMN_SAMPLE_ID,
                                                                                                                               self.COLUMN_PROJECT_ID,
                                                                                                                               self.COLUMN_METAGENOME_ID,
                                                                                                                               self.COLUMN_ORGANIZATION,
                                                                                                                               self.COLUMN_COUNTRY,
                                                                                                                               self.COLUMN_LOCATION,
                                                                                                                               self.COLUMN_LONGITUDE,
                                                                                                                               self.COLUMN_LATITUDE,
                                                                                                                               self.COLUMN_COLLECTION_DATE,
                                                                                                                               self.COLUMN_CREATED_ON,
                                                                                                                               self.COLUMN_ENV_PACKAGE_TYPE,
                                                                                                                               self.COLUMN_FEATURE,
                                                                                                                               self.COLUMN_MATERIAL,
                                                                                                                               self.COLUMN_INVESTIGATION_TYPE,
                                                                                                                               self.COLUMN_SEQ_METHOD,
                                                                                                                               self.COLUMN_SEQUENCE_TYPE,
                                                                                                                               self.COLUMN_PUBLIC,
                                                                                                                               self.COLUMN_DOMAIN_GROUP,
                                                                                                                               self.COLUMN_PHYLUM_GROUP,
                                                                                                                               self.COLUMN_CLASS_GROUP,
                                                                                                                               self.COLUMN_ORDER_GROUP,
                                                                                                                               self.COLUMN_FAMILY_GROUP,
                                                                                                                               self.COLUMN_GENUS_GROUP,
                                                                                                                               self.COLUMN_SPECIES_GROUP,
                                                                                                                               self.COLUMN_STATUS))

            with tqdm(total = len(dict_biosamples_status)) as pbar:
                if detail_csv:
                    for _, line in detail_csv.items():
                        fw.write('%s\n' % line)
                    pbar.update(len(detail_csv))

                for sample_id, sample_detail in dict_biosamples_status.items():
                    status = sample_detail[self.COLUMN_STATUS]

                    # print(sample_id, status)
                    if status is None:
                        _metagenome_id = self.get_metagenome_id(sample_id)

                        _domain_group = ''
                        _phylum_group = ''
                        _class_group = ''
                        _order_group = ''
                        _family_group = ''
                        _genus_group = ''
                        _species_group = ''
                        _status = self.STATUS_OK
                        if _metagenome_id:
                            if _metagenome_id == self.STATUS_TIMEOUT:
                                _status = self.STATUS_TIMEOUT
                                _metagenome_id = ''
                            else:
                                _taxonomy = self.get_species_collection(_metagenome_id)
                                # pprint(_taxonomy)

                                _status = _taxonomy[self.COLUMN_STATUS]
                                if _status == self.STATUS_OK:
                                    _domain_group = '|'.join(_taxonomy[self.RANK_DOMAIN])
                                    _phylum_group = '|'.join(_taxonomy[self.RANK_PHYLUM])
                                    _class_group = '|'.join(_taxonomy[self.RANK_CLASS])
                                    _order_group = '|'.join(_taxonomy[self.RANK_ORDER])
                                    _family_group = '|'.join(_taxonomy[self.RANK_FAMILY])
                                    _genus_group = '|'.join(_taxonomy[self.RANK_GENUS])
                                    _species_group = '|'.join(_taxonomy[self.RANK_SPECIES])

                        fw.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (sample_id,
                                                                                                                                           sample_detail[self.COLUMN_PROJECT_ID],
                                                                                                                                           _metagenome_id,
                                                                                                                                           sample_detail[self.COLUMN_ORGANIZATION],
                                                                                                                                           sample_detail[self.COLUMN_COUNTRY],
                                                                                                                                           sample_detail[self.COLUMN_LOCATION],
                                                                                                                                           sample_detail[self.COLUMN_LONGITUDE],
                                                                                                                                           sample_detail[self.COLUMN_LATITUDE],
                                                                                                                                           sample_detail[self.COLUMN_COLLECTION_DATE],
                                                                                                                                           sample_detail[self.COLUMN_CREATED_ON],
                                                                                                                                           sample_detail[self.COLUMN_ENV_PACKAGE_TYPE],
                                                                                                                                           sample_detail[self.COLUMN_FEATURE],
                                                                                                                                           sample_detail[self.COLUMN_MATERIAL],
                                                                                                                                           sample_detail[self.COLUMN_INVESTIGATION_TYPE],
                                                                                                                                           sample_detail[self.COLUMN_SEQ_METHOD],
                                                                                                                                           sample_detail[self.COLUMN_SEQUENCE_TYPE],
                                                                                                                                           sample_detail[self.COLUMN_PUBLIC],
                                                                                                                                           _domain_group,
                                                                                                                                           _phylum_group,
                                                                                                                                           _class_group,
                                                                                                                                           _order_group,
                                                                                                                                           _family_group,
                                                                                                                                           _genus_group,
                                                                                                                                           _species_group,
                                                                                                                                           _status))
                        pbar.update(1)
        fw.close()

        self.show_print("  File with species name: %s" % self.TAXONOMIC_RANK_FILE_RAW, [self.LOG_FILE])
        self.show_print("", [self.LOG_FILE])

    def get_metagenome_id(self, sample_id):
        # sample_id = 'mgs134870'
        requestURL = self.BASE_SAMPLE_URL.replace('<MGS>', sample_id)

        metagenome_id = None
        try:
            response = requests.get(requestURL, headers = {"Accept": "application/json"}, timeout = self.TIMEOUT)

            if response.ok:
                try:
                    data = json.loads(response.text)
                    metagenome_id = data[self.KEY_METAGENOMES][0][0]
                except Exception as e:
                    _search = re.search('["]mgm[0-9]+[.][0-9]+["]', response.text)
                    if _search:
                        metagenome_id = _search.group()
                        metagenome_id = metagenome_id.replace('"', '')
        except Timeout:
            metagenome_id = self.STATUS_TIMEOUT

        return metagenome_id

    def get_species_collection(self, metagenome_id):

        def append_items(tuplas):
            data = []
            for item in tuplas:
                _item = item[0]
                data.append(_item)
            return data

        # metagenome_id = 'mgm4767073.3'
        # https://api.mg-rast.org/metagenome/mgm4767073.3?verbosity=full
        # https://www.mg-rast.org/mgmain.html?mgpage=overview&metagenome=mgm4767073.3
        requestURL = self.BASE_METAGENOME_URL.replace('<MGM>', metagenome_id)

        taxonomic_data = {}
        status = self.STATUS_OK
        try:
            response = requests.get(requestURL, headers = {"Accept": "application/json"}, timeout = self.TIMEOUT)

            if response.ok:
                try:
                    data = json.loads(response.text)
                    domain_data = data[self.KEY_STATISTICS][self.KEY_TAXONOMY][self.RANK_DOMAIN]
                    phylum_data = data[self.KEY_STATISTICS][self.KEY_TAXONOMY][self.RANK_PHYLUM]
                    class_data = data[self.KEY_STATISTICS][self.KEY_TAXONOMY][self.RANK_CLASS]
                    order_data = data[self.KEY_STATISTICS][self.KEY_TAXONOMY][self.RANK_ORDER]
                    family_data = data[self.KEY_STATISTICS][self.KEY_TAXONOMY][self.RANK_FAMILY]
                    genus_data = data[self.KEY_STATISTICS][self.KEY_TAXONOMY][self.RANK_GENUS]
                    species_data = data[self.KEY_STATISTICS][self.KEY_TAXONOMY][self.RANK_SPECIES]

                    taxonomic_data[self.RANK_DOMAIN] = append_items(domain_data)
                    taxonomic_data[self.RANK_PHYLUM] = append_items(phylum_data)
                    taxonomic_data[self.RANK_CLASS] = append_items(class_data)
                    taxonomic_data[self.RANK_ORDER] = append_items(order_data)
                    taxonomic_data[self.RANK_FAMILY] = append_items(family_data)
                    taxonomic_data[self.RANK_GENUS] = append_items(genus_data)
                    taxonomic_data[self.RANK_SPECIES] = append_items(species_data)
                except Exception as e:
                    status = self.STATUS_ERROR
            else:
                status = self.STATUS_ERROR
        except Timeout:
            status = self.STATUS_TIMEOUT

        taxonomic_data[self.COLUMN_STATUS] = status

        return taxonomic_data

    def separate_types_classifications(self):

        def get_array_item(group):
            _arr = []
            for item in group.split('|'):
                if re.search('^[A-Z][a-z]+$', item.strip()):
                    _arr.append(item.strip())
            return _arr

        self.show_print("Raw taxonomy file: %s" % self.TAXONOMIC_RANK_FILE_RAW, [self.LOG_FILE])
        self.show_print("", [self.LOG_FILE])

        df = pd.read_csv(filepath_or_buffer = self.TAXONOMIC_RANK_FILE_RAW, sep = '\t', header = 0, index_col = False)
        df = df.where(pd.notnull(df), '')

        species = {}
        species_sp = {}
        others = {}
        for _, row in df.iterrows():
            _status = row[self.COLUMN_STATUS]
            _arr_species = row[self.COLUMN_SPECIES_GROUP]

            if _status == self.STATUS_OK:
                _arr_species = _arr_species.split('|')
                for organism in _arr_species:
                    '''
                        Cyanophage NATL1A-7
                        NC10 bacterium 'Dutch sediment'
                        ORF expression vector pSOX
                        uncultured virus
                        uncultured Oceanospirillales bacterium HF0500_09M11
                        uncultured Sphingobacteriales bacterium HF0130_33B19
                        uncultured ammonia-oxidizing archaeon
                        uncultured archaeon
                        uncultured bacterium
                        uncultured crenarchaeote
                        Western X phytoplasma

                        Curvibacter putative symbiont of Hydra magnipapillata
                        Frankia symbiont of Datisca glomerata
                        Wolbachia endosymbiont of Brugia malayi
                        Wolbachia endosymbiont of Cadra cautella
                        Wolbachia endosymbiont of dryinid wasps
                        Wolbachia endosymbiont of Culex pipiens pipiens
                        Wolbachia endosymbiont of Odontotermes sp. BKS-2010
                        Wolbachia endosymbiont of Apoica sp. JKS-343
                        Wolbachia endosymbiont of Chloropidae sp. JKS-375
                        Wolbachia endosymbiont of Nasutitermes sp. W7S3
                        Wolbachia endosymbiont of Drosophila takahashii subgroup
                        Wolbachia endosymbiont of Bryobia spec. V VIDR-2008
                        Rickettsia endosymbiont of Ixodes scapularis
                        Rickettsia endosymbiont of Halyzia 16-guttata orange ladybird
                        Rickettsia endosymbiont of Curculionid weevil species
                        Rickettsia endosymbiont of Aulogymnus balani/skianeuros

                        Chondrilla aff. nucula CHOND
                        Crassula aff. perforata Fishbein 3778
                        Marathrum cf. oxycarpum Cachoeirc 9/1996
                        Cynoglossus cf. arel/macrolepidotus
                        Xenopus (Silurana) tropicalis

                        Candidatus Accumulibacter phosphatis
                        Candidatus Koribacter versatilis
                        Candidatus Puniceispirillum marinum
                        Candidatus Poribacteria sp. WGA-A3

                        Plasmid pT48
                        Plasmid pTiB6S3
                        Plasmid pTiC58
                        Plasmid pTiS4
                        IncP-1beta multiresistance plasmid pB8
                        Red-recombineering helper plasmid RSFGamBet
                        Red-recombineering helper plasmid RSFGamBetkan
                        Red-recombineering helper plasmid RSFRedTER
                        Birmingham IncP-alpha plasmid

                        Vector miniTn7-bar
                        Vector mini-Tn7-Tel
                        Mutation screening vector pSPRX
                        Broad host range expression vector pRK415iq
                        Expression vector pTEV4
                        Expression vector pViet
                        Expression vector pXGWA

                        Flavobacteria bacterium BAL38
                        Flavobacteria bacterium BBFL7
                        Flavobacteria bacterium MS024-2A
                        Flavobacteria bacterium MS024-3C

                        Phage 21
                        Mycobacterium phage Cjw1
                        Mycobacterium phage Omega
                        Mycobacterium phage ScottMcG
                        Mycobacterium phage Spud

                        Vitis hybrid cultivar

                        Simian virus 40
                        Sulfolobus turreted icosahedral virus

                        Bisgaard Taxon 7
                        Cloning Vector pBANF-bar
                        Escherichia Stx1 converting bacteriophage

                        Colaspis nr. flavicornis BMNH#704425
                        Cryptoscenea nr. obscurior SLW-2008
                        Diacantha nr. fenestrata JJG232

                        Tenthredinidae gen. sp. SS-2008
                        Cyrtosiini gen. sp. MDT-2010
                        Bagoini gen. sp. DDM-2009
                        Coffea spp. mixed genomic library
                        Evexus n. sp. JRC-2004-r

                        Pauropodinae gen. sp.
                        Cricetinae gen. sp.
                    '''
                    if organism:
                        _organism = organism.split()

                        # print('>>>', organism)
                        if re.search('^Candidatus[\s]', organism):
                            _organism = _organism[1:]
                            organism = ' '.join(_organism)
                        elif re.search('[\s]of[\s]', organism):
                            _organism = organism.split(' of ')
                            organism = _organism[1]
                            _organism = organism.split()

                        # print('>>>', organism)
                        if len(_organism) >= 2 and re.search('^[A-Z]', organism) and (re.search('(spp|sp)[.]', _organism[1]) or \
                                                                                      re.search('gen[.][\s]sp[.]$', organism) or \
                                                                                      re.search('^(gen|n)[.][\s]sp[.][\s]', ' '.join(_organism[1:]))):
                            if organism not in species_sp:
                                species_sp.update({organism: 1})
                            else:
                                species_sp.update({organism: species_sp[organism] + 1})
                        elif re.search('^[a-z]', organism) or \
                             re.search('Uncultured|[(]|[\']|-like', organism) or \
                             re.search('^[A-Z][a-z]+[\s][A-Z]+', organism) or \
                             re.search('[\s]oral[\s]', organism) or \
                             (re.search('^Plasmid[\s]', organism) or re.search('[\s]plasmid[\s]', organism) or re.search('[\s]plasmid$', organism)) or \
                             (re.search('^Vector[\s]', organism) or re.search('[\s]vector[\s]', organism) or re.search('[\s]vector$', organism)) or \
                             (re.search('^Phage[\s]', organism) or re.search('[\s]phage[\s]', organism) or re.search('[\s]phage$', organism)) or \
                             (re.search('[\s]bacterium[\s]', organism) or re.search('[\s]bacterium$', organism)) or \
                             (re.search('[\s]bacteriophage[\s]', organism) or re.search('[\s]bacteriophage$', organism)) or \
                             (re.search('[\s]hybrid[\s]', organism) or re.search('[\s]hybrid$', organism)) or \
                             (re.search('[\s]virus[\s]', organism) or re.search('[\s]virus$', organism) or re.search('virus[\s]', organism)) or \
                             (len(_organism) >= 2 and (re.search('[A-Z]{2}', _organism[0]) or \
                                                       re.search('[A-Z]{2}', _organism[1]) or \
                                                       re.search('[/]', _organism[1]) or \
                                                       re.search('[0-9]', _organism[1]) or \
                                                       _organism[1].lower() == 'x')) or \
                             (len(_organism) >= 3 and re.search('[/]', _organism[2])) or \
                             len(_organism) == 1:
                            if organism not in others:
                                others.update({organism: 1})
                            else:
                                others.update({organism: others[organism] + 1})
                        else:
                            if re.search('aff[.]', _organism[1]) or \
                               re.search('cf[.]', _organism[1]) or \
                               re.search('nr[.]', _organism[1]):
                                _organism = '%s %s' % (_organism[0], _organism[2])
                            else:
                                _organism = ' '.join(_organism[0:2])

                            if _organism not in species:
                                # print('>>>', _organism)
                                _taxonomic_rank = {}
                                '''
                                _taxonomic_rank[self.RANK_SUPERKINGDOM] = get_array_item(row[self.COLUMN_DOMAIN_GROUP])
                                _taxonomic_rank[self.RANK_PHYLUM] = get_array_item(row[self.COLUMN_PHYLUM_GROUP])
                                _taxonomic_rank[self.RANK_CLASS] = get_array_item(row[self.COLUMN_CLASS_GROUP])
                                _taxonomic_rank[self.RANK_ORDER] = get_array_item(row[self.COLUMN_ORDER_GROUP])
                                _taxonomic_rank[self.RANK_FAMILY] = get_array_item(row[self.COLUMN_FAMILY_GROUP])
                                _taxonomic_rank[self.RANK_GENUS] = get_array_item(row[self.COLUMN_GENUS_GROUP])
                                '''

                                # Only the first id is registered
                                species.update({_organism: {self.KEY_COUNT: 1,
                                                            self.COLUMN_SAMPLE_ID: row[self.COLUMN_SAMPLE_ID],
                                                            self.COLUMN_PROJECT_ID: row[self.COLUMN_PROJECT_ID],
                                                            self.COLUMN_METAGENOME_ID: row[self.COLUMN_METAGENOME_ID],
                                                            self.KEY_TAXONOMY: _taxonomic_rank}})
                            else:
                                current = species[_organism].copy()
                                current[self.KEY_COUNT] = current[self.KEY_COUNT] + 1
                                species.update({_organism: current})

        return species, species_sp, others

    def create_taxonomic_rank_file(self, dict_species, dict_species_sp, dict_others):
        dict_species_status = dict_species.copy()
        for term, data in dict_species.items():
            current = data.copy()
            current[self.COLUMN_STATUS] = None
            dict_species_status.update({term: current})

        # pprint(dict_species_status)
        # exit()

        detail_csv = {}
        if os.path.exists(self.TAXONOMIC_RANK_FILE):
            df = pd.read_csv(filepath_or_buffer = self.TAXONOMIC_RANK_FILE, sep = '\t', header = 0, index_col = False)
            df = df.where(pd.notnull(df), '')

            for index, row in df.iterrows():
                _term = row[self.RANK_SPECIES]

                current = dict_species_status[_term]
                current[self.COLUMN_STATUS] = self.STATUS_REVIEWED
                dict_species_status.update({_term: current})

                _row = '\t'.join(str(value) for value in row.to_numpy())
                detail_csv.update({index + 1: _row})

        # pprint(dict_species_status)
        # pprint(detail_csv)
        # exit()

        bacteria_genus = self.get_genus_bacteria()

        self.show_print("Obtaining taxonomic information of %s organisms from Entrez..." % len(dict_species_status), [self.LOG_FILE])

        with open(self.TAXONOMIC_RANK_FILE, 'w') as fw:
            fw.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (self.COLUMN_SAMPLE_ID,
                                                                                                                               self.COLUMN_PROJECT_ID,
                                                                                                                               self.COLUMN_METAGENOME_ID,
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

            with tqdm(total = len(dict_species_status)) as pbar:
                if detail_csv:
                    for _, line in detail_csv.items():
                        fw.write('%s\n' % line)
                    pbar.update(len(detail_csv))

                for term, detail in dict_species_status.items():
                    status = detail[self.COLUMN_STATUS]

                    # print(term, status)
                    if status is None:
                        # GTDB Bacteria
                        _this_genus = term.split()[0]
                        if _this_genus in bacteria_genus:
                            taxonomic_rank = {}
                            taxonomic_rank.update({self.RANK_SUPERKINGDOM: self.SUPERKINGDOM_BACTERIA})
                            taxonomic_rank.update({self.KEY_OUTPUTMESSAGE: self.STATUS_OK})
                        else:
                            taxonomic_rank = self.get_taxonomic_rank(term)
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

                        fw.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (detail[self.COLUMN_SAMPLE_ID],
                                                                                                                                           detail[self.COLUMN_PROJECT_ID],
                                                                                                                                           detail[self.COLUMN_METAGENOME_ID],
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
            scientific_name_sp = record[0][self.KEY_SCIENTIFIC_NAME]
            taxonomic_rank.update({self.RANK_SPECIES: scientific_name_sp})

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

    # https://github.com/tqdm/tqdm/blob/master/examples/tqdm_wget.py
    def download_file_ftp_pb(self, url, filename):

        class TqdmUpTo(tqdm):
            """Provides `update_to(n)` which uses `tqdm.update(delta_n)`."""
            def update_to(self, b = 1, bsize = 1, tsize = None):
                """
                b    : int, optional
                        Number of blocks transferred so far [default: 1].
                bsize: int, optional
                        Size of each block (in tqdm units) [default: 1].
                tsize: int, optional
                        Total size (in tqdm units). If [default: None] remains unchanged.
                """
                if tsize is not None:
                    self.total = tsize
                self.update(b * bsize - self.n)  # will also set self.n = b * bsize

        try:
            with TqdmUpTo(unit = 'B', unit_scale = True, miniters = 1, desc = url.split('/')[-1]) as t:  # all optional kwargs
                urlretrieve(url, filename = filename, reporthook = t.update_to, data = None)
        except Exception as e:
            return False
        return True

    def get_genus_bacteria(self):
        filename = os.path.basename(self.BASE_GTDB_BACTERIA)
        dtdb_bacteria = os.path.join(self.OUTPUT_PATH, filename)

        is_download = True
        if not self.check_path(dtdb_bacteria):
            self.show_print("Downloading file: %s" % filename, [self.LOG_FILE])

            is_download = self.download_file_ftp_pb(self.BASE_GTDB_BACTERIA, dtdb_bacteria)

            if not is_download:
                self.show_print("ERROR, something went wrong", [self.LOG_FILE])

            self.show_print("", [self.LOG_FILE])

        bacteria_genus = {}
        if is_download:
            with open(dtdb_bacteria, 'r') as fr:
                for line in fr:
                    line = line.strip()
                    line = line.split('\t')[1]
                    _line = line.split(';')

                    for item in _line:
                        if item.startswith('g__'):
                            _genus = item.replace('g__', '')
                            if re.search('^[A-Z][a-z]+$', _genus):
                                if _genus not in bacteria_genus:
                                    bacteria_genus.update({_genus: 1})
                                else:
                                    bacteria_genus.update({_genus: bacteria_genus[_genus] + 1})
            fr.close()

        return bacteria_genus

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

    def validate_taxonomy(self):

        def get_array_item(group):
            _arr = []
            for item in group.split('|'):
                if re.search('^[A-Z][a-z]+$', item.strip()):
                    _arr.append(item.strip())
            return _arr

        self.show_print("Raw taxonomy file: %s" % self.TAXONOMIC_RANK_FILE_RAW, [self.LOG_FILE])
        self.show_print("", [self.LOG_FILE])

        df = pd.read_csv(filepath_or_buffer = self.TAXONOMIC_RANK_FILE_RAW, sep = '\t', header = 0, index_col = False)
        df = df.where(pd.notnull(df), '')

        _cond_status = df[self.COLUMN_STATUS] == self.STATUS_OK
        filtered_df = df[_cond_status]

        # print(filtered_df)
        # print(filtered_df.shape)

        taxonomy_detail = {}
        for _, row in filtered_df.iterrows():
            _taxonomic_rank = {}
            _taxonomic_rank[self.RANK_SUPERKINGDOM] = get_array_item(row[self.COLUMN_DOMAIN_GROUP])
            _taxonomic_rank[self.RANK_PHYLUM] = get_array_item(row[self.COLUMN_PHYLUM_GROUP])
            _taxonomic_rank[self.RANK_CLASS] = get_array_item(row[self.COLUMN_CLASS_GROUP])
            _taxonomic_rank[self.RANK_ORDER] = get_array_item(row[self.COLUMN_ORDER_GROUP])
            _taxonomic_rank[self.RANK_FAMILY] = get_array_item(row[self.COLUMN_FAMILY_GROUP])
            _taxonomic_rank[self.RANK_GENUS] = get_array_item(row[self.COLUMN_GENUS_GROUP])
            taxonomy_detail.update({row[self.COLUMN_SAMPLE_ID]: _taxonomic_rank})

        # pprint(taxonomy_detail)

        self.show_print("Taxonomy file: %s" % self.TAXONOMIC_RANK_FILE, [self.LOG_FILE])
        self.show_print("", [self.LOG_FILE])

        df = pd.read_csv(filepath_or_buffer = self.TAXONOMIC_RANK_FILE, sep = '\t', header = 0, index_col = False)
        df = df.where(pd.notnull(df), '')

        _cond_status = df[self.COLUMN_STATUS] == self.STATUS_OK
        _cond_superkingdom = df[self.RANK_SUPERKINGDOM] == self.SUPERKINGDOM_EUKARYOTA
        _cond_kingdom = df[self.RANK_KINGDOM] == self.KINGDOM_FUNGI

        _cond = _cond_status & _cond_superkingdom & _cond_kingdom
        filtered_df = df[_cond]

        # print(filtered_df)
        # print(filtered_df.shape)

        self.show_print("Taxonomic data validation...", [self.LOG_FILE])

        with open(self.TAXONOMIC_RANK_FILE_VALIDATED, 'w') as fw:
            fw.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (self.COLUMN_SAMPLE_ID,
                                                                                                           self.COLUMN_PROJECT_ID,
                                                                                                           self.COLUMN_METAGENOME_ID,
                                                                                                           self.KEY_COUNT,
                                                                                                           self.RANK_SUPERKINGDOM,
                                                                                                           'Val.',
                                                                                                           self.RANK_KINGDOM,
                                                                                                           self.RANK_PHYLUM,
                                                                                                           'Val.',
                                                                                                           self.RANK_CLASS,
                                                                                                           'Val.',
                                                                                                           self.RANK_ORDER,
                                                                                                           'Val.',
                                                                                                           self.RANK_FAMILY,
                                                                                                           'Val.',
                                                                                                           self.RANK_GENUS,
                                                                                                           'Val.',
                                                                                                           self.RANK_SPECIES_DB,
                                                                                                           self.RANK_SPECIES,
                                                                                                           self.COLUMN_STATUS))

            for _, row in filtered_df.iterrows():
                _sample_id = row[self.COLUMN_SAMPLE_ID]
                _species = row[self.RANK_SPECIES]

                _superkingdom_ncbi = row[self.RANK_SUPERKINGDOM]
                _phylum_ncbi = row[self.RANK_PHYLUM]
                _class_ncbi = row[self.RANK_CLASS]
                _order_ncbi = row[self.RANK_ORDER]
                _family_ncbi = row[self.RANK_FAMILY]
                _genus_ncbi = row[self.RANK_GENUS]

                _superkingdom_mgrast = taxonomy_detail[_sample_id][self.RANK_SUPERKINGDOM]
                _phylum_mgrast = taxonomy_detail[_sample_id][self.RANK_PHYLUM]
                _class_mgrast = taxonomy_detail[_sample_id][self.RANK_CLASS]
                _order_mgrast = taxonomy_detail[_sample_id][self.RANK_ORDER]
                _family_mgrast = taxonomy_detail[_sample_id][self.RANK_FAMILY]
                _genus_mgrast = taxonomy_detail[_sample_id][self.RANK_GENUS]

                _status_superkingdom = self.STATUS_NO_OK
                if _superkingdom_ncbi in _superkingdom_mgrast:
                    _status_superkingdom = self.STATUS_OK

                _status_phylum = self.STATUS_NO_OK
                if _phylum_ncbi in _phylum_mgrast:
                    _status_phylum = self.STATUS_OK

                _status_class = self.STATUS_NO_OK
                if _class_ncbi in _class_mgrast:
                    _status_class = self.STATUS_OK

                _status_order = self.STATUS_NO_OK
                if _order_ncbi in _order_mgrast:
                    _status_order = self.STATUS_OK

                _status_family = self.STATUS_NO_OK
                if _family_ncbi in _family_mgrast:
                    _status_family = self.STATUS_OK

                _status_genus = self.STATUS_NO_OK
                if _genus_ncbi in _genus_mgrast:
                    _status_genus = self.STATUS_OK

                _status = self.STATUS_NO_OK
                if _status_superkingdom == self.STATUS_OK and \
                   _status_phylum == self.STATUS_OK and \
                   _status_class == self.STATUS_OK and \
                   _status_order == self.STATUS_OK and \
                   _status_family == self.STATUS_OK and \
                   _status_genus == self.STATUS_OK:
                    _status = self.STATUS_OK

                fw.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (row[self.COLUMN_SAMPLE_ID],
                                                                                                               row[self.COLUMN_PROJECT_ID],
                                                                                                               row[self.COLUMN_METAGENOME_ID],
                                                                                                               row[self.KEY_COUNT],
                                                                                                               row[self.RANK_SUPERKINGDOM],
                                                                                                               _status_superkingdom,
                                                                                                               row[self.RANK_KINGDOM],
                                                                                                               row[self.RANK_PHYLUM],
                                                                                                               _status_phylum,
                                                                                                               row[self.RANK_CLASS],
                                                                                                               _status_class,
                                                                                                               row[self.RANK_ORDER],
                                                                                                               _status_order,
                                                                                                               row[self.RANK_FAMILY],
                                                                                                               _status_family,
                                                                                                               row[self.RANK_GENUS],
                                                                                                               _status_genus,
                                                                                                               row[self.RANK_SPECIES_DB],
                                                                                                               row[self.RANK_SPECIES],
                                                                                                               _status))
        fw.close()

        self.show_print("  File with taxonomic data validation: %s" % self.TAXONOMIC_RANK_FILE_VALIDATED, [self.LOG_FILE])
        self.show_print("", [self.LOG_FILE])

def main(args):
    try:
        start = oparse.start_time()
        oparse.create_directory(oparse.OUTPUT_PATH)
        oparse.LOG_FILE = os.path.join(oparse.OUTPUT_PATH, oparse.LOG_NAME)
        oparse.EDGES_FILE = os.path.join(oparse.OUTPUT_PATH, oparse.EDGES_FILE)
        oparse.NODES_FILE = os.path.join(oparse.OUTPUT_PATH, oparse.NODES_FILE)
        oparse.TAXONOMIC_RANK_FILE_RAW = os.path.join(oparse.OUTPUT_PATH, oparse.TAXONOMIC_RANK_FILE_RAW)
        oparse.TAXONOMIC_RANK_FILE = os.path.join(oparse.OUTPUT_PATH, oparse.TAXONOMIC_RANK_FILE)
        oparse.TAXONOMIC_RANK_FILE_FILLED = os.path.join(oparse.OUTPUT_PATH, oparse.TAXONOMIC_RANK_FILE_FILLED)
        oparse.TAXONOMIC_RANK_FILE_VALIDATED = os.path.join(oparse.OUTPUT_PATH, oparse.TAXONOMIC_RANK_FILE_VALIDATED)
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

        # oparse.validate_taxonomy()
        # exit()

        # raw_path = 'C:\\Users\\Glen\\Dropbox\\Developmet\\network\\scripts\\mg-rast\\data'
        raw_path = '/home/g13nj45p3r/Dropbox/Developmet/network/scripts/mg-rast/data'

        # Minimus
        # raw_file = 'searchResultsMetadata-min.csv'

        # All
        raw_file = 'searchResultsMetadata.csv'

        raw_file = os.path.join(raw_path, raw_file)

        biosamples = oparse.read_exported_file_mgrast(raw_file)
        # pprint(biosamples)
        oparse.create_taxonomic_rank_file_raw(biosamples)

        dict_species, dict_species_sp, dict_others = oparse.separate_types_classifications()
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
