#!/usr/bin/env python3

"""
This module consists of a modified classes and methods of DGRPy project.
@author: Farid MUSA(mmtechslv)
"""

import pandas as pd
import numpy as np
from csv import reader
from itertools import chain

__author__ = "Farid MUSA"
__copyright__ = "Copyright (C) 2019, DGRPy Project"
__credits__ = ["Farid MUSA"]
__license__ = "GPLv3"
__version__ = "1.0"
__maintainer__ = "Farid MUSA"
__email__ = "farid.musa.h@gmail.com"
__status__ = "Development"


###Following constants should not be modified
#JADVAL_ITS  links internal and external representation of values
JADVAL_ITS = {'lineage':'Consensus Lineage',
'r2Rank': {'d':'Domain','k':'Kingdom','p':'Phylum','c':'Class','o':'Order','f':'Family','g':'Genus','s':'Specie'},
'r2rank': {'d':'domain','k':'kingdom','p':'phylum','c':'class','o':'order','f':'family','g':'genus','s':'specie'},
'rank2r': {'domain':'d','kingdom':'k','phylum':'p','class':'c','order':'o','family':'f','genus':'g','specie':'s'},
'rank2Rank': {'domain':'Domain','kingdom':'Kingdom','phylum':'Phylum','class':'Class','order':'Order','family':'Family','genus':'Genus','specie':'Specie'},
'Rank2r': {'Domain':'d','Kingdom':'k','Phylum':'p','Class':'c','Order':'o','Family':'f','Genus':'g','Specie':'s'}}
JADVAL_MAIN_RANKS = ['d','k','p','c','o','f','g','s']  #JADVAL_MAIN_RANKS list of main taxonomic ranks used in Jadval. Do not modify. Order is important!


class JadvalMain:
    """Jadval is main class of DGRPy and parental class for all others. This class contains definitions, variables, constants and methods shared with other classes."""
    
    @staticmethod
    def generate_lineages_from_taxa(iTaxa,iMissingRank=False,iDesiredRanks=False,iDropRanks=False): #iMissingRank - if True will include rank preffix(such as "s__") even if rank is missing or among iDropRanks; iDesiredRanks - list of desired ranks; iDropRanks - list of undesired ranks that should be removed, this parameter is useless if iMissingRank is set to False    
        if iDesiredRanks and not all(e in JADVAL_MAIN_RANKS for e in iDesiredRanks):
            print('Impossible characters found in iDesiredRanks. Please use: '+','.join(JADVAL_MAIN_RANKS))
            return False
        iDropRanks = iDropRanks if iDropRanks else []
        add_next = (lambda r,t: ((r+'__'+t+'; ') if t != None else r+'__'+'; ')) if iMissingRank else (lambda r,t: ((r+'__'+t+'; ') if t != None else '')) 
        def make_lineage(taxon):
            tmp_lineage = ''
            for rank in make_ranks:
                tmp_lineage = tmp_lineage + add_next(rank,taxon.loc[rank] if (rank not in iDropRanks) else None)   #Add succeeding rank to lineage
            return tmp_lineage[:-2]
        if not iDesiredRanks:
            make_ranks = JADVAL_MAIN_RANKS
            new_lineages = iTaxa.apply(make_lineage,axis=1)
        else:
            make_ranks = [rank for rank in JADVAL_MAIN_RANKS if rank in iDesiredRanks]
            new_lineages = iTaxa.apply(make_lineage,axis=1)
        return new_lineages
    
    #This function find common lineages between iPrimaryLineages and iSecondaryLineages.
    #Result contains dict with {'assigned'} - reserved for lineages that are common to both lineage series
    #{'unassigned'} - reserved for lineages that are not shared between two lineage series
    #{'correlations'} - contains correlation ids between two lineages
    #Both iPrimaryLineages and iSecondaryLineages must be pandas Series with index as taxon ids and values as lineages in format that is common to both of them
    @staticmethod
    def correlate_lineages(iPrimaryLineages,iSecondaryLineages):
        s_linages = iSecondaryLineages[iSecondaryLineages.isin(iPrimaryLineages)].drop_duplicates().reset_index()
        s_linages.columns = ['secondary_id','lineage']
        s_linages = s_linages.set_index('lineage')
        p_lineages = iPrimaryLineages.reset_index()
        p_lineages.columns = ['primary_id','lineage']
        p_lineages = p_lineages.set_index('lineage')
        correlated_lineages = p_lineages.join(s_linages)
        correlations = correlated_lineages.applymap(lambda x: x if pd.notna(x) else None).set_index('primary_id')
        return {'assigned':correlated_lineages.index,'unassigned':iPrimaryLineages[~iPrimaryLineages.isin(iSecondaryLineages)],'correlations':correlations.iloc[:,0]}
    
    #Reads CSV file and returns its content
    @staticmethod
    def _read_csv(iFilename,iSeparator=',', iQuote='"'): 
        file_content= []
        with open(iFilename, 'r') as otu_file:
            otu_file_reader = reader(otu_file, delimiter=iSeparator, quotechar=iQuote)
            for row in otu_file_reader:
                file_content.append(row)
        return file_content
    
   
class JadvalTaxa(JadvalMain):
    """JadvalTaxa is a based on a DGRPy class JadvalOTU
    """
    def __init__(self,iFilename = None):
        self.taxonomy = pd.DataFrame(columns=(['lineage']+JADVAL_MAIN_RANKS)) #Dataframe with OTU taxonomy with rows as OTU ids and columns with various formatted taxonomic information
        if iFilename:
            self.load_taxa(iFilename)
        
    #Following method loads classic OTU table. Classic OTU table must have following structure. First row must contain headers. Columns: 1st for OTU ids; 2nd to last column for OTU reads, and if specified last column for taxonomy consensus lineages
    def load_taxa(self, iFilename, iSeparator=',', iQuote='"'): #iFilename - must be entered; iSeparators - characters that separates data, by default is comma; iQuote - character used for quoting in csv file, by default character is <">; iConsensusLineage - if last column contains consensus lineages for every OTU, default is True; iReversed - if table is reversed, by default OTUs rows of first column and samples are amoung columns of first row
        pre_data  = self._read_csv(iFilename,iSeparator,iQuote)
        last_col = len(pre_data[0])-1
        pre_otu_labels = [str(elem[0]) for elem in pre_data[1:]]
        pre_otu_lineages = [elem[last_col] for elem in pre_data[1:]]
        self.taxonomy = pd.DataFrame(index=pre_otu_labels,columns=self.taxonomy.columns)
        self.__construct_taxonomy_from_lineages(pd.Series(data=pre_otu_lineages,index=pre_otu_labels))
        self.reconstruct_internal_lineages()
        return
    
    #This function initiates self.taxonomy by breaking down lineages and storing ranks separately for further use
    def __construct_taxonomy_from_lineages(self,lineages):
        def allocater(lineage,ref_ranks):
            taxa_dict = {e[0]:e[1] for e in [e.strip().split('__') for e in lineage.split(';') if ('__' in e)]} #Efficiently converts lineage into dictionary
            taxa_dict_allowed = {rank:taxa_dict[rank] for rank in taxa_dict.keys() if rank in ref_ranks} #Drops unallowed ranks and orders ranks based on JADVAL_MAIN_RANKS rank order
            #Following loop sets unavailable ranks to None in order to avoid later problems
            for key in ref_ranks:
                if not (key in taxa_dict_allowed.keys()):
                    taxa_dict_allowed[key] = None
            taxa_list_ordered = [taxa_dict_allowed[rank] for rank in ref_ranks] #Orders ranks based on JADVAL_MAIN_RANKS rank order
            return [lineage]+taxa_list_ordered
        allocater_vectorized = np.vectorize(allocater,excluded=['ref_ranks'],otypes=[list])
        self.taxonomy = pd.DataFrame(index=list(lineages.index),data=list(allocater_vectorized(lineage=list(lineages.values),ref_ranks=JADVAL_MAIN_RANKS)),columns=self.taxonomy.columns)
        self.__fix_missing_taxa()
        return
    
    #This method sets taxons with "" = None. For example, lineages sometimes contain ranks such as (...; g__; s__) and (...; g__;). If consensus lineages are only available taxonomic information then such taxa are basically same and must be fixed.
    def __fix_missing_taxa(self):
        self.taxonomy.loc[:,JADVAL_MAIN_RANKS] = self.taxonomy.loc[:,JADVAL_MAIN_RANKS].applymap(lambda x: None if (x=='') else x)
        return
    
    #This method returns list of available rank levels. Method checks highest and lowest rank available in self.taxonomy that is also present among JADVAL_MAIN_RANKS
    def get_avail_ranks(self):
        return [rank for rank in JADVAL_MAIN_RANKS if self.taxonomy.loc[:,rank].notna().any()]
    
    #This method reconstructs self.taxonomy.loc[:,'lineage'] by initiated taxa
    def reconstruct_internal_lineages(self):
        self.taxonomy.loc[:,'lineage'] = self.generate_lineages(iMissingRank=True,iDesiredRanks=self.get_avail_ranks(),iDropRanks=False)
        return
    
    #This method removes OTUs that do not have any assigned taxonomy from self.taxonomy dataframe and returns OTUs that was removed
    def drop_otus_without_taxa(self):
        otu_indices_drop = self.taxonomy.loc[self.taxonomy.loc[:,JADVAL_MAIN_RANKS].agg(lambda rank: len(''.join(map(lambda x: (str(x or '')),rank))),axis=1) < 1].index
        self.taxonomy = self.taxonomy.drop(otu_indices_drop)
        return otu_indices_drop
    
    
    #This function generate desired lineages
    def generate_lineages(self,iMissingRank=False,iDesiredRanks=False,iDropRanks=False): #iMissingRank - if True will include rank preffix(such as "s__") even if rank is missing or among iDropRanks; iDesiredRanks - list of desired ranks; iDropRanks - list of undesired ranks that should be removed, this parameter is useless if iMissingRank is set to False 
        return self.generate_lineages_from_taxa(iTaxa=self.taxonomy,iMissingRank=iMissingRank,iDesiredRanks=iDesiredRanks,iDropRanks=iDropRanks)
    
    #This method loads JadvalDB and stores truncated database taxonomy and tree for further usage
    def load_database(self,iJadvalDB_obj, iRankTolerance=False):
        correlation_result = self.__derive_taxonomy_db(iJadvalDB_obj, iRankTolerance)
        if not correlation_result:
            return False
        return correlation_result
    
    #This method loads taxonomy database(JadvalDB) and correlates internal OTU ids with database taxon ids with rank tolerance given in iRank
    def __derive_taxonomy_db(self,iJadvalDB, iRankTolerance=False):
        shared_avail_ranks = [rank for rank in self.get_avail_ranks() if rank in iJadvalDB.get_avail_ranks()]
        db_lineages = iJadvalDB.generate_lineages(True,shared_avail_ranks)
        #Here rank lists are generated until it reaches iRank. For example, if iRankTolerance='f' and shared_avail_ranks=['c', 'o', 'f', 'g', 's'] then iTolRanks will contain {'f': ['g', 's'], 'g': ['s']}
        if not iRankTolerance:
           iTolRanks = False
        else:
            if (iRankTolerance in shared_avail_ranks) and ((shared_avail_ranks.index(iRankTolerance)+1)!=len(shared_avail_ranks)):
                tol_index = len(shared_avail_ranks) - shared_avail_ranks.index(iRankTolerance)
                iTolRanks = {shared_avail_ranks[-(e+1)]:shared_avail_ranks[-e:] for e in range(1,tol_index)}
            else:
                return False
        correlation_result = self.generate_taxa_correlations(db_lineages,shared_avail_ranks,iTolRanks)
        return correlation_result
        
    #This function correlate taxon ids between available lineages and iLineages. If iTolRanks is set then correlation of uncorrelated taxa will continue until all ranks are assigned or iTolRanks are exausted.
    def generate_taxa_correlations(self,iLineages,iForRanks,iTolRanks=False):
        self_lineages = self.generate_lineages(True,iForRanks)    
        correlation_results = self.correlate_lineages(self_lineages,iLineages)
        correlations = correlation_results['correlations']
        uncorrelated_lineages = correlation_results['unassigned']
        uncorrelated_otu_ids = uncorrelated_lineages.index
        correlated_lineages = correlation_results['assigned']
        if iTolRanks:
            rank_trials = {'init':{'corr':correlations.to_dict(),'uncorr_ids':uncorrelated_otu_ids,'added':sum(correlations.notna()),'lineages':correlation_results['assigned']}}
            for limiting_rank in [rank for rank in JADVAL_MAIN_RANKS[::-1] if rank in iTolRanks.keys()]:
                self_lineages = self.generate_lineages(True,iForRanks,iTolRanks[limiting_rank])
                unassigned_self_lineages = self_lineages.loc[uncorrelated_otu_ids]
                correlation_results = self.correlate_lineages(unassigned_self_lineages,iLineages)
                new_correlations = correlation_results['correlations']
                correlations.loc[uncorrelated_otu_ids] = new_correlations
                uncorrelated_lineages = correlation_results['unassigned']
                uncorrelated_otu_ids = uncorrelated_lineages.index
                correlated_lineages.append(correlation_results['assigned'])
                rank_trials[limiting_rank] = {'corr':new_correlations.to_dict(),'uncorr_ids':uncorrelated_otu_ids,'added':sum(new_correlations.notna()),'lineages':correlation_results['assigned'],'tolerated_ranks':iTolRanks[limiting_rank]}
        return {'corr_as_series':correlations,'total':sum(correlations.notna()),'rank_trials':rank_trials if iTolRanks else None}

class JadvalTaxonomy(JadvalMain):
    """JadvalTaxonomy DGRPy class contain methods and variables shared amoung JadvalTaxonomyGreengenes, JadvalTaxonomySILVA, JadvalTaxonomyNCBI, JadvalTaxonomyOTL, JadvalTaxonomyRDP """
    def __init__(self):
        super().__init__()
        self.taxonomy_df = pd.DataFrame(columns=(['lineage']+JADVAL_MAIN_RANKS)) #Dataframe with OTU taxonomy with rows as OTU ids and columns with various formatted taxonomic information

    
    #This function generate desired lineages
    def generate_lineages(self,iMissingRank=False,iDesiredRanks=False,iDropRanks=False): #iMissingRank - if True will include rank preffix(such as "s__") even if rank is missing or among iDropRanks; iDesiredRanks - list of desired ranks; iDropRanks - list of undesired ranks that should be removed, this parameter is useless if iMissingRank is set to False 
        return self.generate_lineages_from_taxa(iTaxa=self.taxonomy_df,iMissingRank=iMissingRank,iDesiredRanks=iDesiredRanks,iDropRanks=iDropRanks)
    
    #This method returns list of available rank levels. Method checks highest and lowest rank available in self.taxonomy_df that is also present among JADVAL_MAIN_RANKS
    def get_avail_ranks(self):
        return [rank for rank in JADVAL_MAIN_RANKS if self.taxonomy_df.loc[:,rank].notna().any()]
    
    #This method sets taxons with "" = None. For example, lineages sometimes contain ranks such as (...; g__; s__) and (...; g__;). If consensus lineages are only available taxonomic information then such taxa are basically same and must be fixed.
    def _fix_missing_taxa(self):
        self.taxonomy_df.loc[:,JADVAL_MAIN_RANKS] = self.taxonomy_df.loc[:,JADVAL_MAIN_RANKS].applymap(lambda x: None if (x=='') else x)
        return
    
    #This method reconstructs self.taxonomy_df.loc[:,'lineage'] by initiated taxa
    def reconstruct_internal_lineages(self):
        self.taxonomy_df.loc[:,'lineage'] = self.generate_lineages(iMissingRank=True,iDesiredRanks=self.get_avail_ranks(),iDropRanks=False)
        return 

class JadvalTaxonomyGreengenes(JadvalTaxonomy):
    """JadvalTaxonomyGreengenes DGRPy class that is responsible for loading, maintaining, and manipulating GreenGenes Taxonomy"""
    def __init__(self,iTaxaMapFile=None):
        super().__init__()
        if iTaxaMapFile:
            self.load_taxonomy_database(iTaxaMapFile)
    
    #Following method loads taxonomic database
    def load_taxonomy_database(self,iTaxaMapFile):       
        self.__load_taxonomy_map(iTaxaMapFile)
        self.__init_internal_taxonomy()
        self.reconstruct_internal_lineages()
        return
    
    #Following method loads greengeens taxonomic database
    def __load_taxonomy_map(self,iTaxaMapFile): #iDatabaseName - Name of database; iTreeFile - Tree file in format given in iTreeFormat, by default is newick; iTaxaMapFile - file with ID - Taxonomy data, by default values are separated by \t
        tmp_taxa_map = []
        with open(iTaxaMapFile, 'r') as map_file:
            taxa_map_file = reader(map_file, delimiter='\t')
            for taxa in taxa_map_file:
                tmp_taxa_map.append(taxa)
        self.taxonomy_map = pd.Series(data=[e[1] for e in tmp_taxa_map],index=[e[0] for e in tmp_taxa_map])
        return
    
    #This function initiates self.taxonomy_df by breaking down lineages and storing ranks separately (Vectorization is not an option. Reoptimization is necessary!)
    def __init_internal_taxonomy(self):
        def allocater(lineage,ref_ranks):
            taxa_dict = {e[0]:e[1] for e in [e.strip().split('__') for e in lineage.split(';') if ('__' in e)]} #Efficiently converts lineage into dictionary
            taxa_dict_allowed = {rank:taxa_dict[rank] for rank in taxa_dict.keys() if rank in ref_ranks} #Drops unallowed ranks and orders ranks based on JADVAL_MAIN_RANKS rank order
            #Following loop sets unavailable ranks to None in order to avoid later problems
            for key in ref_ranks:
                if not (key in taxa_dict_allowed.keys()):
                    taxa_dict_allowed[key] = None
            taxa_list_ordered = [taxa_dict_allowed[rank] for rank in ref_ranks] #Orders ranks based on JADVAL_MAIN_RANKS rank order
            return [lineage]+taxa_list_ordered
        allocater_vectorized = np.vectorize(allocater,excluded=['ref_ranks'],otypes=[list])
        self.taxonomy_df = pd.DataFrame(index=list(self.taxonomy_map.index),data=list(allocater_vectorized(lineage=list(self.taxonomy_map.values),ref_ranks=JADVAL_MAIN_RANKS)),columns=self.taxonomy_df.columns)
        self._fix_missing_taxa()
        return True

class JadvalTaxonomySILVA(JadvalTaxonomy):
    """JadvalTaxonomySILVA DGRPy class that is responsible for loading, maintaining, and manipulating SILVA Taxonomy"""
    def __init__(self,iTaxaMapFile=None):
        super().__init__()
        if iTaxaMapFile:
            self.load_taxonomy_database(iTaxaMapFile)
    
    #Following method loads taxonomic database
    def load_taxonomy_database(self, iTaxaMapFile):       
        self.__load_taxonomy_map(iTaxaMapFile)
      
        self.__init_internal_taxonomy()   
        self.reconstruct_internal_lineages()
        return
    
    #Following method loads greengeens taxonomic database
    def __load_taxonomy_map(self,iTaxaMapFile): #iTaxaMapFile - file with ID - Taxonomy data, by default values are separated by \t
        tmp_taxa_map = []
        with open(iTaxaMapFile, 'r') as map_file:
            taxa_map_file = reader(map_file, delimiter='\t')
            for taxon in taxa_map_file:
                tmp_taxa_map.append(taxon)
        self.taxonomy_map = pd.DataFrame(data=[[e[0], e[2]] for e in tmp_taxa_map],index=[e[1] for e in tmp_taxa_map],columns=['lineage','level'])
        return
    
    #This function initiates self.taxonomy_df by breaking down lineages and storing ranks separately (This function is NOT OPTIMIZED. Cython can by considered)
    def __init_internal_taxonomy(self):
        SILVA_RANKS = JADVAL_MAIN_RANKS[:-1]
        tmp_taxonomy_df = pd.DataFrame(columns=(['lineage']+SILVA_RANKS),index=self.taxonomy_map.index,data=None)
        rank = SILVA_RANKS[0]
        tmp_taxonomy_df.loc[:,rank] = self.taxonomy_map.loc[self.taxonomy_map.loc[:,'level']==JADVAL_ITS['r2rank'][rank],:].apply(lambda row: row['lineage'][:-1].split(';')[0],axis=1)
        for rank in SILVA_RANKS[1:]:    
            tmp_taxonomy_df.loc[:,rank] = self.taxonomy_map.loc[self.taxonomy_map.loc[:,'level']==JADVAL_ITS['r2rank'][rank],:].apply(lambda row: row['lineage'][:-1].split(';'),axis=1)
        def reassign(taxons,r):
            new_taxa = {r:taxons.pop()}
            for r_i in SILVA_RANKS[:SILVA_RANKS.index(r)]:
                for t_i in range(len(taxons)):
                    if taxons[t_i] in tmp_taxonomy_df.loc[tmp_taxonomy_df.loc[:,r_i].notna(),r_i].values:
                        new_taxa[r_i] = taxons[t_i]
                        taxons.pop(t_i)
                        break;
            return pd.Series(new_taxa)
        for rank in SILVA_RANKS[1:]:    
            tmp_taxonomy_df.update(tmp_taxonomy_df.loc[tmp_taxonomy_df.loc[:,rank].notna(),rank].apply(reassign,r=rank))
        tmp_taxonomy_df.insert(loc=len(tmp_taxonomy_df.columns),column='s',value=None)
        tmp_taxonomy_df = tmp_taxonomy_df.applymap(lambda x: None if pd.isna(x) else x)
        self.taxonomy_df = tmp_taxonomy_df
        self._fix_missing_taxa()
        return
    
