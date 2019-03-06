from . import _jadval as jadval
from . import _tm_constants as TMC
import csv
import os


def error_pass(code):
    if type(code)==str:
        if code not in TMC.ERROR_CODES.keys():
            return True
        else:
            return False
    elif type(code)==type:
        if issubclass(code,BaseException):
            return False
        else:
            return True
    else:
        return True
    return True

def error_codes():
    return TMC.ERROR_CODES.keys()

def error_echo(code):
    return TMC.ERROR_CODES[code] if code in TMC.ERROR_CODES.keys() else False

def check_avail_db_type(db_type):
    return True if db_type in TMC.AVAIL_DATABASES else 'TPE-1'

def assert_greengenes_taxa_map(map_file):
    check = os.path.isfile(map_file) and os.stat(map_file).st_size > 0
    return True if check else 'TPE-3'

def assert_silva_taxa_map(map_file):
    check = os.path.isfile(map_file) and os.stat(map_file).st_size > 0
    return True if check else 'TPE-4'

def assert_lineages(lineage_file):
    check = os.path.isfile(lineage_file) and os.stat(lineage_file).st_size > 0
    return True if check else 'TPE-5'

def read_db(database_type,map_file):
    error_code = False
    if database_type == 'greengenes':
        file_pass = assert_greengenes_taxa_map(map_file)
        if error_pass(file_pass):
            return jadval.Greengenes(map_file) #Loading Greengeens database
        else:
            error_code = file_pass
    elif database_type == 'silva':
        file_pass = assert_silva_taxa_map(map_file)
        if error_pass(file_pass):
            return jadval.SILVA(map_file) #Loading Greengeens database
        else:
            error_code = file_pass
    else:
        error_code = 'TPE-2'
    return error_code
                
def read_lineages(filename):
    error_code = False
    file_pass = assert_lineages(filename)
    if error_pass(file_pass):
        jadval_with_lineages = jadval.Taxa(filename)
        removed_taxa = jadval_with_lineages.drop_otus_without_taxa()
        return {'lineages':jadval_with_lineages,'removed':removed_taxa}
    else:
        error_code = file_pass
    return error_code

def lineage_correlator(db_type,map_file,lineage_file):
    error_code = False
    db_avail = check_avail_db_type(db_type)
    if error_pass(db_avail):
        ref_taxa = read_db(db_type,map_file)
        if error_pass(ref_taxa):
            target_taxa = read_lineages(lineage_file)
            if error_pass(target_taxa):
                ids_without_taxa = list(target_taxa['removed'])
                target_taxa = target_taxa['lineages']
                correlations = target_taxa.load_database(ref_taxa)
                return {'removed_taxa':ids_without_taxa,'correlations':correlations['corr_as_series'],'total_assignments':correlations['total']}
            else:
                error_code = target_taxa
        else:
            error_code = ref_taxa
    else:
        error_code = db_avail
    return error_code

def csv_output(filename,ref_type,correlations):
    corr_data = [[target_id,ref_id] for target_id, ref_id in correlations.items()]
    corr_data.sort(key=lambda x:x[0])
    headers = ['Target ID',ref_type+' ID']
    try:
        with open(filename, 'w') as out_file:
            csv_writer = csv.writer(out_file)
            csv_writer.writerow(headers)
            csv_writer.writerows(corr_data)
    except Exception as e:
        return e
    return
    
