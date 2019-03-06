from . import _tm_main

def error_handler(error_code):
    if error_code in _tm_main.error_codes():
        return _tm_main.error_echo(error_code)
    else:
        return 'Unknown Error. ERROR_CODE: %s' %error_code

def assert_inputs(db_type,ref_path,target_path,output_path):
    error = False
    error = _tm_main.check_avail_db_type(db_type)
    if _tm_main.error_pass(error):
        if db_type == 'greengenes':
            error = _tm_main.assert_greengenes_taxa_map(ref_path)
        elif db_type == 'silva':
            error = _tm_main.assert_silva_taxa_map(ref_path)
        
        if _tm_main.error_pass(error):
            error = _tm_main.assert_lineages(target_path)
            if _tm_main.error_pass(error):
                error = _tm_main.os.path.isfile(output_path)
                if _tm_main.error_pass(error):
                    return True
    return error

def evaluate_correlations(db_type,ref_path,target_path,output_path):
    echo_msgs = False
    error = False
    assert_pass = assert_inputs(db_type,ref_path,target_path,output_path)
    if assert_pass:
        result = _tm_main.lineage_correlator(db_type,ref_path,target_path)
        if _tm_main.error_pass(result):
            echo_msgs = []
            correlations = result['correlations']
            echo_msgs.append('Total correlations: %d' %result['total_assignments'])
            if len(result['removed_taxa'])>0:
                echo_msgs.append('Missing taxa IDs: [%s]' %(','.join(result['removed_taxa'])))
            file_gen = _tm_main.csv_output(output_path,db_type,correlations.loc[correlations.notna()])
            if file_gen == True:
                echo_msgs.append('Generating %s' %output_path)
            else:
                error = file_gen
        else:
            error = result
    else:
        error = assert_pass
        
    return echo_msgs, error

def get_greengenes(ref_path=False,target_path=False,output_path=False):
    return evaluate_correlations('greengenes',ref_path,target_path,output_path)

def get_silva(ref_path=False,target_path=False,output_path=False):
    return evaluate_correlations('silva',ref_path,target_path,output_path)