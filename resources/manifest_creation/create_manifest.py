#!/usr/bin/env python3

# python standard library
import os
import collections
import argparse
import getpass

# general purpose third party packages
import numpy as np
nnz = np.count_nonzero
import pandas as pd
import petl as etl
etl.config.display_index_header = True
import MySQLdb

# Passwords
# Use usernames and passwords from the pathdev-passwords repo.
# ENA submission database (host:shap, db:subtrack)
try:
    subtrack_username = os.environ['SUBTRACK_USER']
except KeyError:
    subtrack_username = input('subtrack username: ')
try:
    subtrack_password = os.environ['SUBTRACK_PASSWORD']
except KeyError:
    subtrack_password = getpass.getpass('subtrack password: ')
# MLWH (iRODS) database (host:mlwh-db-ro, db:mlwarehouse)
try:
    mlwh_username = os.environ['MLWH_USER']
except KeyError:
    mlwh_username = input('mlwh username: ')
try:
    mlwh_password = os.environ['MLWH_PASSWORD']
except KeyError:
    mlwh_password = getpass.getpass('mlwh password: ')


def set_irods_name(row):
    # First calculate prefix using nonhuman dependent on id_run and species (taxon)
    # Then calculate tag string (empty string if no tag)
    # Then calculate file extension
    # Then concat these three strings
    if (
        (row['taxon_id'] in [7165, 7173, 30066, 62324]) & # Anopheles gambiae, arabiensis, merus and funestus respectively
        (row['id_run'] <= 6750)
    ):
        prefix = "%d_%d_nonhuman" % (row['id_run'], row['position'])
    elif (
        ~(row['taxon_id'] in [7165, 7173, 30066, 62324]) & # Anopheles gambiae, arabiensis, merus and funestus respectively
        (row['id_run'] <= 7100)
    ):
        prefix = "%d_%d_nonhuman" % (row['id_run'], row['position'])
    else:
        prefix = "%d_%d" % (row['id_run'], row['position'])
        
    if np.isnan(row['tag_index']):
        tag_string = ''
    else:
        tag_string = '#%d' % row['tag_index']
    
    if row['id_run'] <= 10300:
        file_extension = '.bam'
    else:
        file_extension = '.cram'
    
    irods_filename = prefix + tag_string + file_extension
    return irods_filename


def set_sample_id(row):
    if row['sample'] is None:
        return 'DS_%s' % row['name']
    else:
        return 'DS_%s' % row['sample']


def is_valid_alfresco_code(s):
    try: 
        int(s[0:4])
        return True
    except ValueError:
        return False


def create_build_manifest(
    taxon_ids,
    studies_to_exclude,
    miseq_runs_to_allow,
    input_dir,
    mlwh_missing_exceptions_fn,
    mlwh_study_exceptions_fn,
    mlwh_sample_exceptions_fn,
    mlwh_taxon_exceptions_fn,
    sequencescape_alfresco_study_mappings_fn,
    output_columns
):
    """Create a DataFrame to be used as a build manifest.
    
    A build manifest here is a list of all the "lanelets" that need to be
    included in a build. The output will typically be written to a
    tab-delimited file, either to use as input to vr-pipe, or perhaps as
    input to a vr-track DB.

    Args:
        taxon_ids (list of int): The taxons of the build (P. falciparum is '5833',
            P. vivax is '5855').
        sequencescape_alfresco_study_mappings_fn (str): filename of
            mappings from sequencscape to alfresco study mappings

    Returns:
        pd.DataFrame: The build manifest.

    """
    
    # Read in taxon exceptions file
    df_mlwh_taxon_exceptions = pd.read_csv(mlwh_taxon_exceptions_fn, sep='\t',
                                           dtype={'tag_index': str, 'taxon_id': int},
                                           index_col='irods_filename')
    df_mlwh_taxon_exceptions['tag_index'].fillna('', inplace=True)
    # Identify which samples in exceptions file match the taxon of current build, and which don't
    df_mlwh_exceptions_this_taxon = df_mlwh_taxon_exceptions.loc[df_mlwh_taxon_exceptions['taxon_id'].isin(taxon_ids)]
    df_mlwh_exceptions_other_taxa = df_mlwh_taxon_exceptions.loc[~(df_mlwh_taxon_exceptions['taxon_id'].isin(taxon_ids))]

    # Read in sample exceptions file
    df_mlwh_sample_exceptions = pd.read_csv(mlwh_sample_exceptions_fn, sep='\t', index_col='irods_filename')

    # Read in study exceptions file
    df_mlwh_study_exceptions = pd.read_csv(mlwh_study_exceptions_fn, sep='\t', index_col='irods_filename')

    # Read in missing exceptions file
    df_mlwh_missing_exceptions = pd.read_csv(mlwh_missing_exceptions_fn, sep='\t', index_col='irods_filename')
    df_mlwh_missing_exceptions['derivative_sample_id'] = 'DS_' + df_mlwh_missing_exceptions['sample_id']
    df_mlwh_missing_exceptions = df_mlwh_missing_exceptions.loc[
        df_mlwh_missing_exceptions['taxon_id'].isin(taxon_ids)
    ]

    # Read in SequenceScape-Alfresco study mappings
    df_sequencescape_alfresco_study_mappings = \
        pd.read_csv(sequencescape_alfresco_study_mappings_fn, sep='\t', index_col='seqscape_study_id')
    sequencescape_alfresco_study_mappings_dict = \
        df_sequencescape_alfresco_study_mappings.loc[:, 'alfresco_study_code'].to_dict()
    
    # Read in data from mlwh matching this taxon, plus samples from exceptions file matching this taxon
    conn = MySQLdb.connect(
        host='mlwh-db-ro',
        user=mlwh_username,
        password=mlwh_password,
        db='mlwarehouse',
        port=3435
    )
    
    sql_query = 'SELECT \
        study.name as study_name, \
        study.id_study_lims as study_lims, \
        sample.supplier_name as sample, \
        sample.name, \
        sample.sanger_sample_id, \
        sample.taxon_id, \
        iseq_product_metrics.id_run, \
        iseq_product_metrics.position, \
        iseq_product_metrics.tag_index, \
        iseq_product_metrics.num_reads, \
        iseq_product_metrics.human_percent_mapped, \
        iseq_run_lane_metrics.instrument_name, \
        iseq_run_lane_metrics.instrument_model, \
        iseq_run_lane_metrics.paired_read, \
        iseq_run_lane_metrics.qc_complete, \
        iseq_flowcell.manual_qc, \
        iseq_flowcell.requested_insert_size_from, \
        iseq_flowcell.requested_insert_size_to, \
        iseq_flowcell.forward_read_length, \
        iseq_run_status_dict.description \
    FROM \
        study, \
        iseq_flowcell, \
        sample, \
        iseq_product_metrics, \
        iseq_run_status, \
        iseq_run_lane_metrics, \
        iseq_run_status_dict \
    WHERE \
        study.id_study_tmp = iseq_flowcell.id_study_tmp and \
        iseq_flowcell.id_sample_tmp = sample.id_sample_tmp and \
        iseq_flowcell.manual_qc = 1 and \
        iseq_product_metrics.id_iseq_flowcell_tmp = iseq_flowcell.id_iseq_flowcell_tmp and \
        iseq_run_status.id_run = iseq_product_metrics.id_run and \
        iseq_product_metrics.id_run = iseq_run_lane_metrics.id_run and \
        iseq_product_metrics.position = iseq_run_lane_metrics.position and \
        iseq_run_status.iscurrent = 1 and \
        ( ( iseq_run_lane_metrics.instrument_model != "MiSeq" ) or (iseq_product_metrics.id_run in (%s)) ) and \
        iseq_run_status.id_run_status_dict = iseq_run_status_dict.id_run_status_dict and \
        study.faculty_sponsor = "Dominic Kwiatkowski" and \
        ( ( sample.taxon_id in (%s) ) or ' % (
            ', '.join([str(x) for x in miseq_runs_to_allow]),
            ', '.join([str(x) for x in taxon_ids])
        )
    sql_query = sql_query + ' or '.join(
        df_mlwh_taxon_exceptions.loc[df_mlwh_taxon_exceptions['taxon_id'].isin(taxon_ids)].apply(
            lambda x: '(iseq_product_metrics.id_run="%s" and \
                      iseq_product_metrics.position="%s" and \
                      iseq_product_metrics.tag_index="%s")' % (x['id_run'], x['position'], x['tag_index']),
            1
        )
    )
    sql_query = sql_query + ')'
    df_return = pd.read_sql(sql_query, conn)

    # Replace missing taxon_id with -1 (can't have missing int values in pandas Series)
    df_return['taxon_id'] = df_return['taxon_id'].fillna(-1).astype('int32')

    # Determine file name in iRods
    df_return['irods_filename'] = df_return.apply(set_irods_name, 1)
    df_return.set_index('irods_filename', inplace=True)

    # Determine alfresco study and change any incorrect studies
    df_return['alfresco_study_code'] = df_return['study_lims'].apply(
        lambda x: sequencescape_alfresco_study_mappings_dict[int(x)] if int(x) in
                  sequencescape_alfresco_study_mappings_dict else ''
    )
    empty_alfresco_study_code = (df_return['alfresco_study_code'] == '')
    df_return.loc[empty_alfresco_study_code, 'alfresco_study_code'] = \
        df_return.loc[empty_alfresco_study_code, 'study_name']
    study_exception_indexes = df_mlwh_study_exceptions.index[
        df_mlwh_study_exceptions.index.isin(df_return.index)
    ]
    df_return.loc[study_exception_indexes, 'alfresco_study_code'] = \
        df_mlwh_study_exceptions.loc[study_exception_indexes, 'alfresco_study_code']

    # Change any incorrect taxon IDs
    taxon_exception_indexes = df_mlwh_exceptions_this_taxon.index[
        df_mlwh_exceptions_this_taxon.index.to_series().isin(df_return.index)
    ]
    df_return.loc[taxon_exception_indexes, 'taxon_id'] = df_mlwh_exceptions_this_taxon.loc[taxon_exception_indexes, 'taxon_id']

    # Determine sample ID and change any incorrect sample IDs
    df_return['derivative_sample_id'] = df_return.apply(set_sample_id, 1)
    sample_exception_indexes = df_mlwh_sample_exceptions.index[
        df_mlwh_sample_exceptions.index.isin(df_return.index)
    ]    
    df_return.loc[sample_exception_indexes, 'derivative_sample_id'] = 'DS_' + df_mlwh_sample_exceptions.loc[sample_exception_indexes, 'sample_id']
    df_return.drop(['sample', 'name'], axis=1, inplace=True)
    
    # Remove any samples that are not in this taxon
    df_return = df_return.loc[~df_return.index.isin(df_mlwh_exceptions_other_taxa.index)]
    
    # Merge in missing exceptions
    df_return = df_return.append(df_mlwh_missing_exceptions.drop('study_group', 1), sort=False)
    
    # Remove any samples that don't have an alfresco study
    invalid_codes = ~df_return['alfresco_study_code'].apply(is_valid_alfresco_code)
    if np.count_nonzero(invalid_codes) > 0:
        print('Removing %d files with invalid alfresco_study_code:' % np.count_nonzero(invalid_codes))
        print(df_return['alfresco_study_code'].loc[invalid_codes].value_counts())
        df_return = df_return.loc[~invalid_codes]

    # Create other columns
    df_return.rename(columns={'alfresco_study_code': 'study', 'paired_read': 'paired'}, inplace=True)
    df_return['lane'] = df_return.index.to_series().apply(lambda x: x.split('.')[0])
    df_return['path'] = df_return.apply(lambda x: "%s/%s" % (input_dir, x.name), 1)
    df_return['reads'] = df_return['num_reads'].fillna(-1).astype(int)
    df_return['irods_path'] = df_return.apply(lambda x: "/seq/%d/%s" % (x['id_run'], x.name), 1)
    df_return['sample'] = df_return['derivative_sample_id'].str.replace('DS_', '')

    # Read in data from subtrack, most importantly to get run accessions
    # Note that here we are doing this differently from previous notebooks,
    # we are pulling in both .bam and .srf files if they exist
    # rather than setting an id_run threshold to guess which runs have .srf and which have .bam
    sql_query = "\
    SELECT \
        files.file_name as subtrack_filename, \
        files.bytes as subtrack_files_bytes, \
        files.timestamp as subtrack_files_timestamp, \
        files.public_date as subtrack_files_public_date, \
        submission.id as subtrack_submission_id, \
        submission.created as subtrack_submission_created, \
        submission.release_date as subtrack_submission_release_date, \
        submission.timestamp as subtrack_submission_timestamp, \
        submission.ext_db as subtrack_submission_ext_db, \
        submission.ebi_sample_acc, \
        submission.ebi_exp_acc, \
        submission.ebi_study_acc, \
        submission.ebi_sub_acc, \
        submission.ebi_run_acc \
    FROM \
        submission, \
        files \
    WHERE \
        files.sub_id = submission.id AND \
        (files.file_name in (%s) or files.file_name in (%s))\
    " % (
        "'" + "', '".join(df_return.index) + "'",
        "'" + "', '".join(df_return.index.to_series().apply(lambda x: x.replace('.bam', '.srf'))) + "'"
    )

    conn = MySQLdb.connect(
        host='shap',
        user=subtrack_username,
        password=subtrack_password,
        db='subtrack',
        port=3303
    )

    df_run_accessions = pd.read_sql(sql_query, conn)
    df_run_accessions['lane'] = df_run_accessions['subtrack_filename'].apply(lambda x: x.split('.')[0])
    df_run_accessions.set_index('lane', inplace=True)
    
    # Merge in subtrack data
    df_return = df_return.join(df_run_accessions, on='lane', how='left')

    # Remove unwanted studies
    df_return = df_return.loc[~df_return['study'].isin(studies_to_exclude)]

    # Remove lanes with zero reads
    df_return = df_return.loc[df_return['reads'] != 0]
    
    # Remove control samples (see emails from Sonia 03/05/2018 12:20 and from Richard 09/10/2018 13:26)
#     df_return = df_return.loc[df_return['sample'].str.slice(0, 1) != 'C']
    df_return = df_return.loc[(df_return['sample'].str.slice(0, 1) != 'C') & (df_return['sample'] != 'control')]
    
    # Restrict to certain columns
    if output_columns is not None:
        df_return = (
            df_return[output_columns]
            .sort_values(['study', 'sample', 'id_run', 'position', 'tag_index'])
        )
    
    print(df_return.shape)

    return df_return


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Script to create the pf6.4 manifest")
    parser.add_argument('-o', '--output', dest='output_fn',
                        help='path to the output manifest file [manifest.txt]', default='manifest.txt')
    parser.add_argument('-p', '--previous_imports_dir', dest='previous_imports_dir',
                        help='path to a directory containing previous manifests', default='.')
    args = parser.parse_args()

    # Inputs
    mlwh_missing_exceptions_fn = '<INSERT PATH HERE>/SIMS/meta/mlwh/mlwh_missing_exceptions.txt'
    mlwh_study_exceptions_fn = '<INSERT PATH HERE>/SIMS/meta/mlwh/mlwh_study_exceptions.txt'
    mlwh_sample_exceptions_fn = '<INSERT PATH HERE>/SIMS/meta/mlwh/mlwh_sample_exceptions.txt'
    mlwh_taxon_exceptions_fn = '<INSERT PATH HERE>/SIMS/meta/mlwh/mlwh_taxon_exceptions.txt'
    sequencescape_alfresco_study_mappings_fn = '<INSERT PATH HERE>/SIMS/meta/mlwh/sequencescape_alfresco_study_mappings.txt'

    # IDs of relevant taxa to include
    # taxon_ids=[5833, 5855, 5858, 5821, 7165, 7173, 30066], # Pf, Pv, Pm, Pb, Ag, A. arabiensis, A. merus
    taxon_ids = [5833, 36329, 5847, 57267, 137071]  # Pf, 3D7, V1, Dd2, HB3

    # Exclude studies
    # 1089 is excluded as this is R&D samples. 1204 and 1176 are excluded as these are CP1 samples. 1175 is excluded as
    # some Pf R&D samples were incorrectly tagged with this study. 1157 is a Pv study - there are two suspected Pf
    # samples in this study, but we need further investigation, and possible study change both in ROMA and in
    # Sequencescape/mlwh/ENA/iRODS before we can include
    studies_to_exclude = ['1089-R&amp;D-GB-ALCOCK', '1204-PF-GM-CP1',
                          '1176-PF-KE-CP1', '1175-VO-KH-STLAURENT',
                          '1157-PV-MULTI-PRICE']

    # Allow MiSeq runs
    # Two Miseq runs within 24 samples on each from study 1095-PF-TZ-ISHENGOMA that were included in Pf 6.0
    miseq_runs_to_allow = [13809, 13810]

    # Manifest columns
    output_columns = ['path', 'study', 'sample', 'lane', 'reads', 'paired', 'irods_path', 'sanger_sample_id',
                      'taxon_id', 'study_lims', 'study_name', 'id_run', 'position', 'tag_index', 'qc_complete',
                      'manual_qc', 'description', 'instrument_name', 'instrument_model', 'forward_read_length',
                      'requested_insert_size_from', 'requested_insert_size_to', 'human_percent_mapped',
                      'subtrack_filename', 'subtrack_files_bytes', 'ebi_run_acc']

    # Redundant variables
    input_dir = '<INSERT PATH HERE>/input'

    # Call function to create manifest
    df_assay_data = create_build_manifest(
        taxon_ids,
        studies_to_exclude,
        miseq_runs_to_allow,
        input_dir,
        mlwh_missing_exceptions_fn,
        mlwh_study_exceptions_fn,
        mlwh_sample_exceptions_fn,
        mlwh_taxon_exceptions_fn,
        sequencescape_alfresco_study_mappings_fn,
        output_columns
    )
    print(df_assay_data.head())

    # Remove samples SPT43042 (failed BQSR) and RCN12090 (failed mapping due to metadata problem)
    df_assay_data = df_assay_data.loc[~df_assay_data['sample'].isin(['SPT43042', 'RCN12090'])]


    # Sanity check on some key runs
    print()
    print('Running sanity checks...')
    # 30617 is a run for which QC has been on hold for a long time
    print(df_assay_data.loc[df_assay_data['id_run'] == 30617,
                            ['study', 'sample', 'lane', 'reads', 'study_name', 'ebi_run_acc']].shape)
    # 30854 is a run for which QC has only passed today. Should have 5 lanes
    print(df_assay_data.loc[df_assay_data['id_run'] == 30854,
                            ['study', 'sample', 'lane', 'reads', 'study_name', 'ebi_run_acc']].shape)
    # 31115 is a run for which QC has only passed today. Should have 6 lanes as lanes 2 and 4 were failed
    print(df_assay_data.loc[df_assay_data['id_run'] == 31115,
                            ['study', 'sample', 'lane', 'reads', 'study_name', 'ebi_run_acc']].shape)
    # 31240 is the run that Sonia previously wanted me to look at (see https://github.com/malariagen/parasite-ops/issues/113)
    print(df_assay_data.loc[df_assay_data['id_run'] == 31240,
                            ['study', 'sample', 'lane', 'reads', 'study_name', 'ebi_run_acc']].shape)

    print()
    print('Running additional checks...')
    print('Check which runs have missing ENA run accessions')
    print('Number of runs missing ENA run accessions:',
          df_assay_data.loc[df_assay_data['ebi_run_acc'].isnull()].groupby(['study', 'id_run']).size())
    print('Check number of runs in 1180-PF-TRAC2DONDORP with run accessions:')
    # Are there any other runs from 1180-PF-TRAC2DONDORP that have run accessions?
    df_assay_data.loc[df_assay_data['study'] == '1180-PF-TRAC2DONDORP'].groupby(['study', 'id_run']).size()
    # Are there any other runs from 1180-PF-TRAC2DONDORP that have run accessions?
    df_assay_data.loc[df_assay_data['id_run'].isin([29491, 29502, 29591])].groupby(['study', 'id_run']).size()


    # Determine files used in previous builds
    previous_imports = collections.OrderedDict()
    previous_imports['Pf6.0'] = 'pf_60_import'
    previous_imports['Pf6.1'] = 'pf_61_import'
    previous_imports['Pf6.2'] = 'pf_62_import'
    previous_imports['Pf6.3'] = 'pf_63_import'
    df_import = collections.OrderedDict()

    print('Determining samples used in previous builds:')
    for previous_import in previous_imports:
        print(previous_import)
        import_fn = f'{args.previous_imports_dir}/{previous_imports[previous_import]}.txt'
        df_import[previous_import] = pd.read_csv(import_fn, sep='\t')
        df_import[previous_import]['irods_filename'] = df_import[previous_import]['path'].apply(lambda x: x.split('/')[-1])
        df_import[previous_import].set_index('irods_filename', inplace=True)
        print(df_import[previous_import].shape)
        df_assay_data.loc[df_assay_data.index.to_series().isin(df_import[previous_import].index), 'previous_release'] = previous_import

    # Did any samples get new lanes since the previous builds?
    samples_with_new_lanes_since = collections.OrderedDict()

    print('Check whether any samples get new lanes since the previous builds:')
    for previous_import in previous_imports:
        samples_with_new_lanes_since[previous_import] = df_assay_data.loc[
            df_assay_data['previous_release'].isnull() & df_assay_data['sample'].isin(df_import[previous_import]['sample']),
            'sample'
        ]
        print(samples_with_new_lanes_since[previous_import].shape[0], "samples with new lanes since", previous_import)

    print('Check how many new samples there (since when?):')
    # How many new samples are there? Which studies are these in?
    print(df_assay_data.loc[df_assay_data['previous_release'].isnull()].groupby(['study', 'sample']).size().shape)



    # Manually fix name of study 1180-PF-TRAC2DNODORP
    # We see above that this was correctly named in Pf 6.1, but newer samples appear to have been named incorrectly.
    # I previously flagged this to Sonia, and it has been changed from 1180-PF-TRAC2DNODORP to 1180-PF-TRAC2DONDORP,
    # but to be consistent with Alfresco and previous builds, it needs to be 1180-PF-TRAC2-DONDORP.
    trac_studies = ['1044-PF-KH-FAIRHURST', '1052-PF-TRAC-WHITE', '1180-PF-TRAC2-DONDORP', '1180-PF-TRAC2DONDORP',
                    '1195-PF-TRAC2-DONDORP', '1196-PF-TRAC2-FAIRHURST']
    df_assay_data.loc[df_assay_data['study'] == '1180-PF-TRAC2DONDORP', 'study'] = '1180-PF-TRAC2-DONDORP'
    # Sanity check the above worked as expected
    print('Check that study 1180-PF-TRAC2DNODORP has been manually changed to 1180-PF-TRAC2-DONDORP:')
    print(
        df_assay_data.fillna('???').loc[df_assay_data['study'].isin(trac_studies)]
        .groupby(['study', 'previous_release', 'sample']).size().groupby(['study', 'previous_release'])
        .size()
    )

    # Ensure manifest includes only new samples, and those with new lanes
    samples_in_manifest = df_assay_data.loc[df_assay_data['previous_release'].isnull(), 'sample'].unique()
    print('Check number of new samples (since previous releases) in manifest:')
    print(samples_in_manifest.shape)
    df_irods_manifest = df_assay_data.loc[df_assay_data['sample'].isin(samples_in_manifest)]
    print(df_irods_manifest.shape)
    print('Print new samples summary:')
    print(df_irods_manifest.loc[df_irods_manifest['sample'].duplicated(keep=False),
                                ['study', 'sample', 'lane', 'reads', 'ebi_run_acc', 'previous_release']])
    print(df_irods_manifest.fillna('-').loc[df_irods_manifest['sample'].duplicated(keep=False), ].groupby(['study', 'id_run','previous_release']).size())

    # Write out file
    df_irods_manifest.to_csv(args.output_fn, sep='\t', index=False)
    print(f'Created manifest at {os.path.abspath(args.output_fn)}')

