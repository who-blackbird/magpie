#!/usr/bin/env python

import re
import configparser
import argparse
import pandas as pd
import yaml
from itertools import chain

from classes.my_classes import Overlap


def main():
    """
    Parse arguments
    """
    parser = argparse.ArgumentParser(description='Calculate db overlaps')
    parser.add_argument('--config', help='Config file')
    parser.add_argument('--input', help='Input file')
    parser.add_argument('--output', help='Output file')
    parser.add_argument('--columns', help='Columns file')
    args = parser.parse_args()

    """
    Read config file
    """

    config = configparser.ConfigParser()
    config.read(args.config)

    """
    Read input file and define variables
    """
    my_dict_data = {}

    # input_file_path = config.get('DEFAULT', 'input')
    input_file_path = args.input
    input_file = open(input_file_path)
    line = input_file.readline()
    header = line.split('\t')

    overlapping = Overlap()
    overlap_threshold = float(config.get('DEFAULT', 'overlap_threshold'))  # 0.7
    bnd_window = float(config.get('DEFAULT', 'bnd_window'))  # 200

    """
    Get list of fields to be extracted
    """

    unique_dbs_str = config.get('DEFAULT', 'unique_dbs')
    unique_dbs = [x.strip() for x in unique_dbs_str.split(',')]  # ['gnomad', 'ngc38', 'wgs13k', 'ngclov']

    dbs_list = []
    for e in unique_dbs:
        e_list = list(filter(lambda x: e in x and not x.startswith('left') and not x.startswith('right')
                             , header))
        dbs_list.append(e_list)

    dbs = list(chain.from_iterable(dbs_list))  # eg: ['gnomad_s', 'gnomad_e', 'gnomad_ov_id'...]

    fields_list = [['CHROM', 'POS', 'END', 'SVTYPE'], dbs]
    fields = list(chain.from_iterable(fields_list))  # unlist

    """
    Build list of indices, since this may change. You can stick with the column names from above
    """

    index = []
    for f in fields:
        index.append(header.index(f))

    """
    Build dict elements
    """
    dic_element = 0
    for line in input_file:
        information = line.split('\t')

        my_dict_data[dic_element] = {}
        my_dict_data[dic_element]['CHROM'] = information[header.index('CHROM')]
        my_dict_data[dic_element]['POS'] = information[header.index('POS')]
        my_dict_data[dic_element]['END'] = information[header.index('END')]
        my_dict_data[dic_element]['SVTYPE'] = information[header.index('SVTYPE')]
        my_dict_data[dic_element]['dbs'] = {}

        # Let's first initialize all the dictionaries to avoid overwriting them
        for db in unique_dbs:
            if information[header.index(db + '_s')] != "":
                my_dict_data[dic_element]['dbs'][db] = {}
                my_start_split = information[header.index(db + '_s')].split(',')
                id_counter = 0
                for start in my_start_split:
                    my_dict_data[dic_element]['dbs'][db][id_counter] = {}
                    id_counter = id_counter + 1

        for db in dbs:
            if information[header.index(db)] != "":

                # split to get db name
                my_split_db_str = db.split("_")
                db_name = my_split_db_str[0]

                my_db_element = my_split_db_str[1]

                if my_db_element == 's':
                    my_start_split = information[header.index(db)].split(',')
                    id_counter = 0
                    for start in my_start_split:
                        my_dict_data[dic_element]['dbs'][db_name][id_counter].update({'START': start})
                        id_counter = id_counter + 1

                elif my_db_element == 'e':

                    my_end_split = information[header.index(db)].split(',')
                    id_counter = 0
                    for end in my_end_split:
                        my_dict_data[dic_element]['dbs'][db_name][id_counter].update({'END': end})
                        id_counter = id_counter + 1

                else:
                    # wgs13k has different ID
                    db_id = ''
                    if db_name == 'wgs13k':
                        db_id = my_split_db_str[1]
                    else:
                        db_id = my_split_db_str[2]

                    my_id_split = information[header.index(db)].split(',')

                    id_counter = 0
                    for id in my_id_split:
                        my_values_split = id.split('|')

                        db_type = ""
                        SN_value = ""

                        if db == 'gnomad_ov_id':
                            db_type = my_values_split[1]
                            SN_value = my_values_split[2]
                        else:
                            db_type = my_values_split[0]
                            SN_value = re.sub('SN=', '', my_values_split[1])

                        my_dict_data[dic_element]['dbs'][db_name][id_counter].update({'SN': SN_value})
                        my_dict_data[dic_element]['dbs'][db_name][id_counter].update({'db_type': db_type})

                        id_counter = id_counter + 1

        dic_element = dic_element + 1

    """
    Read dataframe in pandas
    """
    df_pd = pd.read_csv(input_file_path, delimiter='\t')

    """
    For each db, calculate overlaps and add to pandas dataframe
    """
    for db in unique_dbs:  # for each db

        db_results = []

        for dic_element in my_dict_data.keys():  # for each row

            dic_element_value = []

            query_start = my_dict_data[dic_element]['POS']
            query_end = my_dict_data[dic_element]['END']
            query_type = my_dict_data[dic_element]['SVTYPE']

            if query_type == 'TRA' or query_type == 'BND':
                query_start = int(query_start) - bnd_window
                query_end = int(query_start) + bnd_window

            if query_start == query_end:
                query_end = int(query_end) + 1

            # print("Query: " + str(dic_element), query_start, query_end, query_type)

            if db in my_dict_data[dic_element]['dbs'].keys():

                for id_counter in my_dict_data[dic_element]['dbs'][db].keys():
                    db_start = my_dict_data[dic_element]['dbs'][db][id_counter]['START']
                    db_end = my_dict_data[dic_element]['dbs'][db][id_counter]['END']
                    sn = my_dict_data[dic_element]['dbs'][db][id_counter]['SN']
                    db_type = my_dict_data[dic_element]['dbs'][db][id_counter]['db_type']

                    # print("\t", db, id_counter, db_start, db_end, sn, db_type)
                    if db_type == query_type:
                        MRPO = overlapping.reciprocal_percentage_overlap([query_start, query_end], [db_start, db_end])
                        # try:
                        #     MRPO = overlapping.reciprocal_percentage_overlap([query_start, query_end], [db_start, db_end])
                        # except Exception:
                        #     print("\t", query_start, query_end, id_counter, db_start, db_end, sn, db_type)

                        if MRPO >= overlap_threshold:
                            dic_element_value.append(sn)

            else:
                dic_element_value = ['0']

            db_results.append(';'.join(dic_element_value))

        db_results_ref = [x if x != '' else '0' for x in db_results]

        df_pd[db] = db_results_ref

    """
    Define columns to export
    """
    columns_path = args.columns

    with open(columns_path, 'r') as colnames:
        cols_list = yaml.safe_load(colnames)
        fixed_cols = cols_list['fixed']
        samples_cols = cols_list['samples']

    cols = list(chain.from_iterable([fixed_cols, samples_cols]))

    """
    Melt data frame
    """
    outfile = df_pd[cols]

    outfile_melt = pd.melt(outfile, id_vars=fixed_cols,
                           value_vars=samples_cols,
                           var_name='SAMPLE', value_name='GT')

    outfile_melt['GT'] = outfile_melt['GT'].str.split(":", n=1, expand=True)

    """
    Filter out ref genotypes
    """
    outfile_alt = outfile_melt.loc[-outfile_melt['GT'].isin(['0/0', './.', './0', '0/.', '0'])]

    cols_ordered = list(chain.from_iterable([['SAMPLE', 'GT'], list(fixed_cols)]))

    outfile_out = outfile_alt[cols_ordered]

    """
    Save to file
    """
    out_file_path = args.output

    outfile_out.to_csv(out_file_path, sep='\t')


if __name__ == '__main__':
    main()