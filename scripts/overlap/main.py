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
    parser.add_argument('--verbose', help='increase output verbosity', action='store_true')
    args = parser.parse_args()

    """
    Read config file
    """
    if args.verbose:
        print("Reading config file")
    config = configparser.ConfigParser()
    config.read(args.config)

    """
    Read input file and define variables
    """
    if args.verbose:
        print("Reading input file and defining variables")

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
    if args.verbose:
        print("Getting list of fields to be extracted")

    unique_dbs_str = config.get('DEFAULT', 'unique_dbs')
    unique_dbs = [x.strip() for x in unique_dbs_str.split(',')]  # ['gnomad', 'ngc38', 'wgs13k', 'ngclov', 'bkl']

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
    if args.verbose:
        print("Building list of indices")

    index = []
    for f in fields:
        index.append(header.index(f))

    """
    Build dict elements
    """
    if args.verbose:
        print("Building dictionary elements")

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
                my_split = information[header.index(db)].split(',')

                id_counter = 0

                if my_db_element == 's':
                    for start in my_split:
                        my_dict_data[dic_element]['dbs'][db_name][id_counter].update({'START': start})
                        id_counter = id_counter + 1

                elif my_db_element == 'e':
                    for end in my_split:
                        my_dict_data[dic_element]['dbs'][db_name][id_counter].update({'END': end})
                        id_counter = id_counter + 1

                else:
                    # wgs13k has different ID
                    # db_id = ''
                    # if db_name in ['wgs13k', 'bkl']:
                    # if db_name in ['wgs13k']:
                    #     db_id = my_split_db_str[1]
                    # else:
                    #     db_id = my_split_db_str[2]

                    # my_id_split = information[header.index(db)].split(',')

                    for id in my_split:
                        my_values_split = id.split('|')

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
    if args.verbose:
        print("Reading data frame in pandas")

    # df_pd = pd.read_csv(input_file_path, delimiter='\t')

    tmplist = []
    # pool = mp.Pool(4)  # TO DO multiprocessing

    for chunk in pd.read_csv(input_file_path, delimiter='\t', chunksize=2000):
        # c = pool.apply_async(process_frame, [chunk])
        tmplist.append(chunk)

    df_pd = pd.concat(tmplist, axis=0)
    del tmplist


    """
    For each db, calculate overlaps and add to pandas dataframe
    """
    if args.verbose:
        print("Calculating overlaps")

    for db in unique_dbs:  # for each db

        db_results = []

        for dic_element in my_dict_data.keys():  # for each row

            dic_element_value = []

            query_start = my_dict_data[dic_element]['POS']
            query_end = my_dict_data[dic_element]['END']
            query_type = my_dict_data[dic_element]['SVTYPE']

            if query_type == 'TRA' or query_type == 'INS':
                query_start = int(query_start) - bnd_window
                query_end = int(query_start) + bnd_window

            if query_start == query_end:
                query_end = int(query_end) + 1

            # print("Query: " + str(dic_element), query_start, query_end, query_type)

            if db in my_dict_data[dic_element]['dbs'].keys():

                for id_counter in my_dict_data[dic_element]['dbs'][db].keys():
                    db_start = my_dict_data[dic_element]['dbs'][db][id_counter]['START']
                    db_end = my_dict_data[dic_element]['dbs'][db][id_counter]['END']

                    if db == 'bkl':
                        try:
                            ovl = overlapping.overlap([int(query_start), int(query_end)],
                                                      [int(db_start), int(db_end)])
                            db_size = (int(ovl[1]) - int(ovl[0])) + 1
                        except Exception:
                            db_size = 0

                        dic_element_value.append(db_size)

                    else:
                        sn = my_dict_data[dic_element]['dbs'][db][id_counter]['SN']
                        db_type = my_dict_data[dic_element]['dbs'][db][id_counter]['db_type']

                        if db_type == query_type:
                            MRPO = overlapping.reciprocal_percentage_overlap([int(query_start), int(query_end)],
                                                                             [int(db_start), int(db_end)])

                            if query_type == 'TRA' or query_type == 'INS':
                                dic_element_value.append(sn)
                            else:
                                if MRPO >= overlap_threshold:
                                    dic_element_value.append(sn)

            else:
                dic_element_value = ['0']

            # print("Query: " + str(dic_element_value), query_start, query_end, query_type)

            ##If bkl sum all values, else return comma separated
            if db == 'bkl':
                bkl_lst = list(map(int, dic_element_value))
                db_results.append(sum(bkl_lst))

            else:
                db_results.append(';'.join(dic_element_value))

        db_results_ref = [x if x != '' else '0' for x in db_results]

        df_pd[db] = db_results_ref

    """
    Define columns to export
    """
    if args.verbose:
        print("Defining columns to export")

    columns_path = args.columns

    with open(columns_path, 'r') as colnames:
        cols_list = yaml.safe_load(colnames)
        fixed_cols = cols_list['fixed']
        samples_cols = cols_list['samples']

    cols = list(chain.from_iterable([fixed_cols, samples_cols]))

    """
    Melt data frame
    """
    if args.verbose:
        print("Melting data frame")

    outfile = df_pd[cols]

    outfile_melt = pd.melt(outfile, id_vars=fixed_cols,
                           value_vars=samples_cols,
                           var_name='SAMPLE', value_name='GT')

    """
    Filtering and reformatting
    """
    if args.verbose:
        print("Filtering and reformatting")

    # outfile_alt = outfile_melt.loc[-outfile_melt['GT'].isin(['0/0', './.', './0', '0/.', '0'])]
    outfile_alt = outfile_melt.loc[outfile_melt['GT'] != './.:NaN:0:0,0:--:NaN:NaN:NaN:NAN:NAN:NAN']

    ##add column for sample ID
    outfile_alt['ID'] = outfile_alt.GT.str.split(':').str[7]
    # outfile_alt['ID_sample'] = outfile_alt.GT.str.split(':').str[7]

    ##keep only GT from GT column
    outfile_alt['GT'] = outfile_alt.GT.str.split(':').str[0]
    # outfile_alt.replace(to_replace='GT', value=outfile_alt.GT.str.split(':').str[0]) #not working, to fix

    ##reorder by input columns file
    cols_ordered = list(chain.from_iterable([['SAMPLE', 'GT'], list(fixed_cols)]))
    # cols_ordered = list(chain.from_iterable([['SAMPLE', 'GT', 'ID_SAMPLE'], list(fixed_cols)]))
    outfile_out = outfile_alt[cols_ordered]

    """
    Save to file
    """
    if args.verbose:
        print("Saving to file")

    out_file_path = args.output
    out_file_melt_path = re.sub('tab', 'melt.tab', out_file_path)

    outfile.to_csv(out_file_path, sep='\t', index=False)
    outfile_out.to_csv(out_file_melt_path, sep='\t', index=False)



if __name__ == '__main__':
    main()