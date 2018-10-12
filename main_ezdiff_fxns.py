import os
import numpy as np
import csv
from csv import reader, writer
import pandas
import sys


data = sys.argv[1]
output_dir = sys.argv[2]


vars_of_interest = ['wmacck', 'wmrtk', 'wmsdk', 'g2acck', 'g2avrtk', 'g2sdrtk']


def ezdiff(MRT, VRT, pc, s=1.0):
    def logit(p): return numpy.log(p / (1 - p))
    if pc == 1:
        pc = .99

    r = (logit(pc) * (((pc**2) * logit(pc)) - pc * logit(pc) + pc - 0.5)) / VRT
    v = numpy.sign(pc - 0.5) * s * (r)**0.25
    a = (s**2 * logit(pc)) / v
    y = (-1 * v * a) / (s**2)
    MDT = (a / (2 * v)) * ((1 - numpy.exp(y)) / (1 + numpy.exp(y)))
    t = MRT - MDT
    v = v * .1
    a = a * .1
    return([v, a, t])


def ezdiff_parser(data, output_dir, gen_numbers):
    input_data_vals=extract_data_from_csv_into_dict(data)
    with open(data, 'U') as file:
        input_data_vals = reader(file)
        # zip iterates over the csv file, *operator
        input_data_vals = map(list, zip(*input_data_vals))
        new_dict = dict((rows[0], rows[1:]) for rows in (input_data_vals))

    key_val = []
    for key in new_dict.keys():
        for value in new_dict[key]:
            key_val.append([key, value])

    new_dict = dict()
    for variable in vars_of_interest:
        for val in new_dict[variable]:
            if ' ' in val:
                print 'empty'
            else:
                new_dict[variable] = []

    for key_val_idx, key_val_value in enumerate(key_val):
        if key_val_value[0] in new_dict.keys():
            new_dict[key_val_value[0]].append(key_val_value[1])

    df = pandas.DataFrame.from_dict(new_dict, orient='columns', dtype=int)

    if not os.path.exists(os.path.join(output_dir, 'parsed_output_dir.csv')):
        with open(os.path.join(output_dir, 'parsed_output_dir.csv'), 'wb') as csvfile:
            writer = csv.writer(csvfile, delimiter=',')
            df.to_csv(os.path.join(output_dir, 'parsed_output_dir.csv'),
                      header=True, index=True)

    else:
        print 'data has been parsed'

    if gen_numbers == 1:

        csv_for_new_vals = pandas.read_csv(
            os.path.join(output_dir, 'parsed_output_dir.csv'))

        wmrtk_conv = csv_for_new_vals['wmrtk_conv'].tolist()
        wmsdk_conv = csv_for_new_vals['wmsdk_conv'].tolist()
        wmacck_conv = csv_for_new_vals['wmacck_conv'].tolist()
        g2avrtk_conv = csv_for_new_vals['g2avrtk_conv'].tolist()
        g2sdrtk_conv = csv_for_new_vals['g2sdrtk_conv'].tolist()
        g2acck_conv = csv_for_new_vals['g2acck_conv'].tolist()

        with open(os.path.join(output_dir, 'ez_drift_diffusion_data.csv'), 'wb') as csvfile:
            writer = csv.writer(csvfile, delimiter=',')
            writer.writerow(
                ('CPT_drift_rate', 'CPT_boundary_separation', 'CPT_nondecision_time'))
        for idx, value in enumerate(wmrtk_conv):
            [v, a, t] = ezdiff(
                wmrtk_conv[idx], wmsdk_conv[idx], wmacck_conv[idx])
            writer.writerow((v, a, t))
        writer.writerow(
            ('GNG_drift_rate', 'GNG_boundary_separation', 'GNG_nondecision_time'))
        for idx, g2 in enumerate(g2avrtk_conv):
            [v, a, t] = ezdiff(g2avrtk_conv[idx],
                               g2sdrtk_conv[idx], g2acck_conv[idx])
            writer.writerow((v, a, t))
    elif gen_numbers == 0:
        print "data has been parsed"


ezdiff_parser(data, output_dir, sys.argv[3])
