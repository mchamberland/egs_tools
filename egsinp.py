import sys
from os import listdir
from os.path import join, isfile, splitext, exists


def generate_transformations_from_dwells():
    all_files = [f for f in listdir('.') if isfile(join('.', f))]
    dwell_files = [f for f in all_files if splitext(f)[1] == '.dwell']
    start = ':start transformation:\n\ttranslation = '
    stop = '\n:stop transformation:\n\n'

    for f in dwell_files:
        the_transformations = ''
        if not exists(f):
            sys.exit('Could not open file.')

        with open(f, 'r') as file:
            all_values = [line.split() for line in file]
            transformation = [', '.join(values[0:3]) for values in all_values]
            the_weights = ', '.join(values[3] for values in all_values)

        for t in transformation:
            the_transformations += start + ''.join(t) + stop

        transf_filename = splitext(f)[0] + '.transf'
        with open(transf_filename, 'wt') as the_transf_file:
            the_transf_file.write(the_transformations)

        weight_filename = splitext(f)[0] + '.weights'
        with open(weight_filename, 'wt') as weight_file:
            weight_file.write(the_weights)
