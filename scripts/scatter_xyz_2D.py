import argparse
from os.path import splitext, basename
import csv
import numpy as np
import matplotlib.pyplot as plt

AXES_DICT = {'x': 0, 'y': 1, 'z': 2}

parser = argparse.ArgumentParser(description='Produce a 2D scatter plot from a list of 3D (x,y,z) coordinates.')

parser.add_argument('the_file',
                    help='The text file containing the coordinates of the points to plot.')

parser.add_argument('-a', '--axes', dest='axes', type=str, nargs=2, default='x y',
                    help='The 2 coordinates to plot (e.g., "x y", "y z", or "x z".')

parser.add_argument('-c', '--compare', dest='compare', type=str,
                    help='A text file containing a second dataset to compare to.')


args = parser.parse_args()

ax1 = args.axes[0]
ax2 = args.axes[1]

with open(args.the_file) as f:
    dialect = csv.Sniffer().sniff(f.read(1024))
    f.seek(0)
    tmp = list(csv.reader(f, dialect, quoting=csv.QUOTE_NONNUMERIC))
    the_coordinates = np.array(list(map(list, zip(*tmp))))

plt.scatter(the_coordinates[AXES_DICT[ax1]], the_coordinates[AXES_DICT[ax2]], label=args.the_file)
plt.xlabel(ax1)
plt.ylabel(ax2)
plt.title(args.the_file)

if args.compare:
    with open(args.compare) as f:
        dialect = csv.Sniffer().sniff(f.read(1024))
        f.seek(0)
        tmp = list(csv.reader(f, dialect, quoting=csv.QUOTE_NONNUMERIC))
        the_comparison = np.array(list(map(list, zip(*tmp))))

    plt.scatter(the_comparison[AXES_DICT[ax1]], the_comparison[AXES_DICT[ax2]], c='r', label=args.compare)
    plt.legend(loc='upper left')

plt.savefig('{}.png'.format(splitext(basename(args.the_file))[0]))
plt.show()
