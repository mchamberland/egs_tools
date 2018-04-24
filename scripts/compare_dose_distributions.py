#!/home/mchamber/anaconda3/envs/py3/bin/python
import argparse
import dose_distribution.dose3d as dd
import dose_distribution.manip as dman
import matplotlib.pyplot as plt
from scipy.stats import norm

parser = argparse.ArgumentParser(description='Compare 2 3ddose distributions.')

parser.add_argument('dose_reference',
                    help='The reference dose.')
parser.add_argument('dose_compare',
                    help='The dose that is being compared to the reference.')

parser.add_argument('-t', '--threshold', dest='threshold', type=float, default=0,
                    help='The dose threshold for the comparison. The threshold can be a percentage of the maximum dose'
                         'in the reference distribution or can be an absolute number.')

parser.add_argument('-u', '--units', dest='units', default='%',
                    help='The units of the dose threshold. "%" (default) or "absolute".')

parser.add_argument('-n', '--nbins', dest='nbins', type=int, default=30,
                    help='The number of bins to use in the histogram. Defaults to 30.')


args = parser.parse_args()


dose_reference = dd.DoseDistribution(args.dose_reference)
dose_compare = dd.DoseDistribution(args.dose_compare)

t = dman.calculate_normalized_differences(dose_compare, dose_reference, args.threshold, args.units)
mean, stdev = norm.fit(t)
plt.hist(t, bins=30)
plt.title(r'$\mathrm{Histogram\ of\ t:}\ \mu=%.3f,\ \sigma=%.3f$' % (mean, stdev))
plt.show()
