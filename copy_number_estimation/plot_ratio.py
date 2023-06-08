'''
Plot the SUN ratio of each individual
'''

import matplotlib.pyplot as plt
import numpy as np
import argparse

# Read SUN fractions of each SUN set in th order: NOTCH2NLA/B, NOTCH2NLC, NOTCH2NLR, unique NOTCH2NLR, NOTCH2
parser = argparse.ArgumentParser()
parser.add_argument('AB')
parser.add_argument('C')
parser.add_argument('R')
parser.add_argument('uniR')
parser.add_argument('N2')
args = parser.parse_args()

# Get sample number
sample = args.AB.lstrip('all_A_B_R48_only_').rstrip('_pos_all.data')

plt.style.use("seaborn")

# Load data into variables
x1, y1 = np.loadtxt(args.AB, delimiter=' ', unpack=True)
x2, y2 = np.loadtxt(args.C, delimiter=' ', unpack=True)
x3, y3 = np.loadtxt(args.R, delimiter=' ', unpack=True)
x4, y4 = np.loadtxt(args.uniR, delimiter=' ', unpack=True)
x5, y5 = np.loadtxt(args.N2, delimiter=' ', unpack=True)

# Plot
fig, (ax1, ax2, ax3, ax4, ax5) = plt.subplots(5, sharex=True, sharey=True)  # Five subplots
fig.suptitle(sample + ' SUN Depth')  # Figure title
ax1.bar(x1, y1, width=20, ec='#2C80C0')  # Bar plot
ax2.bar(x2, y2, width=20, ec='#2C80C0')
ax3.bar(x3, y3, width=20, ec='#2C80C0')
ax4.bar(x4, y4, width=20, ec='#2C80C0')
ax5.bar(x5, y5, width=20, ec='#2C80C0')

ax1.set_title('NOTCH2NLA/B', fontsize=7.5, x=0.5, y=0.95)  # Titles of each bar plot
ax2.set_title("NOTCH2NLC", fontsize=7.5, x=0.5, y=0.95)
ax3.set_title("NOTCH2NLR", fontsize=7.5, x=0.5, y=0.95)
ax4.set_title("Unique 6 NOTCH2NLR", fontsize=7.5, x=0.5, y=0.95)
ax5.set_title("NOTCH2", fontsize=7.5, x=0.5, y=0.95)
plt.yticks(np.arange(0, 1.0, 0.2))
plt.xticks(np.arange(0, 81000, 10000))
fig.text(0.5, 0.04, 'Alignment Position to GRCh38 NOTCH2', ha='center', va='center')
fig.text(0.06, 0.5, 'SUN Fraction of Reads', ha='center', va='center', rotation='vertical')

# plt.show()
plt.savefig(sample + '_SUN_ratio', dpi = 600)