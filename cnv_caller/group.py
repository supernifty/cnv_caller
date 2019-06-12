#!/usr/bin/env python
'''
  
'''

import argparse
import collections
import csv
import logging
import math
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
matplotlib.rcParams.update({'font.size': 18})
import matplotlib.patches as mpatches

WIDTH, HEIGHT = 18, 12

COLOURS = {
  '-2': '#ff0000', 
  '-1': '#ffa500',
  '0':  '#00ff00',
  '+1': '#ee82ee',
  '+2': '#0000ff'
}

LABELS = [
  ('Hom Deletion', '-2'),
  ('Het Deletion', '-1'),
  ('No change', '0'),
  ('Het copy', '+1'),
  ('2+ copies', '+2')
]

import version

def add_cnv(out, current, calls, counts, minlen, minbases, mincount):
  if current['status'] != None and current['status'] != '0' and current['finish'] - current['start'] >= minlen and current['bases'] >= minbases and current['count'] >= mincount:
    sys.stdout.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(current['chrom'], current['start'], current['finish'], current['count'], current['bases'], current['finish'] - current['start'], current['status']))
    calls.append({'start': current['start'], 'finish': current['finish'], 'status': current['status']})
    counts[current['status']] += current['finish'] - current['start']

def plot_chrom(plot, chrom, regions, xs, ys):
  filename = '{}.{}.png'.format(plot, chrom)
  plt.figure(figsize=(WIDTH, HEIGHT))
  for key in xs.keys():
    plt.scatter(xs[key], ys[key], color=COLOURS[key], alpha=0.4, s=4)
  plt.title('CNVs for chromosome {}'.format(chrom))
  for call in regions:
    plt.axvspan(call['start'], call['finish'], facecolor=COLOURS[call['status']], alpha=0.5)

  patches = []
  for label in LABELS:
    patches.append(mpatches.Patch(color=COLOURS[label[1]], label=label[0]))
    
  plt.legend(patches, [label[0] for label in LABELS], loc='lower left', ncol=5, labelspacing=0., framealpha=0.7 )
  plt.tight_layout()
  plt.savefig(filename)
  plt.close()
  logging.info('plotted %s', filename)

def main(fh, minlen, minbases, mincount, plot):
  logging.info('reading from stdin...')
  current = {'start': None, 'finish': None, 'chrom': None, 'status': None, 'bases': 0, 'count': 0}
  counts = collections.defaultdict(int)

  # chrom start finish tumour normal ratio outcome
  xs = collections.defaultdict(list)
  ys = collections.defaultdict(list)
  regions = []
  last_chrom = None

  sys.stdout.write('chrom\tstart\tfinish\tsegments\tbases\tlength\tstatus\n')
  for row in csv.DictReader(fh, delimiter='\t'):
    row['start'] = int(row['start'])
    row['finish'] = int(row['finish'])
    counts['considered'] += row['finish'] - row['start']
    if row['outcome'] != current['status'] or row['chrom'] != current['chrom']: # new status
      add_cnv(sys.stdout, current, regions, counts, minlen, minbases, mincount)
      current = {'chrom': row['chrom'], 'start': row['start'], 'finish': row['finish'], 'status': row['outcome'], 'count': 1, 'bases': row['finish'] - row['start']}
    else: # same status
      current['finish'] = row['finish'] 
      current['count'] += 1
      current['bases'] += row['finish'] - row['start']

    if plot is not None and last_chrom is not None and row['chrom'] != last_chrom:
      plot_chrom(plot, last_chrom, regions, xs, ys)
      xs = collections.defaultdict(list)
      ys = collections.defaultdict(list)
      regions = []

    xs[row['outcome']].append((row['start'] + row['finish']) / 2)
    ys[row['outcome']].append(max(-4, min(4, float(row['ratio'])))) # limit to -4 to +4

    last_chrom = row['chrom']

  # check the last one
  add_cnv(sys.stdout, current, regions, counts, minlen, minbases, mincount)
  if plot is not None:
    plot_chrom(plot, last_chrom, regions, xs, ys)
      
  logging.info('done. %s', counts)

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Group and plot CNVs')
  parser.add_argument('--version', action='version', version=version.PROGRAM_VERSION)
  parser.add_argument('--minlen', required=False, default=1000, help='min length of abnormal CNV')
  parser.add_argument('--minbases', required=False, default=1000, help='min bases considered of abnormal CNV')
  parser.add_argument('--mincount', required=False, default=3, help='min bases considered of abnormal CNV')
  parser.add_argument('--plot', required=False, help='prefix for plots')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)

  main(sys.stdin, args.minlen, args.minbases, args.mincount, args.plot)

