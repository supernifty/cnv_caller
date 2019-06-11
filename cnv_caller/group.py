#!/usr/bin/env python
'''
  
'''

import argparse
import collections
import csv
import logging
import math
import sys

import version

def main(fh, minlen, minbases):
  logging.info('reading from stdin...')
  current = {'start': None, 'finish': None, 'chrom': None, 'status': None, 'bases': 0, 'count': 0}
  counts = collections.defaultdict(int)

  # chrom start finish tumour normal ratio outcome
  sys.stdout.write('chrom\tstart\tfinish\tsegments\tbases\tlength\tstatus\n')
  for row in csv.DictReader(fh, delimiter='\t'):
    row['start'] = int(row['start'])
    row['finish'] = int(row['finish'])
    counts['considered'] += row['finish'] - row['start']
    if row['outcome'] != current['status'] or row['chrom'] != current['chrom']: # new status
      if current['status'] != None and current['status'] != '0' and current['finish'] - current['start'] >= minlen and current['bases'] >= minbases:
        sys.stdout.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(current['chrom'], current['start'], current['finish'], current['count'], current['bases'], current['finish'] - current['start'], current['status']))
        counts[current['status']] += current['finish'] - current['start']
      current = {'chrom': row['chrom'], 'start': row['start'], 'finish': row['finish'], 'status': row['outcome'], 'count': 1, 'bases': row['finish'] - row['start']}
    else: # same status
      current['finish'] = row['finish'] 
      current['count'] += 1
      current['bases'] += row['finish'] - row['start']

  # check the last one
  if current['status'] != None and current['finish'] - current['start'] > minlen:
    sys.stdout.write('{}\t{}\t{}\t{}\t{}\n'.format(current['chrom'], current['start'], current['finish'], current['count'], current['bases'], current['finish'] - current['start'], current['status']))
      
  logging.info('done. %s', counts)

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Group and plot CNVs')
  parser.add_argument('--version', action='version', version=version.PROGRAM_VERSION)
  parser.add_argument('--minlen', required=False, default=10000, help='min length of abnormal CNV')
  parser.add_argument('--minbases', required=False, default=1000, help='min bases considered of abnormal CNV')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)

  main(sys.stdin, args.minlen, args.minbases)

