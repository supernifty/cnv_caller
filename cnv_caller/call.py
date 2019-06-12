#!/usr/bin/env python
'''
  
'''

import argparse
import logging
import math
import sys
import version

import pysam

def main(tumour, normal, bed, minlen):
  logging.info('starting...')

  tumour_file = pysam.AlignmentFile(tumour, "rb")
  normal_file = pysam.AlignmentFile(normal, "rb")

  tumour_total_aligned = 0
  normal_total_aligned = 0
  regions = []

  sys.stdout.write('chrom\tstart\tfinish\ttumour\tnormal\tratio\toutcome\n')

  for region_idx, region in enumerate(open(bed, 'r')):
    fields = region.strip('\n').split('\t')
    if len(fields) < 3:
      continue
    chrom, start, finish = fields[0], int(fields[1]), int(fields[2])
    tumour_region_aligned = 0
    result = {'chrom': chrom, 'start': start, 'finish': finish, 'tumour': 0, 'normal': 0}
    for x in tumour_file.pileup(chrom, start, finish):
      aligned = x.get_num_aligned()
      tumour_total_aligned += aligned
      result['tumour'] += aligned
      
    for x in normal_file.pileup(chrom, start, finish):
      aligned = x.get_num_aligned()
      normal_total_aligned += aligned
      result['normal'] += aligned

    if region_idx < 10 or region_idx % 1000 == 0:
      logging.debug('%i regions added...', region_idx)

    ratio = (result['tumour'] + 2) / (result['normal'] + 1)
    regions.append(result)

  total_ratio = (tumour_total_aligned + 2) / (normal_total_aligned + 1)

  logging.debug('tumour reads: %i normal reads: %i ratio: %.2f', tumour_total_aligned, normal_total_aligned, tumour_total_aligned / normal_total_aligned)

  # assess each region
  start = {'type': None, 'chrom': None, 'pos': None}
  for result in regions:
    if result['tumour'] == 0:
      sys.stdout.write('{}\t{}\t{}\t{}\t{}\t{:.2f}\t{}\n'.format(result['chrom'], result['start'], result['finish'], result['tumour'], result['normal'], -999, '-2'))
    elif result['normal'] == 0:
      sys.stdout.write('{}\t{}\t{}\t{}\t{}\t{:.2f}\t{}\n'.format(result['chrom'], result['start'], result['finish'], result['tumour'], result['normal'], 999, '+2'))
    else:
      ratio = math.log(result['tumour'] / result['normal'] / total_ratio, 2)
      if ratio > 0.8:
        outcome = '+2'
      elif ratio > 0.3:
        outcome = '+1'
      elif ratio > -0.6:
        outcome = '0'
      elif ratio > -2:
        outcome = '-1'
      else:
        outcome = '-2'
      sys.stdout.write('{}\t{}\t{}\t{}\t{}\t{:.2f}\t{}\n'.format(result['chrom'], result['start'], result['finish'], result['tumour'], result['normal'], ratio, outcome))

  logging.info('done')

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Assess MSI')
  parser.add_argument('--version', action='version', version=version.PROGRAM_VERSION)
  parser.add_argument('--tumour', required=True, help='tumour bam')
  parser.add_argument('--normal', required=True, help='normal bam')
  parser.add_argument('--bed', required=True, help='bed file of regions')
  parser.add_argument('--minlen', required=False, default=10000, help='min length of LOH')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)

  main(args.tumour, args.normal, args.bed, args.minlen)

