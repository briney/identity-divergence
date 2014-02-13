#!/usr/bin/python
# filename: iden_div.py


###########################################################################
#
# Copyright (c) 2013 Bryan Briney.  All rights reserved.
#
# @version: 1.0.0
# @author: Bryan Briney
# @license: MIT (http://opensource.org/licenses/MIT) 
#
###########################################################################


import os
import argparse
from pymongo import MongoClient
from Bio import pairwise2, SeqIO
from multiprocessing import Pool, cpu_count
import matplotlib as mpl
mpl.use('pdf')
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser("For a MongoDB collection, plots the germline divergence against the sequence identity to a given 'subject' sequence.")
parser.add_argument('-d', '--database', dest='db', required=True, help="Name of the MongoDB database to query. Required")
parser.add_argument('-c', '--collection', dest='collection', default=None, help="Name of the MongoDB collection to query. If not provided, all collections in the given database will be processed iteratively.")
parser.add_argument('-o', '--output', dest='output', required=True, help="Output directory figure files. The figure file(s) will be 'output/<db>_<collection>_<standard>.pdf'. Required")
parser.add_argument('-i', '--ip', dest='ip', default='localhost', help="The IP address for the MongoDB server.  Defaults to 'localhost'.")
parser.add_argument('-p', '--port', dest='port', default=27017, help="The port for the MongoDB server.  Defaults to '27017'.")
parser.add_argument('-s', '--standard', dest='standard', required=True, help='Path to a file containing the standard sequence(s) for which identity/divergence will be calculated, in FASTA format. All sequences in the standard file will iteratively processed. Required')
parser.add_argument('-x', '--chain', dest='chain', default='heavy', choices=['heavy', 'kappa', 'lambda', 'light'], help="The chain type of the subject sequence.  Options are 'heavy', 'kappa', 'lambda' and 'light'.  Default is 'heavy'.")
parser.add_argument('-n', '--no_update', dest='no_update', action='store_true', default=False, help="Does not update the MongoDB with iden_div info. Can save some time if the idenentity calculations aren't needed again.")
args = parser.parse_args()


def get_standards():
	standards = []
	for s in SeqIO.parse(open(args.standard, 'r'), 'fasta'):
		standards.append([s.id, str(s.seq)])
	return standards

def get_collections():
	if args.collection:
		conn = MongoClient(args.ip, args.port)
		db   = conn[args.db]
		subjects = db.collection_names()
		subjects.remove('system.indexes')
		return sorted(subjects)
	return [args.collection,]

def get_chain():
	if args.chain == 'light': 
		return ['kappa', 'lambda']
	return [args.chain,]

def query(collection):
	conn = MongoClient(args.ip, args.port)
	db = conn[args.db]
	coll = db[collection]
	chain = get_chain()
	print_query_info()
	results = coll.find({'chain': {'$in': chain}},{'_id': 0, 'seq_id': 1, 'nt_identity.v': 1, 'vdj_aa': 1}).limit(25000)
	output = []
	for r in results:
		output.append([r['nt_identity']['v'], r['vdj_aa'], r['seq_id']])
	return output

def update_db(standard, scores, collection):
	conn = MongoClient(args.ip, args.port, max_pool_size=1000)
	db = conn[args.db]
	coll = db[collection]
	print_update_info()
	for score in scores:
		coll.find_and_modify(query={'seq_id': score[2]}, update={'$set': {'iden_div': {standard.lower(): float(score[1])}}})
	print_done()

def identity(standard, seqs):
	global scores
	scores = []
	print_single_standard(standard)
	pool = Pool(processes=cpu_count())
	for seq in seqs:
		pool.apply_async(do_alignment, args=(seq,standard[1]), callback=log_result)
	pool.close()
	pool.join()
	print_done()
	return standard[0]

def do_alignment(seq, standard):
	identity = seq[0]
	sequence = seq[1]
	seq_id = seq[2]
	score = pairwise2.align.globalxx(sequence, standard, one_alignment_only=1, score_only=1)
	norm_score = 100 * float(score) / max(len(sequence), len(standard))
	output = [identity, norm_score, seq_id]
	return output

def log_result(result):
	scores.append(result)

def make_figure(standard_id, scores, collection):
	print_fig_info()
	fig_file = os.path.join(args.output, '{0}_{1}_{2}.pdf'.format(args.db, collection, standard_id))
	x = [100.0 - s[0] for s in scores]
	y = [s[1] for s in scores]
	xmin = min(x)
	xmax = max(x)
	ymin = min(y)
	# ymax = max(y)
	# plot params
	plt.subplots_adjust(hspace=0.95)
	plt.subplot(111)
	plt.hexbin(x, y, bins='log', cmap=mpl.cm.jet, mincnt=2, gridsize=100)
	plt.title(standard_id, fontsize=18)
	# set and label axes
	plt.axis([xmin-2, xmax+2, ymin-2, 102])
	# plt.gca().invert_xaxis()
	plt.xlabel('Germline divergence')
	plt.ylabel('{0} identity'.format(standard_id))
	# make and label the colorbar
	cb = plt.colorbar()
	cb.set_label('Sequence count (log10)', labelpad=10)
	# save figure and close
	plt.savefig(fig_file)
	plt.close()
	print_done()

def print_standards_info(standards):
	print ''
	print ''
	print 'Found {} standard sequence(s):'.format(len(standards))
	print ', '.join([s[0] for s in standards])

def print_collections_info(collections):
	print ''
	print 'Found {} collection(s):'.format(len(collections))
	print ', '.join(collections)

def print_single_standard(standard):
	print ''
	print 'Standard ID: {}'.format(standard[0])
	print 'Calculating pairwise identities...'

def print_single_collection(collection):
	print ''
	print ''
	print '----------------------------------------'
	print 'Collection: {}'.format(collection)
	print '----------------------------------------'
	print ''

def print_query_info():
	print ''
	print 'Querying for comparison sequences...'

def print_fig_info():
	print ''
	print 'Making the identity/divergence figure...'

def print_update_info():
	print ''
	print 'Updating the MongoDB database with identity scores...'

def print_done():
	print 'Done.'

def main():
	standards = get_standards()
	print_standards_info(standards)
	collections = get_collections()
	print_collections_info(collections)
	for collection in collections:
		print_single_collection(collection)
		seqs = query(collection)
		for standard in standards:
			standard_id = identity(standard, seqs)
			make_figure(standard_id, scores, collection)
			if not args.no_update:
				update_db(standard_id, scores, collection)


if __name__ == '__main__':
	scores = []
	main()