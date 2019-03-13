#! /usr/bin/env python
import sys
import pathogenprofiler as pp
import argparse
import json
import tbprofiler as tbp
import os


def reprofile(args):
	old_results = json.load(open(args.json))
	conf_file = pp.filecheck(args.db)
	conf = json.load(open(conf_file))
	new_results = old_results.copy()
	variant_dump = []
	for var in old_results["dr_variants"]:
		del var["drug"]
		var["gene_id"] = var["locus_tag"]
		variant_dump.append(var)
	for var in old_results["other_variants"]:
		variant_dump.append(var)
		var["gene_id"] = var["locus_tag"]
	new_results["variants"] = variant_dump
	del new_results["other_variants"]
	del new_results["dr_variants"]
	new_results = pp.db_compare(db_file=conf["json_db"],mutations=new_results)
	tbp.reformat_annotations(new_results)
	json.dump(new_results,open("%s.results.json"%args.prefix,"w"))


def main_lineage(args):
	conf_file = pp.filecheck(tbp._ROOT+"/../"+args.db+".config.json")
	conf = json.load(open(conf_file))
	pp.filecheck(args.bcf)
	bcf = pp.bcf(args.bcf)
	mutations = bcf.get_bed_gt(conf["barcode"],conf["ref"])
	results = {}
	results["barcode"] = pp.barcode(mutations,conf["barcode"])
	print(results["barcode"])
	tbp.barcode2lineage(results)
	if args.prefix:
		outfile = "%s.lineage.%s" % (args.prefix,args.outfmt)
		O = open(outfile,"w")
		if args.outfmt=="json":
			json.dump(results["lineage"],O)
		elif args.outfmt=="txt":
			O.write(tbp.text.lineagejson2text(results["lineage"]))
		O.close()

if __name__=="__main__":
	parser = argparse.ArgumentParser(description='TBProfiler Utils',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	subparsers = parser.add_subparsers(help="Task to perform")

	parser_sub = subparsers.add_parser('reprofile', help='Create a phylogeny based on SNPs', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser_sub.add_argument('json',default=None,type=str,help='Sample prefix')
	parser_sub.add_argument('--prefix','-p',default="tbprofiler",help='Sample prefix')
	parser_sub.add_argument('--db',default='drdb_v2',help='Full path to mutation database json file to use (default: TBProfiler panel)')
	parser_sub.set_defaults(func=reprofile)


	parser_sub = subparsers.add_parser('gbcf_lineage', help='Create a phylogeny based on SNPs', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser_sub.add_argument('bcf',default=None,type=str,help='Sample prefix')
	parser_sub.add_argument('--prefix','-p',default="tbprofiler",help='Sample prefix')
	parser_sub.add_argument('--outfmt',default='json',choices=["json","txt"],type=str,help="Output format")
	parser_sub.add_argument('--db',default='drdb_v2',help='Full path to mutation database json file to use (default: TBProfiler panel)')
	parser_sub.set_defaults(func=main_lineage)

	args = parser.parse_args()
	args.func(args)
