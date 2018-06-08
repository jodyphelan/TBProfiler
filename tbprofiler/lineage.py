from __future__ import division
import sys
import math
import variant_calling
from collections import defaultdict

def stdev(arr):
	mean = sum(arr)/len(arr)
	return math.sqrt(sum([(x-mean)**2 for x in arr])/len(arr))

def lineage(self):
	calls = variant_calling.htsbox_calls(self,self.params["lineage_bed"])
	lin_meta_dict = {}
	lin_pos_dict = defaultdict(dict)
	for l in open(self.params["lineage_bed"]):
		arr = l.rstrip().split("\t")
		lin_meta_dict[arr[3]] = arr[6:]
		lin_pos_dict[arr[0]][arr[1]] = {"ref":arr[4],"alt":arr[5],"lin":arr[3]}

	lin_support = defaultdict(list)
	sample_lin = set()
	for chrom in lin_pos_dict:
		for pos in lin_pos_dict[chrom]:
			if pos not in calls[chrom]: continue
			if set([c[0] for c in calls[chrom][pos]])==set(["N"]): continue
			if lin_pos_dict[chrom][pos]["alt"] in [x[0] for x in calls[chrom][pos]]:
				lin_support[lin_pos_dict[chrom][pos]["lin"]].append((sum([x[2] for x in calls[chrom][pos] if x[0]==lin_pos_dict[chrom][pos]["alt"]]),sum([x[2] for x in calls[chrom][pos] if x[0]!=lin_pos_dict[chrom][pos]["alt"]])))
			else:
				lin_support[lin_pos_dict[chrom][pos]["lin"]].append((0,sum([x[2] for x in calls[chrom][pos] if x[0]!=lin_pos_dict[chrom][pos]["alt"]])))
	if self.params["verbose"]==2:
		sys.stderr.write("%s\n" % lin_support)
	lin_frac = defaultdict(float)
	for l in lin_support:
		print stdev([x[0]/(x[0]+x[1]) for x in lin_support[l]])
		if stdev([x[0]/(x[0]+x[1]) for x in lin_support[l]])>0.15: continue
		lin_pos_reads = sum([x[0] for x in lin_support[l]])
		lin_neg_reads = sum([x[1] for x in lin_support[l]])
		lf = lin_pos_reads/(lin_pos_reads+lin_neg_reads)
		if lf<0.05:continue
		lin_frac[l] = lf
		sample_lin.add(l)
	lin_str = ""
	lin_list = []
	for l in sorted(list(sample_lin)):
		lin_list.append({"lin":l,"frac":lin_frac[l],"family":lin_meta_dict[l][0],"spoligotype":lin_meta_dict[l][1],"rd":lin_meta_dict[l][2]})
		lin_str+= "%s\t%s\t%s\t%s\t%s\n" % (l,lin_frac[l],lin_meta_dict[l][0],lin_meta_dict[l][1],lin_meta_dict[l][2])
	return lin_list
