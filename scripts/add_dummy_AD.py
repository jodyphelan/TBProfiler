#! /usr/bin/env python
import sys

for l in sys.stdin:
	if l[0]=="#":
		if "#CHROM" in l:
			sys.stdout.write('##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths (high-quality bases)">\n')
		sys.stdout.write(l)
	else:
		row = l.rstrip().split()
		if "AD" in row[8]:
			sys.stdout.write(l)
		else:
			row[8]+=":AD"
			row[9]+=":0,100"
			sys.stdout.write("%s\n" % "\t".join(row))
