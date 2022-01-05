#! /usr/bin/env python3
import sys
import re

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
			r = re.search("AF=([0-9\.]+)",row[7])
			if r:
				alt = int(100*float(r.group(1)))
				ref = int(100-alt)
				row[9]+=f":{ref},{alt}"
			else:
				row[9]+=":0,100"
			sys.stdout.write("%s\n" % "\t".join(row))
