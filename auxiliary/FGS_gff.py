#!/usr/bin/env python

import sys

if len(sys.argv) < 3:
	print sys.argv[0], "FGS-output new-GFF-file"
	sys.exit()

infile = open(sys.argv[1], "r")
outfile = open(sys.argv[2], "w")

print >>outfile, "##gff-version 3"
for aline in infile:
	aline = aline.strip()
	if aline[0] == '>':
		seqid = aline[1:]
	else:
		subs = aline.split()
		id = seqid + "_" + subs[0] + "_" + subs[1] + "_" + subs[2]
		des = [seqid, "FGS", "CDS", subs[0], subs[1], ".", subs[2], str(int(subs[3]) - 1), "ID=" + id + ";product=predicted protein"] 
		print >>outfile, "\t".join(des)
outfile.close()
infile.close()
		

