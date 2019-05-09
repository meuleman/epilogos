#!/usr/bin/env python
import sys
import csv
import errno

try: # necessary because the input will be the output from bgzip
    parser = csv.reader(sys.stdin, delimiter='\t')
    linenum=0;
    for line in parser:
        linenum = linenum + 1
        site = line[0:3]
        scores = line[3:]
        data = sorted(zip(scores,list(range(1,len(scores)+1))), key = lambda tup: float(tup[0]))
        sys.stdout.write(site[0]+'\t'+site[1]+'\t'+site[2]+'\tid:'+str(linenum)+',qcat:[ ')
        sys.stdout.write('['+data[0][0]+','+str(data[0][1])+']')
        for i in range(1,len(scores)):
            sys.stdout.write(', ['+data[i][0]+','+str(data[i][1])+']')
        sys.stdout.write(' ]\n')
except IOError as e: # python 2.X
    if e.errno != errno.EPIPE:
        raise
except BrokenPipeError as e: # python 3.X
    if e.errno != errno.EPIPE:
        raise
exit()
