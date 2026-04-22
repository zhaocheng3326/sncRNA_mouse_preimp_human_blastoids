#!/usr/bin/python
import sys
import re
for line in sys.stdin:
    if line.startswith("@"):
        print(line,end="")
    else:
        sam=line.strip("\n").split("\t")
        rl=len(sam[9])
        if (sam[2] == "*"):
            continue
        for tag in sam[11:len(sam)]:
            md=sam[12].replace("MD:Z:","")
            if len(re.findall("NM:i:",tag)) > 0:
                tag_pos=tag.index("NM:i:")
                num_mismatch=int(tag[tag_pos+5:])
                if (md.startswith("0")):
                    num_mismatch=num_mismatch-1
                elif md.endswith("A0") | md.endswith("T0") | md.endswith("G0") | md.endswith("C0")  :
                    num_mismatch=num_mismatch-1
                if (rl <= 20 and num_mismatch <1):
                    print(line,end="")
                elif (rl <= 40 and rl >20 and num_mismatch <2):
                    print(line,end="")
                elif (rl > 40 and num_mismatch <3):
                    print(line,end="")
                continue


