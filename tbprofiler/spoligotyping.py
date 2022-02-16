import os 


def counts2spoligotype(counts,cutoff=None):
    spacers = []
    if cutoff==None:
        cutoff = min([10,max([x["count"] for x in counts])*0.2])
    for k in counts:
        spacers.append("1" if k['count']>cutoff else "0")
    
    octal = []
    for i in range(0,40,3):
        tmp = "".join([spacers[i],spacers[i+1],spacers[i+2]])
        if tmp=="000":octal.append("0")
        elif tmp=="001":octal.append("1")
        elif tmp=="010":octal.append("2")
        elif tmp=="011":octal.append("3")
        elif tmp=="100":octal.append("4")
        elif tmp=="101":octal.append("5")
        elif tmp=="110":octal.append("6")
        elif tmp=="111":octal.append("7")
        # else:ps.log("Don't know what to do with %s" % tmp,ext=T)

    octal.append("0" if spacers[42]=="0" else "1")
    sitvit_str = "".join(["n" if x=="1" else "o" for x in spacers])
    binary_str = "".join(spacers)
    octal_str = "".join(octal)
    return binary_str,octal_str
