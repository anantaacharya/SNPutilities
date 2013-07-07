#ananta Acharya
#ananta@uga.edu

#FormatChange
#compatible formats
#AA,AT,TT---Diploid
#A/A A/T T/T--DiploidSep
#A, A/T, T--hetero
#A,W,T--IUPAC
#0,1,2--012
#0,1 for each SNP---binary or 01

#The following are supoorted for output only
#Beagle
#binary by gen
#1,2,3,4 for ACGT


def main():
    print "Format change"
    print """#compatible formats
        1. AA,AT,TT---Diploid
        2. A/A A/T T/T--DiploidSep
        3. A, A/T, T--hetero
        4. A,W,T--IUPAC
        5. 0,1,2--012
        6. 0,1 for each SNP---binary or 01
        7. Beagle (A/T in subsequent lines)
        8. binary by gen (binary in beagle)
        9. 1,2,3,4 for ACGT
        """
    fromFormat=raw_input("Change format from: \n1st col: Marker Names, 1st row: Genotype Names\n 2nd col single alleles for binary, A/T for 012")
    toFormat=raw_input("Change format to: ")
    inf=raw_input("input file location")
    sep=raw_input("column separator")
    if sep=="tab":
        sep="\t"
    elif sep=="space":
        sep=" "
    if toFormat in ["2","3"]:
        sep2=raw_input("Allele separator within genotype ie / for A/T, : for A:T ")
    missing=raw_input("missing Data")
    mono=raw_input("remove monomorphic? y/n")
    outf=raw_input("out file location")
    
    if fromFormat=="1":
        
        header,SNPlist,IUPACMatrix=Dip2IUPAC(inf,sep,missing)
    elif fromFormat=="2":
        sep2=raw_input("Allele separator within genotype ie / for A/T, : for A:T: ")
        header,SNPlist,IUPACMatrix=DipSep2IUPAC(inf,sep,sep2,missing)
    elif fromFormat=="3":
        sep2=raw_input("Allele separator within genotype ie / for A/T, : for A:T: ")
        header,SNPlist,IUPACMatrix=Hetero2IUPAC(inf,sep,sep2, missing)
    elif fromFormat=="4":
        
        header,SNPlist,IUPACMatrix=IUPAC2IUPACread(inf,sep,missing)
    elif fromFormat=="5":
        codes=raw_input("codes: (0,1,2), oe 1,2,3, or -1, 0, 1")
        flag=raw_input("Are alleles in 2nd column?: ")
        
        if flag in ["y","yes", "Y","YES", "Yes"]:
            flag=True
        else:
            flag=False
        codes=codes.strip().split(",")
        header,SNPlist,IUPACMatrix=f0122IUPAC(inf,sep,flag,missing, codes)
    elif fromFormat=="6":
        print "Expects single SNP in second column and SNP name on first column"
        header,SNPlist,IUPACMatrix = binary2IUPAC(inf,sep,missing)
    elif fromFormat=="7":
        
        header,SNPlist,IUPACMatrix = Beagle2IUPAC(inf,sep,missing)
    print "Reading finished\n removing missing and/or monomorphic"    
    monof=False
    if "y" in mono:
        monof=True
    header,SNPlist,IUPACMatrix=removeMissing(header,SNPlist,IUPACMatrix, monof)
    print "now writing" +str(len(SNPlist))+" SNPs"
    if toFormat=="1":
        IUPAC2Dip(header,SNPlist,IUPACMatrix,sep,outf)
        
    elif toFormat=="2":
        
        IUPAC2DipSep(header,SNPlist,IUPACMatrix,sep,sep2,outf)

    elif toFormat=="3":
        
        IUPAC2Hetero(header,SNPlist,IUPACMatrix,sep,sep2,outf)

    elif toFormat=="4":
        IUPAC2IUPACwrite(header,SNPlist,IUPACMatrix,sep,outf)

    elif toFormat=="5":
        IUPAC2012(header,SNPlist,IUPACMatrix,sep,outf)

    elif toFormat=="6":
        IUPAC2binary(header,SNPlist,IUPACMatrix,sep,outf)
    elif toFormat=="7":
        IUPAC2Beagle(header,SNPlist,IUPACMatrix,sep,outf)
    elif toFormat=="8":
        IUPAC2binarybygen(header,SNPlist,IUPACMatrix,sep,outf)
    elif toFormat=="9":
        IUPAC21234(header,SNPlist,IUPACMatrix,sep,outf)

def Dip2IUPAC(inf,sep,missing):
    with open(inf,"r") as infile:
        IUPACMatrix=dict()
        header=infile.next().strip().split(sep)[1::]
        header=[geno.replace('"','') for geno in header]
        
        SNPlist=list()
        for line in infile:
           linelist=line.strip().split(sep)
           SNP=linelist[0].replace('"','')
           SNPlist.append(SNP)
           eachline=linelist[1::]
           for index, eachgenotype in enumerate(eachline):
               eachgenotype=eachgenotype.replace('"','')
               if missing in eachgenotype:
                   IUPACMatrix[SNP,header[index]]="NA"
                   
               else:
                   IUPACMatrix[SNP,header[index]]=iupac([eachgenotype[0],eachgenotype[1]])
                
            

    return header,SNPlist,IUPACMatrix


def DipSep2IUPAC(inf,sep,sep2, missing):
    with open(inf,"r") as infile:
        IUPACMatrix=dict()
        header=infile.next().strip().split(sep)[1::]
        SNPlist=list()
        
        for line in infile:
           linelist=line.strip().split(sep)
           SNP=linelist[0]
           SNPlist.append(SNP)
           eachline=linelist[1::]
           for index, eachgenotype in enumerate(eachline):
               if missing in eachgenotype:
                   IUPACMatrix[SNP,header[index]]="NA"
                   
               else:
               
                    IUPACMatrix[SNP,header[index]]=iupac(eachgenotype.strip().split(sep2))
                
            

    return header,SNPlist,IUPACMatrix

    

def Hetero2IUPAC(inf,sep, sep2, missing):
    with open(inf,"r") as infile:
        IUPACMatrix=dict()
        header=infile.next().strip().split(sep)[1::]
        SNPlist=list()
        for i,each in enumerate(header):
            print i,each
        for line in infile:
           linelist=line.strip().split(sep)
           SNP=linelist[0]
           SNPlist.append(SNP)
           eachline=linelist[1:len(header)+1]
           for index, eachgenotype in enumerate(eachline):
               
               if missing in eachgenotype:
                   IUPACMatrix[SNP,header[index]]="NA"
                   
               else:
               
                    IUPACMatrix[SNP,header[index]]=iupac(eachgenotype.strip().split(sep2))
    return header,SNPlist,IUPACMatrix


def Beagle2IUPAC(inf,sep, missing):
    with open(inf,"r") as infile:
        IUPACMatrix=dict()
        header=infile.next().strip().split(sep)[2::2]
        SNPlist=list()
        for i,each in enumerate(header):
            print i,each
        for line in infile:
           linelist=line.strip().split(sep)
           SNP=linelist[1]
           SNPlist.append(SNP)
           eachline=linelist[2:2*len(header)+2]
           pairlist=[[eachline[i],eachline[i+1]] for i in range(0,len(eachline)-1,2)]
           for index, eachgenotype in enumerate(pairlist):
               
               if missing in eachgenotype:
                   IUPACMatrix[SNP,header[index]]="NA"
                   
               else:
               
                    IUPACMatrix[SNP,header[index]]=iupac(eachgenotype)
    return header,SNPlist,IUPACMatrix

def IUPAC2IUPACread(inf,sep, missing):
    with open(inf,"r") as infile:
        IUPACMatrix=dict()
        header=infile.next().strip().split(sep)[1::]
        SNPlist=list()
        
        for line in infile:
           linelist=line.strip().split(sep)
           SNP=linelist[0]
           SNPlist.append(SNP)
           eachline=linelist[1::]
           for index, eachgenotype in enumerate(eachline):
               if missing in eachgenotype:
                   IUPACMatrix[SNP,header[index]]="NA"
                   
               else:
                   IUPACMatrix[SNP,header[index]]=eachgenotype
    return header,SNPlist,IUPACMatrix

def f0122IUPAC(inf,sep,flag, missing, codes):
    with open(inf,"r") as infile:
        IUPACMatrix=dict()
        
        SNPlist=list()
        if flag:
            header=infile.next().strip().split(sep)[2::]
            for line in infile:
               linelist=line.strip().split(sep)
               SNP=linelist[0]
               SNPlist.append(SNP)
               alleles=linelist[1].strip().split("/")
               snp1, snp2=alleles[0].upper(), alleles[1].upper()
               
               eachline=linelist[2::]
               for index, eachgenotype in enumerate(eachline):
                   if missing in eachgenotype:
                       IUPACMatrix[SNP,header[index]]="NA"
                   
                   
                   elif eachgenotype==codes[0]:
                       
                       IUPACMatrix[SNP,header[index]]=snp1
                   elif  eachgenotype==codes[2]:
                       
                       IUPACMatrix[SNP,header[index]]=snp2
                   elif eachgenotype==codes[1]:
                       IUPACMatrix[SNP,header[index]]=iupac([snp1,snp2])
                   else:
                       IUPACMatrix[SNP,header[index]]="NA"
                   
        else:
            header=infile.next().strip().split(sep)[1::]
            snp1, snp2="A","T"
            for line in infile:
               linelist=line.strip().split(sep)
               SNP=linelist[0]
           
               SNPlist.append(SNP)
               
               eachline=linelist[1::]
               for index, eachgenotype in enumerate(eachline):
                   if missing in eachgenotype:
                       IUPACMatrix[SNP,header[index]]="NA"
                   
                   
                   elif eachgenotype==codes[0]:
                       
                       IUPACMatrix[SNP,header[index]]=snp1
                   elif  eachgenotype==codes[2]:
                       
                       IUPACMatrix[SNP,header[index]]=snp2
                   elif eachgenotype==codes[1]:
                        IUPACMatrix[SNP,header[index]]=iupac([snp1,snp2])
                   else:
                       IUPACMatrix[SNP,header[index]]="NA"
            
    return header,SNPlist,IUPACMatrix

def binary2IUPAC(inf,sep, missing):
    with open(inf,"r") as infile:
        IUPACMatrix=dict()
        header=infile.next().strip().split(sep)[2::]
        SNPlist=list()
        SNPalleles=dict()
        binaryMatrix=dict()
        for line in infile:
           linelist=line.strip().split(sep)
           SNP=linelist[0]
           if SNP not in SNPlist:
               SNPlist.append(SNP)
           SNPalleles.setdefault(SNP,[]).append(linelist[1].upper())
           
           
           eachline=linelist[2::]
           for index, eachgenotype in enumerate(eachline):
               if missing in eachgenotype:
                   binaryMatrix[SNP,linelist[1],header[index]]="NA"
                   
               else:
                   binaryMatrix[SNP,linelist[1],header[index]]=eachgenotype
        for SNP in SNPlist:
            alleles= SNPalleles[SNP]
            for index, geno in enumerate(header):
                if (binaryMatrix[SNP,alleles[0],header[index]]=="NA") and (binaryMatrix[SNP,alleles[1],header[index]]=="NA"):
                    IUPACMatrix[SNP,header[index]]="NA"
                elif binaryMatrix[SNP,alleles[0],header[index]]=="1":
                    if binaryMatrix[SNP,alleles[1],header[index]]=="1":
                
                        IUPACMatrix[SNP,header[index]]=iupac([alleles[0],alleles[1]])
                    else:
                        IUPACMatrix[SNP,header[index]]=alleles[0]
                elif binaryMatrix[SNP,alleles[1],header[index]]=="1":
                    IUPACMatrix[SNP,header[index]]=alleles[1]
                

                
    return header,SNPlist,IUPACMatrix

def IUPAC2Dip(header,SNPlist,IUPACMatrix,sep,outf):
    towrite="SNP"+sep
    towrite+=sep.join(header)
    towrite+="\n"
    for SNP in SNPlist:
        eachline=SNP
        for geno in header:
            
            gen=IUPACMatrix[SNP,geno]
            if gen=="NA":
                eachline+=sep+"NA"
            else:
                gen=reverse_iupac(gen)
                if len(gen)==1:
                    
                    eachline+=sep+gen+gen
                else:
                    eachline+=sep+"".join(gen)
        eachline+="\n"
        towrite+=eachline
    with open(outf,"a") as outfile:
        outfile.write(towrite)
    print "file written in %s" %outf        

    
def IUPAC2DipSep(header,SNPlist,IUPACMatrix,sep,sep2,outf):
    
    towrite="SNP"+sep
    towrite+=sep.join(header)
    towrite+="\n"
    for SNP in SNPlist:
        eachline=SNP
        for geno in header:
            gen=IUPACMatrix[SNP,geno]
            if gen=="NA":
                eachline+=sep+"NA"
            else:
                gen=reverse_iupac(gen)
                if len(gen)==1:
                    
                    eachline+=sep+gen+sep2+gen
                else:
                    eachline+=sep+sep2.join(gen)
        eachline+="\n"
        towrite+=eachline
    with open(outf,"a") as outfile:
        outfile.write(towrite)
    print "file written in %s" %outf


def IUPAC2Beagle(header,SNPlist,IUPACMatrix,sep,outf):
    sep="\t"
    towrite="I"+sep+"SNP"+sep
    towrite+=sep.join([each+sep+each for each in header])
    towrite+="\n"
    for SNP in SNPlist:
        eachline="M"+sep+SNP
        for geno in header:
            gen=IUPACMatrix[SNP,geno]
            if gen=="NA":
                eachline+=sep+"NA"+sep+"NA"
            else:
                gen=reverse_iupac(gen)
                if len(gen)==0:
                    eachline+=sep+"NA"+sep+"NA"
                elif len(gen)==1:
                    
                    eachline+=sep+gen+sep+gen
                else:
                    eachline+=sep+gen[0]+sep+gen[1]
        eachline+="\n"
        towrite+=eachline
    with open(outf,"a") as outfile:
        outfile.write(towrite)
    print "file written in %s" %outf   
def IUPAC2Hetero(header,SNPlist,IUPACMatrix,sep,sep2,outf):
    towrite="SNP"+sep
    towrite+=sep.join(header)
    towrite+="\n"
    for SNP in SNPlist:
        eachline=SNP
        for geno in header:
            gen=IUPACMatrix[SNP,geno]
            if gen=="NA":
                eachline+=sep+"NA"
            else:
                gen=reverse_iupac(gen)
                    
                if len(gen)==1:
                    
                    eachline+=sep+gen
                else:
                    eachline+=sep+sep2.join(gen)
        eachline+="\n"
        towrite+=eachline
    with open(outf,"a") as outfile:
        outfile.write(towrite)
    print "file written in %s" %outf   

def IUPAC2IUPACwrite(header,SNPlist,IUPACMatrix,sep,outf):
    towrite="SNP"+sep
    towrite+=sep.join(header)
    towrite+="\n"
    for SNP in SNPlist:
        eachline=SNP
        for geno in header:
            eachline+=sep+IUPACMatrix[SNP,geno]
        eachline+="\n"
        towrite+=eachline
    with open(outf,"a") as outfile:
        outfile.write(towrite)
    print "file written in %s" %outf

def IUPAC2012(header,SNPlist,IUPACMatrix,sep,outf):
    towrite="SNP"+sep+"alleles"+sep
    towrite+=sep.join(header)
    towrite+="\n"
    for SNP in SNPlist:
        
        SNPs=set()
        for geno in header:
            
            gen=reverse_iupac(IUPACMatrix[SNP,geno])
            SNPs|=set(gen)
        SNPslist=list(SNPs)
        
        if len(SNPslist)>1:
            eachline=SNP+sep+SNPslist[0].upper()+"/"+SNPslist[1].lower()
        
            for geno in header:
                gen=IUPACMatrix[SNP,geno]
                if gen=="NA":
                    eachline+=sep+"NA"
                else:
                    
                    if gen in ["W","R","M","K","S","Y"]:
                        
                        eachline+=sep+"1"
                    elif gen==SNPslist[0]:
                        eachline+=sep+"0"
                    elif gen==SNPslist[1]:
                        eachline+=sep+"2"
                    else:
                        eachline+=sep+"NA"
            eachline+="\n"
            towrite+=eachline
        elif len(SNPslist)>0:
            eachline=SNP+sep+SNPslist[0].upper()+"/"+SNPslist[0].lower()
        
            for geno in header:
                gen=IUPACMatrix[SNP,geno]
                if gen==SNPslist[0]:
                    eachline+=sep+"0"
                else:

                    eachline+=sep+"NA"

            eachline+="\n"
            towrite+=eachline
    with open(outf,"a") as outfile:
        outfile.write(towrite)
    print "file written in %s" %outf       

def IUPAC2binary(header,SNPlist,IUPACMatrix,sep, outf):


    towrite="SNP"+sep+"allele"+sep
    towrite+=sep.join(header)
    towrite+="\n"
    for SNP in SNPlist:
        
        SNPs=set()
        for geno in header:
            
            gen=reverse_iupac(IUPACMatrix[SNP,geno])
            SNPs|=set(gen)
        SNPs=list(SNPs)
        
        eachline=SNP+sep+SNPs[0]
        for geno in header:
            gen=IUPACMatrix[SNP,geno]
            if gen=="NA":
                eachline+=sep+"NA"
            else:
                
                if gen in ["W","R","M","K","S","Y"]:
                    
                    eachline+=sep+"1"
                elif gen==SNPs[0]:
                    eachline+=sep+"1"
                elif gen==SNPs[1]:
                    eachline+=sep+"0"
        eachline+="\n"
        if len(SNPs)>1:
            eachline+=SNP+sep+SNPs[1]
            for geno in header:
                gen=IUPACMatrix[SNP,geno]
                if gen=="NA":
                    eachline+=sep+"NA"
                else:
                    
                    if gen in["W","R","M","K","S","Y"]:
                        
                        eachline+=sep+"1"
                    elif gen==SNPs[0]:
                        eachline+=sep+"0"
                    elif gen==SNPs[1]:
                        eachline+=sep+"1"
            eachline+="\n"
        towrite+=eachline
    with open(outf,"a") as outfile:
        outfile.write(towrite)
    print "file written in %s" %outf   

def IUPAC21234(header,SNPlist,IUPACMatrix,sep, outf):
    sep="\t"

    towrite=sep*2
    towrite+=(sep*2).join(SNPlist)
    towrite+="\n"
    snpcodedict={"A":"1","C":"2","G":"3","T":"4"}
    pops=list()
    popdict=dict()
    for geno in header:
        geno=geno.replace('"','')
        pops.append(geno.split("_")[0])
    pops=list(set(pops))
    x=100
    for i,pop in enumerate(pops):
        if "AMP" in pop:
            idd=pop[3::]
        else:
            idd=x
            x+=1
        popdict[pop]=str(idd)
    
    for geno in header:
        
        
        eachline=geno+sep+popdict[geno.replace('"','').split("_")[0]]
       # print geno
        for SNP in SNPlist:
          #  print SNP
            gen=reverse_iupacNA(IUPACMatrix[SNP,geno])
           # print IUPACMatrix[SNP,geno],gen
          #  raw_input()
            if "NA" in gen:
                eachline+=sep+"-9"+sep+"-9"
            else:
                
                if len(gen)==1:
                    
                    eachline+=sep+snpcodedict[gen[0]]+sep+snpcodedict[gen[0]]
                else:
                   # try:
                    
                    eachline+=sep+snpcodedict[gen[0]]+sep+snpcodedict[gen[1]]
                   # except:
                    #    print IUPACMatrix[SNP,geno],gen
                    #    raw_input


   
        eachline+="\n"
        
        towrite+=eachline
    with open(outf,"a") as outfile:
        outfile.write(towrite)
    print "file written in %s" %outf   




def IUPAC2binarybygen(header,SNPlist,IUPACMatrix,sep, outf):


    towrite="SNP"+sep
    towrite+=sep.join([each+sep+each for each in header])
    towrite+="\n"
    for SNP in SNPlist:
        
        SNPs=set()
        for geno in header:
            
            gen=reverse_iupac(IUPACMatrix[SNP,geno])
            SNPs|=set(gen)
        SNPs=list(SNPs)
        eachline=SNP
        for geno in header:
            gen=IUPACMatrix[SNP,geno]
            if gen=="NA":
                eachline+=sep+"NA"+sep+"NA"
            else:
                
                if gen in ["W","R","M","K","S","Y"]:
                    
                    eachline+=sep+"1"+sep+"1"
                elif gen==SNPs[0]:
                    eachline+=sep+"1"+sep+"0"
                elif gen==SNPs[1]:
                    eachline+=sep+"0"+sep+"1"
        eachline+="\n"
##        eachline+=SNP+sep+SNPs[1]
##        for geno in header:
##            gen=IUPACMatrix[SNP,geno]
##            if gen=="NA":
##                eachline+=sep+"NA"
##            else:
##                
##                if gen in["W","R","M","K","S","Y"]:
##                    
##                    eachline+=sep+"1"
##                elif gen==SNPs[0]:
##                    eachline+=sep+"0"
##                elif gen==SNPs[1]:
##                    eachline+=sep+"1"
##        eachline+="\n"
        towrite+=eachline
    with open(outf,"a") as outfile:
        outfile.write(towrite)
    print "file written in %s" %outf   

def removeMissing(header,SNPlist,IUPACMatrix, mono):
    newSNPlist=list()
    snpinfo={"miss":0,"mono":0}
    if mono:
        for SNP in SNPlist:
            eachline=SNP
            gencall=list()
            for geno in header:
                gencall.append(IUPACMatrix[SNP,geno])
            gencallU=list(set(gencall))
            if len(gencallU)==1:
                if "NA" in gencallU:
                    snpinfo["miss"]=snpinfo.get("miss")+1
                else:
                    snpinfo["mono"]=snpinfo.get("mono")+1
            elif len(gencallU)==2:
                if "NA" in gencallU:
                    snpinfo["mono"]=snpinfo.get("mono")+1
                else:
                    newSNPlist.append(SNP)
            else:
                newSNPlist.append(SNP)
        print "%s SNPs read \n%s SNPs with all missing and \n%s monomorphic SNPs removed " %(len(SNPlist), snpinfo["miss"], snpinfo["mono"])

    else:
        for SNP in SNPlist:
            eachline=SNP
            gencall=list()
            for geno in header:
                gencall.append(IUPACMatrix[SNP,geno])
            gencallU=list(set(gencall))
            if len(gencallU)==1:
                if "NA" in gencallU:
                    snpinfo["miss"]=snpinfo.get("miss")+1
                else:
                    newSNPlist.append(SNP)
           
            else:
                newSNPlist.append(SNP)
        print "%s SNPs read \n%s SNPs with all missing removed " %(len(SNPlist), snpinfo["miss"])    
    return header, newSNPlist, IUPACMatrix        



def iupac(listt):

    if len(set(listt))==1:
        return listt[0]

    elif "A" in listt and "T" in listt:
        return "W"
    elif "A" in listt and "C" in listt:
        return "M"
    elif "A" in listt and "G" in listt:
        return "R"
    elif "C" in listt and "T" in listt:
        return "Y"
    elif "C" in listt and "G" in listt:
        return "S"
    elif "G" in listt and "T" in listt:
        return "K"
    else:
        return "N"


def reverse_iupac(letter):

    if letter in ["A","C","G","T"]:
        return letter

    elif letter=="W":
        return ["A","T"]
    elif letter=="M":

        return ["A","C"]
    elif letter=="R":

        return ["A","G"]
    elif letter=="Y":

        return ["C","T"]
    elif letter=="S":

        return ["C","G"]
    elif letter=="K":

        return ["G","T"]
    else:
        return []


def reverse_iupacNA(letter):

    if letter in ["A","C","G","T"]:
        return letter

    elif letter=="W":
        return ["A","T"]
    elif letter=="M":

        return ["A","C"]
    elif letter=="R":

        return ["A","G"]
    elif letter=="Y":

        return ["C","T"]
    elif letter=="S":

        return ["C","G"]
    elif letter=="K":

        return ["G","T"]
    else:
        return ["NA"]



main()
