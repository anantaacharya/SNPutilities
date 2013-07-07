#ananta acharya
#ananta@uga,edu
#add barcode in front of fasta or fastq file
#use as interactive, open with idle, input values for prompts


def main():
    p=params()
    if p[0]=="2":
        fqtofq(*p[1::])
    elif p[0]=="1":
        fastatofq(*p[1::])
    elif p[0]=="3":
        fastatofasta(*p[1::])
    else:
        print "not valid parameter"

def params():
    formatt=raw_input("change format from, to :options are \n1. fasta,fastq \n2. fastq,fastq(fq) \n3. fasta, fasta: ").strip()
    infile=raw_input("input file").strip().strip('"')
    outfile=raw_input("output file").strip().strip('"')
    

    seqlen=int(raw_input("length of final sequence , if short reads detected A will be padded to end: ").strip())


    bc=raw_input("sequence to add in front, ie barcode + some enzyme cutsite, or only barcode : ").strip().upper()

    return(formatt, infile, outfile, seqlen, bc)
    


def fqtofq(infile, outfile, seqlen, bc):
    wbuffer=""
    with open(infile, "r") as inf:
        linecount=1        
        for line in inf:
            
            if linecount==1:
                wbuffer+=line
                linecount=2
            elif linecount==2:
                seq=line.strip().upper()
                tlen=len(seq)+len(bc)
                if tlen<seqlen:
                    seqtowrite=bc+seq+"A"*(seqlen-tlen)+"\n"
                else:
                    seqtowrite=bc+seq[0:seqlen-len(bc)]+"\n"
                wbuffer+=seqtowrite
                linecount=3
            elif linecount==3:
                linecount=4
                wbuffer+=line
            elif linecount==4:
                linecount=1
                qual=line.strip()
                tlen=len(seq)+len(bc)
                if tlen<seqlen:
                    qualtowrite="F"*len(bc)+qual+"F"*(seqlen-tlen)+"\n"
                else:
                    qualtowrite="F"*len(bc)+qual[0:seqlen-len(bc)]+"\n"
                wbuffer+=qualtowrite
    with open(outfile,"w") as outf:
        outf.write(wbuffer)


def fastatofq(infile, outfile, seqlen, bc):
    wbuffer=""
    with open(infile, "r") as inf:
        linecount=1        
        for line in inf:
            
            if linecount==1:
                linecount=2
                header="@somemachine"+line.strip()[1::]+"\n"
                
                
            elif linecount==2:
                linecount=1
                seq=line.strip().upper()
                tlen=len(seq)+len(bc)
                if tlen<seqlen:
                    seqtowrite=bc+seq+"A"*(seqlen-tlen)+"\n"
                else:
                    seqtowrite=bc+seq[0:seqlen-len(bc)]+"\n"
                thirdline="+\n"
                qualline="F"*seqlen+"\n"
                
                wbuffer+=header+seqtowrite+thirdline+qualline
                
           
    with open(outfile,"w") as outf:
        outf.write(wbuffer)               
            
    
def fastatofasta(infile, outfile, seqlen, bc):
    wbuffer=""
    with open(infile, "r") as inf:
        linecount=1        
        for line in inf:
            
            if linecount==1:
                linecount=2
                header=line
                
                
            elif linecount==2:
                linecount=1
                seq=line.strip().upper()
                tlen=len(seq)+len(bc)
                if tlen<seqlen:
                    seqtowrite=bc+seq+"A"*(seqlen-tlen)+"\n"
                else:
                    seqtowrite=bc+seq[0:seqlen-len(bc)]+"\n"
                                
                wbuffer+=header+seqtowrite
                
           
    with open(outfile,"w") as outf:
        outf.write(wbuffer)           

                
main()    
