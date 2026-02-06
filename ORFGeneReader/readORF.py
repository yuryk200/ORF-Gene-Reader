import sys
import tkinter
from tkinter import filedialog

def read():
    
    #Make tkinter object to open file explorer to select fasta file
    root = tkinter.Tk()
    root.wm_withdraw()
    
    FileName = filedialog.askopenfilename(filetypes = [('All files','*.*')])
    root.destroy()
    
    FileContents = []
    
    # open the file in reading text mode using error checking
    try:
        
        # open the file with the name stored in FileName in read-only mode and assign the file object to Fp1
        Fp1 = open(FileName, 'r')
        sequence = Fp1.readline()
        
        Data = Fp1.read()
        
    except IOError:
        print("error unable to read file or file does not exist!!!")
        print("Exiting the program")
        Fp1.close()
        sys.exit(1)
    
    # split the contents of the file at each newline character and store in ListSeq
    ListSeq = Data.split('\n')
    DnaSeq = ('').join(ListSeq) 
    
    # append the value of sequence and DnaSeq to the end of the list FileContents
    FileContents.append(sequence)
    FileContents.append(DnaSeq)
    
    return FileContents


def Compliments(DnaSeq):
    
    ComplimentSeq = ''
    
    # loop through the characters in the DNA sequence and get the reverse complement
    for index in range(0, len(DnaSeq)):
        if DnaSeq[index] == 'T':
            ComplimentSeq += 'A'
        if DnaSeq[index] == 'A':
            ComplimentSeq += 'T'
        if DnaSeq[index] == 'C':
            ComplimentSeq += 'G'
        if DnaSeq[index] == 'G':
            ComplimentSeq += 'C'

    print("The primary sequence 5' to 3' \n", DnaSeq)
    print("\nThe compliment sequence 3' to 5' \n",ComplimentSeq)
    print("\nThe reverse compliment 5' to 3' \n", ComplimentSeq[::-1])

    input("")

    return ComplimentSeq


def Translate(DnaSeq, RFNumber):

    AminoAcidList = []
    AminoAcidSeq = ''

    # defining a dictionary of codons and their corresponding amino acids
    CodonTable = {
       'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
       'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
       'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
       'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
       'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
       'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
       'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
       'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
       'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
       'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
       'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
       'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
       'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
       'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
       'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
       'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W',
       }
    
    print("\n***************************************** reading frame number {:d} ********************************************".format(RFNumber+1))
    DnaSequenceRF = DnaSeq[RFNumber:len(DnaSeq)]
    print("the DNA seq is: \n")
    print(DnaSequenceRF)

    # iterating over every third base starting from the given reading frame number
    for n in range(RFNumber, len(DnaSeq), 3):
        
        # selects the next three bases as a codon
        codon = DnaSeq[n:n+3]
        
        if codon in CodonTable:
            
            # adding the corresponding amino acid to the amino acid sequence
            AminoAcid = CodonTable[codon]
            AminoAcidSeq += AminoAcid
    
    return AminoAcidSeq
       

def main():
    
    # gets the fasta file
    FileContents = read()
    DesLine = FileContents[0]
    DnaSeq = FileContents[1]
    
    # get reverse complement of DNA sequence
    Compseq = Compliments(DnaSeq)
    
    print("**************** The primary strand is: ")
    print(DnaSeq)
    print("**************** the reverse compliment strand is: ")
    print(Compseq)

    print("\n********************** All the codons of the PRIMARY STRAND ************************************")
    
    for RFNumber in range(0, 3):
        
        # translate DNA sequence into amino acid sequence
        AAseq1 = Translate(DnaSeq, RFNumber)
        print("\nthe amino acid sequence of RF {:d} is: \n".format(RFNumber + 1))
        print(AAseq1) 

    print("\n********************** All the codons of the REVERSE COMPLIMENTARY STRAND ************************************")

    for RFNumber in range(0, 3):
        
        # translate reverse complement of DNA sequence into amino acid sequence
        AAseq2 = Translate(Compseq[::-1],  RFNumber)
        print("\nthe compliment amino acid sequence printed right to left of RF {:d} is: \n".format(RFNumber - 3))
        print(AAseq2[::-1])
       
    return 0



main()
        
        
        
    