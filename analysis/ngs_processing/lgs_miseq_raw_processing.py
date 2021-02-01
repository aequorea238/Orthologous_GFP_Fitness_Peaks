import numpy as np

def revcomp(seq):
    """Returns the reverse complement of a DNA sequence"""
    revcomp_seq = ''
    compdict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N':'N'}
    for nt in seq[::-1]:
        revcomp_seq += compdict[nt]
    return revcomp_seq


def get_barcode(seq, start_pattern, bc_length):
    """ Generic barcode-extractor for barcodes flanked by a constant sequence"""
    if seq.count(start_pattern) >= 1:
        barcode = seq[seq.index(start_pattern)+len(start_pattern):seq.index(start_pattern)+len(start_pattern)+bc_length]
        return barcode

    
class MiSeq_Read:
    """ This class exists for the preliminary processing of individual reads obtained from MiSeq'd samples. 
    Each Read has its associated forward and reverse reads, corresponding qualities, barcode, etc."""
    
    def __init__(self, name, r1_seq, r1_quality, r2_seq, r2_quality):
        self.name = name
        self.r1_seq = r1_seq
        self.r1_quality = r1_quality
        self.r2_seq = r2_seq
        self.r2_quality = r2_quality

        self.fastq_quality_code = '!"#$%&\'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~'
        
        
    def trim_extras(self, r1_trim, r2_trim):
        """ This will remove everything before and after the gene/barcode,
        as well as trim the corresponding r1&r2 qualities accordingly. """
        
        self.r1_quality = self.r1_quality[self.r1_seq.index(r1_trim)+len(r1_trim):]
        self.r1_seq = self.r1_seq[self.r1_seq.index(r1_trim)+len(r1_trim):]
        self.r2_quality = self.r2_quality[self.r2_seq.index(r2_trim)+len(r2_trim):]
        self.r2_seq = self.r2_seq[self.r2_seq.index(r2_trim)+len(r2_trim):]
        
        self.trimmed = True
        
        
    def revcomp_read2(self):
        """ This creates the reverse complement of read2 and its corresponding quality"""
        
        r2_rc = revcomp(self.r2_seq)
        self.r2_rc_seq = r2_rc
        self.r2_rc_quality = self.r2_quality[::-1]
        
        
    def get_bc(self, half):
        """ This locates the barcode based on constant sequences immediately flanking it, which should be 
        identical in all clones. The argument should be a string, 'bc_on_r1' or 'bc_on_r2" depending on which read 
        the barcode is on."""
                
        if half == 'n':            
            self.bc_seq = self.r1_seq[0:20]
            self.bc_quality = self.r1_quality[0:20]
            
            
        elif half == 'c':
            self.bc_seq = revcomp(self.r2_seq[0:20])
            self.bc_quality = self.r2_quality[0:20]
            self.bc_quality = self.bc_quality[::-1]
           
      
        

class MiSeq_SubSample:
    """ This class exists for preliminary processing of each individual sample sent for MiSeq.
    Namely: discard the reads without barcodes, trim away the N-shifts, and merge together the read1+read2 into one molecule.
    It takes the fastq file as input. Gene ID must be specified, as well as N- or C-terminal half."""
    
    def __init__(self, fastq, gene, half):
        self.reads = []
        self.gene = gene
        self.half = half
        self.reads_are_trimmed = False
        
        with open(fastq, 'r') as f:
            while True:
                x = f.readline().strip()
                name = x
                x = f.readline().strip()
                r1_seq = x
                x = f.readline()
                x = f.readline().strip()
                r1_quality = x
                x = f.readline()
                x = f.readline().strip()
                r2_seq = x
                x = f.readline()
                x = f.readline().strip()
                r2_quality = x
                self.reads.append(MiSeq_Read(name=name, r1_seq=r1_seq, r1_quality=r1_quality,
                                             r2_seq=r2_seq, r2_quality=r2_quality))
                if not x:
                    break
                
        print ('This sample has a total of %s paired-end reads.' % len(self.reads))
        
        # Different genes/halves have different sequences flanking the barcode
        if self.gene == 'amac':
            if self.half == 'c':
                self.r1_pattern = 'CGACACACTG'
                self.r2_pattern = 'GTCTCAACCT'
            elif self.half == 'n':
                self.r1_pattern = 'CCTAGGAACA'
                self.r2_pattern = 'GCTCGATGCG'
        
        if self.gene == 'cgre':
            if self.half == 'c':
                self.r1_pattern = 'GCCAGTACGA'
                self.r2_pattern = 'GTCTCAACCT'
            elif self.half == 'n':
                self.r1_pattern = 'CCTAGGAACA'
                self.r2_pattern = 'CCGTTTGACT'
                
        if self.gene == 'pplu':
            if self.half == 'c':
                self.r1_pattern = 'GAGGACGGCG'
                self.r2_pattern = 'GTCTCAACCT'
            elif self.half == 'n':
                self.r1_pattern = 'CCTAGGAACA'
                self.r2_pattern = 'CCTTGAAGTC'
                
        
        
    def remove_barcodeless(self):
        """ This will check whether the constant sequences which surround the barcode are present.
        If they are not, then the read is removed from the library."""

        self.reads = [read for read in self.reads if read.r2_seq.count(self.r2_pattern) >= 1 and read.r1_seq.count(self.r1_pattern) >= 1]
        print('There are now %s paired-end reads remaining, after removing the ones without a barcode.' % len(self.reads))
            
        
        
    def trim_extras(self):
        """ This trims away primer sequences and N-shifts, and ensures reads meet the expected minimal length afterwards."""
        for read in self.reads:
            read.trim_extras(self.r1_pattern, self.r2_pattern)
        self.reads_are_trimmed = True
        
        self.reads = [read for read in self.reads if len(read.r1_seq)>=270 and len(read.r2_seq)>=265]
        print("Primer/non-gene/non-bc sequences have been trimmed from this sample's reads. %s reads remain after filtering away weirdly truncated reads." % len(self.reads))
        
    
        
    def get_barcodes(self):
        """ This extracts the barcode from each read."""
        self.bc_counts = {}
        
        if self.reads_are_trimmed == False:
            print('WARNING! Trim me first or no barcodes for you.')

        else:
            for read in self.reads:
                read.get_bc(self.half)
                    
                if read.bc_seq in self.bc_counts and len(read.bc_seq) == 20:
                    self.bc_counts[read.bc_seq] += 1
                elif len(read.bc_seq) == 20:
                    self.bc_counts[read.bc_seq] = 1
                   
            self.reads = [read for read in self.reads if len(read.bc_seq)==20]
            print('Barcodes have been extracted for this sample. There are %s unique barcodes in %s reads.' % (len(self.bc_counts), len(self.reads)))
            
            

            
def get_consensus(keepers_dict, n=5, x=0.8):
    """ Takes a dictionary for input, where each key is a unique barcode and its value is the list of all reads associated with that barcode. Barcodes are discarded if they don't meet the 'n' minimal number of reads. For any given position, a fraction of reads greater than 'x' must agree, or the posision is deemed ambiguous. Returns a dictionary, where keys are barcodes and values are the corresponding consensus sequences. This process should be done independently for forward and reverse reads of N- and C-terminal halves of each gene."""
    consensus_dict = {bc:'' for bc in keepers_dict}
    
    for bc in keepers_dict:
        counts = {'A':[0]*n, 'T':[0]*n, 'C':[0]*n, 'G':[0]*n, 'N':[0]*n}

        for sq in keepers_dict[bc]:
            for i,nt in enumerate(sq[:n]):
                counts[nt][i] += 1
        
        consensus = ''
        for i in range (0,n):
            maxnt = counts[max(counts, key=lambda key: counts[key][i])][i]
            a = counts['A'][i]
            c = counts['C'][i]
            t = counts['T'][i]
            g = counts['G'][i]
            
            if a + c + t + g + counts['N'][i] < len(keepers_dict[bc]):
                break
            elif x=='r2':
                if maxnt==a==c or maxnt==a==t or maxnt==a==g or maxnt==t==g or maxnt==t==c or maxnt==c==g:
                    consensus += 'N'
                else:
                    consensus += max(counts, key=lambda key: counts[key][i])
            elif maxnt < len(keepers_dict[bc])*x:
                consensus += 'N'
            else:
                consensus += max(counts, key=lambda key: counts[key][i])

        consensus_dict[bc] = consensus
        
    return consensus_dict


def r1r2_merger(data, bcs, r1, r2, overlap):
    """This merges the forward and reverse reads for a given gene half. Takes as input a dataframe where the indeces are barcodes, the forward and reverse reads (or gene halves) are separate columns, and 'overlap' is the expected length of sequence overlap between the parts to be merged, which should match."""
    print('Initially, we have %s barcodes associated to their r1 and r2 consensus sequences.' % len(bcs))
    merged = {}
    conflicting = {}
    
    for bc in bcs:
        if data.loc[bc][r2][:overlap] in data.loc[bc][r1]:
            merged[bc] = data.loc[bc][r1][:data.loc[bc][r1].index(data.loc[bc][r2][:overlap])] + data.loc[bc][r2]
        
        
        #if data.loc[bc][r1][len(data.loc[bc][r1])-overlap:] == data.loc[bc][r2][:overlap]:
         #   merged[bc] = data.loc[bc][r1] + data.loc[bc][r2][overlap:]
        else:
            conflicting[bc] = [data.loc[bc][r1], data.loc[bc][r2]]
    
    print('R1 and R2 were successfully merged for %s barcodes.' % len(merged))
    return [merged, conflicting]