class HiSeq_Read:
    """ This class describes individual reads from HiSeq'd samples, with their corresponding primary and secondary barcodes, read quality, etc."""
    def __init__(self, name, sq, quality):
        self.name = name
        self.sq = sq
        self.quality = quality
        
        self.fastq_quality_code = '!"#$%&\'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~'    

        
    def get_barcodes_1and2(self, pattern, bc1_length, bc2_length):
        """ Extracts the primary (20N) and secondary (10N) barcodes from the read, preceding and following the constant region between the two."""
        bc1 = self.sq[self.sq.index(pattern)-bc1_length:self.sq.index(pattern)]
        self.bc1 = bc1
        self.bc1_quality = self.quality[self.sq.index(pattern)-bc1_length:self.sq.index(pattern)]
        
        bc2 = self.sq[self.sq.index(pattern)+len(pattern):self.sq.index(pattern)+len(pattern)+bc2_length]
        self.bc2 = bc2
        self.bc2_quality = self.quality[self.sq.index(pattern)+len(pattern):self.sq.index(pattern)+len(pattern)+bc2_length]


class Gate:
    """ This class exists to process individual HiSeq samples originating from a specific sorted gate. It is populated by HiSeq Reads and associated library-wide statistics (number of unique barcodes, etc.)"""
    def __init__(self, fastq):
        self.reads = []
        with open(fastq, 'r') as f:
            while True:
                name = f.readline().strip()
                sq = f.readline().strip()
                quality = f.readline()
                quality = f.readline().strip()
                                
                if not quality:
                    break
                    
                self.reads.append(HiSeq_Read(name, sq, quality))
                
        print ('This sample has %s reads.' % len(self.reads))
        
        
    def filter_barcodeless(self, pattern='AGGTGCTAG'):
        """ This discards reads lacking the constant sequence between primary and secondary barcodes"""
        self.reads = [read for read in self.reads if pattern in read.sq]
        print ('There are %s reads left after removing the ones without the constant inter-bc sequence.' % len(self.reads))
        
        
    def get_barcodes(self, pattern, bc1_length, bc2_length):
        """Extracts primary and secondary barcodes from all reads"""
        for read in self.reads:
            read.get_barcodes_1and2(pattern, bc1_length, bc2_length)
        print ('Primary and secondary barcodes have been extracted.')
        
          
    def get_unique_1bcs(self):
        self.reads = [read for read in self.reads if read.bc1 != '' and read.bc2 != '']
        print ('%s reads are left after removing the ones with no 1bc despite having constant sq.' % len(self.reads))

        self.unique_1bcs = list(set([read.bc1 for read in self.reads]))
        print ('There are %s unique primary barcodes in this sample.' % len(self.unique_1bcs))
        

    def get_bc1_counts(self):
        """ Counts the number of times a primary barcode appears in the Gate sample """
        self.bc1_totalcounts = {bc:0 for bc in self.unique_1bcs}
        for read in self.reads:
            self.bc1_totalcounts[read.bc1] += 1
            
        print ('Number of times a given primary barcode appears have been counted.')
        
        
    def get_bc1_bc2_dict(self):
        """ Associates each primary barcode to all its secondaries """
        self.bc1_to_bc2_all = {read.bc1:[] for read in self.reads}
        for read in self.reads:
            self.bc1_to_bc2_all[read.bc1].append(read.bc2)
        print ('Primary barcodes have been linked to their secondary barcodes.')