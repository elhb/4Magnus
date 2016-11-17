def comp(sequence):
    '''Function that takes a sequence and complements it'''
    
    # define translation dictionary
    complement = {
                'A':'T',
                'T':'A',
                'C':'G',
                'G':'C',
                'N':'N',
                'R':'Y',
                'Y':'R',
                'K':'M',
                'M':'K',
                'B':'V',
                'V':'B',
                'D':'H',
                'H':'D',
                }
    
    # translate each nucleotide in the sequence
    complementary_sequence = "".join(
                        [ complement.get( nucleotide.upper(), '') for nucleotide in sequence ]
                    )
    
    return complementary_sequence

def revcomp(sequence):
    '''Function that takes a sequence and reversecomplements it'''
    
    # use the comp function to get complement
    complementary = comp(sequence)
    
    # reverse the sequence
    reverse_complementary = complementary[::-1]
    
    return reverse_complementary

def convert_fastqs(in_file_1,in_file_2):
    
    """
    This function reads a pair of gzipped fastqfiles and yields strings on fasta format
    the fasta format consist of the fastq read header and then the read one sequence followed by a stretch of N bases and then the reverse complement of the read two sequence
    the function is totally unavare of base calling qualities so some kind of filtering based on this might be needed later on to mask or remove low quality base calls
    """
    
    # a looop that reads all lines in the infiles
    while True:
        
        
        try:
            # reading the read one file
            header_1   = in_file_1.readline().rstrip()
            sequence_1 = in_file_1.readline().rstrip()
            junk_1     = in_file_1.readline().rstrip()
            quality_1  = in_file_1.readline().rstrip()

            # reading the read two file
            header_2   = in_file_2.readline().rstrip()
            sequence_2 = in_file_2.readline().rstrip()
            junk_2     = in_file_2.readline().rstrip()
            quality_2  = in_file_2.readline().rstrip()
            
            # check if we reached the end of the file
            if '' in [header_1,sequence_1,junk_1,quality_1,header_2,sequence_2,junk_2,quality_2]: raise EOFError

        # if we reached the end of file break the loop
        except EOFError: break

        # use this many N bases between the read sequences
        length_of_N_stretch = 20
    
        # make the long sequence that consist of both reads
        long_sequence = sequence_1 + 'N'.join('' for i in xrange(length_of_N_stretch)) + revcomp(sequence_2)

        # yield the long sequence on fasta format and continue with next iteration in the loop
        yield '>'+header_1[1:]+'\n'+long_sequence

#
# imports
#
import sys
import gzip

#
# open the infiles
#
in_file_1 = gzip.open(sys.argv[1])
in_file_2 = gzip.open(sys.argv[2])

#
# convert the infiles and print to stdout
#
for long_sequence in convert_fastqs(in_file_1,in_file_2): print long_sequence