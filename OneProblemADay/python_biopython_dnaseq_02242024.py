##	credit: https://www.geeksforgeeks.org/dna-protein-python-3/
##	https://biopython.org/docs/dev/api/Bio.SeqUtils.html
##	https://blog.adnansiddiqi.me/python-for-bioinformatics-getting-started-with-sequence-analysis-in-python/

from Bio import SeqIO
#from Bio.Alphabet import generic_dna


from Bio.SeqUtils import gc_fraction


def firstSeqFun(seqIN):
  idx = 0
  for sequence in SeqIO.parse(seqIN,'fasta'):
    print(sequence.id)
    print(sequence.seq)
    print('No of nucleotides: {}'.format(len(sequence)))
    idx+=1

  print('Total Sequence found {}'.format(idx))

def GCcontent(seqIN):
  for sequence in SeqIO.parse(seqIN,'fasta'):
    print(sequence.id)
    print(sequence.seq)

    print(f"GC content of {sequence.id} : {gc_fraction(sequence.seq):.2f}")
    print(len(sequence.seq))

if __name__ == '__main__':
  #firstSeqFun('ls_orchid.fasta')
  firstSeqFun('Atg9b_sequenced_in_placenta.fa')
  GCcontent('Atg9b_sequenced_in_placenta.fa')
#  protein = translate('Atg9b_sequenced_in_placenta.fa')

