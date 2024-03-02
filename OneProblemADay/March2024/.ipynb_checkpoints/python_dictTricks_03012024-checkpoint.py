##	credit: https://stackoverflow.com/questions/16819222/how-do-i-return-dictionary-keys-as-a-list-in-python



Problem: 
##	Today I am working on modifying some files, and it turns out that I need to use the powerful data struction is python: dictionary
##	In particular, I want to get the keys of a dictionary in list format:


##	How do I return dictionary keys as a list in Python?



Answers: 

##	It turns out that there are so many tricks to get this:

newdict = {1:0, 2:0, 3:0}
newdict.keys()



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

