from Bio import SeqIO

out = open('/home/stephen/Институт биоинформатики/Научный проет (Семестр 1)/HIV FASTQ/Untitled Folder/Revcomped/SRR5888729_revcomp.fasta', 'w')

for record in SeqIO.parse('/home/stephen/Институт биоинформатики/Научный проет (Семестр 1)/HIV FASTQ/Untitled Folder/Default/SRR5888729.fasta', "fasta"):
    out.write( ">"+record.id+"\n"+str(record.seq.reverse_complement())+"\n" )
out.close()