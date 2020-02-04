from Bio import SeqIO

out = open('/home/stephen/Институт биоинформатики/Научный проет (Семестр 1)/HIV FASTQ/Untitled Folder/SRR5888726_filtered.fasta', 'w')

for record in SeqIO.parse('/home/stephen/Институт биоинформатики/Научный проет (Семестр 1)/HIV FASTQ/Untitled Folder/SRR5888726_revcomp.fasta', "fasta"):
    if record.id != "SRR5888726:::230216:0:0:0:":
        out.write( ">"+record.id+"\n"+str(record.seq)+"\n" )
out.close()