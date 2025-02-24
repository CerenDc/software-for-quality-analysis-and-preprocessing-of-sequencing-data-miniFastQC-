from Bio import SeqIO
import gzip
import matplotlib.pyplot as plt
import numpy as np
import os

os.getcwd()
records1 = SeqIO.parse('SRR800768_1_sub.fastq', 'fastq')
records2 = SeqIO.parse('SRR800768_2_sub.fastq', 'fastq')

for record in records1:
  seq = record.seq
  phred_scores = record.letter_annotations['phred_quality']
  print(seq)
  print(phred_scores)
  break

records1 = SeqIO.parse('SRR800768_1_sub.fastq', 'fastq')
count = 0   #Question1
for record in records1:
  count+=1

print(count)

count = 0 
for record in records2:
  count+=1

print(count)

for (rec1, rec2) in zip(records1, records2):
  if rec1.name != rec2.name:
    print("ERROR") 

L1 = [] #Question2
for rec in records1:
  L1.append(len(rec.seq))

set(L1)
print(L1)

L2 = [] 
for rec in records2:
  L2.append(len(rec.seq))

set(L2)
print(L2)

#Question3
Dico1 = {'A': [0] * 101, 'T': [0] * 101, 'C': [0] * 101, 'G': [0] * 101}

for record in SeqIO.parse('SRR800768_1_sub.fastq', "fastq"):
    seq = str(record.seq)
    for i, base in enumerate(seq[:101]):  
        if base in Dico1:
            Dico1[base][i] += 1

FreqA1 = Dico1['A']
FreqT1 = Dico1['T']
FreqC1 = Dico1['C']
FreqG1 = Dico1['G']

for i in range(101):
   FreqA1[i] = (FreqA1[i]/100000)*100
   FreqT1[i] = (FreqT1[i]/100000)*100
   FreqC1[i] = (FreqC1[i]/100000)*100
   FreqG1[i] = (FreqG1[i]/100000)*100
   
x = []
for i in range(101):
   x.append(i)

plt.plot(x, FreqA1, color='red', label='%A')
plt.plot(x, FreqT1, color='green', label='%T')
plt.plot(x, FreqC1, color='blue', label='%C')
plt.plot(x, FreqG1, color='black', label='%G')

plt.xlabel('Position in Read')
plt.ylabel('Base Percentage (%)')
plt.title("Base Composition per Read Position for SRR800768_1_sub.fastq")
plt.legend()
plt.show()

Dico2 = {'A': [0] * 101, 'T': [0] * 101, 'C': [0] * 101, 'G': [0] * 101}

for record in SeqIO.parse('SRR800768_2_sub.fastq', "fastq"):
    seq = str(record.seq)
    for i, base in enumerate(seq[:101]):  
        if base in Dico2:
            Dico2[base][i] += 1

FreqA2 = Dico2['A']
FreqT2 = Dico2['T']
FreqC2 = Dico2['C']
FreqG2 = Dico2['G']

for i in range(101):
   FreqA2[i] = (FreqA2[i]/100000)*100
   FreqT2[i] = (FreqT2[i]/100000)*100
   FreqC2[i] = (FreqC2[i]/100000)*100
   FreqG2[i] = (FreqG2[i]/100000)*100
   
x = []
for i in range(101):
   x.append(i)

plt.plot(x, FreqA2, color='red', label='%A')
plt.plot(x, FreqT2, color='green', label='%T')
plt.plot(x, FreqC2, color='blue', label='%C')
plt.plot(x, FreqG2, color='black', label='%G')

plt.xlabel('Position in Read')
plt.ylabel('Base Percentage (%)')
plt.title("Base Composition per Read Position for SRR800768_2_sub.fastq")
plt.legend()
plt.show()

#Question 5 

countA1 = 0
countT1 = 0
countC1 = 0
countG1 = 0

for i in range(101):
   countA1 = countA1 + Dico1['A'][i]
   countT1 = countT1 + Dico1['T'][i]
   countC1 = countC1 + Dico1['C'][i]
   countG1 = countG1 + Dico1['G'][i]

gc_content1 = ((countG1 + countC1) / (countA1 + countT1 + countG1 + countC1)) * 100 

countA2 = 0
countT2 = 0
countC2 = 0
countG2 = 0

for i in range(101):
   countA2 = countA2 + Dico2['A'][i]
   countT2 = countT2 + Dico2['T'][i]
   countC2 = countC2 + Dico2['C'][i]
   countG2 = countG2 + Dico2['G'][i]

gc_content2 = ((countG2 + countC2) / (countA2 + countT2 + countG2 + countC2)) * 100 

#Question 6

Dico_phred1 = {i: [] for i in range(101)} 

for record in SeqIO.parse("SRR800768_1_sub.fastq", "fastq"):
    phred_scores = record.letter_annotations["phred_quality"]
    for i, score in enumerate(phred_scores[:101]):
        Dico_phred1[i].append(score)

q1_1 = [np.percentile(Dico_phred1[i], 25) if Dico_phred1[i] else 0 for i in range(101)]
median1 = [np.median(Dico_phred1[i]) if Dico_phred1[i] else 0 for i in range(101)]
q3_1 = [np.percentile(Dico_phred1[i], 75) if Dico_phred1[i] else 0 for i in range(101)]

plt.plot(x, median1, color="purple", label="Median")
plt.plot(x, q1_1, color="blue", label="Q1")
plt.plot(x, q3_1, color="red", label="Q3")

plt.xlabel("Position dans le read")
plt.ylabel("Score Phred")
plt.title("Phred Score Distribution per Read Position for SRR800768_1_sub.fastq")
plt.legend()
plt.show()

Dico_phred2 = {i: [] for i in range(101)}

for record in SeqIO.parse("SRR800768_2_sub.fastq", "fastq"):
    phred_scores = record.letter_annotations["phred_quality"]
    for i, score in enumerate(phred_scores[:101]):
        Dico_phred2[i].append(score)

q1_2 = [np.percentile(Dico_phred2[i], 25) if Dico_phred2[i] else 0 for i in range(101)]
median2 = [np.median(Dico_phred2[i]) if Dico_phred2[i] else 0 for i in range(101)]
q3_2 = [np.percentile(Dico_phred2[i], 75) if Dico_phred2[i] else 0 for i in range(101)]

plt.plot(x, median2, color="purple", label="Median")
plt.plot(x, q1_2, color="blue", label="Q1")
plt.plot(x, q3_2, color="red", label="Q3")

plt.xlabel("Position dans le read")
plt.ylabel("Score Phred")
plt.title("Phred Score Distribution per Read Position for SRR800768_2_sub.fastq")
plt.legend()
plt.show() 

#Question 8 

def trim_low_quality_bases(record, quality_threshold=20):
    """ Trims low-quality bases (Phred score < threshold) from both ends of a read. """
    phred_scores = record.letter_annotations["phred_quality"]
    start = 0
    while start < len(phred_scores) and phred_scores[start] < quality_threshold:
        start += 1
    end = len(phred_scores) - 1
    while end > start and phred_scores[end] < quality_threshold:
        end -= 1
    return record[start:end+1]

def clean_paired_reads(records1, records2, output_paired1, output_paired2, output_single1, output_single2, quality_threshold=20, min_length=30):
    """ Cleans paired-end reads by trimming low-quality bases and filtering short reads. """
    paired1, paired2 = [], [] 
    single1, single2 = [], []  
    records1 = list(SeqIO.parse(records1, "fastq"))
    records2 = list(SeqIO.parse(records2, "fastq"))

    for rec1, rec2 in zip(records1, records2):
        trimmed_rec1 = trim_low_quality_bases(rec1, quality_threshold)
        trimmed_rec2 = trim_low_quality_bases(rec2, quality_threshold)
        len1, len2 = len(trimmed_rec1.seq), len(trimmed_rec2.seq)
        if len1 >= min_length and len2 >= min_length:
            paired1.append(trimmed_rec1)
            paired2.append(trimmed_rec2)
        elif len1 >= min_length:  
            single1.append(trimmed_rec1)
        elif len2 >= min_length:  
            single2.append(trimmed_rec2)

    with open(output_paired1, 'w') as f1, open(output_paired2, 'w') as f2:
        SeqIO.write(paired1, f1, "fastq")
        SeqIO.write(paired2, f2, "fastq")
    with open(output_single1, 'w') as f1, open(output_single2, 'w') as f2:
        SeqIO.write(single1, f1, "fastq")
        SeqIO.write(single2, f2, "fastq")

    print(f"Processed {records1} & {records2}:")
    print(f"  {len(paired1)} paired reads saved to {output_paired1}")
    print(f"  {len(paired2)} paired reads saved to {output_paired2}")
    print(f"  {len(single1)} single reads saved to {output_single1}")
    print(f"  {len(single2)} single reads saved to {output_single2}")

records1 = "SRR800768_1_sub.fastq"
records2 = "SRR800768_2_sub.fastq"

output_paired1 = "SRR800768_1_sub_clean.fastq"
output_paired2 = "SRR800768_2_sub_clean.fastq"
output_single1 = "SRR800768_1_sub_sing_clean.fastq"
output_single2 = "SRR800768_2_sub_sing_clean.fastq"

clean_paired_reads(records1, records2, output_paired1, output_paired2, output_single1, output_single2)