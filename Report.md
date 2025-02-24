# Report - Exploration of paired fastq files

## Dowloading the files

1. The sequencing run was performed using the Illumina platform.
2. The instrument model used is Illumina HiSeq 2000.
3. The type of library used for sequencing is paired-end. This information can be found in the ‘Library Layout’ section. 
4. The type of molecules that have been sequenced are genomic DNA, as indicated by the library source. If the sequencing had been for RNA, the source would have been labelled as Transcript.
5. The extension of these files is .fastq.gz. This indicates that the files are in FASTQ format. The .gz part signifies that the files are compressed using Gzip to reduce their size.
6. Yes, these files should be different in content. The pair of reads comes from the same DNA fragment, but they correspond to opposite ends of the fragment. This means that read 1 and read 2 contain complementary but distinct sequences. The sequencing process involves inversion, so even though the two reads come from the same fragment, they provide different information.
7. Each file is 380 MB in size, but we are not able to open and read them directly due to their size. For the tutorial, we use smaller subfiles of 22 MB to facilitate analysis and processing.

## Analyzing the files

1. There are 100,000 reads in each of the two files.
2. The read length is 101 for each file.
3. ![Plot 1 - Reads per base](/Users/dincceren/Desktop/ProjetFatsqc/Reads_per_base1)
   ![Plot 2 - Reads per base](/Users/dincceren/Desktop/ProjetFatsqc/Reads_per_base2)

4. We can see that the quantities of A and T are more abundant than those of C and G, with peaks in certain places. 
5. We estimated the GC-content of the Saccharomyces cerevisiae genome at 37.2%, which is slightly lower than the value of 38.3% reported by the data available on the Harvard BioNumbers site ([source](https://bionumbers.hms.harvard.edu/bionumber.aspx?s=n&v=1&id=102126)).
6. ![Plot 1 - Phred_Score_Distribution](/Users/dincceren/Desktop/ProjetFatsqc/Phred_Score_Distribution1)
   ![Plot 2 - Phred_Score_Distribution](/Users/dincceren/Desktop/ProjetFatsqc/Phred_Score_Distribution2)

7. No, the quality of the bases is not uniform throughout the playback. Base quality is lower at the end of the read.  
At the start of sequencing, the results show better performance, and this decreases as the bases are read, indicating lower performance. This may be due to a number of factors, such as reduced polymerase efficiency during DNA strand elongation, enzyme fatigue as the chain is extended, and instrumental limitations in detecting bases at the end of the sequence.