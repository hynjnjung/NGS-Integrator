# NGS-Integrator
NGS Integrator Tool (NGS-Integrator) is a Java-based tool that integrates multiple genome-wide NGS data via the minimum Bayes' Factor from density of sequence reads mapped on the genome. To integrate multiple NGS data, homogeneous replicates or 
heterogeneous NGS data can be integrated into the single data track based on the complements of the minimum Bayes factor (*cMBF*)
(range: 0-1) calculated from signal-to-background noise ratios as a function of genomic position.

This is **NOT** a peak-calling algorithm as it does not specifically locate peaks 
(though you can use a peak-calling tools on the output!).

### Installation
You can download the .jar file from the Epithelial Systems Biology Laboratory website [here](https://hpcwebapps.cit.nih.gov/ESBL/NGS-Integrator/). Move the .jar file to your working folder, and follow the run steps below!

**Requires: Java 1.7 or higher**


## Part 1: Calculation of the complement of minimum Bayes' Factor (Calculator)
The first NGS-Integrator program **calculates** the cMBF for each position for later integration. 
It first splits a given genome/chromosome file into a directory of separate chromosomes to facilitate parallelizing computation. 
For each chromosome file, the program uses a "moving window" to calculate the cMBF.

### Input file
NGS-Integrator takes as input a *sorted* BED file of sequence read counts (generated standard tools, such as *samtools* (depth), *bedtools* (coverage, genomecov), etc.)
The input file (bed file format) should contain the following (whitespace/tab separated) columns, and each region should be adjacent within the chromosome.
`chr  start end read_count`

For example,

```
chr1	0	10	0
chr1	10	20	1
chr1	20	30	1
```

The commandline options for changing the parameters are as below:
```
-h,--help
-i,--input <arg>         [req] input file path, must be a bed file
-m,--medMult <arg>       [opt] median multiple, cannot be 0                     default: 1
-o,--output <arg>        [opt] output file base name                            default: "out"
-w,--windowSize <arg>    [req] window size for calculating cMBF, must be        default: 10000 (bp)
                          a multiple of interval size
-z,--defaultZero <arg>   [opt] default number to replace zero                   default: 0.5
```

### Run
`java -jar NGS-Integrator.jar -i [Input.bed] -o [out_BaseName] -w [10000] [other options]`

### Output
This part produces a directory of chromosome-separated BED files of the calculated cMBFs.
The files have the following tab-separated format, `chr start end cMBF` (see the [UCSC bedGraph Format](https://genome.ucsc.edu/goldenpath/help/bedgraph.html)).
The directory also includes a final concatenated BED file of all chromosomes in the order of the given file.



## Part Two: Integration (Integrator)
The second part of the NGS-Integrator program **integrates** the calculated cMBFs for the different tracks 
by multiplying the cMBF at each position.

### Input file
Inputs should be Part One generated files, or have the same format (`chr start end cMBF`) (BED files).
An output filename also needs to included as the **first** input.
The inputs should be ordered as follows:
`chr start end cMBF`

### Run
`java -cp NGS-Integrator.jar Integrator [Integration_Output_File_Name.bed] [file1.bed] [file2.bed] ...`

### Output
The output will be a BED file of the same format (`chr start end cMBF`) and of the specified name.



## Things to Consider
- Window size: Larger is more accurate, but will take longer and requires more memory. 
- Median multiple: Generally, a larger median multiple will reduce the background level, but may also reduce signal.
- Default zero: 
- If you have multiple separated regions of the same chromosome, please place them into separate files to run the Calculator (Part One).
- NGS-Integrator requires open memory (generally at least 8GB, depending on window size). Running NGS-Integrator when low on space may cause errors such as `Improper file formatting` error.


## Questions/Issues?
Please contact us [here](https://esbl.nhlbi.nih.gov/contact.html) for any questions or bug reporting.

