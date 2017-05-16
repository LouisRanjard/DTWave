# DTWave
Dynamic Time Warping for spectrogram Alignment and AVErage sequence computation


## Introduction

DTWave is a tool for analysing sound sequences, it allows to 
- measure a pairwise distance between two or more sounds,
- compute an average sequence from a pair of sounds.

DTWave computes an edit distance between a pair of sounds which is proportional to the number of operation to transform one sequence into the other one (see references).

Currently, DTWave reads and writes HTK file format (http://htk.eng.cam.ac.uk/) but future releases will be format independent.



## Examples

DTWave is run from the command line, so first open a Terminal and cd (change directory) to the DTWave directory.

Two wav files are provided and have been encoded using HTK (http://htk.eng.cam.ac.uk/) and the config file "config" with the following commands:  
HCopy -C ./htk_config/config ./tui1.wav ./test_files/tui1.htk  
HCopy -C ./htk_config/config ./tui2.wav ./test_files/tui2.htk  
Or using a shellscript command:  
for i in ./*.wav; do htkname=\`echo $i | sed s:.wav:.htk:g\`; HCopy -C ./htk_config/config $i $htkname; done  

Here are some examples of the usage of DTWave, the command is given first ($) and the program output follows (>).

##### Measure a distance between two encoded sound files:  
$./dtwave_unix64 /fullpath_to/test_files/tui1.htk /fullpath_to/test_files/tui2.htk  
read htk file /fullpath_to/test_files/tui1.htk  
read htk file /fullpath_to/test_files/tui2.htk  
17.334454  
Using compression/expansion:  
$ ./dtwave_unix64 -c /fullpath_to/test_files/tui1.htk /fullpath_to/test_files/tui2.htk  
read htk file /fullpath_to/test_files/tui1.htk  
read htk file /fullpath_to/test_files/tui2.htk  
9.962879

##### Compute an average sequence of two encoded sound files:  
the option -w allows to give a weight toward one of the sequence for computing the average. The value must be between 0 and 1, for example -w 0.2 will pull the average toward the second file "tui2.htk", here the average is aimed to be at equal distance from each file (note: because of the way the average is computed it will not fall exactly at the same distance from the two input files, see reference papers):  
$./dtwave_unix64 -w .5 /fullpath_to/test_files/tui1.htk /fullpath_to/test_files/tui2.htk /fullpath_to/test_files/average.htk  
read htk file /fullpath_to/test_files/tui1.htk  
read htk file /fullpath_to/test_files/tui2.htk  
17.334454

##### Measure the distances between the original files and the average:  
$./dtwave_unix64 /fullpath_to/test_files/tui1.htk /fullpath_to/test_files/average.htk  
read htk file /fullpath_to/test_files/tui1.htk  
read htk file /fullpath_to/test_files/average.htk  
8.476674  
$./dtwave_unix64 /fullpath_to/test_files/tui2.htk /fullpath_to/test_files/average.htk  
read htk file /fullpath_to/test_files/tui2.htk  
read htk file /fullpath_to/test_files/average.htk  
9.424361 

##### Analyse a directory containing several encoded sound files (compute a pairwise distance matrix and an average sequence using compression/expansion):   
The program option -a can be ignored and the program will only compute the pairwise distance matrix.  
$./dtwave_unix64 -a /fullpath_to/test_files/  
3 files  
read htk file tui1.htk  
read htk file average.htk  
read htk file tui1.htk  
read htk file tui2.htk  
read htk file average.htk  
read htk file tui2.htk  
8.476674 17.334454 9.424397 

The average sequence is saved as the "average" file in the same directory.  
The pairwise distance matrix is saved in the file "pairwise_distances.csv" that should contains:  
tui1.htk, , 8.476674, 17.334454,   
average.htk, , , 9.424397,   
tui2.htk, , , ,   		



##  Notes

DTWave cannot read HTK compressed encoded sound files (config file with SAVECOMPRESSED=T).



##  Contacts and copyright

DTWave is provided as is, with no warranty. The developer and any other associated parties, provide no warranty for the use of DTWave. While every effort is taken to provide a reliable and trustworthy application, the authors, and associated parties, accept no liability for any damages, loss, or inconvenience caused by the use of DTWave. Comments, suggestions, queries, and bug reports are encouraged though - please feel free to email the author (louis.ranjard@gmail.com).



## References:

To cite the program, please refer to this paper:  
**Louis Ranjard, Michael G. Anderson, Matt J. Rayner, Robert B. Payne, Ian McLean, James V. Briskie, Howard A. Ross, Dianne H. Brunton, Sarah M. N. Woolley, Mark E. Hauber, Bioacoustic distances between the begging calls of brood parasites and their host species: a comparison of metrics and techniques, Behavioral Ecology and Sociobiology, 2010. http://dx.doi.org/10.1007/s00265-010-1065-2**

A technical description of the method is available in:  
**Louis Ranjard, Howard A. Ross, Unsupervised bird song syllable classification using evolving neural networks, Journal of the Acoustical Society of America 123(6):4358-4368, 2008**



--
Copyright (c) 2008-2017, Louis Ranjard

