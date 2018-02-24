# seq-analysis
A set of Perl modules designed to analyze a set of genetic sequences and return meta information on it. This program contains two different classes for analyzing sequences, Align.pm and SeqAnalysis.pm.

***Align.pm*** This is a module which implements the Smith/Waterman algorithm, as provided by Dr. Chun Liang. The algorithm is designed to take one sequence and attempt to find similar structures in a target sequence by matching up the bases. It allows for mismatches and gaps in the target sequence or the query sequence, and is very useful for tracking mutations between different species.

***SeqAnalysis:*** This class is designed to take a single sequence in the constructor, after which a list of methods can be called to perform various sorts of analysis. The functions that can be called are as follows:
* ***codonUsage:*** Creates a histogram of all of the codons that appear in the sequence being analyzed.
* ***detectPolyaSignal:*** Just as the name implies, the method uses regex pattern matching to detect a Poly-A tail in the sequence.
* ***dinucleotideFrequency:*** Create a histogram of all of the codons that appear in the sequence being analyzed.
* ***detectEnzyme:*** This is a method with potentially useful applications, though it's implemented statically rather than dynamically. Detect a list of predefined restriction enzymes within the sequence using regex pattern matching. Supports IUPAC codes
* ***gcContentSeqLength:*** GC content proportion in a sequence is often deterministic of other qualities. This method returns GC coutnt and seq length separately, as per the specifications of the assignment
* ***nucleotideCounter:*** As per the name, counts nucleotides in the sequence, returns them as separate values rather than a hash (as per specifications of the assignment).
* ***printWithSpacer:*** Print the sequence with wrapping adn spacing, with numbered rows and columns for ease of viewing
* ***detectMotifs/detectMotifsWithLabels:*** Search out and return the start and end positions of specified subsequences in a hash
* ***revCompl*** Returns the reverse compliment of a sequence, which is useful for the transcription process, a common operation in bioinformatics.

***mullencr_a4.pl:*** In this case, we apply the alignment algorithm to map a query sequence specified by the user onto an entire FASTA format sequence file, also specified by the user. If a match is found, the sequence analysis functions are called and information displayed to the user. The results of this script can be found in alignment.txt. To produce similar output, the following command can be run: `perl mullencr_a4.pl GeneSeq.fa GAAGGTTCCG 1 -1 -2 2 all`. The format for this command is: `perl mullencr_a4.pl <FASTA_format_sequence_file> <query_sequence> <match_score> <mismatch_score> <gap_score> <max_gaps> <all|one>`. The numerical parameters are all configurations for the algorithm, and determine how to score the alignments. The final parameter specifies wether or not to print all matches found, or only the highest scoring match.
