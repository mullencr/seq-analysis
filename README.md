# seq-analysis
A set of Perl modules designed to analyze a set of genetic sequences and return meta information on it. This program contains two different classes for analyzing sequences, GeneList.pm and SeqAnalysis.pm. The code mullencr_a3.pl invokes these classes as a demonstration.

***GeneList:*** This class is designed to extract gene distribution data from two different files, each with lists of genes and their expressions in certain groups. The doOutput method can be invoked to print a list of the genes that are unique to each file, as well as a list of the genes that appear in both and their comparative expressions.

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
