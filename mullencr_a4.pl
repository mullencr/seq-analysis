#! /usr/bin/perl
use strict;
use Bio::SeqIO;
use Align;
use SeqAnalysis;

my $gene_file = shift;
my $query_seq = shift;
my $match = shift;
my $mismatch = shift;
my $gap = shift;
my $max_error = shift;
my $mode = shift;
open(FOUT, ">alignment.txt");

# Use bio::seqIO to loop through the gene file and work out the alignment info.
my $stream = Bio::SeqIO->new(-file => $gene_file, -format => 'fasta');
while (my $seq = $stream->next_seq()) {
	# ====================================================== Run alignments ===================================================
	print "\n>".($seq->display_id())." ".($seq->desc())."\n\n";
	my $min_len = length($query_seq) - $max_error; 
	my $align = Align->new(-target_name=>$seq->display_id(), -seq1=>$query_seq, -seq2=>$seq->seq(), -match=>$match, -mismatch=>$mismatch, -gap=>$gap, -min_map_len=>$min_len, -max_error=>$max_error, -mode=>$mode);
	
	$align->printSeqWithSpacer();
	
	print "\n>".($seq->display_id())." ".($seq->desc())."\n\n";
	$align->getAlignment();
	
	# Write to the file. If I try to do this in the getAlignment method, it overwrites the file.
	select(FOUT);
	print "\n>".($seq->display_id())." ".($seq->desc())."\n\n";
	$align->getAlignment();
	select(STDOUT);
	print "\n";
	
	# =================================================== Run Seq Analysis ==========================================
	my $seqan = SeqAnalysis->new(-seq_name=>$seq->display_id(), -sequence=>$seq->seq());
    # Functionality 1
    my ($a, $t, $g, $c, $n) = $seqan->nucleotideCounter();
    print "[1] Nucleotide Counts: A=$a T=$t G=$g C=$c Other=$n\n";
    # Functionality 2 & 3
    my ($GCcontent, $SeqLength) = $seqan->gcContentSeqLength();
    my $GCcontentPercentage = ($GCcontent / $SeqLength);
    print "[2] GC Content: $GCcontentPercentage\n";
    print "[3] Sequence length: $SeqLength\n";
    # Funcitonality 4
    $seqan->detectEnzyme();
    # Functionality 5
    print "[5] Dinucleotide Frequency (%)\n\n";
    my %dinucleotide_hash = ();
    init_dhash(\%dinucleotide_hash);
    $seqan->dinucleotideFrequency(\%dinucleotide_hash);
    my $count = 0;
	foreach my $dinuc (sort keys %dinucleotide_hash) {
        $count++;
        print substr("[$dinuc]=$dinucleotide_hash{$dinuc}", 0, 9)." ";
        if($count % 4 == 0) {
                print "\n";
        }
    }
    # Functionality 6
    $seqan->detectPolyaSignal();
    # Funcitonality 7
    print "[7] Codon Usage:\n\n";
    my %codon_dist = ();
    $seqan->codonUsage(\%codon_dist);
    $count = 0;
    foreach my $codon (sort keys %codon_dist) {
        $count++;
        print substr("[$codon]=$codon_dist{$codon}", 0, 11)." ";
        if($count % 4 == 0) {
                print "\n";
        }
    }
    print "\n\n";
    # Functionality 8
    my $woLabel_ref = $seqan->detectMotifs('GAATCC', 'GAATGG', 'GAACCCC');
    my %woLabel = %{$woLabel_ref};
    my $withLabel_ref = $seqan->detectMotifsWithLables('A1'=>'GAATCC', 'A2'=>'GAATGG', 'A3'=>'GAACCCC');
    my %withLabel = %{$withLabel_ref};
    print "[8] Detection of motifs without labels: ('GAATCC', 'GAATGG', 'GAACCCC')\n\nMotif\tStart\tEnd\n";
    foreach my $motif (sort keys %woLabel) {
            my ($start, $end) = split(/:/, $woLabel{$motif});
            print "$motif\t$start\t$end\n";
    }
    print "\n[9] Detection of motifs with labels: ('A1'=>'GAATCC', 'A2'=>'GAATGG', 'A3'=>'GAACCCC')\n\nLabel\tMotif\tStart\tEnd\n";
    foreach my $name (sort keys %withLabel) {
            my ($start, $end) = split(/:/, $withLabel{$name});
            my ($label, $motif) = split(/:/, $name);
            print "$label\t$motif\t$start\t$end\n";
    }
}

sub init_dhash {
        my $hash_ref = shift;
        my @bases = ('A', 'T', 'G', 'C');
        for(my $i = 0; $i < 4; $i++) {
                for(my $c = 0; $c < 4; $c++) {
                        my $dinucleotide = $bases[$i].$bases[$c];
                        $$hash_ref{$dinucleotide} = 0;
                }
        }
}
