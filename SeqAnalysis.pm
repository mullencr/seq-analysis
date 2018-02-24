package SeqAnalysis;
use strict;
# I have this for the polya tail to make sure I have no less than one mismatch.
use String::Similarity;

# Constructor for the seqAnalysis class.
# Pass a sequence to this class and the methods can return analysis on them.
sub new {
	my $class = shift;
	my $self = bless {}, $class;
	
	my %args = @_;
	foreach my $key (keys %args) {
		my $value = $args{$key};
		$self->{$key}=$value;
	}
	return $self;
}

# Create a histogram of the different codons that appear in the 
sub codonUsage {
	my $self = shift;
	my $seq = $self->{'-sequence'};
	my $rev = $self->revCompl();
	my $codon_distribution_ref = shift;
	my $codon_count = 0;
	for(my $i = 0; $i <= length($seq)-3; $i++) {
		$codon_count++;
		my $codon = substr($seq, $i, 3);
		if(!exists($$codon_distribution_ref{$codon})) {
			$$codon_distribution_ref{$codon} = 0;
		}
		$$codon_distribution_ref{$codon} += 1;
	}
	for(my $i = 0; $i <= length($rev)-3; $i++) {
		$codon_count++;
		my $codon = substr($rev, $i, 3);
		if(!exists($$codon_distribution_ref{$codon})) {
			$$codon_distribution_ref{$codon} = 0;
		}
		$$codon_distribution_ref{$codon} += 1;
	}
	my %codons = %{$codon_distribution_ref};
	foreach my $codon (sort keys %codons) {
		$$codon_distribution_ref{$codon} /= $codon_count;
	}
}

# Detect a poly-a tail
sub detectPolyaSignal {
	my $self = shift;
	my $seq = $self->{'-sequence'};
	print "\n[6] Detection of poly(A) signal (AATAAA):\n\n";
	print "No.\tStart\tEnd\tSignal\n";
	my $count = 0;
	for (my $i = 0; $i <= length($seq)-6; $i++) {
		my $target = substr($seq, $i, 6);
		# 5/6 or single mismatch is 0.8333 repeating. 
		if(($target =~ m/^A.T..A$/) && (similarity($target, "AATAAA")) > 0.8) {
			$count++;
			my $end = $i + 6;
			print "$count\t$i\t$end\t$target\n";
		}
	}
	print "\n";
}

# Create a histogram of all of the dinucleotides.
sub dinucleotideFrequency {
	my $self = shift;
	my $seq = $self->{'-sequence'}; 
	my $dinucleotide_hash_ref = shift;
	my %dinucleotide_hash = %{$dinucleotide_hash_ref};
	# For each spot in the sequence, look at each nucleotide.
	my $d_count = 0;
	for(my $i = 0; $i <= length($seq)-2; $i++) {
		$d_count++;
		foreach my $dinuc (sort keys %dinucleotide_hash) {
			if(substr($seq, $i, 2) eq $dinuc) {
				$$dinucleotide_hash_ref{$dinuc} += 1;
			}
		}
	}
	foreach my $dinuc (sort keys %dinucleotide_hash) {
		$$dinucleotide_hash_ref{$dinuc} /= $d_count;
	}
}

# Detect a selection of restriction enzymes within the sequence.
sub detectEnzyme {
	my $self = shift;
	my $seq = $self->{'-sequence'};
	print "[4] Restriction Sites: \n\n";
	print "Name\tPos\tSeq\tIUPAC\tALT\n";
	for(my $i = 0; $i <= length($seq)-8; $i++) {
		my $target = substr($seq, $i, 8);
		if(substr($target, 0, 7) eq 'ACAAGGG') {
			my $tmp = substr($target, 0, 7);
			print "BclI\t$i\t$tmp\tACAAGGG\tACAAGGG\n";
		}
		if($target =~ m/GG[AG]GCA[CT]T/) {
			print "BfmI\t$i\t$target\tGGRGCAYT\tGG[AG]GCA[CT]T\n";
		}
		if($target =~ m/AGG[ACGT]TTTA/) {
			print "EcoRI\t$i\t$target\tAGGNTTTA\tAGG[ACGT]TTTA\n";
		}
		# We do the substring here because if we just matched, we might be off by one or two indeces.
		if(substr($target, 0, 6) =~ m/TGGCC[CT]/) {
			my $tmp = substr($target, 0, 6);
			print "Cac8I\t$i\t$tmp\tTGGCCY\tTGGCC[CT]\n";
		}
	}
	print "\n";
}

# Return values to calculate gc proportions
sub gcContentSeqLength {
	my $self = shift;
	my $seq = $self->{'-sequence'};
	my ($GCcontent, $SeqLength) = (0, 0);
	my ($a, $t, $g, $c, $n) = $self->nucleotideCounter();
	$GCcontent = $g + $c;
	$SeqLength = length($seq);
	return ($GCcontent, $SeqLength);
}

# Straightforward
sub nucleotideCounter {
	my $self = shift;
	my $seq = $self->{'-sequence'};
	my ($a, $t, $g, $c, $n) = (0, 0, 0, 0, 0);
	$a = ($seq =~ s/A/a/g);
	$t = ($seq =~ s/T/t/g);
	$g = ($seq =~ s/G/g/g);
	$c = ($seq =~ s/C/c/g);
	$n = ($seq =~ s/N/n/g);
	if($n == "") {
		$n = "0";
	}
	return ($a, $t, $g, $c, $n);
}

sub printWithSpacer {
	my $self = shift;
	my $gene_name = $self->{'-seq_name'};
	my $gene_seq = $self->{'-sequence'};
	my $space = " ";
	my $line_len_num = 100;
	print ">$gene_name\n\n";
	headers($space, $line_len_num);
	# Working with seqArray as a return value allows us to use the same output function for both file types.
	seqOutput($space, $line_len_num, $gene_seq);
}

sub detectMotifs {
	my $self = shift;
	my $seq = $self->{'-sequence'};
	my @motifs = @_;
	my %hash = ();
	foreach my $motif (@motifs) {
		my $start = -1;
		for(my $i = 0; $i <= length($seq)-length($motif); $i++) {
			if(substr($seq, $i, length($motif)) eq $motif) {
				$start = $i;
				last;
			}
		}
		my $end = $start + length($motif);
		if($start != -1) {
			$hash{$motif} = "$start:$end";
		} else {
			$hash{$motif} = "0:0";
		}
	}
	return \%hash;
}

sub detectMotifsWithLables {
	my $self = shift;
	my $seq = $self->{'-sequence'};
	my %args = @_;
	my %output = ();
	foreach my $label (sort keys %args) {
		my $motif = $args{$label};
		my $start = -1;
		for(my $i = 0; $i <= length($seq)-length($motif); $i++) {
			if(substr($seq, $i, length($motif)) eq $motif) {
				$start = $i;
				last;
			}
		}
		my $end = $start + length($motif);
		if ($start != -1) {
			$output{"$label:$motif"} = "$start:$end";
		} else {
			$output{"$label:$motif"} = "0:0";
		}
	}
	return \%output;
}

# Return reverse compliment
sub revCompl {
  my $self=shift;
  my $seq=$self->{'-sequence'};
  $seq=~tr/gatcGATC/ctagCTAG/;
  return $seq;
}

#========================= Assignment 01 Funcitons ========================
sub headers {
	my $space = shift;
	my $line_len_num = shift;
	my $up_to = $line_len_num/10;
	my $line1 = "".(" " x 5)."";
	for(my $i = 1; $i <= $up_to; $i += 1) {
		$line1 = $line1.$space.(" " x 9)."$i";
	}
	print $line1."\n";
	
	my $line2 = "Line ";
	for(my $i = 0; $i < $up_to; $i += 1) {
		$line2 = $line2.$space."1234567890";
	}
	print $line2."\n";
	
	my $line3 = (" "x5).$space;
	for(my $i = 0; $i < $up_to; $i += 1) {
		$line3 = $line3."__________";
		if ($space eq " ") { 
			$line3 = $line3."_"; 
		}
	}
	print $line3."\n";
}

# This function loops through and writes out all of the lines in their proper format.
sub seqOutput {
	# First we need to unpack the parameters
	my $space = shift;
	my $line_len_num = shift;
	my $final_seq = shift;
	# Now we just need to print the sequence out.
	# The following loop should loop through 4 times if the length were 359 (359/100 = 3.59). I tested this
	for (my $i = 0; $i < (length($final_seq)/$line_len_num); $i += 1) {
		# Print the initial line start. 
		# We can't just do 3 spaces, it would misalign line 10 or anything else that is two-digits.
		print "".(" "x(4 - length($i + 1))).($i + 1)."|"; # For some reason I needed the empty quotes at the beginning, otherwise it doesn't work.
		# Print the sequence in chunks. We need to do a check each time to see if there is anything left.
		for (my $c = 0; $c < ($line_len_num/10); $c += 1) {
			# Normally this loops through 5 or 10 times.
			# Our current index is the number of lines before us * line length + number of blocks before us on this line * 10.
			# With 0 based counting, the number of lines before us should be equal to $i, and the number of blocks before us should be $c
			my $curr_ind = ($i * $line_len_num) + ($c * 10);
			print $space.substr($final_seq, $curr_ind, 10);
		}
		print "\n"; # end the line
	}
	print "\n";
}

1;