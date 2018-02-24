package Align;
use strict;
use Match;

# Attributes:
#my $AlignObject=Align->new(-seq1=>$seqbase1,
 #						  -seq2=>$seqbase2,
  #                                    -match=>$match,
   #                                   -mismatch=>$mismatch,
    #                                  -gap=>$gap,
#  -min_map_len=>$min_map_len,
#  -max_error=>$max_error,
#  -mode=>$mode # mode=all or one
#);   
      
#$AlignObject->printSeqWithSpacer();
#$AlignObject->getAlignment();

sub new {
  my $class = shift;
  my %args=@_;
  my $self = bless {}, $class;
  foreach my $key (keys %args)
  {
    my $value=$args{$key};
    $self->{$key}=$value;
  }
  
  # Create an array of match objects, it can be empty.
  my @matches = ();
  $self->{'matches_array'} = \@matches;
  
  # Populate a class variable with all of the scores from the alignment.
  $self->getScoreMatrix();
  # Get the maximum score, store it in a variable.
  my $score = $self->getMaxScore();
  # If we don't get any alignments that meet our criteria, try it again with the next lowest score.
  my $num_good = 0;
  while($num_good <= 0 && $score >= 0) {
  	$self->getPairs($score);
  	$num_good += $self->traceBack();
  	if($num_good == 0) {
  		$score--;
  	}
  }
  $self->{'highest_found_score'} = $score;
  
  return $self;
}

sub getScoreMatrix {
  my $self=shift;
  my $seq1=$self->{'-seq1'};
  my $seq2=$self->{'-seq2'};
  # scoring scheme
  my $MATCH    =  $self->{'-match'};
  my $MISMATCH = $self->{'-mismatch'};
  my $GAP      = $self->{'-gap'};
  # initialization
  my @matrix;
  $matrix[0][0]{score}   = 0;
  $matrix[0][0]{pointer} = "none";
  for(my $j = 1; $j <= length($seq1); $j++) {
    $matrix[0][$j]{score}   = 0;
    $matrix[0][$j]{pointer} = "none";
  }
  for (my $i = 1; $i <= length($seq2); $i++) {
    $matrix[$i][0]{score}   = 0;
    $matrix[$i][0]{pointer} = "none";
  }

  for(my $i = 1; $i <= length($seq2); $i++) {
    for(my $j = 1; $j <= length($seq1); $j++) {
        my ($diagonal_score, $left_score, $up_score);
        
        # calculate match score
        my $letter1 = substr($seq1, $j-1, 1);
        my $letter2 = substr($seq2, $i-1, 1);       
        if ($letter1 eq $letter2) {
            $diagonal_score = $matrix[$i-1][$j-1]{score} + $MATCH;
        }
        else {
            $diagonal_score = $matrix[$i-1][$j-1]{score} + $MISMATCH;
        }
        
        # calculate gap scores
        $up_score   = $matrix[$i-1][$j]{score} + $GAP;
        $left_score = $matrix[$i][$j-1]{score} + $GAP;
        
        if ($diagonal_score <= 0 and $up_score <= 0 and $left_score <= 0) {
            $matrix[$i][$j]{score}   = 0;
            $matrix[$i][$j]{pointer} = "none";
            next; # terminate this iteration of the loop
        }
        
        # choose best score
        if ($diagonal_score >= $up_score) {
            if ($diagonal_score >= $left_score) {
                $matrix[$i][$j]{score}   = $diagonal_score;
                $matrix[$i][$j]{pointer} = "diagonal";
            }
            else {
                $matrix[$i][$j]{score}   = $left_score;
                $matrix[$i][$j]{pointer} = "left";
            }
        } else {
            if ($up_score >= $left_score) {
                $matrix[$i][$j]{score}   = $up_score;
                $matrix[$i][$j]{pointer} = "up";
            }
            else {
                $matrix[$i][$j]{score}   = $left_score;
                $matrix[$i][$j]{pointer} = "left";
            }
        }
    }
  }
  $self->{'score_matrix'} = \@matrix;
}

sub getMaxScore {
  my $self = shift;
  my $seq1=$self->{'-seq1'};
  my $seq2=$self->{'-seq2'};
  my @matrix = @{ $self->{'score_matrix'} };
  my $max_score = 0;
  
  for(my $i = 1; $i <= length($seq2); $i++) {
    for(my $j = 1; $j <= length($seq1); $j++) {
        # set maximum score
        if ($matrix[$i][$j]{score} > $max_score) {
            $max_score = $matrix[$i][$j]{score};
        }
    }
  }
  return $max_score;
}

sub getPairs {
  my $self = shift;
  my $max_score = shift;
  my $seq1=$self->{'-seq1'};
  my $seq2=$self->{'-seq2'};
  my @matrix = @{ $self->{'score_matrix'} };
  my @coordinates = ();
  
  for(my $i = 1; $i <= length($seq2); $i++) {
  	for(my $j = 1; $j <= length($seq1); $j++) {
        # Check for a match
        if ($matrix[$i][$j]{score} == $max_score) {
            push(@coordinates, "$i:$j");
        }
    }
  }
  $self->{'coordinates'} = \@coordinates;
}

sub traceBack {
  my $self=shift;
  my $seq1=$self->{'-seq1'};
  my $seq2=$self->{'-seq2'};
  my $min_len = $self->{'-min_map_len'};
  my $max_error = $self->{'-max_error'};
  my @matrix = @{ $self->{'score_matrix'} };
  my $coordinates_ref = $self->{'coordinates'};
  my @coordinates = @{$coordinates_ref};
  my $num_good = 0;
  
  if($#coordinates < 0) {
  	return 0;
  }
  foreach my $location (@coordinates) {
	  # trace-back
	  my $align1 = "";
	  my $align2 = "";
	  my $nucl1 = "";
	  my $nucl2 = "";
	  my $verti = "";
	  my $cigar = "";
  	  
	  my ($i, $j) = split(/:/, $location);
	  my $end = $i;
	  while (1) {
	    last if $matrix[$i][$j]{pointer} eq "none";
	    
	    if ($matrix[$i][$j]{pointer} eq "diagonal") {
	        $align1 .= substr($seq1, $j-1, 1);
	        $align2 .= substr($seq2, $i-1, 1);
	        $nucl1=substr($seq1, $j-1, 1);
	        $nucl2=substr($seq2, $i-1, 1);
	        if($nucl1 eq $nucl2){
	        	$verti.="|";
	        	$cigar.="m";
	        }else{
	        	$verti.=" ";
	        	$cigar.="M";
	        }
	        $nucl1="";
	        $nucl2="";
	        $i--; $j--;
	    }
	    elsif ($matrix[$i][$j]{pointer} eq "left") {
	        $align1 .= substr($seq1, $j-1, 1);
	        $align2 .= "-";
	        $verti.=" ";
	        $cigar.="D";
	        $j--;
	    }
	    elsif ($matrix[$i][$j]{pointer} eq "up") {
	        $align1 .= "-";
	        $align2 .= substr($seq2, $i-1, 1);
	        $verti.=" ";
	        $cigar.="I";
	        $i--;
	    }   
	  }
	  $align1 = reverse $align1;
	  $align2 = reverse $align2;
	  $verti = reverse $verti;
	  $cigar = reverse $cigar;
	  my $start = $end - length($align1);
	  
	  # Check the match before we add it to the matches array.
	  my $good = 1;
	  # min_map_len
	  if(length($align1) < $min_len) {
	  	$good = 0;
	  }
	  # both ends valid
	  if(substr($verti, 0, 1) eq " " || substr($verti, length($verti)-1, 1) eq " ") {
	  	$good = 0;
	  }
	  # max errors not exceeded
	  my $err_count = 0;
	  for(my $i = 0; $i < length($verti); $i++) {
	  	if(substr($verti, $i, 1) ne "|") {
	  		$err_count++;
	  		if($err_count > $max_error) {
	  			$good = 0;
	  			last;
	  		}
	  	}
	  }
	  
	  if($good) {
		  my $match = Match->new(-align1 => $align1, -align2 => $align2, -verti => $verti, -cigar => $cigar, -start => $start, -end => $end);
		  push($self->{'matches_array'}, $match);
		  $num_good++;
	  }
  }
  return $num_good;
}

# ============================== Printing ======================================

sub printSeqWithSpacer {
	my $self = shift;
	$self->headers(" ", 100);
	$self->seqOutput(" ", 100);
}

# This function is designed to handle my headers, and is called by the output function.
sub headers {
		my $self = shift;
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
		my $self = shift;
        # First we need to unpack the parameters
        my $space = shift;
        my $line_len_num = shift;
        my $final_seq = $self->{'-seq2'};
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

# ====================== Print Aligns ===============================
sub getAlignment() {
	my $self = shift;
	my $matches_ref = $self->{'matches_array'};
	my @matches = @{$matches_ref};
	print "[Scoring schema]: match=$self->{'-match'}, mismatch=$self->{'-mismatch'}, gap=$self->{'-gap'}\n";
	print "[Search Target]: $self->{'-seq1'}\n";
	print "[Minimum target length]: ".(length($self->{'-seq1'}) - $self->{'-max_error'})."\n";
	print "[Maximum target length]: ".(length($self->{'-seq1'}) + $self->{'-max_error'})."\n";
	print "[Maximum allowed error bases]: ".($self->{'-max_error'})."\n";
	print "[The highest alignment score]: ".($self->{'highest_found_score'})."\n";
	print "[The alignments with the highest score]: ".($#matches+1);
	print "\n\n";
	if($#matches+1 <= 0) {
		print "No valid alignment is found for your criteria!\n\n";
	} else {
		for(my $i = 0; $i < (($self->{'-mode'} eq "one") ? 1 : $#matches+1); $i++) {
			my $match = $matches[$i];
			my $temp_cigar = $match->{'-cigar'};
			my $m = ($temp_cigar =~ tr/m/./);
			my $M = ($temp_cigar =~ tr/M/./);
			my $I = ($temp_cigar =~ tr/I/./);
			my $D = ($temp_cigar =~ tr/D/./); 
			print "[".($i+1)."] Seq$i vs. $self->{'-target_name'}\n\n";
			print " $match->{'-start'} $match->{'-align2'} $match->{'-end'} ";
			print "[Mapped-length=".(length($temp_cigar))." m=$m I=$I D=$D M=$M]\n";
			print " ".(" "x(length($match->{'-start'} + 1)))." $match->{'-verti'} ".(" "x(length($match->{'-end'} + 1)))." \n";
			print " ".(" "x(length($match->{'-start'} + 1)))." $match->{'-align1'} ".(" "x(length($match->{'-end'} + 1)))." \n";
			print " ".(" "x(length($match->{'-start'} + 1)))." $match->{'-cigar'} ".(" "x(length($match->{'-end'} + 1)))." \n\n";
		}
	}
}

1;