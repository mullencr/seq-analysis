package Match;
use strict;
sub new {
  my $class = shift;
  my %args=@_;
  my $self = bless {}, $class;
  foreach my $key (keys %args)
  {
    my $value=$args{$key};
    $self->{$key}=$value;
  }
  return $self;
}
1;