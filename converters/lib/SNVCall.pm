#!/usr/bin/perl

=head2 SNVCall.pm

 Module for using with maf-producing script

=cut

# SNVCall.pm 

package SNVCall;   

use constant DEBUG=>1;

sub new {
 my $class  = shift;  
 my $ch     = shift;
 my $start  = shift;
 my $ref    = shift;
 my $alt    = shift;
 $alt =~s/,.*//;
  
 $start++; # Actual mutated base is one base downstream

 # it is given by the user as an argument
 my $self = {chrom   => $ch,
             start   => $start,
             end     => $start,
             refa    => uc($ref),
             tuma    => uc($alt),
             # Next set is incrementally initialized
             tumorid => "NA",
             normid  => "NA",
             context => "NA",
             
             ref_f1r2 => 0,
             ref_f2r1 => 0,
             alt_f1r2 => 0,
             alt_f2r1 => 0,
             reads    => {},
             type     => "SNP",
             poxoq    => "0.00"
             };        # the internal structure we'll use to represent
 bless( $self, $class );    # make $self an object of class $class

 return $self;
}

=head3 Add read
 add_read needs read id and bit tag (needed for finding out the orientation)
=cut

sub add_read {
 my $self    = shift;
 my $read_id = shift;
 return if ($self->{reads}->{$read_id}); # Do not allow duplicates
 $self->{reads}->{$read_id}++;
}

=head3 Has read?
 do we have this read?
=cut

sub has_read {
 # reads will arrive with mate info)
 my $self    = shift;
 my $read    = shift;
 my $tag     = shift;

 if ($read !~m!(/[1,2]{1})$!) {
  my $mate  = $tag & 64 ? 1 : 2;
  $read.="/$mate";
 }

 return $self->{reads}->{$read} ? 1 : 0;
}

=head3 Update read
 given read id, tag and sequence find the orientation
 of this read and update our counts
=cut

sub update_read {
 # reads will arrive without mate info (or not)
 my $self    = shift;
 my $read    = shift;
 my $tag     = shift;
 my $chr     = shift;
 my $start   = shift;
 my $cig     = shift;
 my $seq     = shift; 

 if ($tag & 256) { # Get rid of secondary alignments
   $self->{reads}->{$read} = undef;
   return;
 }

 if ($read !~m!(/[1,2]{1})$!) {
  my $mate  = $tag & 64 ? 1 : 2;
  $read.="/$mate";
 }

 # Get the base from read that is in the SNP
 if ($self->{start} <= $start) {
   # This should not be frequent
   print STDERR "There is an error, read does not overlap the variant call\n" if DEBUG;
   $self->{reads}->{$read} = undef;
   return;
 }

 my $offset = $self->{start} - $start - 1; # Adjust for 0-based index
 # Checking for soft-clipped sequence
 if ($cig =~/^(\d+)S/) {
     $offset+=$1;  
 }

 my $seq_m = &parse_cig($seq, $cig);

 if ($offset >= length($seq_m)) {
   $self->{reads}->{$read} = undef;
   return;
 }

 my $readbase = uc(substr $seq_m, $offset, 1);
 if ($readbase ne $self->{refa} && $readbase ne $self->{tuma}) { 
   print STDERR "Readbase [$readbase] of read [$seq_m] with CIGAR [$cig] DOESN't match either REF [$self->{refa}] OR ALT [$self->{tuma}] allele\n" if DEBUG;
   # Discard these reads since they do not support either ref or alt
   $self->{reads}->{$read} = undef;
   return;  
 } 

 if ($tag & 64) {       # First read
   print STDERR "First in pair\n" if DEBUG;
   if ($tag & 16 && !($tag & 32)) {
      $readbase eq $self->{refa} ? $self->{ref_f2r1}++ : $self->{alt_f2r1}++;
   } elsif ($tag & 32 || !($tag & 16)) {
     $readbase eq $self->{refa} ? $self->{ref_f1r2}++ : $self->{alt_f1r2}++
   }
 } elsif ($tag & 128) { # Second read
   print STDERR "Second in pair\n" if DEBUG;
   if ($tag & 16 && !($tag & 32)) {
     $readbase eq $self->{refa} ? $self->{ref_f1r2}++ : $self->{alt_f1r2}++;
   } elsif ($tag & 32 || !($tag & 16)) {
     $readbase eq $self->{refa} ? $self->{ref_f2r1}++ : $self->{alt_f2r1}++
   }
 } 
}

=head3
 set_tumor_id - set the name of tumor
 get_tumor_id - get the name of tumor
=cut

sub set_tumor_id {
 my $self = shift;
 my $id   = shift;
 $self->{tumorid} = $id;
}

sub get_tumor_id {
 my $self = shift;
 return $self->{tumorid} eq "NA" ? "TUMOR" : $self->{tumorid};
}

=head3
 set_normal_id - set the name of ref sample
 get_normal_id - get the name of ref sample
=cut

sub set_normal_id {
 my $self = shift;
 my $id   = shift;
 $self->{normid} = $id;
}

sub get_normal_id {
 my $self = shift;
 return $self->{normid} eq "NA" ? "NORMAL" : $self->{normid};
}

=head3
 set_poxoq - set the poxo-q value
 get_poxoq - get the poxo-q value
=cut

sub set_poxoq {
 my $self   = shift;
 my $poxo   = shift;
 $self->{poxoq} = $poxo;
}

sub get_poxoq {
 my $self = shift;
 return $self->{poxoq};
}


=head3
 set_context - set the nucleotide context string
 get_context - get the nucleotide context string
=cut

sub set_context {
 my $self = shift;
 my $cont   = shift;
 $self->{context} = uc($cont);
}

sub get_context {
 my $self = shift;
 my $nucs = shift;
 unless($nuc){ return $self->{context}; }

 #Use $nuc to extract middle $nuc nucleotides
 if ($nuc < length($self->{context})) {
  return substr $self->{context}, (length($self->{context}) - $nuc)/2, $nuc;
 }

 return $self->{context};
}

=head3 
 foxog function returns Foxog, coefficient described in original paper's methods
=cut

sub foxog {
 my $self = shift;

 if ($self->{refa} eq "C" || $self->{refa} eq "A") {
   return $self->{alt_f1r2} > 0 || $self->{alt_f2r1} > 0 ? sprintf '%.2f',$self->{alt_f2r1}/($self->{alt_f1r2} + $self->{alt_f2r1}) : 0;
 } else { # G or T
   return $self->{alt_f1r2} > 0 || $self->{alt_f2r1} > 0 ? sprintf '%.2f',$self->{alt_f1r2}/($self->{alt_f1r2} + $self->{alt_f2r1}) : 0;
 }
}

=head3 to_sting
 Return a tab-delimited set of data for this object
=cut

sub to_string {
 my $self = shift;
 return join("\t",($self->{chrom},$self->{start},$self->{end},$self->{refa},$self->{tuma},
                   $self->{tumorid},$self->{normid},$self->{context},
                   $self->{alt_f1r2},$self->{alt_f2r1},$self->{ref_f1r2},$self->{ref_f2r1},
                   $self->foxog(),$self->{type}),$self->get_poxoq());
}

=head3 parse_cig
 A 'private function for modifying split sequences (with mismatches etc)

 M alignment match (can be a sequence match or mismatch)
 I insertion to the reference
 D deletion from the reference
 N skipped region from the reference
 S soft clipping (clipped sequences present in SEQ)
 H hard clipping (clipped sequences NOT present in SEQ)
 P padding (silent deletion from padded reference)
 = sequence match
 X sequence mismatch
=cut

sub parse_cig {

 my $sm;
 my ($s,$c) = @_;
 my @chars = split("",$s);

 CHUNK:
 while($c =~/(\d+)(\D{1})/g) {
  my $count = $1;
  my $op    = $2;
  # We fully handle only M I D N = and X flags
  # and do nothing about the others (just adding back bases)

  # Insertions:
  if ($op eq "I") {
   map {shift @chars;} (1..$count);
   next CHUNK;
  }

  # Deletions or splits:
  if ($op eq "D" || $op eq "N") {
   map {$sm.="N";} (1..$count);
   next CHUNK;
  }

  # Anything else:
  for (1..$count) {
    my $next_char = @chars > 0 ? shift @chars : "N";
    $sm.=$next_char;
  }
 }

 return $sm;
}

1;
