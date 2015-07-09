package Sam::Seq;

use warnings;
use strict;

use overload '""' => \&string;

use List::Util;
use Storable 'dclone';

use Verbose;

use Sam::Parser;
use Sam::Alignment 0.10 ':flags';

use Fastq::Seq;
use Fasta::Seq;

use constant {
    # :) scaling constant for frequency to phred conversion. The higher
    # the value, the more trust is put into frequency (proovread-1.01: 50)
    PROOVREAD_CONSTANT => 120,
};

our $VERSION = '1.0.0';




=head1 NAME

Sam::Seq.pm

=head1 DESCRIPTION

Class for handling sam reference sequences and its aligned reads.

=cut

=head1 SYNOPSIS

  use Sam::Seq;

  my $sp = Sam::Parser->new(fh => \*SAM);
  my %ss;

  while (my %h = $sp->next_header_line('SQ')) {
      $ss{$h{SN}} = Sam::Seq->new(id=>$h{SN}, len => $h{LN});
  }

  while (my $aln = $sp->next_aln) {
      $ss{$aln->rname}->add_aln($aln);
  }

=cut


##------------------------------------------------------------------------##


=head1 Class ATTRIBUTES

=cut

=head2 $V

Verbose messages are handled using the Verbose.pm module. To
 customize verbose message behaviour, overwrite the attribute with
 another Verbose object created with the Verbose module.

=cut

our $V = Verbose->new();

=head2 $BinSize [20]

=cut

our $BinSize = 20;

=head2 $MaxCoverage [50]

=cut

our $MaxCoverage = 50;

=head2 $PhredOffset [33]

=cut

our $PhredOffset = 33;

=head2 InDelTaboo [0.1]

Trim reads to prevent insertions/deletions within the first/last
 InDelTaboo fraction of the read. N=0 deactivates the feature.

=cut

our $InDelTaboo = 0.1;

=head2 InDelTabooLength [undef]

Trim reads to prevent insertions/deletions within the first/last
 InDelTabooLength bps of the read. N=0 deactivates the feature. If defined,
 superceeds relative InDelTaboo

=cut

our $InDelTabooLength = undef;

=head2 $Trim

Boolean. Deactivate trimming completely, including leading/trailing indels
 and InDelTaboo.

=cut

our $Trim = 1;

=head2 $MaxInsLength

Default 0, which deactivates the feature.

=cut

our $MaxInsLength = 0;

=head2 $FallbackPhred

Default 1. Used if alignment w/o qual is used in quality context

=cut

our $FallbackPhred = 1;

=head2 $RepCoverage

Default = 0, which deactivates the feature. If set triggers
filter_rep_region_alns for each region with higher coverage.

=cut

our $RepCoverage = 0;


=head2 $MinScore/$MinNScore/$MinNCScore

Minimum score/nscore/ncscore cutoffs for accoring filter_by_n/c/score()
functions;

=cut

our $MinScore = undef;
our $MinNScore = undef;
our $MinNCScore = undef;

# DEPRECATED
#=head2 %Freqs2phreds
#
#p = sqrt(f*120);
#
#=cut
#
# our %Freqs2phreds;
# @Freqs2phreds{0..31} = (map{int((($_/50)**(1/2)*50)+0.5)}(0..39));
# @Freqs2phreds{32..100} = (40)x69;
# our %Phreds2freqs = reverse  %Freqs2phreds;

=head1 Class METHODS

=cut

=head2 BinSize

Get/Set $Sam::Seq::BinSize. Default 20.

=cut

sub BinSize{
	my ($class, $size) = @_;
	$BinSize = $size if defined $size;
	return $BinSize;
}

=head2 MaxCoverage

Get/Set $Sam::Seq::MaxCoverage. Default 50.

=cut

sub MaxCoverage{
	my ($class, $cov) = @_;
	$MaxCoverage = $cov if defined $cov;
	return $MaxCoverage;
}

=head2 PhredOffset

Get/Set $Sam::Seq::PhredOffset. Default 33.

=cut

sub PhredOffset{
	my ($class, $offset) = @_;
	$PhredOffset = $offset if defined $offset;
	return $PhredOffset;
}

=head2 RepCoverage

Get/Set $Sam::Seq::RepCoverage. Default 0.

=cut

sub RepCoverage{
	my ($class, $cov) = @_;
	$RepCoverage = $cov if defined $cov;
	return $RepCoverage;
}

=head2 MinScore/MinNScore/MinNCScore

Get/Set $Sam::Seq::MinScore/MinNScore/MinNCScore. Default undef.

=cut

sub MinScore{
	my ($class, $score, $force) = @_;
	$MinScore = $score if defined($score) || $force;
	return $MinScore;
}

sub MinNScore{
	my ($class, $score, $force) = @_;
	$MinNScore = $score if defined($score) || $force;
	return $MinNScore;
}

sub MinNCScore{
	my ($class, $score, $force) = @_;
	$MinNCScore = $score if defined($score) || $force;
	return $MinNCScore;
}

=head2 InDelTaboo

Get/Set $Sam::Seq::InDelTaboo. Default 10.

=cut

sub InDelTaboo{
	my ($class, $indeltaboo) = @_;
	$InDelTaboo = $indeltaboo if defined $indeltaboo;
	return $InDelTaboo;
}

=head2 InDelTabooLength

Get/Set $Sam::Seq::InDelTabooLength. Default undef.

=cut

sub InDelTabooLength{
	my ($class, $indeltaboolength, $force) = @_;
	$InDelTabooLength = $indeltaboolength if defined $indeltaboolength || $force;
	return $InDelTabooLength;
}

=head2 Trim

Get/Set $Sam::Seq::Trim. Default TRUE.

=cut

sub Trim{
	my ($class, $trim) = @_;
	$Trim = $trim if defined $trim;
	return $Trim;
}

=head2 MaxInsLength

Get/Set $Sam::Seq::MaxInsLength. Default 0, which deactivates the feature.

=cut

sub MaxInsLength{
	my ($class, $max_ins_length) = @_;
	$MaxInsLength = $max_ins_length if defined $max_ins_length;
	return $MaxInsLength;
}

=head2 FallbackPhred

Get/Set $Sam::Seq::FallbackPhred. Default 1.

=cut

sub FallbackPhred{
	my ($class, $fbp) = @_;
	$FallbackPhred = $fbp if defined $fbp;
	return $FallbackPhred;
}

=head2 Freqs2phreds

Convert a LIST of frequencies to a LIST of phreds.

=cut

sub Freqs2phreds{
	my $class = shift;
	return map{
            my $p = int( sqrt( $_ * PROOVREAD_CONSTANT ) +.5);
            $p > 40 ? 40 : $p;
        }@_;
}


=head2 Phreds2freqs

Convert a LIST of Phreds to a LIST of frequencies.

=cut

sub Phreds2freqs{
	my $class = shift;
	return map{
            int((( $_**2 / PROOVREAD_CONSTANT ) * 100)+ .5) /100
	}@_;
}


=head2 Qual2phreds

Return the phred values of the quality string accorting to specified offset.

=cut

sub Qual2phreds{
    my $class = shift;
    return map{$_-=$PhredOffset}unpack("W*", $_[0]);
}


=head2 Hx

Takes a reference to an ARRAY of counts, converts the counts to probabilities
 and computes and returns the shannon entropy to describe its composition.
 Omits undef values.

NOTE: Not a Class Method

  # R
  Hx = function(x){
  	p = x/sum(x);
  	-(sum(sapply(p, function(pi){pi*log2(pi)})))
  }

=cut

sub Hx{
	my ($col) = @_;
	my $total = 0;
	my @states = grep{$_}@$col;
	$total += $_ for @states;
	my @Px = map{$_/$total}@states;
	my $Hx;
	$Hx -= $_ for map{$_ * (log($_)/log(2)) }@Px;
	return $Hx;
}


=head2 Trace2cigar

Compress a trace string (MMMDMMIIMMDMM) to a cigar string (3M1D2M2I2M1D2M)

=cut

sub Trace2cigar{
	my ($class, $trace) = @_;
        chomp($trace); # just to be safe

	my $cigar = '';
        my $spos = 0;

        while($trace =~ /(\w)(?!\g{1})/g){
            $cigar.= pos($trace)-$spos.$1;
            $spos = pos($trace);
        }

        # DEPRECATED: cannot handle submatches >32k as the entire
        # submatch is captured
        #	while($trace =~ m/(\w)(\g{1}*)/g){
        #		$cigar .= length($1.$2).$1;
        #	}

	return $cigar;
}


=head2 State_matrix

=cut

sub State_matrix{
        my $self = shift;
	my %p = (
		alns => [$self->alns],
		'length' => $self->len,
		states => $self->{_states},
		matrix => undef,
                ignore_coords => undef,
                qual_weighted => 0,
                use_ref_qual => 0,
		@_
	);

	die unless defined ($p{alns} and $p{'length'}) || $p{matrix};

	# state matrix
	my @S = $p{matrix}
		? @{$p{matrix}}
		: map{[]}1..$p{'length'};	# init state matrix

	# predefined states
	my %states = %{$p{states}};

        # add ref to state matrix
        if($p{use_ref_qual} && $self->ref){
            my @seq = split (//, $self->ref->seq);
            my @freqs = Sam::Seq->Phreds2freqs($self->ref->phreds);

            for(my $i=0; $i<@seq; $i++){
                # never add 0, if nothing more matches, a 0 count might be introduced
                next unless $freqs[$i];
                #ie Dumper($self->{_states}, "@seq", $self->ref) unless defined $self->{_states}{$seq[$i]};
                ($S[$i][$self->{_states}{$seq[$i]}])+= $freqs[$i];
            }
        }

        foreach my $aln (@{$p{alns}}){

		###################
		### prepare aln ###
		# get read seq
		my $seq = $aln->seq;
		my $orig_seq_length = length($seq);
		next unless $orig_seq_length > 50;

                my $qua = $aln->qual;

                if($qua eq "*"){ # replace missing qual with fallback qual
                    my $q = Fastq::Seq->Phreds2Char( [$FallbackPhred] , $PhredOffset );
                    $qua = $q x $orig_seq_length;
                }

		# get read cigar, eg 80M2D3M1IM4
		my @cigar = split(/(\d+)/,$aln->cigar);
		shift @cigar;

                if($cigar[1] eq 'S'){
                    # just move on in query, do nothing else
                    $seq = substr($seq, $cigar[0]);
                    $qua = substr($qua, $cigar[0]);
                    shift @cigar;
                    shift @cigar;
                }
                if($cigar[-1] eq 'S'){
                    $seq = substr($seq, 0, -$cigar[-2]);
                    $qua = substr($qua, 0, -$cigar[-2]);
                    pop @cigar;
                    pop @cigar;
                }
                if($cigar[1] eq 'H'){
                    shift @cigar;
                    shift @cigar;
                }
                if($cigar[-1] eq 'H'){
                    pop @cigar;
                    pop @cigar;
                }


		$V->exit("Empty Cigar") unless @cigar;

		# reference position
		my $rpos = $aln->pos-1;

		if($Trim){
			##################
			### InDelTaboo ###
			# this also removes leading/trailing InDels regardless of InDelTaboo
                    # trim head
			my $mc = 0;
			my $dc = 0;
			my $ic = 0;
			my $indeltaboolength = $InDelTabooLength ? $InDelTabooLength : int($orig_seq_length * $InDelTaboo + 0.5);
			for(my $i=0; $i<@cigar;$i+=2){
				if($cigar[$i+1] eq 'M'){
					if($mc + $ic + $cigar[$i] > $indeltaboolength){
						if($i){# there was something before this match
							# only cut before this match
							# trim cigar
							splice(@cigar, 0, $i);
							# adjust rpos
							$rpos+= ($mc+$dc);
							# trim seq
							substr($seq, 0, $mc+$ic, '');
							substr($qua, 0, $mc+$ic, '');
						}
						last;
					}
					$mc+=$cigar[$i];
				}elsif($cigar[$i+1] eq 'D'){
					$dc+= $cigar[$i];
				}elsif($cigar[$i+1] eq 'I'){
					$ic+= $cigar[$i];
				}else{
					$V->exit("Unknown Cigar '".$cigar[$i+1]."'");
				}
			}

			# have to have kept at least 50 bps and 70% of original read length
			#  to consider read for state matrix
			next if length($seq) < 50  || (length($seq)/$orig_seq_length) < 0.7;

			# trim tail
			my $tail=0;
			for(my $i=$#cigar-1; $i;$i-=2){
				if($cigar[$i+1] eq 'M'){
					$tail+=$cigar[$i];
					if($tail > $indeltaboolength){
						if($i < $#cigar-1){# there is after this match
							# only cut before this match
							my $tail_cut = $tail-$cigar[$i];
							# trim cigar
							splice(@cigar, -($#cigar-($i+1)));
							# trim seq
							substr($seq, -$tail_cut, $tail_cut, '');
							substr($qua, -$tail_cut, $tail_cut, '');
						}
						last;
					}
				}elsif($cigar[$i+1] eq 'D'){ # ignore leading deletions, but adjust rpos

				}elsif($cigar[$i+1] eq 'I'){
					$tail+=$cigar[$i];
				}else{
					$V->exit("Unknown Cigar '".$cigar[$i+1]."'");
				}
			}

			# have to have kept at least 50 bps and 70% of original read length
			#  to consider read for state matrix
			next if length($seq) < 50  || (length($seq)/$orig_seq_length) < 0.7;
		}


		#######################
		### cigar to states ###
		my @states;
                my @squals;

		# cigar counter, increment by 2 to capture count and type of cigar (10,M) (3,I) (5,D) ...

		my $qpos = 0; # usually 0, >0 for cigar S
		for (my $i=0; $i<@cigar;$i+=2) {
                    if ($cigar[$i+1] eq 'M') {
                        push @states, split(//,substr($seq,$qpos,$cigar[$i]));
                        push @squals, split(//,substr($qua,$qpos,$cigar[$i])) if $p{qual_weighted};
                        $qpos += $cigar[$i];
                    } elsif ($cigar[$i+1] eq 'D') {
                        push @states, ('-') x $cigar[$i];
                        if($p{qual_weighted}){
                            my $qbefore = $qpos > 1 ? substr($qua,$qpos-1,1) : substr($qua,$qpos,1);
                            my $qafter  = $qpos < length($qua) ? substr($qua,$qpos,1) : substr($qua,$qpos-1,1);
                            my $q = $qbefore lt $qafter ? $qbefore : $qafter;
                            push @squals, ($q) x $cigar[$i];
                        };
                    } elsif ($cigar[$i+1] eq 'I') {
                        if ($i) {
                            # append to prev state
                            #print STDERR "@cigar\n" unless @states;
                            if ($states[$#states] eq '-') {
                                # some mappers, e.g. bowtie2 produce 1D1I instead of
                                # mismatchas (1M), as it is cheaper. This needs to be
                                # corrected to a MM
                                $states[$#states] = substr($seq,$qpos,$cigar[$i]);
                                $squals[$#states] = substr($qua,$qpos,$cigar[$i]) if $p{qual_weighted};

                            } else {
                                $states[$#states] .= substr($seq,$qpos,$cigar[$i]);
                                $squals[$#states] .= substr($qua,$qpos,$cigar[$i]) if $p{qual_weighted};
                            }
                        } else {
                            $states[0] = substr($seq,$qpos,$cigar[$i]);
                            $squals[0] = substr($qua,$qpos,$cigar[$i]) if $p{qual_weighted};
                        }
                        $qpos += $cigar[$i];
                    } else {
                        $V->exit("Unknown Cigar '".$cigar[$i+1]."'");
                    }
		}


		########################
		### states to matrix ###

		for(my $i=0; $i<@states; $i++) {
                    if ($p{ignore_coords} && _is_in_range($rpos, $p{ignore_coords})) {
                        $rpos++;
                        next;
                    }

                    my $state = $states[$i];

                    if (length ($state) > 1 && ! exists $states{$state}) {
                        $states{$state} = scalar keys %states;
                    }

                    if($p{qual_weighted}){
                        my $squal = $squals[$i];
                        my @phreds = Sam::Seq->Qual2phreds($squal);
                        my @freqs = Sam::Seq->Phreds2freqs(@phreds);
                        #print STDOUT "$squal @phreds @freqs\n";
                        my $freq = List::Util::min(@freqs);
                        ($S[$rpos][$states{$state}])+= $freq; # match/gap states always exist
                    }else{
                        ($S[$rpos][$states{$state}])++; # match/gap states always exist
                    }
                    $rpos++;
		}

	}

	return \@S, \%states;

}


##------------------------------------------------------------------------##

=head1 Constructor METHOD

=head2 new

Create a Sam::Seq object. Takes key => value representation.
Returns a Sam::Seq object.

  ## defaults
  id => undef,
  len => undef,           # length of the reference sequence, required
  ref => undef,           # reference seq object, Fasta::Seq or Fastq::Seq
  con => undef,           # consensus seq, Fastq::Seq (+cov)
  ref_merge_regions => [],# regions in the ref seq, to be included in
                          #  consensus calling, ARRAY of Tuples (start, offset)
  max_coverage => 50,     # assumed maximum coverage value for consensus
                          #  quality value calculation
  is => undef,            # see Sam::Parser->is()

  bin_size => $Sam::Seq::BinSize,
  bin_max_coverage => $Sam::Seq::BinMaxCoverage,
  phred_offset => $Sam::Seq::PhredOffset,


=cut

# alias for cloning
*clone = \&new;

sub new{
	my $proto = shift;
	my $self;
        my $class;

        if ($class = ref $proto) { # clone
            $proto->{_is} = undef; # can't clone code
            $self = bless(dclone($proto), $class);
            $self->is($self->{is}) if $self->{is};
            return $self;
        }

        $class = $proto;

	$self = {
		## defaults
		id => undef,
		len => undef,			# length of the reference sequence, required
		ref => undef,			# reference seq object, Fasta::Seq or Fastq::Seq
		con => undef,			# consensus seq, Fastq::Seq (+cov)
		max_coverage => $MaxCoverage,
		bin_size => $BinSize,
		bin_max_bases => $BinSize * $MaxCoverage,
		phred_offset => $PhredOffset,
		is => undef,
		## custom overwrites
		@_,
		## protected
		_is => sub{return 1},	# don't set directly, use is()
		_alns => {},			# reads or idx pos of reads in sam
		_bin_alns => undef,		# position bins containing aln refs
		_bin_scores => undef,	# position bins containing nscores of alns in _bin_alns
		_bin_bases => undef,
		_bin_lengths => undef,
		_aln_idc => 0,			#
		_state_matrix => [],	# consensus state matrix
		_states => {			# consensus states
			A => 0,
			T => 1,
			G => 2,
			C => 3,
			'-' => 4,
			N => 5,
			# .. complex states, dynamically added
		}

	};


	bless $self, $class;

	# prepare is test routine
	$self->is($self->{is}) if $self->{is};#

	# init bins
	$self->_init_read_bins;
	return $self;
}







##------------------------------------------------------------------------##

=head1 Public METHODS

=cut

=head2 add_aln_by_score

Add a sam alignment object of a mapped illumina read to the object, based on
 position and score. Takes a file position as second object. If provided this
 position is stored instead of the actual alignment, the position is assumed to
 be a index pointing to the alignment in the file indicated by C<< $sam_seq->sam
 >>. Returns internal alignment id (>0) if aln has been added, 0 if the
 alignment did not make the threshold or undef, if alignment couldn't be
 assessed, e.g. because it had no score/was unaligned.

  $sam_seq->add_aln_by_score($aln);
  $sam_seq->add_aln_by_score($aln, tell($sam_fh);


=cut

sub add_aln_by_score{
	my ($self, $aln) = @_;

	my $bin = $self->bin($aln);

	my $ncscore = $aln->ncscore;
        return undef unless defined $ncscore;

	# if bin_bases are full, check if new ncscore is good enough
	if( $self->{_bin_bases}[$bin] > $self->{bin_max_bases} ){
            # ignore scores, that are too low
            if( $ncscore <= $self->{_bin_scores}[$bin][-1] ){
                return 0;
            }else{ # sufficient score
                my $iid = $self->{_bin_alns}[$bin][-1];
                $self->remove_aln_by_iid($iid);
            }
	}

        my $bases = $aln->length;
        $self->{_bin_bases}[$bin] += $bases;
	my $id = $self->add_aln($aln);

	# set score/id at the right place
	my $i = @{$self->{_bin_scores}[$bin]} - 1;
	$i-- while $i >= 0 && $ncscore > $self->{_bin_scores}[$bin][$i];
	# store new  score and _id of aln at correct position
	splice(@{$self->{_bin_scores}[$bin]}, $i+1, 0, $ncscore);
	splice(@{$self->{_bin_alns}[$bin]}, $i+1, 0, $id);
	splice(@{$self->{_bin_lengths}[$bin]}, $i+1, 0, $bases);

	return $id;
}

=head2 add_aln

Add a Sam::Alignment to Sam::Seq.

=cut


sub add_aln{
	my ($self, $aln) = @_;

	$self->{_alns}{++$self->{_aln_idc}} = $aln;

	return $self->{_aln_idc};
}


=head2 remove_aln_by_iid

Remove a Sam::Alignment from Sam::Seq (including position bins, if it exits).
 Returns the removed Sam::Alignment or undef.

=cut

sub remove_aln_by_iid{
	my ($self, $id) = @_;

	my $aln = delete $self->{_alns}{$id};
	defined $aln || return;

	my $bin = $self->bin($aln);
        my $ba = $self->{_bin_alns}[$bin];
        if (@$ba) {

            my $idx = List::Util::first {$ba->[$_] == $id} 0..@$ba-1;
            defined $idx || return;

            splice(@{$self->{_bin_scores}[$bin]}, $idx, 1);
            splice(@{$self->{_bin_alns}[$bin]}, $idx, 1);
            my $rm_bases = splice(@{$self->{_bin_lengths}[$bin]}, $idx, 1);
            $self->{_bin_bases}[$bin] -= $rm_bases;
        }

	return $aln;
}


=head2 haplo_consensus

=cut

sub haplo_consensus{
	my $self = shift;
        my %p = (
                 hcrs => [],
                 ignore_coords => undef,
                 qual_weighted => 0,
                 use_ref_qual => 0,
                 @_
                );

	$self->_init_state_matrix(
                                  ignore_coords => $p{ignore_coords},
                                  qual_weighted => $p{qual_weighted},
                                  use_ref_qual => $p{use_ref_qual},
                                 );

	$self->_add_pre_calc_fq(@{$p{hcrs}}) if $p{hcrs};

        # get variants
        $self->variants(reuse_matrix => 1);
        $self->penalize_variants();
        my $hc = $self->haplo_coverage(reuse_matrix => 1);
        $self->filter_by_coverage($hc) if $hc;

        # recall
	$self->_init_state_matrix(
                                  ignore_coords => $p{ignore_coords},
                                  qual_weighted => $p{qual_weighted},
                                  use_ref_qual => $p{use_ref_qual},
                                 );
        $self->variants(
            min_prob => 0.2,
            min_freq => 2,
        );
	$self->_haplo_consensus;

	return $self->{con};

}


=head2 consensus

Calculate and the consensus sequence from state matrix. Returns a Fastq::Seq
 object, with an additional entry C<< $seq->{cov} >>, which contains phred
 like ascii coded string, representing the per base coverages, offset 33.

=cut

sub consensus{
	my $self = shift;
        my %p = (
                 hcrs => [],
                 ignore_coords => undef,
                 qual_weighted => 0,
                 use_ref_qual => 0,
                 @_
                );

	$self->_init_state_matrix(
                                  ignore_coords => $p{ignore_coords},
                                  qual_weighted => $p{qual_weighted},
                                  use_ref_qual => $p{use_ref_qual},
                                 );

	$self->_add_pre_calc_fq(@{$p{hcrs}}) if $p{hcrs};
	$self->_consensus;
	return $self->{con};
}


=head2 coverage

Calculate and return LIST of accurate per base coverages.
By default the state matrix is calculated, wether or not it has
 been computed before. Set first parameter to TRUE to prevent
 unneccessary recalculation.

=cut

sub coverage{
	my ($self, $reuse_matrix) = (@_,0);
	# calculate from _state_matrix, not the fastest way but accurate and
	#  already implemented :)
	# compute state_matrix if required/wanted
        $self->_init_state_matrix() if (!$self->{_state_matrix} || !$reuse_matrix);

	my @covs;
	foreach my $col(@{$self->{_state_matrix}}){
		if($col){
			my $s = 0;
			$s += $_ for (grep{$_}@$col);
			push @covs, $s;
		}else{
			push @covs, 0;
		}
	};
	return @covs;
}

=head2 chimera

By default the state matrix is calculated, wether or not it has
 been computed before. Set first parameter to TRUE to prevent
 unneccessary recalculation.

=cut

sub chimera{
	my ($self, $reuse_matrix) = (@_,0);
	# compute state_matrix if required/wanted
        $self->_init_state_matrix() if (!$self->{_state_matrix} || !$reuse_matrix);

	my @bin_bases = @{$self->{_bin_bases}};
	return unless @bin_bases > 20; # need at least 20 bins to make sense

	# low coverage bin: bin_bases[bin] << bin_max_bases
	# threshold 10% (educated guess)
	my $bin_threshold = $self->{bin_max_bases}/5 +1;

	my $lcov_bin_count = 0;
	my @lcov_bin_idxs = ();
	# find lcov_bins
	# skip terminal bins
	for (my $i=5; $i<@bin_bases-5; $i++){
		if($bin_bases[$i] <= $bin_threshold){
			$lcov_bin_count++
		}elsif($lcov_bin_count){
			# require at least 2 consecutive low cov bins to trigger chimera check
			# yet only consider local coverage drops (<5 bins)
			push @lcov_bin_idxs, [($i-$lcov_bin_count,$i-1)] if ($lcov_bin_count >= 1 && $lcov_bin_count < 5);
			$lcov_bin_count = 0; # reset
		}
	}

	my @coords;
	# get the state matrix columns
	foreach my $lcov_bin_idxs (@lcov_bin_idxs){
		my $mat_from = ($lcov_bin_idxs->[0]-1) * $self->{bin_size};
		my $mat_to = ($lcov_bin_idxs->[1]+2) * $self->{bin_size} -1;

		# uncovered columns in lcov are zero quality anyways
		next if grep{!@$_}@{$self->{_state_matrix}}[$mat_from .. $mat_to];

		#heterozygosity proof -> left hand and right hand matrix separately
		my $fl = $lcov_bin_idxs->[0]-4;
		my $tr = $lcov_bin_idxs->[1]+5;
		my $delta = int(($tr - $fl - 1)/2);
		my $tl = $fl + $delta;
		my $fr = $tr - $delta;

		my @alns_l_by_bin = $self->alns_by_bins($fl, $tl);

		my @alns_l;
		foreach(@alns_l_by_bin){
			push @alns_l, @$_;
		};
		my @alns_r_by_bin = $self->alns_by_bins($fr, $tr);
		my @alns_r;
		foreach(@alns_r_by_bin){
			push @alns_r, @$_;
		};

		my ($mat_l) = $self->State_matrix(
			alns => \@alns_l,
			#states => $self->{_states},
			#'length' => $self->len,
		);

		my ($mat_r) = $self->State_matrix(
			alns => \@alns_r,
			#states => $self->{_states},
			#'length' => $self->len,
		);

		my @mat_r = @{$mat_r}[$mat_from .. $mat_to];
		my @mat_l = @{$mat_l}[$mat_from .. $mat_to];

		my @hx_delta;
		for(my $i=0; $i< @mat_r; $i++){
			# only consider columns, where both sides contribute
			next unless @{$mat_r[$i]} && @{$mat_l[$i]};

			my $hx_r = Hx($mat_r[$i]);
			my $hx_l = Hx($mat_l[$i]);
			my $hx_gt = $hx_r > $hx_l ? $hx_r : $hx_l;

			# combined column
			my @col;
			my ($max) = sort{$b <=>$a}(scalar @{$mat_l[$i]}, scalar @{$mat_r[$i]});

			for(my $j=0; $j<$max; $j++){
				if($mat_l[$i][$j] && $mat_r[$i][$j]){
					push @col, $mat_l[$i][$j] + $mat_r[$i][$j]
				}elsif($mat_l[$i][$j]){
					push @col, $mat_l[$i][$j];
				}else{
					push @col, $mat_r[$i][$j];
				}
			}

			# delta of combined entropy and greater entropy of both single ones
			push @hx_delta, Hx(\@col) - $hx_gt;

		}

                # Hx > 0.7:
                #  4:1 = 0.72
                #  5:1 = 0.65
                #  8:2 = 0.72
                # 12:3 = 0.72

		push @coords, {
			col_range => [$mat_from + $self->{bin_size}, $mat_to - $self->{bin_size}],
			hx => \@hx_delta,
			score => (grep{$_> 0.7}@hx_delta) / @hx_delta,  # number of + columns normalized to total n columns
		}if @hx_delta;

		#TODO: self chimera
	}

	return @coords;

}


=head2 filter_by_score/filter_by_nscore/filter_by_ncscore

Filter alignments by score/nscore/ncscore. See _ncscore() for details on ncscore
computation.

=cut

sub filter_by_score{
    my $self = shift;
    foreach ($self->aln_iids) {
        my $aln = $self->aln_by_iid($_);
        my $score = $aln->score;
        $self->remove_aln_by_iid($_) if ! defined($score) || $score < $MinScore;
    }
}

sub filter_by_nscore{
    my $self = shift;
    foreach ($self->aln_iids) {
        my $aln = $self->aln_by_iid($_);
        my $nscore = $aln->nscore;
        $self->remove_aln_by_iid($_) if ! defined($nscore) || $nscore < $MinNScore;
    }
}

sub filter_by_ncscore{
    my $self = shift;
    foreach ($self->aln_iids) {
        my $aln = $self->aln_by_iid($_);
        my $ncscore = $aln->ncscore;
        $self->remove_aln_by_iid($_) if ! defined($ncscore) || $ncscore < $MinNCScore;

    }
}

=head2 filter_rep_region_alns

Filter alignments that mostly align to repetitive regions - regions with high
amount of stacked short local alignments.

  convert:

  --        -            --
  --        -            ---
  --        -             ---
  -----  ------------- ----------
  ________________________________

  to:

  -----  ------------- ----------
  ________________________________

=cut

sub filter_rep_region_alns{
    my ($self, $reuse_matrix) = (@_,0);
    $self->_init_state_matrix() if (!$self->{_state_matrix} || !$reuse_matrix);

    # get repetitive regions
    my @cov = $self->coverage(1);
    my $cmax = $RepCoverage;
    my $high = 0;
    my @rwin;

    for (my $i=0; $i<@cov; $i++) {
        if ($cov[$i] < $cmax) {
            if ($high) {
                $high = 0;
                $rwin[$#rwin][1] = $i - $rwin[$#rwin][0];
            }
        } else {
            unless ($high) {
                $high = 1;
                push @rwin, [$i];
            }
        }
    }
    $rwin[$#rwin][1] = @cov - $rwin[$#rwin][0] if $high;


    # filter rep alns
    if (@rwin) {
        # extend repetitive regions by 150bp at each side
        @rwin = map{
            $_->[0]-=150;
            $_->[1]+=300;
            $_;
        }@rwin;

        # check sequence boundaries
        if ( $rwin[0][0] < 0 ){
            $rwin[0][1]+= $rwin[0][0];
            $rwin[0][0] = 0;
        };
        if ( my $too_long = @cov - ($rwin[$#rwin][0] + $rwin[$#rwin][1]) < 0 ){
            $rwin[$#rwin][1]-= $too_long;
        }

        # filter alns
        foreach my $id ($self->aln_iids) {
            my $aln = $self->aln_by_iid($id);
            $self->remove_aln_by_iid($id) if _is_in_range([$aln->pos, $aln->length], \@rwin);
        }
    }
}

sub filter_contained_alns{
    my ($self) = @_;
    my $alns = $self->{_alns};

    # sort idx by coords-length, descending
    my @iids = keys %$alns;
    my @coords = map{[ $alns->{$_}->pos, $alns->{$_}->length ]} @iids;
    my @scores = map{$alns->{$_}->score} @iids;

    my @idx = sort{$coords[$b][1] <=> $coords[$a][1]}(0..$#iids);

    @iids = @iids[@idx];
    @coords = @coords[@idx];
    @scores = @scores[@idx];

    # filter alns
    while (@iids > 1) {
        my $iid = pop @iids;
        my $coo = pop @coords;
        if ($coo->[1] < 21) { # handle very short hits
            $coo->[0] += int($coo->[1]/2);
            $coo->[1] = 1;
        }elsif ($coo->[1] < 21) { # adjust short hits (+-10bp)
            $coo->[0]+=10;
            $coo->[1]-=20;
        }else { # adjust long hits (+- 10%);
            my $ad = int($coo->[1] * 0.1);
            $coo->[0]+=$ad;
            $coo->[1]-=(2*$ad);
        }

        if (_is_in_range($coo, \@coords)) {
            if ($coo->[1] > $coords[$#coords][1]-40) { # hits of almost identical length
                # compare by score
                my $i = @coords;
                if ($scores[$i] > $scores[$i-1]) { # exchange popped and last in queue
                    my $iid_restore = $iid;
                    $iid = pop @iids;
                    pop @coords;
                    push @iids, $iid_restore;
                    push @coords, $coo;
                }
            }
            $self->remove_aln_by_iid($iid);
        }
    }
}


=head2 filter_by_coverage

By default, add_aln_by_score will fill bins up to max_bin_bases, which is
directly controlled by $MaxCoverage. However, on data sets with uneven coverage,
e.g. meta-genomes or transcriptomic data, it might be necessary to adjust the
coverage cutoff of a particular read and filter surplus alignments.

=cut

sub filter_by_coverage{
    my ($self, $cov) = (@_);
    die "coverage required" unless $cov;

    if ($cov < $self->{max_coverage}) {
        $self->{max_coverage} = $cov;
        $self->{max_bin_bases} = $self->{max_coverage} * $self->BinSize;

        for (my $i=0; $i<@{$self->{_bin_alns}}; $i++) {
            my $iids = $self->{_bin_alns}[$i];
            my $lens = $self->{_bin_lengths}[$i];
            next unless @$iids;

            my $alnc = @$iids;
            my $rmc = 0;
            while (1) {
                last if @$lens <2;
                my $tl = 0;
                $tl += $_ for @$lens;
                last if $tl <= $self->{max_bin_bases};
                $self->remove_aln_by_iid($iids->[-1]);
                $rmc++;
            }
        }
    }
}


=head2 penalize_variants

=cut

sub penalize_variants{
    my $self = shift;
    die "Variants not yet called!\n" unless ref $self->{vars} eq 'ARRAY' && @{$self->{vars}};
    my $var_c_tot=0;
    my @ref_seq = split(//, $self->ref->seq);

    my $other;
    foreach my $aid ($self->aln_iids()) {
        $other++;
        my $aln = $self->aln_by_iid($aid);
        my $f = $aln->pos -1;
        my $seq = $aln->seq_aligned;
        my $l = length($seq);

        my $var_c = 0;
        for ( my $i=0; $i<$l; $i++) {
            my $v = $self->{vars}[$f+$i];
            #TODO
            # only consider true variants (>1 state) and purely ATGC states
            next if scalar @$v < 2 || grep{length($_) > 1 || $_ =~ /[^ATGC]/}@$v;
            $var_c++ if substr($seq, $i, 1) ne $ref_seq[$f+$i];
        }

        if ($var_c) {
            $aln = $self->remove_aln_by_iid($aid);
            next if $other%2; # remove every second snp aln
            my $score = $aln->score;
            $score -= ($var_c * -60);
            $aln->score($score);
            $self->add_aln_by_score($aln);
        }
        $var_c_tot+=$var_c;
    }
    #print STDERR "penalized $var_c_tot SNVs in alignments\n";
}


=head2 haplo_coverage

This method tries to find discriminating SNPs in the state matrix, favouring the
(underrepresented) most likely haplotype. It returns the estimated coverage of
this haplotype.

=cut

sub haplo_coverage{
    my ($self) = @_;

    my @hpl_cov;

    # compute all variants
    $self->variants(min_freq => 4);
    my @ref_seq = split(//, $self->ref->seq);

    for ( my $i=0; $i<@{$self->{vars}}; $i++) {
        my $v = $self->{vars}[$i];

        # only consider true variants (>1 state) and purely ATGC states
        next if scalar @$v < 2 || grep{length($_) > 1 || $_ =~ /[^ATGC]/}@$v;

        my $r = $ref_seq[$i];
        my $fc;
        for (my $j=0;$j<@$v;$j++) { if ($v->[$j] eq $r){
            $fc = $self->{freqs}[$i][$j];
            last;
        }}
        push @hpl_cov, $fc;
    }

    return unless @hpl_cov;
    my $hpl_cov = (sort{$a<=>$b}@hpl_cov)[int($#hpl_cov *.75)]; # 75% quantile

    # significance of hpl_cov
    my $high_cov = grep{defined $_ && $_ >= $hpl_cov * 1.5 }@{$self->{covs}};
    my $df = $high_cov ? @hpl_cov / $high_cov : 0;

    #print STDERR "$self->{id}\t$hpl_cov\t$high_cov\t$self->{len}\t$df\n";
    return $df > 0.00015 ? $hpl_cov : undef;
}

##------------------------------------------------------------------------##

=head1 Accessor METHODS

=cut

=head2 id

Get/Set the file handle.

=cut

sub id{
	my ($self, $id) = @_;
	$self->{id} = $id if $id;
	return $self->{id};
}

=head2 len

Get/Set the length.

=cut

sub len{
	my ($self, $length) = @_;
	$self->{len} = $length if $length;
	return $self->{len};
}

=head2 ref

Get/Set the reference sequence Fasta::Seq or Fastq::Seq object.

=cut

sub ref{
	my ($self, $ref) = @_;
	$self->{ref} = $ref if $ref;
	return $self->{ref};
}

=head2 con

Get the consensus sequence, a Fastq::Seq object with an additional entry
 C<< $seq->{cov} >>, which contains an phred like ascii coded string,
 representing the per base coverages, offset 33. Calls the consensus method
 in case no con object is present.

=cut

sub con{
	my ($self) = @_;
	$self->consensus unless $self->{con};
	return $self->{con};
}


=head2 next_aln

Returns the next Sam::Alignment of Sam::Seq. The alns are internally stored
 in a hash and retrieved using each. Therefore the method behaves like each:
 It returns an empty list or undef, respectively, once after all alignments
 have been returned. It returns alignments in apparently random order,
 consistent as long as no alignments are added or removed.

C< get_alns > returns in identical order as long as no alignments
 are added or removed.


=cut

sub next_aln{
	my ($self) = @_;
	my $aln;
	while(1){
		if($self->{sam}){
			my $pos = (each %{$self->{_alns}})[1];
			return undef unless defined $pos;
			$aln = $self->{sam}->aln_by_pos($pos);	#

		}else{
			$aln = (each %{$self->{_alns}})[1];
		}
		return undef unless defined($aln);
		last if &{$self->{_is}}($aln);
	}
	return $aln
}


=head2 alns(<SORTED_BY_POS>)

Returns number of Sam::Alignments of the Sam::Seq in scalar context, a list
 of all Sam::Alignments in list context. Order is identical to C<next_aln> as
 long as no alignments are added or removed or sorted by position, if first
 argument to call is TRUE.

  @alns = $pb->alns();     # LIST of alignments
  $num_alns = $pb->alns();
    # faster than getting list and using it in scalar context
  @alns_sorted_by_pos = $pb->alns(1);
    # LIST of alignments, sorted by pos, slower than without sorting


=cut

sub alns{
	my ($self, $sorted_by_pos) = @_;
	wantarray || return scalar keys %{$self->{_alns}};
        # return objects from _aln
        if($sorted_by_pos){
            return sort{ $a->{'pos'} <=> $b->{'pos'} } values %{$self->{_alns}};
        }else{
            return values %{$self->{_alns}};
        }
}


# TODO: aln_iids

=head2 aln_iids(<SORTED_BY_POS>)

Returns a list of all Sam::Alignments internal IDS (iid). Default order is
 identical to C<next_aln> as long as no alignments are added or removed. If
 first argument TRUE, iids are ordered by alignment position.

=cut

sub aln_iids{
	my ($self, $sorted_by_pos) = @_;
	wantarray || return scalar keys %{$self->{_alns}};
        # return objects from _aln
        if($sorted_by_pos){
            return sort{ $self->{_alns}{$a}{pos} <=> $self->{_alns}{$b}{pos} } keys %{$self->{_alns}};
        }else{
            return keys %{$self->{_alns}};
        }
}


=head2 aln_by_iid

Return alignment by iid.

=cut

sub aln_by_iid{
    my ($self, $iid) = @_;
    die __PACKAGE__."->aln_by_iid: IID ($iid) does not exist\n" unless exists $self->{_alns}{$iid};
    return $self->{_alns}{$iid};
}

=head2 alns_by_bins

Returns list of bins, each containing list of Sam::Alignments, decendingly
 ordered by their score. Takes FROM (default 0) and TO (default last bin)
 as optional arguments.

=cut

sub alns_by_bins{
	my ($self, $from, $to) = (@_, 0);
	my $alns = $self->{_bin_alns};
	$to = @$alns-1 unless defined $to;
	my @bins;
	for(my $i=$from; $i <= $to; $i++ ){
		push @bins, [map{ $self->{_alns}{$_}; }@{$alns->[$i]}];
	};


	return @bins;
}


=head2 bin

Compute the bin of given alignment. The bin is determined based on the
 position of the center of the read (read_start + read_length/2) and the
 currently set bin size. Returns the bin (INT).

=cut

sub bin{
	my ($self,$aln) = @_;
	return int(( $aln->pos + ( $aln->length/2 )) / $self->{bin_size})
}


=head2 is

Get/Set conditions that determine which alignments are returned by the
 next methods. Takes either a reference to a list of property bitmasks
 or a code reference to a customized test which returns 1 and 0 respectively.
 For details on bitmasks see L<Sam::Alignment>.

The test routine is executed with the parameters C<$parser_obj, $aln_obj>
 and for C<next_pair()> additionally with C< $aln_obj2 >.

  # parser returning only BAD_QUALITY alns
  my $sp = Sam::Parser->new(
  	is => [Sam::Alignment->BAD_QUALITY]
  );

  # customized parser that only returns reads with a GC content > 70%.
  my $sp = Sam::Parser->new(
  	is => sub{
  	my ($self, $aln) = @_;
  	return ($aln->seq =~ tr/GC//) / length($aln->seq) > .7 ? 1 : 0;
  })



=cut

sub is{
	my ($self, $is) = @_;
	if($is){
		if(CORE::ref($is) eq 'ARRAY'){
			$self->{_is} = eval 'sub{$_[0]->is('.join(', ', @$is).')}';
		}elsif(CORE::ref($is) eq 'CODE'){
			$self->{_is} = $is;
		}else{
			die (((caller 0)[3])." neither ARRAY nor CODE reference given!\n");
		}
	}
	return $self->{_is};
}


=head2 string

Stringify Sam::Seq object to SAM string.

  print $ss->string();
  print $ss->string(header => 1,sorted => 1);

=cut

sub string{
    my $self = shift;
    my %p = (header => 0, sorted => 0);

    # overload adds (undef, '') to the string call, which crashes @_ to hash
    if (@_ && defined($_[0])) {
        %p = (%p, @_);
    }

    my $string = "";
    if ($p{header}) {
        $string.="\@SQ\tSN:".$self->id."\tLN:".$self->len."\n";
    }
    $string.= $_ for $self->alns($p{sorted});
    return $string;
}

##------------------------------------------------------------------------##

=head1 Private Methods

=head2 _init_read_bins

Initialize bin structure required for storage of mapped read scores and ids.

=cut

sub _init_read_bins{
	my $self = shift;
	my $last_bin = int($self->len / $self->{bin_size});
	$self->{_bin_scores} = [map{[]}(0..$last_bin)];
	$self->{_bin_alns} = [map{[]}(0..$last_bin)];
	$self->{_bin_lengths} = [map{[]}(0..$last_bin)];
	$self->{_bin_bases} = [(0) x ($last_bin+1)];
}

=head2 _init_state_matrix

=cut

sub _init_state_matrix{
	my $self = shift;

	my %p = (
                 append_matrix => 0,
                 ignore_coords => undef,
                 qual_weighted => 0,
                 use_ref_qual => 0,
                 @_
                );

        my ($S, $states) = $self->State_matrix(
                                               ignore_coords => $p{ignore_coords},
                                               qual_weighted => $p{qual_weighted},
                                               use_ref_qual => $p{use_ref_qual},
                                               $p{append_matrix}
                                               ? ('matrix' => $self->{_state_matrix})
                                               : (),
                                              );

	# return state matrix
	$self->{_state_matrix} = $S;
	$self->{_states} = $states;

	return $self;

}

=head2 _add_pre_calc_fq

Add partial, already corrected sequences information (FASTQ) to the
 state_matrix, to include them in consensus.

=cut

sub _add_pre_calc_fq{
	my $self = shift;
	my @S = @{$self->{_state_matrix}};
	foreach my $coords (@_){
		my ($seq) = $self->ref->substr_seq($coords);
		my @seq = split (//, $seq->seq);
		my @freqs = Sam::Seq->Phreds2freqs($seq->phreds);

		for(my $i=0; $i<length($seq->seq); $i++){
			# never add 0, if nothing more matches, a 0 count might be introduced
			next unless $freqs[$i];
			($S[$i+$coords->[0]][$self->{_states}{$seq[$i]}])+=	$freqs[$i];
		}
	}
	$self->{_state_matrix} = \@S;
	return $self;
}


=head2 _haplo_consensus

=cut

sub _haplo_consensus{
    my $self = shift;
    my %states_rev = reverse %{$self->{_states}}; # works since values are also unique
    my $seq = '';
    my @freqs;
    my $trace;

    my $vcovs = $self->{covs};
    my $vvars = $self->{vars};
    my $vfreqs = $self->{freqs};

    for (my $i=0; $i<$self->len; $i++){
        # uncovered col
        unless ($vcovs->[$i]){
            $seq.= $self->{ref} ? substr($self->{ref}{seq}, $i, 1) : 'n';
            push @freqs, 0;
            #$trace.='M';
            next;
        }
        # TODO MaxInsertSize

        my $j = 0;
        if ($self->{ref} && @{$vvars->[$i]} >1) { # SNVs
            my $r = substr($self->{ref}{seq}, $i, 1);
            for (my $k=0; $k<@{$vvars->[$i]};$k++) { # loop through SNVs
                if ($vvars->[$i][$k] eq $r){ # find first ref matching state
                    $j=$k;
                    last;
                }
            }
        }

        #print STDERR "$i\t$vvars->[$i][$j] : @{$vvars->[$i]} @{$vfreqs->[$i]}\n" if $j;
        $seq.= $vvars->[$i][$j];
        push @freqs, $vfreqs->[$i][$j];


    }

    $self->{con} = Fastq::Seq->new(
        '@'.$self->{id},
        $seq,
        '+',
        Fastq::Seq->Phreds2Char( [Sam::Seq->Freqs2phreds(@freqs)] , $self->{phred_offset} ),
        cov => Fastq::Seq->Phreds2Char([@freqs], $self->{phred_offset}),
        phred_offset => $self->{phred_offset},
        #trace => $trace,
        #cigar => Sam::Seq->Trace2cigar($trace),
    );
}


=head2 _consensus

=cut

sub _consensus{
	my $self = shift;
	my %states_rev = reverse %{$self->{_states}}; # works since values are also unique
	my $seq = '';
	my @freqs;
	my $trace;
	my $col_c = -1;

	foreach my $col (@{$self->{_state_matrix}}){
		$col_c++;
		# uncovered col
		unless (scalar @$col){
			$seq.= $self->{ref} ? substr($self->{ref}{seq}, $col_c, 1) : 'n';
			push @freqs, 0;
			$trace.='M';
			next;
		}

		my $idx=undef;
		my $max_freq=0;
		my $i;
		my $cov;

		# majority vote
		for($i=0; $i<@$col; $i++){
			# get all defined states
			if (defined(my $freq = $col->[$i])){
				# add state frequency to coverage
				$cov+=$freq;
				# check if more frequent than previous states
				if($freq > $max_freq){
					# exception: long prominent inserts (>3) are very ugly,
					# probable mapping artefacts caused by cheap gap costs
					# compared to mismatch which might lead to long gaps
					# at read ends close to error rich regions. These long
					# gaps should not be considered.
					# Insert state has to have idx > 4 (not A,T,G,C or -)
					# for check to make sense
					next if($MaxInsLength && $i > 4 && length $states_rev{$i} > $MaxInsLength);

					$max_freq = $freq;
					$idx = $i;
				};
			}
		};

		# check $max_freq, necessary due to long gap exception
		unless ($max_freq){
			$seq.= $self->{ref} ? substr($self->{ref}{seq}, $col_c, 1) : 'n';
			push @freqs, 0;
			$trace.='M';
			next;
		}

		# insertion on reference
		if ($idx == 4){
			$trace.='I';
			next
		};

		# get most prominent state
		my $con = $states_rev{$idx};
		#use Data::Dumper;
		#printf "%s\t%s %s\n", $self->{ref}->id(), $con, $idx if $con =~ /-/;
		$seq.= $con;
		push @freqs, ($max_freq) x length($con);
		$trace.= length($con) == 1 ? 'M' : 'M'.('D'x (length($con)-1));
	}

	# compress trace to cigar



	$self->{con} = Fastq::Seq->new(
		'@'.$self->{id},
		$seq,
		'+',
		Fastq::Seq->Phreds2Char( [Sam::Seq->Freqs2phreds(@freqs)] , $self->{phred_offset} ),
		cov => Fastq::Seq->Phreds2Char([@freqs], $self->{phred_offset}),
		phred_offset => $self->{phred_offset},
		trace => $trace,
		cigar => Sam::Seq->Trace2cigar($trace),
	);

	return $self;
}


=head2 variants

By default the state matrix is calculated, wether or not it has been computed
 before. Set C<reuse_matrix => 1> to prevent unneccessary recalculation.

=cut

*_variants = \&variants; # backward comp

sub variants{
    my $self = shift;
    die __PACKAGE__."->verbose: uneven number of options: @_" if @_%2;

    my %p = (
        min_prob => 0,
        min_freq => 4,
        reuse_matrix => 0,
        @_
    );

    # compute state_matrix if required/wanted
    $self->_init_state_matrix() if (!$self->{_state_matrix} || !$p{reuse_matrix});

    my @seq;
    my %states_rev = reverse %{$self->{_states}}; # works since values are also unique
    $self->{covs} = [];
    $self->{vars} = [];
    $self->{freqs} = [];
    $self->{probs} = [];

    foreach my $col (@{$self->{_state_matrix}}) {
        # cov
        unless($col){
            push @{$self->{covs}}, 0;
            push @{$self->{vars}}, ['?'];
            push @{$self->{freqs}},[''];
            push @{$self->{probs}},[''];
            next;
        }

        my $cov;
        my %vars;
        # variants
        for (my $i=0; $i<@$col; $i++) {
            if (defined(my $v = $col->[$i])) {
                $cov+= $v;
                $vars{$states_rev{$i}} = $v;
            }
        }

        push @{$self->{covs}}, $cov;
        my @vars = sort{$vars{$b} <=> $vars{$a}}keys %vars;
        my @freqs = @vars{@vars};
        my @probs = map{$_/$cov}@freqs;
        my $k = @vars;
        if ($p{min_freq} ){
            $k = grep{$_>= $p{min_freq}}@freqs;
        }

        if($p{min_prob}){
            my $kp = grep{$_>= $p{min_prob}}@probs;
            $k = $kp if $kp < $k;
        }

        $k--;
        push @{$self->{vars}}, $k >= 0 ? [@vars[0..$k]] : [$vars[0]];
        push @{$self->{freqs}}, $k >= 0 ? [@freqs[0..$k]] : [$freqs[0]];
        push @{$self->{probs}}, $k >= 0 ? [@probs[0..$k]] : [$probs[0]];
    }
    # rel freq

    return $self;
}


sub _phred_Hx{
	my ($self, $col) = @_;
	my $total = 0;
	my @states = grep{$_}@$col;
	$total += $_ for @states;
	my @Px = map{$_/$total}@states;
	my $Hx;
	$Hx -= $_ for map{$_ * log $_}@Px;
	$Hx = 1 if $Hx > 1;
	$Hx = 0 if $Hx < 0;
	# prevent phreds >40 if actual coverage > $MaxCoverage
	# TODO: just a sloppy fix...
	# TODO !!!
	$total = $self->{max_coverage} if $total > $self->{max_coverage};
	return chr(
		int(
		((1-$Hx)/($self->{max_coverage}/$total))
		*30)							# phred range - phred minimum
		+ 10							# phred minimum
		+ $self->{phred_offset});  		# phred offset
}



=head2 _stddev

Takes a reference to a list of values and returns mean and stddev of the
 list. Additionally takes the mean of the list as second parameter, in it
 has been calculated before to speed up computation.

=cut

sub _stddev{
	my($values, $mean1) = (@_);
	#Prevent division by 0 error in case you get junk data
	return undef unless scalar @$values;

	# calculate mean unless given
	unless(defined $mean1){
		# Step 1, find the mean of the numbers
		my $total1 = 0;
		$total1 += $_  for @$values;
		my $mean1 = $total1 / (scalar @$values);
	}


	# find the mean of the squares of the differences
	# between each number and the mean
	my $total2 = 0;
	$total2 += ($mean1-$_)**2 for @$values;
	my $mean2 = $total2 / (scalar @$values);

	# standard deviation is the square root of the
	# above mean
	return sqrt($mean2);
}

=head2 _is_in_range

Test if a value/range lies within a given set of ranges. ranges are expected in
[OFFSET, LENGTH] format.

  _is_in_range(5, [[0, 3], [4,7]])
  _is_in_range([2,2], [[0, 3], [4,7]])

=cut

sub _is_in_range{
    my ($c, $ranges) = @_;
    die __PACKAGE__."::_is_in_range: requires exactly to arguments: VALUE or RANGE[OFFSET, LENGTH],  RANGES[[OFFSET, LENGTH][OFFSET, LENGTH]]" unless @_ == 2;

    if (CORE::ref $c eq "ARRAY") {
        my $c1 = $c->[0];
        my $c2 = $c->[0] + $c->[1]-1;
        for my $r (@$ranges){
            if (
                ($c1 >= $r->[0] && $c1 < $r->[0] + $r->[1]) &&
                ($c2 >= $r->[0] && $c2 < $r->[0] + $r->[1])
            ){
                return 1;
            }
        }
    }elsif (! CORE::ref $c) {
        for my $r (@$ranges){
            return 1 if $c >= $r->[0] && $c < $r->[0] + $r->[1];
        }
    }else {
        die __PACKAGE__."::_is_in_range: first arguments needs to be SCALAR or ARRAY ref";
    }
    return 0;
}

##------------------------------------------------------------------------##

=head1 AUTHOR

Thomas Hackl S<thomas.hackl@uni-wuerzburg.de>

=cut



1;
