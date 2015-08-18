#!/usr/bin/env perl
use warnings;
use strict;

use Data::Dumper;

use constant {
    NUC => 0,
    POS => 1,
    HPL => 2,
    SCORE => 3,
    HPL_IDX => 4,
    MIN_FRAG_OVL => 2,
};

my $v1 = "AT"x30;
my $v2 = $v1;
$v2 =~ tr/ATGC/TACG/; # simply complement

srand(5);

my $cov = 5;
my $f_len = 6;
my $n = length($v1)*$cov/$f_len;

my $frags = sim_frags(ref => $v1, length => $f_len, n => $n);
$frags = filter_frags(frags => $frags, min_length => 3);
$frags = assign_frag_states(frags => $frags, ref => $v1);
print $v1,"\n";
print_frags(frags => $frags, ref => $v1);

# use longest frag as seed frag
my $c=0;
while (@$frags) {
    last if $c++ > 2;
    my $frag_longest_idx = (sort{length($frags->[$b][0]) <=> length($frags->[$a][0])}(0..$#$frags))[0];
    print_frags(frags => $frags, ref => $v1);

    print @$frags." $frag_longest_idx\n";

    my ($hpl, $idx_from, $idx_to) = compute_local_haplotype( frags => $frags, seed_idx => $frag_longest_idx );
    print_haplotypes(hpl => $hpl, pos => $frags->[$frag_longest_idx][POS], ref => $v1);
    print_frags(frags => [@{$frags}[$idx_from .. $idx_to]], ref => $v1);

    my $len = $idx_to - $idx_from;
    splice(@$frags, $idx_from, $len);

}


##----------------------------------------------------------------------------##

sub compute_local_haplotype{
    my %p = (@_);
    my $idx = $p{seed_idx};
    my $frags = $p{frags};


    my @hpl = ($frags->[$idx][2], $frags->[$idx][2], [], []);
    $hpl[1] =~ tr/10/01/;
    my $pos = $frags->[$idx][1];
    my $f_max_len = $frags->[$idx][2];

    # use max ref similar (most "1") complement as h[0];
    if ($hpl[0] =~ tr/1// < $hpl[1] =~ tr/1//){
        @hpl = ($hpl[1], $hpl[0]);
    }

    # limit search space of potentially overlapping frags
    # I know length of longest frag - frags farther away cannot overlap
    # start of fragment range: [h_start - longest_frag + min_ovl, h_end - min_ovl]

    # only use frags with at least 2bp overlap
    my ($xp1, $xp2) = extender_range_pre(hpl => \@hpl, pos => $pos, idx => $idx);
    my ($xs1, $xs2) = extender_range_suf(hpl => \@hpl, pos => $pos, idx => $idx, frags => $frags);

    my @xfrags = (); # extendy frags
    push @xfrags, @{$frags}[$xp1..$xp2, $xs1..$xs2];

    print_haplotypes(hpl => \@hpl, pos => $pos, ref => $v1);
    #print_frags(frags => \@xfrags, ref => $v1);

    # end init
    
    my $c = 0;
    while (@xfrags) {
        #last if ++$c > 12;
        score_frags_vs_haplotype(frags => \@xfrags, hpl => \@hpl, pos => $pos);

        # loop to frags, look for best extender, at the same time assign containees
        my $max_score = 0;
        my $max_score_idx;
        my @keep_idx;
        for (my $i=$#xfrags; $i; $i--) {
            if (
                $xfrags[$i][POS] >= $pos &&
                    $xfrags[$i][POS]+length($xfrags[$i][HPL]) <= $pos + length ($hpl[0])
                ) {             # containee
                push @{$hpl[$xfrags[$i][POS]+2]}, $xfrags[$i];
            } else {
                if ($xfrags[$i][SCORE] > $max_score) {
                    $max_score = $xfrags[$i][SCORE];
                    push @keep_idx, $max_score_idx if defined($max_score_idx);
                    $max_score_idx = $i;
                } else {
                    push @keep_idx, $i;
                }
            }
        }

        my $is_pre = 0;
        if ( defined $max_score_idx ) {
            my $xfrag = $xfrags[$max_score_idx];
            push @{$hpl[$xfrags[$max_score_idx][POS]+2]}, $xfrag;

            # DONE: assign max score frag and extend hpl
            my $shift = $xfrag->[POS] - $pos;
            if ( $shift <0  ) { # prepend
                $is_pre++;
                my $pre = substr($xfrag->[HPL], 0, abs($shift));
                my $pre_c = $pre;
                $pre_c =~ tr/10/01/;
                $hpl[$xfrag->[HPL_IDX]].= $pre; 
                $hpl[! $xfrag->[HPL_IDX]].= $pre_c;

                $pos = $xfrag->[POS];
                $idx = $max_score_idx;
            } else {            # append
                my $shift = length($xfrag->[HPL]) - length($hpl[0]) + $shift;
                my $suf = substr($xfrag->[HPL], -$shift);
                my $suf_c = $suf;
                $suf_c =~ tr/10/01/;
                $hpl[$xfrag->[HPL_IDX]].= $suf; 
                $hpl[! $xfrag->[HPL_IDX]].= $suf_c;
            }
        }
        
        # DONE: remove containees and assigned max score frag, keep rest
        @xfrags = @xfrags[@keep_idx];
        
        # TODO: update extender_range and xfrags
        if ( $is_pre ) {
            ($xp1, $xp2) = extender_range_pre(hpl => \@hpl, pos => $pos, idx => $xp1);
            push @xfrags, @{$frags}[$xp1..$xp2];
        } else {
            ($xs1, $xs2) = extender_range_suf(hpl => \@hpl, pos => $pos, idx => $xs2, frags => $frags);
            push @xfrags, @{$frags}[$xs1..$xs2];
        }
        
        #print_haplotypes(hpl => \@hpl, pos => $pos, ref => $v1);
        #print_frags(frags => \@xfrags, ref => $v1);
    }

    return (\@hpl, $xp1, $xs2);
}



sub extender_range_suf{
    my %p = (@_); # pos, idx, hpl, frags
    my $pos_max = $p{pos} + length($p{hpl}[0]) - MIN_FRAG_OVL;
    my $idx_suf_min = $p{idx} + 1;
    $idx_suf_min = $#$frags if $idx_suf_min > $#$frags;
    my $idx_suf_max = $idx_suf_min;
    while($idx_suf_max < @{$p{frags}} && $frags->[$idx_suf_max][POS] < $pos_max){ $idx_suf_max++ };
    $idx_suf_max--;
    return ($idx_suf_min, $idx_suf_max);
}

sub extender_range_pre{
    my %p = (@_); # pos, idx, hpl
    my $pos_min = $p{pos} - length($p{hpl}[0]) + MIN_FRAG_OVL;
    my $idx_pre_max = $p{idx} - 1;
    $idx_pre_max = 0 if $idx_pre_max < 0;
    my $idx_pre_min = $idx_pre_max;
    while($idx_pre_min && $frags->[$idx_pre_min][POS] > $pos_min){ $idx_pre_min-- };
    $idx_pre_min++;
    return ($idx_pre_min, $idx_pre_max);
}


sub score_frags_vs_haplotype{
    my %p = (@_);
    foreach (@{$p{frags}}) {
        my $shift = $p{pos} - $_->[POS];
        my ($m1, $m2);
        if ($shift < 0){ # pre
            my $s1 = substr($_->[HPL], 0, abs($shift)) ^ $p{hpl}[0];
            $m1 = $s1 =~ tr/\0//; # count matches
            my $s2 = substr($_->[HPL], 0, abs($shift)) ^ $p{hpl}[1];
            $m2 = $s2 =~ tr/\0//; # count matches
        }else { # suf
            my $s1 = substr($p{hpl}[0], 0, abs($shift)) ^ $_->[HPL];
            $m1 = $s1 =~ tr/\0//; # count matches
            my $s2 = substr($p{hpl}[1], 0, abs($shift)) ^ $_->[HPL];
            $m2 = $s2 =~ tr/\0//; # count matches
        }

        if($m2 < $m1){
            $_->[SCORE] =  $m1;
            $_->[HPL_IDX] = 0;
        }else {
            $_->[SCORE] = $m2;
            $_->[HPL_IDX] = 1;
        } 
    }
}

sub print_haplotypes{
    my %p = (@_);
    my $suf = length($p{ref}) - ($p{pos} + length($p{hpl}[0]));
    print "HAPLOTYPES:\n";
    print "-" x $p{pos} , $p{hpl}[0], "-" x $suf,"\n";
    print "-" x $p{pos} , $p{hpl}[1], "-" x $suf,"\n";
};
   
sub print_frags{
    my %p = (@_);
    my $max_len = length($p{ref});
    print "FRAGS\n";
    foreach (@{$p{frags}}) {
        my ($f, $p, $h) = @$_;

        my $suf = $max_len - ($p + length($h));
        print "-" x $p , $h, "-" x $suf,"\n";
    }
}    

sub sim_frags{
    my %p = (@_);
    my $max = length($p{ref})+$p{length}-1;
    my @p = sort{$a <=> $b} map{ int(rand($max)) - $p{length} +1}(1..$p{n});
    my @f;
    for (my $i=0; $i<@p; $i++) {
        my ($p, $l) = ($p[$i], $p{length});
        if ( $p < 0 ) { # short fragments at start
            $l += $p; $p=0
        }
        if ( (my $ovh = $p + $l - length($p{ref})) > 0) { # short frags at end
            $l-=$ovh;
        }
        push @f, [substr($p{ref}, $p, $l), $p];
        $f[-1][0] =~ tr/ATGC/TACG/ if int(rand(2));
    }
    return \@f;
}

sub assign_frag_states{
    my %p = (@_);
    foreach (@{$p{frags}}) {
        my ($f, $p) = @$_;
        
        my $dx = substr($p{ref}, $p, length($f)) ^ $f;
        $dx =~ tr/\0/0/c;
        $dx =~ tr/\0/1/;
        $_->[2] = $dx;
    };
    return $frags;
}

sub filter_frags{
    my %p = (@_);
    return [grep{length($_->[0]) >2}@{$p{frags}}];
}
