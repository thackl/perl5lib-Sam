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
    FRAG_MIN_OVL => 2,
    FRAG_MIN_LEN => 3,
    FRAG_MIN_MAX_OVL_FRAC => .75
};

print "SRAND: ", srand(4263044362),"\n";
# SEEDs with errors:
# FIXED:
# 4263044362 -- switch with 0.5 err/frag
#   fixed by removing non-overlapping pre frags from x_pot
# 2586335731 -- req. refinement - successfully refined
#   11111011111111001111111111111111111111111111111111---1111111
#         *       * *
#   11111001111111100111111111111111111111111111111111---1111111
# 3602735325 -- req. refinement - success
# 3272643173 -- switch - prefix scoring bug
# 706007940  -- switch - dunno

my $v1 = "AT"x60;
my $cov = 5;
my $f_len = 8;
my $n = length($v1)*$cov/$f_len;

my $frags_all = sim_frags(ref => $v1, length => $f_len, n => $n);

substr($v1, 5, 1, "A");
substr($v1, 6, 1, "T");
substr($v1, 15, 1, "A");
substr($v1, 16, 1, "T");
substr($v1, 35, 1, "A");
substr($v1, 45, 1, "A");

# run through frags block-wise
# filter by length, ...
# push @block_frags if overlapping
my $frags = [];
my $max_end;
my $max_len;
my $max_len_i;
my $i;
my @haplotypes;
my $support_matrix = {0 => [], 1 => []};

# frag blocks;
foreach my $frag (@$frags_all) {
    next unless length($frag->{seq}) > FRAG_MIN_LEN;

    if ( @$frags && $max_end - $frag->{pos} > FRAG_MIN_OVL-2 ) { # >3bp overlap
        $i++;
        $max_end = $frag->{end} if $frag->{end} > $max_end;
        if ( length($frag->{seq}) > $max_len ){
            $max_len = length($frag->{seq});
            $max_len_i = $i;
        }
        push @$frags, $frag;
        next;
    }

    if (@$frags) {
        push @haplotypes, compute_local_haplotype(frags => $frags, seed_idx => $max_len_i);
    }

    $i = 0;
    $frags = [$frag];
    $max_end = $frag->{end};
    $max_len = length($frag->{seq});
    $max_len_i = $i;
}

# in case we exited for loop on short frag
if (@$frags) {
    push @haplotypes, compute_local_haplotype(frags => $frags, seed_idx => $max_len_i);
}

#print_frags(frags => \@haplotypes, use_hpl => 1, label => "haplotypes");

# merge into global haplotype
my $global_haplotype = '-' x length($v1);
foreach ( @haplotypes ) {
    substr($global_haplotype, $_->{pos}, length($_->{hpl}), $_->{hpl});
}

print $global_haplotype,"\n";

# refine by majority vote of supporters
my $refine = '';
#my $_0 = '';
#my $_1 = '';
for ( my $i=0; $i<length($v1); $i++ ) {
    my $s0 = $support_matrix->{0}[$i] // 0;
    my $s1 = $support_matrix->{1}[$i] // 0;
    my $h = substr($global_haplotype, $i, 1);

    if ($h) {
        if ( $s1 < $s0 ) {
            $refine.= '*';
            substr($global_haplotype, $i, 1, "0");
        }else {
            $refine.= ' ';
        }
    } else {
        if ( $s0 < $s1 ) {
            $refine .= '*';
            substr($global_haplotype, $i, 1, "1");
        } else {
            $refine .= ' ';
        }
    }
}

print $refine,"\n";
print $global_haplotype,"\n";

#print Dumper($support_matrix);


##----------------------------------------------------------------------------##

sub compute_local_haplotype{
    my %p = (@_);
    my $frags = $p{frags};

    #print "\ncompute_local_haplotype(seed_idx => $p{seed_idx})\n";

    add_frag_haplotypes(frags => $frags, ref => $v1);
    print_frags(frags => $frags, max_len => length($v1), use_hpl => 1, label => "[local frags]");

    # seed
    my $hpl = {%{$frags->[$p{seed_idx}]}};
    $hpl->{hplc} = $hpl->{hpl};
    $hpl->{hplc} =~ tr/10/01/;
    $hpl->{1} = [];
    $hpl->{0} = [];
    my $f_max_len = length($hpl->{seq});

    # use max ref similar (most "1") complement as h[0];
    if ($hpl->{hpl} =~ tr/1// < $hpl->{hplc} =~ tr/1//){
        ($hpl->{hpl}, $hpl->{hplc}) = ($hpl->{hplc}, $hpl->{hpl});
    }

    if (@$frags == 1) { # short-circuit singletons
        return $hpl;
    }

    print_frags(frags => [$hpl], max_len => length($v1), use_hpl => 1, label => "[seed frag]");

    extend_suf(hpl => $hpl, frags => $frags, last_idx => $p{seed_idx});
    #print_frags(frags => [$hpl], max_len => length($v1), use_hpl => 1); #, label => "[sufx frag]");

    my $max_look_behind = length($frags->[$p{seed_idx}]{seq}) - FRAG_MIN_OVL;
    extend_pre(hpl => $hpl, frags => $frags, first_idx => $p{seed_idx}, max_look_behind => $max_look_behind);
    #print_frags(frags => [$hpl], max_len => length($v1), use_hpl => 1); #, label => "[prex frag]");

    return $hpl;
}



sub extend_pre{
    my %p = (reuse_frags => [], @_);
    my $hpl = $p{hpl};
    my $f = $p{frags};

    # look back x frags for best ovl
    my $i;
    my $max_ovl = 0;
    my $min_ovl_i = $p{first_idx};
    for ($i=$p{first_idx}-1; $i>=0; $i--) {
        # can't use ovl in look-behind
        # --0000000------ < valid but missed
        # --1111---------
        # --11111-------- < no ovl
        # -------1111111-
        # my $ovl =  $f->[$i]->{end} - $hpl->{pos} + 1; # min ovl
        # last if $ovl < FRAG_MIN_OVL;

        my $ovl = $hpl->{end} - $f->[$i]->{pos} + 1; # min ovl
        unless ( $ovl < FRAG_MIN_OVL ){ # seen a frag with valid ovl
            $min_ovl_i = $i;
        }

        last if $f->[$i]{pos} < $hpl->{pos} - $p{max_look_behind};
    }
    $i = $min_ovl_i;

    my $x_all = $p{reuse_frags};
    if ( $i < $p{first_idx}) { # new extender
        unshift @$x_all, @{$frags}[$i .. $p{first_idx}-1];
    }
    return unless @$x_all;

    # assign scores to all frags
    score_and_flip_frags(frags => $x_all, hpl => $hpl);

    # potential extenders
    my $x_pot = [];
    my $x_pre = []; # look-behind hits with too short overlap

    # remove containees, record support/coverage
    foreach ( @$x_all ) {
        if ( $_->{pos} < $hpl->{pos} ) {
            if ( $_->{end} - $hpl->{pos} + 1 >= FRAG_MIN_OVL ){
                push @$x_pot, $_;
            }else {
                push @$x_pre, $_;
            }
        } else {
            add_frag_to_support_matrix($_);
        }
    }
    #print_frags(frags => $x_con, use_hpl => 1, label => "[containee frags]");
    return unless @$x_pot;

    # best extender (highest score)
    my $x_best;
    ($x_best, @$x_pot) = @{$x_pot}[sort { $x_pot->[$b]{score} <=> $x_pot->[$a]{score} || length($x_pot->[$b]{seq}) <=> length($x_pot->[$a]{seq}) }(0..$#$x_pot)];

    push @$x_pot, @$x_pre;
    print_frags(frags => [$x_best], use_hpl => 1, use_scores => 1, label => "[pre-xtend best]", max_len => length($v1));
    print_frags(frags => $x_pot , use_hpl => 1, use_scores => 1, label => "[pre-xtend sub]", max_len => length($v1));

    add_frag_to_support_matrix($x_best);

    my $x_len = $hpl->{pos} - $x_best->{pos};  # extension length
    my $x_hpl = substr($x_best->{hpl}, 0, $x_len);
    my $x_hplc = $x_hpl;
    $x_hplc =~ tr/10/01/;

    $hpl->{pos} = $x_best->{pos}; # set new pos
    $hpl->{hpl} = $x_hpl.$hpl->{hpl};
    $hpl->{hplc} = $x_hplc.$hpl->{hplc};


    #print_frags(frags => [$hpl], use_hpl => 1, use_scores => 1, label => "[extended hpl]", max_len => length($v1));

    extend_pre(hpl => $hpl, frags => $p{frags}, first_idx => $i, reuse_frags => $x_pot, max_look_behind => $p{max_look_behind});
}

sub extend_suf{
    my %p = (reuse_frags => [], @_);
    my $hpl = $p{hpl};
    my $f = $p{frags};

    # look at the next x frags for best ovl
    my $i;
    my $max_ovl = 0;
    for ($i=$p{last_idx}+1; $i<@$f; $i++) {
        # filter prior:
        # next if length($f[$i]->{seq}) < FRAG_MIN_LEN; # ignore shorties

        my $ovl = $hpl->{end} - $f->[$i]{pos} + 1; # min ovl
        last if $ovl < FRAG_MIN_OVL;

        # probably insignificant pros causes trouble with containees
        #  111111
        #  11111      # strong containee
        #     111111  # week extendee
        #if ($ovl >= $max_ovl ){
        #    $max_ovl = $ovl
        #}else {
        #    last if $ovl < FRAG_MIN_MAX_OVL_FRAC * $max_ovl;
        #}
    }
    $i--;

    my $x_all = $p{reuse_frags};
    if ( $i > $p{last_idx}) { # new extenders
        push @$x_all, @{$frags}[$p{last_idx}+1 .. $i ];
    }
    return unless @$x_all;

    score_and_flip_frags(frags => $x_all, hpl => $hpl);

    # potential extenders
    my $x_pot = [];
    my $x_con = [];

    # remove containees
    foreach ( @$x_all ) {
        if ( $_->{end} > $hpl->{end} ){
            push @$x_pot, $_;
        }else {
            add_frag_to_support_matrix($_);
        }
    }

    #print_frags(frags => $x_con, use_hpl => 1, label => "[containee frags]");
    return unless @$x_pot;

    my $x_best;
    ($x_best, @$x_pot) = @{$x_pot}[sort { $x_pot->[$b]{score} <=> $x_pot->[$a]{score} || length($x_pot->[$b]{seq}) <=> length($x_pot->[$a]{seq}) }(0..$#$x_pot)];

    add_frag_to_support_matrix($x_best);

    #print_frags(frags => [$x_best], use_hpl => 1, use_scores => 1, label => "[suf-xtend best]", max_len => length($v1));
    #print_frags(frags => $x_pot, use_hpl => 1, use_scores => 1, label => "[suf-xtend sub]", max_len => length($v1));

    my $x_len = $x_best->{end} - $hpl->{end}; # extension length
    my $x_hpl = substr($x_best->{hpl}, -$x_len);
    my $x_hplc = $x_hpl;
    $x_hplc =~ tr/10/01/;

    $hpl->{end} = $x_best->{end}; # set new end
    $hpl->{hpl}.= $x_hpl;
    $hpl->{hplc}.= $x_hplc;

    #print_frags(frags => [$hpl], use_hpl => 1, use_scores => 1, label => "[extended hpl]", max_len => length($v1));

    extend_suf(hpl => $hpl, frags => $p{frags}, last_idx => $i, reuse_frags => $x_pot);
}


sub score_and_flip_frags{
    my %p = (@_);
    foreach (@{$p{frags}}) {
        my $shift = $p{hpl}{pos} - $_->{pos};
        my ($m1, $m2);
        if ($shift > 0){ # pre
            my $ovl = $_->{end} - $p{hpl}{pos}+1;
            #print "prefix\n",substr($_->{hpl}, -$ovl),"\n",$p{hpl}{hpl},"\n";
            my $s1 = substr($_->{hpl}, -$ovl) ^ $p{hpl}{hpl};
            $m1 = $s1 =~ tr/\0//; # count matches
            my $s2 = substr($_->{hpl}, -$ovl) ^ $p{hpl}{hplc};
            $m2 = $s2 =~ tr/\0//; # count matches
        } else { # suf or containee
            #print "suffix\n",substr($p{hpl}{hpl}, abs($shift)),"\n",$_->{hpl},"\n";
            my $s1 = substr($p{hpl}{hpl}, abs($shift)) ^ $_->{hpl};
            $m1 = $s1 =~ tr/\0//; # count matches
            my $s2 = substr($p{hpl}{hplc}, abs($shift)) ^ $_->{hpl};
            $m2 = $s2 =~ tr/\0//; # count matches
        }

        if($m2 < $m1){
            $_->{score} =  $m1;
            $_->{flipped} = 0;
        }else {
            $_->{score} = $m2;
            $_->{hpl} =~ tr/10/01/;
            $_->{flipped} = 1;
        }
    }
}

sub add_frag_to_support_matrix{
    my ($f) = (@_);
    my $k = $f->{pos};
    foreach my $state ( split(//, $f->{hpl})) {
        $support_matrix->{$state}[$k]++;
        $k++;
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
    print $p{label},"\n" if $p{label};
    foreach (@{$p{frags}}) {
        my $k = $p{use_hpl} ? "hpl" : "seq";
        my $suf = $p{max_len} ? $p{max_len} - ($_->{pos} + length($_->{$k})) : 0;
        my $score = "";
        if ( $p{use_scores} && exists $_->{score} ) {
            $score = " [$_->{score}]". ($_->{flipped} ? " !" : "");
        }
        print "-" x $_->{pos} , $_->{$k}, "-" x $suf,$score,"\n";
    }
}

sub sim_frags{
    my %p = (@_);
    my $max = length($p{ref})+$p{length}-1;
    my @p = sort{$a <=> $b} map{ int(rand($max)) - $p{length}+1}(0..$p{n});
    my @f;
    for (my $i=0; $i<@p; $i++) {
        my ($p, $l) = ($p[$i], int(rand($p{length}-2))+3 );

        if ( $p < 0 ) { # short fragments at start
            $l += $p; $p=0
        }
        if ( (my $ovh = $p + $l - length($p{ref})) > 0) { # short frags at end
            $l-=$ovh;
        }
        next if $l<1;
        my $f = {
            seq => substr($p{ref}, $p, $l),
            pos => $p
        };
        $f->{seq} =~ tr/ATGC/TACG/ if int(rand(2)); # randomize strand
        $f->{end} = $p + length($f->{seq}) -1;
        # 1 error / 50 % frags
        if (int(rand(2))) {
            my $pos = int(rand(length($f->{seq})));
            my $org = substr($f->{seq}, $pos, 1);
            my $err = (grep{$_ ne $org}qw(A T G C))[int(rand(3))];
            substr($f->{seq}, $pos, 1, $err);
        }
        push @f, $f;
    }
    return \@f;
}

sub add_frag_haplotypes{
    my %p = (@_);
    foreach (@{$p{frags}}) {
        my $dx = substr($p{ref}, $_->{pos}, length($_->{seq})) ^ $_->{seq};
        $dx =~ tr/\0/0/c;
        $dx =~ tr/\0/1/;
        $_->{hpl} = $dx;

        #print " "x $_->{pos}, substr($p{ref}, $_->{pos}, length($_->{seq})),"\n",
        #    " "x $_->{pos}, $_->{seq},"\n",
        #    " "x $_->{pos}, $_->{hpl},"\n";
    };
    return $frags;
}

sub filter_frags{
    my %p = (@_);
    return [grep{length($_->[0]) >2}@{$p{frags}}];
}
