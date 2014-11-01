#!/usr/bin/perl -w
# benchmark.pl --- 
# Author: Thomas Hackl <thackl@TiAkAl>
# Created: 31 Oct 2014
# Version: 0.01
package bm;
use warnings;
use strict;

use Data::Dumper;
use Benchmark qw(:all);


my $aln = {pos => 1234, seq => "ATGC"x25};
my $bin_size = 20;
bless $aln,"bm" ;

#$self->{opt} = "NM:i:2	MD:Z:2^G94	AS:i:474	XS:i:470";

cmpthese(-1, {
              A => sub{
                  my ($self) = @_;
                  return int(( $aln->pos + ( length($aln->seq)/2 )) / $bin_size)
              },
              B => sub{
                  my ($self) = @_;
                  return int(( $aln->pos_get + ( length($aln->seq_get)/2 )) / $bin_size)
              },
              C => sub{
                  my ($self) = @_;
                  return int(( $aln->{pos} + ( length($aln->{seq})/2 )) / $bin_size)
              }

});

sub bm::pos{
    my ($self, $pos, $force) = @_;
    if(defined $pos || $force){
        $self->{pos} = $pos;
    }
    return $self->{pos};
}

sub bm::seq{
    my ($self, $seq, $force) = @_;
    if(defined $seq || $force){
        $self->{seq} = $seq;
    }
    return $self->{seq};
}

sub bm::pos_get{
    my ($self) = @_;
    return $self->{pos};
}

sub bm::seq_get{
    my ($self) = @_;
    return $self->{seq};
}


=pod

cmpthese(-1, {
              A => sub{
                  while($self->{opt} =~ /(\w\w):(\w):([^\t]+)/g){
                      $self->{A}{$1} = [$2, $3];
                  };
              },
              B => sub{
                  foreach(split("\t", $self->{opt})){
                      my ($k, @v) = split(":", $_);
                      $self->{B}{$k} = \@v;
                  };
              },
              C => sub{
                  my @opt = split(/[\t:]/, $self->{opt});
                  for(my $i=0;$i<@opt;$i+=3){
                      $self->{C}{$opt[$i]} = [$opt[$i+1], $opt[$i+2]];
                  };
              },
              D => sub{
                  foreach(split("\t", $self->{opt})){
                      $self->{D}{substr($_, 0, 2)} = [substr($_, 3, 1), substr($_, 5)];
                  };
              },
              E => sub{
                  foreach(split("\t", $self->{opt})){
                      my ($k, $s1, $t, $s2, $v) = unpack("A2 A A A A*", $_);
                      $self->{E}{$k} = [$t, $v];
                  };
              },
              F => sub{
                  foreach(split("\t", $self->{opt})){
                      $self->{F}{unpack("A2", $_)} = [unpack("x3 A x A*", $_)];
                  };
              },
             });

=cut
