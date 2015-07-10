package Sam::Parser;

use warnings;
use strict;

use Sam::Alignment qw(:flags);

our $VERSION = '1.3.0';

=head1 NAME

Sam::Parser.pm

=head1 DESCRIPTION

Parser module for SAM format files.

=cut

=head1 SYNOPSIS

  use Sam::Parser;
  use Sam::Alignment ':flags';

  # SAM from STDIN
  my $sp = Sam::Parser->new();

  # BAM, profits from .bai (faster and more functionality)
  my $sp = Sam::Parser->new(file => "/path/to/file.bam");

  # Fancy
  open(my $fh, 'samtools view file.bam | filter-bam |');
  my $sp = Sam::Parser->new(fh => $fh);

  # parser for file handle and with customized is routine
  # read starts with 'C'
  my $sp = Sam::Parser->new(
    fh => \*SAM,
    is => sub{ substr($_[0]->seq, 0, 1) eq 'C' }
  );

  # header
  print $sp->header;

  # print read ids of all reads with bad quality
  while( my $aln = $sp->next_aln() ){
    print $aln->qname() if $aln->is_bad_quality();
  }

  # seek the begin of the alignments for reparsing
  $sp->seek_alignment_section();

  # reset the 'is' routine
  $sp->is(MAPPED_BOTH);

  # print sequences of read pairs with both reads mapped
  while( my ($aln1, $aln2) = $sp->next_pair() ){
    print $aln1->seq().", ".$aln2->seq()."\n";
  }

=cut

=head1 Constructor METHOD

=head2 new

Initialize a sam parser object. Takes parameters in key => value format.

  fh => \*STDIN,
  file => undef,
  is => undef,
  mode => '<',   # read,
                 # '+>': read+write (clobber file first)
                 # '+<': read+write (append)
                 # '>' : write (clobber file first)
                 # '>>': write (append)
=back

=cut

my @ATTR_SCALAR = qw(file fh is mode
                     samtools_path samtools region
                     _header_fh _is_bai _idxstats_fh _idxstats
                );

my %SELF;
@SELF{@ATTR_SCALAR} = (undef) x scalar @ATTR_SCALAR;

sub new{
	my $class = shift;

	my $self = {
            %SELF,
            # defaults
            mode => '<',
            samtools => 'samtools',
            # overwrite defaults
            @_,
            # protected
            _line_buffer => undef,
            _aln_section => undef,
            _is => undef,
            _is_bai => undef,
            _idxstats_fh => undef,
	};

	bless $self, $class;

	# open file in read/write mode
        die "Either file or fh required\n" if ($self->file && $self->fh);
        $self->fh(\*STDIN) if (!$self->file && !$self->fh);

        $self->samtools(join("/", grep{$_}($self->samtools_path, $self->samtools)));

        $self->file2fh if $self->file;
        $self->cache_header unless $self->_is_bai;

	# prepare _is test routine
	$self->is($self->{is}) if $self->{is};

	return $self;
}

sub DESTROY{
    my $self = shift;
    foreach my $fh (qw(fh _header_fh _idxstats_fh)) {
        close $self->$fh if $self->$fh;
    }
}








############################################################################


=head1 Object METHODS

=cut

=head2 fh, mode, samtools, samtools_path

Get/set ...

=cut

sub _init_accessors{
    no strict 'refs';

    # generate simple accessors closure style
    foreach my $attr ( @ATTR_SCALAR ) {
        next if $_[0]->can($attr); # don't overwrite explicitly created subs
        *{__PACKAGE__ . "::$attr"} = sub {
            $_[0]->{$attr} = $_[1] if @_ == 2;
            return $_[0]->{$attr};
        }
    }
}


=head2 next_aln

Loop through sam file and return next 'Sam::Alignment' object (meeting the
 'is' criteria if specified).

=cut

sub next_aln{
	my ($self) = @_;
	my $sam = readline($self->{fh});
        return unless defined $sam; # eof

        my $aln = Sam::Alignment->new($sam);
        return $aln if !$self->{_is} or &{$self->{_is}}($aln);
	return;
}

=head2 next_seq

Return next Sam::Seq object. Only works with BAM.

=cut

sub next_seq{
    my ($self) = @_;
    die (((caller 0)[3]).": only works on indexed BAM files\n") unless $self->_is_bai;

    my ($id, $len);
    if ($self->region){
        ($id = $self->region) =~ s/(:([0-9,]+)|:([0-9,]+-[0-9,]+))$//;
        $len = ($self->idxstat($id))[1];
    }else {
        ($id, $len) = $self->next_idxstat;
        ($id, $len) = $self->next_idxstat if $id eq "*"; # dont look at unmapped
    }

    return unless defined($id);
    my $ss = Sam::Seq->new(id => $id, len => $len);

    my $sp = Sam::Parser->new(file => $self->file, region => $id);
    while ( my $aln = $sp->next_aln ) {
        $ss->add_aln_by_score($aln);
    }

    return $ss;
}

=head2 next_idxstat

Read idxstats. Returns LIST (id, len, #reads mapped, #reads unmapped).

=cut

sub next_idxstat{
    my ($self) = @_;
    die (((caller 0)[3]).": only works on indexed BAM files\n") unless $self->_is_bai;
    my $s = readline($self->_idxstats_fh);
    return unless defined $s;
    chomp($s);
    return split("\t", $s)
}


=head2 aln_by_pos

Get an alignment, based on a byte offset position, return 'Sam::Alignment'
 object (if meeting the 'is' criteria if specified) or undef. Also resets
 the filehandle for C<< next_... >> methods.

  $aln = $sp->aln_by_pos($pos);

=cut

sub aln_by_pos{
	my ($self, $seek) = @_;
	my $fh = $self->{fh};

	seek($fh, $seek, 0) || return undef;

	my $l = scalar <$fh>;
	return if $l =~ m/^@/; # header section

	# return sam aln object
	my $aln = Sam::Alignment->new($l);
	return $aln if !$self->{_is} || &{$self->{_is}}($aln);
}

=head2 next_header_line

Parse linewise through sam file header information. The method returns the
 entire line if used in SCALAR context, in LIST context a tag => value list
 corresponding to the TAG:VALUE format of sam header lines, for convenience
 each @CO comment line is returned as CO => <comment> and can be written
 directly to a hash just like the other lines. In LIST context, the C<raw>
 key contains the entire line.

To retrieve a specific header line, provide the corresponding header
 subsection key. If no key is given, any next header line is returned.
 Returns FALSE if no (more) matching header lines are in the file.

The following
 subsection keys can but don't need to be present in a standard sam file.

  #key   meaning (mandatory TAGs of subsection entry)
  @HD => header line (VN)
  @SQ => reference sequence directory (SN, LN)
  @RG => read groups (ID)
  @PG => program (ID)
  @CO => comments

  # print names of all reference sequences from header
  while(%h = $sp->next_header('@SQ')){
  	print $h{'SN'}."\n";
  }


=cut

{ # backward comp
    no warnings 'once';
    *next_header_line = \&next_header;
}

sub next_header{
    my ($self, $search_tag) = (@_, '@');
    my $fh = $self->{_header_fh};

    # loop if line buffer was empty or did not return
    # get next header line
    while ( my $sam = <$fh> ) {
        if (my ($tag, $content) = $sam =~ /^(\@?(?:$search_tag)\w{0,2})\s(.*)/) {
            if (wantarray) {
                return $tag eq '@CO'
                    ? (CO => $content, raw => $sam, tag => $tag)
                    : ($content =~ /(\w\w):(\S+)/g, raw => $sam, tag => $tag);
            } else {
                return $sam;
            }
        }
        next;
    }
    return;
}

=head2 seek_alignment_section

Reset the file handle to the start of the alignment section. Resets
 C<next_aln(), next_pair()> to the first aln in the sam file.

NOTE: this operation does only work on real files, not on STDIN.

=cut

sub seek_alignment_section{
	my ($self) = @_;
	die "".((caller 0)[3]).": Filehandle is a pipe, operation requires real file!" unless -f $self->fh;
	$self->{_line_buffer} = undef;
	seek($self->fh, $self->{_aln_section} ? $self->{_aln_section} : 0 ,0);
	return $self;
}

=head2 seek_header_section

Reset the file handle to the start of the header section (beginning of the
 file). Allows you to reread the header information.

NOTE: this operation does only work on real files, not on STDIN.

=cut


sub seek_header_section{
	my ($self) = @_;
	die "".((caller 0)[3]).": Filehandle is a pipe, operation requires real file!" unless -f $self->fh;
	$self->{_line_buffer} = undef;
	seek($self->fh, 0,0);
	return $self;
}

=head2 seek

Set the filehandle to the specified byte offset. Takes two
optional arguments "POSITION" (0), "WHENCE" (0), see perl "seek" for more.
Returns 'true' on success.

NOTE: this operation does only work on real files, not on STDIN.

=cut

sub seek{
	my ($self, $offset, $whence) = (@_, 0, 0);
	return seek($self->fh, $offset, $whence);
}


=head2 append_aln

Append an alignment to the file, provided as object or string. Returns the
 byte offset position in the file.

NOTE: In case a string is provided, make sure it contains trailing newline
 since no further test is performed.

=cut

sub append_aln{
	my ($self, $aln) = @_;
	my $pos = tell($self->{fh});
	print {$self->{fh}} ref $aln ? $aln->raw : $aln;
	return $pos;
}


=head2 tell

Return the byte offset of the current append filehandle position

=cut

sub tell{
	my ($self) = @_;
	return tell($self->{fh});
}

=head2 append_tell

DEPRECATED: use C<< $fp->tell() >>

Return the byte offset of the current append filehandle position

=cut

sub append_tell{
	shift->tell(@_)
}

=head2 header

=cut

sub header{
    my ($self) = @_;
    return $self->{_header} if defined($self->{_header}); # cached

    my $cmd = $self->samtools." view -H ".$self->file;
    return qx($cmd);
}

############################################################################

=head1 Accessor METHODS

=head2 file

Get/set SAM/BAM file.

=cut

sub file{
    my ($self, $file) = @_;
    if (defined $file) {
        $self->{file} = $file;
        $self->file2fh(); # update filehandle
    }
    return $self->{file};
}

=head2 file2fh

=cut

sub file2fh{
    my ($self) = @_;

    die (((caller 0)[3]).": only read '<' and write '>' mode supported\n") unless $self->mode eq '<' or $self->mode eq '>';
    # write
    if ( $self->mode eq '>' ) {
        my $cmd = $self->samtools." view ".$self->file;
        open($self->{fh}, "|-", $cmd) or die (((caller 0)[3]).": $!\n$cmd\n");
        return $self->{fh};
    }

    -e $self->file.'.bai' ? $self->_is_bai(1) : $self->_is_bai(0);

    # read
    die $self->file." doesn't exist\n" unless -e $self->file;

    if ( $self->_is_bai ){
        my $cmd = $self->samtools." view ".$self->file.($self->region ? " ".$self->region : "");
        open($self->{fh}, "-|", $cmd) or die (((caller 0)[3]).": $!\n$cmd\n");

        my $idxstats_cmd = $self->samtools." idxstats ".$self->file;
        open($self->{_idxstats_fh}, "-|", $idxstats_cmd) or die (((caller 0)[3]).": $!\n$idxstats_cmd\n");

        my $header_cmd = $self->samtools." view -H ".$self->file;
        open($self->{_header_fh}, "-|", $header_cmd) or die (((caller 0)[3]).": $!\n$header_cmd\n");
    }else{ # only one chance to cache header
        my $cmd = $self->samtools." view -h ".$self->file.($self->region ? " ".$self->region : "");
        open($self->{fh}, "-|", $cmd) or die (((caller 0)[3]).": $!\n$cmd\n");
    }

    return $self->fh;
}

=head2 cache_header

Read and cache header from stream. Automatically done if no .bai.

=cut

sub cache_header{
    my ($self) = @_;
    my $fh = $self->fh;
    my $h = '';
    my $c;
    while (1) {
        $c = $fh->getc();
        last if ! defined($c) or $c ne '@';
        $h.= $c.<$fh>;
    }
    $fh->ungetc(ord($c)) if defined $c;       # push back peek

    $self->{_header} = $h;
    open($self->{_header_fh}, '<', \$self->{_header}) or die $!;
}

=head2 region

Get/set region. Only works with BAM files.

=cut

sub region{
    my ($self, $region) = @_;
    if (@_>1) {
        $self->{region} = $region;
        $self->file2fh; # update fh
    }
    return $self->{region};
}

=head2 idxstat

Get idxstat for specific ID. Parses and caches idxstats.

=cut

sub idxstat{
    my ($self, $id) = @_;

    unless ( $self->_idxstats ) {
        my $idxstats_cmd = $self->samtools." idxstats ".$self->file." |";
        open(IDX, $idxstats_cmd) or die (((caller 0)[3]).": $!\n$idxstats_cmd\n");
        $self->{_idxstats} = {};
        while (defined(<IDX>)) {
            chomp;
            my (@s) = split("\t", $_);
            $self->{_idxstats}{$s[0]} = \@s;
        }
        close IDX;
    }

    die "ID \"$id\" doesn't exists in bam index" unless exists $self->_idxstats->{$id};
    return @{$self->_idxstats->{$id}};
}

=head2 is

Get/Set conditions that determine which alignments are returned by the
 next methods. Takes either a reference to a list of property bitmasks
 or a code reference to a customized test which returns 1 and 0 respectively.
 To explicitly deactivate testing, provide a value that evaluates to FALSE.
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

  # deactivate testing
  my $sp->is(0);

=cut

sub is{
	my ($self, $is) = @_;
	if(@_== 2){
		unless($is){
			$self->{_is} = undef;
		}elsif(ref($is) eq 'ARRAY'){
			$self->{_is} = eval 'sub{$_[0]->is('.join(', ', @$is).')}';
		}elsif(ref($is) eq 'CODE'){
			$self->{_is} = $is;
		}else{
			die (((caller 0)[3])." neither ARRAY nor CODE reference given!\n");
		}
	}
	return $self->{_is};
}

# init closure accessors
__PACKAGE__->_init_accessors();

=head1 AUTHOR

Thomas Hackl S<thackl@lim4.de>

=cut



1;
