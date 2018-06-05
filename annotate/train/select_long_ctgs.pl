#!env perl
use strict;
use warnings;
use Bio::SeqIO;
use Getopt::Long;

my $min_len = 100_000;
my $max_len = 1_000_000;

GetOptions( 'm|min:i' => \$min_len,
            'x|max:i' => \$max_len);

my $in = Bio::SeqIO->new(-format => 'fasta',
			 -file   => shift);

my $out = Bio::SeqIO->new(-format => 'fasta');
while ( my $s = $in->next_seq ) {
    my $len = $s->length;
    next unless ($len > $min_len && $len < $max_len);
    $out->write_seq($s);
}
