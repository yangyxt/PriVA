
use strict;
use warnings;
use Bio::DB::BigFile;
use Bio::DB::BigFile::Constants;
use JSON::PP;
my $path = shift @ARGV;
my $wig = Bio::DB::BigFile->bigWigFileOpen($path);
die "Cannot open BigWig\n" unless $wig;
my $list = $wig->chromList();
die "BigWig has no chromosome dictionary\n" unless $list;
my %chroms;
for (my $entry = $list->head; $entry; $entry = $entry->next) {
    my ($name, $size) = ($entry->name, $entry->size);
    die "Invalid BigWig chromosome\n" unless defined($name) && length($name) && $size > 0;
    die "Duplicate BigWig chromosome\n" if exists $chroms{$name};
    $chroms{$name} = 0 + $size;
    # Exercise actual signal/index reads without materializing whole tracks.
    my $window = $size < 4096 ? $size : 4096;
    $wig->bigWigSingleSummary($name, 0, $window, bbiSumMean);
    $wig->bigWigSingleSummary($name, $size - $window, $size, bbiSumMean) if $size > $window;
}
die "Empty BigWig chromosome dictionary\n" unless keys %chroms;
print JSON::PP->new->canonical->encode({chroms => \%chroms,
    reader => 'Bio::DB::BigFile', reader_module => $INC{'Bio/DB/BigFile.pm'}});
