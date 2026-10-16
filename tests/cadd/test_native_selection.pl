#!/usr/bin/env perl
# Calls the installed official VEP methods on synthetic overlap objects.
# This checks adapter mechanics only, not annotation or scientific qualification.
use strict;
use warnings;
use Test::More;
use Scalar::Util qw(refaddr);
use CaddSelection;

{
    package Fixture::Transcript;
    sub get_all_Attributes { return []; }
    sub is_canonical { return $_[0]->{canonical}; }
    sub biotype { return 'protein_coding'; }
    sub translation { return undef; }
    sub length { return $_[0]->{length}; }
}
{
    package Fixture::Allele;
    sub isa { return $_[1] eq 'Bio::EnsEMBL::Variation::BaseTranscriptVariationAllele' || UNIVERSAL::isa(@_); }
    sub transcript { return $_[0]->{transcript}; }
    sub feature { return $_[0]->{transcript}; }
    sub base_variation_feature { return $_[0]->{vf}; }
}
{
    package Fixture::NonTranscript;
    sub base_variation_feature { return $_[0]->{vf}; }
}
{
    package Fixture::Overlap;
    sub get_all_alternate_VariationFeatureOverlapAlleles { return $_[0]->{alleles}; }
}
{
    package Fixture::Variant;
    sub get_all_VariationFeatureOverlaps { return [bless {alleles => $_[0]->{alleles}}, 'Fixture::Overlap']; }
    sub allele_string { return $_[0]->{allele_string} || 'A/G'; }
}
sub allele {
    my ($gene, $canonical, $length) = @_;
    return bless {transcript => bless {_gene_stable_id => $gene, _source_cache => 'Ensembl', canonical => $canonical,
                                      length => $length}, 'Fixture::Transcript'}, 'Fixture::Allele';
}
sub variant {
    my (@alleles) = @_;
    my $vf = bless {alleles => \@alleles}, 'Fixture::Variant';
    $_->{vf} = $vf for @alleles;
    return $vf;
}
my $order = ['canonical', 'length'];
my $plugin = CaddSelection->new({pick_order => $order}, 'cadd_transcript_sources=Ensembl+RefSeq');
my ($short, $canonical, $tie1, $tie2, $other) =
    (allele('GENEA', 0, 1000), allele('GENEA', 1, 100), allele('GENEB', 0, 50),
     allele('GENEB', 0, 50), allele('GENEC', 0, 10));
my $regulatory = bless {}, 'Fixture::NonTranscript';
my $motif = bless {}, 'Fixture::NonTranscript';
my $vf = variant($short, $tie1, $regulatory, $canonical, $tie2, $other, $motif);
my $factory = bless {pick_order => $order}, 'Bio::EnsEMBL::VEP::OutputFactory';
my $direct = $factory->pick_VariationFeatureOverlapAllele_per_gene($vf->{alleles});
my @selected;
my @annotations = map {$plugin->run($_)} @{$vf->{alleles}};
for (my $i = 0; $i < @annotations; $i++) {
    my $ordinal = $annotations[$i]->{CADD_PICK_ORDINAL};
    $selected[$ordinal - 1] = $vf->{alleles}->[$i] if $ordinal;
}
is_deeply([map {refaddr($_)} @selected], [map {refaddr($_)} @$direct], 'ordinals reproduce native per_gene order exactly');
is(scalar(@annotations), 7, 'all original PriVA rows remain');
is($annotations[0]->{CADD_PICK_ORDINAL}, 0, 'native canonical choice excludes longer noncanonical transcript');
ok($annotations[1]->{CADD_PICK_ORDINAL}, 'native stable tie keeps first transcript');
is($annotations[4]->{CADD_PICK_ORDINAL}, 0, 'second tied transcript remains unselected');
is($annotations[2]->{CADD_PICK_ORDINAL}, 1, 'native nontranscript row first');
is($annotations[6]->{CADD_PICK_ORDINAL}, 2, 'native motif row preserved before transcript groups');
is($annotations[0]->{CADD_PICK_COUNT}, 5, 'complete native row count retained');
like($annotations[0]->{CADD_PICK_CONTEXT}, qr/^[a-f0-9]{64}$/, 'source and explicit Perl context recorded');
is_deeply($plugin->{context}->{pick_order}, $order, 'official resolved pick order retained');
is($plugin->{context}->{perl_perturb_keys}, '0', 'hash perturbation explicitly recorded');
eval {CaddSelection->new({pick_order => $order, per_gene => 1})};
like($@, qr/incompatible VEP option per_gene/, 'row filtering option rejected');
my $bad = allele(undef, 0, 1);
variant($bad);
eval {$plugin->run($bad)};
like($@, qr/cached native gene/, 'missing native cached gene cannot trigger remote lookup');
my $multi = allele('GENE', 0, 1);
variant($multi)->{allele_string} = 'A/G/T';
eval {$plugin->run($multi)};
like($@, qr/biallelic/, 'multiallelic representation rejected');
my $ensembl_plugin = CaddSelection->new({pick_order => $order}, 'cadd_transcript_sources=Ensembl');
my ($ens, $ref) = (allele('ENSG', 0, 1), allele('REFG', 1, 100));
$ref->{transcript}->{_source_cache} = 'RefSeq';
variant($ens, $ref);
is($ensembl_plugin->run($ref)->{CADD_PICK_ORDINAL}, 0, 'explicit Ensembl CADD pool retains RefSeq PriVA row with zero ordinal');
is($ensembl_plugin->run($ens)->{CADD_PICK_ORDINAL}, 1, 'explicit Ensembl pool uses the same native picker');
variant($ref);
eval {$ensembl_plugin->run($ref)};
like($@, qr/intergenic representation/, 'RefSeq-only gap cannot silently lose a variant or invent intergenic annotation');
{
    local $ENV{PERL_PERTURB_KEYS} = '2';
    eval {CaddSelection->new({pick_order => $order})};
    like($@, qr/PERL_PERTURB_KEYS=0/, 'unreproducible hash environment rejected');
}
done_testing();
