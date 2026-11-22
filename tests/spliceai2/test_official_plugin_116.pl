#!/usr/bin/env perl
use strict;
use warnings;
use Digest::SHA ();
use File::Spec;
use Getopt::Long qw(GetOptions);
use Test::More;

# Boundary qualification, not full VEP/cache or scientific acceptance. Execute
# unchanged official release 116 parse_data/run with real Ensembl 116 objects.
# Only plugin file registration and tabix candidate retrieval are replaced.
# Default scores are invented. Optional --projected-vcf reads a bounded local
# projection; publisher score strings are never embedded in this repository.
my ($plugin_file, $vep_lib, $projected_vcf);
GetOptions('plugin-file=s' => \$plugin_file, 'vep-lib=s' => \$vep_lib,
           'projected-vcf=s' => \$projected_vcf) or die "Invalid options\n";
die "Require --plugin-file and --vep-lib\n"
  unless defined $plugin_file && defined $vep_lib && !@ARGV;
$plugin_file = File::Spec->rel2abs($plugin_file);
$vep_lib = File::Spec->rel2abs($vep_lib);
my $plugin_sha = 'd08ce7d8c2e6b3229cef638beff145e465e224f0e3fcd43fd49d501ebc2740b1';
my %api_sha = (
  'Bio/EnsEMBL/Variation/VariationFeature.pm' =>
    '1e63a719da9c8b10afd75285e1e30b73b54cf58044ec4be06c8a4228ca494502',
  'Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm' =>
    '68a7d2593faff3484d76a3e54d257a2b7045e2e665466cec40c28a7d87070780',
  'Bio/EnsEMBL/Variation/BaseVariationFeatureOverlapAllele.pm' =>
    '45492183b3f4647888680c08e8d594a6e4ee228f0eacc51145762c4444a90cda',
);

sub digest {
  open my $fh, '<:raw', $_[0] or die "Cannot open $_[0]: $!\n";
  return Digest::SHA->new(256)->addfile($fh)->hexdigest;
}
die "Official release 116 SpliceAI.pm identity mismatch\n" unless digest($plugin_file) eq $plugin_sha;
unshift @INC, $vep_lib;
require $plugin_file;
require Bio::EnsEMBL::ApiVersion;
require Bio::EnsEMBL::Transcript;
require Bio::EnsEMBL::Variation::TranscriptVariation;
require Bio::EnsEMBL::Variation::TranscriptVariationAllele;
die "Require actual Ensembl 116 API\n" unless Bio::EnsEMBL::ApiVersion::software_version() == 116;
for my $module (sort keys %api_sha) {
  die "Actual Ensembl 116 API identity mismatch: $module\n"
    unless defined $INC{$module} && digest($INC{$module}) eq $api_sha{$module};
}

{
  package FixtureSpliceAI116;
  our @ISA = ('SpliceAI');
  sub add_file { push @{$_[0]->{fixture_files}}, $_[1] }
  sub get_data {
    my ($self, @query) = @_;
    $self->{fixture_query} = \@query;
    return $self->{fixture_rows};
  }
}

sub record {
  my ($pos, $ref, $alt, $symbol, $values) = @_;
  die "Record requires eight values\n" unless @$values == 8;
  return join("\t", 'chr1', $pos, '.', $ref, $alt, '.', '.',
    'SpliceAI=' . join('|', $alt, $symbol, @$values));
}
sub plugin {
  my ($records) = @_;
  my $p = FixtureSpliceAI116->new({output_format => 'vcf'},
    'snv=fixture-snv', 'indel=fixture-indel');
  $p->{fixture_rows} = [map { $p->parse_data($_) } @$records];
  return $p;
}
sub tva {
  my (%v) = @_;
  my $vf = Bio::EnsEMBL::Variation::VariationFeature->new_fast({
    chr => ($v{chr} // 'chr1'), start => $v{pos},
    end => ($v{end} // ($v{pos} + length($v{ref}) - 1)),
    strand => ($v{vf_strand} // 1),
    allele_string => join('/', $v{ref}, @{$v{alts} // [$v{alt}]}),
  });
  my $tr = Bio::EnsEMBL::Transcript->new(
    -start => 944203, -end => 959309, -strand => ($v{transcript_strand} // -1));
  $tr->{_gene_symbol} = $v{symbol};
  $tr->{_gene_hgnc} = $v{hgnc};
  my $tv = Bio::EnsEMBL::Variation::TranscriptVariation->new_fast({
    base_variation_feature => $vf, feature => $tr,
  });
  # Keep the real overlap object strongly held for this isolated boundary call.
  return Bio::EnsEMBL::Variation::TranscriptVariationAllele->new_fast({
    base_variation_feature_overlap => $tv, allele_number => ($v{allele_number} // 1),
  }, 1);
}
sub expected {
  my ($symbol, $values) = @_;
  my @names = qw(DS_AG DS_AL DS_DG DS_DL DP_AG DP_AL DP_DG DP_DL);
  my %result = (SpliceAI_pred_SYMBOL => $symbol);
  @result{map { "SpliceAI_pred_$_" } @names} = @$values;
  return \%result;
}
my $edge = ['0.25', '0.04', '0.02', '-0.0', '2', '-7', '4', '0'];
my $other = ['0.7', '0.1', '0.05', '0.01', '8', '-6', '10', '-2'];

subtest 'official 116 implementation and actual API provenance' => sub {
  is(digest($plugin_file), $plugin_sha, 'official release 116 plugin bytes are pinned');
  is(Bio::EnsEMBL::ApiVersion::software_version(), 116, 'actual core API is 116');
  for my $module (sort keys %api_sha) {
    is(digest($INC{$module}), $api_sha{$module}, "actual API bytes: $module");
  }
  is(FixtureSpliceAI116->can('run'), SpliceAI->can('run'), 'run is inherited unchanged');
  is(FixtureSpliceAI116->can('parse_data'), SpliceAI->can('parse_data'), 'parser is inherited unchanged');
  my $allele = tva(pos => 957127, ref => 'GG', alt => 'AA', symbol => 'NOC2L');
  isa_ok($allele, 'Bio::EnsEMBL::Variation::TranscriptVariationAllele');
  isa_ok($allele->variation_feature, 'Bio::EnsEMBL::Variation::VariationFeature');
  is_deeply($allele->base_variation_feature->alt_alleles, ['AA'], 'actual allele API exposes the intact MNV');
  diag("Perl=$^V plugin=$plugin_file sha256=$plugin_sha API=116 lib=$vep_lib");
  diag('Only add_file/get_data are fixtures. No full VEP, cache/reference lookup or real tabix retrieval is certified.');
};

subtest 'intact 957127 GG>AA on the negative model/transcript strand' => sub {
  my $p = plugin([record(957127, 'GG', 'AA', 'NOC2L', $edge)]);
  is_deeply($p->run(tva(pos=>957127, ref=>'GG', alt=>'AA', symbol=>'NOC2L')),
    expected('NOC2L', $edge), 'intact forward-genomic MNV keeps all eight lexical values');
  is_deeply($p->{fixture_query}, ['chr1',957127,957128], 'real plugin queries the complete two-base interval');
  is_deeply($p->run(tva(pos=>957127, ref=>'CC', alt=>'TT', symbol=>'NOC2L')),
    {}, 'model/transcript strand does not reverse-complement source alleles');
  is_deeply($p->run(tva(pos=>957127, ref=>'GG', alt=>'AT', symbol=>'NOC2L')),
    {}, 'wrong intact alternate remains missing');
  is_deeply($p->run(tva(pos=>957127, ref=>'GA', alt=>'AA', symbol=>'NOC2L')),
    {}, 'wrong reference remains missing');
  is_deeply($p->run(tva(pos=>957128, ref=>'GG', alt=>'AA', symbol=>'NOC2L')),
    {}, 'wrong locus remains missing even if supplied as a candidate');
  is_deeply($p->run(tva(pos=>957127, ref=>'GG', alt=>'AA', symbol=>'OTHER')),
    {}, 'wrong gene remains missing');
  is_deeply($p->run(tva(pos=>957127, ref=>'GG', alt=>'AA', symbol=>'noc2l')),
    {}, 'gene symbols match exactly including case');
  is_deeply($p->run(tva(pos=>957127, ref=>'GG', alt=>'AA', hgnc=>'NOC2L')),
    expected('NOC2L', $edge), 'HGNC fallback is retained');
};

subtest 'VF orientation and allele number are independent of transcript strand' => sub {
  my $p = plugin([record(957127,'GG','AA','NOC2L',$edge),
                  record(957127,'GG','TT','NOC2L',$other)]);
  is_deeply($p->run(tva(pos=>957127, ref=>'CC', alt=>'TT', vf_strand=>-1, symbol=>'NOC2L')),
    expected('NOC2L', $edge), 'negative VF orientation alone restores genomic GG>AA');
  is_deeply($p->run(tva(pos=>957127, ref=>'GG', alt=>'AA', transcript_strand=>1, symbol=>'NOC2L')),
    expected('NOC2L', $edge), 'official plugin itself does not validate model/transcript strand');
  is_deeply($p->run(tva(pos=>957127, ref=>'GG', alts=>['AA','TT'], allele_number=>2, symbol=>'NOC2L')),
    expected('NOC2L', $other), 'actual 116 allele_number API selects the second intact allele');
};

subtest 'release 116 requires exact representations unlike historical 113 normalization' => sub {
  my $p = plugin([record(957126,'CGG','CAA','NOC2L',$edge)]);
  is_deeply($p->run(tva(pos=>957127, ref=>'GG', alt=>'AA', symbol=>'NOC2L')),
    {}, 'padded versus minimal MNV is not normalized by 116');
  is_deeply($p->run(tva(pos=>957126, ref=>'CGG', alt=>'CAA', symbol=>'NOC2L')),
    expected('NOC2L', $edge), 'same exact padded representation matches');
  note('GG>AA at 957127 is already minimal. The publisher MNV stays intact; no decomposition, normalization or model-strand flip is applied.');
};

if (defined $projected_vcf) {
  subtest 'bounded actual projected negative-model-strand MNV' => sub {
    open my $fh, '<', $projected_vcf or die "Cannot open projected VCF: $!\n";
    my (@matching, $lines, $bytes);
    while (my $line = <$fh>) {
      die "Projected fixture exceeds 100000 lines/32 MiB\n" if ++$lines > 100000 || ($bytes += length($line)) > 32*1024*1024;
      next if $line =~ /^#/;
      chomp $line;
      my @fields = split /\t/, $line;
      next unless @fields == 8 && $fields[0] =~ /^(?:chr)?1$/ && $fields[1] eq '957127'
        && $fields[3] eq 'GG' && $fields[4] eq 'AA';
      my @info = split /\|/, $fields[7];
      next unless @info == 10 && $info[0] eq 'SpliceAI=AA' && $info[1] eq 'NOC2L';
      push @matching, [$line, [@info[2..9]]];
    }
    close $fh;
    is(scalar @matching, 1, 'exactly one intact NOC2L MNV projection at 957127');
    return unless @matching == 1;
    my ($line, $values) = @{$matching[0]};
    my $p = plugin([$line]);
    is_deeply($p->run(tva(pos=>957127, ref=>'GG', alt=>'AA', symbol=>'NOC2L')),
      expected('NOC2L',$values), 'all eight actual projected values survive the official 116 boundary unchanged');
    is_deeply($p->run(tva(pos=>957127, ref=>'CC', alt=>'TT', symbol=>'NOC2L')),
      {}, 'negative model strand does not flip the actual MNV');
    is_deeply($p->run(tva(pos=>957127, ref=>'GG', alt=>'AA', symbol=>'OTHER')),
      {}, 'actual projected scores do not leak to another gene');
    diag('Actual scores are read from the bounded local artifact. Source-to-projection score provenance is checked by validate_real_mnv.py separately.');
  };
}

done_testing();
