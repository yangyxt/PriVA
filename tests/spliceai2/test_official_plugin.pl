#!/usr/bin/env perl
use strict;
use warnings;
use Digest::SHA ();
use Getopt::Long qw(GetOptions);
use File::Spec;
use Test::More;

# Boundary test, not an end-to-end VEP or cache qualification. All records below
# are invented. Load the installed official plugin and actual Ensembl matcher;
# replace only tabix I/O and the small transcript/VariationFeature object API.
# Usage (with the pinned priva_acmg Perl and its bin directory on PATH):
#   perl tests/spliceai2/test_official_plugin.pl \
#     --plugin-file /path/to/VEP_plugin/SpliceAI.pm \
#     --vep-lib /path/to/share/ensembl-vep-113.2-0
my ($plugin_file, $vep_lib);
GetOptions('plugin-file=s' => \$plugin_file, 'vep-lib=s' => \$vep_lib)
  or die "Invalid options\n";
die "Require --plugin-file and --vep-lib\n"
  unless defined $plugin_file && defined $vep_lib && !@ARGV;
$plugin_file = File::Spec->rel2abs($plugin_file);
$vep_lib = File::Spec->rel2abs($vep_lib);
my $plugin_sha = '500a534c403bcbcc15753718768f0a1f352cc6b93f9119e5eba25745f7205b9a';
my $matcher_sha = 'a85aa34bb337798dd45599270340aa56740f2d53496be1590e610480ff8f6d2e';
sub digest {
  open my $fh, '<:raw', $_[0] or die "Cannot open $_[0]: $!\n";
  return Digest::SHA->new(256)->addfile($fh)->hexdigest;
}
die "Official SpliceAI.pm identity mismatch\n" unless digest($plugin_file) eq $plugin_sha;
unshift @INC, $vep_lib;
require $plugin_file;
my $matcher_file = $INC{'Bio/EnsEMBL/Variation/Utils/Sequence.pm'};
die "Actual Ensembl matcher identity mismatch\n"
  unless defined $matcher_file && digest($matcher_file) eq $matcher_sha;

{
  package FixtureSpliceAI;
  our @ISA = ('SpliceAI');
  sub add_file { push @{$_[0]->{fixture_files}}, $_[1] }
  sub get_data {
    my ($self, @query) = @_;
    $self->{fixture_query} = \@query;
    return $self->{fixture_rows};
  }
  package FixtureVariationFeature;
  sub ref_allele_string { $_[0]->{ref} }
  sub strand { $_[0]->{strand} }
  package FixtureTranscriptAllele;
  sub variation_feature { $_[0]->{vf} }
  sub variation_feature_seq { $_[0]->{alt} }
  sub transcript { $_[0]->{transcript} }
}

sub record {
  my ($pos, $ref, $alt, $symbol, $values) = @_;
  die "Synthetic record must contain eight values\n" unless @$values == 8;
  return join("\t", '1', $pos, '.', $ref, $alt, '.', '.',
    'SpliceAI=' . join('|', $alt, $symbol, @$values));
}
sub plugin {
  my ($records, $cutoff) = @_;
  my @args = ('snv=synthetic-snv', 'indel=synthetic-indel');
  push @args, "cutoff=$cutoff" if defined $cutoff;
  my $p = FixtureSpliceAI->new({output_format => 'vcf'}, @args);
  $p->{fixture_rows} = [map { $p->parse_data($_) } @$records];
  return $p;
}
sub tva {
  my (%v) = @_;
  return bless {
    vf => bless({chr => ($v{chr} // '1'), start => $v{pos},
      end => ($v{end} // $v{pos}), strand => ($v{vf_strand} // 1), ref => $v{ref}},
      'FixtureVariationFeature'),
    alt => $v{alt},
    transcript => {_gene_symbol => $v{symbol}, _gene_hgnc => $v{hgnc},
      strand => ($v{transcript_strand} // 1)},
  }, 'FixtureTranscriptAllele';
}
sub expected {
  my ($symbol, $v, $flag) = @_;
  my @names = qw(DS_AG DS_AL DS_DG DS_DL DP_AG DP_AL DP_DG DP_DL);
  my %r = (SpliceAI_pred_SYMBOL => $symbol);
  @r{map { "SpliceAI_pred_$_" } @names} = @$v;
  $r{SpliceAI_cutoff} = $flag if defined $flag;
  return \%r;
}
my $edge = ['0.25', '0.04', '0.02', '-0.0', '2', '-7', '4', '0'];
my $low  = ['0.01', '0.02', '0.03', '0.04', '-3', '5', '7', '-9'];
my $high = ['0.7', '0.1', '0.05', '0.01', '8', '-6', '10', '-2'];

subtest 'real implementation provenance' => sub {
  is(digest($plugin_file), $plugin_sha, 'plugin bytes match audited release-113 commit 10054085');
  is(digest($matcher_file), $matcher_sha, 'installed VEP 113.2 matcher bytes are pinned');
  is(\&SpliceAI::get_matched_variant_alleles,
     \&Bio::EnsEMBL::Variation::Utils::Sequence::get_matched_variant_alleles,
     'plugin imports the actual Ensembl matcher, without replacement');
  is(FixtureSpliceAI->can('run'), SpliceAI->can('run'), 'run is inherited unchanged');
  is(FixtureSpliceAI->can('parse_data'), SpliceAI->can('parse_data'), 'parser is inherited unchanged');
  diag("Perl=$^V plugin=$plugin_file sha256=$plugin_sha matcher=$matcher_file sha256=$matcher_sha");
  diag('Only file registration/tabix get_data and light VEP API objects are fixtures; no real VEP cache, reference lookup or BGZF index is exercised.');
};

subtest 'one-gene INFO record and all eight field positions' => sub {
  my $line = record(101, 'A', 'C', 'GENE_A', $edge);
  my @cols = split /\t/, $line;
  is(scalar @cols, 8, 'projected synthetic VCF has eight fields');
  unlike($cols[7], qr/[,;]/, 'one SpliceAI value and no extra INFO tags');
  is(scalar(split /\|/, $cols[7]), 10, 'one allele, one symbol and eight numbers');
  my $p = plugin([$line]);
  is_deeply($p->run(tva(pos=>101, ref=>'A', alt=>'C', symbol=>'GENE_A')),
    expected('GENE_A', $edge), 'field order and lexical score/distance values survive real plugin');
  is_deeply($p->{fixture_query}, ['1',101,101], 'actual plugin requests the expected locus');
};

subtest 'two alternate alleles at the same position' => sub {
  my $p = plugin([record(102,'A','C','GENE_A',$edge), record(102,'A','T','GENE_A',$high)]);
  is_deeply($p->run(tva(pos=>102,ref=>'A',alt=>'C',symbol=>'GENE_A')), expected('GENE_A',$edge), 'C selects C values');
  is_deeply($p->run(tva(pos=>102,ref=>'A',alt=>'T',symbol=>'GENE_A')), expected('GENE_A',$high), 'T selects T values');
  is_deeply($p->run(tva(pos=>102,ref=>'A',alt=>'G',symbol=>'GENE_A')), {}, 'absent alternate does not inherit another score');
};

subtest 'gene association and transcript strand boundary' => sub {
  my $p = plugin([record(103,'G','A','GENE_PLUS',$edge), record(103,'G','A','GENE_MINUS',$low)]);
  is_deeply($p->run(tva(pos=>103,ref=>'G',alt=>'A',symbol=>'GENE_PLUS',transcript_strand=>1)), expected('GENE_PLUS',$edge), 'plus transcript receives its exact symbol record');
  is_deeply($p->run(tva(pos=>103,ref=>'G',alt=>'A',symbol=>'GENE_MINUS',transcript_strand=>-1)), expected('GENE_MINUS',$low), 'minus transcript receives its exact symbol record');
  is_deeply($p->run(tva(pos=>103,ref=>'G',alt=>'A',symbol=>'GENE_PLUS',transcript_strand=>-1)), expected('GENE_PLUS',$edge), 'documented limitation: plugin itself does not validate model/transcript strand');
  is_deeply($p->run(tva(pos=>103,ref=>'G',alt=>'A',symbol=>'gene_plus')), {}, 'symbol matching is case-sensitive');
  is_deeply($p->run(tva(pos=>103,ref=>'G',alt=>'A',symbol=>'OTHER')), {}, 'nonmatching gene stays missing');
  is_deeply($p->run(tva(pos=>103,ref=>'G',alt=>'A',hgnc=>'GENE_PLUS')), expected('GENE_PLUS',$edge), 'official HGNC fallback is retained');
  is_deeply($p->run(tva(pos=>103,ref=>'G',alt=>'A',symbol=>'OTHER',hgnc=>'GENE_PLUS')), {}, 'present nonmatching gene symbol does not fall back to HGNC');
};

subtest 'anchored insertion and deletion' => sub {
  my $ins = plugin([record(200,'A','AT','GENE_A',$edge)]);
  is_deeply($ins->run(tva(pos=>200,ref=>'A',alt=>'AT',symbol=>'GENE_A')), expected('GENE_A',$edge), 'anchored insertion uses real allele matching');
  my $del = plugin([record(300,'AT','A','GENE_A',$low)]);
  is_deeply($del->run(tva(pos=>300,end=>301,ref=>'AT',alt=>'A',symbol=>'GENE_A')), expected('GENE_A',$low), 'anchored deletion uses real allele matching');
};

subtest 'minimal allele representations use real normalization' => sub {
  my $ins = plugin([record(200,'A','AT','GENE_A',$edge)]);
  is_deeply($ins->run(tva(pos=>201,end=>200,ref=>'-',alt=>'T',symbol=>'GENE_A')), expected('GENE_A',$edge), 'minimal insertion runs actual VariationFeature to_VCF_record and matcher');
  is_deeply($ins->{fixture_query}, ['1',200,201], 'insertion query bounds are swapped by actual plugin');
  my $del = plugin([record(300,'AT','A','GENE_A',$low)]);
  is_deeply($del->run(tva(pos=>301,ref=>'T',alt=>'-',symbol=>'GENE_A')), expected('GENE_A',$low), 'minimal deletion matches anchored cache allele');
  my $padded = plugin([record(400,'GA','GT','GENE_A',$high)]);
  is_deeply($padded->run(tva(pos=>401,ref=>'A',alt=>'T',symbol=>'GENE_A')), expected('GENE_A',$high), 'padded versus minimal substitution is normalized by actual matcher');
  note('get_data supplies candidates: these checks do not prove shifted candidate retrieval from a real index.');
};

subtest 'VariationFeature orientation is separate from transcript/model strand' => sub {
  my $p = plugin([record(500,'A','C','GENE_A',$edge)]);
  is_deeply($p->run(tva(pos=>500,ref=>'T',alt=>'G',vf_strand=>-1,symbol=>'GENE_A')), expected('GENE_A',$edge), 'actual matcher reverse-complements negative VF alleles');
  is_deeply($p->run(tva(pos=>500,ref=>'T',alt=>'G',vf_strand=>1,symbol=>'GENE_A')), {}, 'same alleles with positive VF strand do not match');
};

subtest 'cutoff is an inclusive flag and never filters score values' => sub {
  my $p = plugin([record(600,'C','T','GENE_A',$edge), record(601,'C','T','GENE_A',$low)], '0.25');
  # The tabix fixture returns both positions; the real matcher rejects the other locus.
  is_deeply($p->run(tva(pos=>600,ref=>'C',alt=>'T',symbol=>'GENE_A')), expected('GENE_A',$edge,'PASS'), 'exact 0.25 boundary passes and retains -0.0/low values');
  is_deeply($p->run(tva(pos=>601,ref=>'C',alt=>'T',symbol=>'GENE_A')), expected('GENE_A',$low,'FAIL'), 'lower values remain present with FAIL flag');
  my $plain = plugin([record(601,'C','T','GENE_A',$low)]);
  is_deeply($plain->run(tva(pos=>601,ref=>'C',alt=>'T',symbol=>'GENE_A')), expected('GENE_A',$low), 'no cutoff keeps every value and adds no flag');
};

subtest 'one projected row per allele and symbol avoids overwrite' => sub {
  my @rows = (record(700,'A','G','GENE_A',$edge), record(700,'A','G','GENE_B',$low));
  my $p = plugin(\@rows);
  is(scalar @{$p->{fixture_rows}}, 2, 'two genes, not three top-k records per gene');
  is_deeply($p->run(tva(pos=>700,ref=>'A',alt=>'G',symbol=>'GENE_A')), expected('GENE_A',$edge), 'first gene survives distinct-symbol record');
  is_deeply($p->run(tva(pos=>700,ref=>'A',alt=>'G',symbol=>'GENE_B')), expected('GENE_B',$low), 'second gene has independent values');
  my $bad = plugin([@rows,record(700,'A','G','GENE_A',$high)]);
  is_deeply($bad->run(tva(pos=>700,ref=>'A',alt=>'G',symbol=>'GENE_A')), expected('GENE_A',$high), 'hazard documented: duplicate same-gene records overwrite; producer must reject them');
};

subtest 'unsupported packed INFO formats are not silently certified' => sub {
  my $line = record(800,'A','C','GENE_A',$edge);
  my $extra = plugin([$line.';EXTRA=1']);
  my $result = $extra->run(tva(pos=>800,ref=>'A',alt=>'C',symbol=>'GENE_A'));
  is($result->{SpliceAI_pred_DP_DL}, '0;EXTRA=1', 'extra INFO contaminates final field: projection must contain only SpliceAI');
  my $packed = plugin([$line.',C|GENE_B|'.join('|',@$low)]);
  is_deeply($packed->run(tva(pos=>800,ref=>'A',alt=>'C',symbol=>'GENE_B')), {}, 'comma-packed second gene is lost: separate VCF records are required');
};

subtest 'absence remains absence' => sub {
  my $p = plugin([]);
  is_deeply($p->run(tva(pos=>900,ref=>'A',alt=>'T',symbol=>'GENE_A')), {}, 'no candidate produces no evidence rather than zeros');
};

done_testing();
