#!/usr/bin/env perl
use strict;
use warnings;
use Digest::SHA ();
use Getopt::Long qw(GetOptions);
use JSON::PP ();

# Exercise the unchanged official constructor/header parser on actual tiny
# BGZF+tabix files. No plugin method, tabix call or returned field is mocked.
my ($plugin_file, $vep_lib, $snv, $indel);
GetOptions('plugin-file=s' => \$plugin_file, 'vep-lib=s' => \$vep_lib,
           'snv=s' => \$snv, 'indel=s' => \$indel) or die "Invalid options\n";
die "Require pinned plugin/API and both indexed files\n"
  unless defined($plugin_file) && defined($vep_lib) && defined($snv) && defined($indel) && !@ARGV;
open my $fh, '<:raw', $plugin_file or die "Cannot read plugin: $!\n";
die "Official release116 plugin identity mismatch\n"
  unless Digest::SHA->new(256)->addfile($fh)->hexdigest eq
    'd08ce7d8c2e6b3229cef638beff145e465e224f0e3fcd43fd49d501ebc2740b1';
unshift @INC, $vep_lib;
require $plugin_file;
my $plugin = SpliceAI->new({output_format => 'vcf'},
  "snv=$snv", "indel=$indel", 'cutoff=0.5');
my $headers = $plugin->get_header_info();
print JSON::PP->new->canonical->encode($headers), "\n";
