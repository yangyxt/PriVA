package CaddSelection;
# Thin adapter around installed official VEP methods; no scoring/ranking port.
# Production use requires independent reference qualification for this exact
# plugin, official VEP source, cache/profile and explicit Perl hash environment.
use strict;
use warnings;
use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepPlugin);
use Bio::EnsEMBL::VEP::OutputFactory;
use Bio::EnsEMBL::VEP::Constants;
use Scalar::Util qw(refaddr);
use Digest::SHA qw(sha256_hex);
use JSON::PP;

sub file_sha {
    my ($path) = @_;
    open(my $handle, '<', $path) or die "CaddSelection cannot read $path: $!";
    binmode($handle);
    return Digest::SHA->new(256)->addfile($handle)->hexdigest;
}

sub selection_context {
    my ($order, $config, $sources) = @_;
    return {
        adapter_id => 'official-vep-native-per-gene-v1',
        plugin_sha256 => file_sha(__FILE__),
        output_factory_sha256 => file_sha($INC{'Bio/EnsEMBL/VEP/OutputFactory.pm'}),
        vep_release => "$Bio::EnsEMBL::VEP::Constants::VEP_VERSION",
        perl_version => "$^V",
        perl_hash_seed => $ENV{PERL_HASH_SEED},
        perl_perturb_keys => $ENV{PERL_PERTURB_KEYS},
        pick_order => [@$order],
        cadd_transcript_sources => [@$sources],
        # Merged/RefSeq pools and transcript-reference correction differ from
        # the original CADD reference defaults. They are explicit profile
        # inputs requiring qualification, never silently treated as equivalent.
        # 24 October 2026 CT: official distance is an ARRAY. Copy its exact
        # ordered values into provenance; stringification records a process
        # address and prevents durable exact-context validation. Config and
        # native selection are untouched; old address contexts remain rejected.
        annotation_options => {map {$_ =>
            $_ eq 'distance' && ref($config->{$_}) eq 'ARRAY' ? [@{$config->{$_}}]
            : defined($config->{$_}) ? "$config->{$_}" : ''}
            qw(assembly cache_version db_version merged refseq use_transcript_ref minimal distance
               regulatory sift polyphen ccds domains numbers canonical total_length)},
    };
}

sub new {
    my $class = shift;
    my $self = $class->SUPER::new(@_);
    my $config = $self->config;
    die "CaddSelection requires explicit PERL_HASH_SEED and PERL_PERTURB_KEYS=0 before Perl starts"
        unless defined($ENV{PERL_HASH_SEED}) && length($ENV{PERL_HASH_SEED})
            && defined($ENV{PERL_PERTURB_KEYS}) && $ENV{PERL_PERTURB_KEYS} eq '0';
    for my $option (qw(pick pick_allele pick_allele_gene per_gene most_severe summary
                       coding_only no_intergenic process_ref_homs individual individual_zyg)) {
        die "CaddSelection requires all overlap rows; incompatible VEP option $option"
            if $config->{$option};
    }
    die 'CaddSelection requires the resolved official VEP pick_order array'
        unless ref($config->{pick_order}) eq 'ARRAY' && @{$config->{pick_order}};
    my $source_parameter = $self->params_to_hash->{cadd_transcript_sources} || '';
    die 'CaddSelection requires explicit cadd_transcript_sources=Ensembl or Ensembl+RefSeq; neither is qualified by default'
        unless $source_parameter eq 'Ensembl' || $source_parameter eq 'Ensembl+RefSeq';
    $self->{cadd_transcript_sources} = [split(/\+/, $source_parameter)];
    # These exact native methods only need pick_order and the already cached
    # objects. No new annotation/database query or second VEP process is used.
    $self->{selector} = bless {pick_order => [@{$config->{pick_order}}]},
                             'Bio::EnsEMBL::VEP::OutputFactory';
    $self->{context} = selection_context($config->{pick_order}, $config, $self->{cadd_transcript_sources});
    $self->{context_sha256} = sha256_hex(JSON::PP->new->canonical->encode($self->{context}));
    return $self;
}

sub feature_types { return ['Transcript', 'RegulatoryFeature', 'MotifFeature', 'Intergenic']; }
sub get_header_info {
    return {
        CADD_PICK_ORDINAL => 'Native official VEP per_gene output ordinal; zero means retained PriVA-only row',
        CADD_PICK_COUNT => 'Number of native per_gene rows before CADD motif aggregation',
        CADD_PICK_CONTEXT => 'SHA256 of native selection source, pick_order and explicit Perl hash environment',
    };
}

sub run {
    my ($self, $vfoa) = @_;
    my $vf = $vfoa->base_variation_feature;
    die 'CaddSelection requires an explicitly biallelic variant'
        unless scalar(split('/', $vf->allele_string)) == 2;
    my $cache = $vf->{_priva_cadd_native_selection};
    unless ($cache) {
        my $all = $self->{selector}->get_all_VariationFeatureOverlapAlleles($vf);
        die 'CaddSelection found no native overlap objects' unless @$all;
        my %sources = map {$_ => 1} @{$self->{cadd_transcript_sources}};
        my @pool;
        for my $object (@$all) {
            # Never let the native method fall back to an online GeneAdaptor.
            die 'CaddSelection requires cached native gene stable IDs'
                if $object->isa('Bio::EnsEMBL::Variation::BaseTranscriptVariationAllele')
                    && !$object->transcript->{_gene_stable_id};
            if ($object->isa('Bio::EnsEMBL::Variation::BaseTranscriptVariationAllele')) {
                my $source = $object->transcript->{_source_cache} || '';
                die 'CaddSelection requires explicit native merged-cache transcript source'
                    unless $source eq 'Ensembl' || $source eq 'RefSeq';
                next unless $sources{$source};
            }
            push @pool, $object;
        }
        # RefSeq-only overlap can leave the Ensembl view empty, while a separate
        # Ensembl-only VEP might emit intergenic. Do not invent that missing row.
        die 'CaddSelection source pool has no native rows; intergenic representation requires qualification'
            unless @pool;
        my $selected = $self->{selector}->pick_VariationFeatureOverlapAllele_per_gene(\@pool);
        my %ordinal = map {refaddr($_) => 0} @$all;
        for (my $index = 0; $index < @$selected; $index++) {
            $ordinal{refaddr($selected->[$index])} = $index + 1;
        }
        $cache = $vf->{_priva_cadd_native_selection} = {
            ordinal => \%ordinal, count => scalar(@$selected), context => $self->{context_sha256}};
    }
    die 'CaddSelection context changed within a variant' unless $cache->{context} eq $self->{context_sha256};
    die 'CaddSelection output object is absent from the native overlap set'
        unless exists($cache->{ordinal}->{refaddr($vfoa)});
    # Never filter a PriVA output row. The separate view adapter validates the
    # complete ordinal set, catching plugin warnings and later plugin filtering.
    return {CADD_PICK_ORDINAL => $cache->{ordinal}->{refaddr($vfoa)},
            CADD_PICK_COUNT => $cache->{count}, CADD_PICK_CONTEXT => $cache->{context}};
}

1;
