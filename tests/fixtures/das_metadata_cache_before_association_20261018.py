# Exact former cache block, source SHA c4c12539f8ad7e1c9fa134687721359d90a58d85dab150b305185e8a94917269.
# Reference-only success-byte baseline; never an installation route.
domain_metadata = None
if am_vcf is not None:
    domain_meta_cache = os.path.join(output_dir, f'domain_metadata.{assembly}.pkl')
    if os.path.exists(domain_meta_cache):
        logger.info(f"[DOMAIN_META] Loading cached metadata from {domain_meta_cache}")
        with open(domain_meta_cache, 'rb') as f:
            domain_metadata = pickle.load(f)
        # Validate cache has prot_length
        sample_meta = next(iter(domain_metadata.values()), {}) if domain_metadata else {}
        if sample_meta.get('prot_length') is None:
            logger.warning(f"[DOMAIN_META] Stale cache (no prot_length), regenerating...")
            domain_metadata = None
    if domain_metadata is None:
        logger.info(f"[DOMAIN_META] Extracting from AM VCF: {am_vcf}")
        domain_metadata = extract_domain_metadata_parallel(am_vcf, threads=threads)
        with open(domain_meta_cache, 'wb') as f:
            pickle.dump(domain_metadata, f)
        logger.info(f"[DOMAIN_META] Cached {len(domain_metadata)} domains to {domain_meta_cache}")
