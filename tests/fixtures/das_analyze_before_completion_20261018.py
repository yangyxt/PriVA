# Exact original analyze_domain_data before explicit completion status.
# Source SHA256: be892d70c7d52e2b7ad320e285d4875fd0275b6e75009ab9d3888d6dc6bb51a8
# Reference-only fixture, never an installation route.
def analyze_domain_data(pickle_file: str,
                       output_dir: str,
                       threads: int = 12,
                       fdr_threshold: float = 0.01,
                       assembly: str = 'hg19',
                       am_vcf: str = None,
                       hcseeker_spots_tsv: str = None,
                       clinvar_vcf: str = None,
                       coldspot_overlap_threshold: float = 0.6,
                       common_af_threshold: float = 0.05,
                       interpro_mapping_pickle: str = None):
    """
    Analyze domain data from the pickle file and apply post-KS filters.

    Args:
        pickle_file: Path to the pickle file containing scores_dict
        output_dir: Directory to save analysis outputs
        threads: Number of threads for parallel processing
        fdr_threshold: FDR threshold for significance after correction
        assembly: Assembly version (hg19 or hg38)
        am_vcf: Path to VEP-annotated AlphaMissense VCF (needed for Filter 1)
        hcseeker_spots_tsv: Path to HC_spots.{assembly}.tsv (None = skip Filter 1)
        clinvar_vcf: Path to clinvar.{assembly}.vep.vcf.gz (None = skip Filter 2)
        coldspot_overlap_threshold: Overlap fraction threshold for Filter 1 (default 0.5)
        common_af_threshold: AF threshold for Filter 3 (default 0.05)
        interpro_mapping_pickle: Selected mapping path (None retains the bundled default)
    """
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Load the data
    with open(pickle_file, 'rb') as f:
        scores_dict = pickle.load(f)

    # Generate transcript-protein position-domain mapping
    transcript_map = generate_transcript_pp_domain_map(scores_dict)
    
    # Save the mapping to both JSON and pickle files
    mapping_json = os.path.join(output_dir, f'transcript_pp_domain_mapping.{assembly}.json')
    mapping_pickle = os.path.join(output_dir, f'transcript_pp_domain_mapping.{assembly}.pkl')
    
    with open(mapping_json, 'w') as f:
        json.dump(transcript_map, f, indent=2)
    with open(mapping_pickle, 'wb') as f:
        pickle.dump(transcript_map, f)
    
    logger.info(f"Saved transcript-pp-domain mapping to {mapping_json} and {mapping_pickle}")
    
    # Continue with existing analysis
    domain_scores = visualize_domain_distribution(scores_dict, output_dir, threads, assembly)
    dm_instance = DomainNormalizer()

    # Standalone callers retain the original bundle; explicit selections never
    # fall back to it. Compression comes from bytes, not the legacy suffix.
    interpro_map_pickle = (os.path.join(os.path.dirname(self_directory), 'data', 'InterPro', 'Interpro_entry_mapping.pkl.gz')
                           if interpro_mapping_pickle is None else interpro_mapping_pickle)
    interpro_map_dict = load_interpro_mapping(interpro_map_pickle)

    functional_map = os.path.join(os.path.dirname(self_directory), 'data', 'InterPro', 'curated_InterPro_func_domains.tsv.gz')
    func_map_df = pd.read_table(functional_map, low_memory=False)
    func_map_dict = dict(zip(func_map_df['IPR_ID'], func_map_df['Molecular_Function_GO_Terms']))
    func_pred_dict = dict(zip(func_map_df['IPR_ID'], func_map_df['Functionality_Assessment']))

    # --- Load/cache domain metadata (needed for ref set filtering and post-KS filters) ---
    domain_metadata = None
    if am_vcf is not None:
        domain_meta_cache = os.path.join(output_dir, f'domain_metadata.{assembly}.pkl')
        domain_metadata = load_or_extract_domain_metadata(
            cache_path=domain_meta_cache, vcf_path=am_vcf, assembly=assembly,
            extractor=lambda: extract_domain_metadata_parallel(am_vcf, threads=threads),
            source_paths=(__file__, os.path.join(self_directory, 'stat_protein_domain_amscores.py')),
            logger=logger)

    # Pick out the reference (curated functional) domains
    ref_domain_scores = {
        domain: domain_s
        for domain, domain_s in domain_scores.items()
        if identify_functional_domain(domain, dm_instance, func_map_dict, func_pred_dict, interpro_map_dict,
                                      domain_metadata=domain_metadata) == "Functional"
    }
    logger.info(f"Found {len(ref_domain_scores)} curated functional domains/sites as references")

    # Prepare the list of query domains (exclude the curated references)
    query_domains = [(d, s) for d, s in domain_scores.items() if d not in ref_domain_scores]
    logger.info(f"Found {len(query_domains)} query domains to test against the composite reference")

    if len(ref_domain_scores) == 0 or len(query_domains) == 0:
        logger.warning("No references or no query domains found; skipping KS analysis.")
        return

    # --- Build a single fixed equal-weight composite reference ONCE ---
    # Choose composite size T: median of query sizes, capped to save memory
    query_sizes = np.array([len(s) for _, s in query_domains], dtype=int)
    if len(query_sizes) == 0:
        logger.warning("No query sizes found; skipping KS analysis.")
        return
    T_default = int(np.max(query_sizes))
    T_cap = 200000
    T = max(50, min(T_default, T_cap))  # keep it reasonable
    seed_global = 42  # make reproducible
    rng_global = np.random.default_rng(seed_global)

    ref_arrays = list(ref_domain_scores.values())
    fixed_composite = build_fixed_equal_weight_composite(ref_arrays, T=T, rng=rng_global)
    logger.info(f"Built fixed equal-weight composite reference of size T={len(fixed_composite)} from {len(ref_arrays)} references")

    # --- Classify ref domains by InterPro type for stratified KS ---
    ref_by_type = {}
    for dp in ref_domain_scores:
        domain_name = dp.split(':', 1)[-1]
        interpro_result = dm_instance.query_interpro_entry_vep_anno(domain_name, interpro_map_dict)
        if interpro_result is None:
            ref_by_type.setdefault(None, []).append(ref_domain_scores[dp])
            continue
        entries = interpro_result.get('interpro_entries')
        itype = entries[0][1] if entries else None
        ref_by_type.setdefault(itype, []).append(ref_domain_scores[dp])

    # Build type-specific composites (only for types with enough refs)
    type_composites = {}
    MIN_REF_FOR_TYPE = 100
    for itype, arrays in ref_by_type.items():
        if itype is None or len(arrays) < MIN_REF_FOR_TYPE:
            continue
        T_type = min(T_cap, sum(len(a) for a in arrays))
        type_composites[itype] = build_fixed_equal_weight_composite(arrays, T=max(50, T_type), rng=rng_global)
    logger.info(f"Built type-specific composites for: {list(type_composites.keys())}")

    # --- Classify query domains by InterPro type ---
    query_type_map = {}
    for d, _ in query_domains:
        domain_name = d.split(':', 1)[-1]
        interpro_result = dm_instance.query_interpro_entry_vep_anno(domain_name, interpro_map_dict)
        if interpro_result is None:
            query_type_map[d] = None
            continue
        entries = interpro_result.get('interpro_entries')
        query_type_map[d] = entries[0][1] if entries else None

    # --- Parallel KS over query domains (type-stratified) ---
    # Each query is tested against its type-specific composite if available,
    # otherwise falls back to the global composite.
    cap = T
    n_tasks = len(query_domains)
    n_cores = min(threads, n_tasks) if n_tasks > 0 else 1
    logger.info(f"Processing {n_tasks} query domains using {n_cores} cores for type-stratified KS tests")

    ks_results = {}
    tasks = []
    for i, (domain_path, scores) in enumerate(query_domains):
        seed = (seed_global + i) & 0xFFFFFFFF
        itype = query_type_map.get(domain_path)
        composite = type_composites.get(itype, fixed_composite)
        task_cap = min(cap, len(composite))
        tasks.append((domain_path, scores, composite, task_cap, seed))

    if n_tasks > 0:
        with mp.Pool(n_cores) as pool:
            for res in pool.imap_unordered(process_domain_ks, tasks):
                domain_path, result = res
                ks_results[domain_path] = result

    if not ks_results:
        logger.warning("No KS results were produced.")
        return

    # --- Collect p-values and apply BH FDR across query domains ---
    domains = list(ks_results.keys())
    raw_pvals = [ks_results[d]['p_value'] for d in domains]

    rejected, pvals_bh, _, _ = multipletests(raw_pvals, alpha=fdr_threshold, method='fdr_bh') # Benjamini-Hochberg procedure

    # --- Write results table ---
    df_rows = []
    for d, adj_p, rej in zip(domains, pvals_bh, rejected):
        row = {
            'domain': d,
            'ks_stat': ks_results[d]['ks_stat'],
            'p_value': ks_results[d]['p_value'],
            'fdr_corrected_p_value': float(adj_p),
            # significant = "more tolerant than composite" at FDR
            'is_more_tolerant': bool(rej),
            'n_query': ks_results[d]['n_query'],
            'n_comp_used': ks_results[d]['n_comp'],
            'composite_T': int(len(fixed_composite)),
            'median_query': ks_results[d]['median_query'],
            'median_composite': ks_results[d]['median_comp']
        }
        df_rows.append(row)

    results_df = pd.DataFrame(df_rows)
    out_tsv = os.path.join(output_dir, f'domain_tolerance_analysis.{assembly}.tsv')
    results_df.to_csv(out_tsv, sep='\t', index=False)
    logger.info(f"Composite KS analysis results saved to {out_tsv}")
    logger.info(f"Completed KS analysis for {len(results_df)} domains; "
                f"found {int(np.sum(rejected))} domains significantly more tolerant than composite (FDR < {fdr_threshold})")

    # --- Define "intolerant" set for downstream: those NOT significantly more tolerant ---
    not_more_tolerant = set(results_df.loc[~results_df['is_more_tolerant'], 'domain'].unique().tolist())
    analysis_intolerant_domains = set(not_more_tolerant)
    n_ks_intolerant = len(analysis_intolerant_domains)
    logger.info(f"[KS] Intolerant set (not-tolerant, before adding refs): {n_ks_intolerant}")

    # --- Log Family-type entries that survived the KS test ---
    n_family_survived = 0
    n_family_total = 0
    for domain_path in domains:
        domain_name = domain_path.split(':', 1)[-1]
        interpro_result = dm_instance.query_interpro_entry_vep_anno(domain_name, interpro_map_dict)
        if interpro_result is not None:
            entries = interpro_result.get('interpro_entries')
            if entries and entries[0][1] == 'Family':
                n_family_total += 1
                if domain_path in analysis_intolerant_domains:
                    n_family_survived += 1
    logger.info(f"[KS] Family-type entries: {n_family_total} total tested, {n_family_survived} survived into intolerant set")

    # Keep a copy of the pre-filter set for the summary TSV
    ks_intolerant_raw = set(analysis_intolerant_domains)

    # --- Post-KS filters (each filter is independent and optional) ---
    filter_summaries = {}
    excluded_sets = {}

    # domain_metadata already loaded above (before ref set building)

    # Filter 1: HCSeeker coldspot overlap
    if hcseeker_spots_tsv is not None:
        if domain_metadata is None:
            logger.warning("[FILTER1:COLDSPOT] Skipped: --am_vcf required for domain AA range extraction")
        else:
            analysis_intolerant_domains, excl, summ = filter_coldspot_overlap(
                analysis_intolerant_domains, domain_metadata, hcseeker_spots_tsv,
                dm_instance, interpro_map_dict,
                coldspot_overlap_threshold)
            filter_summaries['coldspot'] = summ
            excluded_sets['coldspot'] = excl

    # Filter 2: High-confidence benign ClinVar missense
    if clinvar_vcf is not None:
        analysis_intolerant_domains, excl, summ = filter_benign_clinvar(
            analysis_intolerant_domains, clinvar_vcf)
        filter_summaries['benign'] = summ
        excluded_sets['benign'] = excl

    # Filter 3: Common missense (gnomAD AF > threshold) — always runs
    analysis_intolerant_domains, excl, summ = filter_common_missense(
        analysis_intolerant_domains, scores_dict, common_af_threshold)
    filter_summaries['common'] = summ
    excluded_sets['common'] = excl

    # Filter 4: Whole-protein Family classifiers — requires domain_metadata
    if domain_metadata is not None:
        analysis_intolerant_domains, excl, summ = filter_whole_protein_family(
            analysis_intolerant_domains, domain_metadata,
            dm_instance, interpro_map_dict)
        filter_summaries['family'] = summ
        excluded_sets['family'] = excl
    else:
        logger.warning("[FILTER4:FAMILY] Skipped: --am_vcf required")

    # --- Bimodality rescue: recover KS-tolerant GO MF domains with significant bimodality ---
    tolerant_set = set(results_df.loc[results_df['is_more_tolerant'], 'domain'].tolist())
    rescued_domains, rescue_summary = rescue_bimodal_go_mf_domains(
        tolerant_set, domain_scores, fixed_composite,
        dm_instance, interpro_map_dict,
        type_composites=type_composites,
        domain_metadata=domain_metadata,
        threads=threads, fdr_threshold=fdr_threshold)
    analysis_intolerant_domains = analysis_intolerant_domains | rescued_domains

    # --- Rescue pass 2: non-GO-MF Domain/Repeat/Site entries with high median AM ---
    rescued_non_gomf, rescue_non_gomf_summary = rescue_non_gomf_domains(
        tolerant_set, domain_scores, fixed_composite,
        dm_instance, interpro_map_dict,
        type_composites=type_composites,
        threads=threads, fdr_threshold=fdr_threshold, min_median_am=0.75)
    analysis_intolerant_domains = analysis_intolerant_domains | rescued_non_gomf

    # --- Add curated refs AFTER all filters and rescues (refs are never filtered) ---
    analysis_intolerant_domains = analysis_intolerant_domains | set(ref_domain_scores.keys())

    # --- Summary logging ---
    n_f1 = len(excluded_sets.get('coldspot', set()))
    n_f2 = len(excluded_sets.get('benign', set()))
    n_f3 = len(excluded_sets.get('common', set()))
    n_f4 = len(excluded_sets.get('family', set()))
    n_rescued_gomf = len(rescued_domains)
    n_rescued_non_gomf = len(rescued_non_gomf)
    logger.info(f"[SUMMARY] Total domains analyzed: {len(domain_scores)}")
    logger.info(f"[SUMMARY] Curated functional references: {len(ref_domain_scores)}")
    logger.info(f"[SUMMARY] DAS intolerant (post-KS): {n_ks_intolerant}")
    logger.info(f"[SUMMARY] Removed by Filter 1 (coldspot): {n_f1}")
    logger.info(f"[SUMMARY] Removed by Filter 2 (benign): {n_f2}")
    logger.info(f"[SUMMARY] Removed by Filter 3 (common): {n_f3}")
    logger.info(f"[SUMMARY] Removed by Filter 4 (whole-protein family): {n_f4}")
    logger.info(f"[SUMMARY] Rescued GO MF (not bimodal, median AM >= 0.564): {n_rescued_gomf}")
    logger.info(f"[SUMMARY] Rescued non-GO-MF Domain/Repeat/Site (not bimodal, median AM >= 0.85): {n_rescued_non_gomf}")
    logger.info(f"[SUMMARY] Final intolerant domains: {len(analysis_intolerant_domains)}")

    # --- Save filter summary TSV ---
    # Resolve domain metadata for aa coordinates
    dm_for_summary = domain_metadata or {}
    if filter_summaries:
        summary_rows = []
        for dp in sorted(ks_intolerant_raw):
            meta = dm_for_summary.get(dp, {})
            row = {
                'domain_path': dp,
                'gene_symbol': meta.get('gene_symbol', ''),
                'aa_start': meta.get('aa_start', np.nan),
                'aa_end': meta.get('aa_end', np.nan),
                'in_original_intolerant': True,
            }
            if 'coldspot' in filter_summaries:
                s = filter_summaries['coldspot'].get(dp, {})
                row['coldspot_max_frac_of_coldspot'] = s.get('coldspot_max_frac_of_coldspot', np.nan)
                row['coldspot_max_frac_of_domain'] = s.get('coldspot_max_frac_of_domain', np.nan)
                row['coldspot_removed'] = s.get('removed', False)
            if 'benign' in filter_summaries:
                s = filter_summaries['benign'].get(dp, {})
                row['n_benign_missense'] = s.get('n_benign_missense', np.nan)
                row['benign_removed'] = s.get('removed', False)
            if 'common' in filter_summaries:
                s = filter_summaries['common'].get(dp, {})
                # Support both legacy (flat) and new (per-residue) summary fields
                row['n_common_missense'] = s.get('n_fully_common_residues', s.get('n_common_missense', np.nan))
                row['max_af'] = s.get('frac_common', s.get('max_af', np.nan))
                row['common_removed'] = s.get('removed', False)
            if 'family' in filter_summaries:
                s = filter_summaries['family'].get(dp, {})
                row['interpro_type'] = s.get('interpro_type', '')
                row['domain_coverage'] = s.get('coverage', np.nan)
                row['prot_length'] = s.get('prot_length', np.nan)
                row['family_removed'] = s.get('removed', False)
            row['final_intolerant'] = dp in analysis_intolerant_domains
            summary_rows.append(row)

        summary_df = pd.DataFrame(summary_rows)
        summary_tsv = os.path.join(output_dir, f'domain_filter_results.{assembly}.tsv')
        summary_df.to_csv(summary_tsv, sep='\t', index=False)
        logger.info(f"Saved filter summary ({len(summary_df)} domains) to {summary_tsv}")

    # --- Save final intolerant domain set ---
    intolerant_domains_pickle = os.path.join(output_dir, f'all_intolerant_domains.{assembly}.pkl')
    with open(intolerant_domains_pickle, 'wb') as f:
        pickle.dump(analysis_intolerant_domains, f)
    logger.info(f"Saved final intolerant domains (N={len(analysis_intolerant_domains)}) to {intolerant_domains_pickle}")

    # --- Save curated functional reference set (used to build the KS composite) ---
    curated_refs_pickle = os.path.join(output_dir, f'curated_functional_references.{assembly}.pkl')
    curated_refs = set(ref_domain_scores.keys())
    with open(curated_refs_pickle, 'wb') as f:
        pickle.dump(curated_refs, f)
    logger.info(f"Saved curated functional references (N={len(curated_refs)}) to {curated_refs_pickle}")

    # --- Save excluded domain sets per filter ---
    if excluded_sets:
        excluded_pickle = os.path.join(output_dir, f'excluded_domains.{assembly}.pkl')
        with open(excluded_pickle, 'wb') as f:
            pickle.dump(excluded_sets, f)
        logger.info(f"Saved excluded domain sets to {excluded_pickle}")
