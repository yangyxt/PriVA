# GoFCards network dependencies: reliability, cost, and what fails if they are down

This records what the GoFCards build actually depends on over the network, how
those services behave when measured, and what is lost if they are unavailable.

It also corrects two earlier claims of mine that were wrong. Both were made from
single samples taken during a bad window, and both are superseded by the
measurements below.

Producer: `/paedyl01/disk1/yangyxt/PriVA/scripts/build_gofcards_exact_gof_cache.py`
Response cache: `<gofcards_workdir>/gofcards_summary_cache.jsonl`

---

## 1. There are two endpoints, and they behave very differently

I had been describing them as one thing. They are not.

| | Backend table endpoint | Summary endpoint |
| --- | --- | --- |
| Path | `/admin-api/backend/data/hg19/GainFunCards_SNV/geneSymbol/page` | `/admin-api/backend/data/hg19/variantLevel/summary` |
| Method | POST | GET |
| **Calls per build** | **2** (one for SNV, one for Indel) | **2,033** (one per allele) |
| Response size | ~2.4 MB | a few hundred bytes |
| Typical duration | 50–64 seconds | 4–5 seconds |
| Supplies | `AAChange_refGene`, the ANNOVAR protein change | `hg38_start` and `Accession` |

Almost every failure I reported was the **table** endpoint — the one that
transfers 2.4 MB and takes about a minute. The endpoint that needs 2,033 calls
is the small, fast one, and it is considerably better behaved than I said.

---

## 2. My measurements were contaminated by my own sandbox

The environment this work was done in sets a proxy:

```
HTTPS_PROXY=http://127.0.0.1:18081
HTTP_PROXY=http://127.0.0.1:18081
```

Both `curl` and Python's `requests` honour those variables by default. Every
timing and failure I reported earlier went through that proxy, and I did not
account for it.

Measuring the same 20 alleles both ways:

```
DIRECT (proxy bypassed)        20/20 success    median 5.0s    min 2.4s    max 6.4s
THROUGH THE SANDBOX PROXY      10/10 success    median 4.3s
```

**No failures in either configuration.** The proxy does not degrade the small
summary calls at all.

What the proxy *did* break was the 2.4 MB table request: it returned
`502 Bad Gateway` on a request that succeeded in 50 seconds when sent directly.
That is a gateway timeout on a long transfer, and it is a property of this
sandbox, not of GoFCards and not of the pipeline. It is why the build had to be
launched with the proxy variables unset.

---

## 3. Full measurement history

| When | Endpoint | Route | Observation |
| --- | --- | --- | --- |
| 2026-07-27 | table | proxied | completely unreachable, `http=000` after 62 s |
| 2026-07-28 early | summary | proxied via `curl` | 12–14 s per call |
| 2026-07-28 early | summary | proxied via `requests` | 3 of 5 calls returned non-JSON |
| 2026-07-28 early | table | proxied | `502 Bad Gateway`; `200` in 50 s direct |
| 2026-07-28 later | summary | **direct** | **20/20 success, median 5.0 s** |
| 2026-07-28 later | summary | **proxied** | **10/10 success, median 4.3 s** |

The service does have genuinely bad windows. Yesterday's total outage was real,
and this morning's non-JSON responses were real. But it is not persistently
broken, and its steady state is roughly five seconds per call.

### The two claims I got wrong

- **"13 seconds per call."** Measured once, through the proxy, during a slow
  window. The measured median is 5.0 seconds direct and 4.3 seconds proxied.
- **"Returned non-JSON for three of five calls."** True at that moment, and
  false as a description of the service. Twenty consecutive calls succeeded
  afterwards.

Both statements described one bad window and presented it as the norm.

---

## 4. Corrected cost of a first build

```
2,033 calls x 5.0 s median = 2.8 hours
```

My earlier estimates were "roughly an hour" and then "7.3 hours". Both were from
single samples. The first was too optimistic, the second too pessimistic.

This cost is paid **only on a first build**, or after GoFCards publishes a new
release. Once the response cache exists, subsequent builds make zero summary
calls.

---

## 5. What already protects the build

### Per-call failures are handled

Each call gets four attempts with 5, 10 and 15 second backoff:

```python
for attempt in range(1, 5):
    try:
        reply = session.get(SUMMARY_ENDPOINT, headers=API_HEADERS,
                            params=params, timeout=timeout_seconds)
        reply.raise_for_status()
        data = (reply.json() or {}).get("data") or {}
        break
    except Exception as exc:
        if attempt == 4:
            raise RuntimeError(...)
        time.sleep(5 * attempt)
```

Note that `reply.json()` sits **inside** the `try`. A response that returns HTTP
200 but is not JSON — exactly the failure seen this morning — raises and is
retried, rather than silently producing an empty record. That matters: the
alternative would be a variant quietly losing its ClinVar accession with no
error anywhere.

### Progress is never lost

Every response is appended to `gofcards_summary_cache.jsonl` and flushed
immediately, one line per allele:

```json
{"allele_key": "SNV|9|5073770|5073770|G|T", "data": {"hg38_start": 5073770, "Accession": "14662"}}
```

An outage two hours into a first build keeps those two hours of work. A rerun
re-fetches only what is missing. In the most recent build the cache already held
2,033 entries and the fetcher requested exactly the 51 that were absent:

```
[gofcards] summary cache already holds 2033 alleles
[gofcards] summary endpoint: 2084 alleles; 2084 carry a GRCh38 position; 1250 carry a ClinVar accession
```

### The risk that remains

Four retries spanning 30 seconds does not survive a multi-hour outage of the
kind seen on 2026-07-27. In that case the build stops and must be rerun later.
The cache means it resumes rather than restarts.

---

## 6. What is actually lost if the summary endpoint is unavailable

Exactly two things. Everything else in the cache comes from the public workbook,
the two table calls, and local reference data.

### The ClinVar accession for 1,245 variants

This degrades ClinVar linking from two independent confirmations to one. It does
not break it. In the most recent build:

```
variants_annotated                  = 525
confirmed_by_gofcards_variation_id  = 520
confirmed_by_genomic_coordinates    = 525
confirmed_by_both_routes            = 520
```

Coordinate matching found **every** link on its own, and where GoFCards also
supplied an identifier the two agreed **520 out of 520, with zero
disagreements**. Losing the accession costs the second opinion, not the link.

### The GRCh38 position for variants the liftover cannot place

```
GRCh38 coordinates: 2027 from the liftover, 1 recovered from the GoFCards endpoint
```

**One variant** in this build — CRLF2 `loc_X:1314966:A->C_grch37`, which sits in
the pseudoautosomal region at the tip of the X chromosome, where the boundaries
moved between GRCh37 and GRCh38 so the chain file cannot place it.

### What does not depend on it at all

The GRCh37 records, both VEP runs, the liftover, every protein and coding key,
the gene provenance check, the reviewed-mechanism gate, and all 3,154 evidence
entries. None of these touch the summary endpoint.

---

## 7. Summary

This is a slow dependency with occasional multi-hour outages. It costs about
three hours on a first build, is fully cached afterwards, retries transient
failures correctly, and never loses completed work.

If it were unavailable permanently, the cache would lose one GRCh38 coordinate
and a second opinion on ClinVar links that already agree with the first.

That is a real fragility worth knowing about, and a smaller one than I first
described.
