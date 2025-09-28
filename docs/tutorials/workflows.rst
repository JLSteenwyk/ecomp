Workflows
=========

The following mini-workflows mirror the style of the PhyKIT tutorials.  They are
intended as copy/paste starting points; adjust paths to suit your dataset.

Compress an alignment directory
-------------------------------

.. code-block:: bash

    mkdir -p archives
    for fasta in data_large/ascomycota_data/*.fasta; do
        base=$(basename "$fasta" .fasta)
        tree="${fasta}.tre"
        tree_args=()
        [ -f "$tree" ] && tree_args=(--tree "$tree")
        ecomp compress "$fasta" "archives/${base}.ecomp" "${tree_args[@]}" \
            -m "archives/${base}.json"
    done

    # Quick validation pass
    for archive in archives/*.ecomp; do
        ecomp inspect "$archive" --summary
    done

Benchmark canonical codecs vs. eComp
------------------------------------

.. code-block:: bash

    python scripts/compare_compressors.py \
        data_large/ascomycota_data/EOG092D005G.fasta \
        data_large/ascomycota_data/EOG092D00LL.fasta \
        --codecs ecomp gzip bzip2 xz --runner internal \
        --output docs/ascomycota_small_benchmark.json

    # Review the JSON to ensure eComp beats gzip/bzip2; plot later if desired.

Tree-guided ordering sweep
--------------------------

.. code-block:: bash

    python scripts/benchmark_tree_ordering.py \
        --output docs/ascomycota_tree_ordering.json \
        data_large/ascomycota_data/*.fasta

    # Post-process in Python to compute average size deltas
    python - <<'PY'
    import json
    from pathlib import Path
    data = json.loads(Path('docs/ascomycota_tree_ordering.json').read_text())
    groups = {}
    for record in data:
        groups.setdefault(record['alignment'], {})[record['mode']] = record
    gains = [group['tree_provided']['compression_ratio'] - group['baseline']['compression_ratio']
             for group in groups.values() if 'tree_provided' in group]
    print('Average ratio delta (provided tree):', sum(gains)/len(gains))
    PY

Integrate with Python pipelines
-------------------------------

.. code-block:: python

    from evolutionary_compression import read_alignment, compress_alignment

    frame = read_alignment("alignment.fasta")
    with open("alignment.fasta.tre") as handle:
        frame.metadata["tree_newick"] = handle.read()
    result = compress_alignment(frame)
    print("payload bytes", len(result.payload))
    print("metadata", result.metadata)

Use the patterns above as building blocks for more complex workflows (e.g.,
concatenating results, integrating with Snakemake or Airflow).  When in doubt,
start with a small subset of alignments and confirm round-trip behaviour before
launching large jobs.
