Workflows
=========

The following mini-workflows mirror the style of the PhyKIT tutorials. They are
intended as copy/paste starting points; adjust paths to match your dataset.

Compress an alignment directory
-------------------------------

.. code-block:: bash

    SUPP=${EVOCOMP_DATA_ROOT:-../EVOCOMP_MANUSCRIPT}
    mkdir -p archives
    for fasta in "$SUPP/data_large/ascomycota_data"/*.fasta; do
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

    SUPP=${EVOCOMP_DATA_ROOT:-../EVOCOMP_MANUSCRIPT}
    target="$SUPP/data_large/ascomycota_data/EOG092D005G.fasta"
    /usr/bin/time -p ecomp compress "$target" --output "$target.ecomp"
    /usr/bin/time -p gzip -k "$target"
    /usr/bin/time -p bzip2 -k "$target"
    /usr/bin/time -p xz -k "$target"

    jq '.compression_ratio' "$target.ecomp.json"

    # Clean up auxiliary archives
    rm -f "$target".{ecomp,gz,bz2,xz,json}

The CLI emits JSON metadata (when ``-m`` or ``--meta`` is provided) that you can
post-process with Python, ``jq``, or pandas.

Tree-guided ordering sweep
--------------------------

.. code-block:: bash

    SUPP=${EVOCOMP_DATA_ROOT:-../EVOCOMP_MANUSCRIPT}
    target="$SUPP/data_large/ascomycota_data/EOG092D005G.fasta"
    tree="${target}.tre"
    for mode in baseline greedy mst; do
        ECOMP_SEQUENCE_ORDER="$mode" ecomp compress "$target" "${target}.${mode}.ecomp" \
            --tree "$tree" -m "${target}.${mode}.json"
    done

    python - <<'PY'
    import json
    from pathlib import Path
    import os

    supp = Path(os.environ.get("EVOCOMP_DATA_ROOT", "../EVOCOMP_MANUSCRIPT"))
    base = supp / "data_large/ascomycota_data/EOG092D005G.fasta"
    ratios = {}
    for mode in ("baseline", "greedy", "mst"):
        meta = json.loads(Path(f"{base}.{mode}.json").read_text())
        ratios[mode] = meta["compression_ratio"]
    print("Compression ratios by mode:")
    for mode, ratio in ratios.items():
        print(f"  {mode}: {ratio:.3f}")
    PY

    rm -f "$target".{baseline,greedy,mst}.ecomp
    rm -f "$target".{baseline,greedy,mst}.json

Integrate with Python pipelines
-------------------------------

.. code-block:: python

    from evolutionary_compression import read_alignment, compress_alignment

    frame = read_alignment("data/fixtures/small_phylo.fasta")
    with open("data/fixtures/small_phylo.tree") as handle:
        frame.metadata["tree_newick"] = handle.read()
    result = compress_alignment(frame)
    print("payload bytes", len(result.payload))
    print("metadata", result.metadata)

Use the patterns above as building blocks for more intricate workflows (e.g.,
concatenating results, integrating with Snakemake or Airflow). When in doubt,
start with a small subset of alignments and confirm round-trip behaviour before
launching large jobs.
