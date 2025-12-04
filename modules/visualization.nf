process VISUALIZATION {
    label 'process_low'
    publishDir "${params.outdir}/visualization", mode: 'copy'

    input:
    path(summary_csv)

    output:
    path("plots/*.png"), emit: plots

    script:
    """
    pip install seaborn pandas matplotlib > /dev/null
    echo "PATH: \$PATH"
    which visualize_results.py || echo "visualize_results.py not found in PATH"
    
    # Try running with python3 explicitly if found, or assume it's in PATH
    if which visualize_results.py > /dev/null; then
        python3 \$(which visualize_results.py) ${summary_csv} plots
    else
        echo "Script not found, checking current dir"
        ls -l
        # Fallback if bin is not in PATH (should not happen in Nextflow)
        # But we can't easily access projectDir here without passing it.
        # Let's assume it is in PATH or fail with clear message.
        exit 1
    fi
    """
}
