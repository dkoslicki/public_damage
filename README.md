# Damage Bayes

Damage Bayes provides a production-ready command line interface for estimating
ancient DNA damage patterns using the Bayesian approach developed in the
`Damage-Bayes-v1.ipynb` notebook. The tool exposes two commands:

* ``damage-bayes count`` – generate 5′/3′ nucleotide frequency tables from
  FASTQ reads.
* ``damage-bayes analyze`` – run the PyMC model on paired frequency tables and
  produce a PDF report alongside a CSV summary of posterior means.

## Installation

The project follows a standard Python packaging layout and can be installed from
a source checkout:

```bash
python -m pip install .
```

Once published to PyPI or Conda, the same entry point will be available through
``pip install damage-bayes`` or ``conda install damage-bayes``.

## Usage

### Counting terminal bases

```bash
damage-bayes count data/subsample_1.fastq data/subsample_2.fastq \
  --output-dir freq_tables --n-bp 25 --trim 1
```

This command writes ``*_5_end_freq`` and ``*_3_end_freq`` tab-separated files to
``freq_tables``. The FASTQ parser accepts multi-line sequences and ignores any
non-ATCG characters when tabulating counts, reproducing the notebook behaviour
without brittle shell scripts.

### Running the Bayesian model

```bash
damage-bayes analyze --input-dir freq_tables \
  --output-pdf damage_summary.pdf \
  --output-csv damage_summary.csv \
  --model-start 2 --model-end 10 \
  --draws 1000 --tune 1000 --chains 4 --cores 1
```

The analysis command expects matching ``*_5_end_freq`` and ``*_3_end_freq`` files
for each sample. It fits the same PyMC model as the notebook, then generates a
smiley-plot PDF and a CSV table containing posterior means of ``alpha`` and
``beta`` for both C→T (5′) and G→A (3′) transitions.

See ``damage-bayes --help`` for the complete list of options.

## Testing

```bash
pytest
```

The test suite covers FASTQ parsing edge cases, deterministic nucleotide counts
from the provided test data, and smoke tests for the Bayesian model with reduced
sampling iterations.
