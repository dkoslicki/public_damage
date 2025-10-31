from damage_bayes.fastq import read_fastq


def test_multiline_sequences_are_parsed(tmp_path):
    fastq_path = tmp_path / "example.fastq"
    fastq_path.write_text(
        "@read1\n" "ACGT\n" "+\n" "!!!!\n" "@read2\n" "AAAA\nTT\n" "+\n" "!!!!!!\n"
    )

    records = list(read_fastq(str(fastq_path)))
    assert len(records) == 2
    assert records[0].sequence == "ACGT"
    assert records[1].sequence == "AAAATT"
    assert len(records[1].quality) == len(records[1].sequence)
