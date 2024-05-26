from collections import Counter
import gzip
from typing import Callable
from typing import Iterator

import pandas as pd
from pysam import VariantFile
from pysam.libcbcf import VariantRecord, VariantRecordInfo as Info


biotype_ranks = {
    'protein_coding_CDS_not_defined': 1,
    'protein_coding': 2,
    'protein_coding_LoF': 3,
}

canonical_ranks = {
    False: 0,
    True: 1,
}

consequence_ranks = {
    'downstream_gene_variant': -1,
    'intergenic_variant': -1,
    'upstream_gene_variant': -1,
}

impact_ranks = {
    'MODIFIER': 0,
    'LOW': 1,
    'MODERATE': 2,
    'HIGH': 3,
}


def _sort_transcripts_key(transcript: dict):
    return (
        canonical_ranks.get(transcript['canonical'], 0),
        impact_ranks.get(transcript['impact'], 0),
        consequence_ranks.get(transcript['consequence'], 0),
        biotype_ranks.get(transcript['biotype'], 0),
    )


def explore_transcript_key(dataset: list[dict], key: str):
    print(f'---  {key} summary  ---')
    counter = Counter()
    for entry in dataset:
        for transcript in entry['annotations_per_transcript']:
            counter[transcript[key]] += 1
    for key, val in counter.most_common():
        print(f'''{key if key else "''"} --> {val}''')
    print()


def get_info_val(info: Info, key: str, default: float = 0.0, cast_func: Callable = float):
    return cast_func(info.get(key, [default])[0])


def _parse_transcript(transcript: dict, canonical_ids: set[str]):
    try:
        transcript_id, transcript_version = transcript['Feature'].split('.')
    except ValueError:
        transcript_id, transcript_version = transcript['Feature'], ''
    return {
        'amino_acid_change': transcript['Amino_acids'],
        'biotype': transcript['BIOTYPE'],
        'canonical': transcript_id in canonical_ids,
        'consequence': transcript['Consequence'],
        'impact': transcript['IMPACT'],
        'gene_name': transcript['SYMBOL'],
        'ncbi_id': transcript['Gene'],
        'strand': transcript['STRAND'],
        'transcript_id': transcript_id,
        'transcript_version': transcript_version,
    }


def add_transcripts_data(info: dict, variant_new: dict, keys: list[str], canonical_ids: set[str]):
    values_per_transcript = (entry.split('|') for entry in info['vep_output'])
    transcripts = (dict(zip(keys, values)) for values in values_per_transcript)
    transcripts = [_parse_transcript(t, canonical_ids) for t in transcripts]
    transcripts.sort(key=_sort_transcripts_key)
    variant_new.update(**transcripts[-1])


def _parse_variant(variant: VariantRecord, keys: list[str], canonical_ids: set[str]):
    ref, alt = variant.ref, variant.alts[0]
    info = dict(variant.info)
    variant_new = {
        'variant_id': f'{variant.chrom}_{variant.pos}_{ref}_{alt}',
        'ref': ref,
        'alt': alt,
        **{key: val for key, val in info.items() if key != 'vep_output'},
    }
    add_transcripts_data(info, variant_new, keys, canonical_ids)
    return variant_new


def _process_gff_annotation(annotation: str):
    fields = annotation.split(';')
    data = {}
    for field in fields:
        key, val = field.strip().split('=')
        data[key] = val
    return data


def _parse_refseq(gff_file: str):
    canonical_ids = set()
    with gzip.open(gff_file, 'rt') as fin:
        for row in fin:
            if row.startswith('#'):
                continue
            columns = row.split('\t')
            feature_type = columns[2]
            if feature_type != 'mRNA':
                continue
            annotation = columns[8]
            data = _process_gff_annotation(annotation)
            if 'tag' in data:  # the value of data['tag'] (if present) is always 'RefSeq Select'
                transcript_id = data['Name'].split('.')[0]
                canonical_ids.add(transcript_id)
    return canonical_ids


def parse_vcf(vcf_file: str, canonical_ids: set[str]) -> Iterator[dict]:
    with VariantFile(vcf_file, 'r') as reader:
        keys = reader.header.info['vep_output'].description.split(':', maxsplit=1)[1].strip().split('|')
        for variant in reader.fetch():
            yield _parse_variant(variant, keys, canonical_ids)


def main():
    gff_file = 'data/GCF_000001405.25_GRCh37.p13_genomic.gff.gz'
    canonical_ids = _parse_refseq(gff_file)
    vcf_file = 'data/clinvar_2020_from_clingen.vcf.gz'
    data = parse_vcf(vcf_file, canonical_ids)
    df = pd.DataFrame(data)
    df.to_csv('data/clinvar_annotated.tsv.gz', sep='\t', index=False, compression='gzip')


if __name__ == '__main__':
    main()
