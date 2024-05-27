from collections import Counter
import gzip
import os
from typing import Callable
from typing import Iterator

import pandas as pd
import requests
from pysam import VariantFile
from pysam.libcbcf import VariantRecord, VariantRecordInfo as Info
import tqdm


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


def update_gnomad_data(variant_new: dict):
    """As these fields represent counts and frequencies, they can be safely replaced them with zeros when missing."""
    prefix = 'gnomad_'
    for key in ('ac', 'af', 'an', 'hom'):
        full_key = prefix + key
        if full_key not in variant_new:
            variant_new[full_key] = 0


def update_spliceai_data(variant_new: dict):
    max_score = None
    prefix = 'spliceai_'
    for key in ('ag', 'al', 'dg', 'dl'):
        full_key = prefix + key
        if (score := variant_new.get(full_key)) is None:
            continue
        del variant_new[full_key]
        if max_score is None:
            max_score = score
        else:
            max_score = max(max_score, score)
    if max_score is not None:
        variant_new['spliceai'] = max_score


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


def add_transcripts_data(info: Info, variant_new: dict, keys: list[str], canonical_ids: set[str]):
    values_per_transcript = (entry.split('|') for entry in info['vep_output'])
    transcripts = (dict(zip(keys, values)) for values in values_per_transcript)
    transcripts = [_parse_transcript(t, canonical_ids) for t in transcripts]
    transcripts.sort(key=_sort_transcripts_key)
    variant_new.update(**transcripts[-1])


def _parse_variant(variant: VariantRecord, keys: list[str], canonical_ids: set[str]) -> dict:
    ref, alt = variant.ref, variant.alts[0]
    info = variant.info
    variant_new = {
        'variant_id': f'{variant.chrom}_{variant.pos}_{ref}_{alt}',
        'ref': ref,
        'alt': alt,
        **{key: val for key, val in info.items() if key != 'vep_output'},
    }
    add_transcripts_data(info, variant_new, keys, canonical_ids)
    update_gnomad_data(variant_new)
    update_spliceai_data(variant_new)
    return variant_new


def parse_vcf(vcf_file: str, canonical_ids: set[str]) -> Iterator[dict]:
    with VariantFile(vcf_file, 'r') as reader:
        keys = reader.header.info['vep_output'].description.split(':', maxsplit=1)[1].strip().split('|')
        for variant in reader.fetch():
            yield _parse_variant(variant, keys, canonical_ids)


def _process_gff_annotation(annotation: str) -> dict[str, str]:
    fields = annotation.split(';')
    data = {}
    for field in fields:
        key, val = field.strip().split('=')
        data[key] = val
    return data


def _get_canonical_transcript_ids(gff_file: str):
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


def _parse_gene_metrics(metrics_file: str):
    with gzip.open(metrics_file, 'rt') as fin:
        keys = next(fin).strip('\n').split('\t')
        for line in fin:
            values = line.strip('\n').split('\t')
            data = dict(zip(keys, values))
            pass


def _get_remote_file_size(url: str) -> int:
    response = requests.head(url)
    response.raise_for_status()
    file_size = response.headers.get('Content-Length')
    return int(file_size)


def download_file(url: str, output_file: str, *, chunk_size: int = 65536, compress: bool = False, force: bool = False):
    if compress and not output_file.endswith('.gz'):
        output_file = output_file + '.gz'
    if not force and os.path.exists(output_file):
        return output_file
    file_size = _get_remote_file_size(url)
    response = requests.get(url, stream=True)
    response.raise_for_status()
    print(f'Downloading {url.split("/")[-1]} into {output_file}', flush=True)
    if not compress:
        _open = open
    else:
        _open = gzip.open
    progress_bar = tqdm.tqdm(total=file_size, unit='B', unit_scale=True, unit_divisor=1024, miniters=1)
    with _open(output_file, 'wb') as fout, progress_bar:
        for chunk in response.iter_content(chunk_size):
            if chunk:
                fout.write(chunk)
                progress_bar.update(chunk_size)


def download_refseq(force: bool = False):
    url = 'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.25_GRCh37.p13/GCF_000001405.25_GRCh37.p13_genomic.gff.gz'
    output_file = 'output/refseq.gff.gz'
    return download_file(url, output_file, force=force)


def download_gene_metrics(force: bool = False):
    url = 'https://gnomad-public-us-east-1.s3.amazonaws.com/release/4.1/constraint/gnomad.v4.1.constraint_metrics.tsv'
    output_file = 'output/gene_metrics.tsv'
    return download_file(url, output_file, compress=True, force=force)


def main():
    gff_file = 'output/refseq.gff.gz'
    canonical_ids = _get_canonical_transcript_ids(gff_file)
    vcf_file = 'output/clinvar_2020_from_clingen.vcf.gz'
    data = parse_vcf(vcf_file, canonical_ids)
    df = pd.DataFrame(data)
    df.to_csv('output/clinvar_annotated.tsv.gz', sep='\t', index=False, compression='gzip')


if __name__ == '__main__':
    metrics_file = download_gene_metrics()
    _parse_gene_metrics(metrics_file)
    # main()
    # from time import perf_counter_ns
    # start = perf_counter_ns()
    # download_refseq(force=True)
    # end = perf_counter_ns()
    # print(f'{(end - start) / 10**9}:.4f')
