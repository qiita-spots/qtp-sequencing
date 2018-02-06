# -----------------------------------------------------------------------------
# Copyright (c) 2014--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------

from qiita_client import QiitaTypePlugin, QiitaArtifactType

from .validate import validate
from .summary import generate_html_summary

# DEfine the supported artifact types
artifact_types = [
    QiitaArtifactType('SFF', 'Raw SFF files', False, False, True,
                      [('raw_sff', True)]),
    QiitaArtifactType('FASTA_Sanger', 'Raw fasta files from Sanger sequencing',
                      False, False, False, [('raw_fasta', True)]),
    QiitaArtifactType('FASTQ', 'Raw fastq files, with or without paired ends',
                      False, False, True, [('raw_forward_seqs', True),
                                           ('raw_reverse_seqs', False),
                                           ('raw_barcodes', True)]),
    QiitaArtifactType('FASTA', 'Raw fasta files', False, False, False,
                      [('raw_fasta', True), ('raw_qual', False)]),
    QiitaArtifactType('per_sample_FASTQ', 'Raw per sample FASTQ files', False,
                      False, True, [('raw_forward_seqs', True),
                                    ('raw_reverse_seqs', False)]),
    QiitaArtifactType('Demultiplexed',
                      'Demultiplexed and QC sequeunces with QIIME labels',
                      True, True, False, [('preprocessed_fasta', True),
                                          ('preprocessed_fastq', True),
                                          ('preprocessed_demux', False),
                                          ('log', False)]),
]

# Initialize the plugin
plugin = QiitaTypePlugin('Target Gene type', '0.1.0',
                         'Target gene artifact types plugin',
                         validate, generate_html_summary, artifact_types)
