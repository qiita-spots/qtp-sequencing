# -----------------------------------------------------------------------------
# Copyright (c) 2014--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------

from unittest import main
from tempfile import mkdtemp
from os import remove
from os.path import exists, isdir, join, dirname
from shutil import rmtree, copyfile
from json import dumps

from qiita_client.testing import PluginTestCase
from gzip import GzipFile

from qtp_target_gene.fastqc import (generate_fastqc_commands,
                   generate_multiqc_commands, _guess_fastqc_filename, fastqc)


class SummaryTestsNotDemux(PluginTestCase):
    def setUp(self):
        self.out_dir = mkdtemp()
        self._clean_up_files = [self.out_dir]

    def tearDown(self):
        for fp in self._clean_up_files:
            if exists(fp):
                if isdir(fp):
                    rmtree(fp)
                else:
                    remove(fp)

    # MultiQC tests
    def test_generate_multiqc_commands(self):

        exp_cmd = ['multiqc infile.txt --outdir outdir '
                   '--filename MultiQC.html file_dir']

        obs_cmd = generate_multiqc_commands('file_dir', 'out_dir')

        self.assertEqual(obs_cmd, exp_cmd)

    def test_generate_fastqc_commands_fwd_rev(self):
        fd, fp = mkstemp()
        close(fd)
        with open(fp, 'w') as f:
            f.write(MAPPING_FILE)
        self._clean_up_files.append(fp)

        exp = ['mkdir -p output/s1; fastqc --outdir "output/s1" '
               '--noextract fastq/s1.fastq fastq/s1.R2.fastq',
               'mkdir -p output/s2; fastqc --outdir "output/s2" '
               '--noextract fastq/s2.fastq.gz '
               'fastq/s2.R2.fastq.gz',
               'mkdir -p output/s3; fastqc --outdir "output/s3" '
               '--noextract fastq/s3.fastq fastq/s3.R2.fastq']

        obs, samp = generate_fastqc_commands(
            ['fastq/s1.fastq', 'fastq/s2.fastq.gz', 'fastq/s3.fastq'],
            ['fastq/s1.R2.fastq', 'fastq/s2.R2.fastq.gz', 'fastq/s3.R2.fastq'],
            fp, 'output')

        self.assertEqual(obs, exp)

    def test_generate_fastqc_commands_fwd(self):
        fd, fp = mkstemp()
        close(fd)
        with open(fp, 'w') as f:
            f.write(MAPPING_FILE)
        self._clean_up_files.append(fp)

        exp = ['mkdir -p output/s1; fastqc --outdir "output/s1" '
               '--noextract fastq/s1.fastq',
               'mkdir -p output/s2; fastqc --outdir "output/s2" '
               '--noextract fastq/s2.fastq.gz',
               'mkdir -p output/s3; fastqc --outdir "output/s3" '
               '--noextract fastq/s3.fastq']

        obs, samp = generate_fastqc_commands(
                ['fastq/s1.fastq', 'fastq/s2.fastq.gz', 'fastq/s3.fastq'], [],
                fp, 'output')

        self.assertEqual(obs, exp)

    def test_fastqc(self):
        # map file
        fd, fp = mkstemp()
        close(fd)
        with open(fp, 'w') as f:
            f.write(MAPPING_FILE)
        self._clean_up_files.append(fp)

        # generating filepaths
        in_dir = mkdtemp()
        self._clean_up_files.append(in_dir)

        fp1_1 = join(in_dir, 's1_R1.fastq.gz')
        fp1_2 = join(in_dir, 's1_R2.fastq.gz')
        copyfile('support_files/kd_test_1_R1.fastq.gz', fp1_1)
        copyfile('support_files/kd_test_1_R2.fastq.gz', fp1_2)

        jid = self.qclient.post('/apitest/processing_job/', data=data)['job']

        out_dir = mkdtemp()
        self._clean_up_files.append(out_dir)

        success, fps, msg = fastqc(self.qclient, jid,
                                     fp, out_dir)

        # Check all output from fastqc
        self.assertEqual("", msg)
        self.assertTrue(success)
        self.assertEqual(2, len(fps))

        # Verify files exist and have correct names
        exp_fps = [[(join(out_dir, 'multiqc.tar.gz'),
                    'tgz')],
                   [(join(out_dir, 'fastqc.tar.gz'),
                    'tgz')]]

        self.assertItemsEqual(exp_fps, obs_fps)

        for f_a in exp_fps:
            assert exists(f_a[0][0])


MAPPING_FILE = (
    "#SampleID\tplatform\tbarcode\texperiment_design_description\t"
    "library_construction_protocol\tcenter_name\tprimer\trun_prefix\t"
    "instrument_model\tDescription\n"
    "SKB7.640196\tILLUMINA\tA\tA\tA\tANL\tA\ts3\tIllumina MiSeq\tdesc1\n"
    "SKB8.640193\tILLUMINA\tA\tA\tA\tANL\tA\ts1\tIllumina MiSeq\tdesc2\n"
    "SKD8.640184\tILLUMINA\tA\tA\tA\tANL\tA\ts2\tIllumina MiSeq\tdesc3\n"
)

