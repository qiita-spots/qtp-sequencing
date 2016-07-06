# -----------------------------------------------------------------------------
# Copyright (c) 2014--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------

from unittest import TestCase, main
from tempfile import mkdtemp, mkstemp
from os import remove, close, environ
from os.path import exists, isdir, basename, splitext, join
from shutil import rmtree
from json import dumps

from qiita_client import QiitaClient, ArtifactInfo
from h5py import File
from qiita_ware.demux import to_hdf5

from qtp_target_gene.validate import (
    _validate_multiple, _validate_per_sample_FASTQ, _validate_demux_file,
    _validate_demultiplexed, validate)


CLIENT_ID = '19ndkO3oMKsoChjVVWluF7QkxHRfYhTKSFbAVt8IhK7gZgDaO4'
CLIENT_SECRET = ('J7FfQ7CQdOxuKhQAf1eoGgBAE81Ns8Gu3EKaWFm3IO2JKh'
                 'AmmCWZuabe0O5Mp28s1')


class ValidateTests(TestCase):
    @classmethod
    def setUpClass(cls):
        server_cert = environ.get('QIITA_SERVER_CERT', None)
        cls.qclient = QiitaClient("https://localhost:21174", CLIENT_ID,
                                  CLIENT_SECRET, server_cert=server_cert)

    @classmethod
    def tearDownClass(cls):
        cls.qclient.post("/apitest/reset/")

    def setUp(self):
        self._clean_up_files = []

    def tearDown(self):
        for fp in self._clean_up_files:
            if exists(fp):
                if isdir(fp):
                    rmtree(fp)
                else:
                    remove(fp)

    def _create_template_and_job(self, prep_info, files, atype):
        data = {'prep_info': dumps(prep_info),
                'study': 1,
                'data_type': '16S'}
        template = self.qclient.post(
            '/apitest/prep_template/', data=data)['prep']

        parameters = {'template': template,
                      'files': dumps(files),
                      'artifact_type': atype}

        data = {'command': 6,
                'parameters': dumps(parameters),
                'status': 'running'}
        job_id = self.qclient.post(
            '/apitest/processing_job/', data=data)['job']

        return job_id

    def test_validate_multiple(self):
        prep_info = {
            '1.SKB2.640194': {'run_prefix': 'prefix1'},
            '1.SKM4.640180': {'run_prefix': 'prefix1'},
            '1.SKB3.640195': {'run_prefix': 'prefix2'}}
        files = {'raw_forward_seqs': ['/path/to/prefix1.fastq',
                                      '/path/to/prefix2.fastq'],
                 'raw_barcodes': ['/path/to/prefix1_b.fastq',
                                  '/path/to/prefix2_b.fastq']}
        atype = "FASTQ"
        job_id = self._create_template_and_job(prep_info, files, atype)

        obs_success, obs_ainfo, obs_error = _validate_multiple(
            self.qclient, job_id, prep_info, files, atype)

        self.assertTrue(obs_success)
        filepaths = [
            ('/path/to/prefix1_b.fastq', 'raw_barcodes'),
            ('/path/to/prefix2_b.fastq', 'raw_barcodes'),
            ('/path/to/prefix1.fastq', 'raw_forward_seqs'),
            ('/path/to/prefix2.fastq', 'raw_forward_seqs')]
        exp = [ArtifactInfo(None, "FASTQ", filepaths)]
        self.assertEqual(obs_ainfo, exp)
        self.assertEqual(obs_error, "")

    def test_validate_multiple_single_lane(self):
        prep_info = {"1.SKB2.640194": {"not_a_run_prefix": "prefix1"},
                     "1.SKM4.640180": {"not_a_run_prefix": "prefix1"},
                     "1.SKB3.640195": {"not_a_run_prefix": "prefix2"}}
        files = {'raw_forward_seqs': ['/path/to/prefix1.fastq'],
                 'raw_barcodes': ['/path/to/prefix1_b.fastq']}
        atype = "FASTQ"
        job_id = self._create_template_and_job(prep_info, files, atype)

        obs_success, obs_ainfo, obs_error = _validate_multiple(
            self.qclient, job_id, prep_info, files, atype)

        self.assertTrue(obs_success)
        filepaths = [
            ('/path/to/prefix1_b.fastq', 'raw_barcodes'),
            ('/path/to/prefix1.fastq', 'raw_forward_seqs')]
        exp = [ArtifactInfo(None, atype, filepaths)]
        self.assertEqual(obs_ainfo, exp)
        self.assertEqual(obs_error, "")

    def test_validate_multiple_error(self):
        # Filepath type not supported
        prep_info = {"1.SKB2.640194": {"run_prefix": "prefix1"},
                     "1.SKM4.640180": {"run_prefix": "prefix1"},
                     "1.SKB3.640195": {"run_prefix": "prefix2"}}
        files = {'Unknown': ['/path/to/file1.fastq']}
        atype = "FASTQ"
        job_id = self._create_template_and_job(prep_info, files, atype)
        obs_success, obs_ainfo, obs_error = _validate_multiple(
            self.qclient, job_id, prep_info, files, atype)
        self.assertFalse(obs_success)
        self.assertIsNone(obs_ainfo)
        self.assertEqual(obs_error,
                         "Filepath type(s) Unknown not supported by artifact "
                         "type FASTQ. Supported filepath types: raw_barcodes, "
                         "raw_forward_seqs, raw_reverse_seqs")

        # Number of provided files != Num run prefix values
        files = {'raw_forward_seqs': ['/path/to/file1.fastq'],
                 'raw_barcodes': ['/path/to/file1_b.fastq',
                                  '/path/to/file2_b.fastq',
                                  '/path/to/file3_b.fastq']}
        job_id = self._create_template_and_job(prep_info, files, atype)
        obs_success, obs_ainfo, obs_error = _validate_multiple(
            self.qclient, job_id, prep_info, files, atype)
        error = ("Error creating artifact. Offending files:\nraw_forward_seqs:"
                 " The number of provided files (1) doesn't match the number "
                 "of run prefix values in the prep info (2): file1.fastq\n"
                 "raw_barcodes: The number of provided files (3) doesn't "
                 "match the number of run prefix values in the prep info (2): "
                 "file1_b.fastq, file2_b.fastq, file3_b.fastq")
        self.assertFalse(obs_success)
        self.assertIsNone(obs_ainfo)
        self.assertItemsEqual(obs_error.split('\n'), error.split('\n'))

        # File doesn't match any run prefix
        files = {'raw_forward_seqs': ['/path/to/file1.fastq',
                                      '/path/to/prefix2.fastq'],
                 'raw_barcodes': ['/path/to/file1_b.fastq',
                                  '/path/to/prefix2_b.fastq']}
        job_id = self._create_template_and_job(prep_info, files, atype)
        obs_success, obs_ainfo, obs_error = _validate_multiple(
            self.qclient, job_id, prep_info, files, atype)
        error = ("Error creating artifact. Offending files:\nraw_forward_seqs:"
                 " The provided files do not match the run prefix values in "
                 "the prep information: file1.fastq\n"
                 "raw_barcodes: The provided files do not match the run "
                 "prefix values in the prep information: file1_b.fastq")
        self.assertFalse(obs_success)
        self.assertIsNone(obs_ainfo)
        self.assertItemsEqual(obs_error.split('\n'), error.split('\n'))

        # A required filepath type is missing
        files = {'raw_forward_seqs': ['/path/to/prefix1.fastq',
                                      '/path/to/prefix2.fastq'],
                 'raw_reverse_seqs': ['/path/to/prefix1_rev.fastq',
                                      '/path/to/prefix2_rev.fastq']}
        job_id = self._create_template_and_job(prep_info, files, atype)
        obs_success, obs_ainfo, obs_error = _validate_multiple(
            self.qclient, job_id, prep_info, files, atype)
        self.assertFalse(obs_success)
        self.assertIsNone(obs_ainfo)
        self.assertEqual(obs_error,
                         "Missing required filepath type(s): raw_barcodes")

        # No run prefix and more than 1 lane
        prep_info = {"1.SKB2.640194": {"not_a_run_prefix": "prefix1"},
                     "1.SKM4.640180": {"not_a_run_prefix": "prefix1"},
                     "1.SKB3.640195": {"not_a_run_prefix": "prefix2"}}
        files = {'raw_forward_seqs': ['/path/to/prefix1.fastq',
                                      '/path/to/prefix2.fastq'],
                 'raw_barcodes': ['/path/to/prefix1_b.fastq',
                                  '/path/to/prefix2_b.fastq']}
        job_id = self._create_template_and_job(prep_info, files, atype)
        obs_success, obs_ainfo, obs_error = _validate_multiple(
            self.qclient, job_id, prep_info, files, atype)
        error = ("Error creating artifact. Offending files:\nraw_forward_seqs:"
                 " Only one file per type is allowed. Please provide the "
                 "column 'run_prefix' if you need more than one file per "
                 "type: prefix1.fastq, prefix2.fastq\n"
                 "raw_barcodes: Only one file per type is allowed. Please "
                 "provide the column 'run_prefix' if you need more than one "
                 "file per type: prefix1_b.fastq, prefix2_b.fastq")
        self.assertFalse(obs_success)
        self.assertIsNone(obs_ainfo)
        self.assertItemsEqual(obs_error.split('\n'), error.split('\n'))

    def test_validate_per_sample_FASTQ_run_prefix(self):
        prep_info = {"1.SKB2.640194": {"run_prefix": "prefix1"},
                     "1.SKM4.640180": {"run_prefix": "prefix2"},
                     "1.SKB3.640195": {"run_prefix": "prefix3"}}
        files = {'raw_forward_seqs': ['/path/to/prefix1_file.fastq',
                                      '/path/to/prefix2_file.fastq',
                                      '/path/to/prefix3_file.fastq']}
        job_id = self._create_template_and_job(
            prep_info, files, "per_sample_FASTQ")
        obs_success, obs_ainfo, obs_error = _validate_per_sample_FASTQ(
            self.qclient, job_id, prep_info, files)

        self.assertTrue(obs_success)
        filepaths = [('/path/to/prefix1_file.fastq', 'raw_forward_seqs'),
                     ('/path/to/prefix2_file.fastq', 'raw_forward_seqs'),
                     ('/path/to/prefix3_file.fastq', 'raw_forward_seqs')]
        exp = [ArtifactInfo(None, "per_sample_FASTQ", filepaths)]
        self.assertEqual(obs_ainfo, exp)
        self.assertEqual(obs_error, "")

    def test_validate_per_sample_FASTQ(self):
        prep_info = {"1.SKB2.640194": {"not_a_run_prefix": "prefix1"},
                     "1.SKM4.640180": {"not_a_run_prefix": "prefix1"},
                     "1.SKB3.640195": {"not_a_run_prefix": "prefix2"}}
        files = {'raw_forward_seqs': ['/path/to/SKB2.640194_file.fastq',
                                      '/path/to/SKM4.640180_file.fastq',
                                      '/path/to/SKB3.640195_file.fastq']}
        job_id = self._create_template_and_job(
            prep_info, files, "per_sample_FASTQ")
        obs_success, obs_ainfo, obs_error = _validate_per_sample_FASTQ(
            self.qclient, job_id, prep_info, files)
        self.assertTrue(obs_success)
        filepaths = [('/path/to/SKB2.640194_file.fastq', 'raw_forward_seqs'),
                     ('/path/to/SKM4.640180_file.fastq', 'raw_forward_seqs'),
                     ('/path/to/SKB3.640195_file.fastq', 'raw_forward_seqs')]
        exp = [ArtifactInfo(None, "per_sample_FASTQ", filepaths)]
        self.assertEqual(obs_ainfo, exp)
        self.assertEqual(obs_error, "")

    def test_validate_per_sample_FASTQ_error(self):
        # Filepath type not supported
        prep_info = {"1.SKB2.640194": {"run_prefix": "prefix1"},
                     "1.SKM4.640180": {"run_prefix": "prefix2"},
                     "1.SKB3.640195": {"run_prefix": "prefix3"}}
        files = {'Unknown': ['/path/to/file1.fastq']}
        job_id = self._create_template_and_job(
            prep_info, files, "per_sample_FASTQ")
        obs_success, obs_ainfo, obs_error = _validate_per_sample_FASTQ(
            self.qclient, job_id, prep_info, files)
        self.assertFalse(obs_success)
        self.assertIsNone(obs_ainfo)
        self.assertEqual(obs_error,
                         "Filepath type(s) Unknown not supported by artifact "
                         "type per_sample_FASTQ. Supported filepath types: "
                         "raw_forward_seqs, raw_reverse_seqs")

        # Missing raw_forward_seqs
        files = {'raw_reverse_seqs': ['/path/to/file1.fastq']}
        job_id = self._create_template_and_job(
            prep_info, files, "per_sample_FASTQ")
        obs_success, obs_ainfo, obs_error = _validate_per_sample_FASTQ(
            self.qclient, job_id, prep_info, files)
        self.assertFalse(obs_success)
        self.assertIsNone(obs_ainfo)
        self.assertEqual(obs_error,
                         "Missing required filepath type: raw_forward_seqs")

        # Count mismatch
        files = {'raw_forward_seqs': ['/path/to/file1.fastq'],
                 'raw_reverse_seqs': ['/path/to/file1.fastq',
                                      '/path/to/file1.fastq']}
        job_id = self._create_template_and_job(
            prep_info, files, "per_sample_FASTQ")
        obs_success, obs_ainfo, obs_error = _validate_per_sample_FASTQ(
            self.qclient, job_id, prep_info, files)
        self.assertFalse(obs_success)
        self.assertIsNone(obs_ainfo)
        self.assertEqual(obs_error,
                         "The number of provided files doesn't match the "
                         "number of samples (3): 1 raw_forward_seqs, "
                         "2 raw_reverse_seqs (optional, 0 is ok)")

        # Run prefix mismatch
        files = {'raw_forward_seqs': ['/path/to/prefix1_fwd.fastq',
                                      '/path/to/prefix2_fwd.fastq',
                                      '/path/to/Aprefix3_fwd.fastq']}
        job_id = self._create_template_and_job(
            prep_info, files, "per_sample_FASTQ")
        obs_success, obs_ainfo, obs_error = _validate_per_sample_FASTQ(
            self.qclient, job_id, prep_info, files)
        self.assertFalse(obs_success)
        self.assertIsNone(obs_ainfo)
        self.assertEqual(obs_error,
                         "The provided files do not match the run prefix "
                         "values in the prep information. Offending files: "
                         "raw_forward_seqs: Aprefix3_fwd.fastq, "
                         "raw_reverse_seqs: ")

        # Non-unique run-prefix values
        prep_info = {"1.SKB2.640194": {"run_prefix": "prefix1"},
                     "1.SKM4.640180": {"run_prefix": "prefix1"},
                     "1.SKB3.640195": {"run_prefix": "prefix3"}}
        job_id = self._create_template_and_job(
            prep_info, files, "per_sample_FASTQ")
        obs_success, obs_ainfo, obs_error = _validate_per_sample_FASTQ(
            self.qclient, job_id, prep_info, files)
        self.assertFalse(obs_success)
        self.assertIsNone(obs_ainfo)
        self.assertEqual(obs_error,
                         "The values for the column 'run_prefix' are not "
                         "unique for each sample. Repeated values: prefix1 "
                         "(2)")

        # Sample id mismatch
        prep_info = {"1.SKB2.640194": {"not_a_run_prefix": "prefix1"},
                     "1.SKM4.640180": {"not_a_run_prefix": "prefix1"},
                     "1.SKB3.640195": {"not_a_run_prefix": "prefix3"}}
        job_id = self._create_template_and_job(
            prep_info, files, "per_sample_FASTQ")
        obs_success, obs_ainfo, obs_error = _validate_per_sample_FASTQ(
            self.qclient, job_id, prep_info, files)
        self.assertFalse(obs_success)
        self.assertIsNone(obs_ainfo)
        self.assertEqual(obs_error,
                         "The provided files are not prefixed by sample id. "
                         "Please provide the 'run_prefix' column in your prep "
                         "information. Offending files: raw_forward_seqs: "
                         "prefix1_fwd.fastq, prefix2_fwd.fastq, "
                         "Aprefix3_fwd.fastq, raw_reverse_seqs: ")

    def test_create_artifact_demultipelexed_error(self):
        out_dir = mkdtemp()
        self._clean_up_files.append(out_dir)

        # Filepath type not supported
        prep_info = {"1.SKB2.640194": {"run_prefix": "prefix1"},
                     "1.SKM4.640180": {"run_prefix": "prefix2"},
                     "1.SKB3.640195": {"run_prefix": "prefix3"}}
        files = {'Unknown': ['/path/to/file1.fastq']}
        job_id = self._create_template_and_job(
            prep_info, files, "Demultiplexed")
        obs_success, obs_ainfo, obs_error = _validate_demultiplexed(
            self.qclient, job_id, prep_info, files, out_dir)
        self.assertFalse(obs_success)
        self.assertIsNone(obs_ainfo)
        self.assertEqual(obs_error,
                         "Filepath type(s) Unknown not supported by artifact "
                         "type Demultiplexed. Supported filepath types: "
                         "log, preprocessed_demux, preprocessed_fasta, "
                         "preprocessed_fastq")

        # More than a single filepath type
        files = {'preprocessed_fastq': ['/path/to/file1.fastq',
                                        '/path/to/file2.fastq']}
        job_id = self._create_template_and_job(
            prep_info, files, "Demultiplexed")
        obs_success, obs_ainfo, obs_error = _validate_demultiplexed(
            self.qclient, job_id, prep_info, files, out_dir)
        self.assertFalse(obs_success)
        self.assertIsNone(obs_ainfo)
        self.assertEqual(obs_error,
                         "Only one filepath of each file type is supported, "
                         "offending types:\n"
                         "preprocessed_fastq (2): /path/to/file1.fastq, "
                         "/path/to/file2.fastq")

        # demux, fasta and fastq not provided
        files = {'log': ['/path/to/file1.log']}
        job_id = self._create_template_and_job(
            prep_info, files, "Demultiplexed")
        obs_success, obs_ainfo, obs_error = _validate_demultiplexed(
            self.qclient, job_id, prep_info, files, out_dir)
        self.assertFalse(obs_success)
        self.assertIsNone(obs_ainfo)
        self.assertEqual(obs_error,
                         "Either a 'preprocessed_demux', 'preprocessed_fastq' "
                         "or 'preprocessed_fasta' file should be provided.")

    def _generate_files(self, sample_names):
        fd, fastq_fp = mkstemp(suffix=".fastq")
        close(fd)
        with open(fastq_fp, 'w') as f:
            f.write(FASTQ_SEQS.format(**sample_names))

        demux_fp = "%s.demux" % fastq_fp
        with File(demux_fp) as f:
            to_hdf5(fastq_fp, f)

        out_dir = mkdtemp()

        self._clean_up_files.extend([fastq_fp, demux_fp, out_dir])

        return demux_fp, fastq_fp, out_dir

    def test_validate_demux_file(self):
        demux_fp, _, out_dir = self._generate_files({'s1': 's1', 's2': 's2'})
        prep_info = {"1.SKB2.640194": {"run_prefix": "s1"},
                     "1.SKM4.640180": {"run_prefix": "s2"},
                     "1.SKB3.640195": {"run_prefix": "s3"},
                     "1.SKB6.640176": {"run_prefix": "s4"}}
        files = {'preprocessed_demux': demux_fp}
        job_id = self._create_template_and_job(
            prep_info, files, "Demultiplexed")
        obs_success, obs_ainfo, obs_error = _validate_demux_file(
            self.qclient, job_id, prep_info, out_dir, demux_fp)
        self.assertTrue(obs_success)
        name = splitext(basename(demux_fp))[0]
        exp_fastq_fp = join(out_dir, "%s.fastq" % name)
        exp_fasta_fp = join(out_dir, "%s.fasta" % name)
        exp_demux_fp = join(out_dir, basename(demux_fp))
        filepaths = [
            (exp_fastq_fp, 'preprocessed_fastq'),
            (exp_fasta_fp, 'preprocessed_fasta'),
            (exp_demux_fp, 'preprocessed_demux')]
        exp = [ArtifactInfo(None, "Demultiplexed", filepaths)]
        self.assertEqual(obs_ainfo, exp)
        self.assertEqual(obs_error, "")
        with File(exp_demux_fp) as f:
            self.assertItemsEqual(f.keys(), ["1.SKB2.640194", "1.SKM4.640180"])

    def test_validate_demux_file_infer(self):
        demux_fp, _, out_dir = self._generate_files({'s1': 'SKB2.640194',
                                                     's2': 'SKM4.640180'})
        prep_info = {"1.SKB2.640194": {"not_a_run_prefix": "s1"},
                     "1.SKM4.640180": {"not_a_run_prefix": "s2"},
                     "1.SKB3.640195": {"not_a_run_prefix": "s3"},
                     "1.SKB6.640176": {"not_a_run_prefix": "s4"}}
        files = {'preprocessed_demux': [demux_fp]}
        job_id = self._create_template_and_job(
            prep_info, files, "Demultiplexed")
        obs_success, obs_ainfo, obs_error = _validate_demux_file(
            self.qclient, job_id, prep_info, out_dir, demux_fp)
        self.assertTrue(obs_success)
        name = splitext(basename(demux_fp))[0]
        exp_fastq_fp = join(out_dir, "%s.fastq" % name)
        exp_fasta_fp = join(out_dir, "%s.fasta" % name)
        exp_demux_fp = join(out_dir, basename(demux_fp))
        filepaths = [
            (exp_fastq_fp, 'preprocessed_fastq'),
            (exp_fasta_fp, 'preprocessed_fasta'),
            (exp_demux_fp, 'preprocessed_demux')]
        exp = [ArtifactInfo(None, "Demultiplexed", filepaths)]
        self.assertEqual(obs_ainfo, exp)
        self.assertEqual(obs_error, "")
        with File(exp_demux_fp) as f:
            self.assertItemsEqual(f.keys(), ["1.SKB2.640194", "1.SKM4.640180"])

    def test_validate_demux_file_error(self):
        demux_fp, _, out_dir = self._generate_files({'s1': 's1', 's2': 's2'})

        # Run prefix not provided and demux samples do not match
        prep_info = {"1.SKB2.640194": {"not_a_run_prefix": "prefix1"},
                     "1.SKM4.640180": {"not_a_run_prefix": "prefix2"},
                     "1.SKB3.640195": {"not_a_run_prefix": "prefix3"}}
        files = {'preprocessed_demux': [demux_fp]}
        job_id = self._create_template_and_job(
            prep_info, files, "Demultiplexed")
        obs_success, obs_ainfo, obs_error = _validate_demux_file(
            self.qclient, job_id, prep_info, out_dir, demux_fp)
        self.assertFalse(obs_success)
        self.assertIsNone(obs_ainfo)
        self.assertEqual(obs_error,
                         'The sample ids in the demultiplexed files do not '
                         'match the ones in the prep information. Please, '
                         'provide the column "run_prefix" in the prep '
                         'information to map the existing sample ids to the '
                         'prep information sample ids.')

        # Incorrect run prefix column
        prep_info = {"1.SKB2.640194": {"run_prefix": "prefix1"},
                     "1.SKM4.640180": {"run_prefix": "prefix2"},
                     "1.SKB3.640195": {"run_prefix": "prefix3"}}
        job_id = self._create_template_and_job(
            prep_info, files, "Demultiplexed")
        obs_success, obs_ainfo, obs_error = _validate_demux_file(
            self.qclient, job_id, prep_info, out_dir, demux_fp)
        self.assertFalse(obs_success)
        self.assertIsNone(obs_ainfo)
        self.assertEqual(obs_error,
                         'The sample ids in the "run_prefix" columns from the '
                         'prep information do not match the ones in the demux '
                         'file. Please, correct the column "run_prefix" in '
                         'the prep information to map the existing sample ids '
                         'to the prep information sample ids.')

    def test_validate_error(self):
        parameters = {'template': 1,
                      'files': dumps(
                          {"preprocessed_demux": ["/path/file1.demux"]}),
                      'artifact_type': 'UNKNOWN'}

        data = {'command': 6,
                'parameters': dumps(parameters),
                'status': 'running'}
        job_id = self.qclient.post(
            '/apitest/processing_job/', data=data)['job']

        obs_success, obs_ainfo, obs_error = validate(
            self.qclient, job_id, parameters, '')
        self.assertFalse(obs_success)
        self.assertIsNone(obs_ainfo)
        self.assertEqual(obs_error,
                         "Unknown artifact_type UNKNOWN. Supported types: "
                         "'SFF', 'FASTQ', 'FASTA', 'FASTA_Sanger', "
                         "'per_sample_FASTQ', 'Demultiplexed'")

FASTQ_SEQS = """@{s1}_1 orig_bc=abc new_bc=abc bc_diffs=0
xyz
+
ABC
@{s2}_1 orig_bc=abw new_bc=wbc bc_diffs=4
qwe
+
DFG
@{s1}_2 orig_bc=abw new_bc=wbc bc_diffs=4
qwe
+
DEF
"""

if __name__ == '__main__':
    main()
