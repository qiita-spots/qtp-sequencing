# -----------------------------------------------------------------------------
# Copyright (c) 2014--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------

from unittest import main
from tempfile import mkdtemp, mkstemp
from os import remove, close
from os.path import exists, isdir, basename, splitext, join, dirname, abspath
from inspect import currentframe, getfile
from shutil import copyfile
from shutil import rmtree
from json import dumps

from qiita_client import ArtifactInfo
from qiita_client.testing import PluginTestCase
from h5py import File
from qiita_files.demux import to_hdf5

from qtp_sequencing.validate import (
    _validate_multiple, _validate_per_sample_FASTQ, _validate_demux_file,
    _validate_demultiplexed, validate)


class ValidateTests(PluginTestCase):
    def setUp(self):
        self.source_dir = join(dirname(abspath(getfile(currentframe()))),
                               'test_data')
        self.fastq = join(self.source_dir, 'file.fastq')
        self.fastqqz = join(self.source_dir, 'file.fastq.gz')

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
                      'analysis': None,
                      'files': dumps(files),
                      'artifact_type': atype}

        data = {'command': dumps(
            ['Sequencing Data Type', '2022.11', 'Validate']),
                'parameters': dumps(parameters),
                'status': 'running'}
        job_id = self.qclient.post(
            '/apitest/processing_job/', data=data)['job']

        return job_id, parameters

    def test_validate_multiple(self):
        test_dir = mkdtemp()
        self._clean_up_files.append(test_dir)

        copyfile(self.fastq, f'{test_dir}/prefix1.fastq')
        copyfile(self.fastq, f'{test_dir}/prefix2.fastq')
        copyfile(self.fastq, f'{test_dir}/prefix1_b.fastq')
        copyfile(self.fastq, f'{test_dir}/prefix2_b.fastq')

        prep_info = {
            '1.SKB2.640194': {'run_prefix': 'prefix1'},
            '1.SKM4.640180': {'run_prefix': 'prefix1'},
            '1.SKB3.640195': {'run_prefix': 'prefix2'}}
        files = {'raw_forward_seqs': [f'{test_dir}/prefix1.fastq',
                                      f'{test_dir}/prefix2.fastq'],
                 'raw_barcodes': [f'{test_dir}/prefix1_b.fastq',
                                  f'{test_dir}/prefix2_b.fastq']}
        atype = "FASTQ"
        job_id, _ = self._create_template_and_job(prep_info, files, atype)

        obs_success, obs_ainfo, obs_error = _validate_multiple(
            self.qclient, job_id, prep_info, files, atype)

        self.assertEqual(obs_error, "")
        self.assertTrue(obs_success)
        filepaths = [
            (f'{test_dir}/prefix1_b.fastq.gz', 'raw_barcodes'),
            (f'{test_dir}/prefix2_b.fastq.gz', 'raw_barcodes'),
            (f'{test_dir}/prefix1.fastq.gz', 'raw_forward_seqs'),
            (f'{test_dir}/prefix2.fastq.gz', 'raw_forward_seqs')]
        exp = [ArtifactInfo(None, "FASTQ", filepaths)]
        self.assertEqual(obs_ainfo, exp)

    def test_validate_multiple_single_lane(self):
        test_dir = mkdtemp()
        self._clean_up_files.append(test_dir)

        copyfile(self.fastq, f'{test_dir}/prefix1.fastq')
        copyfile(self.fastq, f'{test_dir}/prefix1_b.fastq')

        prep_info = {"1.SKB2.640194": {"not_a_run_prefix": "prefix1"},
                     "1.SKM4.640180": {"not_a_run_prefix": "prefix1"},
                     "1.SKB3.640195": {"not_a_run_prefix": "prefix2"}}
        files = {'raw_forward_seqs': [f'{test_dir}/prefix1.fastq'],
                 'raw_barcodes': [f'{test_dir}/prefix1_b.fastq']}
        atype = "FASTQ"
        job_id, _ = self._create_template_and_job(prep_info, files, atype)

        obs_success, obs_ainfo, obs_error = _validate_multiple(
            self.qclient, job_id, prep_info, files, atype)

        self.assertEqual(obs_error, "")
        self.assertTrue(obs_success)
        filepaths = [
            (f'{test_dir}/prefix1_b.fastq.gz', 'raw_barcodes'),
            (f'{test_dir}/prefix1.fastq.gz', 'raw_forward_seqs')]
        exp = [ArtifactInfo(None, atype, filepaths)]
        self.assertEqual(obs_ainfo, exp)

    def test_validate_multiple_error(self):
        # Filepath type not supported
        prep_info = {"1.SKB2.640194": {"run_prefix": "prefix1"},
                     "1.SKM4.640180": {"run_prefix": "prefix1"},
                     "1.SKB3.640195": {"run_prefix": "prefix2"}}
        files = {'Unknown': ['/path/to/file1.fastq']}
        atype = "FASTQ"
        job_id, _ = self._create_template_and_job(prep_info, files, atype)
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
        job_id, _ = self._create_template_and_job(prep_info, files, atype)
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
        self.assertCountEqual(obs_error.split('\n'), error.split('\n'))

        # File doesn't match any run prefix
        files = {'raw_forward_seqs': ['/path/to/file1.fastq',
                                      '/path/to/prefix2.fastq'],
                 'raw_barcodes': ['/path/to/file1_b.fastq',
                                  '/path/to/prefix2_b.fastq']}
        job_id, _ = self._create_template_and_job(prep_info, files, atype)
        obs_success, obs_ainfo, obs_error = _validate_multiple(
            self.qclient, job_id, prep_info, files, atype)
        error = ("Error creating artifact. Offending files:\nraw_forward_seqs:"
                 " The provided files do not match the run prefix values in "
                 "the prep information: file1.fastq\n"
                 "raw_barcodes: The provided files do not match the run "
                 "prefix values in the prep information: file1_b.fastq")
        self.assertFalse(obs_success)
        self.assertIsNone(obs_ainfo)
        self.assertCountEqual(obs_error.split('\n'), error.split('\n'))

        # A required filepath type is missing
        files = {'raw_forward_seqs': ['/path/to/prefix1.fastq',
                                      '/path/to/prefix2.fastq'],
                 'raw_reverse_seqs': ['/path/to/prefix1_rev.fastq',
                                      '/path/to/prefix2_rev.fastq']}
        job_id, _ = self._create_template_and_job(prep_info, files, atype)
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
        job_id, _ = self._create_template_and_job(prep_info, files, atype)
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
        self.assertCountEqual(obs_error.split('\n'), error.split('\n'))

    def test_validate_SFF(self):
        prep_info = {"1.SKB2.640194": {"run_prefix": "GAX40"},
                     "1.SKM4.640180": {"run_prefix": "GAX40"},
                     "1.SKB3.640195": {"run_prefix": "GAX50"}}
        files = {'raw_sff': ['/path/to/GAX401.sff', '/path/to/GAX402.sff',
                             '/path/to/GAX501.sff']}
        job_id, _ = self._create_template_and_job(prep_info, files, "SFF")
        obs_success, obs_ainfo, obs_error = _validate_multiple(
            self.qclient, job_id, prep_info, files, 'SFF')
        self.assertTrue(obs_success)
        filepaths = [('/path/to/GAX401.sff', 'raw_sff'),
                     ('/path/to/GAX402.sff', 'raw_sff'),
                     ('/path/to/GAX501.sff', 'raw_sff')]
        exp = [ArtifactInfo(None, "SFF", filepaths)]
        self.assertEqual(obs_ainfo, exp)
        self.assertEqual(obs_error, "")

        # let's test a failure
        files = {'raw_sff': ['/path/to/GAX401.sff', '/path/to/GAX402.sff']}
        job_id, _ = self._create_template_and_job(prep_info, files, "SFF")
        obs_success, obs_ainfo, obs_error = _validate_multiple(
            self.qclient, job_id, prep_info, files, 'SFF')
        error = ("Error creating artifact. Offending files:\nraw_sff: The "
                 "following run prefixes in the prep information file do not "
                 "match any file: GAX50")
        self.assertFalse(obs_success)
        self.assertIsNone(obs_ainfo)
        self.assertCountEqual(obs_error, error)

    def test_validate_per_sample_FASTQ_run_prefix(self):
        f1 = join(self.source_dir, 'SKB2.640194_file.fastq')
        f2 = join(self.source_dir, 'SKM4.640180_file.fastq')
        f3 = join(self.source_dir, 'SKB3.640195_file.fastq')
        raw_files = [f1, f2, f3]
        for x in raw_files:
            copyfile(self.fastq, x)
            self._clean_up_files.append(x)
        prep_info = {"1.SKB2.640194": {"run_prefix": "prefix1"},
                     "1.SKM4.640180": {"run_prefix": "prefix2"},
                     "1.SKB3.640195": {"run_prefix": "prefix3"}}
        files = {'raw_forward_seqs': raw_files}
        job_id, _ = self._create_template_and_job(
            prep_info, files, "per_sample_FASTQ")
        obs_success, obs_ainfo, obs_error = _validate_per_sample_FASTQ(
            self.qclient, job_id, prep_info, files)
        self.assertEqual(obs_error, "")
        self.assertTrue(obs_success)
        filepaths = [('%s.gz' % x, 'raw_forward_seqs') for x in raw_files]
        exp = [ArtifactInfo(None, "per_sample_FASTQ", filepaths)]
        self.assertEqual(obs_ainfo, exp)

    def test_validate_per_sample_FASTQ(self):
        prep_info = {"1.SKB2.640194": {"not_a_run_prefix": "prefix1"},
                     "1.SKM4.640180": {"not_a_run_prefix": "prefix1"},
                     "1.SKB3.640195": {"not_a_run_prefix": "prefix2"}}
        files = {'raw_forward_seqs': ['/path/to/SKB2.640194_file.fastq',
                                      '/path/to/SKM4.640180_file.fastq',
                                      '/path/to/SKB3.640195_file.fastq']}
        job_id, _ = self._create_template_and_job(
            prep_info, files, "per_sample_FASTQ")
        obs_success, obs_ainfo, obs_error = _validate_per_sample_FASTQ(
            self.qclient, job_id, prep_info, files, True)
        self.assertEqual(obs_error, 'Some of the files are empty: '
                         'SKB2.640194_file.fastq, SKM4.640180_file.fastq, '
                         'SKB3.640195_file.fastq')
        self.assertFalse(obs_success)

        f1 = join(self.source_dir, 'SKB2.640194_file.fastq')
        f2 = join(self.source_dir, 'SKM4.640180_file.fastq')
        f3 = join(self.source_dir, 'SKB3.640195_file.fastq')
        raw_files = [f1, f2, f3]
        for x in raw_files:
            copyfile(self.fastq, x)
            self._clean_up_files.append(x)

        prep_info = {"1.SKB2.640194": {"not_a_run_prefix": "prefix1"},
                     "1.SKM4.640180": {"not_a_run_prefix": "prefix1"},
                     "1.SKB3.640195": {"not_a_run_prefix": "prefix2"}}
        files = {'raw_forward_seqs': raw_files}
        job_id, _ = self._create_template_and_job(
            prep_info, files, "per_sample_FASTQ")
        obs_success, obs_ainfo, obs_error = _validate_per_sample_FASTQ(
            self.qclient, job_id, prep_info, files)
        self.assertTrue(obs_success)

        filepaths = [('%s.gz' % x, 'raw_forward_seqs') for x in raw_files]
        exp = [ArtifactInfo(None, "per_sample_FASTQ", filepaths)]
        self.assertEqual(obs_ainfo, exp)
        self.assertEqual(obs_error, "")

    def test_validate_per_sample_FASTQ_preprocessed_fastq(self):
        f1 = join(self.source_dir, 'SKB2.640194_file.fastq')
        f2 = join(self.source_dir, 'SKM4.640180_file.fastq')
        f3 = join(self.source_dir, 'SKB3.640195_file.fastq')
        copyfile(self.fastq, f1)
        copyfile(self.fastq, f2)
        copyfile(self.fastq, f3)
        self._clean_up_files.append(f1)
        self._clean_up_files.append(f2)
        self._clean_up_files.append(f3)

        prep_info = {"1.SKB2.640194": {"not_a_run_prefix": "prefix1"},
                     "1.SKM4.640180": {"not_a_run_prefix": "prefix1"},
                     "1.SKB3.640195": {"not_a_run_prefix": "prefix2"}}
        files = {'preprocessed_fastq': [f1, f2, f3]}
        job_id, _ = self._create_template_and_job(
            prep_info, files, "per_sample_FASTQ")
        obs_success, obs_ainfo, obs_error = _validate_per_sample_FASTQ(
            self.qclient, job_id, prep_info, files)
        self.assertTrue(obs_success)
        filepaths = [(f1 + '.gz', 'preprocessed_fastq'),
                     (f2 + '.gz', 'preprocessed_fastq'),
                     (f3 + '.gz', 'preprocessed_fastq')]
        exp = [ArtifactInfo(None, "per_sample_FASTQ", filepaths)]
        self.assertEqual(obs_ainfo, exp)
        self.assertEqual(obs_error, "")
        # making sure the regular fastq files doesn't exist anymore but
        # the gz do
        self.assertFalse(exists(f1))
        self.assertTrue(exists(f1 + '.gz'))

        f1 = join(self.source_dir, 'SKB2.640194_file_R1.fastq')
        f2 = join(self.source_dir, 'SKB2.640194_file_R2.fastq')
        f3 = join(self.source_dir, 'SKB2.640194_file_unmatched_R1.fastq')
        f4 = join(self.source_dir, 'SKB2.640194_file_unmatched_R2.fastq')
        f5 = join(self.source_dir, 'SKM4.640180_file_R1.fastq')
        f6 = join(self.source_dir, 'SKM4.640180_file_R2.fastq')
        f7 = join(self.source_dir, 'SKM4.640180_file_unmatched_R1.fastq')
        f8 = join(self.source_dir, 'SKM4.640180_file_unmatched_R2.fastq')
        f9 = join(self.source_dir, 'SKB3.640195_file_R1.fastq')
        fA = join(self.source_dir, 'SKB3.640195_file_R2.fastq')
        fB = join(self.source_dir, 'SKB3.640195_file_unmatched_R1.fastq')
        fC = join(self.source_dir, 'SKB3.640195_file_unmatched_R2.fastq')
        raw_files = [f1, f2, f3, f4, f5, f6, f7, f8, f9, fA, fB, fC]
        for x in raw_files:
            copyfile(self.fastq, x)
            self._clean_up_files.append(x)

        prep_info = {"1.SKB2.640194": {"not_a_run_prefix": "prefix1"},
                     "1.SKM4.640180": {"not_a_run_prefix": "prefix1"},
                     "1.SKB3.640195": {"not_a_run_prefix": "prefix2"}}
        files = {'preprocessed_fastq': raw_files}
        job_id, _ = self._create_template_and_job(
            prep_info, files, "per_sample_FASTQ")
        obs_success, obs_ainfo, obs_error = _validate_per_sample_FASTQ(
            self.qclient, job_id, prep_info, files)
        self.assertEqual(obs_error, "")
        self.assertTrue(obs_success)
        filepaths = [('%s.gz' % x, 'preprocessed_fastq') for x in raw_files]
        exp = [ArtifactInfo(None, "per_sample_FASTQ", filepaths)]
        self.assertEqual(obs_ainfo, exp)

    def test_validate_per_sample_FASTQ_error(self):
        # Filepath type not supported
        prep_info = {"1.SKB2.640194": {"run_prefix": "prefix1"},
                     "1.SKM4.640180": {"run_prefix": "prefix2"},
                     "1.SKB3.640195": {"run_prefix": "prefix3"}}
        files = {'Unknown': ['/path/to/file1.fastq']}
        job_id, _ = self._create_template_and_job(
            prep_info, files, "per_sample_FASTQ")
        obs_success, obs_ainfo, obs_error = _validate_per_sample_FASTQ(
            self.qclient, job_id, prep_info, files)
        self.assertFalse(obs_success)
        self.assertIsNone(obs_ainfo)
        self.assertEqual(obs_error,
                         "Filepath type(s) Unknown not supported by artifact "
                         "type per_sample_FASTQ. Supported filepath types: "
                         "raw_forward_seqs, raw_reverse_seqs, "
                         "preprocessed_fastq")

        # Missing raw_forward_seqs and preprocessed_fastq
        files = {'raw_reverse_seqs': ['/path/to/file1.fastq']}
        job_id, _ = self._create_template_and_job(
            prep_info, files, "per_sample_FASTQ")
        obs_success, obs_ainfo, obs_error = _validate_per_sample_FASTQ(
            self.qclient, job_id, prep_info, files)
        self.assertFalse(obs_success)
        self.assertIsNone(obs_ainfo)
        self.assertEqual(obs_error,
                         "Missing required filepath type: raw_forward_seqs "
                         "or preprocessed_fastq")

        # Raw forward seqs and preprocessed_fastq
        files = {'raw_forward_seqs': ['/path/to/file1.fastq'],
                 'preprocessed_fastq': ['/path/to/file1.fastq']}
        job_id, _ = self._create_template_and_job(
            prep_info, files, "per_sample_FASTQ")
        obs_success, obs_ainfo, obs_error = _validate_per_sample_FASTQ(
            self.qclient, job_id, prep_info, files)
        self.assertFalse(obs_success)
        self.assertIsNone(obs_ainfo)
        self.assertEqual(obs_error,
                         "If raw_forward_seqs is provided, preprocessed_fastq "
                         "should not be provided")

        # Preprocessed fastq and raw_reverse_seqs
        files = {'raw_reverse_seqs': ['/path/to/file1.fastq'],
                 'preprocessed_fastq': ['/path/to/file1.fastq']}
        job_id, _ = self._create_template_and_job(
            prep_info, files, "per_sample_FASTQ")
        obs_success, obs_ainfo, obs_error = _validate_per_sample_FASTQ(
            self.qclient, job_id, prep_info, files)
        self.assertFalse(obs_success)
        self.assertIsNone(obs_ainfo)
        self.assertEqual(obs_error,
                         "If preprocessed_fastq is provided, raw_reverse_seqs "
                         "should not be provided")

        # Run prefix mismatch
        files = {'raw_forward_seqs': ['/path/to/prefix1_fwd.fastq',
                                      '/path/to/prefix2_fwd.fastq',
                                      '/path/to/Aprefix3_fwd.fastq']}
        job_id, _ = self._create_template_and_job(
            prep_info, files, "per_sample_FASTQ")
        obs_success, obs_ainfo, obs_error = _validate_per_sample_FASTQ(
            self.qclient, job_id, prep_info, files)
        self.assertFalse(obs_success)
        self.assertIsNone(obs_ainfo)
        self.assertEqual(obs_error,
                         'The provided files are not prefixed by sample id or '
                         'do not match the run prefix values in the prep '
                         'information. Offending files:\n raw_forward_seqs: '
                         'Aprefix3_fwd.fastq\nraw_reverse_seqs: ')

        # Non-unique run-prefix values
        prep_info = {"1.SKB2.640194": {"run_prefix": "prefix1"},
                     "1.SKM4.640180": {"run_prefix": "prefix1"},
                     "1.SKB3.640195": {"run_prefix": "prefix3"}}
        job_id, _ = self._create_template_and_job(
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
        job_id, _ = self._create_template_and_job(
            prep_info, files, "per_sample_FASTQ")
        obs_success, obs_ainfo, obs_error = _validate_per_sample_FASTQ(
            self.qclient, job_id, prep_info, files)
        self.assertFalse(obs_success)
        self.assertIsNone(obs_ainfo)
        self.assertEqual(obs_error,
                         'The provided files are not prefixed by sample id. '
                         'Offending files:\n raw_forward_seqs: '
                         'prefix1_fwd.fastq, prefix2_fwd.fastq, '
                         'Aprefix3_fwd.fastq\nraw_reverse_seqs: ')

    def test_create_artifact_demultipelexed_error(self):
        out_dir = mkdtemp()
        self._clean_up_files.append(out_dir)

        # Filepath type not supported
        prep_info = {"1.SKB2.640194": {"run_prefix": "prefix1"},
                     "1.SKM4.640180": {"run_prefix": "prefix2"},
                     "1.SKB3.640195": {"run_prefix": "prefix3"}}
        files = {'Unknown': ['/path/to/file1.fastq']}
        job_id, _ = self._create_template_and_job(
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
        job_id, _ = self._create_template_and_job(
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
        job_id, _ = self._create_template_and_job(
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
        with File(demux_fp, 'w') as f:
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
        job_id, _ = self._create_template_and_job(
            prep_info, files, "Demultiplexed")
        obs_success, obs_ainfo, obs_error = _validate_demux_file(
            self.qclient, job_id, prep_info, out_dir, demux_fp)
        self.assertTrue(obs_success)
        name = splitext(basename(demux_fp))[0]
        exp_fastq_fp = join(out_dir, "%s.fastq.gz" % name)
        exp_fasta_fp = join(out_dir, "%s.fasta.gz" % name)
        exp_demux_fp = join(out_dir, basename(demux_fp))
        filepaths = [
            (exp_fastq_fp, 'preprocessed_fastq'),
            (exp_fasta_fp, 'preprocessed_fasta'),
            (exp_demux_fp, 'preprocessed_demux')]
        exp = [ArtifactInfo(None, "Demultiplexed", filepaths)]
        self.assertEqual(obs_ainfo, exp)
        self.assertEqual(obs_error, "")
        with File(exp_demux_fp) as f:
            self.assertCountEqual(f.keys(), ["1.SKB2.640194", "1.SKM4.640180"])

    def test_validate_demux_file_without_demux(self):
        demux_fp, fastq_fp, out_dir = self._generate_files(
            {'s1': 's1', 's2': 's2'})
        remove(demux_fp)
        prep_info = {"1.SKB2.640194": {"run_prefix": "s1"},
                     "1.SKM4.640180": {"run_prefix": "s2"},
                     "1.SKB3.640195": {"run_prefix": "s3"},
                     "1.SKB6.640176": {"run_prefix": "s4"}}
        files = {'preprocessed_fastq': [fastq_fp]}
        job_id, _ = self._create_template_and_job(
            prep_info, files, "Demultiplexed")
        obs_success, obs_ainfo, obs_error = _validate_demultiplexed(
            self.qclient, job_id, prep_info, files, out_dir)
        self.assertEqual(obs_error, "")
        self.assertTrue(obs_success)

        # we are gonna just check that the demux file looks good and [0] is
        # because this only returns one element in the object list
        demux = [f[0] for f in obs_ainfo[0].files
                 if f[1] == 'preprocessed_demux'][0]
        with File(demux) as f:
            self.assertCountEqual(f.keys(), ["1.SKB2.640194", "1.SKM4.640180"])

    def test_validate_demux_file_infer(self):
        demux_fp, _, out_dir = self._generate_files({'s1': 'SKB2.640194',
                                                     's2': 'SKM4.640180'})
        prep_info = {"1.SKB2.640194": {"not_a_run_prefix": "s1"},
                     "1.SKM4.640180": {"not_a_run_prefix": "s2"},
                     "1.SKB3.640195": {"not_a_run_prefix": "s3"},
                     "1.SKB6.640176": {"not_a_run_prefix": "s4"}}
        files = {'preprocessed_demux': [demux_fp]}
        job_id, _ = self._create_template_and_job(
            prep_info, files, "Demultiplexed")
        obs_success, obs_ainfo, obs_error = _validate_demux_file(
            self.qclient, job_id, prep_info, out_dir, demux_fp)
        self.assertTrue(obs_success)
        name = splitext(basename(demux_fp))[0]
        exp_fastq_fp = join(out_dir, "%s.fastq.gz" % name)
        exp_fasta_fp = join(out_dir, "%s.fasta.gz" % name)
        exp_demux_fp = join(out_dir, basename(demux_fp))
        filepaths = [
            (exp_fastq_fp, 'preprocessed_fastq'),
            (exp_fasta_fp, 'preprocessed_fasta'),
            (exp_demux_fp, 'preprocessed_demux')]
        exp = [ArtifactInfo(None, "Demultiplexed", filepaths)]
        self.assertEqual(obs_ainfo, exp)
        self.assertEqual(obs_error, "")
        with File(exp_demux_fp) as f:
            self.assertCountEqual(f.keys(), ["1.SKB2.640194", "1.SKM4.640180"])

    def test_validate_demux_file_error(self):
        demux_fp, _, out_dir = self._generate_files({'s1': 's1', 's2': 's2'})

        # Run prefix not provided and demux samples do not match
        prep_info = {"1.SKB2.640194": {"not_a_run_prefix": "prefix1"},
                     "1.SKM4.640180": {"not_a_run_prefix": "prefix2"},
                     "1.SKB3.640195": {"not_a_run_prefix": "prefix3"}}
        files = {'preprocessed_demux': [demux_fp]}
        job_id, _ = self._create_template_and_job(
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
        job_id, _ = self._create_template_and_job(
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
                      'analysis': None,
                      'files': dumps(
                          {"preprocessed_demux": ["/path/file1.demux"]}),
                      'artifact_type': 'UNKNOWN'}

        data = {'command': dumps(
            ['Sequencing Data Type', '2022.11', 'Validate']),
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
                         "'per_sample_FASTQ', 'FASTA_preprocessed', "
                         "'Demultiplexed'")

    def test_validate_success(self):
        test_dir = mkdtemp()
        out_dir = mkdtemp()
        self._clean_up_files.append(test_dir)
        self._clean_up_files.append(out_dir)

        copyfile(self.fastq, f'{test_dir}/prefix1.fastq')
        copyfile(self.fastq, f'{test_dir}/prefix2.fastq')
        copyfile(self.fastq, f'{test_dir}/prefix1_b.fastq')
        copyfile(self.fastq, f'{test_dir}/prefix2_b.fastq')

        prep_info = {
            '1.SKB2.640194': {'run_prefix': 'prefix1'},
            '1.SKM4.640180': {'run_prefix': 'prefix1'},
            '1.SKB3.640195': {'run_prefix': 'prefix2'}}
        files = {'raw_forward_seqs': [f'{test_dir}/prefix1.fastq',
                                      f'{test_dir}/prefix2.fastq'],
                 'raw_barcodes': [f'{test_dir}/prefix1_b.fastq',
                                  f'{test_dir}/prefix2_b.fastq']}
        atype = "FASTQ"
        job_id, params = self._create_template_and_job(prep_info, files, atype)

        obs_success, obs_ainfo, obs_error = validate(
            self.qclient, job_id, params, out_dir)
        self.assertEqual(obs_error, '')
        self.assertTrue(obs_success)
        self.assertEqual(len(obs_ainfo), 1)
        exp_files = [
            (f'{test_dir}/prefix1.fastq.gz', 'raw_forward_seqs'),
            (f'{test_dir}/prefix2.fastq.gz', 'raw_forward_seqs'),
            (f'{test_dir}/prefix1_b.fastq.gz', 'raw_barcodes'),
            (f'{test_dir}/prefix2_b.fastq.gz', 'raw_barcodes'),
            (f'{out_dir}/index.html', 'html_summary')]
        self.assertCountEqual(obs_ainfo[0].files, exp_files)

    def test_validate_FASTA_preprocessed(self):
        prep_info = {"1.SKB2.640194": {"run_prefix": "s1"},
                     "1.SKM4.640180": {"run_prefix": "s2"},
                     "1.SKB3.640195": {"run_prefix": "s3"},
                     "1.SKB6.640176": {"run_prefix": "s4"}}

        files = {'preprocessed_fasta': [
            '/path/to/s1_blah.fna', '/path/to/s2_blah.fna',
            '/path/to/s3_blah.fna', '/path/to/s4_blah.fna']}

        out_dir = mkdtemp()
        self._clean_up_files.append(out_dir)

        atype = 'FASTA_preprocessed'
        job_id, _ = self._create_template_and_job(prep_info, files, atype)
        obs_success, obs_ainfo, obs_error = _validate_multiple(
            self.qclient, job_id, prep_info, files, atype, True)
        self.assertEqual(obs_error, "")
        self.assertTrue(obs_success)

        files = [('/path/to/s1_blah.fna.gz', 'preprocessed_fasta'),
                 ('/path/to/s2_blah.fna.gz', 'preprocessed_fasta'),
                 ('/path/to/s3_blah.fna.gz', 'preprocessed_fasta'),
                 ('/path/to/s4_blah.fna.gz', 'preprocessed_fasta')]
        exp = [ArtifactInfo(None, atype, files)]
        self.assertEqual(obs_ainfo, exp)


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
