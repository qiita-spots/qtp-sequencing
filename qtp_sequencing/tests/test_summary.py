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

from qtp_sequencing.summary import (
    generate_html_summary, _summary_demultiplexed, _summary_not_demultiplexed,
    _summary_FASTA_preprocessed)


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

    def test_generate_html_summary_no_demux(self):
        # Create a job in Qiita
        artifact_id = 1
        parameters = {'input_data': artifact_id}
        data = {'command': dumps(['Sequencing Data Type', '2021.07',
                                  'Generate HTML summary']),
                'parameters': dumps(parameters),
                'status': 'running'}
        job_id = self.qclient.post(
            '/apitest/processing_job/', data=data)['job']

        # Qiita will return a filepath, but in the test environment, these
        # files do not exist - create them
        files = self.qclient.get(
            '/qiita_db/artifacts/%s/' % artifact_id)['files']

        bcds_fp = files['raw_barcodes'][0]
        self._clean_up_files.append(bcds_fp)
        with GzipFile(bcds_fp, mode='w', mtime=1) as fh:
            fh.write(BARCODES.encode())
        fwd_fp = files['raw_forward_seqs'][0]
        self._clean_up_files.append(fwd_fp)
        with GzipFile(fwd_fp, mode='w', mtime=1) as fh:
            fh.write(READS.encode())

        # Run the test
        obs_success, obs_ainfo, obs_error = generate_html_summary(
            self.qclient, job_id, parameters, self.out_dir)

        # asserting reply
        self.assertTrue(obs_success)
        self.assertIsNone(obs_ainfo)
        self.assertEqual(obs_error, "")

        # asserting content of html
        res = self.qclient.get("/qiita_db/artifacts/%s/" % artifact_id)
        html_fp = res['files']['html_summary'][0]
        self._clean_up_files.append(html_fp)
        with open(html_fp) as html_f:
            html = html_f.read()

        self.assertEqual(html, EXP_HTML)

    def test_generate_html_summary_demux(self):
        artifact_id = 2
        parameters = {'input_data': artifact_id}
        data = {'command': dumps(['Sequencing Data Type', '2021.07',
                                  'Generate HTML summary']),
                'parameters': dumps(parameters),
                'stuatus': 'running'}
        job_id = self.qclient.post(
            '/apitest/processing_job/', data=data)['job']

        # Qiita will return a filepath, but in the test environment, these
        # files fo not exist create them
        files = self.qclient.get(
            '/qiita_db/artifacts/%s/' % artifact_id)['files']
        demux_fp = files['preprocessed_demux'][0]
        copyfile(join(dirname(__file__), 'test_data', '101_seqs.demux'),
                 demux_fp)
        self._clean_up_files.append(demux_fp)

        obs_success, obs_ainfo, obs_error = generate_html_summary(
            self.qclient, job_id, parameters, self.out_dir)

        # asserting reply
        self.assertTrue(obs_success)
        self.assertIsNone(obs_ainfo)
        self.assertEqual(obs_error, "")

        # asserting content of html
        res = self.qclient.get("/qiita_db/artifacts/%s/" % artifact_id)
        html_fp = res['files']['html_summary'][0]
        self._clean_up_files.append(html_fp)

        with open(html_fp) as html_f:
            html = html_f.read()
        self.assertRegex(html, '\n'.join(EXP_HTML_DEMUX_REGEXP))

    def test_summary_not_demultiplexed_gzipped_no_header(self):
        artifact_type = 'SFF'
        filepaths = {'raw_sff': [join(dirname(__file__), 'test_data',
                                      'Fasting_Example.sff.gz')]}

        obs = _summary_not_demultiplexed(artifact_type, filepaths)
        exp = ['<table border="1" class="dataframe">', '  <thead>',
               '    <tr style="text-align: right;">',
               '      <th>filename</th>', '      <th>md5</th>',
               '      <th>file_type</th>', '    </tr>', '  </thead>',
               '  <tbody>', '    <tr>',
               '      <td>Fasting_Example.sff.gz</td>',
               '      <td>2b636ee4833270da41a39aa5fc0a1622</td>',
               '      <td>raw_sff</td>', '    </tr>', '  </tbody>', '</table>']
        self.assertEqual(obs, '\n'.join(exp))

    def test_summary_not_demultiplexed_gzipped_header(self):
        test_dir = mkdtemp()
        self._clean_up_files.append(test_dir)

        bcds_fp = join(test_dir, 'barcodes.fastq.gz')
        with GzipFile(bcds_fp, mode='w', mtime=1) as fh:
            fh.write(BARCODES.encode())

        fwd_fp = join(test_dir, 'reads.fastq.gz')
        with GzipFile(fwd_fp, mode='w', mtime=1) as fh:
            fh.write(READS.encode())

        artifact_type = 'FASTQ'
        filepaths = {'raw_forward_seqs': [fwd_fp],
                     'raw_barcodes': [bcds_fp]}

        obs = _summary_not_demultiplexed(artifact_type, filepaths)
        exp = ['<table border="1" class="dataframe">', '  <thead>',
               '    <tr style="text-align: right;">',
               '      <th>filename</th>', '      <th>md5</th>',
               '      <th>reads</th>', '      <th>file_type</th>', '    </tr>',
               '  </thead>', '  <tbody>', '    <tr>',
               '      <td>barcodes.fastq.gz</td>',
               '      <td>1c6eefa11d8641a6f853b56801351e5a</td>',
               '      <td>4</td>', '      <td>raw_barcodes</td>', '    </tr>',
               '    <tr>', '      <td>reads.fastq.gz</td>',
               '      <td>f3909bcab34565a5d4b88c300d40bbfc</td>',
               '      <td>4</td>', '      <td>raw_forward_seqs</td>',
               '    </tr>', '  </tbody>', '</table>']
        self.assertCountEqual(obs, '\n'.join(exp))

    def test_summary_not_demultiplexed_not_gzipped_header(self):
        test_dir = mkdtemp()
        self._clean_up_files.append(test_dir)

        bcds_fp = join(test_dir, 'barcodes.fastq')
        with open(bcds_fp, 'w') as fh:
            fh.write(BARCODES)

        fwd_fp = join(test_dir, 'reads.fastq')
        with open(fwd_fp, 'w') as fh:
            fh.write(READS)

        artifact_type = 'FASTQ'
        filepaths = {'raw_forward_seqs': [fwd_fp],
                     'raw_barcodes': [bcds_fp]}

        obs = _summary_not_demultiplexed(artifact_type, filepaths)
        exp = ['<table border="1" class="dataframe">', '  <thead>',
               '    <tr style="text-align: right;">',
               '      <th>filename</th>', '      <th>md5</th>',
               '      <th>reads</th>', '      <th>file_type</th>',
               '    </tr>', '  </thead>', '  <tbody>',
               '    <tr>', '      <td>barcodes.fastq</td>',
               '      <td>d0478899bf5c6ba6019b596690e8c261</td>',
               '      <td>4</td>', '      <td>raw_barcodes</td>',
               '    </tr>', '    <tr>', '      <td>reads.fastq</td>',
               '      <td>97328e860ef506f7b029997b12bf9885</td>',
               '      <td>4</td>', '      <td>raw_forward_seqs</td>',
               '    </tr>',  '  </tbody>', '</table>']

        self.assertCountEqual(obs, '\n'.join(exp))

    def test_summary_not_demultiplexed_empty(self):
        test_dir = mkdtemp()
        self._clean_up_files.append(test_dir)

        fwd_fp = join(test_dir, 'reads.fastq')
        with open(fwd_fp, 'w'):
            pass

        artifact_type = 'per_sample_FASTQ'
        filepaths = {'preprocessed_fastq': [fwd_fp]}
        obs = _summary_not_demultiplexed(artifact_type, filepaths)
        exp = [
            '<table border="1" class="dataframe">', '  <thead>',
            '    <tr style="text-align: right;">', '      <th>filename</th>',
            '      <th>md5</th>', '      <th>reads</th>',
            '      <th>file_type</th>', '    </tr>', '  </thead>', '  <tbody>',
            '    <tr>', '      <td>reads.fastq</td>',
            '      <td>d41d8cd98f00b204e9800998ecf8427e</td>',
            '      <td>0</td>', '      <td>preprocessed_fastq</td>',
            '    </tr>', '  </tbody>', '</table>']
        self.assertCountEqual(obs, '\n'.join(exp))

    def test_summary_demultiplexed(self):
        artifact_type = 'Demultiplexed'
        filepaths = {
            'preprocessed_demux': [join(dirname(__file__), 'test_data',
                                   '101_seqs.demux')],
            'preprocessed_fastq': ['ignored']}

        obs = _summary_demultiplexed(artifact_type, filepaths)
        exp = ['<h3>Features</h3>',
               '<b>Total</b>: 49', '<br/>',
               '<b>Max</b>: 151', '<br/>',
               '<b>Mean</b>: 151', '<br/>',
               '<b>Standard deviation</b>: 151', '<br/>',
               '<b>Median</b>: 0', '<br/>',
               '<img src = "data:image/png;base64,.*"/>']
        self.assertEqual(obs[:-1], exp[:-1])
        self.assertRegex(obs[-1], exp[-1])

    def test_summary_demultiplexed_no_demux(self):
        artifact_type = 'Demultiplexed'
        filepaths = {'preprocessed_fastq': ['ignored']}

        obs = _summary_demultiplexed(artifact_type, filepaths)
        self.assertIsNone(obs)

    def test_summary_FASTA_preprocessed(self):
        artifact_type = 'FASTA_preprocessed'

        indir = join(dirname(__file__), 'test_data')
        input = join(indir, 'scaffolds.fasta.gz')
        fna_files = []
        for rp in ['s1', 's2', 's3', 's4']:
            output = join(indir, '%s.fasta.gz' % rp)
            copyfile(input, output)
            fna_files.append(output)
        files = {'preprocessed_fasta': fna_files}

        obs = _summary_FASTA_preprocessed(
            artifact_type, files, self.out_dir)
        exp = ['<!DOCTYPE html>\n', '<html lang="en">\n', '<head>\n',
               '    <meta http-equiv="Content-Type" content="text/html;'
               'charset=utf-8" >\n']
        self.assertEqual(obs[:4], exp)


READS = """@MISEQ03:123:000000000-A40KM:1:1101:14149:1572 1:N:0:TCCACAGGAGT
GGGGGGTGCCAGCCGCCGCGGTAATACGGGGGGGGCAAGCGTTGTTCGGAATTACTGGGCGTAAAGGGCTCGTAGGCG\
GCCCACTAAGTCAGACGTGAAATCCCTCGGCTTAACCGGGGAACTGCGTCTGATACTGGATGGCTTGAGGTTGGGAGA\
GGGATGCGGAATTCCAGGTGTAGCGGTGAAATGCGTAGATATCTGGAGGAACACCGGTGGCGAAGGCGGCATCCTGGA\
CCAATTCTGACGCTGAG
+
CCCCCCCCCFFFGGGGGGGGGGGGHHHHGGGGGFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF\
FFF-.;FFFFFFFFF9@EFFFFFFFFFFFFFFFFFFF9CFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFBFFFFFFF\
FFFFFFECFFFEFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFBFFFFFFFFFFFFFFFFF;CDEA@FFFFFFFFF\
FFFFFFFFFFFFFFFFF
@MISEQ03:123:000000000-A40KM:1:1101:14170:1596 1:N:0:TCCACAGGAGT
ATGGCGTGCCAGCAGCCGCGGTAATACGGAGGGGGCTAGCGTTGTTCGGAATTACTGGGCGTAAAGCGCACGTAGGCG\
GCTTTGTAAGTTAGAGGTGAAAGCCCGGGGCTCAACTCCGGAACTGCCTTTAAGACTGCATCGCTAGAATTGTGGAGA\
GGTGAGTGGAATTCCGAGTGTAGAGGTGAAATTCGTAGATATTCGGAAGAACACCAGTGGCGAAGGCGACTCACTGGA\
CACATATTGACGCTGAG
+
CCCCCCCCCCFFGGGGGGGGGGGGGHHHGGGGGGGGGHHHGGGGGHHGGGGGHHHHHHHHGGGGHHHGGGGGHHHGHG\
GGGGHHHHHHHHGHGHHHHHHGHHHHGGGGGGHHHHHHHGGGGGHHHHHHHHHGGFGGGGGGGGGGGGGGGGGGGGGG\
GGFGGFFFFFFFFFFFFFFFFFF0BFFFFFFFFFFFFEFFFFFFFFFFFFFFFFFFFFFFFFFFFDFAFCFFFFFFFF\
FFFFFFFFFFBDFFFFF
@MISEQ03:123:000000000-A40KM:1:1101:14740:1607 1:N:0:TCCACAGGAGT
AGTGTGTGCCAGCAGCCGCGGTAATACGTAGGGTGCGAGCGTTAATCGGAATTACTGGGCGTAAAGCGTGCGCAGGCG\
GTTCGTTGTGTCTGCTGTGAAATCCCCGGGCTCAACCTGGGAATGGCAGTGGAAACTGGCGAGCTTGAGTGTGGCAGA\
GGGGGGGGGAATTCCGCGTGTAGCAGTGAAATGCGTAGAGATGCGGAGGAACACCGATGGCGAAGGCAACCCCCTGGG\
ATAATATTTACGCTCAT
+
AABCCFFFFFFFGGGGGGGGGGGGHHHHHHHEGGFG2EEGGGGGGHHGGGGGHGHHHHHHGGGGHHHGGGGGGGGGGG\
GEGGGGHEG?GBGGFHFGFFHHGHHHGGGGCCHHHHHFCGG01GGHGHGGGEFHH/DDHFCCGCGHGAF;B0;DGF9A\
EEGGGF-=C;.FFFF/.-@B9BFB/BB/;BFBB/..9=.9//:/:@---./.BBD-@CFD/=A-::.9AFFFFFCEFF\
./FBB############
@MISEQ03:123:000000000-A40KM:1:1101:14875:1613 1:N:0:TCCACAGGAGT
GGTGGGTGCCAGCCGCCGCGGTAATACAGAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCGCGTAGGTG\
GTTTGTTAAGTTGGATGTGAAATCCCCGGGCTCAACCTGGGAACTGCATTCAAAACTGACAAGCTAGAGTATGGTAGA\
GGGTGGTGGAATTTCCTGTGTAGCGGTGAAATGCGTAGATATAGGAAGGAACACCAGTGGCGAAGGCGACCACCTGGA\
CTGATACTGACACTGAG
+
CCCCCCCCCFFFGGGGGGGGGGGGGHHHHHHGGGGGHHHHGGGGGHHGGGGGHHHHHHHHGGGGHHHGGGGGGGGGHH\
GGGHGHHHHHHHHHHHHHHHHHGHHHGGGGGGHHHHHHHHGGHGGHHGHHHHHHHHHFHHHHHHHHHHHHHHHGHHHG\
HHGEGGDGGFFFGGGFGGGGGGGGGGGFFFFFFFDFFFAFFFFFFFFFFFFFFFFFFFFFFFFFFDFFFFFFFEFFFF\
FFFFFB:FFFFFFFFFF
"""

BARCODES = """@MISEQ03:123:000000000-A40KM:1:1101:14149:1572 1:N:0:TCCACAGGAGT
TCCACAGGAGT
+
CCCCCCCCCFF
@MISEQ03:123:000000000-A40KM:1:1101:14170:1596 1:N:0:TCCACAGGAGT
TCCACAGGAGT
+
CCCCCCCCCCF
@MISEQ03:123:000000000-A40KM:1:1101:14740:1607 1:N:0:TCCACAGGAGT
TCCACAGGAGT
+
AABCCFFFFFF
@MISEQ03:123:000000000-A40KM:1:1101:14875:1613 1:N:0:TCCACAGGAGT
TCCACAGGAGT
+
CCCCCCCCCFF
"""

EXP_HTML = """<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>filename</th>
      <th>md5</th>
      <th>file_type</th>
      <th>reads</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>1_s_G1_L001_sequences_barcodes.fastq.gz</td>
      <td>37fbe374c834e6fb0cb96f986dd3aecb</td>
      <td>raw_barcodes</td>
      <td>4</td>
    </tr>
    <tr>
      <td>1_s_G1_L001_sequences.fastq.gz</td>
      <td>66942f7322839dab50eec931f3926a2d</td>
      <td>raw_forward_seqs</td>
      <td>4</td>
    </tr>
  </tbody>
</table>"""

EXP_HTML_DEMUX_REGEXP = [
    '<h3>Features</h3>',
    '<b>Total</b>: 49', '<br/>',
    '<b>Max</b>: 151', '<br/>',
    '<b>Mean</b>: 151', '<br/>',
    '<b>Standard deviation</b>: 151', '<br/>',
    '<b>Median</b>: 0', '<br/>',
    '<img src = "data:image/png;base64,.*"/>']

if __name__ == '__main__':
    main()
