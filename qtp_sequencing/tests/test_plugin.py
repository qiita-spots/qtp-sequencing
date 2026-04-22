# -----------------------------------------------------------------------------
# Copyright (c) 2014--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------

from unittest import main
from os import remove, makedirs, close, write
from os.path import exists, isdir, dirname
from shutil import rmtree
from tempfile import mkdtemp, mkstemp
from json import dumps
from gzip import GzipFile
from time import sleep

from qiita_client.testing import PluginTestCase

from qtp_sequencing import plugin


class PluginTests(PluginTestCase):
    def setUp(self):
        self.out_dir = mkdtemp()
        self._clean_up_files = [self.out_dir]
        makedirs(self.base_data_dir, exist_ok=True)

    def tearDown(self):
        for fp in self._clean_up_files:
            if exists(fp):
                if isdir(fp):
                    rmtree(fp)
                else:
                    remove(fp)

    def _wait_job(self, job_id):
        for i in range(10):
            status = self.qclient.get_job_info(job_id)['status']
            if status != 'running':
                break
            sleep(1)
        return status

    def test_plugin_summary(self):
        artifact_id = 1
        data = {'command': dumps(['Sequencing Data Type', '2022.11',
                                  'Generate HTML summary']),
                'parameters': dumps({'input_data': artifact_id}),
                'status': 'running'}
        job_id = self.qclient.post(
            '/apitest/processing_job/', data=data)['job']
        # Qiita will return a filepath, but in the test environment, these
        # files do not exist - create them
        files = self.qclient.get(
            '/qiita_db/artifacts/%s/' % artifact_id,
            no_file_fetching=True)['files']

        bcds_fp = files['raw_barcodes'][0]['filepath']
        self._clean_up_files.append(bcds_fp)
        makedirs(dirname(bcds_fp), exist_ok=True)
        with GzipFile(bcds_fp, mode='w', mtime=1) as fh:
            fh.write(BARCODES.encode())
        self.qclient.push_file_to_central(bcds_fp)

        fwd_fp = files['raw_forward_seqs'][0]['filepath']
        self._clean_up_files.append(fwd_fp)
        makedirs(dirname(fwd_fp), exist_ok=True)
        with GzipFile(fwd_fp, mode='w', mtime=1) as fh:
            fh.write(READS.encode())
        self.qclient.push_file_to_central(fwd_fp)

        plugin("https://localhost:21174", job_id, self.out_dir)
        self._wait_job(job_id)
        obs = self.qclient.get_job_info(job_id)
        self.assertEqual(obs['status'], 'success')

    def test_plugin_validate(self):
        f, fp = mkstemp(suffix="prefix1.fastq", dir=self.base_data_dir)
        write(f, READS.encode("utf-8"))
        close(f)
        self.qclient.push_file_to_central(fp)

        f, fp2 = mkstemp(suffix="prefix1_b.fastq", dir=self.base_data_dir)
        write(f, BARCODES.encode("utf-8"))
        close(f)
        self.qclient.push_file_to_central(fp2)

        prep_info = {"1.SKB2.640194": {"not_a_run_prefix": "prefix1"}}
        files = {'raw_forward_seqs': [fp],
                 'raw_barcodes': [fp2]}
        atype = "FASTQ"
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

        plugin("https://localhost:21174", job_id, self.out_dir)
        self._wait_job(job_id)
        obs = self.qclient.get_job_info(job_id)
        self.assertEqual(obs['status'], 'success')

    def test_plugin_error(self):
        fd, log_fp = mkstemp(suffix=".file1.log", dir=self.base_data_dir)
        close(fd)
        parameters = {
            'template': 1,
            'analysis': None,
            'files': dumps(
                {'log': [self.qclient.push_file_to_central(log_fp)]}),
            'artifact_type': "Demultiplexed"}
        data = {'command': dumps(
            ['Sequencing Data Type', '2022.11', 'Validate']),
                'parameters': dumps(parameters),
                'status': 'running'}
        job_id = self.qclient.post(
            '/apitest/processing_job/', data=data)['job']
        plugin("https://localhost:21174", job_id, self.out_dir)
        self._wait_job(job_id)
        obs = self.qclient.get_job_info(job_id)
        self.assertEqual(obs['status'], 'error')


CONFIG_FILE = """
[main]
SERVER_CERT = %s

# Oauth2 plugin configuration
CLIENT_ID = %s
CLIENT_SECRET = %s
"""

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

if __name__ == '__main__':
    main()
