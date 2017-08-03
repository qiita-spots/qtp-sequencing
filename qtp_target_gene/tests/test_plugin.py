# -----------------------------------------------------------------------------
# Copyright (c) 2014--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------

from unittest import main
from os import remove
from os.path import exists, isdir, join
from shutil import rmtree
from tempfile import mkdtemp
from json import dumps
from gzip import GzipFile
from time import sleep

from qiita_client.testing import PluginTestCase

from qtp_target_gene import plugin


class PluginTests(PluginTestCase):
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

    def _wait_job(self, job_id):
        for i in range(10):
            status = self.qclient.get_job_info(job_id)['status']
            if status != 'running':
                break
            sleep(1)
        return status

    def test_plugin_summary(self):
        artifact_id = 1
        data = {'command': dumps(['Target Gene type', '0.1.0',
                                  'Generate HTML summary']),
                'parameters': dumps({'input_data': artifact_id}),
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
            fh.write(BARCODES)
        fwd_fp = files['raw_forward_seqs'][0]
        self._clean_up_files.append(fwd_fp)
        with GzipFile(fwd_fp, mode='w', mtime=1) as fh:
            fh.write(READS)

        plugin("https://localhost:21174", job_id, self.out_dir)
        self._wait_job(job_id)
        obs = self.qclient.get_job_info(job_id)
        self.assertEqual(obs['status'], 'success')

    def test_plugin_validate(self):
        fp = join(self.out_dir, 'prefix1.fastq')
        with open(fp, 'w') as f:
            f.write(READS)
        fp2 = join(self.out_dir, 'prefix1_b.fastq')
        with open(fp2, 'w') as f:
            f.write(BARCODES)
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
                      'files': dumps(files),
                      'artifact_type': atype}

        data = {'command': dumps(['Target Gene type', '0.1.0', 'Validate']),
                'parameters': dumps(parameters),
                'status': 'running'}
        job_id = self.qclient.post(
            '/apitest/processing_job/', data=data)['job']

        plugin("https://localhost:21174", job_id, self.out_dir)
        self._wait_job(job_id)
        obs = self.qclient.get_job_info(job_id)
        self.assertEqual(obs['status'], 'success')

    def test_plugin_error(self):
        parameters = {'template': 1,
                      'files': dumps({'log': ['/path/to/file1.log']}),
                      'artifact_type': "Demultiplexed"}
        data = {'command': dumps(['Target Gene type', '0.1.0', 'Validate']),
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
