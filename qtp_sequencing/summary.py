# -----------------------------------------------------------------------------
# Copyright (c) 2014--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------

from hashlib import md5
from os.path import basename, join, dirname, exists
from base64 import b64encode
from io import BytesIO

from qiita_files.demux import stats as demux_stats
from qiita_client.util import system_call

import matplotlib
import pandas as pd

matplotlib.use('Agg')

import matplotlib.pyplot as plt # noqa

FILEPATH_TYPE_NO_FQTOOLS = ['SFF', 'FASTA_preprocessed']


def _generate_html_summary(artifact_type, filepaths, out_dir):
    """Helper method to generate html_summary

    Parameters
    ----------
    artifact_type : str
        The artifact_type to summarize
    filepaths : [(str, str)]
        A list of string pairs where the first element is the filepath and the
        second is the filepath type
    out_dir : str
        The output folder

    Returns
    -------
    str
        The str representing the artifact html summary
    """
    # we have 2 main cases: Demultiplexed and everything else,
    # splitting on those
    if artifact_type == 'Demultiplexed':
        artifact_information = '\n'.join(_summary_demultiplexed(
            artifact_type, filepaths))
        if artifact_information is None:
            raise ValueError("We couldn't find a demux file in your artifact")
    elif artifact_type == 'FASTA_preprocessed':
        artifact_information = '\n'.join(_summary_FASTA_preprocessed(
            artifact_type, filepaths, out_dir))
    else:
        artifact_information = _summary_not_demultiplexed(
            artifact_type, filepaths)

    return artifact_information


def generate_html_summary(qclient, job_id, parameters, out_dir):
    """Generates the HTML summary of a sequencing type artifact

    Parameters
    ----------
    qclient : qiita_client.QiitaClient
        The Qiita server client
    job_id : str
        The job id
    parameters : dict
        The parameter values to validate and create the artifact
    out_dir : str
        The path to the job's output directory

    Returns
    -------
    bool, None, str
        Whether the job is successful
        Ignored
        The error message, if not successful

    Raises
    ------
    ValueError
        - If there is any error gathering the information from the server
        - If the artifact is 'Demultiplexed' but it doesn't have a demux file
    """
    # Step 1: gather file information from qiita using REST api
    artifact_id = parameters['input_data']
    qclient_url = "/qiita_db/artifacts/%s/" % artifact_id
    artifact_info = qclient.get(qclient_url)

    # 1a. getting the file paths
    filepaths = artifact_info['files']
    # 1.b get the artifact type_info
    artifact_type = artifact_info['type']

    artifact_information = _generate_html_summary(
        artifact_type, filepaths, out_dir)

    of_fp = join(out_dir, "artifact_%d.html" % artifact_id)
    with open(of_fp, 'w') as of:
        of.write(artifact_information)

    # Step 3: add the new file to the artifact using REST api
    success = True
    error_msg = ''
    try:
        qclient.patch(qclient_url, 'add', '/html_summary/', value=of_fp)
    except Exception as e:
        success = False
        error_msg = str(e)

    return success, None, error_msg


def _summary_not_demultiplexed(artifact_type, filepaths):
    """Generates the HTML summary for non Demultiplexed artifacts

    Parameters
    ----------
    artifact_type : str
        The artifact type
    filepaths : [(str, str)]
        A list of string pairs where the first element is the filepath and the
        second is the filepath type

    Returns
    -------
    list
        A list of strings with the html summary
    """
    # loop over each of the fps/fps_type pairs
    artifact_information = []
    errors = []
    df = None
    for fps_type, fps in sorted(filepaths.items()):
        if fps_type in {'html_summary', 'log'}:
            continue
        # Step 2: generate HTML summary
        # md5, from http://stackoverflow.com/a/3431838
        for i, fp in enumerate(fps):
            fn = basename(fp)
            with open(fp, "rb") as f:
                hash_md5 = md5()
                for chunk in iter(lambda: f.read(4096), b""):
                    hash_md5.update(chunk)
            data = {'filename': fn, 'md5': hash_md5.hexdigest(),
                    'file_type': fps_type}

            if artifact_type not in FILEPATH_TYPE_NO_FQTOOLS:
                # check if the validate summary is present
                if i == 0:
                    fdata = f'{dirname(fp)}/qtp-sequencing-validate-data.csv'
                    if exists(fdata):
                        df = pd.read_csv(fdata, index_col=None)

                if df is None:
                    cmd = f'fqtools count {fp}'
                    std_out, std_err, return_value = system_call(cmd)
                    if std_err or return_value != 0:
                        errors.append(f'{fn}: {std_err}')
                        reads = None
                    else:
                        reads = int(std_out)
                else:
                    reads = df[(df.filename == fn) &
                               (df.file_type == fps_type)]
                    # [0] there is only one value
                    reads = reads.reads.values[0]
                data['reads'] = reads

            artifact_information.append(data)

    if errors:
        raise ValueError('Found errors: \n %s' % ''.join(errors))

    df = pd.DataFrame(artifact_information)
    order = ['file_type', 'reads'] if 'reads' in df.columns else ['file_type']
    df.sort_values(order, inplace=True)

    return df.to_html(index=False)


def _summary_demultiplexed(artifact_type, filepaths):
    """Generates the HTML summary for Demultiplexed artifacts

    Parameters
    ----------
    artifact_type : str
        The artifact type
    filepaths : [(str, str)]
        A list of string pairs where the first element is the filepath and the
        second is the filepath type

    Returns
    -------
    list
        A list of strings with the html summary
    """
    # loop over each of the fps/fps_type pairs to find the demux_fp
    demux_fps = filepaths.get('preprocessed_demux', None)
    if demux_fps is None:
        return None

    # If demux_fps exists, it should contain only a single file
    demux_fp = demux_fps[0]

    # generating html summary
    artifact_information = []
    sn, smax, smin, smean, sstd, smedian, shist, shist_edge = demux_stats(
        demux_fp)
    artifact_information.append("<h3>Features</h3>")
    artifact_information.append('<b>Total</b>: %d' % sn)
    artifact_information.append("<br/>")
    artifact_information.append('<b>Max</b>: %d' % smax)
    artifact_information.append("<br/>")
    artifact_information.append('<b>Mean</b>: %d' % smean)
    artifact_information.append("<br/>")
    artifact_information.append('<b>Standard deviation</b>: %d' % sstd)
    artifact_information.append("<br/>")
    artifact_information.append('<b>Median</b>: %d' % smedian)
    artifact_information.append("<br/>")

    # taken from http://stackoverflow.com/a/9141911
    plt.bar(shist_edge[:-1], shist, width=1)
    plt.xlim(min(shist_edge), max(shist_edge))
    plt.xlabel('Sequence Length')
    plt.ylabel('Number of sequences')
    plot = BytesIO()
    plt.savefig(plot, format='png')
    plot.seek(0)
    artifact_information.append(
        '<img src = "data:image/png;base64,{}"/>'.format(
            b64encode(plot.getvalue()).decode('utf-8')))

    return artifact_information


def _summary_FASTA_preprocessed(artifact_type, filepaths, out_dir):
    """Generates the HTML summary for Demultiplexed artifacts

    Parameters
    ----------
    artifact_type : str
        The artifact type
    filepaths : [(str, str)]
        A list of string pairs where the first element is the filepath and the
        second is the filepath type
    out_dir : str
        The output folder

    Returns
    -------
    list
        A list of strings with the html summary
    """
    files = filepaths.get('preprocessed_fasta')
    cmd = f"quast %s -o {out_dir}/quast" % ' '.join(files)
    std_out, std_err, return_value = system_call(cmd)
    if return_value != 0:
        artifact_information = (
            "Std out: %s\nStd err: %s\n\nCommand run was:\n%s" % (
                std_out, std_err, cmd))
    else:
        with open(f'{out_dir}/quast/report.html', 'r') as f:
            artifact_information = f.readlines()

    return artifact_information
