# -----------------------------------------------------------------------------
# Copyright (c) 2014--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------

from hashlib import md5
from gzip import open as gopen
from os.path import basename, join
from urllib import quote
from base64 import b64encode
from StringIO import StringIO


from qiita_ware.demux import stats as demux_stats
import matplotlib.pyplot as plt


FILEPATH_TYPE_TO_NOT_SHOW_HEAD = ['SFF']
LINES_TO_READ_FOR_HEAD = 10


def generate_html_summary(qclient, job_id, parameters, out_dir):
    """Generates the HTML summary of a target gene type artifact

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

    # we have 2 main cases: Demultiplexed and everything else,
    # splitting on those
    if artifact_type == 'Demultiplexed':
        artifact_information = _summary_demultiplexed(
            artifact_type, filepaths)
        if artifact_information is None:
            raise ValueError("We couldn't find a demux file in your artifact")
    else:
        artifact_information = _summary_not_demultiplexed(
            artifact_type, filepaths)

    of_fp = join(out_dir, "artifact_%d.html" % artifact_id)
    with open(of_fp, 'w') as of:
        of.write('\n'.join(artifact_information))

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
    for fps_type, fps in sorted(filepaths.items()):
        # Step 2: generate HTML summary
        # md5, from http://stackoverflow.com/a/3431838
        for fp in fps:
            with open(fp, "rb") as f:
                hash_md5 = md5()
                for chunk in iter(lambda: f.read(4096), b""):
                    hash_md5.update(chunk)

            # getting head of the files
            header = []
            if artifact_type not in FILEPATH_TYPE_TO_NOT_SHOW_HEAD:
                # we need to encapsulate the full for loop because gzip will
                # not raise an error until you try to read
                try:
                    with gopen(fp, 'r') as fin:
                        header = [
                            next(fin) for x in xrange(LINES_TO_READ_FOR_HEAD)]
                except IOError:
                    with open(fp, 'r') as fin:
                        header = [
                            next(fin) for x in xrange(LINES_TO_READ_FOR_HEAD)]
            filename = basename(fp)
            artifact_information.append(
                "<h3>%s (%s)</h3>" % (filename, fps_type))
            artifact_information.append("<b>MD5:</b>: %s</br>" %
                                        hash_md5.hexdigest())
            if header:
                artifact_information.append(
                    "<p style=\"font-family:'Courier New', Courier, monospace;"
                    "font-size:10;\">%s</p><hr/>" % ("<br/>".join(header)))

    return artifact_information


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
    plot = StringIO()
    plt.savefig(plot, format='png')
    plot.seek(0)
    uri = 'data:image/png;base64,' + quote(b64encode(plot.buf))
    artifact_information.append('<img src = "%s"/>' % uri)

    return artifact_information


def generate_fastqc_commands(forward_seqs, reverse_seqs, map_file, out_dir):
    """Generates the FastQC commands

    Parameters
    ----------
    seqs : list of str
        The list of seqs filepaths
    map_file : str
        The path to the mapping file
    out_dir : str
        The job output directory

    Returns
    -------
    list of str
        The FastQC commands

    Raises
    ------

    Notes
    -----
    This is presently reproducing functionality from the kneaddata pipeline
    in qp_shotgun. 
    """

    samples = make_read_pairs_per_sample(forward_seqs, reverse_seqs, map_file)

    cmds = []
    prefixes = []

    fps = []
    for run_prefix, sample, f_fp, r_fp in samples:
        prefixes.append(run_prefix)
        if r_fp is None:
            cmds.append('mkdir -p %s; fastqc --outdir "%s" %s' %
                        (join(out_dir, run_prefix), join(out_dir, run_prefix),
                         f_fp))
            fps.append((f_fp, None))
        else:
            cmds.append('mkdir -p %s; fastqc --outdir "%s" %s %s' %
                        (join(out_dir, run_prefix), join(out_dir, run_prefix),
                         f_fp, r_fp))
            fps.append((f_fp, r_fp))

    return cmds, samples


def generate_multiqc_commands(file_dir, out_dir):
    """Generates the multiqc commands

    Parameters
    ----------
    file_dir : str
        Filepath to directory of QC info
    out_dir : str
        The job output directory

    Returns
    -------
    list of str
        The multiqc commands

    Notes
    -----
    """

    # make input file list

    cmds = []

    cmds.append('multiqc --outdir %s --filename MultiQC.html %s'
                '%s' % (out_dir, file_dir))

    return cmds


def _guess_fastqc_filename(fp):
    f_p = basename(fp)

    exts = ['.fastq', '.fq', '.gz', '.gzip']
    while splitext(f_p)[1] in exts:
        f_p = splitext(f_p)[0]

    return "%s_fastqc.html" % f_p, "%s_fastqc.zip" % f_p



def _run_commands(qclient, job_id, commands, msg):
    for i, cmd in enumerate(commands):
        qclient.update_job_step(job_id, msg % i)
        std_out, std_err, return_value = system_call(cmd)
        if return_value != 0:
            error_msg = ("Error running FastQC:\nStd out: %s\nStd err: %s"
                         "\n\nCommand run was:\n%s"
                         % (std_out, std_err, cmd))
            return False, error_msg

    return True, ""


def fastqc(qclient, job_id, parameters, out_dir):
    """Run FastQC and MultiQC on a sequence artifact

    Parameters
    ----------
    qclient : tgp.qiita_client.QiitaClient
        The Qiita server client
    job_id : str
        The job id
    parameters : dict
        The parameter values to run split libraries
    out_dir : str
        Yhe path to the job's output directory

    Returns
    -------
    bool, list, str
        The results of the job
    """
    # Step 1 get the rest of the information need to run FastQC
    qclient.update_job_step(job_id, "Step 1 of 3: Collecting information")
    artifact_id = parameters['input']

    # removing input from parameters so it's not part of the final command
    del parameters['input']

    # Get the artifact filepath information
    artifact_info = qclient.get("/qiita_db/artifacts/%s/" % artifact_id)
    fps = artifact_info['files']

    # Get the artifact metadata
    prep_info = qclient.get('/qiita_db/prep_template/%s/'
                            % artifact_info['prep_information'][0])
    qiime_map = prep_info['qiime-map']

    # Generate temp dir for exection of QC steps
    out_dir = mkdtemp()

    # Step 2 generating command FastQC
    qclient.update_job_step(job_id, "Step 2 of 3: Generating"
                                    " FastQC command")

    fqc_out_dir = join(out_dir, 'FastQC')
    
    rs = fps['raw_reverse_seqs'] if 'raw_reverse_seqs' in fps else []
    fqc_cmds, samples = generate_fastqc_commands(fps['raw_forward_seqs'],
                                                 rs, qiime_map, fqc_out_dir,
                                                 parameters)

    # Step 3 execute FastQC
    msg = "Step 3 of 3: Executing FastQC job (%d/{0})".format(len(fqc_cmds))
    success, msg = _run_commands(qclient, job_id, fqc_cmds, msg)
    if not success:
        return False, None, msg

    # Generate MultiQC command
    mqc_out_dir = join(out_dir, 'MultiQC')

    qclient.update_job_step(job_id, "Step 4 of 3: Generating"
                                    " MultiQC command")
    
    mqc_cmds, samples = generate_multiqc_commands(mqc_out_dir,
                                                  parameters)

    # Execute MultiQC
    msg = "Step 5 of 3: Executing MultiQC job (%d/{0})".format(len(mqc_cmds))
    success, msg = _run_commands(qclient, job_id, mqc_cmds, msg)
    if not success:
        return False, None, msg

    # Compress QC files
    commands = []

    fqc_cmd = "tar -czvf %s %s" % (join(out_dir,'fastqc.tar.gz'), fqc_out_dir)
    mqc_cmd = "tar -czvf %s %s" % (join(out_dir,'multiqc.tar.gz'), mqc_out_dir)

    msg = "Step 6 of 3: Compressing FastQC and MultiQC output"
    success, msg = _run_commands(qclient, job_id, commands, msg)
    if not success:
        return False, None, msg
    

    # Step 4 generating artifacts

    ainfo = []

    ainfo.extend([
        ArtifactInfo('MultiQC data', 'tgz',
                     [(join(out_dir,'multiqc.tar.gz'), 'tgz')]),
        ArtifactInfo('FastQC data', 'tgz',
                     [(join(out_dir,'fastqc.tar.gz'), 'tgz')])])


    return True, ainfo, ""
