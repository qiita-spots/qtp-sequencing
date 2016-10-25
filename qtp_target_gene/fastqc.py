# -----------------------------------------------------------------------------
# Copyright (c) 2014--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------


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
            cmds.append('mkdir -p %s; fastqc --outdir "%s" --noextract %s' %
                        (join(out_dir, run_prefix), join(out_dir, run_prefix),
                         f_fp))
            fps.append((f_fp, None))
        else:
            cmds.append('mkdir -p %s; fastqc --outdir "%s" --noextract %s %s' %
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


def fastqc(qclient, job_id, map_file, files):
    """Run FastQC and MultiQC on a sequence artifact

    Parameters
    ----------
    qclient : qiita_client.QiitaClient
        The Qiita server client
    filepaths: list of tup
        List of per_sample_fastq filepath artifacts from validate function
    map_file : str
        The path to the mapping file
    Returns
    -------
    bool, list of tup, str
        The results of the job
    """

    # Get forward and (if exists) reverse seq filepaths
    fwd_fps = [x for x, y in filepaths if y == 'raw_forward_seqs']
    rev_fps = [x for x, y in filepaths if y == 'raw_reverse_seqs']

    # Generate temp dir for exection of QC steps
    out_dir = mkdtemp()

    # Generating commands for FastQC
    qclient.update_job_step(job_id, "Step 1 of 6: Generating"
                                    " FastQC command")

    fqc_out_dir = join(out_dir, 'FastQC')

    fqc_cmds, samples = generate_fastqc_commands(fwd_fps, rev_fps, map_file,
                                                 fqc_out_dir)

    # Execute FastQC
    msg = "Step 2 of 6: Executing FastQC job (%d/{0})".format(len(fqc_cmds))
    success, msg = _run_commands(qclient, job_id, fqc_cmds, msg)
    if not success:
        return False, None, msg

    # Generate MultiQC command
    qclient.update_job_step(job_id, "Step 3 of 6: Generating"
                                    " MultiQC command")

    mqc_out_dir = join(out_dir, 'MultiQC')

    mqc_cmds, samples = generate_multiqc_commands(mqc_out_dir,
                                                  fqc_out_dir)

    # Execute MultiQC
    msg = "Step 4 of 6: Executing MultiQC job (%d/{0})".format(len(mqc_cmds))
    success, msg = _run_commands(qclient, job_id, mqc_cmds, msg)
    if not success:
        return False, None, msg

    # Compress QC files
    commands = []

    fqc_cmd = "tar -czvf %s %s" % (join(out_dir,'fastqc.tar.gz'), fqc_out_dir)
    mqc_cmd = "tar -czvf %s %s" % (join(out_dir,'multiqc.tar.gz'), mqc_out_dir)

    msg = "Step 5 of 6: Compressing FastQC and MultiQC output"
    success, msg = _run_commands(qclient, job_id, commands, msg)
    if not success:
        return False, None, msg

    return True, [(join(out_dir,'multiqc.tar.gz'), 'tgz'),
                  (join(out_dir,'fastqc.tar.gz'), 'tgz')], msg
