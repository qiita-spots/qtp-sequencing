# -----------------------------------------------------------------------------
# Copyright (c) 2014--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------

from os.path import join, basename
from tempfile import mkdtemp
from itertools import zip_longest
from qiita_client.util import system_call, get_sample_names_by_run_prefix


def make_read_pairs_per_sample(forward_seqs, reverse_seqs, map_file):
    """Recovers read pairing information

    Parameters
    ----------
    forward_seqs : list of str
        The list of forward seqs filepaths
    reverse_seqs : list of str
        The list of reverse seqs filepaths
    map_file : str
        The path to the mapping file

    Returns
    -------
    samples: list of tup
        list of 3-tuples with run prefix, sample name, fwd read fp, rev read fp

    Raises
    ------
    ValueError
        If the rev is not an empty list and the same length than fwd seqs
        The prefixes of the run_prefix don't match the file names

    Notes
    -----
    At this stage it is required that if reverse sequences are present that all
    samples have both a forward and a reverse sequence. However, the read
    trimming step can sometimes eliminate all reverse reads, especially in low
    coverage samples with poor overall reverse read quality.

    Note that this functionality should be merged into a more central location
    in Qiita, perhaps along with a refactored and more useful fastq artifact
    structure (see https://github.com/qiita-spots/qtp-target-gene/issues/8)

    Currently this code is being borrowed form the kneaddata plugin pull req
    """

    # sort forward seqs
    forward_seqs.sort()

    # check that rev seqs are same len
    if reverse_seqs:
        if len(forward_seqs) != len(reverse_seqs):
            raise ValueError('Your reverse and forward files are of different '
                             'length. Forward: %s. Reverse: %s.' %
                             (', '.join(forward_seqs),
                              ', '.join(reverse_seqs)))
        reverse_seqs.sort()

    # get run prefixes
    # These are prefixes that should match uniquely to forward reads
    # sn_by_rp is dict of samples keyed by run prefixes
    sn_by_rp = get_sample_names_by_run_prefix(map_file)

    # make pairings
    samples = []
    used_prefixes = set()
    for i, (fwd_fp, rev_fp) in enumerate(zip_longest(forward_seqs,
                                                     reverse_seqs)):
        # fwd_fp is the fwd read filepath
        fwd_fn = basename(fwd_fp)

        # iterate over run prefixes and make sure only one matches
        run_prefix = None
        for rp in sn_by_rp:
            if fwd_fn.startswith(rp) and run_prefix is None:
                run_prefix = rp
            elif fwd_fn.startswith(rp) and run_prefix is not None:
                raise ValueError('Multiple run prefixes match this fwd read: '
                                 '%s' % fwd_fn)

        # make sure that we got one matching run prefix:
        if run_prefix is None:
            raise ValueError('No run prefix matching this fwd read: %s'
                             % fwd_fn)

        if run_prefix in used_prefixes:
            raise ValueError('This run prefix matches multiple fwd reads: '
                             '%s' % run_prefix)

        if rev_fp is None:
            samples.append((run_prefix, sn_by_rp[run_prefix], fwd_fp, None))
        else:
            rev_fn = basename(rev_fp)
            # if we have reverse reads, make sure the matching pair also
            # matches the run prefix:
            if not rev_fn.startswith(run_prefix):
                raise ValueError('Reverse read does not match this run prefix.'
                                 '\nRun prefix: %s\nForward read: %s\n'
                                 'Reverse read: %s\n' %
                                 (run_prefix, fwd_fn, rev_fn))

            samples.append((run_prefix, sn_by_rp[run_prefix], fwd_fp,
                            rev_fp))

        used_prefixes.add(run_prefix)

    return(samples)


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


def fastqc(qclient, job_id, map_file, filepaths):
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
    tar_commands = []

    tar_commands.append("tar -czvf %s %s" % (join(out_dir, 'fastqc.tar.gz'),
                        fqc_out_dir))
    tar_commands.append("tar -czvf %s %s" % (join(out_dir, 'multiqc.tar.gz'),
                        mqc_out_dir))

    msg = "Step 5 of 6: Compressing FastQC and MultiQC output"
    success, msg = _run_commands(qclient, job_id, tar_commands, msg)
    if not success:
        return False, None, msg

    return True, [(join(out_dir, 'multiqc.tar.gz'), 'tgz'),
                  (join(out_dir, 'fastqc.tar.gz'), 'tgz')], msg
