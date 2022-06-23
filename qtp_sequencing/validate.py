# -----------------------------------------------------------------------------
# Copyright (c) 2014--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------

from os.path import basename, join, splitext, getsize, dirname
from os import remove
from json import loads
from shutil import copy
from h5py import File
from gzip import open as gopen
from collections import defaultdict

from qiita_client import ArtifactInfo
from qiita_client.util import system_call
from qiita_files.util import open_file
from qiita_files.demux import to_hdf5, to_ascii_file

import pandas as pd

from .summary import FILEPATH_TYPE_NO_FQTOOLS, _generate_html_summary


FILEPATH_TYPE_DICT = {
    'SFF': ({'raw_sff'}, set()),
    'FASTQ': ({'raw_forward_seqs', 'raw_barcodes'}, {'raw_reverse_seqs'}),
    'FASTA': ({'raw_fasta'}, {'raw_qual'}),
    'FASTA_Sanger': ({'raw_fasta'}, set()),
    'FASTA_preprocessed': ({'preprocessed_fasta'}, set()),
}

MUST_GZ = {
    # raw input files: FASTQ, per_sample_FASTQ
    'raw_forward_seqs', 'raw_barcodes', 'raw_reverse_seqs', 'raw_fasta',
    # preprocessed files: demultiplexed, trimmed
    'preprocessed_fastq', 'preprocessed_fasta'}


def _gzip_file(filepath, test=False):
    """gzip the given filepath if needed

    Parameters
    ----------
    filepath : string
        The filepath to verify or compress
    test : bolean
        If True do not compress but change the filename, used for unit testing

    Returns
    -------
    str
        the new gz filepath, None if error
    str
        the error, None if success
    """
    error = None
    return_fp = filepath
    if test:
        return_fp = '%s.gz' % filepath
    else:
        is_gzip = False
        try:
            with gopen(filepath, 'rb') as f:
                f.read(1)
            is_gzip = True
        except (OSError, IOError):
            pass

        if not is_gzip:
            gz_cmd = 'pigz -p 5 -c {0} > {0}.gz'.format(filepath)

            std_out, std_err, return_value = system_call(gz_cmd)
            if return_value != 0 and not test:
                error = ("Std out: %s\nStd err: %s\n\nCommand run was:\n%s"
                         % (std_out, std_err, gz_cmd))
            else:
                # removing non gz file
                remove(filepath)
                return_fp = '%s.gz' % filepath
    return return_fp, error


def _validate_multiple(qclient, job_id, prep_info, files, atype, test=False):
    """Validate and fix a new 'SFF', 'FASTQ', 'FASTA' or 'FASTA_Sanger' artifact

    Parameters
    ----------
    qclient : qiita_client.QiitaClient
        The Qiita server client
    job_id : str
        The job id
    prep_info : dict of {str: dict of {str: str}}
        The prep information keyed by sample id
    files : dict of {str: list of str}
        The files to add to the new artifact, keyed by filepath type
    atype: str
        The type of the artifact
    test: bolean, optional
        If True this is being called by a test

    Returns
    -------
    dict
        The results of the job
    """
    qclient.update_job_step(job_id, "Step 2: Validating '%s' files" % atype)
    req_fp_types, opt_fp_types = FILEPATH_TYPE_DICT[atype]
    all_fp_types = req_fp_types | opt_fp_types

    # Check if there is any filepath type that is not supported
    unsupported_fp_types = set(files) - all_fp_types
    if unsupported_fp_types:
        error_msg = ("Filepath type(s) %s not supported by artifact "
                     "type %s. Supported filepath types: %s"
                     % (', '.join(unsupported_fp_types), atype,
                        ', '.join(sorted(all_fp_types))))
        return False, None, error_msg

    # Check if the run_prefix column is present in the prep info
    offending = {}
    types_seen = set()
    if 'run_prefix' in prep_info[next(iter(prep_info))]:
        # We can potentially have more than one lane in the prep information
        # so check that the provided files are prefixed with the values in
        # the run_prefix column
        run_prefixes = set(v['run_prefix'] for k, v in prep_info.items())
        num_prefixes = len(run_prefixes)

        # Check those filepath types that are required
        for ftype, t_files in files.items():
            # SFF is an special case cause we can have multiple files with
            # the same prefix
            if num_prefixes != len(t_files) and atype != 'SFF':
                offending[ftype] = (
                    "The number of provided files (%d) doesn't match the "
                    "number of run prefix values in the prep info (%d): %s"
                    % (len(t_files), num_prefixes,
                       ', '.join(basename(f) for f in t_files)))
            else:
                rps = []
                fps = []
                for fp in t_files:
                    bn = basename(fp)
                    found = [rp for rp in run_prefixes if bn.startswith(rp)]
                    if found:
                        rps.extend(found)
                    else:
                        fps.append(bn)
                if fps:
                    offending[ftype] = (
                        "The provided files do not match the run prefix "
                        "values in the prep information: %s" % ', '.join(fps))
                else:
                    rps = run_prefixes - set(rps)
                    if rps:
                        offending[ftype] = (
                            "The following run prefixes in the prep "
                            "information file do not match any file: %s"
                            % ', '.join(rps))

            types_seen.add(ftype)
    else:
        # If the run prefix column is not provided, we only allow a single
        # lane, so check that we have a single file for each provided
        # filepath type
        for ftype, t_files in files.items():
            if len(t_files) != 1:
                offending[ftype] = (
                    "Only one file per type is allowed. Please provide the "
                    "column 'run_prefix' if you need more than one file per "
                    "type: %s" % ', '.join(basename(fp) for fp in t_files))

            types_seen.add(ftype)

    # Check that all required filepath types where present
    missing = req_fp_types - types_seen
    if missing:
        error_msg = ("Missing required filepath type(s): %s"
                     % ', '.join(missing))
        return False, None, error_msg

    # Check if there was any offending file
    if offending:
        error_list = ["%s: %s" % (k, v) for k, v in offending.items()]
        error_msg = ("Error creating artifact. Offending files:\n%s"
                     % '\n'.join(error_list))
        return False, None, error_msg

    # Everything is ok
    filepaths = []
    for fps_type, fps in files.items():
        for fp in fps:
            if fps_type in MUST_GZ:
                fp, error_msg = _gzip_file(fp, test)
                if error_msg is not None:
                    return False, None, error_msg
            filepaths.append((fp, fps_type))

    # let's count sequences; this is basically the last check
    errors = []
    artifact_information = []
    if atype not in FILEPATH_TYPE_NO_FQTOOLS:
        for fp, fpt in filepaths:
            cmd = f'fqtools count {fp}'
            std_out, std_err, return_value = system_call(cmd)
            fn = basename(fp)
            if std_err or return_value != 0:
                errors.append(f'{fn}: {std_err}')
            else:
                reads = int(std_out)
                artifact_information.append(
                    {'filename': fn, 'reads': reads, 'file_type': fpt})

        if errors:
            raise ValueError('Found errors: \n %s' % ''.join(errors))
        dname = dirname(fp)
        pd.DataFrame(artifact_information).to_csv(
            f'{dname}/qtp-sequencing-validate-data.csv', index=False)

    return True, [ArtifactInfo(None, atype, filepaths)], ""


def _validate_per_sample_FASTQ(qclient, job_id, prep_info, files, test=False):
    """Validate and fix a new 'per_sample_FASTQ' artifact

    Parameters
    ----------
    qclient : qiita_client.QiitaClient
        The Qiita server client
    job_id : str
        The job id
    prep_info : dict of {str: dict of {str: str}}
        The prep information keyed by sample id
    files : dict of {str: list of str}
        The files to add to the new artifact, keyed by filepath type
    test: bolean, optional
        If True this is being called by a test

    Returns
    -------
    dict
        The results of the job
    """
    qclient.update_job_step(
        job_id, "Step 2: Validating 'per_sample_FASTQ' files")

    samples = list(prep_info.keys())
    samples_count = len(samples)

    # Check if there is any filepath type that is not supported
    unsupported_fp_types = set(files) - {'raw_forward_seqs',
                                         'raw_reverse_seqs',
                                         'preprocessed_fastq'}
    if unsupported_fp_types:
        error_msg = ("Filepath type(s) %s not supported by artifact "
                     "type per_sample_FASTQ. Supported filepath types: "
                     "raw_forward_seqs, raw_reverse_seqs, preprocessed_fastq"
                     % ', '.join(unsupported_fp_types))
        return False, None, error_msg

    if 'raw_forward_seqs' in files:
        if 'preprocessed_fastq' in files:
            error_msg = ("If raw_forward_seqs is provided, preprocessed_fastq "
                         "should not be provided")
            return False, None, error_msg
        read_files = files['raw_forward_seqs']
        read_files_count = len(read_files)
        counts_match = read_files_count == samples_count
    elif 'preprocessed_fastq' in files:
        if 'raw_reverse_seqs' in files:
            error_msg = ("If preprocessed_fastq is provided, raw_reverse_seqs "
                         "should not be provided")
            return False, None, error_msg
        read_files = files['preprocessed_fastq']
        read_files_count = len(read_files)
        # In the preprocessed_fastq case, we either have 1 file per sample
        # or 4 files per sample
        counts_match = ((read_files_count == samples_count) or
                        (read_files_count == 4 * samples_count))
    else:
        error_msg = ("Missing required filepath type: raw_forward_seqs or "
                     "preprocessed_fastq")
        return False, None, error_msg

    # Make sure that we hve the same number of files than samples
    if 'raw_reverse_seqs' in files:
        rev_count = len(files['raw_reverse_seqs'])
        counts_match = counts_match and (rev_count == samples_count)
    else:
        rev_count = 0

    if not counts_match:
        error_msg = ("The number of provided files doesn't match the "
                     "number of samples (%d): %d raw_forward_seqs, "
                     "%d raw_reverse_seqs (optional, 0 is ok)"
                     % (samples_count, read_files_count, rev_count))
        return False, None, error_msg

    def _check_files(run_prefixes, read_files, rev_count, files):
        # Check that the provided files match the run prefixes
        fwd_fail = [basename(fp) for fp in read_files
                    if not basename(fp).startswith(tuple(run_prefixes))]
        if rev_count > 0:
            rev_fail = [basename(fp) for fp in files['raw_reverse_seqs']
                        if not basename(fp).startswith(tuple(run_prefixes))]
        else:
            rev_fail = []
        return fwd_fail, rev_fail

    # first let's check via sample sample_names
    run_prefixes = [sid.split('.', 1)[1] for sid in samples]
    fwd_fail, rev_fail = _check_files(run_prefixes, read_files,
                                      rev_count, files)

    # if that doesn't work, let's test via run_prefix
    run_prefix_present = 'run_prefix' in prep_info[samples[0]]
    if (fwd_fail or rev_fail) and run_prefix_present:
        run_prefixes = [v['run_prefix'] for k, v in prep_info.items()]
        if samples_count != len(set(run_prefixes)):
            repeated = ["%s (%d)" % (p, run_prefixes.count(p))
                        for p in set(run_prefixes)
                        if run_prefixes.count(p) > 1]
            error_msg = ("The values for the column 'run_prefix' are not "
                         "unique for each sample. Repeated values: %s"
                         % ', '.join(repeated))
            return False, None, error_msg

        fwd_fail, rev_fail = _check_files(run_prefixes, read_files,
                                          rev_count, files)

    if fwd_fail or rev_fail:
        error_msg = "The provided files are not prefixed by sample id"
        if run_prefix_present:
            error_msg += (" or do not match the run prefix values in the "
                          "prep information.")
        else:
            error_msg += "."
        error_msg += (" Offending files:\n raw_forward_seqs: %s\n"
                      "raw_reverse_seqs: %s" % (', '.join(fwd_fail),
                                                ', '.join(rev_fail)))
        return False, None, error_msg

    filepaths = []
    empty_files = []
    for fps_type, fps in files.items():
        for fp in fps:
            try:
                fp_size = getsize(fp)
            except OSError:
                fp_size = 0
            # 62 is the size of a gzip empty files that we generate
            if fp_size <= 100:
                empty_files.append(basename(fp))

            if fps_type in MUST_GZ:
                fp, error_msg = _gzip_file(fp, test)
                if error_msg is not None:
                    return False, None, error_msg

            filepaths.append((fp, fps_type))

    if empty_files:
        error_msg = "Some of the files are empty: %s" % ', '.join(empty_files)
        return False, None, error_msg

    return True, [ArtifactInfo(None, 'per_sample_FASTQ', filepaths)], ""


def _validate_demux_file(qclient, job_id, prep_info, out_dir, demux_fp,
                         fastq_fp=None, fasta_fp=None, log_fp=None):
    """Validate and fix a 'demux' file and regenerate fastq and fasta files

    Parameters
    ----------
    qclient : qiita_client.QiitaClient
        The Qiita server client
    job_id : str
        The job id
    prep_info : dict of {str: dict of {str: str}}
        The prep information keyed by sample id
    out_dir : str
        The output directory
    demux_fp : str
        The demux file path
    fastq_fp : str, optional
        The original fastq filepath. If demux is correct, it will not be
        regenerated
    fasta_fp : str, optional
        The original fasta filepath. If demux is correct, it will no be
        regenerated
    log_fp : str, optional
        The original log filepath

    Returns
    -------
    dict
        The results og the job
    """
    pt_sample_ids = set(prep_info)
    with open_file(demux_fp) as f:
        demux_sample_ids = set(f.keys())

    if not pt_sample_ids.issuperset(demux_sample_ids):
        # The demux sample ids are different from the ones in the prep template
        qclient.update_job_step(job_id, "Step 3: Fixing sample ids")
        # Atempt 1: the user provided the run prefix column - in this case the
        # run prefix column holds the sample ids present in the demux file
        if 'run_prefix' in prep_info[next(iter(pt_sample_ids))]:
            id_map = {v['run_prefix']: k for k, v in prep_info.items()}
            if not set(id_map).issuperset(demux_sample_ids):
                error_msg = ('The sample ids in the "run_prefix" columns '
                             'from the prep information do not match the '
                             'ones in the demux file. Please, correct the '
                             'column "run_prefix" in the prep information to '
                             'map the existing sample ids to the prep '
                             'information sample ids.')
                return False, None, error_msg
        else:
            # Attempt 2: the sample ids in the demux table are the same that
            # in the prep template but without the prefix
            prefix = next(iter(pt_sample_ids)).split('.', 1)[0]
            prefixed = set("%s.%s" % (prefix, s) for s in demux_sample_ids)
            if pt_sample_ids.issuperset(prefixed):
                id_map = {s: "%s.%s" % (prefix, s) for s in demux_sample_ids}
            else:
                # There is nothing we can do. The samples in the demux file do
                # not match the ones in the prep template and we can't fix it
                error_msg = ('The sample ids in the demultiplexed files do '
                             'not match the ones in the prep information. '
                             'Please, provide the column "run_prefix" in '
                             'the prep information to map the existing sample'
                             ' ids to the prep information sample ids.')
                return False, None, error_msg

        # Fix the sample ids
        # Do not modify the original demux file, copy it to a new location
        new_demux_fp = join(out_dir, basename(demux_fp))
        # this if is important so we don't regenerate the demux file if the
        # user uploads fastq or fna
        if demux_fp != new_demux_fp:
            copy(demux_fp, new_demux_fp)
            demux_fp = new_demux_fp

        with open_file(demux_fp, 'r+') as f:
            for old in f:
                f.move(old, id_map[old])
        # When we fix, we always generate the FASTQ and FASTA file
        # By setting them to None, below will be generated
        fastq_fp = None
        fasta_fp = None

    # If we didn't fix anything, we only generate the files if they don't
    # already exists
    name = splitext(basename(demux_fp))[0]
    if not fastq_fp:
        fastq_fp = join(out_dir, "%s.fastq" % name)
        to_ascii_file(demux_fp, fastq_fp, out_format='fastq')
        fastq_fp, error_msg = _gzip_file(fastq_fp)
        if error_msg is not None:
            return False, None, error_msg

    if not fasta_fp:
        fasta_fp = join(out_dir, "%s.fasta" % name)
        to_ascii_file(demux_fp, fasta_fp, out_format='fasta')
        fasta_fp, error_msg = _gzip_file(fasta_fp)
        if error_msg is not None:
            return False, None, error_msg

    filepaths = [(fastq_fp, 'preprocessed_fastq'),
                 (fasta_fp, 'preprocessed_fasta'),
                 (demux_fp, 'preprocessed_demux')]
    if log_fp:
        filepaths.append((log_fp, 'log'))
    return True, [ArtifactInfo(None, 'Demultiplexed', filepaths)], ""


def _validate_demultiplexed(qclient, job_id, prep_info, files, out_dir):
    """Validate and fix a new 'Demultiplexed' artifact

    Parameters
    ----------
    qclient : qiita_client.QiitaClient
        The Qiita server client
    job_id : str
        The job id
    prep_info : dict of {str: dict of {str: str}}
        The prep information keyed by sample id
    files : dict of {str: list of str}
        The files to add to the new artifact, keyed by filepath type
    out_dir : str
        The output directory

    Returns
    -------
    dict
        The results of the job
    """
    qclient.update_job_step(job_id, "Step 2: Validating 'Demultiplexed' files")

    supported_fp_types = {'preprocessed_fasta', 'preprocessed_fastq',
                          'preprocessed_demux', 'log'}
    unsupported_fp_types = set(files) - supported_fp_types
    if unsupported_fp_types:
        error_msg = ("Filepath type(s) %s not supported by artifact type "
                     "Demultiplexed. Supported filepath types: %s"
                     % (', '.join(unsupported_fp_types),
                        ', '.join(sorted(supported_fp_types))))
        return False, None, error_msg

    # At most one file of each type can be provided
    offending = set(fp_t for fp_t, fps in files.items() if len(fps) > 1)
    if offending:
        errors = ["%s (%d): %s"
                  % (fp_t, len(files[fp_t]), ', '.join(files[fp_t]))
                  for fp_t in sorted(offending)]
        error_msg = ("Only one filepath of each file type is supported, "
                     "offending types:\n%s" % "; ".join(errors))
        return False, None, error_msg

    # Check which files we have available:
    fasta = (files['preprocessed_fasta'][0]
             if 'preprocessed_fasta' in files else None)
    fastq = (files['preprocessed_fastq'][0]
             if 'preprocessed_fastq' in files else None)
    demux = (files['preprocessed_demux'][0]
             if 'preprocessed_demux' in files else None)
    log = (files['log'][0] if 'log' in files else None)
    if demux:
        # If demux is available, use that one to perform the validation and
        # generate the fasta and fastq from it
        success, a_info, error_msg = _validate_demux_file(
            qclient, job_id, prep_info, out_dir, demux, log_fp=log)
    elif fastq:
        # Generate the demux file from the fastq
        demux = join(out_dir, "%s.demux" % splitext(basename(fastq))[0])
        with File(demux, 'w') as f:
            # to_hdf5 expects a list
            to_hdf5([fastq], f)
        # Validate the demux, providing the original fastq
        success, a_info, error_msg = _validate_demux_file(
            qclient, job_id, prep_info, out_dir, demux, fastq_fp=fastq,
            log_fp=log)
    elif fasta:
        # Generate the demux file from the fasta
        demux = join(out_dir, "%s.demux" % splitext(basename(fasta))[0])
        with File(demux, 'w') as f:
            # to_hdf5 expects a list
            to_hdf5([fasta], f)
        # Validate the demux, providing the original fasta
        success, a_info, error_msg = _validate_demux_file(
            qclient, job_id, prep_info, out_dir, demux, fasta_fp=fasta,
            log_fp=log)
    else:
        error_msg = ("Either a 'preprocessed_demux', 'preprocessed_fastq' or "
                     "'preprocessed_fasta' file should be provided.")
        return False, None, error_msg

    return success, a_info, error_msg


def validate(qclient, job_id, parameters, out_dir):
    """Validae and fix a new artifact

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
    dict
        The results of the job

    Raises
    ------
    ValueError
        If there is any error gathering the information from the server
    """
    prep_id = parameters['template']
    files = loads(parameters['files'])
    a_type = parameters['artifact_type']

    qclient.update_job_step(job_id, "Step 1: Collecting prep information")
    prep_info = qclient.get("/qiita_db/prep_template/%s/data/" % prep_id)
    prep_info = prep_info['data']

    _vm = ['SFF', 'FASTQ', 'FASTA', 'FASTA_Sanger', 'FASTA_preprocessed']
    if a_type in _vm:
        reply = _validate_multiple(qclient, job_id, prep_info, files, a_type)
    elif a_type == 'per_sample_FASTQ':
        reply = _validate_per_sample_FASTQ(qclient, job_id, prep_info, files)
    elif a_type == 'Demultiplexed':
        reply = _validate_demultiplexed(qclient, job_id, prep_info, files,
                                        out_dir)
    else:
        error_msg = ("Unknown artifact_type %s. Supported types: 'SFF', "
                     "'FASTQ', 'FASTA', 'FASTA_Sanger', 'per_sample_FASTQ', "
                     "'FASTA_preprocessed', 'Demultiplexed'" % a_type)
        return False, None, error_msg

    status, artifacts, error_msg = reply
    if not status:
        return reply

    # generating html summary
    files = defaultdict(list)
    # artifacts[0].files: there is only one artifact
    for fp, fpt in artifacts[0].files:
        files[fpt].append(fp)
    artifact_information = _generate_html_summary(a_type, files, out_dir)
    summary_fp = f'{out_dir}/index.html'
    with open(summary_fp, 'w') as fp:
        fp.write(artifact_information)

    # inserting the summary into the artifact
    artifacts[0].files.append((summary_fp, 'html_summary'))

    return status, artifacts, error_msg
