# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
"""Writing out derivative files."""
from nipype.pipeline import engine as pe
from nipype.interfaces import utility as niu

from ...niworkflows.engine.workflows import LiterateWorkflow as Workflow
from ...niworkflows.interfaces.cifti import CiftiNameSource
from ...niworkflows.interfaces.surf import GiftiNameSource

from ...config import DEFAULT_MEMORY_MIN_GB
from ...interfaces import DerivativesDataSink


def init_asl_derivatives_wf(
    bids_root,
    metadata,
    output_dir,
    output_spaces,
    standard_spaces,
    name='asl_derivatives_wf',
):
    """
    Set up a battery of datasinks to store derivatives in the right location.

    Parameters
    ----------
    bids_root : str
        Original BIDS dataset path.
    cifti_output : bool
        Whether the ``--cifti-output`` flag was set.
    freesurfer : bool
        Whether FreeSurfer anatomical processing was run.
    metadata : dict
        Metadata dictionary associated to the asl run.
    output_dir : str
        Where derivatives should be written out to.
    output_spaces : OrderedDict
        List of selected ``--output-spaces``.
    use_aroma : bool
        Whether ``--use-aroma`` flag was set.
    name : str
        This workflow's identifier (default: ``asl_derivatives_wf``).

    """
    from smriprep.workflows.outputs import _bids_relative
    workflow = Workflow(name=name)

    inputnode = pe.Node(niu.IdentityInterface(fields=[
        'asl_mask_std', 'asl_mask_t1', 'asl_std',
        'asl_std_ref', 'asl_t1', 'asl_t1_ref', 'asl_native', 'asl_native_ref',
        'asl_mask_native', 'confounds', 'confounds_metadata', 
        'source_file', 'template']),
        name='inputnode')

    raw_sources = pe.Node(niu.Function(function=_bids_relative), name='raw_sources')
    raw_sources.inputs.bids_root = bids_root

    ds_confounds = pe.Node(DerivativesDataSink(
        base_directory=output_dir, desc='confounds', suffix='regressors'),
        name="ds_confounds", run_without_submitting=True,
        mem_gb=DEFAULT_MEMORY_MIN_GB)
    workflow.connect([
        (inputnode, raw_sources, [('source_file', 'in_files')]),
        (inputnode, ds_confounds, [('source_file', 'source_file'),
                                   ('confounds', 'in_file'),
                                   ('confounds_metadata', 'meta_dict')]),
    ])

    if set(['asl', 'run', 'aslref']).intersection(output_spaces):
        ds_asl_native = pe.Node(
            DerivativesDataSink(base_directory=output_dir, desc='preproc',
                                keep_dtype=True, compress=True, SkullStripped=False,
                                RepetitionTime=metadata.get('RepetitionTime'),
                                TaskName=metadata.get('TaskName')),
            name='ds_asl_native', run_without_submitting=True,
            mem_gb=DEFAULT_MEMORY_MIN_GB)
        ds_asl_native_ref = pe.Node(
            DerivativesDataSink(base_directory=output_dir, suffix='aslref', compress=True),
            name='ds_asl_native_ref', run_without_submitting=True,
            mem_gb=DEFAULT_MEMORY_MIN_GB)
        ds_asl_mask_native = pe.Node(
            DerivativesDataSink(base_directory=output_dir, desc='brain',
                                suffix='mask', compress=True),
            name='ds_asl_mask_native', run_without_submitting=True,
            mem_gb=DEFAULT_MEMORY_MIN_GB)

        workflow.connect([
            (inputnode, ds_asl_native, [('source_file', 'source_file'),
                                         ('asl_native', 'in_file')]),
            (inputnode, ds_asl_native_ref, [('source_file', 'source_file'),
                                             ('asl_native_ref', 'in_file')]),
            (inputnode, ds_asl_mask_native, [('source_file', 'source_file'),
                                              ('asl_mask_native', 'in_file')]),
            (raw_sources, ds_asl_mask_native, [('out', 'RawSources')]),
        ])

    # Resample to T1w space
    if 'T1w' in output_spaces or 'anat' in output_spaces:
        ds_asl_t1 = pe.Node(
            DerivativesDataSink(base_directory=output_dir, space='T1w', desc='preproc',
                                keep_dtype=True, compress=True, SkullStripped=False,
                                RepetitionTime=metadata.get('RepetitionTime'),
                                TaskName=metadata.get('TaskName')),
            name='ds_asl_t1', run_without_submitting=True,
            mem_gb=DEFAULT_MEMORY_MIN_GB)
        ds_asl_t1_ref = pe.Node(
            DerivativesDataSink(base_directory=output_dir, space='T1w',
                                suffix='aslref', compress=True),
            name='ds_asl_t1_ref', run_without_submitting=True,
            mem_gb=DEFAULT_MEMORY_MIN_GB)

        ds_asl_mask_t1 = pe.Node(
            DerivativesDataSink(base_directory=output_dir, space='T1w', desc='brain',
                                suffix='mask', compress=True),
            name='ds_asl_mask_t1', run_without_submitting=True,
            mem_gb=DEFAULT_MEMORY_MIN_GB)
        workflow.connect([
            (inputnode, ds_asl_t1, [('source_file', 'source_file'),
                                     ('asl_t1', 'in_file')]),
            (inputnode, ds_asl_t1_ref, [('source_file', 'source_file'),
                                         ('asl_t1_ref', 'in_file')]),
            (inputnode, ds_asl_mask_t1, [('source_file', 'source_file'),
                                          ('asl_mask_t1', 'in_file')]),
            (raw_sources, ds_asl_mask_t1, [('out', 'RawSources')]),
        ])

    # Resample to template (default: MNI)
    if standard_spaces:
        ds_asl_std = pe.Node(
            DerivativesDataSink(base_directory=output_dir, desc='preproc',
                                keep_dtype=True, compress=True, SkullStripped=False,
                                RepetitionTime=metadata.get('RepetitionTime'),
                                TaskName=metadata.get('TaskName')),
            name='ds_asl_std', run_without_submitting=True,
            mem_gb=DEFAULT_MEMORY_MIN_GB)
        ds_asl_std_ref = pe.Node(
            DerivativesDataSink(base_directory=output_dir, suffix='aslref'),
            name='ds_asl_std_ref', run_without_submitting=True,
            mem_gb=DEFAULT_MEMORY_MIN_GB)

        ds_asl_mask_std = pe.Node(
            DerivativesDataSink(base_directory=output_dir, desc='brain',
                                suffix='mask'),
            name='ds_asl_mask_std', run_without_submitting=True,
            mem_gb=DEFAULT_MEMORY_MIN_GB)
        workflow.connect([
            (inputnode, ds_asl_std, [('source_file', 'source_file'),
                                      ('asl_std', 'in_file'),
                                      ('template', 'space')]),
            (inputnode, ds_asl_std_ref, [('source_file', 'source_file'),
                                          ('asl_std_ref', 'in_file'),
                                          ('template', 'space')]),
            (inputnode, ds_asl_mask_std, [('source_file', 'source_file'),
                                           ('asl_mask_std', 'in_file'),
                                           ('template', 'space')]),
            (raw_sources, ds_asl_mask_std, [('out', 'RawSources')]),
        ])


    return workflow
