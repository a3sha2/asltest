# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
"""Writing out derivative files."""
from nipype.pipeline import engine as pe
from nipype.interfaces import utility as niu

from ...niworkflows.niworkflows.engine.workflows import LiterateWorkflow as Workflow
from ...niworkflows.niworkflows.interfaces.cifti import CiftiNameSource
from ...niworkflows.niworkflows.interfaces.surf import GiftiNameSource

from ...config import DEFAULT_MEMORY_MIN_GB
from ...interfaces import DerivativesDataSink


def init_asl_derivatives_wf(
    bids_root,
    cifti_output,
    freesurfer,
    metadata,
    output_dir,
    output_spaces,
    standard_spaces,
    fslr_density=None,
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
    fslr_density : str, optional
        Density of fsLR surface (32k or 59k)
    name : str
        This workflow's identifier (default: ``func_derivatives_wf``).

    """
    from smriprep.workflows.outputs import _bids_relative
    workflow = Workflow(name=name)

    inputnode = pe.Node(niu.IdentityInterface(fields=[
        'aroma_noise_ics', 'asl_aparc_std', 'asl_aparc_t1', 'asl_aseg_std',
        'asl_aseg_t1', 'asl_cifti', 'asl_mask_std', 'asl_mask_t1', 'asl_std',
        'asl_std_ref', 'asl_t1', 'asl_t1_ref', 'asl_native', 'asl_native_ref',
        'asl_mask_native', 'cifti_variant', 'cifti_metadata', 'cifti_density',
        'confounds', 'confounds_metadata',
        'source_file', 'surfaces', 'template']),
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

    if set(['func', 'run', 'asl', 'aslref', 'sbref']).intersection(output_spaces):
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
        if freesurfer:
            ds_asl_aseg_t1 = pe.Node(DerivativesDataSink(
                base_directory=output_dir, space='T1w', desc='aseg', suffix='dseg'),
                name='ds_asl_aseg_t1', run_without_submitting=True,
                mem_gb=DEFAULT_MEMORY_MIN_GB)
            ds_asl_aparc_t1 = pe.Node(DerivativesDataSink(
                base_directory=output_dir, space='T1w', desc='aparcaseg', suffix='dseg'),
                name='ds_asl_aparc_t1', run_without_submitting=True,
                mem_gb=DEFAULT_MEMORY_MIN_GB)
            workflow.connect([
                (inputnode, ds_asl_aseg_t1, [('source_file', 'source_file'),
                                              ('asl_aseg_t1', 'in_file')]),
                (inputnode, ds_asl_aparc_t1, [('source_file', 'source_file'),
                                               ('asl_aparc_t1', 'in_file')]),
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

        if freesurfer:
            ds_asl_aseg_std = pe.Node(DerivativesDataSink(
                base_directory=output_dir, desc='aseg', suffix='dseg'),
                name='ds_asl_aseg_std', run_without_submitting=True,
                mem_gb=DEFAULT_MEMORY_MIN_GB)
            ds_asl_aparc_std = pe.Node(DerivativesDataSink(
                base_directory=output_dir, desc='aparcaseg', suffix='dseg'),
                name='ds_asl_aparc_std', run_without_submitting=True,
                mem_gb=DEFAULT_MEMORY_MIN_GB)
            workflow.connect([
                (inputnode, ds_asl_aseg_std, [('source_file', 'source_file'),
                                               ('asl_aseg_std', 'in_file'),
                                               ('template', 'space')]),
                (inputnode, ds_asl_aparc_std, [('source_file', 'source_file'),
                                                ('asl_aparc_std', 'in_file'),
                                                ('template', 'space')]),
            ])

    # fsaverage space
    if freesurfer and any(space.startswith('fs') for space in output_spaces.keys()):

        extract_surf_info = pe.MapNode(niu.Function(function=_extract_surf_info),
                                       iterfield=['in_file'], name='extract_surf_info',
                                       mem_gb=DEFAULT_MEMORY_MIN_GB, run_without_submitting=True)
        extract_surf_info.inputs.density = fslr_density

        name_surfs = pe.MapNode(GiftiNameSource(
            pattern=r'(?P<LR>[lr])h.(?P<space>\w+).gii',
            template='space-{space}{den}_hemi-{LR}.func'),
            iterfield=['in_file', 'template_kwargs'], name='name_surfs',
            mem_gb=DEFAULT_MEMORY_MIN_GB, run_without_submitting=True)

        ds_asl_surfs = pe.MapNode(DerivativesDataSink(base_directory=output_dir),
                                   iterfield=['in_file', 'suffix'], name='ds_asl_surfs',
                                   run_without_submitting=True,
                                   mem_gb=DEFAULT_MEMORY_MIN_GB)

        workflow.connect([
            (inputnode, extract_surf_info, [('surfaces', 'in_file')]),
            (inputnode, name_surfs, [('surfaces', 'in_file')]),
            (extract_surf_info, name_surfs, [('out', 'template_kwargs')]),
            (inputnode, ds_asl_surfs, [('source_file', 'source_file'),
                                        ('surfaces', 'in_file')]),
            (name_surfs, ds_asl_surfs, [('out_name', 'suffix')]),
        ])

        # CIFTI output
        if cifti_output:
            name_cifti = pe.MapNode(
                CiftiNameSource(), iterfield=['variant', 'density'], name='name_cifti',
                mem_gb=DEFAULT_MEMORY_MIN_GB, run_without_submitting=True)
            cifti_asls = pe.MapNode(
                DerivativesDataSink(base_directory=output_dir, compress=False),
                iterfield=['in_file', 'suffix'], name='cifti_asls',
                run_without_submitting=True, mem_gb=DEFAULT_MEMORY_MIN_GB)
            cifti_key = pe.MapNode(DerivativesDataSink(
                base_directory=output_dir), iterfield=['in_file', 'suffix'],
                name='cifti_key', run_without_submitting=True,
                mem_gb=DEFAULT_MEMORY_MIN_GB)
            workflow.connect([
                (inputnode, name_cifti, [('cifti_variant', 'variant'),
                                         ('cifti_density', 'density')]),
                (inputnode, cifti_asls, [('asl_cifti', 'in_file'),
                                          ('source_file', 'source_file')]),
                (name_cifti, cifti_asls, [('out_name', 'suffix')]),
                (name_cifti, cifti_key, [('out_name', 'suffix')]),
                (inputnode, cifti_key, [('source_file', 'source_file'),
                                        ('cifti_metadata', 'in_file')]),
            ])

    return workflow


def _extract_surf_info(in_file, density):
    import os
    info = {'den': ''}
    if 'fsLR' in os.path.basename(in_file):
        info['den'] = '_den-{}'.format(density)
    return info
