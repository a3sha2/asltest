# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
"""
Orchestrating the ASL-preprocessing workflow
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autofunction:: init_asl_preproc_wf
.. autofunction:: init_asl_derivatives_wf

"""

import os
from collections import OrderedDict

import nibabel as nb
from nipype import logging

from nipype.interfaces.fsl import Split as FSLSplit
from nipype.pipeline import engine as pe
from nipype.interfaces import utility as niu

from ...niworkflows.engine.workflows import LiterateWorkflow as Workflow
from ...niworkflows.interfaces.cifti import GenerateCifti
from ...niworkflows.interfaces.utils import DictMerge

from ...config import DEFAULT_MEMORY_MIN_GB
from ...utils.meepi import combine_meepi_source

from ...interfaces import DerivativesDataSink
from ...interfaces.reports import FunctionalSummary

# BOLD workflows
from .confounds import init_asl_confs_wf, init_carpetplot_wf
from .hmc import init_asl_hmc_wf
from .stc import init_asl_stc_wf

from .registration import init_asl_t1_trans_wf, init_asl_reg_wf
from .resampling import (
    init_asl_std_trans_wf,
    init_asl_preproc_trans_wf,
)
from .outputs import init_asl_derivatives_wf
from .util import init_asl_reference_wf


LOGGER = logging.getLogger('nipype.workflow')


def init_asl_preproc_wf(
    asl2t1w_dof,
    asl_file,
    debug,
    fmap_bspline,
    fmap_demean,
    force_syn,
    ignore,
    low_mem,
    omp_nthreads,
    output_dir,
    output_spaces,
    t2s_coreg,
    use_bbr,
    use_syn,
    layout=None,
    num_asl=1,
):
    """
    This workflow controls the functional preprocessing stages of *aslprep*.

    Workflow Graph
        .. workflow::
            :graph2use: orig
            :simple_form: yes

            from aslprep.workflows.asl import init_func_preproc_wf
            from collections import namedtuple, OrderedDict
            BIDSLayout = namedtuple('BIDSLayout', ['root'])
            wf = init_func_preproc_wf(
                aroma_melodic_dim=-200,
                bold2t1w_dof=9,
                bold_file='/completely/made/up/path/sub-01_task-nback_bold.nii.gz',
                cifti_output=False,
                debug=False,
                dummy_scans=None,
                err_on_aroma_warn=False,
                fmap_bspline=True,
                fmap_demean=True,
                force_syn=True,
                freesurfer=True,
                ignore=[],
                low_mem=False,
                medial_surface_nan=False,
                omp_nthreads=1,
                output_dir='.',
                output_spaces=OrderedDict([
                    ('MNI152Lin', {}), ('fsaverage', {'density': '10k'}),
                    ('T1w', {}), ('fsnative', {})]),
                regressors_all_comps=False,
                regressors_dvars_th=1.5,
                regressors_fd_th=0.5,
                reportlets_dir='.',
                t2s_coreg=False,
                use_aroma=False,
                use_bbr=True,
                use_syn=True,
                layout=BIDSLayout('.'),
                num_bold=1,
            )

    Parameters
    ----------
    aroma_melodic_dim : int
        Maximum number of components identified by MELODIC within ICA-AROMA
        (default is -200, ie. no limitation).
    bold2t1w_dof : 6, 9 or 12
        Degrees-of-freedom for BOLD-T1w registration
    bold_file : str
        BOLD series NIfTI file
    cifti_output : bool
        Generate bold CIFTI file in output spaces
    debug : bool
        Enable debugging outputs
    dummy_scans : int or None
        Number of volumes to consider as non steady state
    err_on_aroma_warn : bool
        Do not crash on ICA-AROMA errors
    fmap_bspline : bool
        **Experimental**: Fit B-Spline field using least-squares
    fmap_demean : bool
        Demean voxel-shift map during unwarp
    force_syn : bool
        **Temporary**: Always run SyN-based SDC
    freesurfer : bool
        Enable FreeSurfer functional registration (bbregister) and resampling
        BOLD series to FreeSurfer surface meshes.
    ignore : list
        Preprocessing steps to skip (may include "slicetiming", "fieldmaps")
    low_mem : bool
        Write uncompressed .nii files in some cases to reduce memory usage
    medial_surface_nan : bool
        Replace medial wall values with NaNs on functional GIFTI files
    omp_nthreads : int
        Maximum number of threads an individual process may use
    output_dir : str
        Directory in which to save derivatives
    output_spaces : OrderedDict
        Ordered dictionary where keys are TemplateFlow ID strings (e.g. ``MNI152Lin``,
        ``MNI152NLin6Asym``, ``MNI152NLin2009cAsym``, or ``fsLR``) strings designating
        nonstandard references (e.g. ``T1w`` or ``anat``, ``sbref``, ``run``, etc.),
        or paths pointing to custom templates organized in a TemplateFlow-like structure.
        Values of the dictionary aggregate modifiers (e.g. the value for the key ``MNI152Lin``
        could be ``{'resolution': 2}`` if one wants the resampling to be done on the 2mm
        resolution version of the selected template).
    regressors_all_comps
        Return all CompCor component time series instead of the top fraction
    regressors_dvars_th
        Criterion for flagging DVARS outliers
    regressors_fd_th
        Criterion for flagging framewise displacement outliers
    reportlets_dir : str
        Absolute path of a directory in which reportlets will be temporarily stored
    t2s_coreg : bool
        For multiecho EPI, use the calculated T2*-map for T2*-driven coregistration
    use_aroma : bool
        Perform ICA-AROMA on MNI-resampled functional series
    use_bbr : bool or None
        Enable/disable boundary-based registration refinement.
        If ``None``, test BBR result for distortion before accepting.
        When using ``t2s_coreg``, BBR will be enabled by default unless
        explicitly specified otherwise.
    use_syn : bool
        **Experimental**: Enable ANTs SyN-based susceptibility distortion correction (SDC).
        If fieldmaps are present and enabled, this is not run, by default.
    layout : BIDSLayout
        BIDSLayout structure to enable metadata retrieval
    num_bold : int
        Total number of BOLD files that have been set for preprocessing
        (default is 1)

    Inputs
    ------
    bold_file
        BOLD series NIfTI file
    t1w_preproc
        Bias-corrected structural template image
    t1w_brain
        Skull-stripped ``t1w_preproc``
    t1w_mask
        Mask of the skull-stripped template image
    t1w_dseg
        Segmentation of preprocessed structural image, including
        gray-matter (GM), white-matter (WM) and cerebrospinal fluid (CSF)
    t1w_asec
        Segmentation of structural image, done with FreeSurfer.
    t1w_aparc
        Parcellation of structural image, done with FreeSurfer.
    t1w_tpms
        List of tissue probability maps in T1w space
    anat2std_xfm
        ANTs-compatible affine-and-warp transform file
    std2anat_xfm
        ANTs-compatible affine-and-warp transform file (inverse)
    subjects_dir
        FreeSurfer SUBJECTS_DIR
    subject_id
        FreeSurfer subject ID
    t1w2fsnative_xfm
        LTA-style affine matrix translating from T1w to FreeSurfer-conformed subject space
    fsnative2t1w_xfm
        LTA-style affine matrix translating from FreeSurfer-conformed subject space to T1w

    Outputs
    -------
    bold_t1
        BOLD series, resampled to T1w space
    asl_mask_t1
        BOLD series mask in T1w space
    asl_std
        BOLD series, resampled to template space
    asl_mask_std
        BOLD series mask in template space
    confounds
        TSV of confounds
    surfaces
        BOLD series, resampled to FreeSurfer surfaces
    aroma_noise_ics
        Noise components identified by ICA-AROMA
    melodic_mix
        FSL MELODIC mixing matrix
    bold_cifti
        BOLD CIFTI image
    cifti_variant
        combination of target spaces for `bold_cifti`

    See also
    --------
      * :py:func:`~aslprep.workflows.asl.util.init_asl_reference_wf`
      * :py:func:`~aslprep.workflows.asl.stc.init_asl_stc_wf`
      * :py:func:`~aslprep.workflows.asl.hmc.init_asl_hmc_wf`
      * :py:func:`~aslprep.workflows.asl.t2s.init_asl_t2s_wf`
      * :py:func:`~aslprep.workflows.asl.registration.init_asl_t1_trans_wf`
      * :py:func:`~aslprep.workflows.asl.registration.init_asl_reg_wf`
      * :py:func:`~aslprep.workflows.asl.confounds.init_asl_confounds_wf`
      * :py:func:`~aslprep.workflows.asl.confounds.init_ica_aroma_wf`
      * :py:func:`~aslprep.workflows.asl.resampling.init_asl_std_trans_wf`
      * :py:func:`~aslprep.workflows.asl.resampling.init_asl_preproc_trans_wf`
      * :py:func:`~aslprep.workflows.asl.resampling.init_asl_surf_wf`
      * :py:func:`~aslprep.workflows.fieldmap.pepolar.init_pepolar_unwarp_wf`
      * :py:func:`~aslprep.workflows.fieldmap.init_fmap_estimator_wf`
      * :py:func:`~aslprep.workflows.fieldmap.init_sdc_unwarp_wf`
      * :py:func:`~aslprep.workflows.fieldmap.init_nonlinear_sdc_wf`

    """
    from ...config import NONSTANDARD_REFERENCES
    from ..fieldmap.base import init_sdc_wf  # Avoid circular dependency (#1066)

    # Filter out standard spaces to a separate dict
    std_spaces = OrderedDict([
        (key, modifiers) for key, modifiers in output_spaces.items()
        if key not in NONSTANDARD_REFERENCES])
    volume_std_spaces = OrderedDict([
        (key, modifiers) for key, modifiers in std_spaces.items()
        if not key.startswith('fs')])

    ref_file = asl_file
    mem_gb = {'filesize': 1, 'resampled': 1, 'largemem': 1}
    bold_tlen = 10
    multiecho = isinstance(bold_file, list)

    if multiecho:
        tes = [layout.get_metadata(echo)['EchoTime'] for echo in bold_file]
        ref_file = dict(zip(tes, bold_file))[min(tes)]

    if os.path.isfile(ref_file):
        asl_tlen, mem_gb = _create_mem_gb(ref_file)

    wf_name = _get_wf_name(ref_file)
    LOGGER.log(25, ('Creating asl processing workflow for "%s" (%.2f GB / %d TRs). '
                    'Memory resampled/largemem=%.2f/%.2f GB.'),
               ref_file, mem_gb['filesize'], bold_tlen, mem_gb['resampled'], mem_gb['largemem'])

    # For doc building purposes
    if not hasattr(layout, 'parse_file_entities'):
        LOGGER.log(25, 'No valid layout: building empty workflow.')
        metadata = {
            'RepetitionTime': 2.0,
            'SliceTiming': [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9],
            'PhaseEncodingDirection': 'j',
        }
        fmaps = [{
            'suffix': 'phasediff',
            'phasediff': 'sub-03/ses-2/fmap/sub-03_ses-2_run-1_phasediff.nii.gz',
            'magnitude1': 'sub-03/ses-2/fmap/sub-03_ses-2_run-1_magnitude1.nii.gz',
            'magnitude2': 'sub-03/ses-2/fmap/sub-03_ses-2_run-1_magnitude2.nii.gz',
        }]
        run_stc = True
        multiecho = False
    else:
        # Find associated sbref, if possible
        entities = layout.parse_file_entities(ref_file)
        entities['suffix'] = 'sbref'
        entities['extension'] = ['nii', 'nii.gz']  # Overwrite extensions
        files = layout.get(return_type='file', **entities)
        refbase = os.path.basename(ref_file)
        if 'sbref' in ignore:
            LOGGER.info("Single-band reference files ignored.")
        elif files and multiecho:
            LOGGER.warning("Single-band reference found, but not supported in "
                           "multi-echo workflows at this time. Ignoring.")
        elif files:
            sbref_file = files[0]
            sbbase = os.path.basename(sbref_file)
            if len(files) > 1:
                LOGGER.warning(
                    "Multiple single-band reference files found for {}; using "
                    "{}".format(refbase, sbbase))
            else:
                LOGGER.log(25, "Using single-band reference file {}".format(sbbase))
        else:
            LOGGER.log(25, "No single-band-reference found for {}".format(refbase))

        metadata = layout.get_metadata(ref_file)

        # Find fieldmaps. Options: (phase1|phase2|phasediff|epi|fieldmap|syn)
        fmaps = []
        if 'fieldmaps' not in ignore:
            for fmap in layout.get_fieldmap(ref_file, return_list=True):
                if fmap['suffix'] == 'phase':
                    LOGGER.warning("""\
Found phase1/2 type of fieldmaps, which are not currently supported. \
aslprep will discard them for susceptibility distortion correction. \
Please, follow up on this issue at \
.""")
                else:
                    fmap['metadata'] = layout.get_metadata(fmap[fmap['suffix']])
                    fmaps.append(fmap)

        # Run SyN if forced or in the absence of fieldmap correction
        if force_syn or (use_syn and not fmaps):
            fmaps.append({'suffix': 'syn'})

        # Short circuits: (True and True and (False or 'TooShort')) == 'TooShort'
        run_stc = ("SliceTiming" in metadata and
                   'slicetiming' not in ignore and
                   (_get_series_len(ref_file) > 4 or "TooShort"))

    # Check if MEEPI for T2* coregistration target


    # By default, force-bbr for t2s_coreg unless user specifies otherwise

    use_bbr = True

    # Build workflow
    workflow = Workflow(name=wf_name)
    workflow.__desc__ = """

Functional data preprocessing

: For each of the {num_asl} asl runs found per subject (across all
tasks and sessions), the following preprocessing was performed.
""".format(num_asl=num_asl)

    workflow.__postdesc__ = """\
All resamplings can be performed with *a single interpolation
step* by composing all the pertinent transformations (i.e. head-motion
transform matrices, susceptibility distortion correction when available,
and co-registrations to anatomical and output spaces).
Gridded (volumetric) resamplings were performed using `antsApplyTransforms` (ANTs),
configured with Lanczos interpolation to minimize the smoothing
effects of other kernels [@lanczos].
"""

    inputnode = pe.Node(niu.IdentityInterface(
        fields=['asl_file', 'subjects_dir', 'subject_id',
                't1w_preproc', 't1w_brain', 't1w_mask', 't1w_dseg', 't1w_tpms',
                't1w_aseg', 't1w_aparc',
                'anat2std_xfm', 'std2anat_xfm', 'template',
                'joint_anat2std_xfm', 'joint_std2anat_xfm', 'joint_template']),
        name='inputnode')
    inputnode.inputs.asl_file = asl_file

    outputnode = pe.Node(niu.IdentityInterface(
        fields=['asl_t1', 'asl_t1_ref', 'asl_mask_t1', 'asl_std', 'asl_std_ref', 
                'asl_mask_std','asl_native', 'confounds','confounds_metadata']),
        name='outputnode')

    # BOLD buffer: an identity used as a pointer to either the original BOLD
    # or the STC'ed one for further use.
    boldbuffer = pe.Node(niu.IdentityInterface(fields=['asl_file']), name='aslbuffer')

    summary = pe.Node(
        FunctionalSummary(
            slice_timing=run_stc,
            registration=('FSL'),
            registration_dof=asl2t1w_dof,
            pe_direction=metadata.get("PhaseEncodingDirection"),
            tr=metadata.get("RepetitionTime")),
        name='summary', mem_gb=DEFAULT_MEMORY_MIN_GB, run_without_submitting=True)
    summary.inputs.dummy_scans = dummy_scans

    # CIfTI output: currently, we only support fsaverage{5,6}

    asl_derivatives_wf = init_asl_derivatives_wf(
        bids_root=layout.root,
        metadata=metadata,
        output_dir=output_dir,
        output_spaces=output_spaces,
        standard_spaces=list(std_spaces.keys()),
    )

    workflow.connect([
        (outputnode, func_derivatives_wf, [
            ('asl_t1', 'inputnode.asl_t1'),
            ('asl_t1_ref', 'inputnode.asl_t1_ref'),
            ('asl_mask_t1', 'inputnode.asl_mask_t1'),
            ('asl_native', 'inputnode.asl_native'),
            ('confounds', 'inputnode.confounds'),
            ('confounds_metadata', 'inputnode.confounds_metadata'),
        ]),
    ])

    # Generate a tentative boldref
    bold_reference_wf = init_asl_reference_wf(omp_nthreads=omp_nthreads)
    bold_reference_wf.inputs.inputnode.dummy_scans = dummy_scans
    if sbref_file is not None:
        workflow.connect([
            (val_sbref, bold_reference_wf, [('out_file', 'inputnode.sbref_file')]),
        ])

    # Top-level BOLD splitter
    asl_split = pe.Node(FSLSplit(dimension='t'), name='asl_split',
                         mem_gb=mem_gb['filesize'] * 3)

    # HMC on the asl
    asl_hmc_wf = init_asl_hmc_wf(name='asl_hmc_wf',
                                   mem_gb=mem_gb['filesize'],
                                   omp_nthreads=omp_nthreads)

    # calculate ASL registration to T1w
    asl_reg_wf = init_asl_reg_wf(name='asl_reg_wf',
                                   use_bbr=use_bbr,
                                   asl2t1w_dof=asl2t1w_dof,
                                   mem_gb=mem_gb['resampled'],
                                   omp_nthreads=omp_nthreads,
                                   use_compression=False)

    # apply ASL registration to T1w
    asl_t1_trans_wf = init_asl_t1_trans_wf(name='asl_t1_trans_wf',
                                             use_fieldwarp=(bool(fmaps) or use_syn),
                                             multiecho=multiecho,
                                             mem_gb=mem_gb['resampled'],
                                             omp_nthreads=omp_nthreads,
                                             use_compression=False)

    # get confounds
    asl_confounds_wf = init_asl_confs_wf(
        mem_gb=mem_gb['largemem'],
        metadata=metadata,
        name='asl_confounds_wf')
    asl_confounds_wf.get_node('inputnode').inputs.t1_transform_flags = [False]

    # Apply transforms in 1 shot
    # Only use uncompressed output if AROMA is to be run
    asl_asl_trans_wf = init_asl_preproc_trans_wf(
        mem_gb=mem_gb['resampled'],
        omp_nthreads=omp_nthreads,
        use_compression=not low_mem,
        use_fieldwarp=(bool(fmaps) or use_syn),
        name='asl_asl_trans_wf'
    )
    asl_asl_trans_wf.inputs.inputnode.name_source = ref_file

    # SLICE-TIME CORRECTION (or bypass) #############################################
    if run_stc is True:  # bool('TooShort') == True, so check True explicitly
        asl_stc_wf = init_asl_stc_wf(name='asl_stc_wf', metadata=metadata)
        workflow.connect([
            (asl_reference_wf, asl_stc_wf, [
                ('outputnode.skip_vols', 'inputnode.skip_vols')]),
            (asl_stc_wf, aslbuffer, [('outputnode.stc_file', 'asl_file')]),
        ])
        if not multiecho:
            workflow.connect([
                (bold_reference_wf, bold_stc_wf, [
                    ('outputnode.asl_file', 'inputnode.asl_file')])])
        else:  # for meepi, iterate through stc_wf for all workflows
            meepi_echos = boldbuffer.clone(name='meepi_echos')
            meepi_echos.iterables = ('asl_file', asl_file)
            workflow.connect([
                (meepi_echos, bold_stc_wf, [('asl_file', 'inputnode.asl_file')])])
    elif not multiecho:  # STC is too short or False
        # bypass STC from original BOLD to the splitter through boldbuffer
        workflow.connect([
            (asl_reference_wf, aslbuffer, [('outputnode.asl_file', 'asl_file')])])
    else:
        # for meepi, iterate over all meepi echos to boldbuffer
        aslbuffer.iterables = ('asl_file', asl_file)

    # SDC (SUSCEPTIBILITY DISTORTION CORRECTION) or bypass ##########################
    asl_sdc_wf = init_sdc_wf(
        fmaps, metadata, omp_nthreads=omp_nthreads,
        debug=debug, fmap_demean=fmap_demean, fmap_bspline=fmap_bspline)
    # If no standard space is given, use the default for SyN-SDC
    if not volume_std_spaces or 'MNI152NLin2009cAsym' in volume_std_spaces:
        asl_sdc_wf.inputs.inputnode.template = 'MNI152NLin2009cAsym'
    else:
        asl_sdc_wf.inputs.inputnode.template = next(iter(volume_std_spaces))

    if not fmaps:
        LOGGER.warning('SDC: no fieldmaps found or they were ignored (%s).',
                       ref_file)
    elif fmaps[0]['suffix'] == 'syn':
        LOGGER.warning(
            'SDC: no fieldmaps found or they were ignored. '
            'Using EXPERIMENTAL "fieldmap-less SyN" correction '
            'for dataset %s.', ref_file)
    else:
        LOGGER.log(25, 'SDC: fieldmap estimation of type "%s" intended for %s found.',
                   fmaps[0]['suffix'], ref_file)

    # Overwrite ``out_path_base`` of sdcflows' DataSinks
    for node in asl_sdc_wf.list_node_names():
        if node.split('.')[-1].startswith('ds_'):
            asl_sdc_wf.get_node(node).interface.out_path_base = 'aslprep'

    # MULTI-ECHO EPI DATA #############################################
    #  not multicehco for asl

    # MAIN WORKFLOW STRUCTURE #######################################################
    workflow.connect([
        # Generate early reference
        (inputnode, basl_reference_wf, [('asl_file', 'inputnode.asl_file')]),
        # ASL buffer has slice-time corrected if it was run, original otherwise
        (aslbuffer, asl_split, [('asl_file', 'in_file')]),
        # HMC
        (asl_reference_wf, asl_hmc_wf, [
            ('outputnode.raw_ref_image', 'inputnode.raw_ref_image'),
            ('outputnode.asl_file', 'inputnode.asl_file')]),
        (asl_reference_wf, summary, [
            ('outputnode.algo_dummy_scans', 'algo_dummy_scans')]),
        # EPI-T1 registration workflow
        (inputnode, asl_reg_wf, [
            ('t1w_brain', 'inputnode.t1w_brain'),
            # Undefined if --no-freesurfer, but this is safe
            ('subjects_dir', 'inputnode.subjects_dir'),
            ('subject_id', 'inputnode.subject_id'),
        (inputnode, asl_t1_trans_wf, [
            ('asl_file', 'inputnode.name_source'),
            ('t1w_brain', 'inputnode.t1w_brain'),
            ('t1w_mask', 'inputnode.t1w_mask')]),
        # unused if multiecho, but this is safe
        (asl_hmc_wf, asl_t1_trans_wf, [('outputnode.xforms', 'inputnode.hmc_xforms')]),
        (asl_reg_wf, asl_t1_trans_wf, [
            ('outputnode.itk_asl_to_t1', 'inputnode.itk_asl_to_t1')]),
        (asl_t1_trans_wf, outputnode, [('outputnode.asl_t1', 'asl_t1'),
                                        ('outputnode.asl_t1_ref', 'asl_t1_ref')]),
        (asl_reg_wf, summary, [('outputnode.fallback', 'fallback')]),
        # SDC (or pass-through workflow)
        (inputnode, asl_sdc_wf, [
            ('joint_template', 'inputnode.templates'),
            ('joint_std2anat_xfm', 'inputnode.std2anat_xfm')]),
        (inputnode, asl_sdc_wf, [('t1w_brain', 'inputnode.t1_brain')]),
        (asl_reference_wf, asl_sdc_wf, [
            ('outputnode.ref_image', 'inputnode.asl_ref'),
            ('outputnode.ref_image_brain', 'inputnode.asl_ref_brain'),
            ('outputnode.asl_mask', 'inputnode.asl_mask')]),
        # For t2s_coreg, replace EPI-to-T1w registration inputs
        (asl_sdc_wf if not t2s_coreg else asl_t2s_wf, asl_reg_wf, [
            ('outputnode.asl_ref_brain', 'inputnode.ref_asl_brain')]),
        (als_sdc_wf if not t2s_coreg else asl_t2s_wf, asl_t1_trans_wf, [
            ('outputnode.asl_ref_brain', 'inputnode.ref_asl_brain'),
            ('outputnode.asl_mask', 'inputnode.ref_asl_mask')]),
        (asl_sdc_wf, asl_t1_trans_wf, [
            ('outputnode.out_warp', 'inputnode.fieldwarp')]),
        (asl_sdc_wf, asl_asl_trans_wf, [
            ('outputnode.out_warp', 'inputnode.fieldwarp'),
            ('outputnode.asl_mask', 'inputnode.asl_mask')]),
        (asl_sdc_wf, summary, [('outputnode.method', 'distortion_correction')]),
        # Connect bold_confounds_wf
        (inputnode, asl_confounds_wf, [('t1w_tpms', 'inputnode.t1w_tpms'),
                                        ('t1w_mask', 'inputnode.t1w_mask')]),
        (asl_hmc_wf, asl_confounds_wf, [
            ('outputnode.movpar_file', 'inputnode.movpar_file')]),
        (asl_reg_wf, asl_confounds_wf, [
            ('outputnode.itk_t1_to_asl', 'inputnode.t1_asl_xform')]),
        (asl_reference_wf, asl_confounds_wf, [
            ('outputnode.skip_vols', 'inputnode.skip_vols')]),
        (asl_confounds_wf, outputnode, [
            ('outputnode.confounds_file', 'confounds'),
        ]),
        (asl_confounds_wf, outputnode, [
            ('outputnode.confounds_metadata', 'confounds_metadata'),
        ]),
        # Connect asl_asl_trans_wf
        (asl_split, asl_asl_trans_wf, [
            ('out_files', 'inputnode.asl_file')]),
        (asl_hmc_wf, asl_asl_trans_wf, [
            ('outputnode.xforms', 'inputnode.hmc_xforms')]),
        # Summary
        (outputnode, summary, [('confounds', 'confounds_file')]),
    ])

    # for standard EPI data, pass along correct file
    if not multiecho:
        workflow.connect([
            (inputnode, func_derivatives_wf, [
                ('asl_file', 'inputnode.source_file')]),
            (asl_asl_trans_wf, asl_confounds_wf, [
                ('outputnode.asl', 'inputnode.asl'),
                ('outputnode.asl_mask', 'inputnode.asl_mask')]),
            (asl_split, asl_t1_trans_wf, [
                ('out_files', 'inputnode.asl_split')]),
        ])
    else:  # for meepi, create and use optimal combination
        workflow.connect([
            # update name source for optimal combination
            (inputnode, asl_derivatives_wf, [
                (('asl_file', combine_meepi_source), 'inputnode.source_file')]),
            (asl_asl_trans_wf, skullstrip_asl_wf, [
                ('outputnode.asl', 'inputnode.in_file')]),
            (asl_t2s_wf, asl_confounds_wf, [
                ('outputnode.asl', 'inputnode.asl'),
                ('outputnode.asl_mask', 'inputnode.asl_mask')]),
            (asl_t2s_wf, asl_t1_trans_wf, [
                ('outputnode.asl', 'inputnode.asl_split')]),
        ])

    if fmaps:
        from ..fieldmap.unwarp import init_fmap_unwarp_report_wf
        # Report on ASL correction
        fmap_unwarp_report_wf = init_fmap_unwarp_report_wf()
        workflow.connect([
            (inputnode, fmap_unwarp_report_wf, [
                ('t1w_dseg', 'inputnode.in_seg')]),
            (asl_reference_wf, fmap_unwarp_report_wf, [
                ('outputnode.ref_image', 'inputnode.in_pre')]),
            (asl_reg_wf, fmap_unwarp_report_wf, [
                ('outputnode.itk_t1_to_asl', 'inputnode.in_xfm')]),
            (asl_sdc_wf, fmap_unwarp_report_wf, [
                ('outputnode.asl_ref', 'inputnode.in_post')]),
        ])

        # Overwrite ``out_path_base`` of unwarping DataSinks
        for node in fmap_unwarp_report_wf.list_node_names():
            if node.split('.')[-1].startswith('ds_'):
                fmap_unwarp_report_wf.get_node(node).interface.out_path_base = 'aslprep'

        if force_syn and fmaps[0]['suffix'] != 'syn':
            syn_unwarp_report_wf = init_fmap_unwarp_report_wf(
                name='syn_unwarp_report_wf', forcedsyn=True)
            workflow.connect([
                (inputnode, syn_unwarp_report_wf, [
                    ('t1w_dseg', 'inputnode.in_seg')]),
                (asl_reference_wf, syn_unwarp_report_wf, [
                    ('outputnode.ref_image', 'inputnode.in_pre')]),
                (asl_reg_wf, syn_unwarp_report_wf, [
                    ('outputnode.itk_t1_to_asl', 'inputnode.in_xfm')]),
                (asl_sdc_wf, syn_unwarp_report_wf, [
                    ('outputnode.syn_asl_ref', 'inputnode.in_post')]),
            ])

            # Overwrite ``out_path_base`` of unwarping DataSinks
            for node in syn_unwarp_report_wf.list_node_names():
                if node.split('.')[-1].startswith('ds_'):
                    syn_unwarp_report_wf.get_node(node).interface.out_path_base = 'aslprep'

    # Map final ASL mask into T1w space (if required)
    if 'T1w' in output_spaces or 'anat' in output_spaces:
        from ...niworkflows.interfaces.fixes import (
            FixHeaderApplyTransforms as ApplyTransforms
        )

        aslmask_to_t1w = pe.Node(
            ApplyTransforms(interpolation='MultiLabel', float=True),
            name='aslmask_to_t1w', mem_gb=0.1
        )
        workflow.connect([
            (asl_reg_wf, aslmask_to_t1w, [
                ('outputnode.itk_asl_to_t1', 'transforms')]),
            (asl_t1_trans_wf, aslmask_to_t1w, [
                ('outputnode.asl_mask_t1', 'reference_image')]),
            (asl_asl_trans_wf if not multiecho else bold_t2s_wf, boldmask_to_t1w, [
                ('outputnode.asl_mask', 'input_image')]),
            (boldmask_to_t1w, outputnode, [
                ('output_image', 'asl_mask_t1')]),
        ])

    if set(['func', 'run', 'asl', 'aslref']).intersection(output_spaces):
        workflow.connect([
            (asl_asl_trans_wf, outputnode, [
                ('outputnode.asl', 'asl_native')]),
            (asl_asl_trans_wf, asl_derivatives_wf, [
                ('outputnode.asl_ref', 'inputnode.asl_native_ref'),
                ('outputnode.asl_mask', 'inputnode.asl_mask_native')]),
        ])

    if volume_std_spaces:
        # Apply transforms in 1 shot
        # Only use uncompressed output if AROMA is to be run
        asl_std_trans_wf = init_asl_std_trans_wf(
            mem_gb=mem_gb['resampled'],
            omp_nthreads=omp_nthreads,
            standard_spaces=volume_std_spaces,
            name='asl_std_trans_wf',
            use_compression=not low_mem,
            use_fieldwarp=bool(fmaps),
        )
        workflow.connect([
            (inputnode, asl_std_trans_wf, [
                ('joint_template', 'inputnode.templates'),
                ('joint_anat2std_xfm', 'inputnode.anat2std_xfm'),
                ('asl_file', 'inputnode.name_source')]),
            (asl_hmc_wf, asl_std_trans_wf, [
                ('outputnode.xforms', 'inputnode.hmc_xforms')]),
            (asl_reg_wf, asl_std_trans_wf, [
                ('outputnode.itk_asl_to_t1', 'inputnode.itk_asl_to_t1')]),
            (asl_asl_trans_wf if not multiecho else asl_t2s_wf, asl_std_trans_wf, [
                ('outputnode.asl_mask', 'inputnode.asl_mask')]),
            (asl_sdc_wf, asl_std_trans_wf, [
                ('outputnode.out_warp', 'inputnode.fieldwarp')]),
            (asl_std_trans_wf, outputnode, [('outputnode.asl_std', 'asl_std'),
                                             ('outputnode.asl_std_ref', 'asl_std_ref'),
                                             ('outputnode.asl_mask_std', 'asl_mask_std')]),
        ])


        if 'MNI152NLin2009cAsym' in std_spaces:
            carpetplot_wf = init_carpetplot_wf(
                standard_spaces=std_spaces,
                mem_gb=mem_gb['resampled'],
                metadata=metadata,
                name='carpetplot_wf')
            workflow.connect([
                (inputnode, carpetplot_wf, [
                    ('joint_std2anat_xfm', 'inputnode.std2anat_xfm')]),
                (asl_asl_trans_wf if not multiecho else bold_t2s_wf, carpetplot_wf, [
                    ('outputnode.asl', 'inputnode.asl'),
                    ('outputnode.asl_mask', 'inputnode.asl_mask')]),
                (asl_reg_wf, carpetplot_wf, [
                    ('outputnode.itk_t1_to_asl', 'inputnode.t1_asl_xform')]),
                (asl_confounds_wf, carpetplot_wf, [
                    ('outputnode.confounds_file', 'inputnode.confounds_file')]),
            ])

        if not multiecho:
            workflow.connect([
                (asl_split, asl_std_trans_wf, [
                    ('out_files', 'inputnode.asl_split')])
            ])
        else:
            split_opt_comb = bold_split.clone(name='split_opt_comb')
            workflow.connect([
                (asl_t2s_wf, split_opt_comb, [
                    ('outputnode.asl', 'in_file')]),
                (split_opt_comb, asl_std_trans_wf, [
                    ('out_files', 'inputnode.asl_split')
                ])
            ])

        # Artifacts resampled in MNI space can only be sinked if they
        # were actually generated. See #1348.
        # Uses the parameterized outputnode to generate all outputs
        workflow.connect([
            (asl_std_trans_wf, func_derivatives_wf, [
                ('poutputnode.templates', 'inputnode.template'),
                ('poutputnode.asl_std_ref', 'inputnode.asl_std_ref'),
                ('poutputnode.asl_std', 'inputnode.asl_std'),
                ('poutputnode.asl_mask_std', 'inputnode.asl_mask_std'),
            ]),
        ])

            mrg_conf_metadata = pe.Node(niu.Merge(2), name='merge_confound_metadata',
                                        run_without_submitting=True)
            mrg_conf_metadata2 = pe.Node(DictMerge(), name='merge_confound_metadata2',
                                         run_without_submitting=True)
            workflow.disconnect([
                (asl_confounds_wf, outputnode, [
                    ('outputnode.confounds_file', 'confounds'),
                ]),
                (asl_confounds_wf, outputnode, [
                    ('outputnode.confounds_metadata', 'confounds_metadata'),
                ]),
            ])

    # SURFACES ##################################################################################
    #no surface for now

    # REPORTING ############################################################
    ds_report_summary = pe.Node(
        DerivativesDataSink(desc='summary', keep_dtype=True),
        name='ds_report_summary', run_without_submitting=True,
        mem_gb=DEFAULT_MEMORY_MIN_GB)

    ds_report_validation = pe.Node(
        DerivativesDataSink(base_directory=reportlets_dir,
                            desc='validation', keep_dtype=True),
        name='ds_report_validation', run_without_submitting=True,
        mem_gb=DEFAULT_MEMORY_MIN_GB)

    workflow.connect([
        (summary, ds_report_summary, [('out_report', 'in_file')]),
        (bold_reference_wf, ds_report_validation, [
            ('outputnode.validation_report', 'in_file')]),
    ])

    # Fill-in datasinks of reportlets seen so far
    for node in workflow.list_node_names():
        if node.split('.')[-1].startswith('ds_report'):
            workflow.get_node(node).inputs.base_directory = reportlets_dir
            workflow.get_node(node).inputs.source_file = ref_file

    return workflow


def _get_series_len(asl_fname):
    from ...niworkflows.interfaces.registration import _get_vols_to_discard
    img = nb.load(asl_fname)
    if len(img.shape) < 4:
        return 1

    skip_vols = _get_vols_to_discard(img)

    return img.shape[3] - skip_vols


def _create_mem_gb(asl_fname):
    asl_size_gb = os.path.getsize(asl_fname) / (1024**3)
    asl_tlen = nb.load(asl_fname).shape[-1]
    mem_gb = {
        'filesize': asl_size_gb,
        'resampled': asl_size_gb * 4,
        'largemem': asl_size_gb * (max(asl_tlen / 100, 1.0) + 4),
    }

    return asl_tlen, mem_gb


def _get_wf_name(asl_fname):
    """
    Derive the workflow name for supplied BOLD file.

    >>> _get_wf_name('/completely/made/up/path/sub-01_task-nback_bold.nii.gz')
    'func_preproc_task_nback_wf'
    >>> _get_wf_name('/completely/made/up/path/sub-01_task-nback_run-01_echo-1_bold.nii.gz')
    'func_preproc_task_nback_run_01_echo_1_wf'

    """
    from nipype.utils.filemanip import split_filename
    fname = split_filename(asl_fname)[1]
    fname_nosub = '_'.join(fname.split("_")[1:])
    # if 'echo' in fname_nosub:
    #     fname_nosub = '_'.join(fname_nosub.split("_echo-")[:1]) + "_bold"
    name = "asl_preproc_" + fname_nosub.replace(
        ".", "_").replace(" ", "").replace("-", "_").replace("_asl", "_wf")

    return name


def _to_join(in_file, join_file):
    """Join two tsv files if the join_file is not ``None``."""
    from ...niworkflows.interfaces.utils import JoinTSVColumns
    if join_file is None:
        return in_file

    res = JoinTSVColumns(in_file=in_file, join_file=join_file).run()
    return res.outputs.out_file
