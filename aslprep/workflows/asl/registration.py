# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
"""
Registration workflows
++++++++++++++++++++++

.. autofunction:: init_asl_reg_wf
.. autofunction:: init_asl_t1_trans_wf
.. autofunction:: init_bbreg_wf
.. autofunction:: init_fsl_bbr_wf

"""

import os
import os.path as op

import pkg_resources as pkgr

from nipype.pipeline import engine as pe
from nipype import logging
from nipype.interfaces import utility as niu, fsl, c3
from niworkflows.engine.workflows import LiterateWorkflow as Workflow
# See https://github.com/poldracklab/aslprep/issues/768
from ...niworkflows.interfaces.freesurfer import (
    PatchedConcatenateLTA as ConcatenateLTA,
    PatchedBBRegisterRPT as BBRegisterRPT,
    PatchedMRICoregRPT as MRICoregRPT,
    PatchedLTAConvert as LTAConvert)
from ...niworkflows.interfaces.fixes import FixHeaderApplyTransforms as ApplyTransforms
from ...niworkflows.interfaces.images import extract_wm
from ...niworkflows.interfaces.itk import MultiApplyTransforms
from ...niworkflows.interfaces.registration import FLIRTRPT
from ...niworkflows.interfaces.utils import GenerateSamplingReference
from ...niworkflows.interfaces.nilearn import Merge

from ...config import DEFAULT_MEMORY_MIN_GB
from ...interfaces import DerivativesDataSink


LOGGER = logging.getLogger('nipype.workflow')


def init_asl_reg_wf(use_bbr, asl2t1w_dof, mem_gb, omp_nthreads,
                     use_compression=True, write_report=True, name='asl_reg_wf'):
    """
    Build a workflow to run same-subject, asl-to-T1w image-registration.

    Calculates the registration between a reference asl image and T1w-space
    using a boundary-based registration (BBR) cost function.
    If FreeSurfer-based preprocessing is enabled, the ``bbregister`` utility
    is used to align the asl images to the reconstructed subject, and the
    resulting transform is adjusted to target the T1 space.
    If FreeSurfer-based preprocessing is disabled, FSL FLIRT is used with the
    BBR cost function to directly target the T1 space.

    Workflow Graph
        .. workflow::
            :graph2use: orig
            :simple_form: yes

            from aslprep.workflows.asl.registration import init_asl_reg_wf
            wf = init_asl_reg_wf(freesurfer=True,
                                  mem_gb=3,
                                  omp_nthreads=1,
                                  use_bbr=True,
                                  asl2t1w_dof=9)

    Parameters
    ----------
    freesurfer : bool
        Enable FreeSurfer functional registration (bbregister)
    use_bbr : bool or None
        Enable/disable boundary-based registration refinement.
        If ``None``, test BBR result for distortion before accepting.
    asl2t1w_dof : 6, 9 or 12
        Degrees-of-freedom for asl-T1w registration
    mem_gb : float
        Size of asl file in GB
    omp_nthreads : int
        Maximum number of threads an individual process may use
    name : str
        Name of workflow (default: ``asl_reg_wf``)
    use_compression : bool
        Save registered asl series as ``.nii.gz``
    use_fieldwarp : bool
        Include SDC warp in single-shot transform from asl to T1
    write_report : bool
        Whether a reportlet should be stored

    Inputs
    ------
    ref_asl_brain
        Reference image to which asl series is aligned
        If ``fieldwarp == True``, ``ref_asl_brain`` should be unwarped
    t1w_brain
        Skull-stripped ``t1w_preproc``
    t1w_dseg
        Segmentation of preprocessed structural image, including
        gray-matter (GM), white-matter (WM) and cerebrospinal fluid (CSF)
    subjects_dir
        FreeSurfer SUBJECTS_DIR
    subject_id
        FreeSurfer subject ID
    fsnative2t1w_xfm
        LTA-style affine matrix translating from FreeSurfer-conformed subject space to T1w

    Outputs
    -------
    itk_asl_to_t1
        Affine transform from ``ref_asl_brain`` to T1 space (ITK format)
    itk_t1_to_asl
        Affine transform from T1 space to asl space (ITK format)
    fallback
        Boolean indicating whether BBR was rejected (mri_coreg registration returned)

    See also
    --------
      * :py:func:`~aslprep.workflows.asl.registration.init_bbreg_wf`
      * :py:func:`~aslprep.workflows.asl.registration.init_fsl_bbr_wf`

    """
    workflow = Workflow(name=name)
    inputnode = pe.Node(
        niu.IdentityInterface(
            fields=['ref_asl_brain', 't1w_brain', 't1w_dseg',
                    'subjects_dir', 'subject_id']),
        name='inputnode'
    )

    outputnode = pe.Node(
        niu.IdentityInterface(fields=[
            'itk_asl_to_t1', 'itk_t1_to_asl', 'fallback']),
        name='outputnode'
    )
    bbr_wf = init_fsl_bbr_wf(use_bbr=use_bbr, asl2t1w_dof=asl2t1w_dof)

    workflow.connect([
        (inputnode, bbr_wf, [
            ('ref_asl_brain', 'inputnode.in_file'),
            ('subjects_dir', 'inputnode.subjects_dir'),
            ('subject_id', 'inputnode.subject_id'),
            ('t1w_dseg', 'inputnode.t1w_dseg'),
            ('t1w_brain', 'inputnode.t1w_brain')]),
        (bbr_wf, outputnode, [('outputnode.itk_asl_to_t1', 'itk_asl_to_t1'),
                              ('outputnode.itk_t1_to_asl', 'itk_t1_to_asl'),
                              ('outputnode.fallback', 'fallback')]),
    ])

    if write_report:
        ds_report_reg = pe.Node(
            DerivativesDataSink(keep_dtype=True),
            name='ds_report_reg', run_without_submitting=True,
            mem_gb=DEFAULT_MEMORY_MIN_GB)

        def _asl_reg_suffix(fallback):
            if fallback:
                return 'coreg' 
            

        workflow.connect([
            (bbr_wf, ds_report_reg, [
                ('outputnode.out_report', 'in_file'),
                (('outputnode.fallback', _asl_reg_suffix), 'desc')]),
        ])

    return workflow


def init_asl_t1_trans_wf(mem_gb, omp_nthreads, multiecho=False, use_fieldwarp=False,
                          use_compression=True, name='asl_t1_trans_wf'):
    """
    Co-register the reference asl image to T1w-space.

    The workflow uses :abbr:`BBR (boundary-based registration)`.

    Workflow Graph
        .. workflow::
            :graph2use: orig
            :simple_form: yes

            from aslprep.workflows.asl.registration import init_asl_t1_trans_wf
            wf = init_asl_t1_trans_wf(freesurfer=True,
                                       mem_gb=3,
                                       omp_nthreads=1)

    Parameters
    ----------
    freesurfer : bool
        Enable FreeSurfer functional registration (bbregister)
    use_fieldwarp : bool
        Include SDC warp in single-shot transform from asl to T1
    multiecho : bool
        If multiecho data was supplied, HMC already performed
    mem_gb : float
        Size of asl file in GB
    omp_nthreads : int
        Maximum number of threads an individual process may use
    use_compression : bool
        Save registered asl series as ``.nii.gz``
    name : str
        Name of workflow (default: ``asl_reg_wf``)

    Inputs
    ------
    name_source
        asl series NIfTI file
        Used to recover original information lost during processing
    ref_asl_brain
        Reference image to which asl series is aligned
        If ``fieldwarp == True``, ``ref_asl_brain`` should be unwarped
    ref_asl_mask
        Skull-stripping mask of reference image
    t1w_brain
        Skull-stripped bias-corrected structural template image
    t1w_mask
        Mask of the skull-stripped template image
    t1w_aseg
        FreeSurfer's ``aseg.mgz`` atlas projected into the T1w reference
        (only if ``recon-all`` was run).
    t1w_aparc
        FreeSurfer's ``aparc+aseg.mgz`` atlas projected into the T1w reference
        (only if ``recon-all`` was run).
    asl_split
        Individual 3D asl volumes, not motion corrected
    hmc_xforms
        List of affine transforms aligning each volume to ``ref_image`` in ITK format
    itk_asl_to_t1
        Affine transform from ``ref_asl_brain`` to T1 space (ITK format)
    fieldwarp
        a :abbr:`DFM (displacements field map)` in ITK format

    Outputs
    -------
    asl_t1
        Motion-corrected asl series in T1 space
    asl_t1_ref
        Reference, contrast-enhanced summary of the motion-corrected asl series in T1w space
    asl_mask_t1
        asl mask in T1 space
    asl_aseg_t1
        FreeSurfer's ``aseg.mgz`` atlas, in T1w-space at the asl resolution
        (only if ``recon-all`` was run).
    asl_aparc_t1
        FreeSurfer's ``aparc+aseg.mgz`` atlas, in T1w-space at the asl resolution
        (only if ``recon-all`` was run).

    See also
    --------
      * :py:func:`~aslprep.workflows.asl.registration.init_bbreg_wf`
      * :py:func:`~aslprep.workflows.asl.registration.init_fsl_bbr_wf`

    """
    from .util import init_asl_reference_wf
    workflow = Workflow(name=name)
    inputnode = pe.Node(
        niu.IdentityInterface(
            fields=['name_source', 'ref_asl_brain', 'ref_asl_mask',
                    't1w_brain', 't1w_mask', 't1w_aseg', 't1w_aparc',
                    'asl_split', 'fieldwarp', 'hmc_xforms',
                    'itk_asl_to_t1']),
        name='inputnode'
    )

    outputnode = pe.Node(
        niu.IdentityInterface(fields=[
            'asl_t1', 'asl_t1_ref', 'asl_mask_t1']),
        name='outputnode'
    )

    gen_ref = pe.Node(GenerateSamplingReference(), name='gen_ref',
                      mem_gb=0.3)  # 256x256x256 * 64 / 8 ~ 150MB

    mask_t1w_tfm = pe.Node(
        ApplyTransforms(interpolation='MultiLabel', float=True),
        name='mask_t1w_tfm', mem_gb=0.1
    )

    workflow.connect([
        (inputnode, gen_ref, [('ref_asl_brain', 'moving_image'),
                              ('t1w_brain', 'fixed_image'),
                              ('t1w_mask', 'fov_mask')]),
        (inputnode, mask_t1w_tfm, [('ref_asl_mask', 'input_image')]),
        (gen_ref, mask_t1w_tfm, [('out_file', 'reference_image')]),
        (inputnode, mask_t1w_tfm, [('itk_asl_to_t1', 'transforms')]),
        (mask_t1w_tfm, outputnode, [('output_image', 'asl_mask_t1')]),
    ])


    asl_to_t1w_transform = pe.Node(
        MultiApplyTransforms(interpolation="LanczosWindowedSinc", float=True, copy_dtype=True),
        name='asl_to_t1w_transform', mem_gb=mem_gb * 3 * omp_nthreads, n_procs=omp_nthreads)

    # merge 3D volumes into 4D timeseries
    merge = pe.Node(Merge(compress=use_compression), name='merge', mem_gb=mem_gb)

    # Generate a reference on the target T1w space
    gen_final_ref = init_asl_reference_wf(omp_nthreads, pre_mask=True)

    if not multiecho:
        # Merge transforms placing the head motion correction last
        nforms = 2 + int(use_fieldwarp)
        merge_xforms = pe.Node(niu.Merge(nforms), name='merge_xforms',
                               run_without_submitting=True, mem_gb=DEFAULT_MEMORY_MIN_GB)
        if use_fieldwarp:
            workflow.connect([
                (inputnode, merge_xforms, [('fieldwarp', 'in2')])
            ])

        workflow.connect([
            # merge transforms
            (inputnode, merge_xforms, [
                ('hmc_xforms', 'in%d' % nforms),
                ('itk_asl_to_t1', 'in1')]),
            (merge_xforms, asl_to_t1w_transform, [('out', 'transforms')]),
            (inputnode, asl_to_t1w_transform, [('asl_split', 'input_image')]),
        ])

    else:
        from nipype.interfaces.fsl import Split as FSLSplit
        asl_split = pe.Node(FSLSplit(dimension='t'), name='asl_split',
                             mem_gb=DEFAULT_MEMORY_MIN_GB)

        workflow.connect([
            (inputnode, asl_split, [('asl_split', 'in_file')]),
            (asl_split, asl_to_t1w_transform, [('out_files', 'input_image')]),
            (inputnode, asl_to_t1w_transform, [('itk_asl_to_t1', 'transforms')]),
        ])

    workflow.connect([
        (inputnode, merge, [('name_source', 'header_source')]),
        (gen_ref, asl_to_t1w_transform, [('out_file', 'reference_image')]),
        (asl_to_t1w_transform, merge, [('out_files', 'in_files')]),
        (merge, gen_final_ref, [('out_file', 'inputnode.asl_file')]),
        (mask_t1w_tfm, gen_final_ref, [('output_image', 'inputnode.asl_mask')]),
        (merge, outputnode, [('out_file', 'asl_t1')]),
        (gen_final_ref, outputnode, [('outputnode.ref_image', 'asl_t1_ref')]),
    ])

    return workflow



def init_fsl_bbr_wf(use_bbr, asl2t1w_dof, name='fsl_bbr_wf'):
    """
    Build a workflow to run FSL's ``flirt``.

    This workflow uses FSL FLIRT to register a asl image to a T1-weighted
    structural image, using a boundary-based registration (BBR) cost function.
    It is a counterpart to :py:func:`~aslprep.workflows.asl.registration.init_bbreg_wf`,
    which performs the same task using FreeSurfer's ``bbregister``.

    The ``use_bbr`` option permits a high degree of control over registration.
    If ``False``, standard, rigid coregistration will be performed by FLIRT.
    If ``True``, FLIRT-BBR will be seeded with the initial transform found by
    the rigid coregistration.
    If ``None``, after FLIRT-BBR is run, the resulting affine transform
    will be compared to the initial transform found by FLIRT.
    Excessive deviation will result in rejecting the BBR refinement and
    accepting the original, affine registration.

    Workflow Graph
        .. workflow ::
            :graph2use: orig
            :simple_form: yes

            from aslprep.workflows.asl.registration import init_fsl_bbr_wf
            wf = init_fsl_bbr_wf(use_bbr=True, asl2t1w_dof=9)


    Parameters
    ----------
    use_bbr : bool or None
        Enable/disable boundary-based registration refinement.
        If ``None``, test BBR result for distortion before accepting.
    asl2t1w_dof : 6, 9 or 12
        Degrees-of-freedom for asl-T1w registration
    name : str, optional
        Workflow name (default: fsl_bbr_wf)

    Inputs
    ------
    in_file
        Reference asl image to be registered
    t1w_brain
        Skull-stripped T1-weighted structural image
    t1w_dseg
        FAST segmentation of ``t1w_brain``
    fsnative2t1w_xfm
        Unused (see :py:func:`~aslprep.workflows.asl.registration.init_bbreg_wf`)
    subjects_dir
        Unused (see :py:func:`~aslprep.workflows.asl.registration.init_bbreg_wf`)
    subject_id
        Unused (see :py:func:`~aslprep.workflows.asl.registration.init_bbreg_wf`)

    Outputs
    -------
    itk_asl_to_t1
        Affine transform from ``ref_asl_brain`` to T1w space (ITK format)
    itk_t1_to_asl
        Affine transform from T1 space to asl space (ITK format)
    out_report
        Reportlet for assessing registration quality
    fallback
        Boolean indicating whether BBR was rejected (rigid FLIRT registration returned)

    """
    workflow = Workflow(name=name)
    workflow.__desc__ = """\
The asl reference was then co-registered to the T1w reference using
`flirt` [FSL {fsl_ver}, @flirt] with the boundary-based registration [@bbr]
cost-function.
Co-registration was configured with nine degrees of freedom to account
for distortions remaining in the asl reference.
""".format(fsl_ver=FLIRTRPT().version or '<ver>')

    inputnode = pe.Node(
        niu.IdentityInterface([
            'in_file',
            'subjects_dir', 'subject_id',  # BBRegister
            't1w_dseg', 't1w_brain']),  # FLIRT BBR
        name='inputnode')
    outputnode = pe.Node(
        niu.IdentityInterface(['itk_asl_to_t1', 'itk_t1_to_asl', 'out_report', 'fallback']),
        name='outputnode')

    wm_mask = pe.Node(niu.Function(function=extract_wm), name='wm_mask')
    flt_bbr_init = pe.Node(FLIRTRPT(dof=6, generate_report=not use_bbr,
                                    uses_qform=True), name='flt_bbr_init')

    invt_bbr = pe.Node(fsl.ConvertXFM(invert_xfm=True), name='invt_bbr',
                       mem_gb=DEFAULT_MEMORY_MIN_GB)

    #  asl to T1 transform matrix is from fsl, using c3 tools to convert to
    #  something ANTs will like.
    fsl2itk_fwd = pe.Node(c3.C3dAffineTool(fsl2ras=True, itk_transform=True),
                          name='fsl2itk_fwd', mem_gb=DEFAULT_MEMORY_MIN_GB)
    fsl2itk_inv = pe.Node(c3.C3dAffineTool(fsl2ras=True, itk_transform=True),
                          name='fsl2itk_inv', mem_gb=DEFAULT_MEMORY_MIN_GB)

    workflow.connect([
        (inputnode, flt_bbr_init, [('in_file', 'in_file'),
                                   ('t1w_brain', 'reference')]),
        (inputnode, fsl2itk_fwd, [('t1w_brain', 'reference_file'),
                                  ('in_file', 'source_file')]),
        (inputnode, fsl2itk_inv, [('in_file', 'reference_file'),
                                  ('t1w_brain', 'source_file')]),
        (invt_bbr, fsl2itk_inv, [('out_file', 'transform_file')]),
        (fsl2itk_fwd, outputnode, [('itk_transform', 'itk_asl_to_t1')]),
        (fsl2itk_inv, outputnode, [('itk_transform', 'itk_t1_to_asl')]),
    ])

    # Short-circuit workflow building, use rigid registration
    if use_bbr is False:
        workflow.connect([
            (flt_bbr_init, invt_bbr, [('out_matrix_file', 'in_file')]),
            (flt_bbr_init, fsl2itk_fwd, [('out_matrix_file', 'transform_file')]),
            (flt_bbr_init, outputnode, [('out_report', 'out_report')]),
        ])
        outputnode.inputs.fallback = True

        return workflow

    flt_bbr = pe.Node(
        FLIRTRPT(cost_func='bbr', dof=asl2t1w_dof, generate_report=True),
        name='flt_bbr')

    FSLDIR = os.getenv('FSLDIR')
    if FSLDIR:
        flt_bbr.inputs.schedule = op.join(FSLDIR, 'etc/flirtsch/bbr.sch')
    else:
        # Should mostly be hit while building docs
        LOGGER.warning("FSLDIR unset - using packaged BBR schedule")
        flt_bbr.inputs.schedule = pkgr.resource_filename('aslprep', 'data/flirtsch/bbr.sch')

    workflow.connect([
        (inputnode, wm_mask, [('t1w_dseg', 'in_seg')]),
        (inputnode, flt_bbr, [('in_file', 'in_file'),
                              ('t1w_brain', 'reference')]),
        (flt_bbr_init, flt_bbr, [('out_matrix_file', 'in_matrix_file')]),
        (wm_mask, flt_bbr, [('out', 'wm_seg')]),
    ])

    # Short-circuit workflow building, use boundary-based registration
    if use_bbr is True:
        workflow.connect([
            (flt_bbr, invt_bbr, [('out_matrix_file', 'in_file')]),
            (flt_bbr, fsl2itk_fwd, [('out_matrix_file', 'transform_file')]),
            (flt_bbr, outputnode, [('out_report', 'out_report')]),
        ])
        outputnode.inputs.fallback = False

        return workflow

    transforms = pe.Node(niu.Merge(2), run_without_submitting=True, name='transforms')
    reports = pe.Node(niu.Merge(2), run_without_submitting=True, name='reports')

    compare_transforms = pe.Node(niu.Function(function=compare_xforms), name='compare_transforms')

    select_transform = pe.Node(niu.Select(), run_without_submitting=True, name='select_transform')
    select_report = pe.Node(niu.Select(), run_without_submitting=True, name='select_report')

    fsl_to_lta = pe.MapNode(LTAConvert(out_lta=True), iterfield=['in_fsl'],
                            name='fsl_to_lta')

    workflow.connect([
        (flt_bbr, transforms, [('out_matrix_file', 'in1')]),
        (flt_bbr_init, transforms, [('out_matrix_file', 'in2')]),
        # Convert FSL transforms to LTA (RAS2RAS) transforms and compare
        (inputnode, fsl_to_lta, [('in_file', 'source_file'),
                                 ('t1w_brain', 'target_file')]),
        (transforms, fsl_to_lta, [('out', 'in_fsl')]),
        (fsl_to_lta, compare_transforms, [('out_lta', 'lta_list')]),
        (compare_transforms, outputnode, [('out', 'fallback')]),
        # Select output transform
        (transforms, select_transform, [('out', 'inlist')]),
        (compare_transforms, select_transform, [('out', 'index')]),
        (select_transform, invt_bbr, [('out', 'in_file')]),
        (select_transform, fsl2itk_fwd, [('out', 'transform_file')]),
        (flt_bbr, reports, [('out_report', 'in1')]),
        (flt_bbr_init, reports, [('out_report', 'in2')]),
        (reports, select_report, [('out', 'inlist')]),
        (compare_transforms, select_report, [('out', 'index')]),
        (select_report, outputnode, [('out', 'out_report')]),
    ])

    return workflow


def compare_xforms(lta_list, norm_threshold=15):
    """
    Computes a normalized displacement between two affine transforms as the
    maximum overall displacement of the midpoints of the faces of a cube, when
    each transform is applied to the cube.
    This combines displacement resulting from scaling, translation and rotation.

    Although the norm is in mm, in a scaling context, it is not necessarily
    equivalent to that distance in translation.

    We choose a default threshold of 15mm as a rough heuristic.
    Normalized displacement above 20mm showed clear signs of distortion, while
    "good" BBR refinements were frequently below 10mm displaced from the rigid
    transform.
    The 10-20mm range was more ambiguous, and 15mm chosen as a compromise.
    This is open to revisiting in either direction.

    See discussion in
    `GitHub issue #681`_ <https://github.com/poldracklab/aslprep/issues/681>`_
    and the `underlying implementation
    <https://github.com/nipy/nipype/blob/56b7c81eedeeae884ba47c80096a5f66bd9f8116/nipype/algorithms/rapidart.py#L108-L159>`_.

    Parameters
    ----------

      lta_list : list or tuple of str
          the two given affines in LTA format
      norm_threshold : float (default: 15)
          the upper bound limit to the normalized displacement caused by the
          second transform relative to the first

    """
    from ...niworkflows.interfaces.surf import load_transform
    from nipype.algorithms.rapidart import _calc_norm_affine

    bbr_affine = load_transform(lta_list[0])
    fallback_affine = load_transform(lta_list[1])

    norm, _ = _calc_norm_affine([fallback_affine, bbr_affine], use_differences=True)

    return norm[1] > norm_threshold
