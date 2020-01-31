# -*- coding: utf-8 -*-
# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
"""
Head-Motion Estimation and Correction (HMC) of asl images
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autofunction:: init_asl_hmc_wf

"""

from nipype.pipeline import engine as pe
from nipype.interfaces import utility as niu, fsl
from ...niworkflows.niworkflows.engine.workflows import LiterateWorkflow as Workflow
from ...niworkflows.niworkflows.interfaces import NormalizeMotionParams
from ...niworkflows.niworkflows.interfaces.itk import MCFLIRT2ITK
from ...config import DEFAULT_MEMORY_MIN_GB


def init_asl_hmc_wf(mem_gb, omp_nthreads, name='asl_hmc_wf'):
    """
    Build a workflow to estimate head-motion parameters.

    This workflow estimates the motion parameters to perform
    :abbr:`HMC (head motion correction)` over the input
    :abbr:`asl (blood-oxygen-level dependent)` image.

    Workflow Graph
        .. workflow::
            :graph2use: orig
            :simple_form: yes

            from fmriprep.workflows.asl import init_asl_hmc_wf
            wf = init_asl_hmc_wf(
                mem_gb=3,
                omp_nthreads=1)

    Parameters
    ----------
    mem_gb : float
        Size of asl file in GB
    omp_nthreads : int
        Maximum number of threads an individual process may use
    name : str
        Name of workflow (default: ``asl_hmc_wf``)

    Inputs
    ------
    asl_file
        asl series NIfTI file
    raw_ref_image
        Reference image to which asl series is motion corrected

    Outputs
    -------
    xforms
        ITKTransform file aligning each volume to ``ref_image``
    movpar_file
        MCFLIRT motion parameters, normalized to SPM format (X, Y, Z, Rx, Ry, Rz)

    """
    workflow = Workflow(name=name)
    workflow.__desc__ = """\
Head-motion parameters with respect to the asl reference
(transformation matrices, and six corresponding rotation and translation
parameters) are estimated before any spatiotemporal filtering using
`mcflirt` [FSL {fsl_ver}, @mcflirt].
""".format(fsl_ver=fsl.Info().version() or '<ver>')

    inputnode = pe.Node(niu.IdentityInterface(fields=['asl_file', 'raw_ref_image']),
                        name='inputnode')
    outputnode = pe.Node(
        niu.IdentityInterface(fields=['xforms', 'movpar_file']),
        name='outputnode')

    # Head motion correction (hmc)
    mcflirt = pe.Node(fsl.MCFLIRT(save_mats=True, save_plots=True),
                      name='mcflirt', mem_gb=mem_gb * 3)

    fsl2itk = pe.Node(MCFLIRT2ITK(), name='fsl2itk',
                      mem_gb=0.05, n_procs=omp_nthreads)

    normalize_motion = pe.Node(NormalizeMotionParams(format='FSL'),
                               name="normalize_motion",
                               mem_gb=DEFAULT_MEMORY_MIN_GB)

    workflow.connect([
        (inputnode, mcflirt, [('raw_ref_image', 'ref_file'),
                              ('asl_file', 'in_file')]),
        (inputnode, fsl2itk, [('raw_ref_image', 'in_source'),
                              ('raw_ref_image', 'in_reference')]),
        (mcflirt, fsl2itk, [('mat_file', 'in_files')]),
        (mcflirt, normalize_motion, [('par_file', 'in_file')]),
        (fsl2itk, outputnode, [('out_file', 'xforms')]),
        (normalize_motion, outputnode, [('out_file', 'movpar_file')]),
    ])

    return workflow
