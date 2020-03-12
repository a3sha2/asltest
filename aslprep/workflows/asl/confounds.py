# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
"""
Calculate ASL confounds
^^^^^^^^^^^^^^^^^^^^^^^^

.. autofunction:: init_asl_confs_wf

"""
from os import getenv
from nipype.pipeline import engine as pe
from nipype.interfaces import utility as niu, fsl
from nipype.algorithms import confounds as nac

from templateflow.api import get as get_template
from ...niworkflows.engine.workflows import LiterateWorkflow as Workflow
from ...niworkflows.interfaces.confounds import ExpandModel, SpikeRegressors
from ...niworkflows.interfaces.fixes import FixHeaderApplyTransforms as ApplyTransforms
from ...niworkflows.interfaces.images import SignalExtraction
from ...niworkflows.interfaces.masks import ROIsPlot
from ...niworkflows.interfaces.utility import KeySelect



from ...niworkflows.interfaces.utils import (
    TPM2ROI, AddTPMs, AddTSVHeader, TSV2JSON, DictMerge
)

from ...config import DEFAULT_MEMORY_MIN_GB
from ...interfaces import (
    GatherConfounds,ASLSummary, DerivativesDataSink
)


def init_asl_confs_wf(
    mem_gb,
    metadata,
    name="asl_confs_wf",
):
    """
    Build a workflow to generate and write out confounding signals.

    This workflow calculates confounds for a BOLD series, and aggregates them
    into a :abbr:`TSV (tab-separated value)` file, for use as nuisance
    regressors in a :abbr:`GLM (general linear model)`.
    The following confounds are calculated, with column headings in parentheses:

    #. Region-wise average signal (``csf``, ``white_matter``, ``global_signal``)
    #. DVARS - original and standardized variants (``dvars``, ``std_dvars``)
    #. Framewise displacement, based on head-motion parameters
       (``framewise_displacement``)
    #. Temporal CompCor (``t_comp_cor_XX``)
    #. Anatomical CompCor (``a_comp_cor_XX``)
    #. Cosine basis set for high-pass filtering w/ 0.008 Hz cut-off
       (``cosine_XX``)
    #. Non-steady-state volumes (``non_steady_state_XX``)
    #. Estimated head-motion parameters, in mm and rad
       (``trans_x``, ``trans_y``, ``trans_z``, ``rot_x``, ``rot_y``, ``rot_z``)


    Prior to estimating aCompCor and tCompCor, non-steady-state volumes are
    censored and high-pass filtered using a :abbr:`DCT (discrete cosine
    transform)` basis.
    The cosine basis, as well as one regressor per censored volume, are included
    for convenience.

    Workflow Graph
        .. workflow::
            :graph2use: orig
            :simple_form: yes

            from fmriprep.workflows.bold.confounds import init_bold_confs_wf
            wf = init_bold_confs_wf(
                mem_gb=1,
                metadata={},
                regressors_all_comps=False,
                regressors_dvars_th=1.5,
                regressors_fd_th=0.5,
            )

    Parameters
    ----------
    mem_gb : float
        Size of BOLD file in GB - please note that this size
        should be calculated after resamplings that may extend
        the FoV
    metadata : dict
        BIDS metadata for BOLD file
    name : str
        Name of workflow (default: ``bold_confs_wf``)
    regressors_all_comps: bool
        Indicates whether CompCor decompositions should return all
        components instead of the minimal number of components necessary
        to explain 50 percent of the variance in the decomposition mask.
    regressors_dvars_th
        Criterion for flagging DVARS outliers
    regressors_fd_th
        Criterion for flagging framewise displacement outliers

    Inputs
    ------
    bold
        BOLD image, after the prescribed corrections (STC, HMC and SDC)
        when available.
    bold_mask
        BOLD series mask
    movpar_file
        SPM-formatted motion parameters file
    skip_vols
        number of non steady state volumes
    t1w_mask
        Mask of the skull-stripped template image
    t1w_tpms
        List of tissue probability maps in T1w space
    t1_bold_xform
        Affine matrix that maps the T1w space into alignment with
        the native BOLD space

    Outputs
    -------
    confounds_file
        TSV of all aggregated confounds
    rois_report
        Reportlet visualizing white-matter/CSF mask used for aCompCor,
        the ROI for tCompCor and the BOLD brain mask.
    confounds_metadata
        Confounds metadata dictionary.

    """
    workflow = Workflow(name=name)
    workflow.__desc__ = """\
Several confounding time-series were calculated based on the
*preprocessed BOLD*: framewise displacement (FD), DVARS and
three region-wise global signals.
FD and DVARS are calculated for each functional run, both using their
implementations in *Nipype* [following the definitions by @power_fd_dvars].
The three global signals are extracted within the CSF, the WM, and
the whole-brain masks.
Additionally, a set of physiological regressors were extracted to
allow for component-based noise correction [*CompCor*, @compcor].
Principal components are estimated after high-pass filtering the
*preprocessed BOLD* time-series (using a discrete cosine filter with
128s cut-off) for the two *CompCor* variants: temporal (tCompCor)
and anatomical (aCompCor).
tCompCor components are then calculated from the top 5% variable
voxels within a mask covering the subcortical regions.
This subcortical mask is obtained by heavily eroding the brain mask,
which ensures it does not include cortical GM regions.
For aCompCor, components are calculated within the intersection of
the aforementioned mask and the union of CSF and WM masks calculated
in T1w space, after their projection to the native space of each
functional run (using the inverse BOLD-to-T1w transformation). Components
are also calculated separately within the WM and CSF masks.
For each CompCor decomposition, the *k* components with the largest singular
values are retained, such that the retained components' time series are
sufficient to explain 50 percent of variance across the nuisance mask (CSF,
WM, combined, or temporal). The remaining components are dropped from
consideration.
The head-motion estimates calculated in the correction step were also
placed within the corresponding confounds file.
The confound time series derived from head motion estimates and global
signals were expanded with the inclusion of temporal derivatives and
quadratic terms for each [@confounds_satterthwaite_2013].
Frames that exceeded a threshold of  mm FD or  standardised DVARS
were annotated as motion outliers.
"""
    inputnode = pe.Node(niu.IdentityInterface(
        fields=['asl', 'asl_mask', 'movpar_file', 'skip_vols',
                't1w_mask', 't1w_tpms', 't1_bold_xform']),
        name='inputnode')
    outputnode = pe.Node(niu.IdentityInterface(
        fields=['confounds_file', 'confounds_metadata']),
        name='outputnode')

    # Get masks ready in T1w space
    #acc_tpm = pe.Node(AddTPMs(indices=[0, 2]), name='tpms_add_csf_wm')  # acc stands for aCompCor
    #csf_roi = pe.Node(TPM2ROI(erode_mm=0, mask_erode_mm=30), name='csf_roi')
    #wm_roi = pe.Node(TPM2ROI(
        #erode_prop=0.6, mask_erode_prop=0.6**3),  # 0.6 = radius; 0.6^3 = volume
        #name='wm_roi')
    #acc_roi = pe.Node(TPM2ROI(
        #erode_prop=0.6, mask_erode_prop=0.6**3),  # 0.6 = radius; 0.6^3 = volume
        #name='acc_roi')

    # Map ROIs in T1w space into BOLD space
    #csf_tfm = pe.Node(ApplyTransforms(interpolation='NearestNeighbor', float=True),
                      #name='csf_tfm', mem_gb=0.1)
    #wm_tfm = pe.Node(ApplyTransforms(interpolation='NearestNeighbor', float=True),
                     #name='wm_tfm', mem_gb=0.1)
    #acc_tfm = pe.Node(ApplyTransforms(interpolation='NearestNeighbor', float=True),
                      #name='acc_tfm', mem_gb=0.1)
    #tcc_tfm = pe.Node(ApplyTransforms(interpolation='NearestNeighbor', float=True),
                      #name='tcc_tfm', mem_gb=0.1)

    # Ensure ROIs don't go off-limits (reduced FoV)
    #csf_msk = pe.Node(niu.Function(function=_maskroi), name='csf_msk')
    #wm_msk = pe.Node(niu.Function(function=_maskroi), name='wm_msk')

    #acc_msk = pe.Node(niu.Function(function=_maskroi), name='acc_msk')
    #tcc_msk = pe.Node(niu.Function(function=_maskroi), name='tcc_msk')

    # DVARS
    dvars = pe.Node(nac.ComputeDVARS(save_nstd=True, save_std=True, remove_zerovariance=True),
                    name="dvars", mem_gb=mem_gb)

    # Frame displacement
    fdisp = pe.Node(nac.FramewiseDisplacement(parameter_source="SPM"),
                    name="fdisp", mem_gb=mem_gb)

    # a/t-CompCor
    #rg_lbl_cc = pe.Node(niu.Merge(3), name='merge_rois_cc', run_without_submitting=True)

    #tcompcor = pe.Node(
        #TCompCor(components_file='tcompcor.tsv', header_prefix='t_comp_cor_', pre_filter='cosine',
                 #save_pre_filter=True, save_metadata=True, percentile_threshold=.05,
                 #failure_mode='NaN'),
        #name="tcompcor", mem_gb=mem_gb)

    #acompcor = pe.Node(
        #ACompCor(components_file='acompcor.tsv', header_prefix='a_comp_cor_', pre_filter='cosine',
                 #save_pre_filter=True, save_metadata=True, mask_names=['combined', 'CSF', 'WM'],
                 #merge_method='none', failure_mode='NaN'),
        #name="acompcor", mem_gb=mem_gb)

    # Set number of components
    

    # Set TR if present
    
    # Global and segment regressors
    #signals_class_labels = ["csf", "white_matter", "global_signal"]
    #mrg_lbl = pe.Node(niu.Merge(3), name='merge_rois', run_without_submitting=True)
    #signals = pe.Node(SignalExtraction(class_labels=signals_class_labels),
                      #name="signals", mem_gb=mem_gb)

    # Arrange confounds
    add_dvars_header = pe.Node(
        AddTSVHeader(columns=["dvars"]),
        name="add_dvars_header", mem_gb=0.01, run_without_submitting=True)
    add_std_dvars_header = pe.Node(
        AddTSVHeader(columns=["std_dvars"]),
        name="add_std_dvars_header", mem_gb=0.01, run_without_submitting=True)
    add_motion_headers = pe.Node(
        AddTSVHeader(columns=["trans_x", "trans_y", "trans_z", "rot_x", "rot_y", "rot_z"]),
        name="add_motion_headers", mem_gb=0.01, run_without_submitting=True)
    concat = pe.Node(GatherConfounds(), name="concat", mem_gb=0.01, run_without_submitting=True)

    # CompCor metadata
    

    # Expand model to include derivatives and quadratics
    #model_expand = pe.Node(ExpandModel(
        #model_formula='(dd1(rps + wm + csf + gsr))^^2 + others'),
        #name='model_expansion')

    # Add spike regressors
    #spike_regress = pe.Node(SpikeRegressors(
        #fd_thresh=regressors_fd_th,
        #dvars_thresh=regressors_dvars_th),
        #name='spike_regressors')

    # Generate reportlet (ROIs)
    #mrg_compcor = pe.Node(niu.Merge(2), name='merge_compcor', run_without_submitting=True)
    rois_plot = pe.Node(ROIsPlot(colors=['b', 'magenta'], generate_report=True),
                        name='rois_plot', mem_gb=mem_gb)

    #ds_report_bold_rois = pe.Node(
        #DerivativesDataSink(desc='rois', keep_dtype=True),
        #name='ds_report_bold_rois', run_without_submitting=True,
        #mem_gb=DEFAULT_MEMORY_MIN_GB)

    # Generate reportlet (CompCor)
    #mrg_cc_metadata = pe.Node(niu.Merge(2), name='merge_compcor_metadata',
                              #run_without_submitting=True)
    #compcor_plot = pe.Node(
        #CompCorVariancePlot(variance_thresholds=(0.5, 0.7, 0.9),
                            #metadata_sources=['tCompCor', 'aCompCor']),
        #name='compcor_plot')
    #ds_report_compcor = pe.Node(
        #DerivativesDataSink(desc='compcorvar', keep_dtype=True),
        #name='ds_report_compcor', run_without_submitting=True,
        #mem_gb=DEFAULT_MEMORY_MIN_GB)

    # Generate reportlet (Confound correlation)
    #conf_corr_plot = pe.Node(
        #ConfoundsCorrelationPlot(reference_column='global_signal', max_dim=70),
        #name='conf_corr_plot')
    #ds_report_conf_corr = pe.Node(
        #DerivativesDataSink(desc='confoundcorr', keep_dtype=True),
        #name='ds_report_conf_corr', run_without_submitting=True,
        #mem_gb=DEFAULT_MEMORY_MIN_GB)

    #def _pick_csf(files):
        #return files[0]

    #def _pick_wm(files):
        #return files[-1]

    workflow.connect([
        # connect inputnode to each non-anatomical confound node
        (inputnode, dvars, [('asl', 'in_file'),
                            ('asl_mask', 'in_mask')]),
        (inputnode, fdisp, [('movpar_file', 'in_file')]),
        # Collate computed confounds together
        (inputnode, add_motion_headers, [('movpar_file', 'in_file')]),
        (dvars, add_dvars_header, [('out_nstd', 'in_file')]),
        (dvars, add_std_dvars_header, [('out_std', 'in_file')]),
        (fdisp, concat, [('out_file', 'fd')]),
        (add_motion_headers, concat, [('out_file', 'motion')]),
        (add_dvars_header, concat, [('out_file', 'dvars')]),
        (add_std_dvars_header, concat, [('out_file', 'std_dvars')]),
        # Set outputs
        (concat, outputnode, [('out_dict', 'confounds_metadata')]),
        (inputnode, rois_plot, [('asl', 'in_file'),
                                ('asl_mask', 'in_mask')]),
     ])

    return workflow


def init_carpetplot_wf(standard_spaces, mem_gb, metadata, name="asl_carpet_wf"):
    """
    Build a workflow to generate *carpet plots*.

    Resamples the MNI parcellation (ad-hoc parcellation derived from the
    Harvard-Oxford template and others).

    Parameters
    ----------
    mem_gb : float
        Size of BOLD file in GB - please note that this size
        should be calculated after resamplings that may extend
        the FoV
    metadata : dict
        BIDS metadata for BOLD file
    name : str
        Name of workflow (default: ``bold_carpet_wf``)

    Inputs
    ------
    bold
        BOLD image, after the prescribed corrections (STC, HMC and SDC)
        when available.
    bold_mask
        BOLD series mask
    confounds_file
        TSV of all aggregated confounds
    t1_bold_xform
        Affine matrix that maps the T1w space into alignment with
        the native BOLD space
    std2anat_xfm
        ANTs-compatible affine-and-warp transform file

    Outputs
    -------
    out_carpetplot
        Path of the generated SVG file

    """
    inputnode = pe.Node(niu.IdentityInterface(
        fields=['asl', 'asl_mask', 'confounds_file',
                't1_bold_xform', 'std2anat_xfm']),
        name='inputnode')

    outputnode = pe.Node(niu.IdentityInterface(
        fields=['out_carpetplot']), name='outputnode')

    select_std = pe.Node(KeySelect(
        keys=list(standard_spaces.keys()), fields=['std2anat_xfm']),
        name='select_std', run_without_submitting=True)
    select_std.inputs.key = 'MNI152NLin2009cAsym'

    # List transforms
    mrg_xfms = pe.Node(niu.Merge(2), name='mrg_xfms')

    # Warp segmentation into EPI space
    resample_parc = pe.Node(ApplyTransforms(
        float=True,
        input_image=str(get_template(
            'MNI152NLin2009cAsym', resolution=1, desc='carpet',
            suffix='dseg', extension=['.nii', '.nii.gz'])),
        dimension=3, default_value=0, interpolation='MultiLabel'),
        name='resample_parc')

    # Carpetplot and confounds plot
    conf_plot = pe.Node(ASLSummary(
        tr=metadata['RepetitionTime'],
        confounds_list=[
            ('std_dvars', None, 'DVARS'),
            ('framewise_displacement', 'mm', 'FD')]),
        name='conf_plot', mem_gb=mem_gb)
    ds_report_asl_conf = pe.Node(
        DerivativesDataSink(desc='carpetplot', keep_dtype=True),
        name='ds_report_asl_conf', run_without_submitting=True,
        mem_gb=DEFAULT_MEMORY_MIN_GB)

    workflow = Workflow(name=name)
    workflow.connect([
        (inputnode, select_std, [('std2anat_xfm', 'std2anat_xfm')]),
        (inputnode, mrg_xfms, [('t1_asl_xform', 'in1')]),
        (inputnode, resample_parc, [('asl_mask', 'reference_image')]),
        (select_std, mrg_xfms, [('std2anat_xfm', 'in2')]),
        (mrg_xfms, resample_parc, [('out', 'transforms')]),
        # Carpetplot
        (inputnode, conf_plot, [
            ('asl', 'in_func'),
            ('asl_mask', 'in_mask'),
            ('confounds_file', 'confounds_file')]),
        (resample_parc, conf_plot, [('output_image', 'in_segm')]),
        (conf_plot, ds_report_asl_conf, [('out_file', 'in_file')]),
        (conf_plot, outputnode, [('out_file', 'out_carpetplot')]),
    ])
    return workflow


def _add_volumes(bold_file, bold_cut_file, skip_vols):
    """Prepend skip_vols from bold_file onto bold_cut_file."""
    import nibabel as nb
    import numpy as np
    from nipype.utils.filemanip import fname_presuffix

    if skip_vols == 0:
        return bold_cut_file

    bold_img = nb.load(bold_file)
    bold_cut_img = nb.load(bold_cut_file)

    bold_data = np.concatenate((bold_img.dataobj[..., :skip_vols],
                                bold_cut_img.dataobj), axis=3)

    out = fname_presuffix(bold_cut_file, suffix='_addnonsteady')
    bold_img.__class__(bold_data, bold_img.affine, bold_img.header).to_filename(out)

    return out


def _maskroi(in_mask, roi_file):
    import numpy as np
    import nibabel as nb
    from nipype.utils.filemanip import fname_presuffix

    roi = nb.load(roi_file)
    roidata = roi.get_data().astype(np.uint8)
    msk = nb.load(in_mask).get_data().astype(bool)
    roidata[~msk] = 0
    roi.set_data_dtype(np.uint8)

    out = fname_presuffix(roi_file, suffix='_aslmsk')
    roi.__class__(roidata, roi.affine, roi.header).to_filename(out)
    return out
