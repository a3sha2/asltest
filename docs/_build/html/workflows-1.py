from collections import namedtuple
from aslprep.niworkflows.utils.spaces import Reference, SpatialReferences
from aslprep.workflows.base import init_single_subject_wf
BIDSLayout = namedtuple('BIDSLayout', ('root'))
wf = init_single_subject_wf(
    anat_only=False,
    bold2t1w_dof=9,
    cifti_output=False,
    debug=False,
    dummy_scans=None,
    echo_idx=None,
    fmap_bspline=False,
    fmap_demean=True,
    force_syn=True,
    freesurfer=True,
    hires=True,
    ignore=[],
    layout=BIDSLayout('.'),
    longitudinal=False,
    low_mem=False,
    medial_surface_nan=False,
    name='single_subject_wf',
    omp_nthreads=1,
    output_dir='.',
    reportlets_dir='.',
    skull_strip_fixed_seed=False,
    skull_strip_template=Reference('OASIS30ANTs'),
    spaces=SpatialReferences(
        spaces=[('MNI152Lin', {}),
                ('fsaverage', {'density': '10k'}),
                ('T1w', {}),
                ('fsnative', {})],
        checkpoint=True),
    subject_id='test',
    t2s_coreg=False,
    task_id='',
    use_bbr=True,
    use_syn=True,
    bids_filters=None,
    pcasl=pcasl,
)