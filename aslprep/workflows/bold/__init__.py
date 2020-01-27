# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
"""

Pre-processing fMRI - asl signal workflows
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. automodule:: fmriprep.workflows.asl.base
.. automodule:: fmriprep.workflows.asl.util
.. automodule:: fmriprep.workflows.asl.hmc
.. automodule:: fmriprep.workflows.asl.stc
.. automodule:: fmriprep.workflows.asl.t2s
.. automodule:: fmriprep.workflows.asl.registration
.. automodule:: fmriprep.workflows.asl.resampling
.. automodule:: fmriprep.workflows.asl.confounds


"""

from .base import init_asl_preproc_wf
from .util import init_bold_reference_wf
from .hmc import init_asl_hmc_wf
from .stc import init_asl_stc_wf
from .t2s import init_asl_t2s_wf
from .registration import (
    init_asl_t1_trans_wf,
    init_asl_reg_wf,
)
from .resampling import (
    init_asl_std_trans_wf,
    init_asl_surf_wf,
    init_asl_preproc_trans_wf,
)

from .confounds import (
    init_asl_confs_wf
)

__all__ = [
    'init_asl_confs_wf',
    'init_asl_hmc_wf',
    'init_asl_std_trans_wf',
    'init_asl_preproc_trans_wf',
    'init_bold_reference_wf',
    'init_asl_reg_wf',
    'init_asl_stc_wf',
    'init_asl_surf_wf',
    'init_asl_t1_trans_wf',
    'init_asl_t2s_wf',
    'init_asl_preproc_wf'
]
