# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:

# Load modules for compatibility
from niworkflows.interfaces import (
    bids, cifti, freesurfer, images, itk, surf, utils)

from .reports import SubjectSummary, FunctionalSummary, AboutSummary
from .fmap import FieldEnhance, FieldToRadS, FieldToHz, Phasediff2Fieldmap
from .confounds import GatherConfounds, ICAConfounds, ASLSummary
from .multiecho import T2SMap


class DerivativesDataSink(bids.DerivativesDataSink):
    out_path_base = 'aslprep'


__all__ = [
    'bids',
    'images',
    'itk',
    'utils',
    'SubjectSummary',
    'FunctionalSummary',
    'AboutSummary',
    'FieldEnhance',
    'FieldToRadS',
    'FieldToHz',
    'Phasediff2Fieldmap',
    'GatherConfounds',
    'ASLSummary',
    'T2SMap',
    'DerivativesDataSink',
]
