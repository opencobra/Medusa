"""Shared helpers for selecting ensemble members to operate on.

These live in one place so that every entry point agrees on what
``specific_models`` and ``num_models`` mean. Historically they did not: the
deletion functions read ``model.id`` and so required Member objects, while
``ensemble_fva`` concatenated the value into a string and so required member
ids, meaning no single value worked across all three.
"""
from __future__ import absolute_import

from random import sample

from medusa.core.member import Member

# Number of ids listed in an error message before truncating.
MAX_IDS_IN_MESSAGE = 10


def truncated_id_list(ids):
    """Format ids for an error message, capping the number shown."""
    ids = list(ids)
    shown = ids[:MAX_IDS_IN_MESSAGE]
    message = ', '.join(repr(i) for i in shown)
    if len(ids) > len(shown):
        message += ' ... (%i total)' % len(ids)
    return message


def resolve_member_ids(ensemble, specific_models=None, num_models=None):
    """Return the member ids to operate on, as a list of str.

    Parameters
    ----------
    ensemble : medusa.core.ensemble.Ensemble
    specific_models : iterable of str or medusa.core.member.Member, optional
        The members to select, given as member ids or Member objects. A
        single id or Member may be passed directly. If None, selection falls
        back to num_models.
    num_models : int, optional
        Randomly sample this many members. If None (and specific_models is
        also None), every member is selected.

    Returns
    -------
    list of str
    """
    if specific_models is not None and num_models is not None:
        raise ValueError(
            "specific_models and num_models cannot be passed concurrently; "
            "pass one or neither.")

    if specific_models is not None:
        if isinstance(specific_models, (str, Member)):
            specific_models = [specific_models]
        model_list = [member.id if isinstance(member, Member) else str(member)
                      for member in specific_models]
        if not model_list:
            raise ValueError(
                "specific_models is empty. Pass at least one member id, or "
                "leave it as None to use every member of the ensemble.")
        known = {member.id for member in ensemble.members}
        unknown = [member_id for member_id in model_list
                   if member_id not in known]
        if unknown:
            raise KeyError(
                "specific_models refers to ids that are not members of this "
                "ensemble: %s" % truncated_id_list(unknown))
        return model_list

    all_ids = [member.id for member in ensemble.members]
    if num_models is None:
        return all_ids
    if num_models < 1:
        raise ValueError(
            "num_models must be at least 1; got %r." % (num_models,))
    if num_models >= len(all_ids):
        return all_ids
    return sample(all_ids, num_models)
