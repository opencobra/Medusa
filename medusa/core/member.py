
from __future__ import absolute_import

from cobra.core.object import Object

class Member(Object):
    """
    Object representing an individual member (i.e. model) in an ensemble

    Parameters
    ----------
    identifier : string
        The identifier to associate with the member.

    ensemble : medusa.core.ensemble.Ensemble object
        The ensemble that the member belongs to.

    states : dictionary of medusa.core.feature.Feature:component_attribute value
        dictionary of Features mapping to the value of the Feature's
        component_attribute (value type depends on component_attribute type,
        e.g. float for "lb", string for "_gene_reaction_rule") for the member.

    Attributes
    ----------

    """

    def __init__(self,ensemble=None, identifier=None, name=None, states=None):
        Object.__init__(self,identifier,name)
        self.ensemble = ensemble
        self.states = states
