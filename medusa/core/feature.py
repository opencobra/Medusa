
from __future__ import absolute_import

from cobra.core.object import Object

class Feature(Object):
    """
    Feature describing a network component that varies across an ensemble.

    Parameters
    ----------
    identifier : string
        The identifier to associate with the feature. Convention is to append the
        component_attribute to the base_component's id.

    ensemble : medusa.core.ensemble.Ensemble object
        The ensemble that the feature is associated with.

    base_component : cobra.core.reaction.Reaction
        Reference to the Reaction object that the feature describes.

    component_attribute : string
        string indicating the attribute of base_component that the feature
        describes the modification of (e.g. "lb", "ub")

    states : dictionary of string:component_attribute value
        dictionary of model ids mapping to the value of the Feature's
        component_attribute (value type depends on component_attribute type,
        e.g. float for "lb", string for "_gene_reaction_rule")

    Attributes
    ----------

    """

    def __init__(self, identifier=None, name=None,ensemble=None,\
                base_component=None, component_attribute=None, states=None):
        Object.__init__(self,identifier,name)
        self.ensemble = ensemble
        self.base_component = base_component
        self.component_attribute = component_attribute
        self.states = states


    def get_model_state(self,member_id):
        """Get the state of the feature for a particular member
        """
        return self.states[member_id]
