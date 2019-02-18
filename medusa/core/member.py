
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

    def to_model(self):
        """
        Generate a cobra.Model object with the exact state of this member.

        The resulting cobra.Model does not contain any Metabolites, Genes,
        or Reactions that were inactive in the member.

        Returns
        -------
        model : cobra.Model
            The extracted member as a cobrapy model.
        """

        # Set the state of the ensemble.base_model to represent this member
        self.ensemble.set_state(self.id)

        # copy the base model and remove any inactive reactions and
        # associated genes/metabolites
        model = self.ensemble.base_model.copy()
        inactive_rxns = [rxn for rxn in model.reactions if rxn.bounds == (0,0)]
        model.remove_reactions(inactive_rxns, remove_orphans = True)

        return model
