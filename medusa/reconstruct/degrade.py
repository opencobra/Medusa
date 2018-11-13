
from __future__ import absolute_import







# functions for degrading networks to construct ensembles

def degrade_reactions(base_model,num_reactions,num_models=10):
    """
    Removes reactions from an existing COBRA model to generate an ensemble.

    Parameters
    ----------
    base_model: cobra.Model
        Model from which reactions will be removed to create an ensemble.
    num_reactions: int
        The number of reactions to remove to generate each ensemble member.
        Must be smaller than the total number of reactions in the model
    num_models: int
        The number of models to generate by randomly removing num_reactions from
        the base_model. Reactions are removed with replacement.

    Returns
    -------
    Medusa.core.ensemble
        An ensemble
    """
    filler=1
