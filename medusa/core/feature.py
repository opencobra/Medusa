

class Feature(Object):
    """
    Feature describing a network component that varies across an ensemble.

    Parameters
    ----------
    identifier : string
        The identifier to associate with the feature. Default is to append a
        feature count to the reaction ID starting the count at 0 for the first
        feature.

    ensemble : medusa.core.ensemble.Ensemble object
        Either a list of existing Model objects in which case a new Model
        object is instantiated and an ensemble is constructed using the list of
        Models, or None/empty list, in which case an ensemble is created with
        empty attributes.

    base_component : cobra.core.reaction.Reaction

    Attributes
    ----------
    base_model : Model
        A cobra.core.model.Model that contains all variable and invariable components
        of an ensemble.
    members : DictList
        A DictList where the key is the member identifier and the value is a
        medusa.core.member.Member object
    features : DictList
        A DictList where the key is the feature identifier and the value is a
        medusa.core.feature.Feature object
    """

    def __init__(self, identifier=None, ensemble=None, base_component=None,\
                    component_attribute=None, states=None):
