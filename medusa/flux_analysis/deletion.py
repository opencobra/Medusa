from cobra.flux_analysis import single_reaction_deletion, single_gene_deletion

from medusa.flux_analysis._selection import (resolve_member_ids,
                                             truncated_id_list)


def _resolve_specific_genes(ensemble, specific_genes):
    """Validate and normalize the required `specific_genes` argument.

    `specific_genes` is required rather than defaulting to "all genes".
    Previously it defaulted to an empty list, which was forwarded to cobrapy's
    ``gene_list``; cobrapy treats only ``None`` as "every gene" and an empty
    list as "no genes", so the default call deleted nothing and returned an
    empty result for every member without warning.
    """
    if specific_genes is None:
        raise ValueError(
            "specific_genes is required. It previously defaulted to an empty "
            "list, which cobrapy interprets as 'delete no genes', so the "
            "default call silently returned no results. Pass the gene ids to "
            "delete, or every gene with "
            "specific_genes=[gene.id for gene in ensemble.base_model.genes].")

    if isinstance(specific_genes, str):
        # A bare '' is an empty selection, not a gene named ''.
        specific_genes = [specific_genes] if specific_genes else []
    else:
        try:
            specific_genes = list(specific_genes)
        except TypeError:
            raise TypeError(
                "specific_genes must be a gene id, a Gene, or an iterable of "
                "them; got %r." % (specific_genes,)) from None

    if not specific_genes:
        raise ValueError(
            "specific_genes is empty, which would delete no genes and return "
            "an empty result for every member. Pass the gene ids to delete, "
            "or every gene with "
            "specific_genes=[gene.id for gene in ensemble.base_model.genes].")

    gene_ids = []
    for gene in specific_genes:
        gene_id = gene if isinstance(gene, str) else getattr(gene, 'id', None)
        if not gene_id:
            raise ValueError(
                "specific_genes contains an empty or unidentifiable entry: "
                "%r." % (gene,))
        gene_ids.append(gene_id)

    known = {gene.id for gene in ensemble.base_model.genes}
    unknown = [gene_id for gene_id in gene_ids if gene_id not in known]
    if unknown:
        raise KeyError(
            "specific_genes refers to genes that are not in the ensemble's "
            "base model: %s" % truncated_id_list(unknown))

    return gene_ids


def ensemble_single_reaction_deletion(ensemble, num_models=None,
                                        specific_models=None):
    '''
    Performs single reaction deletions on models within an ensemble and
    returns the objective value after optimization with each reaction removed.

    Parameters
    ----------
    ensemble: medusa.core.Ensemble
        The ensemble with which to perform reaction deletions
    num_models: int, optional
        Number of models for which reaction deletions will be performed. The
        number of models indicated will be randomly sampled and reaction
        deletions will be performed on the sampled models. If None, all models
        will be selected (default), or the models specified by specific_models
        will be selected. Cannot be passed concurrently with specific_models.
    specific_models: list of str or medusa.core.member.Member, optional
        The models for which reaction deletions will be performed, given as
        member ids or Member objects. If None, all models will be selected
        (default), or num_models will be randomly sampled and selected.
        Cannot be passed concurrently with num_models.

    Returns
    -------
    dict
        A dict of {member_id: pandas.DataFrame}, where each DataFrame is the
        cobrapy single_reaction_deletion result for that member.
    '''
    model_list = resolve_member_ids(ensemble, specific_models, num_models)

    deletion_results = {}
    with ensemble.base_model:
        for member_id in model_list:
            print('performing deletions for ' + member_id)
            ensemble.set_state(member_id)
            deletion_result = single_reaction_deletion(ensemble.base_model)
            deletion_results[member_id] = deletion_result

    return deletion_results


def ensemble_single_gene_deletion(ensemble, num_models=None,
                                        specific_models=None,
                                        specific_genes=None):
    '''
    Performs single gene deletions on models within an ensemble and
    returns the objective value after optimization with each gene removed.

    Parameters
    ----------
    ensemble: medusa.core.Ensemble
        The ensemble with which to perform gene deletions
    num_models: int, optional
        Number of models for which gene deletions will be performed. The
        number of models indicated will be randomly sampled and gene
        deletions will be performed on the sampled models. If None, all models
        will be selected (default), or the models specified by specific_models
        will be selected. Cannot be passed concurrently with specific_models.
    specific_models: list of str or medusa.core.member.Member, optional
        The models for which gene deletions will be performed, given as member
        ids or Member objects. If None, all models will be selected (default),
        or num_models will be randomly sampled and selected. Cannot be passed
        concurrently with num_models.
    specific_genes: list of str or cobra.core.gene.Gene
        The genes to delete, given as gene ids or Gene objects. A single id or
        Gene may be passed directly. This argument is REQUIRED: passing None,
        an empty list, or an empty string raises ValueError rather than
        silently deleting nothing. To delete every gene, pass
        ``[gene.id for gene in ensemble.base_model.genes]``.

        We recommend identifying genes that are essential in all ensemble
        members first, then excluding those genes from specific_genes. This
        will generally speed up computation.

    Returns
    -------
    dict
        A dict of {member_id: pandas.DataFrame}, where each DataFrame is the
        cobrapy single_gene_deletion result for that member.

    Raises
    ------
    ValueError
        If specific_genes is not given, or is an empty selection.
    KeyError
        If specific_genes names a gene that is not in the base model.
    '''
    gene_ids = _resolve_specific_genes(ensemble, specific_genes)
    model_list = resolve_member_ids(ensemble, specific_models, num_models)

    deletion_results = {}
    with ensemble.base_model:
        for member_id in model_list:
            print('performing deletions for ' + member_id)
            ensemble.set_state(member_id)
            deletion_result = single_gene_deletion(ensemble.base_model,
                                                   gene_ids)
            deletion_results[member_id] = deletion_result

    return deletion_results
