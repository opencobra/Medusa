import pandas as pd

from cobra.test import create_test_model
from cobra.io import load_json_model
from cobra.core import Reaction

from medusa.core.ensemble import Ensemble
from medusa.reconstruct import expand

REACTION_ATTRIBUTES = ['lower_bound', 'upper_bound']
MISSING_ATTRIBUTE_DEFAULT = {'lower_bound':0,'upper_bound':0}



def test_iterative_gapfill_from_binary_phenotypes():
    # load the universal model and a test model
    universal = load_universal_modelseed()
    model = load_modelseed_model('Staphylococcus aureus')

    # Load the biolog composition to be used for gapfilling
    biolog_base_composition = pd.read_csv('./medusa/test/data/biolog_base_composition.csv',sep=',')
    biolog_base_dict = dict(zip(biolog_base_composition['ID'],\
                              [1000 for i in range(0,len(biolog_base_composition['ID']))]))
    biolog_thresholded = pd.read_csv('./medusa/test/data/plata_thresholded.csv',sep='\t',index_col=0)

    # extract the biolog conditions for Staphylococcus aureus
    test_mod_pheno = biolog_thresholded.loc['Staphylococcus aureus']
    test_mod_pheno = list(test_mod_pheno[test_mod_pheno == True].index)

    # check for biolog base components in the model
    add_mets = []
    add_exchanges = []
    for met in list(biolog_base_dict.keys()):
        try:
            model.metabolites.get_by_id(met)
        except:
            print('no '+met)
            add_met = universal.metabolites.get_by_id(met).copy()
            add_mets.append(add_met)

    model.add_metabolites(add_mets)

    for met in list(biolog_base_dict.keys()):
        # Search for exchange reactions
        try:
            model.reactions.get_by_id('EX_'+met)
        except:
            add_met = model.metabolites.get_by_id(met)
            ex_rxn = Reaction('EX_' + met)
            ex_rxn.name = "Exchange reaction for " + met
            ex_rxn.lower_bound = -1000
            ex_rxn.upper_bound = 1000
            ex_rxn.add_metabolites({add_met:-1})
            add_exchanges.append(ex_rxn)

    model.add_reactions(add_exchanges)

    # Find metabolites from the biolog data that are missing in the test model
    # and add them from the universal
    missing_mets = []
    missing_exchanges = []
    media_dicts = {}
    for met_id in test_mod_pheno:
        try:
            model.metabolites.get_by_id(met_id)
        except:
            print(met_id + " was not in model, adding met and exchange reaction")
            met = universal.metabolites.get_by_id(met_id).copy()
            missing_mets.append(met)
            ex_rxn = Reaction('EX_' + met_id)
            ex_rxn.name = "Exchange reaction for " + met_id
            ex_rxn.lower_bound = -1000
            ex_rxn.upper_bound = 1000
            ex_rxn.add_metabolites({met:-1})
            missing_exchanges.append(ex_rxn)
        media_dicts[met_id] = biolog_base_dict.copy()
        media_dicts[met_id] = {'EX_'+k:v for k,v in media_dicts[met_id].items()}
        media_dicts[met_id]['EX_'+met_id] = 1000
    model.add_metabolites(missing_mets)
    model.add_reactions(missing_exchanges)

    # identify transporters for each biolog component in the universal model
    # and pick one that will enable transport in the gapfilling problem.
    transporters_in_universal = [rxn for rxn in universal.reactions if len(rxn.compartments)>1]
    for met in media_dicts.keys():
        metabolite = model.metabolites.get_by_id(met)
        base_met_id = met.split('_')[0]
        rxns_with_metabolite = metabolite.reactions
        transport = False
        for rxn in rxns_with_metabolite:
            metabolites = [met_in_rxn.id for met_in_rxn in rxn.metabolites]
            if (base_met_id+'_e' in metabolites and base_met_id+'_c' in metabolites):
                transport = True

        pick_transporter = {}
        if not transport:
            print("missing transporter for " + metabolite.name)
            for rxn in transporters_in_universal:
                metabolites = [met_in_rxn.id for met_in_rxn in rxn.metabolites]
                if (base_met_id+'_e' in metabolites and base_met_id+'_c' in metabolites):
                    pick_transporter[met] = rxn.id

    # Add the transporters to the model
    transporters_to_add = list(pick_transporter.values())
    transporter_list = []
    for rxn in transporters_to_add:
        transporter_list.append(universal.reactions.get_by_id(rxn).copy())
    model.add_reactions(transporter_list)

    # remove the added transporters from the universal model
    universal.remove_reactions([universal.reactions.get_by_id(rxn) for rxn in transporters_to_add])

    # select a subset of the biolog conditions to perform gapfilling with
    sources = list(media_dicts.keys())
    sub_dict = {sources[0]:media_dicts[sources[0]],
               sources[1]:media_dicts[sources[1]],
               sources[2]:media_dicts[sources[2]],
               sources[3]:media_dicts[sources[3]],
               sources[4]:media_dicts[sources[4]]}

    num_cycles = 5
    lower_bound = 0.05
    flux_cutoff = 1E-10
    ensemble = expand.iterative_gapfill_from_binary_phenotypes(model,universal,sub_dict,num_cycles,\
                                         lower_bound=lower_bound,\
                                         inclusion_threshold=1E-10,\
                                         exchange_reactions=False,\
                                         demand_reactions=False,\
                                         exchange_prefix='EX')

    # the number of models in the ensemble should equal the number of cycles
    # unless a duplicate solution was found; for this test case, this seems
    # to happen ~25% of the time, so we'll loosen the restriction.
    assert len(ensemble.members) > num_cycles/2
    # each member of the ensemble should be able to produce biomass in each
    # biolog condition
    ex_rxns = [rxn for rxn in ensemble.base_model.reactions \
                        if rxn.id.startswith('EX_')]
    for source in sub_dict.keys():
        # close all exchange reactions
        for rxn in ex_rxns:
            rxn.lower_bound = 0
        ensemble.base_model.medium = sub_dict[source]
        for member in ensemble.members:
            ensemble.set_state(member)
            # member should produce the minimum amount of required biomass
            # flux or more
            assert ensemble.base_model.slim_optimize() > lower_bound*0.99


def load_universal_modelseed():
    seed_rxn_table = pd.read_csv('./medusa/test/data/reactions_seed_20180809.tsv',sep='\t')
    seed_rxn_table['id'] = seed_rxn_table['id'] + '_c'
    universal = load_json_model('./medusa/test/data/universal_mundy.json')
    # remove any reactions from the universal that don't have "OK" status
    # in modelSEED (guards against mass and charge-imbalanced reactions)
    ok_ids = list(seed_rxn_table.loc[(seed_rxn_table['status'] == 'OK') | (seed_rxn_table['status'] == 'HB')]['id'])
    remove_rxns = []
    for reaction in universal.reactions:
        if reaction.id not in ok_ids:
            remove_rxns.append(reaction)
    universal.remove_reactions(remove_rxns)
    # remove metabolites from the universal that are no longer present in any
    # reactions.
    mets_in_reactions = []
    for reaction in universal.reactions:
        mets = [met.id for met in reaction.metabolites]
        mets_in_reactions.extend(mets)
    mets_in_reactions = set(mets_in_reactions)

    mets_missing_reactions = []
    for metabolite in universal.metabolites:
        if metabolite.id not in mets_in_reactions:
            mets_missing_reactions.append(metabolite)
    universal.remove_metabolites(mets_missing_reactions)

    universal.repair()
    return universal

def load_modelseed_model(model_name):
    if model_name == 'Staphylococcus aureus':
        model = load_json_model('./medusa/test/data/'+model_name+'.json')
    else:
        raise ValueError('Unsupported model_name provided')
    return model
