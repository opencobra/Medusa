from __future__ import absolute_import

import cobra
import json
import pandas as pd

import cobra.test
from cobra.core import Reaction

from medusa.core.ensemble import Ensemble

from cobra.io import load_json_model
from os.path import abspath, dirname, join

from pickle import load

medusa_directory = abspath(join(dirname(abspath(__file__)), ".."))
data_dir = join(medusa_directory,"test","data","")

def create_test_ensemble(ensemble_name="Staphylococcus aureus"):
    """Returns a previously-generated ensemble for testing
    model_name: str
        One of 'Staphylococcus_aureus_ensemble'
    """
    if ensemble_name == "Staphylococcus aureus":
        with open(join(data_dir, "Staphylococcus_aureus_ensemble.pickle"), 'rb') as infile:
            test_ensemble = load(infile)
    else:
        raise ValueError('ensemble_name does not match one of the test ensembles available')

    return test_ensemble

def create_test_model(model_name="textbook"):
    """Returns a cobra.Model for testing
    model_name: str
        One of ['Staphylococcus aureus'] or any models in cobra.test
    """
    if model_name == "Saureus_seed":
        base_model = cobra.io.load_json_model(join(data_dir, 'Staphylococcus aureus.json'))

    else:
        try:
            base_model = cobra.test.create_test_model(model_name)
        except:
            raise ValueError('model_name does not match one of the test models available')

    return base_model

def load_biolog_plata():
    biolog_base_composition = pd.read_csv(
        join(data_dir,'biolog_base_composition.csv'),sep=',')
    biolog_base_dict = dict(zip(biolog_base_composition['ID'],
                                [1000 for i in range(0,len(biolog_base_composition['ID']))]))
    biolog_thresholded = pd.read_csv(
        join(data_dir,'plata_thresholded.csv'),sep='\t',index_col=0)

    return biolog_base_composition, biolog_base_dict, biolog_thresholded

def load_universal_modelseed():
    seed_rxn_table = pd.read_csv(join(data_dir,'reactions_seed_20180809.tsv'),sep='\t')
    seed_rxn_table['id'] = seed_rxn_table['id'] + '_c'
    universal = load_json_model(join(data_dir,'universal_mundy.json'))
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