from __future__ import absolute_import

import cobra
import json

from medusa.core.ensemble import Ensemble


from os.path import abspath, dirname, join

medusa_directory = abspath(join(dirname(abspath(__file__)), ".."))
data_dir = join(medusa_directory,"test","data","")

def create_test_model(model_name="textbook"):
    """Returns an ensemble of models for testing
    model_name: str
        One of 'ASF356' or 'ASF519'
    """
    if model_name == "ASF356":
        base_model = cobra.io.load_json_model(join(data_dir, 'ASF356_base_model.json'))
        with open(join(data_dir, 'ASF356_base_model_reaction_diffs.json'),'r') as infile:
            reaction_diffs = json.load(infile)

    elif model_name == "ASF519":
        base_model = cobra.io.load_json_model(join(data_dir, 'ASF519_base_model.json'))
        with open(join(data_dir, 'ASF519_base_model_reaction_diffs.json'),'r') as infile:
            reaction_diffs = json.load(infile)

    else:
        base_model = cobra.create_test_model(model_name)



    test_model = medusa.Ensemble(base_id=base_model.id)
    test_model.base_model = base_model
    test_model.reaction_diffs = reaction_diffs
    return test_model
