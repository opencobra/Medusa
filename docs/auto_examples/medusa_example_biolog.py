"""
Generating an ensemble by gapfilling using growth phenotype data
=========================

In this example, we demonstrate how to generate an ensemble by gapfilling a model using growth phenotyping data from biolog growth conditions which contain a single Carbon or Nitrogen source.
"""

# Load the test model for Staphylococcus aureus, originally generated with ModelSEED
import medusa
from medusa.test import create_test_model
model = create_test_model('Saureus_seed')

# Load the biolog data from Plata et al., Nature 2014
from medusa.test import load_biolog_plata
biolog_base_composition, biolog_base_dict, biolog_thresholded = load_biolog_plata()
biolog_base_composition

# extract the data for just Staphylococcus aureus
test_mod_pheno = biolog_thresholded.loc['Staphylococcus aureus']
test_mod_pheno = list(test_mod_pheno[test_mod_pheno == True].index)

# load the universal reaction database
from medusa.test import load_universal_modelseed
from cobra.core import Reaction
universal = load_universal_modelseed()

# check for biolog base components in the model and record
# the metabolites/exchanges that need to be added
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

# now perform gapfilling. This may take 20-30 seconds per cycle.
from medusa.reconstruct.expand import iterative_gapfill_from_binary_phenotypes

# select a subset of the biolog conditions to perform gapfilling with
sources = list(media_dicts.keys())
sub_dict = {sources[0]:media_dicts[sources[0]],
           sources[1]:media_dicts[sources[1]],
           sources[2]:media_dicts[sources[2]],
           sources[3]:media_dicts[sources[3]],
           sources[4]:media_dicts[sources[4]]}

num_cycles = 100
lower_bound = 0.05
flux_cutoff = 1E-10
ensemble = iterative_gapfill_from_binary_phenotypes(model,universal,sub_dict,num_cycles,\
                                     lower_bound=lower_bound,\
                                     inclusion_threshold=1E-10,\
                                     exchange_reactions=False,\
                                     demand_reactions=False,\
                                     exchange_prefix='EX')

print(len(ensemble.members))

# TODO: with the ensemble in hand, generate a histogram of predicted
# flux through biomass to display for the thumbnail.
%matplotlib inline
from medusa.flux_analysis import flux_balance
import matplotlib.pylab as plt
predicted_growth = flux_balance.optimize_ensemble(ensemble, return_flux=['bio1'])
x = predicted_growth['bio1']
fig, ax = plt.subplots()
plt.hist(x=x)
ax.set_xlabel('Biomass flux',size=16)
ax.set_ylabel('Density',size=16)
plt.show()

###############################################################################
# Here is an embedded section of RST text
# ---------------------------------------
#
# We'll use these later to break up the tutorial steps
