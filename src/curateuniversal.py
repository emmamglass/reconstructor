##setting up the universal reactions database
import os
import pickle
from cobra.manipulation.delete import *

script_path = str(os.path.dirname(os.path.realpath(__file__)))

filename = script_path + '/refs/universal.pickle'
with open(filename, 'rb') as f: universal = pickle.load(f)

universal.metabolites.cpd17677_e.formula = 'C17H24O2'
universal.metabolites.cpd17677_c.formula = 'C17H24O2'
universal.metabolites.cpd20970_c.formula = 'C25H39N5O21P3R2'
universal.metabolites.cpd26289_c.formula = 'Cl3Fe'
universal.metabolites.cpd12759_c.formula = 'C5H7NO3R2'
universal.metabolites.cpd21495_c.formula = 'C21H24O3'
universal.metabolites.cpd20345_c.formula = 'C38H49N11O30P4R2'
universal.metabolites.cpd12434_c.formula = 'C25H45NO8R2'
universal.metabolites.cpd25541_c.formula = 'C45H52CoN4O14'
universal.metabolites.cpd12163_c.formula = 'C18H32O16'
universal.metabolites.cpd20971_c.formula = 'C45H57CoN5O13'
universal.metabolites.cpd25094_c.formula = 'C7H4Cl2CuNO2S2'
universal.metabolites.cpd12148_c.formula = 'C6H10O5R2'
universal.metabolites.cpd11851_c.formula = 'C6H10O5'
universal.metabolites.cpd24958_c.formula = 'C42H38FeN5O15'
universal.metabolites.cpd30665_c.formula = 'C25H31ClN3'
universal.metabolites.cpd21089_c.formula = 'C20H24O2'
universal.metabolites.cpd30678_e.formula = 'C8H9Cl3N2S'
universal.metabolites.cpd20944_c.formula = 'C43H69N9O23'
universal.metabolites.cpd17683_c.formula = 'C17H22O3'
universal.metabolites.cpd24567_c.formula = 'PbS'
universal.metabolites.cpd20338_c.formula = 'C19H27N3O17P2R2'
universal.metabolites.cpd11746_e.formula = 'C6H10O5R2'
universal.metabolites.cpd26289_e.formula = 'Cl3Fe'
universal.metabolites.cpd11746_c.formula = 'C6H10O5R2'
universal.metabolites.cpd12039_c.formula = 'HO10P3'
universal.metabolites.cpd21095_c.formula = 'C20H34O7P2'
universal.metabolites.cpd17684_c.formula = 'C18H24O3'
universal.metabolites.cpd20958_c.formula = 'C14H23N6O8'
universal.metabolites.cpd19935_c.formula = 'CH3Hg'
universal.metabolites.cpd20957_c.formula = 'C10H19N5O5'
universal.metabolites.cpd24797_c.formula = 'C42H47N6NiO14'
universal.metabolites.cpd12511_c.formula = 'C24H41O21R'
universal.metabolites.cpd25097_c.formula = 'C7H4ClCuNO3S'
universal.metabolites.cpd20961_c.formula = 'C20H36N10O9'
universal.metabolites.cpd11720_c.formula = 'C5H10O7PR2'
universal.metabolites.cpd21084_c.formula = 'C17H29N3O10R3'
universal.metabolites.cpd17685_c.formula = 'C18H23O4'
universal.metabolites.cpd12077_c.formula = 'C5H8R2'
universal.metabolites.cpd20340_c.formula = 'C28H38N6O24P3R2'
universal.metabolites.cpd25050_c.formula = 'C10H10MoN5O7PS3'
universal.metabolites.cpd21100_c.formula = 'C17H24O3'
universal.metabolites.cpd19503_c.formula = 'C10H10MoN5O8PS2'
universal.metabolites.cpd21353_c.formula = 'C28H41N6O7R5S'
universal.metabolites.cpd25539_c.formula = 'C42H37CoN4O16'
universal.metabolites.cpd24798_c.formula = 'C42H47N6NiO14'
universal.metabolites.cpd17682_c.formula = 'C18H23O4'
universal.metabolites.cpd21288_c.formula = 'C19H22MoN8O15P2S2'
universal.metabolites.cpd17679_c.formula = 'C17H21O4'
universal.metabolites.cpd21101_c.formula = 'C18H26O3'
universal.metabolites.cpd12120_c.formula = 'C6H10O5R2'
universal.metabolites.cpd20943_c.formula = 'C39H65N8O20'
universal.metabolites.cpd21352_c.formula = 'C13H17N6O7R5S'
universal.metabolites.cpd26290_e.formula = 'Cl2Fe'
universal.metabolites.cpd25095_c.formula = 'C8H4Cl5CuNO2S2'
universal.metabolites.cpd30665_e.formula = 'C25H31ClN3'
universal.metabolites.cpd24099_e.formula = 'C71H97CoN18O21P2'
universal.metabolites.cpd11658_c.formula = 'C12H20O10R2'
universal.metabolites.cpd22240_c.formula = 'C72H99CoN18O20P2'
universal.metabolites.cpd20969_c.formula = 'C19H26N3O20P3R2'
universal.metabolites.cpd25093_c.formula = 'C7H4ClCuNO2S2'
universal.metabolites.cpd17682_e.formula = 'C18H23O4'
universal.metabolites.cpd30678_c.formula = 'C8H9Cl3N2S'
universal.metabolites.cpd12234_c.formula = 'C6H11O6R'

#remove mass imbalanced reactions
balanced_rxns = []
charge_only_rxns = []
mass_imbalanced = []

count=0
for rxn in universal.reactions:
	try:
		if not rxn.check_mass_balance(): #Check is the dictionary is empty
			balanced_rxns.append(rxn.id)
		elif (list(rxn.check_mass_balance().keys())[0] == 'charge') & (len(list(rxn.check_mass_balance().keys())) == 1):
			charge_only_rxns.append(rxn.id)
		else:
			mass_imbalanced.append(rxn.id)
	except:
		mass_imbalanced.append(rxn.id)
	count+=1
	print(count)
        
print(len(balanced_rxns))
print(len(charge_only_rxns))
print(len(mass_imbalanced))

count=0
mass_imbalanced_bounds = {}
for rxn_id in mass_imbalanced:
	if rxn_id.startswith("rxn"):
		rxn = universal.reactions.get_by_id(rxn_id)
		universal.remove_reactions(rxn)
		print(rxn)
		count+=1

pickle.dump(universal, open(script_path + '/refs/universal.pickle', 'wb'))


