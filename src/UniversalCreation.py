#!/usr/bin/env python
# coding: utf-8

# In[1]:


import cobra
from copy import deepcopy
import json
import pickle

# building the universal model, add all reactions in together make sure syntax is correct, fist step of model creation

# people can build any universal model that they want, the main story is the optimization method

# you can make any universal creation database that you want. can definitely be improved, can remove all photosynthetic
# reactions that are present in kegg if you want, but my version does not have photosynthetic reactions removed

# minimizing flux is better, since you are minimizing number of enzymes needed to perform a task which is just biologically true?


# In[2]:


universal = cobra.Model('Prokaryote_Universal')
universal.name = 'Prokaryote_Universal'


# In[3]:


# Missing metabolite dictionary
with open('/home/mjenior/Desktop/repos/model_building/refs/compounds.json') as f: compounds = json.load(f)
newCpdDict = {} 

for cpd in compounds:

    cpd_name = cpd['name'].title()
    cpd_formula = cpd['formula']
    
    cpd_id = cpd['id'] + '_c'
    newCpdDict[cpd_id] = cobra.Metabolite(cpd_id, formula=cpd_formula, name=cpd_name, compartment='cytosol')
      
    cpd_id = cpd['id'] + '_e'
    newCpdDict[cpd_id] = cobra.Metabolite(cpd_id, formula=cpd_formula, name=cpd_name, compartment='extracellular')

del compounds


# In[4]:


# Find new reactions
with open('/home/mjenior/Desktop/repos/model_building/refs/reactions.json') as inFile: reactions = json.load(inFile)
currency = set(['cpd00004_c','cpd00003_c','cpd00005_c','cpd00006_c','cpd00002_c','cpd00008_c','cpd00038_c','cpd00031_c','cpd00061_c'])

newReactions = []
for rxn in reactions:
    if['status'] in ['EMPTY', 'CPDFORMERROR']:
        continue
    
    rxn_id = rxn['id'] + '_c'
    name = rxn['name'].title()
    raw_stoich = rxn['stoichiometry'].split(';')
    if raw_stoich == ['']: continue
    raw_stoich = [x.split(':')[0:2] for x in raw_stoich]
    raw_stoich = [[float(x[0]),x[1]+'_c'] for x in raw_stoich]
    
    if rxn['reversibility'] == '>':
            bounds = (0.0, 1000.0)
    elif rxn['reversibility'] == '<':
        bounds = (-1000.0, 0.0)
    else:
        bounds = (-1000.0, 1000.0)
        
    # Cytoplasmic reactions
    if rxn['is_transport'] == 0:
        stoich = {}
        for x in raw_stoich:
            try:
                stoich[newCpdDict[x[1]]] = x[0]
            except:
                stoich[universal.metabolites.get_by_id(x[1])] = x[0]
        
    # Transport reactions
    elif rxn['is_transport'] == 1:
        stoich = {}
        for x in raw_stoich:
            if x[0] < 0.0:
                if x[1] in currency:
                    new_id = x[1]
                else:
                    new_id = x[1].replace('_c', '_e')
                try:
                    stoich[newCpdDict[new_id]] = x[0]
                except:
                    stoich[universal.metabolites.get_by_id(new_id)] = x[0]
            else:
                try:
                    stoich[newCpdDict[x[1]]] = x[0]
                except:
                    stoich[universal.metabolites.get_by_id(x[1])] = x[0]

    # Create reaction object
    new_rxn = cobra.Reaction(rxn_id)
    new_rxn.name = name
    new_rxn.bounds = bounds
    new_rxn.add_metabolites(stoich)
    newReactions.append(new_rxn) 

del newCpdDict


# In[5]:


universal.add_reactions(newReactions)
del newReactions


# In[6]:


universal


# In[7]:


# Add Tom's universal reactions
to_be_rm = ['EX_cpd21851_e','EX_cpd26978_e','EX_cpd35774_e','EX_cpd28247_e','EX_cpd33670_e','EX_cpd32312_e','EX_cpd12556_e','EX_cpd31712_e','EX_cpd37289_e','EX_cpd15276_e','EX_cpd15521_e','EX_cpd15522_e','EX_cpd15523_e','EX_cpd15524_e','EX_cpd15525_e','EX_cpd15526_e','EX_cpd15527_e','EX_cpd15528_e','EX_cpd15529_e','EX_cpd15530_e','EX_cpd15531_e','EX_cpd15532_e','EX_cpd15533_e','EX_cpd15534_e','EX_cpd15535_e','EX_cpd15536_e','EX_cpd15537_e','EX_cpd15538_e','EX_cpd15539_e','EX_cpd15540_e','EX_cpd15541_e','EX_cpd15542_e','EX_cpd15543_e','EX_cpd15544_e','EX_cpd15545_e','EX_cpd15546_e','EX_cpd15547_e','EX_cpd15548_e','EX_cpd15271_e','EX_cpd11825_e','EX_cpd15268_e','EX_cpd15239_e','EX_cpd15277_e','EX_cpd15294_e','EX_cpd11466_e','EX_cpd11468_e','EX_cpd15457_e','EX_cpd15447_e','EX_cpd15582_e','EX_cpd15459_e','EX_cpd15493_e','EX_cpd15432_e','EX_cpd15350_e','EX_cpd15363_e','EX_cpd15362_e','EX_cpd15355_e','EX_cpd15354_e','EX_cpd15358_e','EX_cpd15357_e','EX_cpd15336_e','EX_cpd15337_e','EX_cpd15338_e','EX_cpd15339_e','EX_cpd15340_e','EX_cpd15341_e','EX_cpd15342_e','EX_cpd15343_e','EX_cpd15344_e','EX_cpd15345_e','EX_cpd15346_e','EX_cpd15347_e','EX_cpd15348_e','EX_cpd15349_e','EX_cpd15306_e','EX_cpd15307_e','EX_cpd15308_e','EX_cpd15309_e','EX_cpd15310_e','EX_cpd15311_e','EX_cpd15312_e','EX_cpd15576_e','EX_cpd15448_e','EX_cpd15451_e','EX_cpd15452_e','EX_cpd15208_e','EX_cpd15202_e','EX_cpd15288_e','EX_cpd15198_e','EX_cpd15940_e','EX_cpd03092_e','EX_cpd09429_e','EX_cpd15220_e','EX_cpd15221_e','EX_cpd15224_e','EX_cpd15225_e','EX_cpd15241_e','EX_cpd37153_e','EX_cpd25615_e','EX_cpd27183_e']


# In[8]:


# Correct compartment labels
for cpd in old_universal.metabolites:
    if cpd.compartment == 'c':
        cpd.compartment = 'cytosol'
    elif cpd.compartment == 'e':
        cpd.compartment = 'extracellular'


# In[9]:


# Find missing reactions
newReactions = []
for rxn in old_universal.reactions:
    if rxn.id == 'bio1':
        continue
        
    try:
        test = universal.reactions.get_by_id(rxn.id)
    except:
        newReactions.append(deepcopy(rxn)) 

del old_universal


# In[10]:


universal.add_reactions(newReactions)
del newReactions


# In[11]:


universal


# In[12]:


# Identify and correct missing exchange reactions
for cpd in universal.metabolites:
    if cpd.compartment != 'extracellular': continue
    try:
        test = universal.reactions.get_by_id('EX_' + cpd.id)
    except KeyError:
        universal.add_boundary(cpd, type='exchange', reaction_id='EX_' + cpd.id, lb=-1000.0, ub=1000.0)


# In[13]:


universal


# In[14]:


# Purge or correct bad formula metabolites
matches = ['(', ')', '.', '*']

remove_rxns = []
for cpd in universal.metabolites:
    try:
        if any(x in cpd.formula for x in matches):
            test = 'universal.metabolites.' + cpd.id + '.formula = \'' + cpd.formula + '\''
            print(test)
    except TypeError:
        remove_rxns += list(cpd.reactions)

remove_rxns = list(set(remove_rxns))


# In[15]:


universal.remove_reactions(remove_rxns)
del remove_rxns


# In[3]:


# Correct bad formulas
universal.metabolites.cpd17677_e.formula = 'C17H24O2'
universal.metabolites.cpd17677_c.formula = 'C17H24O2'
universal.metabolites.cpd20970_c.formula = 'C25H186N5O21P3R2'
universal.metabolites.cpd26289_c.formula = 'FeCl3'
universal.metabolites.cpd12759_c.formula = 'C5H7NO3'
universal.metabolites.cpd21495_c.formula = 'C21H24O3'
universal.metabolites.cpd20345_c.formula = 'C38H53N11O30P4R2'
universal.metabolites.cpd12434_c.formula = 'C25H45NO8R'
universal.metabolites.cpd25541_c.formula = 'CoC45H53N4O14'
universal.metabolites.cpd12163_c.formula = 'C12H20O11'
universal.metabolites.cpd20971_c.formula = 'CoC45H57N5O13'
universal.metabolites.cpd25094_c.formula = 'CuClC7H4NO2S2'
universal.metabolites.cpd12148_c.formula = 'C6H10O5'
universal.metabolites.cpd11851_c.formula = 'C6H10O5'
universal.metabolites.cpd24958_c.formula = 'FeC42H41N5O15'
universal.metabolites.cpd30665_c.formula = 'ClC25H40N3'
universal.metabolites.cpd21089_c.formula = 'C20H24O2'
universal.metabolites.cpd30678_e.formula = 'C8H9Cl3N2S'
universal.metabolites.cpd20944_c.formula = 'C43H71N9O23'
universal.metabolites.cpd17683_c.formula = 'C17H22O3'
universal.metabolites.cpd24567_c.formula = 'PbS'
universal.metabolites.cpd20338_c.formula = 'C19H29N3O17P2R2'
universal.metabolites.cpd11746_e.formula = 'C6H10O5'
universal.metabolites.cpd26289_e.formula = 'FeCl3'
universal.metabolites.cpd11746_c.formula = 'C6H10O5'
universal.metabolites.cpd12039_c.formula = 'H2O6P2'
universal.metabolites.cpd21095_c.formula = 'C20H36O7P2'
universal.metabolites.cpd17684_c.formula = 'C18H24O3'
universal.metabolites.cpd20958_c.formula = 'C14H24N6O8'
universal.metabolites.cpd19935_c.formula = 'CH3Hg'
universal.metabolites.cpd20957_c.formula = 'C10H19N5O5'
universal.metabolites.cpd24797_c.formula = 'NiC42H47N6O14'
universal.metabolites.cpd12511_c.formula = 'C24H41O21'
universal.metabolites.cpd25097_c.formula = 'CuClC7H4NO3S'
universal.metabolites.cpd20961_c.formula = 'C20H36N10O9'
universal.metabolites.cpd11720_c.formula = 'C5H11O7P'
universal.metabolites.cpd21084_c.formula = 'C17H28N3O10R'
universal.metabolites.cpd17685_c.formula = 'C18H24O4'
universal.metabolites.cpd12077_c.formula = 'C5H8'
universal.metabolites.cpd20340_c.formula = 'C28H41N6O24P3R2'
universal.metabolites.cpd25050_c.formula = 'MoC10H10N5O7PS3'
universal.metabolites.cpd21100_c.formula = 'C17H24O3'
universal.metabolites.cpd19503_c.formula = 'MoC10H10N5O8PS2'
universal.metabolites.cpd21353_c.formula = 'C28H41N6O7SR5'
universal.metabolites.cpd25539_c.formula = 'CoC42H37N4O16'
universal.metabolites.cpd24798_c.formula = 'NiC42H47N6O14'
universal.metabolites.cpd17682_c.formula = 'C18H24O4'
universal.metabolites.cpd21288_c.formula = 'MoC19H22N8O15P2S2'
universal.metabolites.cpd17679_c.formula = 'C17H22O4'
universal.metabolites.cpd21101_c.formula = 'C18H26O3'
universal.metabolites.cpd12120_c.formula = 'C6H10O5'
universal.metabolites.cpd20943_c.formula = 'C39H66N8O20'
universal.metabolites.cpd21352_c.formula = 'C13H17N6O7SR5'
universal.metabolites.cpd26290_e.formula = 'FeCl2'
universal.metabolites.cpd26290_c.formula = 'FeCl2'
universal.metabolites.cpd25095_c.formula = 'CuC8H4Cl5NO2S2'
universal.metabolites.cpd30665_e.formula = 'ClC25H40N3'
universal.metabolites.cpd24099_e.formula = 'CoC61H86N13O14P'
universal.metabolites.cpd24099_c.formula = 'CoC61H86N13O14P'
universal.metabolites.cpd11658_c.formula = 'C12H20O10'
universal.metabolites.cpd22240_c.formula = 'CoC162H102N18O20P2'
universal.metabolites.cpd20969_c.formula = 'C19H30N3O20P3R2'
universal.metabolites.cpd25093_c.formula = 'CuClC7H4NO2S2'
universal.metabolites.cpd17682_e.formula = 'C18H22O4'
universal.metabolites.cpd30678_c.formula = 'C8H9Cl3N2S'
universal.metabolites.cpd12234_c.formula = 'C6H11O6'


# In[16]:


# Replace biomass functions
biomass = cobra.io.read_sbml_model('/home/mjenior/Desktop/repos/model_building/refs/biomass.sbml')
biomass.remove_reactions([biomass.reactions.rxn10088_c])
biomass_rxns = [x for x in biomass.reactions]
universal.add_reactions(biomass_rxns)

cpd11416_c = universal.metabolites.cpd11416_c
universal.add_boundary(cpd11416_c, type='sink', reaction_id='EX_biomass', lb=0.0, ub=1000.0)


# In[17]:


# Test objectives - Gram positive
universal.objective = 'biomass_GmPos'
universal.slim_optimize()


# In[18]:


# Gram negative
universal.objective = 'biomass_GmNeg'
universal.slim_optimize()


# In[19]:


universal


# In[2]:


iCdR703 = cobra.io.read_sbml_model('/home/mjenior/Desktop/repos/Jenior_Cdifficile_2019/data/reconstructions/iCdR703.sbml')


# In[26]:


transporters = []
exchs = set([x.id for x in iCdR703.exchanges])
for rxn in iCdR703.reactions:
    if rxn.id in exchs: 
        continue
    else:
        compartments = set()
        for x in list(rxn.metabolites): compartments |= set([x.compartment])
        if len(compartments) > 1: transporters.append(rxn.id)


# In[29]:


add_transport = []
total_added = 0
for x in transporters:
    try:
        test = universal.reactions.get_by_id(x)
    except:
        rxn = deepcopy(iCdR703.reactions.get_by_id(x))
        rxn.gene_reaction_rule = ''
        add_transport.append(rxn)
        total_added +=1

universal.add_reactions(add_transport)
print('Transporters added:', total_added)


# In[30]:


add_exchs = []
total_added = 0
for x in iCdR703.exchanges:
    try:
        test = universal.reactions.get_by_id(x.id)
    except:
        add_exchs.append(deepcopy(x))
        total_added +=1
    
universal.add_reactions(add_exchs)
print('Exchanges added:', total_added)


# In[12]:


x = deepcopy(iCdR703.reactions.rxn27496_c)
x.gene_reaction_rule = ''
universal.add_reactions([x])


# In[5]:


universal


# In[2]:


#universal = pickle.load(open('/home/mjenior/Desktop/repos/reconstructor/refs/universal.pickle', 'rb'))


# In[13]:


to_be_rm = ['EX_cpd21851_e','EX_cpd26978_e','EX_cpd35774_e','EX_cpd28247_e','EX_cpd33670_e','EX_cpd32312_e',
            'EX_cpd12556_e','EX_cpd31712_e','EX_cpd37289_e','EX_cpd15276_e','EX_cpd15521_e','EX_cpd15522_e',
            'EX_cpd15523_e','EX_cpd15524_e','EX_cpd15525_e','EX_cpd15526_e','EX_cpd15527_e','EX_cpd15528_e',
            'EX_cpd15529_e','EX_cpd15530_e','EX_cpd15531_e','EX_cpd15532_e','EX_cpd15533_e','EX_cpd15534_e',
            'EX_cpd15535_e','EX_cpd15536_e','EX_cpd15537_e','EX_cpd15538_e','EX_cpd15539_e','EX_cpd15540_e',
            'EX_cpd15541_e','EX_cpd15542_e','EX_cpd15543_e','EX_cpd15544_e','EX_cpd15545_e','EX_cpd15546_e',
            'EX_cpd15547_e','EX_cpd15548_e','EX_cpd15271_e','EX_cpd11825_e','EX_cpd15268_e','EX_cpd15239_e',
            'EX_cpd15277_e','EX_cpd15294_e','EX_cpd11466_e','EX_cpd11468_e','EX_cpd15457_e','EX_cpd15447_e',
            'EX_cpd15582_e','EX_cpd15459_e','EX_cpd15493_e','EX_cpd15432_e','EX_cpd15350_e','EX_cpd15363_e',
            'EX_cpd15362_e','EX_cpd15355_e','EX_cpd15354_e','EX_cpd15358_e','EX_cpd15357_e','EX_cpd15336_e',
            'EX_cpd15337_e','EX_cpd15338_e','EX_cpd15339_e','EX_cpd15340_e','EX_cpd15341_e','EX_cpd15342_e',
            'EX_cpd15343_e','EX_cpd15344_e','EX_cpd15345_e','EX_cpd15346_e','EX_cpd15347_e','EX_cpd15348_e',
            'EX_cpd15349_e','EX_cpd15306_e','EX_cpd15307_e','EX_cpd15308_e','EX_cpd15309_e','EX_cpd15310_e',
            'EX_cpd15311_e','EX_cpd15312_e','EX_cpd15576_e','EX_cpd15448_e','EX_cpd15451_e','EX_cpd15452_e',
            'EX_cpd15208_e','EX_cpd15202_e','EX_cpd15288_e','EX_cpd15198_e','EX_cpd15940_e','EX_cpd03092_e',
            'EX_cpd09429_e','EX_cpd15220_e','EX_cpd15221_e','EX_cpd15224_e','EX_cpd15225_e','EX_cpd15241_e',
            'EX_cpd37153_e','EX_cpd25615_e','EX_cpd27183_e']

remove_rxns = []
for rxn in to_be_rm:
    remove_rxns.append(universal.reactions.get_by_id(rxn))
    remove_rxns += list(list(universal.reactions.get_by_id(rxn).metabolites)[0].reactions)

remove_rxns = list(set(remove_rxns))
universal.remove_reactions(remove_rxns)
print('Total removed:', len(remove_rxns))
del remove_rxns


# In[5]:


to_be_rm = ['EX_cpd03495_e','EX_cpd15511_e','EX_cpd15518_e','EX_cpd15520_e','EX_cpd00286_e']
remove_rxns = []
for rxn in to_be_rm:
    remove_rxns.append(universal.reactions.get_by_id(rxn))
    remove_rxns += list(list(universal.reactions.get_by_id(rxn).metabolites)[0].reactions)

remove_rxns = list(set(remove_rxns))
universal.remove_reactions(remove_rxns)
print('Total removed:', len(remove_rxns))
del remove_rxns


# In[6]:


to_be_rm = ['EX_cpd26014_e','EX_cpd37198_e','EX_cpd37211_e','EX_cpd02229_e','EX_cpd00842_e',
            'EX_cpd02060_e','EX_cpd03113_e','EX_cpd03114_e','EX_cpd03115_e','EX_cpd03116_e',
            'EX_cpd03118_e','EX_cpd03120_e','EX_cpd03122_e','EX_cpd03125_e','EX_cpd03127_e',
            'EX_cpd03129_e','EX_cpd03130_e','EX_cpd15272_e','EX_cpd15274_e','EX_cpd15238_e',
            'EX_cpd35544_e','EX_cpd04233_e']
remove_rxns = []
for rxn in to_be_rm:
    remove_rxns.append(universal.reactions.get_by_id(rxn))
    remove_rxns += list(list(universal.reactions.get_by_id(rxn).metabolites)[0].reactions)

remove_rxns = list(set(remove_rxns))
universal.remove_reactions(remove_rxns)
print('Total removed:', len(remove_rxns))
del remove_rxns


# In[7]:


universal


# In[8]:


pickle.dump(universal, open('/home/mjenior/Desktop/repos/reconstructor/refs/universal.pickle', 'wb'))


# In[ ]:





# In[ ]:




