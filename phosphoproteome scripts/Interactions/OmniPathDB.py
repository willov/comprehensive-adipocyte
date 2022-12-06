
#%%
# Install omnipath and pandas via pip
import omnipath as op
import pandas as pd

#%%
#print(op.interactions.AllInteractions.resources())
interactions = op.interactions.PostTranslational.get(organism="mouse") #Accessed originally 211215
interactions=interactions.rename(columns={"source":"Source", "target":"Target", "n_primary_sources":"Conf"})
interactions.to_csv('OmniPathInteractions.csv')

# %% Get a list of sources used in the interaction list
sources=interactions["sources"]
sources=sources.str.split(';', expand=True).values.flatten()
sources=[s for s in sources if s is not None]

# %%
