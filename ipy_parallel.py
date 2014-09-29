# parallel setup of environment

profile="default" 

from IPython.parallel import *

client = Client(profile=profile)

balanced       = client.load_balanced_view()
direct         = client[:]

balanced.block = False
direct.block   = True

with direct.sync_imports():
    from collections import defaultdict
    import json
    import sys
    
    import numpy as np
    import scipy as sp
    import matplotlib.pyplot as plt
    import pandas as pd
    
    from matplotlib import rcParams
    import matplotlib.cm as cm
    import matplotlib as mpl
    import Bio as bp
    import HTSeq as ht

from IPython.parallel import CompositeError
CompositeError.tb_limit = 1
