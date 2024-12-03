import pandas as pd
import numpy as np
import argparse as ap
from collections import defaultdict



def splicing_altering(df):
    '''
    Basically, we just refer to SpliceVault top-4 events and SpliceAI delta scores.

    The possible splicing altering effects can be categorized into:
    1. Frameshift effect, can be induced by exon skipping, cryptic donor/acceptor activation.
    2. In-frame retention, can be induced by exon skipping, cryptic donor/acceptor activation.
    3. In-frame truncation, can be induced by exon skipping, cryptic donor/acceptor activation.
    '''



def SpliceVault_interpretation(value):
    '''
    The input value is the top-4 SpliceVault events:
    Top1:ES;3-4;12%;Frameshift|Top2:ES;3;0.05%;Frameshift|Top3:CD;+228;0.03%;Frameshift|Top4:CD;+252;0.02%;Frameshift
    '''
    events = value.split('|')
    event_dict = defaultdict(list)
    for event in events:
        event_info = event.split(':')[1]
        event_dict['type'] = event_info.split(';')[0]
        event_dict['position'] = event_info.split(';')[1]
        event_dict['frequency'] = event_info.split(';')[2]
        event_dict['frame_shift'] = event_info.split(';')[3]

    frameshift_frac = len([fr for fr in event_dict['frame_shift'] if fr == "Frameshift"])/4
    if frameshift_frac == 1:
        return "Frameshift"
    
    event_df = pd.DataFrame(event_dict)
    frameshift_freq = event_df.loc[event_df["frame_shift"] == "Frameshift", "frequency"].sum()
    inframe_freq = event_df.loc[event_df["frame_shift"] != "Frameshift", "frequency"].sum()
    

    
    
