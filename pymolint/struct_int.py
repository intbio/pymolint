"""
Get contacts and other interactions between selection in an MDanalysis parsable file (or stream)
"""
from .base import *
import MDAnalysis as mda
import pandas as pd
import numpy as np
__author__="Alexey Shaytan"



class struct2cont:
    """"
    Analyze contacts in 3D structure
    
    """
    
    def __init__(self,struct,selA,selB=None, format=None, d_threshold=3.9,exclude_bonded=False,half_matrix=True):
        #Try to open
        if(isinstance(struct,mda.Universe)):
            self.u=struct
        else:
            self.u=mda.Universe(struct,format=format)
        self.selA_mda=self.u.select_atoms(selA)
        if selB:
            self.selB_mda=self.u.select_atoms(selB)
        else:
            self.selB_mda=self.u.select_atoms('not (bynum %s)'%' '.join([str(i) for i in self.selA_mda.ids]))
        #Let's calculate contacts
        self.contacts=find_contacts(self.selA_mda.positions,self.selA_mda.ids,self.selB_mda.positions,self.selB_mda.ids,\
                                   d_threshold=d_threshold,exclude_bonded=exclude_bonded,half_matrix=half_matrix)
                                                                    
                                                                    
    def get_list(self):
        return self.contacts

    def get_df(self):
        id2name={}
        id2segid={}
        id2resid={}                                                            
              
        id2resid={a.id:a.resid for a in self.selA_mda.atoms}
        id2resid.update({a.id:a.resid for a in self.selB_mda.atoms})
                                                                    
        id2segid={a.id:a.segid for a in self.selA_mda.atoms}
        id2segid.update({a.id:a.segid for a in self.selB_mda.atoms})
                                                                    
        id2name={a.id:a.name for a in self.selA_mda.atoms}
        id2name.update({a.id:a.name for a in self.selB_mda.atoms})
                                                                    
        self.contdf=pd.DataFrame({'A_atom_id':[i[0] for i in self.contacts],'B_atom_id':[i[1] for i in self.contacts],'A_atom_name':[id2name[i[0]] for i in self.contacts],'A_atom_resid':[id2resid[i[0]] for i in self.contacts],'A_atom_segid':[id2segid[i[0]] for i in self.contacts],'B_atom_name':[id2name[i[1]] for i in self.contacts],'B_atom_resid':[id2resid[i[1]] for i in self.contacts],'B_atom_segid':[id2segid[i[1]] for i in self.contacts]})
        return self.contdf
    
    
    def get_num_int_profile(self):
        """
        Gets number of contacts per residue in selection A
        Returns dataframe with three columns 'segid' 'resid' 'num_int' ordered by resid and segid
        """
        self.get_df()
        self.num_int_profile=self.contdf.groupby(['A_atom_segid','A_atom_resid']).size().reset_index()
        self.num_int_profile.columns=['segid','resid','num_int']
        return self.num_int_profile
        
            
            
            
            
            
            
            
                                      