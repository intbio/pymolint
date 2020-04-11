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
        self.d_threshold=d_threshold
        self.exclude_bonded=exclude_bonded
        self.half_matrix=half_matrix
        
        
        self.id2name={}
        self.id2segid={}
        self.id2resid={}                                                            
              
        self.id2resid={a.id:a.resid for a in self.selA_mda.atoms}
        self.id2resid.update({a.id:a.resid for a in self.selB_mda.atoms})
        
        self.id2resname={a.id:a.resname for a in self.selA_mda.atoms}
        self.id2resname.update({a.id:a.resname for a in self.selB_mda.atoms})
                                                                    
        self.id2segid={a.id:a.segid for a in self.selA_mda.atoms}
        self.id2segid.update({a.id:a.segid for a in self.selB_mda.atoms})
                                                                    
        self.id2name={a.id:a.name for a in self.selA_mda.atoms}
        self.id2name.update({a.id:a.name for a in self.selB_mda.atoms})
        
        
        #Let's calculate contacts
        self.contacts=find_contacts(self.selA_mda.positions,self.selA_mda.ids,self.selB_mda.positions,self.selB_mda.ids,\
                                   d_threshold=d_threshold,exclude_bonded=exclude_bonded,half_matrix=half_matrix)
        self.num=len(self.contacts)
                                                    
                                                                    
    def get_list(self):
        return self.contacts

    def get_df(self):
#         print('get_df profile from struct2cont called')
                                                            
        self.contdf=pd.DataFrame({'A_atom_id':[i[0] for i in self.contacts],'A_atom_name':[self.id2name[i[0]] for i in self.contacts],'A_resname':[self.id2resname[i[0]] for i in self.contacts],'A_segid':[self.id2segid[i[0]] for i in self.contacts],'A_resid':[self.id2resid[i[0]] for i in self.contacts],'B_atom_id':[i[1] for i in self.contacts],'B_atom_name':[self.id2name[i[1]] for i in self.contacts],'B_resname':[self.id2resname[i[1]] for i in self.contacts],'B_segid':[self.id2segid[i[1]] for i in self.contacts],'B_resid':[self.id2resid[i[1]] for i in self.contacts],'dist,A':[i[2] for i in self.contacts]})
        return self.contdf
    
    
    
    def get_num_int_profile(self):
        """
        Gets number of contacts per residue in selection A
        Returns dataframe with three columns 'segid' 'resid' 'num_int' ordered by resid and segid
        """
        self.get_df()
        self.num_int_profile=self.contdf.groupby(['A_segid','A_resid']).size().reset_index()
        self.num_int_profile.columns=['segid','resid','num_int']
        return self.num_int_profile
        
            
class trj2cont(struct2cont):
    """
    A class to analyze contacts in a trajectory, should get an MD analysis universe as input
    """
    def __init__(self,u,selA,selB=None, format=None, d_threshold=3.9,exclude_bonded=False,half_matrix=True,time=(None,None)):
        super().__init__(u,selA,selB=selB, d_threshold=d_threshold,exclude_bonded=exclude_bonded,half_matrix=half_matrix)
        self.time=list(time)
        self.contacts_series=[]
        if time[0] is None:
            self.time[0]=0
        if time[1] is None:
            self.time[1]=u.trajectory.n_frames
        
        for ts in u.trajectory[self.time[0]:self.time[1]]:
            self.recalc()
            self.contacts_series.append(self.contacts)
    
    def get_list_series(self):
        return self.contacts_series
    
    def recalc(self):
        self.contacts=find_contacts(self.selA_mda.positions,self.selA_mda.ids,self.selB_mda.positions,self.selB_mda.ids,\
                                   d_threshold=self.d_threshold,exclude_bonded=self.exclude_bonded,half_matrix=self.half_matrix)
        self.num=len(self.contacts)
        
    def get_df(self):
#         print('get_df profile from trj2cont called')
#        # #This way is not efficient but straight forward 
#         df=pd.DataFrame()
#         for t in range(self.time[0],self.time[1]):
#             self.contacts=self.contacts_series[t-self.time[0]]
#             super().get_df()
#             self.contdf['Time']=t
#             df=pd.concat([df, self.contdf])
#         self.contdf=df

#Trying a more efficient way
        self.contdf=pd.DataFrame({'A_atom_id':[i[0] for s in self.contacts_series for i in s],'A_atom_name':[self.id2name[i[0]]  for s in self.contacts_series for i in s],'A_resname':[self.id2resname[i[0]]  for s in self.contacts_series for i in s],'A_segid':[self.id2segid[i[0]]  for s in self.contacts_series for i in s],'A_resid':[self.id2resid[i[0]]  for s in self.contacts_series for i in s],'B_atom_id':[i[1]  for s in self.contacts_series for i in s],'B_atom_name':[self.id2name[i[1]]  for s in self.contacts_series for i in s],'B_resname':[self.id2resname[i[1]]  for s in self.contacts_series for i in s],'B_segid':[self.id2segid[i[1]]  for s in self.contacts_series for i in s],'B_resid':[self.id2resid[i[1]]  for s in self.contacts_series for i in s],'dist':[i[2]  for s in self.contacts_series for i in s],'Time':[frame_num+self.time[0]  for s,frame_num in zip(self.contacts_series,range(len(self.contacts_series))) for i in range(len(s))]})
        return self.contdf

    def get_num_int_profile(self):
        super().get_num_int_profile()
        self.avr_num_int_profile=self.num_int_profile
        self.avr_num_int_profile['num_int']=self.num_int_profile['num_int']/float(self.time[1]-self.time[0])
        return self.avr_num_int_profile


            
            
            
            
                                      