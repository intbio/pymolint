"""
Get contacts and other interactions between selection in an MDanalysis parsable file (or stream)
"""
from .base import *
import MDAnalysis as mda
import pandas as pd
import numpy as np
from tqdm.auto import tqdm
import warnings


__author__="Alexey Shaytan"




class struct2int:
    """"
    Analyze contacts in 3D structure or trajectory
    
    """
    
    def __init__(self,struct,selA,selB=None, format=None,start=None,stop=None,step=None,time=None,
                 d_threshold=3.9,exclude_bonded=False,half_matrix=True,int_type='cont',D_A_thresh=3.5,D_H_A_ang_thresh=30):
        """
        struct - Mdanalysis universe or structure file name
        selA - Mdanalysis selection string
        selB - Mdanalysis selection string 
        contact_type = 'cont', 'hbonds', 'ionic'
        """
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
        self.int_type=int_type
        #cleaner way 
        self.time=slice(start,stop,step)
        # legacy compatibility
        if not (time is None):
            self.time=slice(time)
    
        
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
        
        
        #Let's calculate contacts for structure or trj
        self.contacts_series=[]
        #added to keep frame times 
        self.contacts_series_frame=[]    
        self.contacts_series_time=[]       
       
        #self.num=len(self.contacts) - unused_var
        
        #fixed BUG which lead to an incorrect time striding!  
        if int_type=='cont':
            for ts in tqdm(self.u.trajectory[self.time],unit='steps'):
                cont=find_contacts(self.selA_mda.positions,self.selA_mda.ids,self.selB_mda.positions,self.selB_mda.ids,\
                                       d_threshold=d_threshold,exclude_bonded=exclude_bonded,half_matrix=half_matrix)
                self.contacts_series.append(cont)
                self.contacts_series_frame.append(self.u.trajectory.frame)
                self.contacts_series_time.append(self.u.trajectory.time)
            
                
        elif int_type=='hbonds':
            a_map=[((name in DEFAULT_ACCEPTORS['CHARMM27']) or (name in DEFAULT_DONORS['CHARMM27']))  for name in self.selA_mda.names]
            A_heavy=self.selA_mda[a_map]
            A_heavy_indices=self.selA_mda[a_map].ids
            # lets use MDA topology knowledge and generalize selection to hydrogens belonging to the same residues
            A_hydr=self.selA_mda[a_map].residues.atoms.select_atoms('type H')
            # other way would be to select by mask, benchmark needed
            #A_hydr_coord= self.u.select_atoms('type H and around 2 group selection',selection=self.selA_mda[a_map])
            
            b_map=[((name in DEFAULT_ACCEPTORS['CHARMM27']) or (name in DEFAULT_DONORS['CHARMM27']))  for name in self.selB_mda.names]
            B_heavy=self.selB_mda[b_map]
            B_heavy_indices=self.selB_mda[b_map].ids
            # lets use MDA topology knowledge and generalize selection to hydrogens belonging to the same residues
            B_hydr=self.selB_mda[b_map].residues.atoms.select_atoms('type H')
            # other way would be to select by mask, benchmark needed
            #B_hydr_coord= self.u.select_atoms('type H and around 2 group selection',selection=self.selB_mda[b_map])
            for ts in tqdm(self.u.trajectory[self.time],unit='steps'):
                cont=find_hbonds(A_heavy.positions, A_heavy_indices, A_hydr.positions,
                                 B_heavy.positions, B_heavy_indices, B_hydr.positions,
                                 D_A_thresh=D_A_thresh,D_H_A_ang_thresh=D_H_A_ang_thresh,same=False)
                self.contacts_series.append(cont)
                self.contacts_series_frame.append(self.u.trajectory.frame)
                self.contacts_series_time.append(self.u.trajectory.time)
        self.contacts=self.contacts_series[0]
                
                                                   
                                                                    
    def get_list(self):
        """
        this will return contacts for first frame or structure
        """
        return self.contacts
    
    def get_list_series(self):
        """
        this will return a series of contacts for trj
        """
        return self.contacts_series
    
    def get_df_series(self):
        """
        this will return a dataframe of all contacts for every frame (column Time will be there)
        """

        self.contdf_series=pd.DataFrame({'A_atom_id':[i[0] for s in self.contacts_series for i in s],
                                         'A_atom_name':[self.id2name[i[0]]  for s in self.contacts_series for i in s],
                                         'A_resname':[self.id2resname[i[0]]  for s in self.contacts_series for i in s],
                                         'A_segid':[self.id2segid[i[0]]  for s in self.contacts_series for i in s],
                                         'A_resid':[self.id2resid[i[0]]  for s in self.contacts_series for i in s],
                                         'B_atom_id':[i[1]  for s in self.contacts_series for i in s],
                                         'B_atom_name':[self.id2name[i[1]]  for s in self.contacts_series for i in s],
                                         'B_resname':[self.id2resname[i[1]]  for s in self.contacts_series for i in s],
                                         'B_segid':[self.id2segid[i[1]]  for s in self.contacts_series for i in s],
                                         'B_resid':[self.id2resid[i[1]]  for s in self.contacts_series for i in s],
                                         'dist':[i[2]  for s in self.contacts_series for i in s],
                                         #changed this one to prevent accidental timeshift when using stride
                                         'Frame':np.repeat(self.contacts_series_frame,[len(contacts) for contacts in self.contacts_series]),
                                         'Time':np.repeat(self.contacts_series_time,[len(contacts) for contacts in self.contacts_series])})
                                         #'Time':[frame_num+self.time[0]  for s,frame_num in zip(self.contacts_series,range(len(self.contacts_series))) for i in range(len(s))]})
        if self.int_type=='hbonds':
            self.contdf_series['angle']=[i[3]  for s in self.contacts_series for i in s]
            self.contdf_series['direction']=[i[4]  for s in self.contacts_series for i in s]
            self.contdf_series._metadata=['int_type']
            self.contdf_series.int_type='hbonds'
            
        self.contdf_series._metadata=['int_type']
        self.contdf_series.int_type='contacts'
        return self.contdf_series

    def get_df(self):
        return self.get_df_series()
    
    
    def get_df_avr(self):
        """
        this will return a dataframe with average number of contacs for every pair of atoms
        """
        self.get_df_series()
        if self.int_type=='cont':
            self.df_avr=self.contdf_series.groupby(['A_atom_id','A_atom_name','A_resname',
                                                    'A_segid','A_resid','B_atom_id',
                                                    'B_atom_name','B_resname','B_segid','B_resid']).size().reset_index()
            self.df_avr.columns=['A_atom_id','A_atom_name','A_resname',
                                 'A_segid','A_resid','B_atom_id',
                                 'B_atom_name','B_resname','B_segid'
                                 ,'B_resid','num_int']
        elif self.int_type=='hbonds':
            self.df_avr=self.contdf_series.groupby(['A_atom_id','A_atom_name','A_resname',
                                                    'A_segid','A_resid','B_atom_id',
                                                    'B_atom_name','B_resname','B_segid','B_resid','direction']).size().reset_index()
            
            self.df_avr.columns=['A_atom_id','A_atom_name','A_resname',
                                 'A_segid','A_resid','B_atom_id',
                                 'B_atom_name','B_resname','B_segid'
                                 ,'B_resid','direction','num_int']
        self.df_avr['num_int']=self.df_avr['num_int']/float(len(self.contacts_series_frame))
    
        return self.df_avr

    
    
    
    def get_num_int_profile(self):
        """
        Gets number of contacts per residue in selection A
        Returns dataframe with three columns 'segid' 'resid' 'num_int' ordered by resid and segid
        Returns average along the trajectory
        In case of structure with one frame - total number per residue.
        """
        self.get_df_series()
        self.num_int_profile=self.contdf_series.groupby(['A_segid','A_resid']).size().reset_index()
        self.num_int_profile.columns=['segid','resid','num_int']
        self.num_int_profile['num_int']=self.num_int_profile['num_int']/float(len(self.contacts_series_frame))
        return self.num_int_profile
        
## todo          
class trj2int(struct2int):
    """
    Extends struct2cont to import trajectories
    """
    def __init__(self, struct, trj, selA, **kwargs):
        opened_trj=mda.Universe(struct,trj)#,in_memory=True
        super().__init__(opened_trj, selA, **kwargs)

class trj2cont(struct2int):
    """
    Extends struct2cont to import trajectories
    """
    def __init__(self, struct, trj, selA, **kwargs):
        warnings.warn("struct2cont is deprecated,use trj2int", DeprecationWarning)
        opened_trj=mda.Universe(struct,trj)#,in_memory=True
        super().__init__(opened_trj, selA, **kwargs)

        
class struct2cont(struct2int):
    def __init__(self,*args,**kwargs):
        warnings.warn("struct2cont is deprecated,use struct2int", DeprecationWarning)
        super().__init__(*args,**kwargs)
        
        
        


        
   


            
            
            
            
                                      