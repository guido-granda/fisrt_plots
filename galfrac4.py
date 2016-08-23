from __future__ import print_function
import numpy as np
import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

#redshift list
redshifts=['iz200','iz174','iz156','iz142','iz131','iz113','iz99','iz78']
ivols=['ivol0','ivol10']
direc='./Gonzalez15.VELOCIraptor/'

# set font properties
params = {      'axes.labelsize': 14,
                        'axes.labelweight': 'normal',
                        'axes.titlesize':  16, # Plot Title size
                        'font.size':   14,
                        'legend.fontsize': 14,
                        'xtick.labelsize': 14, # Ticklabels
                        'ytick.labelsize': 14,
                        'font.family' : 'sans-serif',
                        'text.usetex': True} # True for Latex font
                  
plt.rcParams.update(params)


#some parameters
volume_tree=144703.125 #(Mpc/h)^3
box_size= 210.0 # Mpc/h    




for i in range(0,len(redshifts)):
        #Setting collecting arrays
        mcold_t           =np.empty(1)
        mcold_atom_t      =np.empty(1)
        mcold_atom_bulge_t=np.empty(1)
        mcold_burst_t     =np.empty(1)
        mcold_cooling_t   =np.empty(1)
        mcold_major_t     =np.empty(1)
        mcold_minor_t     =np.empty(1)
        mcold_mol_t       =np.empty(1)
        mcold_mol_bulge_t =np.empty(1)
        mcold_recycle_t   =np.empty(1)
        mhalo_t           =np.empty(1)
        mhhalo_t          =np.empty(1)
        mhot_t            =np.empty(1)
        mstars_allburst_t =np.empty(1)
        mstars_burst_t    =np.empty(1)
        mstars_disk_t     =np.empty(1)
        mstars_bulge_t    =np.empty(1)
        vdisk_t           =np.empty(1)
        vhalo_t           =np.empty(1)
        mstars_instab_t   =np.empty(1)

        title='Hot mass fraction vs Stellar mass, all volumes,'+redshifts[i]
        for j in range(0,len(ivols)):
        # set figure
            fileformat = 'png'
            dpi = 300
            fig,axs = plt.subplots()
            direc=direc+redshifts[i]+'/'+ivols[j]+'/galaxies.hdf5'
            ##Reading the data from hdf5 file ####
            with h5py.File(direc,'r') as hf:
                data = hf.get('Output001')
                times =hf.get("Output_Times")
                mcold           =np.array(data.get('mcold')) # Mass of cold gas in the disk of the galaxy         
                mcold_atom      =np.array(data.get('mcold_atom')) #Mass of H gas mass in the disk of the galaxy
                mcold_atom_bulge=np.array(data.get('mcold_atom_bulge'))# Mass of H gas in the bulge of the galaxy
                mcold_burst     =np.array(data.get('mcold_burst'))#Mass of cold gas remaining in ongoing bursts (Msolar/h)
                mcold_cooling   =np.array(data.get('mcold_cooling'))#Cold gas Mass in the galaxy that comes from cooling
                mcold_major     =np.array(data.get('mcold_major'))#Cold gas Mass in the galaxy that comes from major mergers (Msolar/h)
                mcold_minor     =np.array(data.get('mcold_minor'))# Cold gas Mass in the galaxy that comes from minor mergers (Msolar/h).
                mcold_mol       =np.array(data.get('mcold_mol'))#Mass of molecularhydrogen gas mass in the disk of the galaxy (Msolar/h)
                mstars_bulge    =np.array(data.get('mstars_bulge'))#Mass of stars in the bulge (Msolar/h)
                mcold_mol_bulge =np.array(data.get('mcold_mol_bulge'))#Mass of molecular hydrogen gas mass in the bulge of the galaxy (Msolar/h)
                mcold_recycle   =np.array(data.get('mcold_recycle'))# Cold gas Mass in the galaxy that comes from recycling of old stars (Msolar/h).
                mhalo           =np.array(data.get('mhalo'))#Mass of the halo in which the galaxy formed
                mhhalo          =np.array(data.get('mhhalo'))#Mass of the host halo at this output time (Msolar/h)
                mhot            =np.array(data.get('mhot'))#Mass of hot gas in the halo (Msolar/h)
                mstars_allburst =np.array(data.get('mstars_allburst'))#Stellar mass in all bursts (Msolar/h)
                mstars_burst    =np.array(data.get('mstars_burst'))#Stellar mass in ongoing bursts (Msolar/h)
                mstars_disk     =np.array(data.get('mstars_disk'))#Mass of stars in the disk (Msolar/h)
                mstars_instab   =np.array(data.get('mstars_instability'))# mass of stars formed in an instability-triggered burst
                vdisk           =np.array(data.get('vdisk'))
                vhalo           =np.array(data.get('vhalo'))

                # collecting the data of all the volumes
                mcold_t           =np.append(mcold_t,mcold)
                mcold_atom_t      =np.append(mcold_atom_t,mcold_atom)
                mcold_atom_bulge_t=np.append(mcold_atom_bulge_t,mcold_atom_bulge)
                mcold_burst_t     =np.append(mcold_burst_t,mcold_burst)
                mcold_cooling_t   =np.append(mcold_cooling_t,mcold_cooling)
                mcold_major_t     =np.append(mcold_major_t,mcold_major)
                mcold_minor_t     =np.append(mcold_minor_t,mcold_minor)
                mcold_mol_t       =np.append(mcold_mol_t,mcold_mol)
                mcold_mol_bulge_t =np.append(mcold_mol_bulge_t,mcold_mol_bulge)
                mcold_recycle_t   =np.append(mcold_recycle_t,mcold_recycle_t)
                mhalo_t           =np.append(mhalo_t,mhalo)
                mhhalo_t          =np.append(mhhalo_t,mhhalo)
                mhot_t            =np.append(mhot_t,mhot)
                mstars_allburst_t =np.append(mstars_allburst_t,mstars_allburst)
                mstars_bulge_t    =np.append(mstars_bulge_t,mstars_bulge)
                mstars_burst_t    =np.append(mstars_burst_t,mstars_burst)
                mstars_disk_t     =np.append(mstars_disk_t,mstars_disk)
                mstars_instab_t   =np.append(mstars_instab_t,mstars_instab)    
                vdisk_t           =np.append(vdisk_t,vdisk)
                vhalo_t           =np.append(vhalo_t,vhalo)
                print('The dimension of vdisk_t is %2.0f: \n' %(np.shape(vdisk_t)[0]-1))
            direc='./Gonzalez15.VELOCIraptor/'


        
        mcold_t           =mcold_t[1:]
        mcold_atom_t      =mcold_atom_t[1:]
        mcold_atom_bulge_t=mcold_atom_bulge_t[1:]
        mcold_burst_t     =mcold_burst_t[1:]
        mcold_cooling_t   =mcold_cooling_t[1:]
        mcold_major_t     =mcold_major_t[1:]
        mcold_minor_t     =mcold_minor_t[1:]
        mstars_disk_t     =mstars_disk_t[1:]
        mstars_bulge_t    =mstars_bulge_t[1:]
        mhot_t            =mhot_t[1:]      
  
        mcold_total       =mcold_t + mcold_burst_t
        mdisk_total       =mstars_disk_t + mcold_t
        mbulge_total      =mstars_bulge_t+ mcold_burst_t
        mtotal            =mdisk_total+mbulge_total
 
        mstars_total      =mstars_bulge_t+mstars_disk_t
        #mcold_total       =mcold_t+mcold_atom_t+mcold_atom_bulge_t+mcold_burst_t+mcold_cooling_t+mcold_major_t+mcold_minor_t
        #mtotal            =mstars_total+mcold_total 
        id0=(mtotal>0)
        #fcold=mcold_total[id0]/mtotal[id0]
        fhot=mhot_t/mtotal

        axs.plot(mstars_total[id0],fhot[id0],'oy',label=r'$f_{hot}$', markersize=2)
       #r"$f_{cold}\geq 0.5$",markersize =2)
       # axs.plot(vhalo_t[idx2],vdisk_t[idx2],'ob',label=r"$f_{cold}< 0.5$",markersize =2)
        handles, labels = axs.get_legend_handles_labels()
        axs.legend(handles, labels, numpoints=1, loc='upper left') 

        #handles, labels = axs.get_legend_handles_labels()
        #axs.legend(handles, labels, numpoints=1, loc='upper right') #, prop={'size': 'x-small'})
        axs.set_ylim(-0.1,1.1)
        axs.set_xscale('log')
        axs.set_xlabel(r'$M_{*}(M_{\odot}/h)$')
        axs.set_ylabel(r'$f_{hot}$')
        axs.set_title(title)
        axs.grid('True',linestyle='-',which='major' )
        plt.savefig(title,dpi=dpi,format=fileformat)
        

