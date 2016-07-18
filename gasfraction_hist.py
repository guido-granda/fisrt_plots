from __future__ import print_function
import numpy as np
import h5py
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')


#redshift list
redshifts=['iz200','iz174','iz156','iz142','iz131','iz113','iz99','iz78']
#volume list
ivols=['ivol0','ivol10']
nvols=len(ivols)
#directory name 
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
        mcold_t		  =np.empty(1)
        mcold_atom_t      =np.empty(1)
        mcold_atom_bulge_t=np.empty(1)
        mcold_burst_t     =np.empty(1)
        mcold_cooling_t   =np.empty(1)
        mcold_major_t     =np.empty(1)
        mcold_minor_t     =np.empty(1)
        mcold_mol_t       =np.empty(1)
        mstars_bulge_t    =np.empty(1)
        mcold_mol_bulge_t =np.empty(1)
        mcold_recycle_t   =np.empty(1)
        mhalo_t           =np.empty(1)
        mhhalo_t          =np.empty(1)
        mhot_t            =np.empty(1)
        mstars_allburst_t =np.empty(1)
        mstars_bulge_t    =np.empty(1)
        mstars_burst_t    =np.empty(1)
        mstars_disk_t     =np.empty(1)
        vdisk_t           =np.empty(1)
        vhalo_t           =np.empty(1)
        
        title=r'Mcold histogram Mini-Surfs-Velociraptor for gas fraction>0.5 $, '+redshifts[i]
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
                mstars_bulge    =np.array(data.get('mstars_bulge'))#Mass of stars in the bulge (Msolar/h)
                mstars_burst    =np.array(data.get('mstars_burst'))#Stellar mass in ongoing bursts (Msolar/h)
                mstars_disk     =np.array(data.get('mstars_disk'))#Mass of stars in the disk (Msolar/h)
                vdisk=np.array(data.get('vdisk'))
	        vhalo=np.array(data.get('vhalo')) 

                # collecting the data of all the volumes
                mcold_t           =np.append(mcold_t,mcold)
                mcold_atom_t      =np.append(mcold_atom_t,mcold_atom_t)
                mcold_atom_bulge_t=np.append(mcold_atom_bulge_t,mcold_atom_bulge)
                mcold_burst_t     =np.append(mcold_burst_t,mcold_burst)
                mcold_cooling_t   =np.append(mcold_cooling_t,mcold_cooling)
                mcold_major_t     =np.append(mcold_major_t,mcold_major)
                mcold_minor_t     =np.append(mcold_minor_t,mcold_minor)
                mcold_mol_t       =np.append(mcold_mol_t,mcold_mol)
                mstars_bulge_t    =np.append(mstars_bulge_t,mstars_bulge)
                mcold_mol_bulge_t =np.append(mcold_mol_bulge_t,mcold_mol_bulge)
                mcold_recycle_t   =np.append(mcold_recycle_t,mcold_recycle_t)
                mhalo_t           =np.append(mhalo_t,mhalo)
                mhhalo_t          =np.append(mhhalo_t,mhhalo)
                mhot_t            =np.append(mhot_t,mhot)
                mstars_allburst_t =np.append(mstars_allburst_t,mstars_allburst)
                mstars_bulge_t    =np.append(mstars_bulge_t,mstars_bulge)
                mstars_burst_t    =np.append(mstars_burst_t,mstars_burst)
                mstars_disk_t     =np.append(mstars_disk_t,mstars_disk)
                vdisk_t           =np.append(vdisk_t,vdisk)
	        vhalo_t           =np.append(vhalo_t,vhalo) 
            direc='./Gonzalez15.VELOCIraptor/'




        total_mcold   =(mcold_t+mcold_atom_bulge+mcold_cooling+mcold_major+mcold_minor+mcold_mol_bulge+mcold_recycle)[1:]
        total_mass    =total_mcold+mhot[1:]
        fraction      =total_mcold/total_mass
        

#        mcold_t           =
#        mcold_atom_t      =
#        mcold_atom_bulge_t=
#        mcold_burst_t     =
#        mcold_cooling_t   =
#        mcold_major_t     =
#        mcold_minor_t     =
#        mcold_mol_t       =
#        mstars_bulge_t    =
#        mcold_mol_bulge_t =
#        mcold_recycle_t   =
#        mhalo_t           =
#        mhhalo_t          =
#        mhot_t            =
#        mstars_allburst_t =
#        mstars_bulge_t    =
#        mstars_burst_t    =
#        mstars_disk_t     =
#        vdisk_t           =
#        vhalo_t           =
#        mcold_atom_bulge_t=mcold_atom_bulge_t[1:]
#        mcold_mol_bulge_t =mcold_mol_bulge_t[1:]
#        mstars_bulge_t    =mstars_bulge_t[1:]
#        mstars_disk_t     =mstars_disk_t[1:]
#        mchalo_t          =mchalo_t[1:]
#        mhalo_t           =mhalo_t[1:]
#        mhhalo_t          =mhhalo_t[1:]
        vhalo_t           =vhalo_t[1:]
        vdisk_t           =vdisk_t[1:]
#        mbulge_t=mcold_atom_bulge_t+mcold_mol_bulge_t+mstars_bulge_t
#        mtotal  =mbulge_t+mstars_disk_t+mchalo_t+mhalo_t+mhhalo_t
       # ranges
        vhalo_t_min=np.amin(vhalo_t)
        vhalo_t_max=np.amax(vhalo_t)
        vdisk_t_min=np.amin(vdisk_t)
        vdisk_t_max=np.amax(vdisk_t)
        
        print(u'The vhalo goes from %3.5f to %3.5f' %(vhalo_t_min,vhalo_t_max))
        print(u'The vdisk goes from %3.5f to %3.5f \n' %(vdisk_t_min,vdisk_t_max))
 # finding range of values
        n_data1=np.shape(mcold_mol_t)[0]

	n_histogram=int(round(np.sqrt(n_data1)))#/4.0
	left=vdisk_t_min
	right=vdisk_t_max
        idx=(fraction>0.5)#M_Sol/h
        
        dist_bin   =(right-left)/(n_histogram+1)
        l_edge     =left-dist_bin/2.0
        r_edge     =right+dist_bin/2.0
        bin_edges  =np.arange(l_edge,r_edge+dist_bin,dist_bin)
        bin_centers=bin_edges+dist_bin
        bin_centers=bin_centers[1:]
        binwidths  =np.diff(bin_edges)
	################# histogram Vdisk############################3
        hist1, bin_edges1= np.histogram(vdisk_t[idx], bin_edges)
        hist2, bin_edges2= np.histogram(vdisk_t, bin_edges)
	idx1=(hist1>0) 
        idx2=(hist2>0)

        #idx2=(hist1[idx] >0) 
	# over the volume
        hist1=hist1/(volume_tree*nvol)
        hist2=hist2/(volume_tree*nvol)

	#error1=np.sqrt(hist1)/volume_tree
	# correct for log(vmax)
        hist1=hist1/np.diff(bin_edges)
        hist2=hist2/np.diff(bin_edges)

        axs.plot(bin_centers[idx1],hist1[idx1], label=r'$Mcold molecular >10^{7} $', marker='.', linestyle='-', markersize=4, c='r')		  
        axs.plot(bin_centers[idx2],hist2[idx2], label=r'Mcold molecular all ', marker='.', linestyle='-', markersize=4, c='k')
        handles, labels = axs.get_legend_handles_labels()
        axs.legend(handles, labels, numpoints=1, loc='upper right') #, prop={'size': 'x-small'})
        #axs.set_ylim(1E-1,1E1)
        axs.set_xscale('log')
        axs.set_yscale('log')
        axs.set_xlabel(r'$V_{disk}[Km/s]$)')
        axs.set_ylabel(r'$dn/d V_{disk} [h^{3}Mpc^{-3}skm^{-1}]$')
        axs.set_title(title)
        plt.savefig(title+'.png',dpi=dpi,format=fileformat)
        

