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
       # mcold_atom_t=np.empty(1)
       # mcold_mol_t  =np.empty(1)
        mcold_atom_bulge_t=np.empty(1)
        mcold_mol_bulge_t =np.empty(1)
        mstars_bulge_t    =np.empty(1)
        mstars_disk_t     =np.empty(1)
        mchalo_t          =np.empty(1)
        mhalo_t           =np.empty(1)
        mhhalo_t          =np.empty(1)
        vdisk_t=np.empty(1)
        vhalo_t=np.empty(1)
        
        title=r'Mcold histogram Mini-Surfs-Velociraptor for $M_{mol}>10^{7}$, '+redshifts[i]
        for j in range(0,len(ivols)):
        # set figure
            fileformat = 'png'
            dpi = 300
            fig,axs = plt.subplots()
            direc=direc+redshifts[i]+'/'+ivols[j]+'/galaxies.hdf5'
            ##Reading the data from hdf5 file ####
            with h5py.File(direc,'r') as hf:
                data = hf.get('Output001')
               #  np_data = np.array(data)
               #  print('Shape of the array Output001: \n', np_data.shape)
               #  print('List of items of the output: \n',data.items())
                times =hf.get("Output_Times")
               # vdisk=np.array(data.get('vdisk'))
               # vhalo=np.array(data.get('vhalo'))
                mcold_atom_bulge=np.array(data.get('mcold_atom_bulge'))
                mcold_mol_bulge =np.array(data.get('mcold_mol_bulge'))
                mstars_bulge    =np.array(data.get('mstars_bulge'))
                mstars_disk     =np.array(data.get('mstars_disk'))
                mchalo          =np.array(data.get('mchalo'))
                mhalo           =np.array(data.get('mhalo'))
                mhhalo          =np.array(data.get('mhhalo'))
                vdisk=np.array(data.get('vdisk'))
	        vhalo=np.array(data.get('vhalo')) 
     #volume 144703.125
    # Checking Vdisk data and Vhalo data
     #   n_data1=vdisk.shape[0]# they are both the same number
     #   n_data2=vhalo.shape[0]# check above
                mcold_atom_bulge_t=np.append(mcold_atom_bulge_t,mcold_atom_bulge)
                mcold_mol_bulge_t =np.append(mcold_mol_bulge_t,mcold_mol_bulge)
                mstars_bulge_t    =np.append(mstars_bulge_t,mstars_bulge)
                mstars_disk_t     =np.append(mstars_disk_t,mstars_disk)
                mchalo_t          =np.append(mchalo_t,mchalo)
                mhalo_t           =np.append(mhalo_t,mhalo)
                mhhalo_t          =np.append(mhhalo_t,mhhalo)
                vhalo_t=np.append(vhalo_t,vhalo)
                vdisk_t=np.append(vdisk_t,vdisk)
            direc='./Gonzalez15.VELOCIraptor/'
        
        mcold_atom_bulge_t=mcold_atom_bulge_t[1:]
        mcold_mol_bulge_t =mcold_mol_bulge_t[1:]
        mstars_bulge_t    =mstars_bulge_t[1:]
        mstars_disk_t     =mstars_disk_t[1:]
        mchalo_t          =mchalo_t[1:]
        mhalo_t           =mhalo_t[1:]
        mhhalo_t          =mhhalo_t[1:]
        vhalo_t           =vhalo_t[1:]
        vdisk_t           =vdisk_t[1:]
        mbulge_t=mcold_atom_bulge_t+mcold_mol_bulge_t+mstars_bulge_t
        mtotal  =mbulge_t+mstars_disk_t+mchalo_t+mhalo_t+mhhalo_t
       # ranges
       # mcold_atom_min=np.amin(mcold_atom_t)
       # mcold_atom_max=np.amax(mcold_atom_t)
      #  mcold_mol_min=np.amin(mcold_mol_t)
      #  mcold_mol_max=np.amax(mcold_mol_t)
        vhalo_t_min=np.amin(vhalo_t)
        vhalo_t_max=np.amax(vhalo_t)
        vdisk_t_min=np.amin(vdisk_t)
        vdisk_t_max=np.amax(vdisk_t)
        
        #print(u'The atomic mcold goes from %3.5f to %3.5f' %(mcold_atom_min,mcold_atom_max))
        #print(u'The molecular mcold goes from %3.5f to %3.5f' %(mcold_mol_min,mcold_mol_max))
        print(u'The vhalo goes from %3.5f to %3.5f' %(vhalo_t_min,vhalo_t_max))
        print(u'The vdisk goes from %3.5f to %3.5f' %(vdisk_t_min,vdisk_t_max))
 # finding range of values
        n_data1=np.shape(mcold_mol_t)[0]

	n_histogram=int(round(np.sqrt(n_data1)))#/4.0
	left=vdisk_t_min
	right=vdisk_t_max
        idx=(mcold_mol_t>1.0e7)#M_Sol/h
        
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
        

