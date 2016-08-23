from __future__ import print_function
import numpy as np
import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


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
dm_mass = 5.97e9 # Msun/h
halo_mass_min=100*dm_mass # halo mass treshold



for i in range(0,len(redshifts)):
        #Setting collecting arrays
        rdisk_t		  =np.empty(1)
        vdisk_t           =np.empty(1)        
        title=r'Distribution of sizes, '+redshifts[i]
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
                rdisk      =np.array(data.get('rdisk')) # Mass of cold gas in the disk of the galaxy             
                vdisk      =np.array(data.get('vdisk')) #
                # collecting the data of all the volumes
                rdisk_t    =np.append(rdisk_t,rdisk)
                vdisk_t    =np.append(vdisk_t,vdisk)
                print('The dimension of vdisk_t is %2.0f: \n' %(np.shape(vdisk_t)[0]-1))
            direc='./Gonzalez15.VELOCIraptor/'

        rdisk_t =rdisk_t[1:]
        vdisk_t =vdisk_t[1:]

       # ranges
        rdisk_t_min=np.amin(rdisk_t)
        rdisk_t_max=np.amax(rdisk_t)
        
        print(u'The rdisk goes from %3.5f to %3.5f' %(rdisk_t_min,rdisk_t_max))
        # finding histogram 
        n_data1=np.shape(rdisk_t)[0]

	n_histogram=int(round(np.sqrt(n_data1)))
	left  =rdisk_t_min
	right =rdisk_t_max
       # Vdisk criteria selction 
        idb1=(vdisk_t>500 )
        
        dist_bin   =(right-left)/(n_histogram+1)
        l_edge     =left-dist_bin/2.0
        r_edge     =right+dist_bin/2.0
        bin_edges  =np.arange(l_edge,r_edge+dist_bin,dist_bin)
        bin_centers=bin_edges+dist_bin
        bin_centers=bin_centers[1:]
        binwidths  =np.diff(bin_edges)
	################# histogram Vdisk############################
        hist1, bin_edges1= np.histogram(rdisk_t[idb1], bin_edges)
        hist2, bin_edges2= np.histogram(rdisk_t, bin_edges)

	#############################################################
        #idx1=(hist1>0) 
        #idx2=(hist2>0)
        #idxh= idx1 & idx2
	# over the volume
        #  hist1=hist1/(volume_tree*nvols)
	# correct for dv 
        # hist1=hist1/np.diff(bin_edges)
        # plots
        #axs.set_xrange([1.0e1,1.0e3])
        axs.plot(bin_centers,hist1,label=r'$V_{disk} >  500 (km/s) $', marker='.', linestyle='-', markersize=4, c='r')
	axs.plot(bin_centers,hist2,label=r'total', marker='.', linestyle='-', markersize=4, c='b')
	  
        handles, labels = axs.get_legend_handles_labels()
        axs.legend(handles, labels, numpoints=1, loc='upper right') #, prop={'size': 'x-small'})
        axs.set_ylim([-1.0,1.2*np.amax(hist2)])
        axs.set_xlabel(r'$r_{disk}[Mpc/h]$)')
        axs.set_ylabel(r'$N$')
        axs.set_title(title)
        plt.savefig(title+'.png',dpi=dpi,format='png')
