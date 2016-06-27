from __future__ import print_function
import numpy as np
import h5py
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')


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
        mcold_t=np.empty(1)
        title=r'Mcold histogram Mini-Surfs-Velociraptor, '+redshifts[i]
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
                mcold_atom=np.array(data.get('mcold_atom'))
     #volume 144703.125
    # Checking Vdisk data and Vhalo data
     #   n_data1=vdisk.shape[0]# they are both the same number
     #   n_data2=vhalo.shape[0]# check above
                mcold_t=np.append(mcold_t,mcold_atom)
            direc='./Gonzalez15.VELOCIraptor/'
        
        mcold_t=mcold_t[1:]
        mcold_min=np.amin(mcold_t)
        mcold_max=np.amax(mcold_t)
        print(u'The mcold goes from %3.5f to %3.5f' %(mcold_min,mcold_max))
 # finding range of values
        n_data1=np.shape(mcold_t)[0]

	n_histogram=int(round(np.sqrt(n_data1)))#/4.0
	left=np.log10(mcold_min+1.0)
	right=np.log10(mcold_max)

	dist_bin   =(right-left)/(n_histogram+1)
	l_edge     =left-dist_bin/2.0
	r_edge     =right+dist_bin/2.0
	bin_edges  =np.arange(l_edge,r_edge+dist_bin,dist_bin)
	bin_centers=bin_edges+dist_bin
	bin_centers=bin_centers[1:]
        binwidths  =np.diff(10**bin_edges)
	################# histogram Vdisk############################3
	hist1, bin_edges1= np.histogram(np.log10(mcold_t+1.0), bin_edges)
	idx1=(hist1 >0)
	# over the volume
	hist1=hist1/volume_tree
	#error1=np.sqrt(hist1)/volume_tree
	# correct for log(vmax)
	hist1=np.log10(hist1[idx1]/np.diff(np.log(bin_edges)))
        p1_1 = axs.plot(bin_centers[idx1],hist1, label='Mcold', marker='.', linestyle='-', markersize=4, c='k')		  
        #handles, labels = axs.get_legend_handles_labels()
        #axs.legend(handles, labels, numpoints=1, loc='upper right') #, prop={'size': 'x-small'})
        #axs.set_ylim(1E-1,1E1)
        axs.set_xlabel(r'$log(M_{cold}[M_\odot/h])$)')
        axs.set_ylabel(r'$log(dn/dln M_{cold}/h^{3}Mpc^{-3}$')
        axs.set_title(title)
        plt.savefig(title+'.png',dpi=dpi,format=fileformat)
        

