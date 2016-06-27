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
                        'text.fontsize':   14,
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
     vdisk_total=np.empty(1)
     vhalo_total=np.empty(1)
     title=r'Mini-Surfs-Velociraptor, '+redshifts[i]
     for j in range(0,len(ivols)):
		
	# set figure
	fileformat = 'png'
	dpi = 300
	fig = plt.figure()

	axs1 = plt.subplot2grid((3,1), (0,0), rowspan=2)
	axs2 = plt.subplot2grid((3,1), (2,0), sharex=axs1)

	axs1.set_xscale('log')
	axs1.set_yscale('log')

		#axs1.axes.get_xaxis().set_ticks([])

	axs2.set_xscale('log')
	axs2.set_yscale('log')

	direc=direc+redshifts[i]+'/'+ivols[j]+'/galaxies.hdf5'
	##Reading the data from hdf5 file ####
        with h5py.File(direc,'r') as hf:
         data = hf.get('Output001')
		 #  np_data = np.array(data)
		 #  print('Shape of the array Output001: \n', np_data.shape)
		 #  print('List of items of the output: \n',data.items())
	 times =hf.get("Output_Times")
	 vdisk=np.array(data.get('vdisk'))
	 vhalo=np.array(data.get('vhalo'))

		 #volume 144703.125
		# Checking Vdisk data and Vhalo data
		 #   n_data1=vdisk.shape[0]# they are both the same number
		 #   n_data2=vhalo.shape[0]# check above
	 np.append(vdisk_total,vdisk)
	 np.append(vhalo_total,vhalo)
  vdisk_t=vdisk_total[1::]
  vhalo_t=vhalo_total[1::]
  vdisk_min=np.amin(vdisk_t)
  vdisk_max=np.amax(vdisk_t)
  print(u'The Vdisk goes from %3.5f to %3.5f' %(vdisk_min,vdisk_max))
  vhalo_min=np.amin(vhalo_t)
  vhalo_max=np.amax(vhalo_t)
  print(u'The Vhalo goes from %3.5f to %3.5f' %(vhalo_min,vhalo_max))
# finding range of values
  minimum=np.minimum(vdisk_min,vhalo_min)
  maximum=np.maximum(vdisk_max,vhalo_max)


  n_histogram=int(round(np.sqrt(n_data1)))#/4.0
	    left=minimum
	    right=maximum

	    dist_bin   =(right-left)/(n_histogram+1)
	    l_edge     =minimum-dist_bin/2.0
	    r_edge     =maximum+dist_bin/2.0
	    bin_edges  =np.arange(l_edge,r_edge+dist_bin,dist_bin)
	    bin_centers=bin_edges+dist_bin
	    bin_centers=bin_centers[1:]

	 ################# histogram Vdisk############################3
	    hist1, bin_edges1= np.histogram(vdisk, bin_edges)
	    idx1=(hist1 >0)
	 # over the volume
	    hist1=hist1/volume_tree
	    error1=np.sqrt(hist1)/volume_tree
	 # correct for dV 
	    hist1=hist1/np.diff(bin_edges)


	 ################## hist vhalom ###############################
	    hist2, bin_edges2= np.histogram(vhalo, bin_edges)
	    idx2=(hist2 >0)
	 # over the volume
	    hist2=hist2/volume_tree
	    error2=np.sqrt(hist2)/volume_tree
	 # correct for log(vmax)
	    hist2=hist2/np.diff(bin_edges)

	 # joining restrictions
	    idx3= idx1 & idx2
	 # ratio of velocities
	    ratio=hist1[idx3]/hist2[idx3]
	    error_ratio=ratio*(error1[idx3]/hist1[idx3]+error2[idx3]/hist2[idx3])

	########################################

	    axs1.set_ylim(1E-8,1E-2)
	    axs2.set_xlim(1E1,2E3)
	 #top plot
	     
	    p1_1 = axs1.plot(bin_centers[idx1],hist1[idx1], label='Vdisk', marker='.', linestyle='-', markersize=4, c='k')
	    p1_2 = axs1.plot(bin_centers[idx2],hist2[idx2], label='Vhalo', marker='.', linestyle='-', markersize=4, c='r' )
	    p2_1 = axs2.plot(bin_centers[idx3],ratio,marker='.',linestyle='-',markersize='4',c='b')
	    handles, labels = axs1.get_legend_handles_labels()
	    axs1.legend(handles, labels, numpoints=1, loc='upper right') #, prop={'size': 'x-small'})
	    axs2.set_ylim(1E-1,1E1)
	    axs2.set_xlabel(r'$V_{disk,halo}[km/s]$)')
	    axs2.set_ylabel(r'ratio Vdisk/Vhalo')
	    axs1.axvline(x=20,ymin=0,ymax=1,color='k',linestyle='--')
	    axs2.axvline(x=20,ymin=0,ymax=1,color='k',linestyle='--')
	    axs2.axhline(y=2,xmin=0,xmax=1,color='k',linestyle='--')
	    axs2.axhline(y=1,xmin=0,xmax=1,color='k',linestyle='--')
	    axs2.axhline(y=0.5,xmin=0,xmax=1,color='k',linestyle='--')
	    axs1.set_ylabel(r'$dn/dV_{disk,halo} [h^3 Mpc^{-3}]$')
	    axs1.set_title(title)
	    fig.savefig(titlei+'.png',dpi=dpi,format=fileformat)
	    direc='./Gonzalez15.VELOCIraptor/'
	# botton plot 


