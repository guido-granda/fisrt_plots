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
        vdisk_total=np.empty(1)
        vhalo_total=np.empty(1)
        title=r'Vdisk vs Vhalo Mini-Surfs-Velociraptor, '+redshifts[i]
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
                vdisk=np.array(data.get('vdisk'))
                vhalo=np.array(data.get('vhalo'))

     #volume 144703.125
    # Checking Vdisk data and Vhalo data
     #   n_data1=vdisk.shape[0]# they are both the same number
     #   n_data2=vhalo.shape[0]# check above
                vdisk_total=np.append(vdisk_total,vdisk)
                vhalo_total=np.append(vhalo_total,vhalo)
            direc='./Gonzalez15.VELOCIraptor/'
        
        vdisk_t=vdisk_total[1:]
        vhalo_t=vhalo_total[1:]
        vdisk_min=np.amin(vdisk_t)
        vdisk_max=np.amax(vdisk_t)
        print(u'The Vdisk goes from %3.5f to %3.5f' %(vdisk_min,vdisk_max))
        vhalo_min=np.amin(vhalo_t)
        vhalo_max=np.amax(vhalo_t)
        print(u'The Vhalo goes from %3.5f to %3.5f' %(vhalo_min,vhalo_max))
 # finding range of values
        minimum=np.minimum(vdisk_min,vhalo_min)
        maximum=np.maximum(vdisk_max,vhalo_max)
        n_data1=np.shape(vdisk_t)[0]

	     
        axs.plot(vhalo_t,vdisk_t, 'ro',markersize =2)
        #handles, labels = axs.get_legend_handles_labels()
        #axs.legend(handles, labels, numpoints=1, loc='upper right') #, prop={'size': 'x-small'})
        #axs.set_ylim(1E-1,1E1)
        axs.set_xlabel(r'$V_{halo}[km/s]$)')
        axs.set_ylabel(r'$V_{disk}[km/s]$')
        axs.set_title(title)
        plt.savefig(title+'.png',dpi=dpi,format=fileformat)
        

