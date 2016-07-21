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
        vdisk_t=np.empty(1)
        vhalo_t=np.empty(1)
        mstars_bulge_t=np.empty(1)
        mstars_disk_t=np.empty(1)

        title=r'Vdisk vs Vhalo Mini-Surfs-Velociraptor, '+redshifts[i]
        for j in range(0,len(ivols)):
        # set figure
            fileformat = 'png'
            dpi = 300
            fig,axs = plt.subplots()
            direc=direc+redshifts[i]+'/'+ivols[j]+'/galaxies.hdf5'
            ##Reading the data from hdf5 file ####
            with h5py.File(direc,'r') as hf:
                data         = hf.get('Output001')
                times        =hf.get("Output_Times")
                vdisk        =np.array(data.get('vdisk'))
                vhalo        =np.array(data.get('vhalo'))
                mstars_disk  =np.array(data.get('mstars_disk'))
                mstars_bulge =np.array(data.get('mstars_bulge'))
                


                vdisk_t      =np.append(vdisk_t,vdisk)
                vhalo_t      =np.append(vhalo_t,vhalo)
                mstars_disk_t=np.append(mstars_disk_t,mstars_disk)
                mstars_bulge_t=np.append(mstars_bulge_t,mstars_bulge)
            direc='./Gonzalez15.VELOCIraptor/'


        
        vdisk_t       =vdisk_t[1:]
        vhalo_t       =vhalo_t[1:]
        mstars_disk_t =mstars_disk_t[1:]
        mstars_bulge_t=mstars_bulge_t[1:]
        mstars_total  = mstars_disk_t + mstars_bulge_t
        id0= (mstars_total>0)
        fbulge=mstars_bulge_t[id0]/mstars_total[id0]
        vhalo_t           =vhalo_t[id0]
        vdisk_t           =vdisk_t[id0]
       


        idx1=(fbulge>=0.5)
        idx2=(fbulge<0.5)

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

	     
        axs.plot(vhalo_t[idx1],vdisk_t[idx1],'or',label=r"$f_bulge\geq 0.5$",markersize =2)
        axs.plot(vhalo_t[idx2],vdisk_t[idx2],'ob',label=r"$f_bulge< 0.5$",markersize =2)
        handles, labels = axs.get_legend_handles_labels()
        axs.legend(handles, labels, numpoints=1, loc='upper right') 

        #handles, labels = axs.get_legend_handles_labels()
        #axs.legend(handles, labels, numpoints=1, loc='upper right') #, prop={'size': 'x-small'})
        #axs.set_ylim(1E-1,1E1)
        axs.set_xlabel(r'$V_{halo}[km/s]$)')
        axs.set_ylabel(r'$V_{disk}[km/s]$')
        axs.set_title(title)
        plt.savefig(title+'.png',dpi=dpi,format=fileformat)
        

