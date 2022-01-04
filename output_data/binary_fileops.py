
"""
File for extracting data from binary generated
by wave packet dynamics and processing

"""

import os
import sys
import numpy as np
from numpy import linalg as LA
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.animation as animation


def get_data_from_bin():
    filptr = open(sys.argv[1], 'r')
    filedata = filptr.readlines()
    path = filedata[0].strip()
    bin_data = np.fromfile(path, dtype=float)  # loading data into array
    total_filesize = int(len(bin_data))     # total no of elements in the file
    # checking whether plotting is to be one/two dimensional
    dimension = int(filedata[1].split()[0])
    # important for 2d, if surface plot/contour plot
    plot_type = filedata[1].split()[1]
    # whether propagation is done numerically/analytically
    propagate_type = filedata[1].split()[2]
    length_each_wp = int(filedata[2].strip())  # size of data for each frame
    ##------------------------------##
    ## In numerical wave packet propagation, for each multidimensional grid 
    ## point,real,imaginary and absolute values are saved leading to 3 times 
    ## the size ##
    if propagate_type == "numeric":
        # total no of frames in the animation
        num_frames = int(total_filesize/(3*length_each_wp))
    else:
        num_frames = int(total_filesize/length_each_wp)
    fps = float(filedata[10].strip())
    stepsize = float(filedata[11].strip())
    timestep_arr = []       ## array for storing times for different frames 
    ###------------------###
    if len(filedata[9].split()) == num_frames:
        for i in range(num_frames):
            timestep_arr.append(stepsize*float(filedata[9].split()[i]))
    else:
        for i in range(num_frames):
            timestep_arr.append(i*stepsize*float(filedata[9].split()[0]))
    ###------------------###
    total_time = timestep_arr[len(timestep_arr)-1]  ## total time for dynamics
    save_location = filedata[6].strip()     ## location where animation is to be saved 
    ani_filename = filedata[7].strip()      ## animation filename
    ##-------------------##
    ## accumulating data for all frames ##
    if propagate_type != "numeric":
        full_data = []
        for i in range(num_frames):
            cur_wp_val = []
            start_id = i*length_each_wp
            end_id = (i+1)*length_each_wp
            for j in range(start_id, end_id):
                cur_wp_val.append(bin_data[j])
            full_data.append(cur_wp_val)
    else:
        full_data,full_data_chk_w_anal = [],[]
        full_data_real,full_data_real_chk_w_anal = [],[]
        full_data_imag,full_data_imag_chk_w_anal = [],[]
        full_data_prob_density,full_data_prob_density_chk_w_anal = [],[]
        for i in range(num_frames):
            cur_wp_val = []
            cur_wp_val_real = []
            cur_wp_val_imag = []
            cur_wp_val_sqr = []
            start_id = i*3*length_each_wp
            end_id = (i+1)*3*length_each_wp
            for j in range(start_id, end_id, 3):
                cur_wp_val.append(bin_data[j+2])
                cur_wp_val_real.append(bin_data[j])
                cur_wp_val_imag.append(bin_data[j+1])
                cur_wp_val_sqr.append(bin_data[j+2]**2)
            if i%5 !=0:
                full_data.append(cur_wp_val)
                full_data_real.append(cur_wp_val_real)
                full_data_imag.append(cur_wp_val_imag)
                full_data_prob_density.append(cur_wp_val_sqr)
            elif i%5 ==0:
                full_data.append(cur_wp_val)
                full_data_real.append(cur_wp_val_real)
                full_data_imag.append(cur_wp_val_imag)
                full_data_prob_density.append(cur_wp_val_sqr)
                full_data_chk_w_anal.append(cur_wp_val)
                full_data_real_chk_w_anal.append(cur_wp_val_real)
                full_data_imag_chk_w_anal.append(cur_wp_val_imag)
                full_data_prob_density_chk_w_anal.append(cur_wp_val_sqr)
#                np.savetxt("wp_1_2_abs_chk_w_anal.dat",full_data_chk_w_anal,fmt='%20.12e')
#               np.savetxt("wp_1_2_real_chk_w_anal.dat",full_data_real_chk_w_anal,fmt='%20.12e')
#                np.savetxt("wp_1_2_imag_chk_w_anal.dat",full_data_imag_chk_w_anal,fmt='%20.12e')
#                np.savetxt("wp_1_2_sqr_chk_w_anal.dat",full_data_prob_density_chk_w_anal,fmt='%20.12e')
#        sys.exit()
    ##--------------------##
    if dimension == 1:
        X_l = float(filedata[3].strip())
        X_u = float(filedata[4].strip())
        X_numpts = float(filedata[5].split()[0])
        delta_X = (X_u-X_l)/(X_numpts-1)
        X_grid = [X_l+delta_X*k for k in range(X_numpts)]
        if len(X_grid) != length_each_wp:
            print("Data type is not correct")
            sys.exit()
        else:
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.set_xlabel("Time",fontsize=15)
            ax.set_ylabel("$C_{t}$",fontsize=15)
            ax.set_title("Time-correlation function",fontsize = 15)
            ax.xaxis.set_tick_params(labelsize=20)
            ax.yaxis.set_tick_params(labelsize=20)
            ax.tick_params(direction='out',length=6,width=2,color='black')
            fig.tight_layout()
            Writer = animation.writers['ffmpeg']
            writer = Writer(fps=fps,metadata=dict(artist='Me',bitrate=1800))
            ims = []
            for i in range(num_frames):
                ims.append(ax.plot(X_grid,full_data[i],linestyle = '-',color = 'teal',marker = 'o', label = str(timestep_arr[i])+"  fs"))
                ax.legend(loc='best')
            im_ani = animation.ArtistAnimation(fig,ims,interval = total_time/fps,repeat=False,blit=True)
            im_ani.save(save_location+ani_filename,writer=writer)
    elif dimension == 2:
        X_l = float(filedata[3].split()[0])
        X_u = float(filedata[3].split()[1])
        X_label = filedata[3].split()[2] 
        Y_label = filedata[4].split()[2] 
        Y_l = float(filedata[4].split()[0])
        Y_u = float(filedata[4].split()[1])
        X_numpts = int(filedata[5].split()[0])
        Y_numpts = int(filedata[5].split()[1])
        delta_X, delta_Y = (X_u-X_l)/(X_numpts-1), (Y_u-Y_l)/(Y_numpts-1)
        X_grid = [X_l+k*delta_X for k in range(X_numpts)]
        Y_grid = [Y_l+k*delta_Y for k in range(Y_numpts)]
        X,Y = np.meshgrid(X_grid,Y_grid)
        full_data_reshape = np.reshape(full_data,(num_frames,X_numpts,Y_numpts))
        print(len(full_data_reshape))
        Writer = animation.writers['ffmpeg']
        writer = Writer(fps, metadata=dict(artist='Me'), bitrate=1800)
        if plot_type == "surface":
            fig = plt.figure()
            ax = fig.add_subplot(111,projection='3d')
            plot = [ax.plot_surface(X,Y,full_data_reshape[0,:,:].T,color = 'turquoise',alpha=0.5,rstride=1,cstride=1)]
            ax.set_title(str(timestep_arr[0])+"  fs",loc = 'left')
            ax.set_xlabel(X_label,fontsize=15)
            ax.set_ylabel(Y_label,fontsize=15)
#            ax.set_zlabel("$|\psi|^{2}(t)$")
            ax.set_zlabel("$\psi_{imag}(t)$")
            ax.set_zlim(-0.2,0.2)
            def update_frame(frame_no,full_data_reshape,plot):
                plot[0].remove()
                plot[0] = ax.plot_surface(X,Y,full_data_reshape[frame_no,:,:].T,color = 'turquoise',alpha=0.5,rstride=1,cstride=1)
                ax.set_title(str(timestep_arr[frame_no])+"  fs",loc='left')
                ax.set_xlabel(X_label)
                ax.set_ylabel(Y_label)
                ax.set_zlabel("$\psi_{imag}(t)$")
                ax.set_zlim(-0.2,0.2)
            ani = animation.FuncAnimation(fig,update_frame,frames=num_frames,fargs=(full_data_reshape,plot),blit=False,repeat=False)
            ani.save(save_location+ani_filename,writer=writer)
        elif plot_type == "contour":
            fig,ax = plt.subplots()
            ax.set_xlabel(X_label,fontsize=20)
            ax.set_ylabel(Y_label,fontsize=20)
            ax.set_title(str(timestep_arr[0])+"  fs",loc = 'left')
            norm = mpl.colors.Normalize(vmin=-0.30, vmax=0.30)
            m = plt.cm.ScalarMappable(cmap=cm.RdBu)
            m.set_array(full_data_reshape)
            m.set_clim(-0.30,0.30)
            plt.colorbar(m,ticks = [-0.30,-0.20,-0.10,0.0,0.10,0.20,0.30])
            def update_frame(frame_no):
                cs = ax.contour(X,Y,full_data_reshape[frame_no,:,:].T,200,linewidth= 20,cmap = cm.RdBu,norm = norm)
                ax.set_xlabel(X_label,fontsize=15)
                ax.set_ylabel(Y_label,fontsize=15)
                ax.set_title(str(timestep_arr[frame_no])+"  fs",loc = 'left')
            ani = animation.FuncAnimation(fig,update_frame,frames=num_frames,blit=False,repeat=False)
            ani.save(save_location+ani_filename,writer=writer)


get_data_from_bin()


