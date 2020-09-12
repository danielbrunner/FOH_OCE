"""gleiche wie immer nur mit einer gauss bell"""



import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap, BoundaryNorm
from shapely.geometry import LineString, Point
import matplotlib.colors as colors
from _utils import format_paper

format_paper()






x_0=0
y_0=0




v_p_W=1450        # velocity of acoustic waves in water [m/s]
v_p_C=3540
v_s_C=v_p_C/np.sqrt(3)
rho_W=1000      # kg/m3a
rho_C=2500      # kg/m3

########################## pressure at the sea surface
T=np.linspace(1,25,250)

tt=np.linspace(0,20,40)


radius=np.linspace(0,800000,4000)




ocean_depths=np.linspace(6000,20000,250)

omegas=2*np.pi/T
thetas=np.linspace(5.1,21.6,10)
heights=np.linspace(0,6000,10)

#plt.ion()

#
# # TODO was ist mit diesen plots .....
# ###############################################
# fff=[4]
#
# fig, axs = plt.subplots(1, 1, sharex=True, sharey=True)
# fig.set_size_inches(15, 10)
# for ff in fff:
#
#     max_value = []
#     max_period = []
#     for ocean_depth in ocean_depths:
#         transversal=np.load("/home/djamel/PHD_projects/FOH_ICE/data_flat_sea_floor/transversal_normal_data/height_ocean_depth_frequency_/transversal/transversal_depth: " + str(int(ocean_depth)) + ".npy")
#         # normal=np.load("/home/djamel/PHD_projects/FOH_ICE/data_flat_sea_floor/transversal_normal_data/height_ocean_depth_frequency_/normal/normal_depth: " + str(int(ocean_depth)) + ".npy")
#         # max_P=np.load("/home/djamel/PHD_projects/FOH_ICE/data_flat_sea_floor/transversal_normal_data/height_ocean_depth_frequency_/max_P/max_P_depth: " + str(int(ocean_depth)) + ".npy")
#
#
#         max_value=np.append(max_value,np.max(np.abs(transversal[-1, ff, :])))
#         max_period=np.append(max_period,np.argmax(np.abs(transversal[-1, ff, :])))
#
#     max_period=max_period.astype(int)
#
#     max_T=T[max_period]
#
#
#
#     points = np.array([max_T, ocean_depths]).T.reshape(-1, 1, 2)
#     segments = np.concatenate([points[:-1], points[1:]], axis=1)
#
#
#     norm = plt.Normalize(max_value.min(), max_value.max())
#     lc = LineCollection(segments, cmap='hot', norm=norm)
#     lc.set_array(max_value)
#     lc.set_linewidth(5)
#     line = axs.add_collection(lc)
#     #fig.colorbar(line, ax=axs)
#     axs.set_xlim(max_T.min()-2,max_T.max()+2)
#     axs.set_ylim(ocean_depths.min()-1000,ocean_depths.max()+1000)
#     plt.xlabel("T [s]")
#     plt.ylabel("ocean depth [m]")
#     plt.legend("")
# plt.show()
# #########################################################################
#
#
#
# #plt.subplot(131)
#     #plt.imshow(np.flipud(transversal[:,2,:]), extent=(T[-1], T[0], heights[0], heights[-1]), aspect='auto', cmap='gnuplot')
#     #plt.plot(T,transversal[-1,5,:])
#     # plt.subplot(132)
#     # plt.imshow(np.flipud(normal[:,2,:]), extent=(T[-1], T[0], heights[0], heights[-1]), aspect='auto', cmap='gnuplot')
#     # plt.subplot(133)
#     # plt.imshow(np.flipud(max_P[:,2,:]), extent=(T[-1], T[0], heights[0], heights[-1]), aspect='auto', cmap='gnuplot')
#     # #plt.colorbar()
# #plt.show()
#
#     #plt.pause(0.1)
#
#
transversal_tot = np.zeros((250, 250),dtype=complex)
normal_tot = np.zeros((250, 250),dtype=complex)
pressure_tot = np.zeros((250, 250))
R_tot = np.zeros((250, 2250))
kks= np.zeros((250, 250))

sigma=50000
# ampl=100
# gauss=ocean_depths[40]-ampl*np.exp(-(radius-radius[-1]/2)**2/sigma**2)




for ttt, theta_0 in enumerate(thetas):
    print(theta_0)





    for hhe, height in enumerate(heights):
        print(height)

        working_dir = "/home/djamel/PHD_projects/FOH_ICE/gauss_hill_data"
        working_dir = working_dir + "/theta_0_{}_height_{}".format(np.round(theta_0,1),int(height))
        if os.path.exists(working_dir):
            print("exist")
            continue

        if not os.path.exists(working_dir):
            os.mkdir(working_dir)



        for odd,ocean_depth in enumerate(ocean_depths):
            print(odd)




            #profile_upward=np.linspace(ocean_depth,ocean_depth-height,2000)
            #profile_downward=np.linspace(ocean_depth-height,ocean_depth,2000)

            profile=ocean_depth-height*np.exp(-(radius-radius[-1]/5)**2/sigma**2)
            alphas = np.degrees(np.arctan(np.diff(profile) / np.diff(radius)))  # angle in degrees

            sea_surface = LineString([(0,0),(radius[-1],0)])
            sea_bottom  = LineString(zip(radius,np.abs(profile)))

            # alpha_up = -np.degrees(np.arctan((profile[1] - profile[0]) / (radius[1] - radius[0])))
            # alpha_down = -np.degrees(np.arctan((profile[-1] - profile[-2]) / (radius[1] - radius[0])))



            for iii,omega in enumerate(omegas):       # TODO loop fur zeit
                #print(T[17])

                n_t=35

                p=np.sin(np.radians(theta_0))/v_p_W
                q=np.cos(np.radians(theta_0))/v_p_W

                Pressure=np.exp(1j*omega*(p*radius+q*profile-tt[n_t]))      # TODO was machen wir mit der zeit

                x_sum = [0]
                y_sum = [0]

                distance=0
                R=1
                theta = theta_0
                int_pt = Point(0,0)

                distance_tot=[]

                for kk in range(100):


                    int_pt_0=int_pt

                    downward_ray = LineString([int_pt,(int_pt.xy[0][0]+80000*np.tan(np.radians(theta)), int_pt.xy[1][0]+80000)])

                    int_pt=downward_ray.intersection(sea_bottom)
                    if int_pt.is_empty:
                        break

                    distance+=int_pt.distance(int_pt_0)
                    x_sum = np.append(x_sum, int_pt.xy[0][0])
                    y_sum = np.append(y_sum, int_pt.xy[1][0])
                    int_pt_0 = int_pt


                    if int(np.ceil(int_pt.xy[0][0] / (radius[1] - radius[0]))) > 3890:
                        continue

                    alpha = alphas[int(np.ceil(int_pt.xy[0][0] / (radius[1] - radius[0])))]
                    # plt.plot(radius[0:2600]/1000,alphas[0:2600])
                    # plt.ylabel("angle [°]")
                    # plt.xlabel("distance [km]")
                    # plt.show()
                    theta = theta - alpha




                    uppward_ray = LineString([int_pt,(int_pt.xy[0][0]+80000*np.tan(np.radians(theta)),int_pt.xy[1][0]-80000)])
                    int_pt = uppward_ray.intersection(sea_surface)
                    if int_pt.is_empty:
                        break

                    ####x_sum,_sum brauchen wir nur fürs plotten der rays
                    x_sum = np.append(x_sum, int_pt.xy[0][0])
                    y_sum = np.append(y_sum, int_pt.xy[1][0])
                    #########################

                    ################## reflection coeff from gualtieri
                    r_1=rho_C*v_p_C*(1-2*p**2*v_s_C**2)**2*np.cos(np.radians(theta))
                    r_2=4*v_s_C**3*p**2*rho_C*np.sqrt(1-p**2*v_p_C**2)*np.sqrt(1-p**2*v_s_C**2)*np.cos(np.radians(theta))
                    r_3=rho_W*v_p_W*np.sqrt(1-p**2*v_p_C**2)
                    R*=(r_1+r_2-r_3)/(r_1+r_2+r_3)
                    if kk==0:
                        R_tot[odd, iii] = R
                    if np.isnan(R):
                        break
                    ############################################

                    distance += int_pt.distance(int_pt_0)
                    phase=distance/v_p_W
                    p=np.sin(np.radians(theta))/v_p_W
                    q=np.cos(np.radians(theta))/v_p_W

                    Pressure+=R*np.exp(1j*omega*(phase+p*radius+q*profile-tt[n_t]))
                    distance_tot=np.append(distance_tot,alpha)
                kks[odd,iii]=kk



                # points = np.array([radius / 1000, profile]).T.reshape(-1, 1, 2)
                # segments = np.concatenate([points[:-1], points[1:]], axis=1)
                #
                # fig, axs = plt.subplots(1, 1, sharex=True, sharey=True)
                # fig.set_size_inches(15, 10)
                # plt.fill_between(radius / 1000, np.ones(4000) * np.max(np.abs(profile)) + 1, np.abs(profile), linewidth=5,
                #              color="grey")
                # plt.fill_between(radius / 1000, np.abs(profile), np.zeros(4000) + 1 - 50, linewidth=5, color="blue",
                #              alpha=0.2)
                # plt.plot(x_sum / 1000, y_sum - 90, color="green")
                # plt.plot([0, radius[-1] / 1000], [0 - 60, 0 - 60], linewidth=5, color="blue")
                # norm = plt.Normalize(np.real(Pressure).min(), np.real(Pressure).max())
                # lc = LineCollection(segments, cmap='hot', norm=norm)
                # lc.set_array(np.real(Pressure))
                # lc.set_linewidth(10)
                # line = axs.add_collection(lc)
                # fig.colorbar(line, ax=axs)
                # axs.set_xlim(50, 300)
                # axs.set_ylim(ocean_depth, -100)
                # plt.xlabel("distance [km]")
                # plt.ylabel("ocean depth [m]")
                # plt.show()

                normal=np.sum(np.cos(np.abs(np.radians(alphas[0:3999])))*Pressure[0:3999]*(radius[1]-radius[0]))
                transversal=np.sum(np.sin(np.abs(np.radians(alphas[0:3999])))*Pressure[0:3999])*(radius[1]-radius[0])
                transversal_tot[odd,iii]=transversal
                normal_tot[odd,iii]=normal
                pressure_tot[odd,iii]=np.max(np.abs(Pressure))
        np.save(working_dir+"/ampl_factor",pressure_tot)
        np.save(working_dir+"/normal",normal_tot)
        np.save(working_dir+"/transversal",transversal_tot)
        np.save(working_dir+"/kk_tot",kks)
        np.save(working_dir+"/R_tot",R_tot)






#
#
#                 if iii==0:
#                     points = np.array([radius / 1000, profile]).T.reshape(-1, 1, 2)
#                     segments = np.concatenate([points[:-1], points[1:]], axis=1)
#                     fig, axs = plt.subplots(1, 1, sharex=True, sharey=True)
#                     fig.set_size_inches(15, 10)
#                     plt.fill_between(radius / 1000, np.ones(4000) * np.max(np.abs(profile)) + 1, np.abs(profile), linewidth=5,
#                                      color="grey")
#                     plt.fill_between(radius / 1000, np.abs(profile), np.zeros(4000) + 1 - 50, linewidth=5, color="blue",
#                                      alpha=0.2)
#                     plt.plot(x_sum / 1000, y_sum - 90, color="green")
#                     plt.plot([0, radius[-1] / 1000], [0 - 60, 0 - 60], linewidth=5, color="blue")
#                     norm = plt.Normalize(np.real(Pressure).min(), np.real(Pressure).max())
#                     lc = LineCollection(segments, cmap='hot', norm=norm)
#                     lc.set_array(np.real(Pressure))
#                     lc.set_linewidth(10)
#                     line = axs.add_collection(lc)
#                     fig.colorbar(line, ax=axs)
#                     axs.set_xlim(radius.min() / 1000, radius[-1] / 1000)
#                     axs.set_ylim(ocean_depth, -100)
#                     plt.xlabel("distance [km]")
#                     plt.ylabel("ocean depth [m]")
#                     plt.savefig("/home/djamel/Desktop/depth_vs_period/depth_vs_period_normal_"+str(np.round(ocean_depth))+"_"+str(np.round(height))+".png", dpi=150)
#                     plt.close()


        # plt.imshow(np.flipud(pressure_tot), extent=(T[0], T[-1],ocean_depths[0], ocean_depths[-1]), aspect='auto', cmap='gnuplot')
        # plt.xlabel("T [s]")
        # plt.ylabel("height [m]")
        # plt.colorbar()
        # plt.show()
        #
        # plt.imshow(np.abs(np.flipud(kks)), extent=(T[0], T[-1],ocean_depths[0], ocean_depths[-1]), aspect='auto', cmap='gnuplot')
        # plt.xlabel("T [s]")
        # plt.ylabel("height [m]")
        # plt.colorbar()
        # # plt.savefig(
        # #     "/home/djamel/Desktop/depth_vs_period/depth_vs_period_kk_" + str(np.round(theta_0)) + "_" + str(
        # #         np.round(height)) + ".png", dpi=150)
        # plt.show()
        # #plt.close()
        #
        # plt.imshow(np.flipud(R_tot), extent=(T[0], T[-1],ocean_depths[0], ocean_depths[-1]), aspect='auto', cmap='gnuplot')
        # plt.xlabel("T [s]")
        # plt.ylabel("height [m]")
        # plt.colorbar()
        # # plt.savefig(
        # #     "/home/djamel/Desktop/depth_vs_period/depth_vs_period_R_" + str(np.round(theta_0)) + "_" + str(
        # #         np.round(height)) + ".png", dpi=150)
        # plt.show()
        # #plt.close()
        # #
        # #
        # #
        # plt.imshow(np.abs(np.flipud(transversal_tot)), extent=(T[0], T[-1],ocean_depths[0], ocean_depths[-1]), aspect='auto', cmap='gnuplot')
        # plt.xlabel("T [s]")
        # plt.ylabel("height [m]")
        # plt.colorbar()
        # # plt.savefig(
        # #     "/home/djamel/Desktop/depth_vs_period/depth_vs_period_tranversal_" + str(np.round(theta_0)) + "_" + str(
        # #         np.round(height)) + ".png", dpi=150)
        # plt.show()
        # plt.close()
        #
        #
        # plt.figure( figsize=(10, 8))
        # plt.imshow(np.abs(np.flipud(normal_tot)), extent=(T[0], T[-1],heights[0], heights[-1]), aspect='auto', cmap='gnuplot')
        # plt.xlabel("T [s]")
        # plt.ylabel("height [m]")
        # plt.colorbar()
        # # plt.savefig(
        # #     "/home/djamel/Desktop/depth_vs_period/depth_vs_period_tranversal_" + str(np.round(theta_0)) + "_" + str(
        # #         np.round(height)) + ".png", dpi=150)
        # plt.show()
        # plt.close()

# # plt.imshow(np.abs(np.flipud(trans_over_normal)), extent=(heights[0], heights[-1],ocean_depths[0], ocean_depths[-1]), aspect='auto', cmap='gnuplot')
# # plt.xlabel("height [m]")
# # plt.ylabel("depth [m]")
# # plt.colorbar()
# # plt.show()
# #
#
#
#     # np.save("/home/djamel/PHD_projects/FOH_ICE/data_flat_sea_floor/transversal_normal_data/max_P_depth: " + str(
#     #     int(ocean_depth)) + ".npy", transversal_tot)
#     # np.save("/home/djamel/PHD_projects/FOH_ICE/data_flat_sea_floor/transversal_normal_data/transversal_depth: "+str(int(ocean_depth))+".npy",transversal_tot)
#     # np.save("/home/djamel/PHD_projects/FOH_ICE/data_flat_sea_floor/transversal_normal_data/normal_depth: "+str(int(ocean_depth))+".npy",normal_tot)
#     # np.save("/home/djamel/PHD_projects/FOH_ICE/data_flat_sea_floor/transversal_normal_data/iter_step_depth: "+str(int(ocean_depth))+".npy",iter_tot)
#
#
#
# # plt.imshow(transversal_diff_angles ,extent=(T[0],T[-1],thetas[0],thetas[-1]), aspect='auto',cmap='gnuplot')
# # plt.colorbar()
# # plt.show()
# # plt.imshow(normal_diff_angles ,extent=(T[0],T[-1],thetas[0],thetas[-1]), aspect='auto',cmap='gnuplot')
# # plt.colorbar()
# # plt.show()
# # plt.plot(T,np.sum(normal_diff_angles,axis=0),"r")
# # plt.show()
# # plt.plot(T,np.sum(transversal_diff_angles,axis=0),"b")
# # plt.show()
# # plt.subplot(211)
# # plt.plot(T,normal)
# #plt.subplot(212)
# #plt.plot(T,transversal)
# # plt.subplot(313)
# # plt.plot(T,np.abs(transversal)/np.abs(normal))
# # plt.axis([1,12,0,0.5])
#
# #plt.show()
