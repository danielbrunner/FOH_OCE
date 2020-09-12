"""here i want to caluclate the pressure at the bottom of a hill, to check for resonances in dependence of the take off angle etc; hier speichere ich auch
die maximum von den spektren"""


# hier werde ich eine karte des meeres laden, eine source hintun und berechnen wie weit das acoustic field
# reflektiert werden kann in abhängigkeit von theta....

from scipy.interpolate import interp1d
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap, BoundaryNorm
from shapely.geometry import LineString, Point
import matplotlib.colors as colors
from _utils import format_paper
import os
format_paper()



########### staff that i need for plotting
#plt.ion() ## Note this

x_label = [0, 200, 400, 600, 800, 1000]
lat = ['34°W', '32°W', '30°W', '28°W', '26°W', '24°W']
y_label = [0, 200, 400, 600, 800]
lon = ['36°N', '34°N', '32°N', '30°N']
##################################



###################################### map
elevation = Dataset("/home/djamel/PHD_projects/FOH_ICE/sea_bottom_data/GEBCO_2020_30_Jul_2020_f4a81978d75e"
                 "/gebco_2020_n36.0_s27.0_w-34.0_e-22.0.nc", "r").variables['elevation']



dist_x=415
dist_y=415
dist_in_y=np.linspace(0,elevation.shape[0],elevation.shape[0])*dist_y
dist_in_x=np.linspace(0,elevation.shape[1],elevation.shape[1])*dist_x

elevation=np.array(elevation)

####################################################



theta_0=4.0                                 # angle that wave propagates (degree)
radius=np.linspace(0,350000,500)
x_0=dist_in_x[900]
y_0=dist_in_y[1050]

print(x_0)
print(y_0)

angle=0.0            # hier in radians
v_p_W=1450        # velocity of acoustic waves in water [m/s]
v_p_C=3540
v_s_C=v_p_C/np.sqrt(3)
rho_W=1000      # kg/m3a
rho_C=2500      # kg/m3

########################## pressure at the sea surface
T=np.linspace(1,12,80)
omegas=2*np.pi/T
tt=np.linspace(0,20,100)


################################################################
################################################################




x=x_0+np.cos(angle)*radius
y=y_0+np.sin(angle)*radius
n_x=x/dist_x
n_y=y/dist_y


###############################################
######################### creating folder

#define the name of the directory to be created
path = "/home/djamel/Documents/PAPER_SCHREIBEN/FOH_OCE/plots_realistic_reflection_on_seafloor/acustic_field_on_the_sea_floor/x_0: "+str(int(x_0))+\
       " y_0: "+str(int(y_0))+" theta_0: "+str(theta_0)

try:
    os.mkdir(path)
except OSError:
    print ("Creation of the directory %s failed" % path)
else:
    print ("Successfully created the directory %s " % path)


####################################


profile=[]
for nn in range(len(n_x)):
    profile=np.append(profile,elevation[int(n_y[nn]),int(n_x[nn])])
profile_interpol = interp1d(radius, profile, kind='cubic')
radius=np.linspace(0,350000,5000)       # interpolation
profile=profile_interpol(radius)

#############################################



alphas=np.degrees(np.arctan(np.diff(profile)/(radius[1]-radius[0])))                          # angle in degrees

# plt.subplot(211)
# plt.plot(radius,profile)
# plt.subplot(212)
# plt.plot(radius[1::],alphas)
# plt.show()

sea_surface = LineString([(0,0),(1500000,0)])
#sea_bottom = LineString(zip(radius,np.abs(profile)))
sea_bottom = LineString([(0,4000),(1500000,5000)])
alpha=np.degrees(np.arctan(-1000/1500000))
#
# plt.subplot(211)
# plt.plot(radius,profile,"o")
# plt.subplot(212)
# plt.plot(radius[0:-1],alphas,"o")
# plt.show()


######################
#


maximum_of_spectrum=[]
for iii,omega in enumerate(omegas):       # TODO loop fur zeit

    # figure = plt.gcf()  # get current figure
    # figure.set_size_inches(15, 10)
    n_t=20


    p=np.sin(np.radians(theta_0))/v_p_W
    q=np.cos(np.radians(theta_0))/v_p_W

    Pressure=np.cos(omega*(p*radius+q*np.abs(profile)-tt[n_t]))      # TODO was machen wir mit der zeit

    x_sum = [0]
    y_sum = [0]



    distance=0
    R=1
    theta = theta_0
    int_pt = Point(0, 0)


    for kk in range(140):
        int_pt_0=int_pt

        downward_ray = LineString([int_pt,(int_pt.xy[0][0]+50000*np.tan(np.radians(theta)), int_pt.xy[1][0]+50000)])


        int_pt=downward_ray.intersection(sea_bottom)

        distance+=int_pt.distance(int_pt_0)
        x_sum = np.append(x_sum, int_pt.xy[0][0])
        y_sum = np.append(y_sum, int_pt.xy[1][0])
        int_pt_0 = int_pt

        #alpha=alphas[int(int_pt.xy[0][0]/(radius[1]-radius[0]))]
        theta=theta-alpha



        uppward_ray = LineString([int_pt,(int_pt.xy[0][0]+30000*np.tan(np.radians(theta)),int_pt.xy[1][0]-30000)])
        int_pt = uppward_ray.intersection(sea_surface)

        ####x_sum,_sum brauchen wir nur fürs plotten der rays
        x_sum = np.append(x_sum, int_pt.xy[0][0])
        y_sum = np.append(y_sum, int_pt.xy[1][0])
        #########################



        ################## reflection coeff from gualtieri
        r_1=rho_C*v_p_C*(1-2*p**2*v_s_C**2)**2*np.cos(np.radians(theta))
        r_2=4*v_s_C**3*p**2*rho_C*np.sqrt(1-p**2*v_p_C**2)*np.sqrt(1-p**2*v_s_C**2)*np.cos(np.radians(theta))
        r_3=rho_W*v_p_W*np.sqrt(1-p**2*v_p_C**2)
        R*=(r_1+r_2-r_3)/(r_1+r_2+r_3)

        ############################################

        distance += int_pt.distance(int_pt_0)
        phase=distance*omega/v_p_W
        p=np.sin(np.radians(theta))/v_p_W
        q=np.cos(np.radians(theta))/v_p_W
        Pressure+=R*np.cos(phase+omega*(p*radius+q*np.abs(profile)-tt[n_t]))

    maximum_of_spectrum=np.append(maximum_of_spectrum,np.max(Pressure))


    plt.rcParams['axes.facecolor']='white'
    points = np.array([radius / 1000, np.abs(profile)]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    #
    fig, axs = plt.subplots(1, 1, sharex=True, sharey=True)
    fig.set_size_inches(15, 10)
    # plt.fill_between(radius / 1000, np.ones(5000) * np.max(np.abs(profile)) + 1, np.abs(profile), linewidth=5,color="grey")
    # plt.fill_between(radius / 1000, np.abs(profile), np.zeros(5000) + 1-50 , linewidth=5, color="blue", alpha=0.2)
    plt.plot(x_sum / 1000, y_sum-40, color="green")
    plt.plot([0, radius[-1] / 1000], [0-50, 0-50], linewidth=5, color="blue")
    norm = plt.Normalize(Pressure.min(), Pressure.max())
    lc = LineCollection(segments, cmap='hot', norm=norm)
    lc.set_array(Pressure)
    lc.set_linewidth(6)
    line = axs.add_collection(lc)
    fig.colorbar(line, ax=axs)
    axs.set_xlim(radius.min() / 1000, radius.max() / 1000)
    axs.set_ylim(np.abs(profile).max(), -100)
    plt.xlabel("distance [km]")
    plt.ylabel("ocean depth [m]")
    plt.show()


    # stri=path+"/ocean_bottom_profile(flat part)_loc: x_0:"+str(int(x_0))+" y_0:_"+str(int(y_0))+"_period: "\
    #      +str(np.round(T[iii],2))+" theta_0: "+str(theta_0)+".png"
    # plt.savefig(stri)

    stri_max="/home/djamel/PHD_projects/FOH_ICE/maximum_of_the_pressure_amplification_factor/"+"/maximum_of_ocean_bottom_profile(flat part)_loc: x_0:"+str(int(x_0))\
         +" theta_0: "+str(np.round(theta_0,2))
    np.save(stri_max,maximum_of_spectrum)
# #
#aa=np.load("/home/djamel/PHD_projects/FOH_ICE/maximum_of_the_pressure_amplification_factor/maximum_of_ocean_bottom_profile(flat part)_loc: x_0:373629 theta_0: 19.5.npy")
# bb=np.load("/home/djamel/PHD_projects/FOH_ICE/maximum_of_the_pressure_amplification_factor/maximum_of_ocean_bottom_profile(flat part)_loc: x_0:373629 theta_0: 2.95.npy")
#plt.plot(aa,"g")
# plt.plot(bb,"m")
#plt.show()


plt.plot(maximum_of_spectrum)
plt.show()

# plt.subplot(211)
# plt.plot(radius/1000,profile,linewidth=3,color="grey")
# plt.axis([0,80,-5000,-3500])
# plt.ylabel("ocean depth [m]")


# plt.subplot(212)
# plt.plot(x_sum[1::2]/1000,theta_plot,linewidth=3,color="k")
# plt.ylabel("take off angle [°]")
# plt.axis([0,40,-10,10])


#plt.show()

# ax1=plt.subplot(111)
# ax1.imshow(elevation, extent=(0,dist_in_x[-1]/1000,0,dist_in_y[-1]/1000), aspect='auto',cmap='seismic')
# #plt.colorbar()
#
#
# ax1.plot(x/1000,np.abs(dist_in_y[-1]-y)/1000,linewidth=3,color="k")
# ax1.plot(x_0/1000,np.abs(dist_in_y[-1]-y_0)/1000,"yo", markersize=8)
#
# plt.gca().invert_yaxis()
# #plt.rcParams["axes.grid"] = False
# ax1.set_xticks(x_label)
# ax1.set_xticklabels(lat, minor=False, rotation=45)
# ax1.set_yticks(y_label)
# ax1.set_yticklabels(lon, minor=False, rotation=45)
# #
#
# plt.show()





# TODO theta verschieben
# TODO wenn man H verschiebt dann ändert sich OCE??
# TODO wenn man freq verschiebt .....??
# TODO wenn man strahl verschiebt?? also anderes hill oder -----übereiander plotten für verschiedene hills;
#  plot mit verschiedenen stahls und OCE plots übereinander


# TODO ocean bottom map das zeig wo die amplifaction stärker ist
# TODO sum over many thetas to see the resulting amplification on each side of the hill
