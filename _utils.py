# mit diesem script will ich die first arrival berchnen
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
#from scipy import signal


def phases(d_r,dr_1,num_r):


    # dr_1=7600   #distance first receiver
    # # d_r=500    # distance receivers
    # # num_r=85    # number of receivers

    vp_1=3398.0   # velocity erste layer
    vp_2=4016       # velocity zweite layer
    vp_3=5009       # velocity  dritte layer
    vp_4 = 5818  # velocity  dritte layer

    vs_1=2274   # velocity erste layer
    vs_2=2723       # velocity zweite layer
    vs_3=3393       # velocity  dritte layer
    vs_4 =3964  # velocity  dritte layer


    h1=1000          # thickness of the first and second layer
    h2=1000             # thickness of the third layer
    h3=2000  # thickness of the third layer

    #########
    ts=0# sollte eigentlich 0 sein ---> habe es so eingestellt das first arrival uebereinstimmt
    ########

    dist_r=np.zeros((num_r,1))      # distance between receivers
    for ii in range(0,num_r):
        dist_r[ii,0]=dr_1+d_r*ii


    t_dir_p=np.zeros((num_r,1))
    t_dir_s=np.zeros((num_r,1))

    for ii in range(0,num_r):
        t_dir_p[ii,0]=ts+dist_r[ii,0]/vp_1
        t_dir_s[ii,0]=ts+dist_r[ii,0]/vs_1


    # reflected phase
    t_refl_p=np.zeros((num_r,1))
    t_refl_s=np.zeros((num_r,1))

    for ii in range(0,num_r):
        t_refl_p[ii,0]=ts+2.0/vp_1*np.sqrt((dist_r[ii,0]/2.0)**2.0+h1**2.0)
        t_refl_s[ii,0]=ts+2.0/vs_1*np.sqrt((dist_r[ii,0]/2.0)**2.0+h1**2.0)



    # refrected wave at the first layer


    t_refr_p_1=np.zeros((num_r,1))
    t_refr_s_1=np.zeros((num_r,1))

    for ii in range(0,num_r):
        t_refr_p_1[ii,0]=ts+2*h1*np.sqrt((1/vp_1)**2-(1/vp_2)**2)+dist_r[ii]/vp_2
        t_refr_s_1[ii,0]=ts+2*h1*np.sqrt((1/vs_1)**2-(1/vs_2)**2)+dist_r[ii]/vs_2




    # refrected wave at the second layer

    t_refr_p_2=np.zeros((num_r,1))
    t_refr_s_2=np.zeros((num_r,1))

    for ii in range(0,num_r):
        t_refr_p_2[ii,0]=ts+2*h1*np.sqrt((1/vp_1)**2-(1/vp_2)**2)+2*h2*np.sqrt((1/vp_2)**2-(1/vp_3)**2)+dist_r[ii]/vp_3
        t_refr_s_2[ii,0]=ts+2*h1*np.sqrt((1/vs_1)**2-(1/vs_2)**2)+2*h2*np.sqrt((1/vs_2)**2-(1/vs_3)**2)+dist_r[ii]/vs_3






    # refrected wave at the third layer

    t_refr_p_3=np.zeros((num_r,1))
    t_refr_s_3=np.zeros((num_r,1))

    for ii in range(0,num_r):
        t_refr_p_3[ii,0]=ts+2*h1*np.sqrt((1/vp_1)**2-(1/vp_2)**2)+2*h2*np.sqrt((1/vp_2)**2-(1/vp_3)**2)+2*h3*np.sqrt((1/vp_3)**2-(1/vp_4)**2)+dist_r[ii]/vp_4
        t_refr_s_3[ii,0]=ts+2*h1*np.sqrt((1/vs_1)**2-(1/vs_2)**2)+2*h2*np.sqrt((1/vs_2)**2-(1/vs_3)**2)+2*h3*np.sqrt((1/vs_3)**2-(1/vs_4)**2)+dist_r[ii]/vs_4



    return dist_r,t_dir_p,t_dir_s,t_refl_p,t_refl_s,t_refr_p_1,t_refr_s_1, t_refr_p_2,t_refr_s_2,t_refr_p_3,t_refr_s_3



def max_2D_array(A):
    '''
    find maximum in a 2D array
    '''
    len_x=len(A[0,:])
    ma = np.zeros((1, len_x))
    for kk in range(0, len_x):
        ma[0, kk] = max(A[:,kk])
    maa=max(ma[0,:])
    return maa



def _theo_disp(swi,pl=1):
    'diese funktion ist nur fuer dispersion methode zu verwenden'

    green_line = mlines.Line2D([], [], color='tan', markersize=10, linewidth=2, label='Fundamental Mode')
    red_line = mlines.Line2D([], [], color='wheat', markersize=10, linewidth=2, label='Overtones')
    black_line = mlines.Line2D([], [], color='c', markersize=10, linewidth=2, label='2. overtone')
    cyan_line = mlines.Line2D([], [], color='c', markersize=10, linewidth=2, label='3. overtone')
    magenta_line = mlines.Line2D([], [], color='c', markersize=10, linewidth=2, label='4. overtone')
    yellow_line = mlines.Line2D([], [], color='c', markersize=10, linewidth=2, label='5. overtone')

    prov=[green_line, red_line, black_line, cyan_line, magenta_line,yellow_line]

    file = open("/home/djamel/PHD_projects/force_on_hill/scripting/GEOPSY_theo_disp/four_layer_model_5_mode.disp", "r")
    line = file.readlines()

    if swi=='R':
        swi=0
    elif swi=='T':
        swi=1



    b=[]
    ll=0
    for ii in range(0,len(line)):

        if line[ii][0]=='#':
            if line[ii+1][0] != '#':
                b.append(ll)
            elif line[ii+1][0] == '#':
                if line[ii+2][0] == '#':
                    continue
                elif line[ii+2][0] != '#':
                    b.append(ll)
        elif line[ii][0]!='#':
            ll=ll+1
        last=ii
    b.append(last)

    b=[x-1 for x in b[1:]]
    b[0]=0
    b=np.reshape(b,(2,len(b)/2))






    data = np.loadtxt('/home/djamel/PHD_projects/force_on_hill/scripting/GEOPSY_theo_disp/four_layer_model_5_mode.disp', skiprows=8)



    if pl==2:       # bedinotgung zum plotten
        for ii in range(0,len(b[0,:])-1):
            if ii==0:

                #plt.plot(data[b[swi,ii]+1:b[swi,ii+1]-3, 0], 1 / data[b[swi,ii]+1:b[swi,ii+1]-3, 1]/1000.0, linewidth=4,color='chocolate', alpha=0.9)
                plt.plot(data[b[swi,ii]+1:b[swi,ii+1]-3, 0], 1 / data[b[swi,ii]+1:b[swi,ii+1]-3, 1]/1000.0, linewidth=6,color='chocolate', alpha=0.9)
            if ii != 0:
                #plt.plot(data[b[swi, ii] + 1:b[swi, ii + 1] - 3, 0], 1 / data[b[swi, ii] + 1:b[swi, ii + 1] - 3, 1]/1000.0,linewidth=4, color='wheat', alpha=0.9)
                plt.plot(data[b[swi, ii] + 1:b[swi, ii + 1] - 3, 0], 1 / data[b[swi, ii] + 1:b[swi, ii + 1] - 3, 1]/1000.0,linewidth=6, color='lightgrey', alpha=0.9)
            ###############!!!!!!! -2 bei array fuer 20Hz getan weil sonst kurve udberlappt vorher war -2 nicht da!!!!!

        #plt.legend(handles=prov[0:2], loc=3, fontsize=18)      #legend
    return data,b





def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return idx


def _mode_max(DISP,freq,v,swi):
    '''hier wird die maximum values der einzelnen moden gefunden und die freq und v dieses ortes zurueckgegegen'''
    file = open("/home/djamel/GEOPSY/gpdc/four_layer_model_5_mode.disp", "r")
    line = file.readlines()
    #swi='T'
    if swi == 'R':
        swi = 0
    elif swi == 'T':
        swi = 1

    b = []
    ll = 0
    for ii in range(0, len(line)):

        if line[ii][0] == '#':

            if line[ii + 1][0] != '#':

                b.append(ll)

            elif line[ii + 1][0] == '#':

                if line[ii + 2][0] == '#':

                    continue
                elif line[ii + 2][0] != '#':
                    b.append(ll)
        elif line[ii][0] != '#':
            ll = ll + 1
        last = ii
    b.append(last)

    b = [x - 1 for x in b[1:]]
    b[0] = 0
    b = np.reshape(b, (2, len(b) / 2))

    data = np.loadtxt('/home/djamel/GEOPSY/gpdc/four_layer_model_5_mode.disp', skiprows=8)

    fr=[]
    vv=[]
    ma=0.0
    maa=[]
    int_h = 1
    int_l=3

    incr_a=[17,15,10,4]
    incr_e=[70,50,41,33]
    ll=0

    for kk in range(0,4):
        fr=[]
        vv=[]
        for jj in range(incr_a[ll],incr_e[ll],1):
            for ii in range(b[swi,kk]+1,b[swi,kk+1]-2):
                if ii==b[swi,kk]+jj:
                    ww=[data[ii, 0], 1 / data[ii, 1]]
                    fr.append(find_nearest(freq,ww[0]))
                    vv.append(find_nearest(v, ww[1]))


            for ii in range(0, len(vv)):
                prov=max_2D_array(abs(DISP[fr[ii]-int_h:fr[ii]+int_h,vv[ii]-int_l:vv[ii]+int_l]))
                if prov>ma:
                    ma=prov
        maa=np.append(maa,ma)
        ma = 0.0

        # for ii in range(0,len(vv)):
        #     plt.imshow(abs(DISP[fr[ii]-int_h:fr[ii]+int_h,vv[ii]-int_l:vv[ii]+int_l]).T,
        #                aspect='auto', extent=(freq[fr[ii]-int_l], freq[fr[ii]+int_l], v[vv[ii]+int_h], v[vv[ii]-int_h]))



        #plt.plot(data[b[swi, kk] + 1:b[swi, kk + 1] - 2, 0], 1 / data[b[swi, kk] + 1:b[swi, kk + 1] - 2, 1], linewidth=3)


        ll=ll+1


    #plt.show()

    return maa


def filter(data,N,Wn):
    b, a = signal.butter(N, Wn, 'low')
    #data_filt=np.zeros_like(data)
    # for ii in range(0,len(data[0,:,0])):
    #     for jj in range(0, len(data[0, 0, :])):
    data_filt = signal.filtfilt(b, a, data)

    #data_filt[:,ii,jj] = [[signal.filtfilt(b, a, data[:,ii,jj]) for ii in range(0,len(data[0,:,0]))] for jj in range(0,len(data[0,0,:]))]
    return data_filt

def gauss_window(data,t,shift,shift2):
    # wird verwendet um einzelne moden aus den seismgraommen zu "filtern"
    window = signal.gaussian(len(t), std=12)
    window = np.roll(window, shift2)
    window=np.roll(window, shift)
    data_win=window*data
    return data_win


def format_paper():
    ### colors and font size for paper
    plt.style.use('ggplot')
    print(plt.style.available)

    SMALL_SIZE = 16
    MEDIUM_SIZE = 25
    BIGGER_SIZE = 20
    plt.rc('font', size=BIGGER_SIZE)  # controls default text sizes
    plt.rc('axes', titlesize=BIGGER_SIZE)  # fontsize of the axes title
    plt.rc('axes', labelsize=BIGGER_SIZE)  # fontsize of the x and y labels
    plt.rc('xtick', labelsize=BIGGER_SIZE)  # fontsize of the tick labels
    plt.rc('ytick', labelsize=BIGGER_SIZE)  # fontsize of the tick labels
    plt.rc('legend', fontsize=BIGGER_SIZE-1)  # legend fontsize
    # plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title






def max_location(A):
    # find maximum location of a 2D array
    # A_ 2D array

    hor=len(A[0,:])
    ver=len(A[:,0])

    maxx=np.amax(A)
    max_loc=np.argmax(A)
    col=np.mod(max_loc, hor)
    row=max_loc/hor

    return maxx, col, row




# a = np.zeros((5,8))
# a[3,2]=6
# print(a)
# print(max_location(a))

# b=np.zeros((2,2))
# print(a)
# print(np.amax(a))           # Maximum of the flattened array
# print(np.argmax(a))
# print(np.mod(32,3))

# vp=[5.5,6,6.5,7.5]
# vs=[3.17,3.46,3.75,4.33]
# rho=[2400,2600,2750,3000]
#
# print(rho)