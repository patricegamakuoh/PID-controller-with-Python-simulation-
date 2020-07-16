import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

# animate plolspace?
animate=True # True / False

########## vehicle dynamics define ####### 
def vehicle(v,t,u,load):
    # inpulspace
    #  v    = vehicle velocity (m/s)
    #  t    = time (sec)
    #  u    = gas pedal position (-50% to 100%)
    #  load = passenger load + cargo (kg)
    Cd = 0.24    # drag coefficient
    rho = 1.225  # air density (kg/m^3)
    A = 5.0      # cross-sectional area (m^2)
    Fp = 30      # thrust parameter (N/%pedal)
    m = 500      # vehicle mass (kg)

    # calculate derivative of the velocity
    dv_dt = (1.0/(m+load)) * (Fp*u - 0.5*rho*Cd*A*v**2)
    return dv_dt
    
######## simulation #######################
tf = 300.0                    # final time for simulation
nsteps = 951                  # number of time steps
delta_t = tf/(nsteps-1)       # how long is each time step?
lspace = np.linspace(0,tf,nsteps) # linearly spaced time vector

######### simulate step test operation ###########
step = np.zeros(nsteps) # u = valve % open
step[11:] = 50.0        # step up pedal position
load = 200.0            #passenger(s) + cargo load kg

########### velocity initial condition ###########
v0 = 0.0
vs = np.zeros(nsteps)

############## PI Parameters ##################### 
ubias=0.0                                        #
Kc =1.0/1.2*0.5  #K gain constand                #
tauI=40.0                                        #                   
sum_int=0.0                                      #
reset_windup=0.0 #Anti-Reset Windup              #
##################################################

######## for storing the resullspace ################
es=np.zeros(nsteps)
ies=np.zeros(nsteps)
sps = np.zeros(nsteps)

sp = 25.0 # set point

################## simulate with ODEINT ##########
for i in range(nsteps-1):
    #change the set point to any value within the range 
    # if i ==200:
    #     sp=0
    # if i ==30000:
    #     sp=17
    # if i ==150:
    #     sp=20
    if i ==400:
        sp=15

############ PI Loop controller #################                                                                   
    sps[i+1] = sp                               #          
    error =sp-v0                                # 
    es[i+1]=error                               #
    sum_int=sum_int+error+delta_t               #  
    # u=ubias+Kc+error+Kc/tauI*sum_int          #   
    reset_windup = reset_windup+Kc/tauI*error   #       
    u=Kc+error+reset_windup                     #  
#################################################

 ######## clip inpulspace to -50% to 100%  ##########
    if u >= 100.0:
        u = 100.0
        sum_int = sum_int-error * delta_t
    if u <= -50.0:
        u = -50.0
        sum_int = sum_int-error * delta_t
    ies[i+1]=sum_int
    step[i+1]=u
    v = odeint(vehicle,v0,[0,delta_t],args=(u,load))
    v0 = v[-1]   # take the last value
    vs[i+1] = v0 # store the velocity for plotting

#################### plot resullspace ###################
plt.figure()
plt.subplot(2,2,1)
plt.plot(lspace,vs,'b-',linewidth=3)
plt.plot(lspace,sps,'k--',linewidth=2)
plt.ylabel('Velocity (m/s)')
plt.legend(['Velocity','Set Point'],loc='best')
plt.subplot(2,2,2)
plt.plot(lspace,step,'r--',linewidth=3)
plt.ylabel('Gas Pedal')    
plt.legend(['Gas Pedal (%)'],loc='best')
plt.subplot(2,2,3)
plt.plot(lspace,es,'b--',linewidth=3)
plt.legend(['Error (SP-PV)'])
plt.xlabel('Time (sec)')   
plt.subplot(2,2,4)
plt.plot(lspace,ies,'k--',linewidth=3)
plt.legend(['Intergral of Error'])
plt.xlabel('Time (sec)')   
plt.show()

