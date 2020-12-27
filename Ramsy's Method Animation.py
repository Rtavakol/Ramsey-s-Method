import numpy as np
from scipy import integrate
from scipy import special

from matplotlib import pyplot as plt
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import proj3d
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.colors import cnames
from matplotlib import animation

N_trajectories = 2
phase = 0.0 # show initial phase of B1 field
phase_of_2nd_B1_pulse = 0
field_freq = 200.0
Larmor_freq_H = 200.0
Nuclei_name = 'H'
max_range = 4000
cyc = 2.0
cyc1 = 2.0; cyc2 = 0
time_cyc1 = cyc1 * 1.0 / Larmor_freq_H
time= cyc * 1.0 / Larmor_freq_H
time_final = time + 0.005
dt = time_final / max_range
def nuclear_name(name):
    global gamma
    if name == 'H':
        gamma = 267513000.0 #rad/s/T
    if name == 'F':
        gamma = 251662000.0
    return gamma

nuclear_name('H')

fac = 1 # for rotating field we have a factor of 0.5 but for linearly polarized we have a factor of 1.
B_1 = fac * np.pi/2.0/time_cyc1/267513000.0
#dB_1 = (fac * np.pi/2.0/time_cyc1/267513000.0)*0.05
B_1_initial = fac * np.pi/2.0/time_cyc1/267513000.0

B_0 = 2.0 * np.pi * Larmor_freq_H/267513000.0

class Arrow3D(FancyArrowPatch):
    def __init__(self, xs, ys, zs, *args, **kwargs):
        FancyArrowPatch.__init__(self, (0,0), (0,0), *args, **kwargs)
        self._verts3d = xs, ys, zs

    def draw(self, renderer):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, renderer.M)
        self.set_positions((xs[0],ys[0]),(xs[1],ys[1]))
        FancyArrowPatch.draw(self, renderer)

def lorentz_deriv(y,t0):
    global time, field_freq, B_1, B_0, gamma, dt, phase, cyc1, phase_of_2nd_B1_pulse

    w = gamma * B_1
    w_l = gamma * B_0
    w1 = 2.0 * np.pi * field_freq
    """Compute the time-derivative of a Lorentz system."""
    if y[3] > time:
        w = 0.0
        if (y[3]-time)< dt/6000:
            print (y, np.arctan2(y[0],y[1])*180/np.pi)

    if y[3] < (cyc1 /field_freq):
        phase_of_B1_pulse = 0.0

    if y[3] > (cyc1 /field_freq):
        phase_of_B1_pulse = (phase_of_2nd_B1_pulse * np.pi/180.0) - (w1 * cyc1 / field_freq)

        # with liniealy polarized magnetic field with sin wt
    #return [w_l * y[1], w * y[2] *np.sin(w1 * y[3]) - w_l * y[0], - w * y[1] *np.sin(w1 * y[3]), 1.0]

    # with liniealy polarized magnetic field with cos wt
    #return [w_l * y[1] - w * y[2]*np.cos(w1*y[3])*np.sin(phase), w * y[2] *np.cos(w1 * y[3])*np.cos(phase) - w_l * y[0],
    #      y[0]* w * np.cos(w1*y[3])*np.sin(phase) - w * y[1] *np.cos(w1 * y[3])*np.cos(phase), 1.0]


            # Right hand rotating (which is on resonant)
    #return [w_l * y[1] + w * y[2] * np.sin(w1 * y[3] + phase_of_B1_pulse ),
     #       w * y[2]  * np.cos(w1 * y[3] + phase_of_B1_pulse) - w_l * y[0],
      #    - y[0] * w * np.sin(w1 * y[3] + phase_of_B1_pulse) - w * y[1] * np.cos(w1 * y[3] + phase_of_B1_pulse) , 1.0]

    # left hand rotating (which is off resonant)
    #return [w_l * y[1] - w * y[2] * np.sin(w1 * y[3]),
    #        w * y[2] * np.cos(w1 * y[3]) - w_l * y[0],
    #        + y[0] * w * np.sin(w1 * y[3]) - w * y[1] * np.cos(w1 * y[3]), 1.0]

    # Right hand rotating (which is on resonant) with sin wt x^ + cos wt y^
    return [w_l * y[1] - w * y[2] * np.cos(w1 * y[3]),
           w * y[2]  * np.sin(w1 * y[3]) - w_l * y[0],
          + y[0] * w * np.cos(w1 * y[3]) - w * y[1] * np.sin(w1 * y[3]) , 1.0]

def plot_Mz_vs_B1(A,B,C):
    global phase_of_2nd_B1_pulse, cyc1, cyc2
    #B1_form = 'Linearly Polarized Sin wt, ' + str(cyc) ; B1_name_to_put_on_plot = r'Linearly Polarized Sin($\omega_{\rm L}$t)'
    #B1_form =  'Linearly Polarized Cos wt, ' + str(cyc) ; B1_name_to_put_on_plot = r'Linearly Polarized Cos($\omega_{\rm L}$t)'
    B1_form =  'Circularly Polarized B1, ' + str(cyc1) + 'cyc(x) + ,' + str(cyc2) + 'cyc(x+' + str(phase_of_2nd_B1_pulse) + ')' + str(cyc)
    B1_name_to_put_on_plot = 'Circulrary Polarized, ' + str(cyc1) +' cyc(x)' + str(cyc2) + 'cyc(x+' + str(phase_of_2nd_B1_pulse) + ')'
    file_name = 'MzH and MzF vs B1, ' + B1_form + ', '
    fig = plt.figure()
    plt.xticks([(0.75 + w * 0.05) for w in range(11)])
    plt.yticks([(-0.8 + v * 0.1 ) for v in range(17)])
    plt.xlabel(r'Normalized B$_{1}$', fontsize=20)
    plt.ylabel(r'M$_{\rm z}$', fontsize=20, rotation=90)
    plt.plot(A,B, linestyle ='-' , label = r'M$_{\rm zH}$')
    plt.plot(A,C, linestyle = '--' , label = r'M$_{\rm zF}$')
    plt.grid(b=True, which='both', color='0.65', linestyle='-')
    plt.legend(bbox_to_anchor=(0.97, 0.97), loc=1, borderaxespad=0.)
    plt.figtext(0.5, 0.92, B1_name_to_put_on_plot, fontsize=15, color='k', ha='center')
    #plt.savefig(filename = file_name + '.png', format='png', dpi=500)

    text_to_save = np.transpose(np.asarray([A,B,C]))
    #np.savetxt(file_name + '.txt', text_to_save)



# Choose random starting points, uniformly distributed from -15 to 15
#np.random.seed(1)
x0 = [[0.0, 0.0, 1.0, 0.0 ],[1.0,0.0,0.0,0.0]] # -15 + 30 * np.random.random((N_trajectories, 3))

t = np.linspace(0, time_final , max_range) # it was (0, 0.04, 2000)

M_z_H = []; M_z_F = []; x_axis_B1 = []
for i in range(0,1,1):
    B_1 =  B_1_initial #+ (i-5) * dB_1
    x_axis_B1.append(B_1/B_1_initial)
    nuclear_name('H')
    a = integrate.odeint(lorentz_deriv, x0[0], t)
    M_z_H.append(a[-1][2])
    print(gamma*B_1/np.pi/2.0)
    nuclear_name('F')
    b = integrate.odeint(lorentz_deriv, x0[1], t)
    #print(gamma*B_1/np.pi/2.0)
    M_z_F.append(b[-1][2])

x_t = np.asarray([a,b]) # this may show false animation be careful because it is after for command which is changing B1 value
#plot_Mz_vs_B1(x_axis_B1, M_z_H, M_z_F)





a_time = []; a_x = []; a_y = []; a_z = []; b_x = []; b_y = []; b_z = []
for i in range(0,len(a),1):
    if a[i][3] < (cyc)/Larmor_freq_H :
        a_x.append(a[i][0]); a_y.append(a[i][1]) ; a_z.append(a[i][2])
        b_x.append(b[i][0]); b_y.append(b[i][1]) ; b_z.append(b[i][2])



# Set up figure & 3D axis for animation

fig = plt.figure()
ax = fig.add_axes([0, 0, 1, 1], projection='3d')
ax.set_aspect("auto")

#ax.plot(a_x, a_y, a_z, linestyle = '--' ,linewidth=3, c='k')
#ax.plot(b_x, b_y, b_z,linewidth=3, c='r')

#ax.plot([0, 0], [0, 0], [0, b_z[-1]],linewidth=3, c='r')
subtitle = 'Perfect pi/2 Rotation'

#subtitle = 'Dressing Parameters'  + '\n' r'Y$_{n}$=' + str(dress_y_n) + '   '  + r'X$_{n}$=' + \
 #          str(dress_x_n) + ', ' + r'Y$_{3}$=' + str(dress_y_3) + '   '  + r'X$_{3}$=' + str(dress_x_3)  + '\n' 'f$_{d}$=' + str(int(dress_freq)) + '  ' + 'f$_{L}$=' + \
  #         str(int(Larmor_freq)) + '   ' + 'f$_{mod}$=' + str(int(mod_larmor_fre))

fig.suptitle(subtitle, fontsize=14, fontweight='bold')
plt.figtext(0.26, 0.96, r'$\bullet$' + '  ' + 'Proton', fontsize=20, color='b', ha ='right')
plt.figtext(0.185, 0.92, r'$\bullet$' + '  ' + 'Fluorine', fontsize=20, color='r', ha ='right')
#ax.text(1.1, 0, 0, 'X', (0,1,0),fontsize=20,color='blue')
#ax.text(0, 1.1, 0, 'Y', (0,1,0),fontsize=20,color='green')
#ax.text(0, 0, 1.1, 'Z', (0,1,0),fontsize=20,color='red')
a1 = Arrow3D([0,0],[0,0],[0,1], mutation_scale=20, lw=1, arrowstyle="-|>", color="r") # to plot arrow
a2 = Arrow3D([0,1],[0,0],[0,0], mutation_scale=20, lw=1, arrowstyle="-|>", color="b") # to plot arrow
a3 = Arrow3D([0,0],[0,1],[0,0], mutation_scale=20, lw=1, arrowstyle="-|>", color="g") # to plot arrow

a4 = Arrow3D([0,0],[0,0],[0,b_z[-1]], mutation_scale=20, lw=1, arrowstyle="-|>", color="r") # to plot arrow
#a4 = Arrow3D([0,x_t[0,0]],[0,0],[0,0], mutation_scale=20, lw=5, arrowstyle="-|>", color="k")
ax.add_artist(a1)
ax.add_artist(a2)
ax.add_artist(a3)
ax.add_artist(a4)
xc, yc, zc = [], [], []
for i in range(0,100,1):
    xc.append(np.cos(i*np.pi/50))
    yc.append(np.sin(i*np.pi/50))
    zc.append(0.0)
ax.plot(xc, yc, zc,linestyle='--', linewidth=1 , c='k')
#ax.plot(zc, yc, xc,linestyle='--', linewidth=1 , c='k')
#ax.plot(xc, zc, yc,linestyle='--', linewidth=1 , c='k')
ax.plot([1,-1], [0,0], [0,0],linestyle='--', linewidth=1 , c='b')
ax.plot([0,0], [1,-1], [0,0],linestyle='--', linewidth=1 , c='g')
ax.plot([0,0], [0,0], [1,-1],linestyle='--', linewidth=1 , c='r')
u = np.linspace(0, 2 * np.pi, 100)
v = np.linspace(0, np.pi, 100)

x = np.outer(np.cos(u), np.sin(v))
y = np.outer(np.sin(u), np.sin(v))
z = np.outer(np.ones(np.size(u)), np.cos(v))
ax.plot_surface(x, y, z, rstride=5, cstride=5, linewidth=0 , color='w', alpha = 0.5)
#ax.axis('off')



# choose a different color for each trajectory
colors = plt.cm.jet(np.linspace(0, 1, N_trajectories))

# set up lines and points
lines = sum([ax.plot([], [], [], '-', c=c) for c in colors], [])
pts = sum([ax.plot([], [], [], 'o', c=c)for c in colors], [])

# prepare the axes limits
ax.set_xlim((-1.1, 1.1))
ax.set_ylim((-1.1, 1.1))
ax.set_zlim((-1.1, 1.1))
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
# set point-of-view: specified by (altitude degrees, azimuth degrees)
ax.view_init(20, 110)

# initialization function: plot the background of each frame
def init():
    for line, pt in zip(lines, pts):
        line.set_data([], [])
        line.set_3d_properties([])

        pt.set_data([], [])
        pt.set_3d_properties([])
    return lines + pts

# animation function.  This will be called sequentially with the frame number
def animate(i):
    global field_freq, cyc, dt
    # we'll step two time-steps per frame.  This leads to nice results.
    i = (50 * i) % x_t.shape[1]

    for line, pt, xi in zip(lines, pts, x_t):
        x, y, z, tt = xi[:i].T
        line.set_data(x, y)
        line.set_3d_properties(z)



        pt.set_data(x[-1:], y[-1:])
        pt.set_3d_properties(z[-1:])

    num = float(i) * dt * 1000.0
    num1 = '{:04.3f}'.format(num)
    lab = str(num1) + 'msec'
    if num / 1000.0 > cyc * (1.0 / field_freq):
        lab = lab + '\n' r'B$_{1}$ is OFF'
        time_color = 'black'
    else:
        lab = lab + '\n' r'B$_{1}$ is ON'
        time_color = 'red'
    fig.suptitle('Time = ' + lab, fontsize=14, fontweight='bold', color = time_color)


    ax.view_init(30, -0.18 * i) # -1.435 fpr 2000 data point and fL = 200 Hz
    fig.canvas.draw()
    return lines + pts

# instantiate the animator.
anim = animation.FuncAnimation(fig, animate, init_func=init, frames=2000, interval=30, blit=True) # it was 800

# Save as mp4. This requires mplayer or ffmpeg to be installed
#filename = 'Perfect pi/2 Rotation' + '.mp4'
#anim.save(filename, fps=15, extra_args=['-vcodec', 'libx264'])

plt.show()
