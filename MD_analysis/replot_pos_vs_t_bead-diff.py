import matplotlib.pyplot as plt                     # To plot

damp     = 100
psigma   = 1
spacing  = 10
confignr = 1000

endlocation     = '/home/kine/Projects_PhD/P2_PolymerMD/Planar_brush/Diffusion_bead_near_grid/Spacing%i/damp%i_diffseedLgv/Brush/Sigma_bead_' % (spacing,damp)+str(psigma) + '/'

infilename_xoft = endlocation+'xoft_config%i.txt' % confignr # (times_single[i],xthis))
infilename_yoft = endlocation+'yoft_config%i.txt' % confignr
infilename_zoft = endlocation+'zoft_config%i.txt' % confignr

infile_xoft = open(infilename_xoft,'r')
infile_yoft = open(infilename_yoft,'r')
infile_zoft = open(infilename_zoft,'r')

xpos = []
ypos = []
zpos = []
step = []

xlines = infile_xoft.readlines()
ylines = infile_yoft.readlines()
zlines = infile_zoft.readlines()

for i in range(len(xlines)):
    xwords = xlines[i].split()
    ywords = ylines[i].split()
    zwords = zlines[i].split()
    if len(xwords)>0: # I think there is a trailing line
        step.append(int(xwords[0]))
        xpos.append(float(xwords[1]))
        ypos.append(float(ywords[1]))
        zpos.append(float(zwords[1]))

# Could plot xy, xz, yz too

plt.figure(figsize=(6,5))
plt.plot(step, xpos)
#plt.plot(step, xpos, ',')
plt.xlabel(r't')
plt.ylabel(r'x')
plt.title('x(t)')
plt.tight_layout()
plt.show()

plt.figure(figsize=(6,5))
plt.plot(step, ypos)
#plt.plot(step, ypos, ',')
plt.xlabel(r't')
plt.ylabel(r'y')
plt.title('y(t)')
plt.tight_layout()
plt.show()

plt.figure(figsize=(6,5))
plt.plot(step, zpos)
#plt.plot(step, zpos, ',')
plt.xlabel(r't')
plt.ylabel(r'z')
plt.title('z(t)')
plt.tight_layout()
plt.show()
