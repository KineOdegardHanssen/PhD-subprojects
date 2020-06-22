import matplotlib.pyplot as plt

d     = 1.5
sigma = 1

# The alpha parameter makes the color transparent
circle1 = plt.Circle((0, 0), sigma, color='r', alpha=0.7)
circle2 = plt.Circle((0, d), sigma, color='r', alpha=0.7)
circle3 = plt.Circle((0, 2*d), sigma, color='r', alpha=0.7)
circle4 = plt.Circle((d, 0), sigma, color='r', alpha=0.7)
circle5 = plt.Circle((d, d), sigma, color='r', alpha=0.7)
circle6 = plt.Circle((d, 2*d), sigma, color='r', alpha=0.7)
circle7 = plt.Circle((2*d, 0), sigma, color='r', alpha=0.7)
circle8 = plt.Circle((2*d, d), sigma, color='r', alpha=0.7)
circle9 = plt.Circle((2*d, 2*d), sigma, color='r', alpha=0.7)

fig, ax = plt.subplots()

plt.xlim(-d/2.,2.5*d)
plt.ylim(-d/2.,2.5*d)

plt.grid(linestyle='--')

ax.set_aspect(1)

ax.add_artist(circle1)
ax.add_artist(circle2)
ax.add_artist(circle3)
ax.add_artist(circle4)
ax.add_artist(circle5)
ax.add_artist(circle6)
ax.add_artist(circle7)
ax.add_artist(circle8)
ax.add_artist(circle9)

plt.title('Brush, d=%.2f, sigma=%i' % (d,sigma), fontsize=8)

plt.savefig("brush_d%.2f_sigma%i.png" % (d,sigma), bbox_inches='tight')

plt.show()