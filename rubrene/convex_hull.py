# coding: utf-8
rub = io.read("rubrene.cif")
get_ipython().magic(u'cd rubrene/')
rub = io.read("rubrene.cif")
view(rub)
frag1 = io.read("dimers/dimer_4_0_d-7.19A.xyz")
frag1.edit
frag1.edit()
frag1
from scipy.spatial import ConvexHull
hull = ConvexHull(frag1.get_positions())
hull.vertices
hull.neighbors
hull.vertices
pos = frag1.get_positions()
pos[hull.vertices]
pos[1]
hull_points = pos[hull.vertices]
com_frag1 = frag1.get_center_of_mass()
com_frag1
zeroed_hull = hull_points - com_frag1
zeroed_hull
hull2 = ConvexHull(frag1.get_positions(), qhull_options="FA")
hull2
hull2.equations
hull2.max_bound
zeroed_hull
zeroed_2d = zeroed_hull[:,[0,1]]
zeroed_2d
fig = scipy.spatial.convex_hull_plot_2d(zeroed_2d)
import scipy
fig = scipy.spatial.convex_hull_plot_2d(zeroed_2d)
fig = scipy.spatial.convex_hull_plot_2d(zeroed_2d)
zeroed_2d
hull
hull.points
pos2d = pos[:,[0,1]]
hull2d = ConvexHull(pos2d)
frag1.get_center_of_mass()
com = frag1.get_center_of_mass()
frag1.center()
frag1.get_center_of_mass()
view(frag1)
frag1.get_center_of_mass()
hull2d = ConvexHull(frag1.get_positions()[:,[0,1]])
view(frag1)
fig = scipy.spatial.convex_hull_plot_2d(hull2d)
fig.show()
fig.savefig("test.pdf")
hull2d = ConvexHull(frag1.get_positions()[:,[0,1]]*5)
fig = scipy.spatial.convex_hull_plot_2d(hull2d)
fig.savefig("test2.pdf")
