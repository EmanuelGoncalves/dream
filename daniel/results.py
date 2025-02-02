__author__ = 'daniel'

from matplotlib.pyplot import plot, legend, show


#Sub challenge 1

linreg_zmin = [0.08, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5]
linreg_res = [0.17815, 0.18179, 0.17934, 0.17471, 0.16295, 0.12616, 0.04444]

ridge_zmin = [0.08, 0.1, 0.2, 0.3, 0.4, 0.5]
ridge_res = [0.17824, 0.1821, 0.1742, 0.16444, 0.12717, 0.057409]

par_zmin = [0.08, 0.1, 0.2, 0.3, 0.4, 0.5]
par_res = [0.17992, 0.18398, 0.1766, 0.16847, 0.13991, 0.049336]

svm_zmin = [0.08, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5]
svm_res = [0.17705, 0.18047, 0.18108, 0.17869, 0.17507, 0.16436, 0.16531]

#other
bayrdg_zmin = [0.08, 0.1, 0.2]
bayrdg_res = [0.17824, 0.18209, 0.17419]
ard = {0.3: 0.1171}
logreg = {0.4: 0.12662}
tree = {0.2: 0.053013}


plot(linreg_zmin, linreg_res, 'o-',
     ridge_zmin, ridge_res, 'o-',
     par_zmin, par_res, 'o-',
     svm_zmin, svm_res, 'o-'
     )

legend(['linreg', 'ridge', 'par', 'svm'], 'lower left')

show()


#Sub challenge 2
