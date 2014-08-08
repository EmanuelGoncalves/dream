__author__ = 'daniel'

from matplotlib.pyplot import plot, legend, show

linreg_zmin = [0.08, 0.1, 0.15, 0.2]
linreg_res = [0.17815, 0.18179, 0.17934, 0.17471]

ridge_zmin = [0.08, 0.1, 0.2]
ridge_res = [0.17824, 0.1821, 0.1742]

par_zmin = [0.08, 0.1, 0.2]
par_res = [0.17992, 0.18398, 0.1766]

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

legend(['par', 'ridge', 'bayrdg', 'linreg', 'svm'])

show()
