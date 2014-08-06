import numpy as np
import matplotlib.pyplot as plt
from sklearn.linear_model import ElasticNet, enet_path
from sklearn.metrics import r2_score

n_samples, n_features = 50, 200
X = np.random.randn(n_samples, n_features)
coef = 3 * np.random.randn(n_features)
inds = np.arange(n_features)
np.random.shuffle(inds)
coef[inds[10:]] = 0  # sparsify coef
y = np.dot(X, coef)

# add noise
y += 0.01 * np.random.normal((n_samples,))

# Split data in train set and test set
n_samples = X.shape[0]
X_train, y_train = X[:n_samples / 2], y[:n_samples / 2]
X_test, y_test = X[n_samples / 2:], y[n_samples / 2:]

en = ElasticNet(alpha=0.1, l1_ratio=0.7)

y_pred_enet = en.fit(X_train, y_train).predict(X_test)
r2_score_enet = r2_score(y_test, y_pred_enet)

print(en)
print("r^2 on test data : %f" % r2_score_enet)

plt.plot(en.coef_, label='Elastic net coefficients')
plt.plot(coef, '--', label='original coefficients')
plt.legend(loc='best')
plt.show()