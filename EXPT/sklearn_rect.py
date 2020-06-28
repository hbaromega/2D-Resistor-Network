# imports
import numpy as np
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_squared_error, r2_score

# Read data-set
#np.random.seed(0)
#x = np.random.rand(100, 1)
#y = 2 + 3 * x + np.random.rand(100, 1)

#datafile='set_hexagonal.dat'
datafile='set_rectangular.dat'
data = np.loadtxt(datafile)
y = data[:,1]
x = data[:,0]

# Scikit needs reshaping as if it's 2D (m x n) structure (n=1 here, column vector): 
y = y.reshape(-1,1)
x = x.reshape(-1,1)
print('y=',y)
print('x=',x)

# sckit-learn implementation

# Model initialization
regression_model = LinearRegression()
# Fit the data(train the model)
regression_model.fit(x, y)
# Predict
y_pred = regression_model.predict(x)

# model evaluation
rmse = mean_squared_error(y, y_pred)
r2 = r2_score(y, y_pred)

# printing values
print('Slope:' ,regression_model.coef_)
print('Intercept:', regression_model.intercept_)
print('Root mean squared error: ', rmse)
print('R2 score: ', r2)

# plotting values

# data points
plt.scatter(x, y, s=10)
plt.xlabel('x')
plt.ylabel('y')

# pred values
plt.plot(x, y_pred, color='r')
plt.show()


