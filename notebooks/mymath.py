
# coding: utf-8

# In[4]:

n=5
import numpy as np

x = np.arange(n)
y = np.zeros(n)
for i in range(n):
    y[i] = x[i]*x[i] 
z = x*x
print(x)
print(y)
print(z)
    


# In[5]:

import matplotlib.pylab as plt
plt.plot(x,y)
plt.plot(x,z)
plt.show()


# In[ ]:



