import numpy as np
import matplotlib.pyplot as plt
import time


from tri_find_Reff_routine import *


def main():
  start = time.time()
  V0=1.0; R=1.0
  

  test_mode=0
  if test_mode==1:
     Lx,Ly = map(int, input('Enter Lx & Ly\n').split())
     #Ly = int(input('Enter Ly\n'))
     #V0 = float(input('Enter Circuit bias, V0\n'))
     #R = float(input('Enter resistance, R\n'))
     Reff=Effective_resistance(V0,R,Lx,Ly)
     print('Reff for Lx={} and Ly={}:'.format(Lx,Ly),Reff)
     exit() 

  print("V0=",V0, "R=",R)
  L_list = [10, 20, 50, 100, 200] 
  #L_list = [10, 20] 
  for Ly in L_list:
    print('Parameter Ly=',Ly)  
    f0=open('Reff_vs_Lx_fixed_Ly{}_tri.dat'.format(Ly), 'w')
    f1=open('Rratio_vs_Lx_fixed_Ly{}_tri.dat'.format(Ly), 'w')
    print('f0=',f0)
    for Lx in range(2,201):  # Breaks down for Ly=1, hence starting from Ly=2
         print("Lx=",Lx) 
         Reff=Effective_resistance(V0,R,Lx,Ly)
         print("Effective resistance =", Reff, " Anaytically expected",2*R*Lx/float(Ly))
         print(Lx, Reff, file=f0)
         Rratio=Reff*Ly/float(R*Lx)
         Rratio='{0:.3f}'.format(Rratio)
         print(Lx, Rratio, file=f1)
    f0.close() 
    f1.close() 


  end = time.time()
  hr, rem = divmod(end-start, 3600)
  min, sec = divmod(rem, 60)
  print("Execution time = {:0>2} h {:0>2} m {:05.2f} s".format(int(hr),int(min),sec))


  import plot_tri_Reff_vs_Lx_fixed_Ly 
  import plot_tri_Rratio_vs_Lx_fixed_Ly 

  


  
if __name__ == '__main__':
    main()

print("Done!")


