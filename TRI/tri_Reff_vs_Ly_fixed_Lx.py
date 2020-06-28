import numpy as np
import matplotlib.pyplot as plt
import time


from tri_find_Reff_routine import *


def main():
  start = time.time()
  V0=1.0; R=1.0
  

  test_mode=1
  if test_mode==0:
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
  for Lx in L_list:
    print('Parameter Lx=',Lx)  
    f0=open('Reff_vs_Ly_fixed_Lx{}_tri.dat'.format(Lx), 'w')
    f1=open('Rratio_vs_Ly_fixed_Lx{}_tri.dat'.format(Lx), 'w')
    print('f0=',f0)
    for Ly in range(2,201):  # Breaks down for Ly=1, hence starting from Ly=2
         print("Ly=",Ly) 
         Reff=Effective_resistance(V0,R,Lx,Ly)
         print("Effective resistance =", Reff, " Anaytically expected",2*R*Lx/float(Ly))
         print(Ly, Reff, 1/Ly, file=f0)
         Rratio=Reff*Ly/float(R*Lx)
         Rratio='{0:.3f}'.format(Rratio)
         print(Ly, Rratio, file=f1)
    f0.close() 
    f1.close() 


  end = time.time()
  hr, rem = divmod(end-start, 3600)
  min, sec = divmod(rem, 60)
  print("Execution time = {:0>2} h {:0>2} m {:05.2f} s".format(int(hr),int(min),sec))


  import plot_tri_Reff_vs_Ly_fixed_Lx 
  import plot_tri_Rratio_vs_Ly_fixed_Lx 

  


  
if __name__ == '__main__':
    main()

print("Done!")


