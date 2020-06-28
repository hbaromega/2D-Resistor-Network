import numpy as np
from scipy.sparse.linalg import spsolve
from scipy.sparse import csr_matrix as cs
import matplotlib.pyplot as plt
import time



def main():
  start = time.time()
  #Lx,Ly = map(int, input('Enter Lx & Ly\n').split())
  #Ly = int(input('Enter Ly\n'))
  #V0 = float(input('Enter Circuit bias, V0\n'))
  #R = float(input('Enter resistance, R\n'))
  #Lx=100
  #Ly=100
  
  V0=1.0; R=1.0
  print("V0=",V0, "R=",R)
  L_list = [10, 20, 50, 100, 200] 
  #L_list = [10, 20] 
  for Lx in L_list:
    print('Parameter Lx=',Lx)  
    f0=open('Reff_vs_Ly_fixed_Lx{}_rect.dat'.format(Lx), 'w')
    f1=open('Rratio_vs_Ly_fixed_Lx{}_rect.dat'.format(Lx), 'w')
    print('f0=',f0)
    for Ly in range(2,201):  # Breaks down for Ly=1, hence starting from Ly=2
         print("Ly=",Ly)
         Reff=Effective_resistance(V0,R,Lx,Ly)
         print("Effective resistance =", Reff, " Anaytically expected",2*R*Lx/float(Ly))
         print(Ly,Reff,1/Ly, file=f0)
         Rratio=Reff*Ly/float(R*Lx)
         Rratio='{0:.3f}'.format(Rratio)
         print(Ly, Rratio, file=f1)

    f0.close() 
    f1.close() 


  end = time.time()
  hr, rem = divmod(end-start, 3600)
  min, sec = divmod(rem, 60)
  print("Execution time = {:0>2} h {:0>2} m {:05.2f} s".format(int(hr),int(min),sec))


  import plot_rect_Reff_vs_Ly_fixed_Lx 
  import plot_rect_Rratio_vs_Ly_fixed_Lx 


def Effective_resistance(V0,R,Lx,Ly):

    

    R_eff_list=[]
    Lx_list=[]
    R_anaLytical_list=[]
  
    Rp=R
    #G=np.zeros((Lx*Ly,Lx*Ly), float)####conductance matrix
    I=np.zeros(Lx*Ly, float)
    row_G=[]
    col_G=[]
    data_G=[]

    # Mapping: (i,j) --> jLx + i; e.g. (0,0)->0; (0,1)->Lx
				# Original mapping (j-1)Lx+i 

    ## Left bottom corner point (i=0,j=0) # V-indices: 0,0->0; 1,0->1; 1,0->Lx
    #					 # [Natural index sys: 1,1->1; 2,1->2; 1,2->Lx+1]            
    row_G.append(0)
    row_G.append(0)
    row_G.append(0)
    col_G.append(0)
    col_G.append(1)
    col_G.append(Lx)  # relates to V_{12} [V_{01} in Python indexing sys] 
    data_G.append( 1/R+1/(R+Rp) + 1/(R+Rp) ) 
    data_G.append(-1/(R+Rp))
    data_G.append(-1/(R+Rp))
    I[0]=V0/R

    ## Left top corner point (i=0, j=1 to Ly-2) 
    row_G.append((Ly-1)*Lx)
    row_G.append((Ly-1)*Lx)
    row_G.append((Ly-1)*Lx)
    col_G.append((Ly-1)*Lx)
    col_G.append((Ly-1)*Lx+1)
    col_G.append((Ly-2)*Lx)
    data_G.append(1/R+1/(R+Rp) + 1/(R+Rp))
    data_G.append(-1/(R+Rp))
    data_G.append(-1/(R+Rp))

    I[(Ly-1)*Lx]=V0/R

    ## Right top corner point (i=Lx-1, j=Ly-1)
    row_G.append(Lx*Ly-1)
    row_G.append(Lx*Ly-1)
    row_G.append(Lx*Ly-1)
    col_G.append(Lx*Ly-1)
    col_G.append(Lx*Ly-2)
    col_G.append(Lx*Ly-Lx-1)
    data_G.append(1/R+1/(R+Rp)+1/(R+Rp))
    data_G.append(-1/(R+Rp))
    data_G.append(-1/(R+Rp))

    # Right bottom corner point
    row_G.append(Lx-1)
    row_G.append(Lx-1)
    row_G.append(Lx-1)
    col_G.append(Lx-1)
    col_G.append(Lx-2)
    col_G.append(2*Lx-1)
    data_G.append(1/R+1/(R+Rp)+1/(R+Rp))
    data_G.append(-1/(R+Rp))
    data_G.append(-1/(R+Rp))


    # NOTE: In Python's range(1,X), loop runs from 1 to X-1, not 1 to X

    for j in range(1,Ly-1): # j runs from 1 to Ly-2 
        # Left non-corner edge point (i=1)
                       # V-indices:        0,j->jLx;       1,j->jLx+1;     0,j-1->(j-1)Lx,   0,j+1->(j+1)Lx
                       # [Natural indices: 1,j->(j-1)Lx+1; 2,j->(j-1)Lx+2; 1,j-1->(j-2)Lx+1; 1,j+1->jLx+1]
        row_G.append(j*Lx)
        row_G.append(j*Lx)
        row_G.append(j*Lx)
        row_G.append(j*Lx)
        col_G.append(j*Lx)
        col_G.append(j*Lx+1)
        col_G.append((j-1)*Lx)
        col_G.append((j+1)*Lx)   # Earlier code had a wrong sign!
        data_G.append(1/R+ 1/(R+Rp)+1/(R+Rp)+1/(R+Rp))
        data_G.append(-1/(R+Rp))
        data_G.append(-1/(R+Rp))
        data_G.append(-1/(R+Rp))
        I[j*Lx]=1/R

        ## Right non corner edge point
        row_G.append(j*Lx+Lx-1)
        row_G.append(j*Lx+Lx-1)
        row_G.append(j*Lx+Lx-1)
        row_G.append(j*Lx+Lx-1)
        col_G.append(j*Lx+Lx-1)
        col_G.append(j*Lx+Lx-2)
        col_G.append((j-1)*Lx+Lx-1)
        col_G.append((j+1)*Lx+Lx-1)
        data_G.append(1/R + 1/(R+Rp)+1/(R+Rp)+1/(R+Rp))
        data_G.append(-1/(R+Rp))
        data_G.append(-1/(R+Rp))
        data_G.append(-1/(R+Rp))

    for i in range(1,Lx-1):
        # Bottom non-corner edge point
        row_G.append(i)
        row_G.append(i)
        row_G.append(i)
        row_G.append(i)
        col_G.append(i)
        col_G.append(i-1)
        col_G.append(i+1)
        col_G.append(Lx+i)
        data_G.append(1/(R+Rp)+ 1/(R+Rp) +1/(R+Rp))
        data_G.append(-1/(R+Rp))
        data_G.append(-1/(R+Rp))
        data_G.append(-1/(R+Rp))

        ## Top non-corner edge point
        row_G.append((Ly-1)*Lx+i)
        row_G.append((Ly-1)*Lx+i)
        row_G.append((Ly-1)*Lx+i)
        row_G.append((Ly-1)*Lx+i)
        col_G.append((Ly-1)*Lx+i) 
        col_G.append((Ly-1)*Lx+i-1)
        col_G.append((Ly-1)*Lx+i+1)
        col_G.append((Ly-2)*Lx+i) 
        data_G.append(1/(R+Rp) +1/(R+Rp) +1/(R+Rp))
        data_G.append(-1/(R+Rp))
        data_G.append(-1/(R+Rp))
        data_G.append(-1/(R+Rp))
 
        ## Non border inner points
    for i in range(1,Lx-1):
      for j in range(1,Ly-1):
        row_G.append(j*Lx+i)
        row_G.append(j*Lx+i)
        row_G.append(j*Lx+i)
        row_G.append(j*Lx+i)
        row_G.append(j*Lx+i)
        col_G.append(j*Lx+i)
        col_G.append(j*Lx+i-1)
        col_G.append(j*Lx+i+1)
        col_G.append((j-1)*Lx+i)
        col_G.append((j+1)*Lx+i)
        data_G.append(1/(R+Rp) +1/(R+Rp)+1/(R+Rp) +1/(R+Rp))
        data_G.append(-1/(R+Rp))
        data_G.append(-1/(R+Rp))
        data_G.append(-1/(R+Rp))
        data_G.append(-1/(R+Rp))

    G=cs((data_G, (row_G,col_G)), shape=(Lx*Ly,Lx*Ly), dtype=float)
                   ###converting the sparse matrix to compressed sparse column format
    V=spsolve(G, I)###solving the matrix

    I_out=0.0
    #####for calculating net current
    for j in range(0,Ly):
        I_out += V[j*Lx+Lx-1]/R

    
    Reff=V0/I_out # Effective resistance 
    print("V0=",V0,"Lx=",Lx,"Ly=",Ly)
    return Reff

  


if __name__ == '__main__':
    main()

print("Done!")


