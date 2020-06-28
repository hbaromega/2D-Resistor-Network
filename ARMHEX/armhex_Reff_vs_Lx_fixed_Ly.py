import numpy as np
from scipy.sparse.linalg import spsolve
from scipy.sparse import csr_matrix 
import matplotlib.pyplot as plt
import time



def main():
  start = time.time()
    
  V0=1.0; R=1.0
  print("V0=",V0, "R=",R)

  test_mode=1
  if test_mode==1:
     Lx,Ly = map(int, input('Enter Lx & Ly\n').split())
     #Ly = int(input('Enter Ly\n'))
     #V0 = float(input('Enter Circuit bias, V0\n'))
     #R = float(input('Enter resistance, R\n'))
     Reff=Effective_resistance(V0,R,Lx,Ly)
     print('Reff for Lx={} and Ly={}:'.format(Lx,Ly),Reff)
     exit() 



  L_list = [10, 20, 50, 100, 200] 
  #L_list = [10, 20] 
  for Ly in L_list:
    print('Parameter Ly=',Ly)  
    f0=open('Reff_vs_Lx_fixed_Ly{}_armhex.dat'.format(Ly), 'w')
    f1=open('Rratio_vs_Lx_fixed_Ly{}_armhex.dat'.format(Ly), 'w')
    print('f0=',f0)
    for Lx in range(2,201,2):  # Breaks down for Ly=1, hence starting from Ly=2
         print("Lx=",Lx)
         Reff=Effective_resistance(V0,R,Lx,Ly)
         print("Effective resistance =", Reff, " Anaytically expected",2*R*Lx/float(Ly))
         print(Lx, Reff, file=f0)
         print(Lx, Reff*Ly/float(R*Lx), file=f1)
    f0.close() 
    f1.close() 


  end = time.time()
  hr, rem = divmod(end-start, 3600)
  min, sec = divmod(rem, 60)
  print("Execution time = {:0>2} h {:0>2} m {:05.2f} s".format(int(hr),int(min),sec))


  import plot_armhex_Reff_vs_Lx_fixed_Ly 
  import plot_armhex_Rratio_vs_Lx_fixed_Ly 


def Effective_resistance(V0,R,Lx,Ly):



  #NL=2*Lx*Ly+2*Lx*(Ly+1)
  #NL=Lx*(2*Ly+1)+1
  NL=Lx*(2*Ly+1)

  print('Lx=',Lx,' Ly=',Ly,' NL=',NL)

  I=np.zeros(NL, float)
  row_G=[]
  col_G=[]
  data_G=[]
  nz=0

  # Mapping:
  # (i,j,k) -> i + jLx + k LxLy [Original mapping: i + (j-1)Lx + k LxLy  ]



  # Left border points (i=0, j=0 to Ly-1, k=0) : M-type 
  for j in range(Ly):
    row_G.append(j*Lx)
    row_G.append(j*Lx)
    row_G.append(j*Lx)
    col_G.append(j*Lx)
    col_G.append((j+1)*Lx+Lx*Ly)
    col_G.append(j*Lx+Lx*Ly)
    data_G.append(2/R)
    data_G.append(-1/(2*R))
    data_G.append(-1/(2*R))
    I[j*Lx]=V0/R
    nz+=3

  ## Right border points (i=Lx-1, j=1 to Ly-1, k=0) : M-type
  for j in range(Ly):
        row_G.append(Lx-1+j*Lx)
        row_G.append(Lx-1+j*Lx)
        row_G.append(Lx-1+j*Lx)
        col_G.append(Lx-1+j*Lx)
        col_G.append(Lx-1+(j+1)*Lx+Lx*Ly)
        col_G.append(Lx-1+j*Lx+Lx*Ly)
        data_G.append(2/R)
        data_G.append(-1/(2*R))
        data_G.append(-1/(2*R))
        nz+=3



  ## Top border points (i=0 to Lx-1, j=Ly-1, k=1) : S-type
  for i in range(Lx):
    if (i%2==0):
        row_G.append(i+Ly*Lx+Lx*Ly)
        row_G.append(i+Ly*Lx+Lx*Ly)
        row_G.append(i+Ly*Lx+Lx*Ly)
        col_G.append(i+Ly*Lx+Lx*Ly)
        col_G.append(i+1+Ly*Lx+Lx*Ly)
        col_G.append(i+(Ly-1)*Lx)
        data_G.append(1/R)
        data_G.append(-1/(2*R))
        data_G.append(-1/(2*R))
        nz+=3
    elif (i%2!=0):
        row_G.append(i+Ly*Lx+Lx*Ly)
        row_G.append(i+Ly*Lx+Lx*Ly)
        row_G.append(i+Ly*Lx+Lx*Ly)
        col_G.append(i+Ly*Lx+Lx*Ly)
        col_G.append(i-1+Ly*Lx+Lx*Ly)
        col_G.append(i+(Ly-1)*Lx)
        data_G.append(1/R)
        data_G.append(-1/(2*R))
        data_G.append(-1/(2*R))
        nz+=3
        
        
  ## Bottom border points (i=0 to Lx-1, j=0, k=1) : S-type
  for i in range(Lx):
    if (i%2==0):
        row_G.append(i+Lx*Ly)
        row_G.append(i+Lx*Ly)
        row_G.append(i+Lx*Ly)
        col_G.append(i+Lx*Ly)
        col_G.append(i+1+Lx*Ly)
        col_G.append(i)
        data_G.append(1/R)
        data_G.append(-1/(2*R))
        data_G.append(-1/(2*R))
        nz+=3
    elif (i%2!=0):
        row_G.append(i+Lx*Ly)
        row_G.append(i+Lx*Ly)
        row_G.append(i+Lx*Ly)
        col_G.append(i+Lx*Ly)
        col_G.append(i-1+Lx*Ly)
        col_G.append(i)
        data_G.append(1/R)
        data_G.append(-1/(2*R))
        data_G.append(-1/(2*R))
        nz+=3
    
        
  
         
  
  ### For inner points
  # M-type (k=0) Mapping: (i,j) -> jLx + i [Original mapping:  (i,j) -> (j-1)Lx + i]
  for j in range(Ly): # Original range: 1 to Ly
    for i in range(1,Lx-1):  # Original range: 2 to Lx-1
        if (i%2==0):
            row_G.append(i+j*Lx)
            row_G.append(i+j*Lx)
            row_G.append(i+j*Lx)
            row_G.append(i+j*Lx)
            col_G.append(i+j*Lx)
            col_G.append(i-1+j*Lx)
            col_G.append(i+j*Lx+Lx*Ly)
            col_G.append(i+(j+1)*Lx+Lx*Ly)
            data_G.append(3/(2*R))
            data_G.append(-1/(2*R))
            data_G.append(-1/(2*R))
            data_G.append(-1/(2*R))
            nz+=4
        elif (i%2!=0):
            row_G.append(i+j*Lx)
            row_G.append(i+j*Lx)
            row_G.append(i+j*Lx)
            row_G.append(i+j*Lx)
            col_G.append(i+j*Lx)
            col_G.append(i+1+j*Lx)
            col_G.append(i+j*Lx+Lx*Ly)
            col_G.append(i+(j+1)*Lx+Lx*Ly)
            data_G.append(3/(2*R))
            data_G.append(-1/(2*R))
            data_G.append(-1/(2*R))
            data_G.append(-1/(2*R))
            nz+=4

  # S-type (k=1) 
   #         Mapping: (i,j) -> i+jLx + LxLy  
   #         [Original mapping:  (i,j) -> i + (j-1)Lx + LxLy]
  for j in range(1,Ly): # Original range: 2 to Ly
    for i in range(Lx+1): # Orginal range: 1 to Lx
        if (i%2==0):
            row_G.append(i+j*Lx+Lx*Ly)
            row_G.append(i+j*Lx+Lx*Ly)
            row_G.append(i+j*Lx+Lx*Ly)
            row_G.append(i+j*Lx+Lx*Ly)
            col_G.append(i+j*Lx+Lx*Ly)
            col_G.append(i+j*Lx)
            col_G.append(i+(j-1)*Lx)
            col_G.append(i+1+j*Lx+Lx*Ly)
            data_G.append(3/(2*R))
            data_G.append(-1/(2*R))
            data_G.append(-1/(2*R))
            data_G.append(-1/(2*R))
            nz+=4
        elif (i%2!=0):
            row_G.append(i+j*Lx+Lx*Ly)
            row_G.append(i+j*Lx+Lx*Ly)
            row_G.append(i+j*Lx+Lx*Ly)
            row_G.append(i+j*Lx+Lx*Ly)
            col_G.append(i+j*Lx+Lx*Ly)
            col_G.append(i+j*Lx)
            col_G.append(i+(j-1)*Lx)
            col_G.append(i-1+j*Lx+Lx*Ly)
            data_G.append(3/(2*R))
            data_G.append(-1/(2*R))
            data_G.append(-1/(2*R))
            data_G.append(-1/(2*R))
            nz+=4


  print("NL=",NL)

  ''' commenting out below
  values=np.zeros(nz,float)
  columns=np.zeros(nz, int)
  rows=np.zeros(nz, int)
  rowIndex=np.zeros(NL+1,int)
  rowIndex[0]=0
  for i in range (nz):
       values[i]=data_G[i] # Stores values of nz elements
       columns[i]=col_G[i]
       rows[i]=row_G[i]
       rw=row_G[i] # The row-coordinate
       rowIndex[rw+1]=rowIndex[rw+1]+1 # counts nonzero elements in given row
       #rowIndex[rw+1] += 1
  for i in range (1,NL+1):
       rowIndex[i]=rowIndex[i]+rowIndex[i-1]
       #print("i,rowIndex[i]=",i,rowIndex[i])

  #print("rowIndex[Lx], N_nz =",rowIndex[Lx], nz)
  commented out above '''


  print('Sizes of row_G, col_G, and data_G:')
  print(len(row_G),len(col_G),len(data_G))
  sz=len(data_G)
  print('NL=',NL)



  #G=csr_matrix((data_G, (row_G,col_G)), shape=(Lx*Ly,Lx*Ly), dtype=float)# converting the sparse matrix to compressed sparse column format
  #G = csr_matrix((values, columns, rowIndex), shape=(NL, NL))
  G=csr_matrix((data_G,(row_G,col_G)), shape=(NL,NL),dtype=float)
  #print (G)
  #G = csr_matrix((data_G, col_G, rowIndex), shape=(NL, NL))
  #print (G)
  #print(G)
  #G = csr_matrix((values, (columns, rows)), shape=(NL, NL))


  
  V=spsolve(G, I) # Solving the matrix eq. G V = I  [Ax=b format]
  #print(V)
  #print(I)





  # Calculating net current
  I_out=0.0
  for j in range(Ly):
        #I_out+=V[2*Lx-1+j*2*Lx]/R
        I_out+=V[j*Lx+Lx-1]/R
        #print (I_out)


  Reff=V0/I_out # Effective resistance
  print("V0=",V0,"Lx=",Lx,"Ly=",Ly)#, "p=",p_latt)
  print("Effective resistance =", Reff)
  return Reff
    


if __name__ == '__main__':
    main()

print("Done!")
