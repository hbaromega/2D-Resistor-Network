import numpy as np
from scipy.sparse.linalg import spsolve
from scipy.sparse import csr_matrix 
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

  test_mode=1
  if test_mode==1:
     Lx,Ly = map(int, input('Enter Lx & Ly\n').split())
     #Ly = int(input('Enter Ly\n'))
     #V0 = float(input('Enter Circuit bias, V0\n'))
     #R = float(input('Enter resistance, R\n'))
     Reff=Effective_resistance(V0,R,Lx,Ly)
     print('Reff for Lx={} and Ly={}:'.format(Lx,Ly),Reff)
     exit() 




  #L_list = [10, 20, 50, 100, 200] 
  L_list = [2, 3, 4, 5, 6, 7]
  for Lx in L_list:
    print('Parameter Lx=',Lx)  
    f0=open('Reff_vs_Ly_fixed_Lx{}_zighex.dat'.format(Lx), 'w')
    f1=open('Rratio_vs_Ly_fixed_Lx{}_zighex.dat'.format(Lx), 'w')
    print('f0=',f0)
    for Ly in range(2,11):  # Breaks down for Ly=1, hence starting from Ly=2
         print("Ly=",Ly)
         Reff=Effective_resistance(V0,R,Lx,Ly)
         print("Effective resistance =", Reff, " Anaytically expected",2*R*Lx/float(Ly))
         print(Ly, Reff, 1/Ly, file=f0)
         print(Ly, Reff*Ly/float(R*Lx), file=f1)
    f0.close() 
    f1.close() 


  end = time.time()
  hr, rem = divmod(end-start, 3600)
  min, sec = divmod(rem, 60)
  print("Execution time = {:0>2} h {:0>2} m {:05.2f} s".format(int(hr),int(min),sec))


  import plot_zighex_Reff_vs_Ly_fixed_Lx 
  import plot_zighex_Rratio_vs_Ly_fixed_Lx 


def Effective_resistance(V0,R,Lx,Ly):



  #NL=2*Lx*Ly+2*Lx*(Ly+1)
  #NL=Lx*(2*Ly+1)+1
  NL=2*Lx*Ly
  print('V0=',V0,' R=',R)
  print('Lx=',Lx,' Ly=',Ly,' NL=',NL)

  I=np.zeros(NL, float)
  row_G=[]
  col_G=[]
  data_G=[]
  nz=0
 
  # Conductance factors
  oc=1.0/R  
  tc=0.5/R #1.0/(2.0*R)
  fc=0.25/R #1.0/(4.0*R)

  # Mapping:
  # (i,j,k) -> i + jLx + k LxLy                    [in natural indexing : i + (j-1)Lx + k LxLy]

  LL=Lx*Ly

  ## Bottom left corner point (type O): i = 0, j = 0; k = 1
  # P = 0 + 0*Lx + 1*LL = LL
  row_G.append(LL)
  row_G.append(LL)
  row_G.append(LL)
  col_G.append(LL)
  col_G.append(1+LL) # (1,0,1) pt
  col_G.append(0) # (0,0,0) pt
  data_G.append(fc+tc)
  data_G.append(-fc)
  data_G.append(-tc)
  nz+=3 
 
  ## Bottom right corner point (type O): i = Lx-1, j = 0; k = 1
  # P = Lx-1 + 0*Lx + 1*LL = Lx-1+LL
  row_G.append(Lx-1+LL)
  row_G.append(Lx-1+LL)
  row_G.append(Lx-1+LL)
  col_G.append(Lx-1+LL)
  col_G.append(Lx-2+LL) # (Lx-2,0,1) pt
  col_G.append(Lx-1) # (Lx-1,0,0) pt
  data_G.append(fc+tc)
  data_G.append(-fc)
  data_G.append(-tc)
  nz+=3 


  ## Top left corner point (type I): i = 0, j = Ly-1; k = 0.
  # P = 0 + (Ly-1)*Lx + 0*LL = (Ly-1)*Lx
  row_G.append((Ly-1)*Lx)
  row_G.append((Ly-1)*Lx)
  row_G.append((Ly-1)*Lx)
  col_G.append((Ly-1)*Lx)
  col_G.append(1+(Ly-1)*Lx) # (1,Ly-1,0) pt
  #col_G.append(Lx-2+LL) # (Lx-2,0,1) pt
  col_G.append((Ly-1)*Lx+LL) # (0,Ly-1,1) pt
  data_G.append(oc+fc+tc)
  data_G.append(-fc)
  data_G.append(-tc) 
  I[(Ly-1)*Lx]=V0/R
  nz+=3


  ## Top right corner point (type I): i = Lx-1, j = Ly-1; k = 0
  # P = Lx-1 + (Ly-1)*Lx + 0*LL =  Lx-1+(Ly-1)*Lx
  row_G.append( Lx-1+(Ly-1)*Lx )
  row_G.append( Lx-1+(Ly-1)*Lx )
  row_G.append( Lx-1+(Ly-1)*Lx )
  col_G.append( Lx-1+(Ly-1)*Lx )
  col_G.append( Lx-2+(Ly-1)*Lx ) # (Lx-2,Ly-1,0) pt 
  col_G.append( Lx-1+(Ly-1)*Lx+LL ) # (Lx-1,Ly-1,1) pt 
  data_G.append(oc+fc+tc)
  data_G.append(-fc)
  data_G.append(-tc) 
  nz+=3
 

  ## Bottom non-corner border points (type O): i = 1 to Lx − 2, j = 0; k = 1.
  # P = i + 0*Lx + LL = i+LL
  for i in range(1,Lx-1):
      row_G.append(i+LL)
      row_G.append(i+LL)
      row_G.append(i+LL)
      row_G.append(i+LL)
      col_G.append(i+LL)
      col_G.append(i-1+LL) # (i-1,0,1) pt
      col_G.append(i+1+LL) # (i+1,0,1) pt
      col_G.append(i) # (i,0,0) pt
      data_G.append(2*fc+tc)
      data_G.append(-fc)
      data_G.append(-fc)
      data_G.append(-tc)
      nz+=4

  ## Top non-corner border points (type I): i = 1 to Lx − 2, j = Ly-1; k = 0
  # P = i + (Ly-1)*Lx + 0*LL = i+(Ly-1)*Lx
  for i in range(1,Lx-1):
      row_G.append(i+(Ly-1)*Lx)
      row_G.append(i+(Ly-1)*Lx)
      row_G.append(i+(Ly-1)*Lx)
      row_G.append(i+(Ly-1)*Lx)
      col_G.append(i+(Ly-1)*Lx)
      col_G.append(i-1+(Ly-1)*Lx) # (i-1,Ly-1,0) pt
      col_G.append(i+1+(Ly-1)*Lx) # (i+1,Ly-1,0) pt
      col_G.append(i+(Ly-1)*Lx+LL) # (i,Ly-1,1) pt
      data_G.append(2*fc+tc)
      data_G.append(-fc)
      data_G.append(-fc)
      data_G.append(-tc)
      nz+=4

  
  ## Left non-corner border points (type I): i = 0, j = 0 to Ly-2; k = 0
  # P = 0 + j*Lx + 0 = j*Lx
  for j in range(Ly-1):
      if (j%2==0): 
          row_G.append(j*Lx)
          row_G.append(j*Lx)
          row_G.append(j*Lx)
          col_G.append(j*Lx)
          col_G.append(j*Lx+LL) # (0,j,1) pt
          col_G.append((j+1)*Lx+LL) # (0,j+1,1) pt
          data_G.append(oc+2*tc)
          data_G.append(-tc)
          data_G.append(-tc)
          nz+=3
      else:
          row_G.append(j*Lx)
          row_G.append(j*Lx)
          row_G.append(j*Lx)
          row_G.append(j*Lx)
          col_G.append(j*Lx)
          col_G.append(j*Lx+LL) # (0,j,1) pt
          col_G.append((j+1)*Lx+LL) # (0,j+1,1) pt
          col_G.append(1+(j+1)*Lx+LL) # (1,j+1,1) pt
          data_G.append(oc+3*tc)
          data_G.append(-tc)
          data_G.append(-tc)
          data_G.append(-tc)
          nz+=4
      I[j*Lx]=V0/R
     
  ## Left non-corner border points (type O): i = 0, j = 1 to Ly-1; k = 1  
  # P = 0 + j*Lx + 1*LL = j*Lx+LL
  for j in range(1,Ly):
      if (j%2==0): 
          row_G.append(j*Lx+LL)
          row_G.append(j*Lx+LL)
          row_G.append(j*Lx+LL)
          col_G.append(j*Lx+LL) 
          col_G.append(j*Lx) # (0,j,0) pt
          col_G.append((j-1)*Lx) # (0,j-1,0) pt
          data_G.append(oc)
          data_G.append(-tc)
          data_G.append(-tc)
          nz+=3
      else:
          row_G.append(j*Lx+LL)
          row_G.append(j*Lx+LL)
          row_G.append(j*Lx+LL)
          row_G.append(j*Lx+LL)
          col_G.append(j*Lx+LL) 
          col_G.append(j*Lx) # (0,j,0) pt
          col_G.append((j-1)*Lx) # (0,j-1,0) pt
          col_G.append(1+(j-1)*Lx) # (1,j-1,0) pt
          data_G.append(3*tc)
          data_G.append(-tc)
          data_G.append(-tc)
          data_G.append(-tc)
          nz+=4
      
  ## Right non-corner border points (type I): i = Lx-1, j = 0 to Ly-2; k = 0 
  # P = Lx-1 + j*Lx + 0 = Lx-1
  for j in range(Ly-1):
      if (j%2==0): 
          row_G.append(Lx-1+j*Lx)
          row_G.append(Lx-1+j*Lx)
          row_G.append(Lx-1+j*Lx)
          row_G.append(Lx-1+j*Lx)
          col_G.append(Lx-1+j*Lx)
          col_G.append(Lx-1+j*Lx+LL) # (Lx-1,j,1) pt
          col_G.append(Lx-1+(j+1)*Lx+LL) # (Lx-1,j+1,1) pt
          col_G.append(Lx-2+(j+1)*Lx+LL) # (Lx-2,j+1,1) pt
          data_G.append(oc+3*tc)
          data_G.append(-tc)
          data_G.append(-tc)
          data_G.append(-tc)
          nz+=4
      else:
          row_G.append(Lx-1+j*Lx)
          row_G.append(Lx-1+j*Lx)
          row_G.append(Lx-1+j*Lx)
          col_G.append(Lx-1+j*Lx)
          col_G.append(Lx-1+j*Lx+LL) # (Lx-1,j,1) pt
          col_G.append(Lx-1+(j+1)*Lx+LL) # (Lx-1,j+1,1) pt
          data_G.append(2*oc)
          data_G.append(-tc)
          data_G.append(-tc)
          nz+=3

  ## Right non-corner border points (type O): i = Lx-1, j = 1 to Ly-1; k = 1
  # P = Lx-1 + j*Lx + LL = Lx-1+j*Lx+LL
  for j in range(1,Ly):
          row_G.append(Lx-1+j*Lx+LL)
          row_G.append(Lx-1+j*Lx+LL)
          row_G.append(Lx-1+j*Lx+LL)
          col_G.append(Lx-1+j*Lx+LL)
          col_G.append(Lx-1+j*Lx) # (Lx-1,j,0) pt
          col_G.append(Lx-1+(j-1)*Lx) # (Lx-1,j-1,0) pt
          data_G.append(oc)
          data_G.append(-tc)
          data_G.append(-tc)
          nz+=3
     
         
  ## Inner points (type I): i = 1 to Lx − 2, j = 0 to Ly − 2; k = 0.
  # P = i+j*Lx
  for j in range(Ly-1): 
    for i in range(1,Lx-1):  
            row_G.append(i+j*Lx)
            row_G.append(i+j*Lx)
            row_G.append(i+j*Lx)
            row_G.append(i+j*Lx)
            col_G.append(i+j*Lx)
            col_G.append(i+j*Lx+LL) # (i,j,1) pt
            col_G.append(i+(j+1)*Lx+LL) # (i,j+1,1) pt
            col_G.append(i+1+j*Lx+LL) # (i+1,j,1) pt 
            data_G.append(3*tc)
            data_G.append(-tc)
            data_G.append(-tc)
            data_G.append(-tc)
            nz+=4
  
  ## Inner points (type O): i = 1 to Lx − 2, j = 1 to Ly − 1; k = 1
  # P = i+j*Lx+LL
  for j in range(1,Ly): 
    for i in range(1,Lx-1):  
            row_G.append(i+j*Lx+LL)
            row_G.append(i+j*Lx+LL)
            row_G.append(i+j*Lx+LL)
            row_G.append(i+j*Lx+LL)
            col_G.append(i+j*Lx+LL)
            col_G.append(i+j*Lx) # (i,j,0) pt
            col_G.append(i+(j-1)*Lx) # (i,j-1,0) pt
            col_G.append(i+1+j*Lx) # (i+1,j,0) pt 
            data_G.append(3*tc)
            data_G.append(-tc)
            data_G.append(-tc)
            data_G.append(-tc)
            nz+=4


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

  

  #print('row_G=',row_G)
  #print('data_G=',data_G)

  print('Sizes of row_G, col_G, and data_G:')
  print(len(row_G),len(col_G),len(data_G))
  sz=len(data_G)
  print("NL=",NL,' nz=',nz)

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
