import numpy as np   
from scipy.sparse import csr_matrix 
from scipy.sparse.linalg import spsolve



def Effective_resistance(V0,R,Lx,Ly):
   
 


    Rp=R
    Rd=R # diagonal resistance 
    #Rd=100000 # diagonal resistance 
    NL=Lx*Ly
    print("Lx=",Lx, "Ly=",Ly, "NL=",NL)


    I=np.zeros(Lx*Ly, float)
    row_G=[]
    col_G=[]
    data_G=[]
    nz=0

    # Linear mapping: (i,j) -> i + jLx 


    ## Left bottom corner point (i=0,j=0) 
    #! No diag res                                                
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

    ## Left top corner point (i=0, j=Ly-1) 
    #! Diag res at (1,Ly-2)                                [(2,Ly-1) in natural indexing]
    row_G.append((Ly-1)*Lx)
    row_G.append((Ly-1)*Lx)
    row_G.append((Ly-1)*Lx)
    row_G.append((Ly-1)*Lx)
    col_G.append((Ly-1)*Lx)
    col_G.append((Ly-1)*Lx+1)
    col_G.append((Ly-2)*Lx)
    col_G.append(1+(Ly-2)*Lx)
    data_G.append(1/R+1/(R+Rp) + 1/(R+Rp)+1/(R+Rd)) # 2.5/R when Rd=Rp=R
    data_G.append(-1/(R+Rp))
    data_G.append(-1/(R+Rp))
    data_G.append(-1/(R+Rd))

    I[(Ly-1)*Lx]=V0/R

    # Right bottom corner point (i=Lx-1, j=0)
    #! Diag res at  (Lx-2,1)                      [(Lx-1,2) in natural indexing]
    row_G.append(Lx-1)
    row_G.append(Lx-1)
    row_G.append(Lx-1)
    row_G.append(Lx-1)
    col_G.append(Lx-1)
    col_G.append(Lx-2)
    col_G.append(2*Lx-1)
    col_G.append(2*Lx-2) # Lx-2 + 1*Lx according to i+j*Lx mapping
    data_G.append(1/R+1/(R+Rp)+1/(R+Rp)+1/(R+Rd))
    data_G.append(-1/(R+Rp))
    data_G.append(-1/(R+Rp))
    data_G.append(-1/(R+Rd))


    ## Right top corner point (i=Lx-1, j=Ly-1)
    #! No diag res 
    row_G.append(Lx*Ly-1)
    row_G.append(Lx*Ly-1)
    row_G.append(Lx*Ly-1)
    col_G.append(Lx*Ly-1)
    col_G.append(Lx*Ly-2)
    col_G.append(Lx*Ly-Lx-1)
    data_G.append(1/R+1/(R+Rp)+1/(R+Rp))
    data_G.append(-1/(R+Rp))
    data_G.append(-1/(R+Rp))



    ## Edge points ---->
    for j in range(1,Ly-1): # j runs from 1 to Ly-2 

        ## Left non-corner edge point (i=0)
        #! Diag res at (1,j-1)                       [(2,j-1) in natural indexing] 
        row_G.append(j*Lx)
        row_G.append(j*Lx)
        row_G.append(j*Lx)
        row_G.append(j*Lx)
        row_G.append(j*Lx)
        col_G.append(j*Lx)
        col_G.append(j*Lx+1)
        col_G.append((j-1)*Lx)
        col_G.append((j+1)*Lx)   
        col_G.append(1+(j-1)*Lx) # Diag res contributes 
        data_G.append(1/R+ 1/(R+Rp)+1/(R+Rp)+1/(R+Rp)+1/(R+Rd))
        data_G.append(-1/(R+Rp))
        data_G.append(-1/(R+Rp))
        data_G.append(-1/(R+Rp))
        data_G.append(-1/(R+Rd)) # Diag res contributes
        I[j*Lx]=1/R

        ## Right non-corner edge point (i=Lx-1)
        #! Diag res at (Lx-2,j+1)                        [(Lx-1,j+1) in natural indexing] 
        row_G.append(j*Lx+Lx-1)
        row_G.append(j*Lx+Lx-1)
        row_G.append(j*Lx+Lx-1)
        row_G.append(j*Lx+Lx-1)
        row_G.append(j*Lx+Lx-1)
        col_G.append(j*Lx+Lx-1)
        col_G.append(j*Lx+Lx-2)
        col_G.append((j-1)*Lx+Lx-1)
        col_G.append((j+1)*Lx+Lx-1)
        col_G.append((j+1)*Lx+Lx-2)
        data_G.append(1/R + 1/(R+Rp)+1/(R+Rp)+1/(R+Rp)+1/(R+Rd)) # Diag res contributes
        data_G.append(-1/(R+Rp))
        data_G.append(-1/(R+Rp))
        data_G.append(-1/(R+Rp))
        data_G.append(-1/(R+Rd))  # Diag res contributes

    for i in range(1,Lx-1):
        ## Bottom non-corner edge point (j=0)
        #! Diag res at (i-1,1)                     [(i-1,2) in natural indexing]
        row_G.append(i)
        row_G.append(i)
        row_G.append(i)
        row_G.append(i)
        row_G.append(i)
        col_G.append(i)
        col_G.append(i-1)
        col_G.append(i+1)
        col_G.append(Lx+i)
        col_G.append(Lx+i-1)
        data_G.append(1/(R+Rp)+ 1/(R+Rp) +1/(R+Rp)+1/(R+Rd))  # Diag res contributes 
        data_G.append(-1/(R+Rp))
        data_G.append(-1/(R+Rp))
        data_G.append(-1/(R+Rp))
        data_G.append(-1/(R+Rd))  # Diag res contributes

        ## Top non-corner edge point (j=Ly-1)
        #! Diag res at (i+1,Ly-2)                     [(i+1,Ly-1) in natural indexing]
        row_G.append((Ly-1)*Lx+i)
        row_G.append((Ly-1)*Lx+i)
        row_G.append((Ly-1)*Lx+i)
        row_G.append((Ly-1)*Lx+i)
        row_G.append((Ly-1)*Lx+i)
        col_G.append((Ly-1)*Lx+i) 
        col_G.append((Ly-1)*Lx+i-1)
        col_G.append((Ly-1)*Lx+i+1)
        col_G.append((Ly-2)*Lx+i) 
        col_G.append((Ly-2)*Lx+i+1) 
        data_G.append(1/(R+Rp) +1/(R+Rp) +1/(R+Rp)+1/(R+Rd))  # Diag res contributes
        data_G.append(-1/(R+Rp))
        data_G.append(-1/(R+Rp))
        data_G.append(-1/(R+Rp))
        data_G.append(-1/(R+Rd))  # Diag res contributes
 
    ## Non border inner points
    for i in range(1,Lx-1):
      for j in range(1,Ly-1):
        #! Diag res at (i-1,j+1) and (i+1,j-1)
        row_G.append(j*Lx+i)
        row_G.append(j*Lx+i)
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
        col_G.append((j+1)*Lx+i-1)
        col_G.append((j-1)*Lx+i+1)
        data_G.append(1/(R+Rp) +1/(R+Rp)+1/(R+Rp) +1/(R+Rp)+2/(R+Rd))  # Diag res contribute
        data_G.append(-1/(R+Rp))
        data_G.append(-1/(R+Rp))
        data_G.append(-1/(R+Rp))
        data_G.append(-1/(R+Rp))
        data_G.append(-1/(R+Rd))  # Diag res contributes
        data_G.append(-1/(R+Rd))  # Diag res contributes


    print('Sizes --')
    print(len(row_G),len(col_G),len(data_G))

    G=csr_matrix((data_G, (row_G,col_G)), shape=(Lx*Ly,Lx*Ly), dtype=float)
                   ###converting the sparse matrix to compressed sparse column format
    V=spsolve(G, I)###solving the matrix

    I_out=0.0
    #####for calculating net current
    for j in range(0,Ly):
        I_out += V[j*Lx+Lx-1]/R

    
    Reff=V0/I_out # Effective resistance 
    print("V0=",V0,"Lx=",Lx,"Ly=",Ly)
    return Reff


