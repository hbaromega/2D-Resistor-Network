def Effective_resistance(V0,R,Lx,Ly):



  NL=Lx*Ly
  print("Lx=",Lx, "Ly=",Ly, "NL=",NL)
  '''m=int(p_latt*NL)
  Q=np.random.choice(NL,m,replace=False)
  def R(n,Q,R_M,R_I):
    if n in Q: return R_M
    else: return R_I'''


  I=np.zeros(Lx*Ly, float)
  row_G=[]
  col_G=[]
  data_G=[]
  nz=0


  for P in range(NL):
    py = int(P/Lx)             # Note this is an integer division
    px = P - py*Lx

    # LEFT BORDER
    if px==0:
      I[P]=V0/R#(P,Q,R_M,R_I)

      if py==0:

        ## Left bottom corner point (i=0,j=0)
        row_G.append(P)
        row_G.append(P)
        row_G.append(P)
        row_G.append(P)
        col_G.append(0)
        col_G.append(1)
        col_G.append(Lx)
        col_G.append(1+Lx)
        data_G.append(2.5/R)
        data_G.append(-0.5/R)#(-1/(R(0,Q,R_M,R_I)+R(1,Q,R_M,R_I)))
        data_G.append(-0.5/R)#(-1/(R(0,Q,R_M,R_I)+R(Lx,Q,R_M,R_I)))
        data_G.append(-0.5/R)
        nz=nz+3

      elif py==Ly-1:

        ## Left top corner point (i=0, j=Ly-1)
        row_G.append(P)
        row_G.append(P)
        row_G.append(P)
        #row_G.append(P)
        #col_G.append( (Ly-1)*Lx )  # (0,Ly-1) --> 0+(Ly-1)Lx=P
        col_G.append(P)  # (0,Ly-1) --> 0+(Ly-1)Lx=P
        #col_G.append( (Ly-1)*Lx+1 ) # (1,Ly-1)-->(Ly-1)Lx+1=P+1
        col_G.append(P+1) # (1,Ly-1)-->(Ly-1)Lx+1=P+1
        #col_G.append((Ly-2)*Lx) # (0,Ly-2) -->  0+(Ly-2)Lx=P-Lx
        col_G.append(P-Lx) # (0,Ly-2) -->  0+(Ly-2)Lx=P-Lx
        #col_G.append(1+(Ly-2)*Lx)
        data_G.append(2/R)#(1/R(P,Q,R_M,R_I)+ 1/(R(P,Q,R_M,R_I)+R(P+1,Q,R_M,R_I)) + 1/(R(P,Q,R_M,R_I)+R(P-Lx,Q,R_M,R_I)) )
        data_G.append(-0.5/R)#( -1/(R(P,Q,R_M,R_I)+R(P+1,Q,R_M,R_I)) )
        data_G.append(-0.5/R)#( -1/(R(P,Q,R_M,R_I)+R(P-Lx,Q,R_M,R_I)) )
        #data_G.append(-0.5/Rd)
        I[(Ly-1)*Lx]=V0/R#((Ly-1)*Lx,Q,R_M,R_I)
        nz=nz+4

      else:

        ## Left non-corner edge  (i=0, j running)
        row_G.append(py*Lx)
        row_G.append(py*Lx)
        row_G.append(py*Lx)
        row_G.append(py*Lx)
        row_G.append(py*Lx)
        col_G.append(py*Lx)
        col_G.append(py*Lx+1)
        col_G.append((py-1)*Lx)
        col_G.append((py+1)*Lx)
        col_G.append(1+(py+1)*Lx)
        data_G.append(2/R)#((0.5/Rd)+1/R(py*Lx,Q,R_M,R_I) +1/(R(py*Lx,Q,R_M,R_I)+R(py*Lx+1,Q,R_M,R_I)) +1/(R(py*Lx,Q,R_M,R_I)+R((py-1)*Lx,Q,R_M,R_I))+ 1/(R(py*Lx,Q,R_M,R_I)+R((py-1)*Lx,Q,R_M,R_I)))
        data_G.append(-0.5/R)#(-1/(R(py*Lx,Q,R_M,R_I)+R(py*Lx+1,Q,R_M,R_I)))
        data_G.append(-0.5/R)#(-1/(R(py*Lx,Q,R_M,R_I)+R((py-1)*Lx,Q,R_M,R_I)))
        data_G.append(-0.5/R)#(-1/(R(py*Lx,Q,R_M,R_I)+R((py+1)*Lx,Q,R_M,R_I)))
        data_G.append(-0.5/R)
        I[py*Lx]=V0/R#(py*Lx,Q,R_M,R_I)
        nz=nz+5


    # RIGHT BORDER
    elif px==Lx-1:

      if py==Ly-1:
        ## Right top corner point (i=Lx-1, j=Ly-1)
        row_G.append(Lx*Ly-1)
        row_G.append(Lx*Ly-1)
        row_G.append(Lx*Ly-1)
        row_G.append(Lx*Ly-1)
        col_G.append(Lx*Ly-1)
        col_G.append(Lx*Ly-2)
        col_G.append(Lx*Ly-Lx-1)
        col_G.append(Lx-2+(Ly-1)*Lx)
        data_G.append(2.5/R)#(0.5/Rd+1/R(Lx*Ly-1,Q,R_M,R_I) + 1/(R(Lx*Ly-1,Q,R_M,R_I)+R(Lx*Ly-2,Q,R_M,R_I)) + 1/(R(Lx*Ly-1,Q,R_M,R_I)+R(Lx*Ly-Lx-1,Q,R_M,R_I)))
        data_G.append(-0.5/R)#(-1/(R(Lx*Ly-1,Q,R_M,R_I)+R(Lx*Ly-2,Q,R_M,R_I)))
        data_G.append(-0.5/R)#(-1/(R(Lx*Ly-1,Q,R_M,R_I)+R(Lx*Ly-Lx-1,Q,R_M,R_I)))
        data_G.append(-0.5/R)
        nz=nz+4

      elif py==0:
        ## Right bottom corner point (i=Lx-1, j=0)
        row_G.append(Lx-1)
        row_G.append(Lx-1)
        row_G.append(Lx-1)
        col_G.append(Lx-1)
        col_G.append(Lx-2)
        col_G.append(2*Lx-1)
        data_G.append(2/R)#( 1/R(Lx-1,Q,R_M,R_I) +1/(R(Lx-1,Q,R_M,R_I)+R(Lx-2,Q,R_M,R_I)) +1/(R(Lx-1,Q,R_M,R_I)+R(2*Lx-1,Q,R_M,R_I)) )
        data_G.append(-0.5/R)#(-1/(R(Lx-1,Q,R_M,R_I)+R(Lx-2,Q,R_M,R_I)))
        data_G.append(-0.5/R)#(-1/(R(Lx-1,Q,R_M,R_I)+R(2*Lx-1,Q,R_M,R_I)))
        nz=nz+3

      else:
        ## Right non corner edge (i=Lx-1, j running)
        row_G.append(py*Lx+Lx-1)
        row_G.append(py*Lx+Lx-1)
        row_G.append(py*Lx+Lx-1)
        row_G.append(py*Lx+Lx-1)
        row_G.append(py*Lx+Lx-1)
        col_G.append(py*Lx+Lx-1)
        col_G.append(py*Lx+Lx-2)
        col_G.append((py-1)*Lx+Lx-1)
        col_G.append((py+1)*Lx+Lx-1)
        col_G.append((py-1)*Lx+Lx-2)
        data_G.append(2/R)#(0.5/Rd+1/R(py*Lx+Lx-1,Q,R_M,R_I) +1/(R(py*Lx+Lx-1,Q,R_M,R_I)+R(py*Lx+Lx-2,Q,R_M,R_I)) +1/(R(py*Lx+Lx-1,Q,R_M,R_I)+R((py-1)*Lx+Lx-1,Q,R_M,R_I)) +1/(R(py*Lx+Lx-1,Q,R_M,R_I)+R((py+1)*Lx+Lx-1,Q,R_M,R_I)))
        data_G.append(-0.5/R)#(-1/(R(py*Lx+Lx-1,Q,R_M,R_I)+R(py*Lx+Lx-2,Q,R_M,R_I)))
        data_G.append(-0.5/R)#(-1/(R(py*Lx+Lx-1,Q,R_M,R_I)+R((py-1)*Lx+Lx-1,Q,R_M,R_I)))
        data_G.append(-0.5/R)#(-1/(R(py*Lx+Lx-1,Q,R_M,R_I)+R((py+1)*Lx+Lx-1,Q,R_M,R_I)))
        data_G.append(-0.5/R)
        nz=nz+5

    elif px>=1 and px<=Lx-2 and py==0: # Bottom non-corner edge
        row_G.append(px)
        row_G.append(px)
        row_G.append(px)
        row_G.append(px)
        row_G.append(px)
        col_G.append(px)
        col_G.append(px-1)
        col_G.append(px+1)
        col_G.append(Lx+px)
        col_G.append(px+1+Lx)
        data_G.append(2/R)#(0.5/Rd+1/(R(px,Q,R_M,R_I)+R(Lx+px,Q,R_M,R_I)) +1/(R(px,Q,R_M,R_I)+R(px-1,Q,R_M,R_I)) +1/(R(px,Q,R_M,R_I)+R(px+1,Q,R_M,R_I)) )
        data_G.append(-0.5/R)#(-1/(R(px,Q,R_M,R_I)+R(px-1,Q,R_M,R_I)))
        data_G.append(-0.5/R)#(-1/(R(px,Q,R_M,R_I)+R(px+1,Q,R_M,R_I)))
        data_G.append(-0.5/R)#(-1/(R(px,Q,R_M,R_I)+R(Lx+px,Q,R_M,R_I)))
        data_G.append(-0.5/R)
        nz=nz+5

    elif px>=1 and px<=Lx-2 and py==Ly-1: # Top non-corner edge
        row_G.append((Ly-1)*Lx+px)
        row_G.append((Ly-1)*Lx+px)
        row_G.append((Ly-1)*Lx+px)
        row_G.append((Ly-1)*Lx+px)
        row_G.append((Ly-1)*Lx+px)
        col_G.append((Ly-1)*Lx+px)
        col_G.append((Ly-1)*Lx+px-1)
        col_G.append((Ly-1)*Lx+px+1)
        col_G.append((Ly-2)*Lx+px)
        col_G.append((Ly-2)*Lx+px-1)
        data_G.append(2/R)#(0.5/Rd+1/(R((Ly-1)*Lx+px,Q,R_M,R_I)+R((Ly-2)*Lx+px,Q,R_M,R_I)) +1/(R((Ly-1)*Lx+px,Q,R_M,R_I)+R((Ly-1)*Lx+px-1,Q,R_M,R_I)) +1/(R((Ly-1)*Lx+px,Q,R_M,R_I)+R((Ly-1)*Lx+px+1,Q,R_M,R_I)) )
        data_G.append(-0.5/R)#(-1/(R((Ly-1)*Lx+px,Q,R_M,R_I)+R((Ly-1)*Lx+px-1,Q,R_M,R_I)))
        data_G.append(-0.5/R)#(-1/(R((Ly-1)*Lx+px,Q,R_M,R_I)+R((Ly-1)*Lx+px+1,Q,R_M,R_I)))
        data_G.append(-0.5/R)#(-1/(R((Ly-1)*Lx+px,Q,R_M,R_I)+R((Ly-2)*Lx+px,Q,R_M,R_I)))
        data_G.append(-0.5/R)
        nz=nz+5
    elif px>=1 and px<=Lx-2 and py>=1 and py<=Ly-2: # Non-Border inner points
        row_G.append(P)
        row_G.append(P)
        row_G.append(P)
        row_G.append(P)
        row_G.append(P)
        row_G.append(P)
        row_G.append(P)
        col_G.append(P)
        col_G.append(P-1)
        col_G.append(P+1)
        col_G.append(P-Lx)
        col_G.append(P+Lx)
        col_G.append(P-1-Lx)
        col_G.append(P+1+Lx)
        data_G.append(3/R)#(1/Rd+1/(R(P,Q,R_M,R_I)+R(P-Lx,Q,R_M,R_I)) +1/(R(P,Q,R_M,R_I)+R(P-1,Q,R_M,R_I)) +1/(R(P,Q,R_M,R_I)+R(P+1,Q,R_M,R_I))+1/(R(P,Q,R_M,R_I)+R(P+Lx,Q,R_M,R_I)) )
        data_G.append(-0.5/R)#(-1/(R(P,Q,R_M,R_I)+R(P-1,Q,R_M,R_I)))
        data_G.append(-0.5/R)#(-1/(R(P,Q,R_M,R_I)+R(P+1,Q,R_M,R_I)))
        data_G.append(-0.5/R)#(-1/(R(P,Q,R_M,R_I)+R(P-Lx,Q,R_M,R_I)))
        data_G.append(-0.5/R)#(-1/(R(P,Q,R_M,R_I)+R(P+Lx,Q,R_M,R_I)))
        data_G.append(-0.5/R)
        data_G.append(-0.5/R)
        nz=nz+7



  print("Number of nonzero elements =",nz) # +1 coz indexing starts from 0
  #print(" Expected: 5*Lx*Ly - 2*Lx-2*Ly =", 5*Lx*Ly - 2*Lx-2*Ly)

  values=np.zeros(nz,float)
  columns=np.zeros(nz, int)
  rows=np.zeros(nz, int)
  rowIndex=np.zeros(NL+1,int)

  print("NL=",NL)
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

  print("rowIndex[Lx], N_nz =",rowIndex[Lx], nz)




  G=csr_matrix((data_G, (row_G,col_G)), shape=(Lx*Ly,Lx*Ly), dtype=float)# converting the sparse matrix to compressed sparse column format
  #G = csr_matrix((values, columns, rowIndex), shape=(NL, NL))
  #G = csr_matrix((data_G, col_G, rowIndex), shape=(NL, NL))
  #G = csr_matrix((values, (columns, rows)), shape=(NL, NL))


  V=spsolve(G, I) # Solving the matrix eq. G V = I  [Ax=b format]
  V1=np.zeros((Lx,Ly),float)
  for iter in range(Lx*Ly):
    j=int(iter/Lx)
    i=iter-j*Lx
    V1[i,j]=V[iter]




  # Calculating net current
  px=0
  I_in=0.0
  for py in range(Ly):
                 P = py*Lx + px
                 I_in += (V0-V[P])/R#(P,Q,R_M,R_I)
  #print("Total incoming current=",I_in)

  px=Lx-1
  I_out=0.0
  for py in range(Ly):
                 P = py*Lx + px
                 I_out += V[P]/R#(P,Q,R_M,R_I)
  '''# Calculating net current

  px=0
  I_in=0.0
  for py in range(Ly):
                 P = py*Lx + px
                 I_in += (V0-V[P])/R(P,Q,R_M,R_I)
  print("Total incoming current=",I_in)

  px=Lx-1
  I_out=0.0
  for py in range(Ly):
                 P = py*Lx + px
                 I_out += V[P]/R(P,Q,R_M,R_I)
  print("Total outgoing current=",I_out)

  print()'''


  Reff=V0/I_out # Effective resistance
  print("V0=",V0,"Lx=",Lx,"Ly=",Ly)
  print("Effective resistance =", Reff)
  #print("Expected analytical value =",2*R*Lx/float(Ly))
  return Reff
  
if __name__ == '__main__':
    main()

print("Done!")



