#!/usr/bin/env python
##
##   dV_0p.py plots (V_0 - V_p), planar averaged.
##   V_0 is the total DFT potential of the neutral defect
##   V_p is the pristine DFT potential of the same size cell
##   
import sys, string
import numpy as np


def read_input(file_name):
  cell_dim = 1.0
  factor = None
  charge = -1
  plt_dir = 'a1'
  fp = open(file_name,'r')
  lines = fp.readlines()
  for il in range(len(lines)):
    if "file_neutral" in lines[il]:
      w = lines[il].split("=")
      if len(w) < 2 or len(w) > 3:
        print("ERROR while parsing input file: %s, line: %d"%(file_name,il))
        sys.exit()
      file_n = w[1].split()[0]
    if "file_pristine" in lines[il]:
      w = lines[il].split("=")
      if len(w) < 2 or len(w) > 3:
        print("ERROR while parsing input file: %s, line: %d"%(file_name,il))
        sys.exit()
      file_p = w[1].split()[0]
    if "file_type" in lines[il]:
      w = lines[il].split("=")
      if len(w) < 2 or len(w) > 3:
        print("ERROR while parsing input file: %s, line: %d"%(file_name,il))
        sys.exit()
      file_type = w[1].split()[0]
    if "plt_dir" in lines[il]:
      w = lines[il].split("=")
      if len(w) < 2 or len(w) > 3:
        print("ERROR while parsing input file: %s, line: %d"%(file_name,il))
        sys.exit()
      plt_dir = w[1].split()[0]
    if "factor" in lines[il]:
      w = lines[il].split("=")
      if len(w) < 2 or len(w) > 3:
        print("ERROR while parsing input file: %s, line: %d"%(file_name,il))
        sys.exit()
      factor = w[1].split()[0]
    if "charge" in lines[il]:
      w = lines[il].split("=")
      if len(w) < 2 or len(w) > 3:
        print("ERROR while parsing input file: %s, line: %d"%(file_name,il))
        sys.exit()
      charge = eval(w[1])
  return file_n,file_p,file_type,plt_dir,factor,charge

def write2file(file_name,A,v_a):
  fp = open(file_name,'w')
  if len(A) != len(v_a):
    print("Error: len(A) != len(v_a)")
  for i in range(len(A)):
    fp.write("%4.3f %4.8f\n"%(A[i],v_a[i]))
  fp.close()

def pl_avg_a3(vol,a1_dim,a2_dim,a3_dim,step_l,factor):
  A3 = []
  vol_a3 = np.zeros((a3_dim))
  for k in range(a3_dim):
    Sum1 = 0.
    for i in range(a1_dim):
      for j in range(a2_dim):
        Sum1 = Sum1 + vol[i][j][k]
    vol_a3[k] = Sum1/(a2_dim*a1_dim)
    A3.append(k*step_l)
  if factor == "Ryd":
    vol_a3 = vol_a3*rydberg
  elif factor == "Hartree":
    vol_a3 = vol_a3*hartree
  return vol_a3, np.array(A3)

def pl_avg_a1(vol,a1_dim,a2_dim,a3_dim,step_l,factor):
  A1 = []
  vol_a1 = np.zeros((a1_dim))
  for i in range(a1_dim):
    Sum1 = 0.
    for j in range(a2_dim):
      for k in range(a3_dim):
        Sum1 = Sum1 + vol[i][j][k]
    vol_a1[i] = Sum1/(a2_dim*a3_dim)
    A1.append(i*step_l)
  if factor == "Ryd":
    vol_a1 = vol_a2*rydberg
  elif factor == "Hartree":
    vol_a1 = vol_a2*hartree
  return vol_a1,np.array(A1)

def pl_avg_a2(vol,a1_dim,a2_dim,a3_dim,step_l,factor):
  A2 = []
  vol_a2 = np.zeros((a2_dim))
  for j in range(a2_dim):
    Sum1 = 0.
    for i in range(a1_dim):
      for k in range(a3_dim):
        Sum1 = Sum1 + vol[i][j][k]
    vol_a2[j] = Sum1/(a1_dim*a3_dim)
    A2.append(j*step_l)
  if factor == "Ryd":
    vol_a2 = vol_a2*rydberg
  elif factor == "Hartree":
    vol_a2 = vol_a2*hartree
  return vol_a2,np.array(A2)

def xsf_read(file):
  fp = open(file,"r")
  lines = fp.readlines()
  primvec = []
  primcoord = []
  grid = []
  vol = []
  origin = []
  for i in range(len(lines)):
    if "PRIMVEC" in lines[i]:
      for j in range(3):
        w = lines[i+j+1].split()
        primvec.append([eval(w[0]),eval(w[1]),eval(w[2])])
    if "PRIMCOORD" in lines[i]:
      w = lines[i+1].split()
      na = eval(w[0])
      for j in range(na):
        w = lines[i+j+2].split()
        primcoord.append([w[0],eval(w[1]),eval(w[2]),eval(w[3])])
    if "DATAGRID_3D_" in lines[i]:
      w = lines[i+1].split()
      grid = [eval(w[0]),eval(w[1]),eval(w[2])]
      w = lines[i+2].split()
      origin = [eval(w[0]),eval(w[1]),eval(w[2])]
      # Skip the next 3 lines
      a1_index = 0
      a2_index = 0
      z_index = 0
      vol = np.zeros((grid[0],grid[1],grid[2]))
      for j in range(6,len(lines)):
        if "END_DATAGRID" in lines[i+j]:
          break
        words = lines[i+j].split()
        words = list(filter(bool, words))
        for w in words:
          vol[a1_index][a2_index][z_index] = eval(w)
          a1_index = a1_index + 1
          if a1_index == grid[0]:
            a2_index = a2_index+1
            a1_index = 0
            if a2_index == grid[1]:
              z_index = z_index + 1
              a2_index = 0
  primvec = np.array(primvec)
  step = np.array([primvec[0]/grid[0],primvec[1]/grid[1],primvec[2]/grid[2]]) 
  return na, primcoord, grid, origin, step, vol
  

def cub_read(file):
   ierr = 0
   na = 0
   aspecies = []
   acharge = []
   aposition = []
   grid = []
   origin = []
   step = []
   vol = [[[]]]
   try:
      h = open(file, 'r')
   except:
      ierr = 1
   if ierr == 0:
      for i in range(2):
         s = h.readline()
      s = h.readline()
      t = s.split()
      na = int(t[0])
      origin = []
      for j in range(3):
         origin.append(float(t[j + 1]))
      grid = []
      step = []
      for i in range(3):
         s = h.readline()
         t = s.split()
         grid.append(int(t[0]))
         step.append([])
         for j in range(3):
            step[i].append(float(t[j + 1]))
      for i in range(na):
         s = h.readline()
         t = s.split()
         aspecies.append(int(t[0]))
         acharge.append(float(t[1]))
         aposition.append([])
         for j in range(3):
            aposition[i].append(float(t[j + 2]))
      n = int(grid[0] * grid[1] * ((grid[2] - 1) / 6 + 1))
      i = 0
      j = 0
      k = 0
      for m in range(n):
         s = h.readline()
         t = s.split()
         for l in range(6):
            if k < grid[2]:
               vol[i][j].append(float(t[l]))
               k += 1
         if k == grid[2]:
            k = 0
            j += 1
            if j < grid[1]:
               vol[i].append([])
            else:
               k = 0
               j = 0
               i += 1
               if i < grid[0]:
                  vol.append([[]])
      h.close()
   return ierr, na, aspecies, acharge, aposition, grid, origin, step, vol


if len(sys.argv) == 1:
  print("Please provide input file after call: dV_0p.py input_file")
  sys.exit()
# Read the input file
file_n,file_p,file_type,dir,factor,charge = read_input(sys.argv[1])
if dir != "a1" and dir != "a2" and dir != "a3":
  print("Please specify plt_dir in the input file. It takes a1/a2/a3")
  sys.exit()
if file_type != "cube" and file_type != "xsf":
  print("Please specify file_type in the input file. It takes cube/xsf")

# Read the DFT potential files
if file_type == "cube":
  ierr, na, aspecies, acharge, aposition, grid, origin, step, vol_n = cub_read(file_n)
  ierr, na, aspecies, acharge, aposition, grid, origin, step, vol_p = cub_read(file_p)
elif file_type == "xsf":
  na, primcoord, grid, origin, step, vol_n = xsf_read(file_n)
  na, primcoord, grid, origin, step, vol_p = xsf_read(file_p)

# Compute the difference
vol = (np.array(vol_n) - np.array(vol_p))*charge

print("Shape of the data in the file:", np.shape(vol))
if dir == "a1":
  step_l = np.sqrt(np.dot(step[0],step[0]))
  vol_a1, A1 = pl_avg_a1(vol,grid[0],grid[1],grid[2], step_l,factor)
  write2file("pa_dv0p_a1.plot",A1,vol_a1)
if dir == "a2":
  step_l = np.sqrt(np.dot(step[1],step[1]))
  vol_a2, A2 = pl_avg_a2(vol,grid[0],grid[1],grid[2], step_l,factor)
  write2file("pa_dv0p_a2.plot",A2,vol_a2)
if dir == "a3":
  step_l = np.sqrt(np.dot(step[2],step[2]))
  vol_a3, A3 = pl_avg_a3(vol,grid[0],grid[1],grid[2], step_l,factor)
  write2file("pa_dv0p_a3.plot",A3,vol_a3)
