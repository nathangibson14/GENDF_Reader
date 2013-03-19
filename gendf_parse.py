#!/usr/bin/env python

from __future__ import division

import sys
import re
import struct
import math
import h5py
import os

#===============================================================================
# define parseline command for fixed width formats
widths = (11, 11, 11, 11, 11, 11, 4, 2, 3, 5)
formats = ''.join(['{0}s'.format(f) for f in widths])
parseline = struct.Struct(formats).unpack_from

#===============================================================================
def parse_next_line(fh):
    # reads in next line in GENDF file, splits the line based on
    # fixed width columns, and evaluates the columns
    
    line = fh.readline()
    fields = parseline(line)
    
    # replace bad scientific notation with good scientific notation      
    fields = [re.sub(r'(?<=\d)\-','e-',f) for f in fields]
    fields = [re.sub(r'(?<=\d)\+','e+',f) for f in fields]
    
    # convert from string to number
    tmp = []
    for f in fields:
        try:
            tmp.append(eval(f))
        except SyntaxError:
            tmp.append(None)
    fields = tmp
    
    return fields

#===============================================================================
def get_fields(fh,N):
    # get the next N fields from fh
    # function used to deal with the fact the fields might be on
    # multiple lines
    
    fields = []
    for i in range(int(math.ceil(N/6))):
        tmp = parse_next_line(fh)[0:6]
        [fields.append(f) for f in tmp]
        
    return fields[0:N]

#===============================================================================
def read_header(fh,header):        
    
    # title line
    fh.readline()
            
    # CONT record
    fields = parse_next_line(fh)
    header['zaid'] = fields[0]
    header['awr'] = fields[1]
    header['n_dilutions'] = fields[3]
    
    # LIST record
    fields = parse_next_line(fh)
    header['temp'] = fields[0]
    header['NGN'] = fields[2]
    header['NGG'] = fields[3]
    header['NW'] = fields[4]
    
    # Read remainder of LIST
    fields = get_fields(fh,header['NW'])
    i=1
    header['sig0'] = fields[i:i+header['n_dilutions']]
    i += header['n_dilutions']
    header['E'] = fields[i:i+header['NGN']+1]
    i += header['NGN']
    header['E_gamma'] = fields[i+header['NGG']+1]
    
    # Read and check buffer line
    fields = parse_next_line(fh)
    if sum(fields[7:])  != 0:
        raise ValueError
        
    
#===============================================================================
def read_next_record(fh,header,h5file):    
    
    # CONT record
    fields = parse_next_line(fh)
    zaid = fields[0]
    awr = fields[1]
    n_legendre = fields[2]
    n_dilutions = fields[3]
    NGN = fields[5]
    
    # check for end of file
    if fields[6] == 0:
        return 0
    
    # get MF, MT numbers
    MF = fields[7]
    MT = fields[8]
    
    ## only works for MF3 right now
    #if MF != 3:
        #return 0
     
    data = {}
    data['group'] = []
    data['sigma'] = []    
    if MF == 6:
        data['min_group'] = []
    
    
    # read first line in LIST record
    fields = parse_next_line(fh)
    IG1 = fields[5]
    IG = IG1    
    
    while True:
        IG_old = IG
        
        # first line in record            
        IG2LO = fields[3]
        NW = fields[4]
        IG = fields[5]
        
        # get remaining lines in record
        fields = get_fields(fh,NW)
        
        # calculate position after fluxes
        x = n_legendre*n_dilutions
                
        # grab sigma
        data['group'].append(IG)
        data['sigma'].append(fields[x:])                
        if MF == 6:
            data['min_group'].append(IG2LO)
        
        # read next line
        fields = parse_next_line(fh)
        
        # check for end of record
        if fields[8] == 0:
            break
    
    if MF==3:
        MF3toH5(data,h5file,MT,n_legendre,n_dilutions,NGN)
    elif MF==6:
        MF6toH5(data,h5file,MT,n_legendre,n_dilutions,NGN)
    
    print MF, MT    
        
    return 1
    
#===============================================================================
def MF3toH5(data,h5file,MT,n_legendre,n_dilutions,NGN):
    
    # create HDF5 datasets
    grp = h5file["MF3"]
    grp = grp.create_group("MT"+str(MT))
    for i in range(n_legendre):
        P = grp.create_group("P"+str(i))
        for j in range(n_dilutions):
            P.create_dataset("dil"+str(j),(len(data['group']),1), \
                [('group','i'),('sigma','f')])
    
    # put data in HDF5, reverse group numbering          
    for g in range(len(data['group'])):
        y = 0
        for j in range(n_dilutions):
            for i in range(n_legendre):
                grp["P"+str(i)]["dil"+str(j)][-g-1] = \
                    (NGN-data['group'][g]+1,data['sigma'][g][y])
                y+=1

#===============================================================================
def MF6toH5(data,h5file,MT,n_legendre,n_dilutions,NGN):

    # create HDF5 datasets
    grp = h5file["MF6"]
    grp = grp.create_group("MT"+str(MT))
    for i in range(n_legendre):
        P = grp.create_group("P"+str(i))
        for j in range(n_dilutions):
            P.create_group("dil"+str(j))
            P["dil"+str(j)].create_dataset("meta",(len(data['group']),1), \
                [('group','i'),('ijj','i'),('njj','i')])
                
    njj = []
                
    # put ijj,njj data in HDF5
    for g in range(len(data['group'])):
        njj.append(int(len(data['sigma'][g])/(n_legendre*n_dilutions)))
        for j in range(n_dilutions):
            for i in range(n_legendre):
                grp["P"+str(i)]["dil"+str(j)]["meta"][-g-1] = (
                    NGN-data['group'][g]+1,
                    NGN-data['min_group'][g]-njj[-1]+2,
                    njj[-1] )

    # put sigma data in HDF5
    for i in range(n_legendre):
        for j in range(n_dilutions):
            grp["P"+str(i)]["dil"+str(j)].create_dataset("sigma",(sum(njj),1),'f')
    
    x = 0
    for g in reversed(range(len(data["group"]))):
        y = 0
        for gg in range(njj[g]):
            for j in range(n_dilutions):
                for i in range(n_legendre):
                    #~ print g, y
                    grp["P"+str(i)]["dil"+str(j)]["sigma"][x,0] = \
                        data['sigma'][g][y]
                    y+=1
            x += 1

#===============================================================================
def main():

    if len(sys.argv)==2:
        h5filename = sys.argv[1]+'.h5'
    else:
        h5filename = sys.argv[2]+'.h5'

    # delete HDF5 file if it already exists
    try:
        os.remove(h5filename)
    except:
        pass

    # create HDF5 file
    h5file = h5py.File(h5filename,'w-')
    
    # create MF groups
    h5file.create_group("MF3")
    h5file.create_group("MF5")
    h5file.create_group("MF6")
    
    header = {}

    with open(sys.argv[1],'r') as fh:

        read_header(fh,header)
      
        while True:
            ierr = read_next_record(fh,header,h5file)
            if ierr == 0:
                break

#===============================================================================
if __name__ == "__main__":
    main()
