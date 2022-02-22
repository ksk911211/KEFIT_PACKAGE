#!/usr/bin/env python
""" 
A module for data interfacing to MDSplus
"""
from __future__   import print_function, unicode_literals, absolute_import

__all__=['mds'];
__version__='1.0';

import ctypes as _C
from MDSplus import Connection
#from MDSplus._mdsshr import _load_library, MdsException 
from MDSplus._mdsshr import MdsshrException

#ConnectToMds=_load_library('MdsIpShr').ConnectToMds
#ConnectToMds.argtypes=[_C.c_char_p]
#DisconnectFromMds = _load_library('MdsIpShr').DisconnectFromMds
#DisconnectFromMds.argtypes = [_C.c_int]

class _Connection(Connection):
    """
    Updating 'Connection' class in 'MDSplus' to manange the disConnection when it
    becomes unused anymore
    - modified at Nov 2014 from the one designed by D. K. Oh at Aug 2012
    """
    __version__=__version__;
    def __del__(self):
        self.closeConnection()
        
    def closeConnection(self):
        if self.socket != -1:
             if False: #DisconnectFromMds(self.socket) == 0: 
                raise Exception("Error in disconnection")
             else:
                self.socket = -1

#--------------------------------------------------------------------------------#
import subprocess as sub
class mds():
   """
   A modified mdsplus routine by adopting the 'Connection' class modified by 
   D.K. Oh at Aug 2012 from the original one designed by Y.M. Jeon
   
   Example)
   >> from MDS import *           # Import required model 'mds'
   >> g=mds('kstar',5947);        # Connect to mdsplus server and open a specified tree
                                  # For another mdsplus server, 
                                  # >> g=mds('kstar',5947,server='172.17.100.200:8003');
   >> t,v=ip=g.get('\\pcrc03');
   >> plot(t,v);
   
   HISTORY
   - A class 'Connection' updated by adopting the one designed by D.K. Oh at Aug 2012
   - get() function modified to handle a numeric (constant) data         2012-1115 YMJ
   """
   __version__=__version__;
   
#   def __init__(self,tree=None,shot=None,server='172.17.250.21:8005'):
#   def __init__(self,tree=None,shot=None,server='localhost:8005'):
   def __init__(self,tree=None,shot=None,server='172.17.100.200:8005'):
       """ 
       Create an instance of MDSplus.Connection class
         tree     : tree name. Default='kstar'
         shot     : shot number.
         server   : server name. If None, then it is automatically configured
                    172.17.100.200:8005 (iKSTAR), 172.17.100.200:8300 (Control-room)
       Here either 'tree' or 'shot is None, then just do connect to the MDSplus
       server, while it does connect and open the tree for given shot number 
       if both are given correctly.
       """
       if(server): self.server=server;
       else:
          ans=sub.check_output(["/sbin/route","-n"]);
          ans=ans.split("\n");
          for line in ans:
              if('UG ' in line):
                 str0=line.split();
                 ipgate=str0[1];
          if('172.17' in ipgate): self.server='172.17.100.200:8300'; # in control-room
          else:                   self.server='172.17.100.200:8005'; # ikstar
          
       try:
          self.__G__=_Connection(self.server);
          if(tree is not None) and (shot is not None):
             self.open(tree,shot);
       except MdsshrException:
          raise MdsshrException("Error in connection to %s" %(server));
       except:
          raise Exception("Unknown error in connection to %s" %(server));
   
   def __enter__(self): return self;
   
   def __exit__(self,exc_type,exc_val,exc_tb): 
       self.close();
       self.disconnect();

   def open(self,tree,shot):
       """ 
       Open a tree for given shot
         tree    : tree name. Default is 'kstar'
         shot    : shot number
       """
       self.__G__.openTree(tree,shot);
       self.tree=tree;
       self.shot=shot;
       
   def close(self,tree=None,shot=None):
       if tree is None: tree=self.tree;
       if shot is None: shot=self.shot;
       if (shot is not None) and (tree is not None):
          try:
             self.__G__.closeTree(tree,shot);
          except:
             raise MdsshrException("Error in close(): unknown error");
       
   def disconnect(self): self.__G__.closeConnection();
   
   def get(self,sig):
       """
       Get the data for given 'SIG' and return the data in a form of (time, value)
         sig     : a statement for execution. Ex) '\\pcrc03'
       """
       from scipy import array
       try:
          v=self.__G__.get(sig).data();
       except: # no data available
          t=array([]);
          v=array([]);
          return t,v;
          
       try:
          t=self.__G__.get('dim_of(%s)' %(sig)).data();
       except: # no time data available. So it could be a constant or some others
          v=array([v]); # assuming one single numeric for a while
          t=array([]);
       
       return t,v;
   
   def get_data(self,sig):
       """ it equals with 'mds.get(sig).data()' """
       return self.__G__.get(sig).data();
   
   def time(self):
       try:
          t_time = self.__G__.get('\T0_STR').data()
       except:
          t_time = None            
          raise MdsshrException("Error in getting time-info of the shot '%s'" %(self.shot));
       return t_time

#--------------------------------------------------------------------------------#
if __name__=='__main__':
    """  """
    from matplotlib.pyplot import plot, xlim, show
    
    print("----------------------------------------------------------------------");
    print("   Example for how to use it");
    print("----------------------------------------------------------------------");
    print(">>> from MDS import mds");
    print(">>> g=mds('kstar',7821);");
    g=mds('kstar',7821);
    print(">>> Da=g.get('\\tor_ha11')");
    Da=g.get('\\tor_ha11');
    print(">>> t,v=Da");
    print(">>> plot(t,v); xlim(0.0,9.0);");
    t,v = Da;
    plot(t,v); 
    xlim(0.0,9.0);
    show();
    
    
