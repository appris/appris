#!/usr/bin/python

import os, sys     # low level handling, such as command line stuff
import string      # string methods available
import re          # regular expressions
import getopt      # comand line argument handling
from low import *  # custom functions, written by myself


# =============================================================================  
def show_help( ):
  """ displays the program parameter list and usage information """
  stdout( "usage: " + sys.argv[0] + " -f <path>" )
  stdout( " " )
  stdout( " option    description" )
  stdout( " -h        help (this text here)" )
  stdout( " -f        uniprot dat file" )
  stdout( " " )
  sys.exit(1)

# =============================================================================
def handle_arguments():
  """ verifies the presence of all necessary arguments and returns the data dir """

  if len ( sys.argv ) == 1:
    stderr( "no arguments provided." )
    show_help()  
  
  try: # check for the right arguments
    keys, values = getopt.getopt( sys.argv[1:], "hf:" )
  except getopt.GetoptError:
    stderr( "invalid arguments provided." )
    show_help()

  args = {}
  for key, value in keys:
    if key == '-f': args['datfile'] = value
    
  for key in ['datfile']:
    if key.endswith("file"):
      if not args_file_exists(args, key): show_help()
    elif key.endswith("dir"):
      if not args_dir_exists(args, key): show_help()
  return args

# =============================================================================  
def parse_until_doubleslash(fo):
  hash, end = {}, False
  line = fo.readline().strip()
  while not line.startswith("//"):
    if len(line) == 0:
      end = True
      break
    if len(line.split(" ", 1)[0]) != 2:
      key = "SEQ"
      value = line.strip().replace(" ", "")
    else:
      cols =  [e.strip() for e in line.split(" ", 1)]
      if len(cols) != 2: 
        line = fo.readline().strip()
        continue
      key, value = [e.strip() for e in line.split(" ", 1)]
    if not hash.has_key(key): hash[key] = ""
    if key != "SEQ" and len(hash[key]) > 0 and hash[key][-1] != " " and not value.startswith(" "): hash[key] += " "
    hash[key] += value
    line = fo.readline().strip()
  return hash, end
  
# =============================================================================
# === MAIN ====================================================================
# =============================================================================
def main( args ):
  fo = open(args['datfile'])
  while 1:
    hash, end = parse_until_doubleslash(fo)
    if end: break
    print ">" + hash["ID"].split()[0] + "|" + hash["AC"] + " " + hash["OC"]
    print hash["SEQ"]
  fo.close()

# =============================================================================
args = handle_arguments()
main( args )

