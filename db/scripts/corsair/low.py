# collection of generic self-defined functions

import os, sys 				# low level handling, such as command line stuff
import string					# string methods available
import re							# regular expressions
import math						# math functions
import time						# time functions
import anydbm					# index databases (file hash)
from Bio import SeqIO # biopython stuff, to parse fasta files for instance
import warnings
warnings.filterwarnings('ignore', category=DeprecationWarning)


# =============================================================================
def args_file_exists(hash, key):
  if not hash.has_key(key):
    stderr( "file argument \"" + key + "\" missing." )
    return False
  elif not file_exists( hash[key] ):
    stderr( key + " file \"" + hash[key] + "\" does not exist." )
    return False
  return True

# =============================================================================
def args_dir_exists(hash, key):
  if not hash.has_key(key):
    stderr( "dir argument \"" + key + "\"missing." )
    return False
  elif not dir_exists( hash[key] ):
    stderr( key + " dir \"" + hash[key] + "\" does not exist." )
    return False
  return True

# =============================================================================
def overlap_between(pair1, pair2):
  """ checks whether there is an overlap between two start/stop position pairs """
  # ensure that start < stop position
  pair1.sort()
  start1, stop1 = pair1
  pair2.sort()
  start2, stop2 = pair2

  if stop1 < start2: return 0
  elif stop2 < start1: return 0
  else: return 1


# =============================================================================
def humanize_time(secs):
  mins, secs = divmod(secs, 60)
  hours, mins = divmod(mins, 60)
  days, hours = divmod(hours, 24)
  return '%01dd %02d:%02d:%02d' % (days, hours, mins, secs)

# =============================================================================
def sort_by_value(dict):
  """ Returns the keys of dictionary d sorted by their values """
  items=dict.items()
  backitems=[ [v[1],v[0]] for v in items]
  backitems.sort()
  return [ backitems[i][1] for i in range(0,len(backitems))]

# =============================================================================
def convert_time_to_hours( seconds ):
	""" reads float seconds and transforms it into a string of hh:mm:ss """
	hours = math.floor( 1.0*seconds / 3600 )
	seconds = seconds % 3600
	minutes = math.floor( 1.0*seconds / 60 )
	seconds = math.floor(seconds % 60)
	L = []
	L.append( add_leading_zeroes(hours, 2) )
	L.append( add_leading_zeroes(minutes, 2) )
	L.append( add_leading_zeroes(seconds, 2) )
	return string.join( L, ':' )

# =============================================================================
def infomsg( text ):
	""" writes an info message to std.err """
	ctt = time.localtime()
	currenttime = "%s-%s-%s %s:%s:%s" %(ctt[0], ctt[1], ctt[2], ctt[3], add_leading_zeroes(ctt[4],2), add_leading_zeroes(ctt[5],2))
	sys.stderr.write( "<"+currenttime+"> [INFO] " + text + "\n")

# =============================================================================
def info( text ):
	""" writes an info message to std.err using \r, thus is refreshable """
	sys.stderr.write( "\r" + text )

# =============================================================================
def stderr( text ):
	""" just a wrapper for sys.stderr.write() """
	ctt = time.localtime()
	currenttime = "%s-%s-%s %s:%s:%s" %(ctt[0], ctt[1], ctt[2], ctt[3], add_leading_zeroes(ctt[4],2), add_leading_zeroes(ctt[5],2))
	sys.stderr.write( "<"+currenttime+"> ERROR: " + text + "\n")

# =============================================================================	
def stdout( text ):
	""" just a wrapper for sys.stdout.write() """
	sys.stdout.write( text  + "\n")

# =============================================================================
def max( list ):
	""" returns the max value of the list """
	list.sort()
	list.reverse()
	return list[ 0 ]

# =============================================================================
def min( list ):
	""" returns the min value of the list """
	list.sort()
	return list[ 0 ]


# =============================================================================
def split_string_into_lines( string, colwidth ):
	""" splits a string into lines with <colwidth> characters in one line """
	i = 0
	splitstring = ''
	iterations = math.ceil( 1.0 * len( string ) / colwidth )
	for i in range( int ( iterations ) ):
		start = i*colwidth
		end = ((i+1)*colwidth)-1
		splitstring += string[start:end] + '\n'
	if splitstring.endswith( '\n' ):
		splitstring = splitstring.rstrip( '\n' )
	return splitstring


# =============================================================================
def get_global_path( path ):
	if path.startswith( '/' ):
		globalpath = path
	else:
		if path.startswith( './' ):
			path = path[ 2 : ]
		if path == '.':
			path = ''
		curdir = os.getcwd()
		while path.startswith( '../' ):
			curdir = curdir [ : curdir.rindex( '/' ) ]
			path = path[ 3 : ]
		globalpath = curdir + '/' + path
	
	return globalpath


# =============================================================================
def file_exists( path ):
	""" checks whether or not the path points to a valid file """
	if os.path.exists( path ) and os.path.isfile( path ):
		return 1
	else:
		return 0


# =============================================================================
def dir_exists( path ):
	""" checks whether or not the path points to a valid directory """
	if os.path.exists( path ) and os.path.isdir( path ):
		return 1
	else:
		return 0


# =============================================================================
def get_complementary_sequence( s1 ):
	""" returns the complementary sequence for s1. """
	complement = { 'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N' }
	
	s1 = s1.upper()
	s2 = ''
	for i in range( len ( s1 ) ):
		s2 += complement.get( s1[i] )
	
	return s2

# =============================================================================
def get_number_of_nucleotides( string ):
	""" returns the number of nucleotides ([ATCGN-]) within the string. """
	string = string.upper()
	return string.count( 'A' ) + string.count( 'T' ) + string.count( 'C' ) + string.count( 'G' ) + string.count( 'N' ) + string.count( '-' )

# =============================================================================
def write_to_file( file, string ):
	""" writes the given string to the specified file. """
	fw = open( file, 'w' )
	fw.write( string )
	fw.flush()
	fw.close()
	
# =============================================================================
def read_from_file( file ):
	""" reads the content of a file and returns it as a string."""
	fo = open( file, 'r' )
	contents = fo.read()
	fo.close()
	return contents

# =============================================================================
def send_email( addr_from, addr_to, subject, message ):
	""" sends an email with the specific options using sendEmail.pl """
	os.system( 'sendEmail.pl -q -f ' + addr_from + ' -t ' + addr_to + ' -u ' + subject + ' -m ' + message )
	

# =============================================================================
def get_hash_column_from_file( file, column_index ):
  """
  reads a column from a file and returns the values as a hash (key=content, value=numberofappearances).
  if index out of range, this function returns None
  """
  hash = {}
  fo = open( file, 'r' )
  for line in fo:
    line = line.replace( '\n', '' )
    col = line.split()
    if column_index >= len(col):
      stderr( "invalid indexing when retrieving column #%s in file %s" %(column_index,file) )
      return None
    if not hash.has_key( col[column_index] ):
     	hash[ col[column_index] ] = 1
    else:
    	hash[ col[column_index] ] += 1
  fo.close()
  return hash

# =============================================================================
def get_column_from_file( file, column_index ):
  """
  reads a column from a file and returns the values as a list.
  if index out of range, this function returns None
  """
  values = []
  fo = open( file, 'r' )
  for line in fo:
    line = line.replace( '\n', '' )
    col = line.split()
    if column_index >= len(col):
    	stderr( "invalid indexing when retrieving column #%s in file %s" %(column_index,file) )
    	return None
    values.append( col[column_index] )
  fo.close()
  return values
 
 # =============================================================================
def add_leading_zeroes( number, digits ):
 	"""
 	adds leading zeroes to fill up the digits.
 	"""
 	number = int( number )
 	string = str(number)
 	i = 1
 	while i < digits:
 		if number < 10 ** i:
 			string = '0' + string
 		i += 1
 			
 	return string

# =============================================================================
def get_basename( pathandfilename, realbase=0 ):
	"""
	returns the base name of a string containing path and filename
	"""
	path, filename = os.path.split(pathandfilename)
	base, ext = os.path.splitext(filename)
	if realbase:
		if base.find('.') != 0:
			base = base[:base.index('.')]
	return base

# =============================================================================
def sysout( touple ):
	"""
	test routine for updatable text output
	"""
	sys.stdout.write(touple)
	
# =============================================================================
def run_clustalw( input_file, output_format, align_out_file, tree_out_file ):
	"""
	creates an MSA and a tree using clustalw with the specified settings:
	input_file: file containing the sequences to align
	output_format: format in which to produce the MSA.
	fa* = fasta
	phy* = phylip
	default: clustal
	align_out_file: alignment output file
	tree_out_file: tree output file
	
	"""
	output_format = output_format.lower()
	pin, pout = os.popen2("clustalw")
	# load input file
	pin.write( "1\n" + input_file + "\n" )
  # specify output format
	pin.write( "2\n9\n" )
	if output_format.startswith( 'fa' ):
		pin.write( "1\nF\n" )
	elif output_format.startswith( 'phy' ):
		pin.write( "4\n1\n" )
	# else: use default (clustal)
	# align!
	pin.write( "0\n" )
	pin.write( align_out_file + "\n" )
	pin.write( "\n\n\n\n" )
	# generate tree for the MSA
	pin.write( "4\n4\n" )
	pin.write( tree_out_file + "\n" )
	# exit
	pin.write( "\n\nx\n" )
	pin.flush()
	pin.close()
	pout.close()
	
# =============================================================================
def run_tcoffee( input_file, output_format, outfile ):
	"""
	creates an MSA using T-COFFEE with the specified settings:
	input_file: file containing the sequences to align
	output_format: format in which to produce the MSA. (clustalw_aln fasta_aln phylip score_html score_pdf)
	"""
	output_format = output_format.lower()
	os.system("t_coffee -in %s -in Mlalign_id_pair -in Mslow_pair -output %s -outfile %s -cache=no &> /dev/null" %(input_file, output_format, outfile) )
	
# =============================================================================
def run_tranalign( nt_file, aa_file, out_file ):
	"""
	TRANALIGN: Align nucleic coding regions given the aligned proteins.
	nt_file: raw nucleotide sequences
  aa_file: aligned amino acid sequences
  out_file: output file for the nucleotide MSA
	"""
	ok = os.system( "tranalign -asequence " + nt_file 
								  + " -bsequence " + aa_file + " -outseq " + out_file )


# =============================================================================
def catch_bash_cmd_output(command):
	"""
	executes a command and returns the stdout and stderr of that.
	"""
	sin, sout, serr = os.popen3( command )
	sin.close()
	souttext = sout.read()
	serrtext = serr.read()
	sout.close()
	serr.close()
	return souttext, serrtext

# =============================================================================
def index_database_dbm( dbfiles, outfile, type='fasta' ):
	DBM = anydbm.open( outfile, 'c' )
	if type == 'fasta':
		for db in dbfiles:
			handle = open(db)
			for seq_record in SeqIO.parse(handle, "fasta"): DBM[ seq_record.id ] = seq_record.seq.tostring()
			handle.close()
	if type == 'description':
		for db in dbfiles:
			handle = open(db)
			for seq_record in SeqIO.parse(handle, "fasta"): 
				DBM[ seq_record.id ] = seq_record.description[ seq_record.description.index(' ')+1 : ]
			handle.close()
	if type == 'annotation':
		for db in dbfiles:
			handle = open(db)
			for seq_record in SeqIO.parse(handle, "fasta"): 
				DBM[ seq_record.id ] = seq_record.annotations
			handle.close()
	if type == 'name':
		for db in dbfiles:
			handle = open(db)
			for seq_record in SeqIO.parse(handle, "fasta"): 
				DBM[ seq_record.id ] = seq_record.name
			handle.close()
	if type == 'features':
		for db in dbfiles:
			handle = open(db)
			for seq_record in SeqIO.parse(handle, "fasta"): 
				DBM[ seq_record.id ] = seq_record.features
			handle.close()
			
			name
	DBM.close()
	return outfile
