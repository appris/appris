import sys
import re
import pickle

""" PROCESSING OF ALIGN-TREE FILES FOR SLR """

'''
Create a hash with the specie alignment
'''
def convertFastaToHash(fname):    
	try:
		lines=open(fname,'r').readlines()
	except:
		sys.stderr.write('Can not open the file '+fname+'\n')
		sys.exit()
	seqs={}
	name=''
	seq=''
	for line in lines:
		if line[0]=='>': # assuming fasta file
		   if name != '':
			   seqs[name]=seq
			   seq=''
		   name=line[1:].strip()
		else:
		   seq+=''.join(line.split())
	if name != '':
		seqs[name]=seq
	return seqs

'''
Return a list with the position of residues (without gaps)
'''
def ToPrint(stringa):
	'''ritorna una lista con le posizioni corrispondenti ai gap'''
	lista=[]
	pos=0
	for i in range(len(stringa)):
		if stringa[pos]!='-':
			lista.append(pos)
		pos+=1
	return lista
	
def CheckStop(id,stringa):
	STOPCODON=['UAA','UAG','UAR','TAA','TAG','TAR','TGA','UGA','uaa','uag','uar','taa','tag','tar','tga','uga']
	p = re.compile('(((U|T)A(A|G|R))|((T|U)GA))')
	temp=''
	out=''	
	for pos in range(0,len(stringa),3):
		i=0
		while i<=2:
			if pos+i<len(stringa):
				temp=temp+stringa[pos+i]
			i+=1
		#if not p.match(temp):
		if temp not in STOPCODON:
			out=out+temp
		else:
			if id.startswith(name):
				sys.stderr.write(id+':codon:'+`pos`+'\n')
			out=out+'---'
		temp=''
	return out

'''
Trees should be specified in Newick format (a la PAML). The position
of a sequence in the tree should be specified by its number in the
alignment, not by name as it common in PAML.
'''
def TreeMod(tree, name):
	file=open(tree, 'r')
	line=file.read()
	line=line.lower()
	line=line.replace('\n','')
	line=line.replace('\s','')
	num_species=1
	for id in seqs.keys():
		if id.startswith(name):
			human_id=num_species
		lower_id=id.lower()
		#sys.stderr.write('treemod: id=>'+str(lower_id)+' '+str(num_species)+'\n') # trace log    
		line = line.replace(','+str(lower_id)+':',','+str(num_species)+':')
		line = line.replace('('+str(lower_id)+':','('+str(num_species)+':')
		num_species=num_species+1
		
	return (human_id,num_species,line)


def GetRestSpecies(human_id,num_species,treeFile,alignFile,tree2,num_nucleotide,seqfile_body):

	p = re.compile('(\w+):')
	iterator = p.findall(tree2)

	p2 = re.compile('\d+')

	rest_seqfile_body=''
	for specie in iterator:        
		if not p2.match(specie):          
			rest_seqfile_body+=specie+'\n'+'-'*num_nucleotide+'\n'
			#sys.stderr.write('getrest: id=>'+str(specie)+' '+str(num_species)+'\n') # trace log   
			tree2 = tree2.replace(','+specie+':',','+str(num_species)+':')
			tree2 = tree2.replace('('+specie+':','('+str(num_species)+':')
			num_species=num_species+1

	new=open(treeFile, 'w')
	new.write(' '+str(num_species-1)+' '+str(human_id)+'\n'+tree2)
	new.close()
	
	new=open(alignFile, 'w')
	new.write(str(num_species-1)+' '+str(num_nucleotide)+'\n\n')
	new.write(seqfile_body)
	new.write(rest_seqfile_body)
	new.close()        

	
"-------------------MAIN--------------------"
'''
Return alignment and tree files ready for SLR input

	python pslr.py <Human ID> <Input Alignment file> <Input Tree Align file> <Output Alignment file> <Output Tree Align file>

For example:
	python pslr.py hg19 align.faa tree.nh out_align.faa out_tree.nh
'''

name=sys.argv[1]
alig=sys.argv[2]
tree=sys.argv[3]
outAlign=sys.argv[4]
outTree=sys.argv[5]
seqs=convertFastaToHash(alig)
human_key=name

''' fill with gaps the difference of sequence '''
for id in seqs.keys():
	if id.startswith(name):
		human_key=id
	else:
		if human_key!='':
			if len(seqs[human_key])>len(seqs[id]):
				diff=len(seqs[human_key])-len(seqs[id])
				for times in range(diff):
					seqs[id]=seqs[id]+'-'


''' Return a list with the position of residues (without gaps) '''
seqpos=ToPrint(seqs[human_key])

''' Print the new format of alignment '''
seqfile_body=''
num_nucleotide=0
header=True
num_nucleotide=0
seqfile_body=''
for id in seqs.keys():
	out=''
	for pos in seqpos:
		out=out+seqs[id][pos]
	if header:
		num_nucleotide=len(out)
		header=False
	seqfile_body+=id+'\n'
	seqfile_body+=CheckStop(id,out).lower()+'\n'

(human_id,num_species,tree2)=TreeMod(tree, name)
num_species=GetRestSpecies(human_id,num_species,outTree,outAlign,tree2,num_nucleotide,seqfile_body)

