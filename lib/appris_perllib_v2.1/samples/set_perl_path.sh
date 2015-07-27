# setting up perl library path 


# MAIN CONFIGURATION TO CHECKOUT-ROOT
##########################################

# symlink to workdir

export APPRIS_HOME='/Users/jmrodriguez/projects/Encode/release_7'

export PROTEO_HOME='/Users/jmrodriguez/projects/Encode/proteo'
 
#
# symlinkk to used appris-checkouts ( diff. checkouts below ) 
#

export APPRIS_LIB=$APPRIS_HOME/lib/Perl
export APPRIS_CONFIG=$APPRIS_HOME/conf  
export APPRIS_SCRIPTS=$APPRIS_HOME/scripts

export PROTEO_SCRIPTS=$PROTEO_HOME/scripts


echo " "
echo "Checking pathes...."
echo " " 

path_to_check=(  
	'APPRIS_HOME'   'APPRIS_LIB'    'APPRIS_CONFIG'  'APPRIS_SCRIPTS'
	'PROTEO_HOME'   'PROTEO_SCRIPTS'
)
vars_to_check=(
	$APPRIS_HOME    $APPRIS_LIB     $APPRIS_PIPELINE   $APPRIS_SCRIPTS
	$PROTEO_HOME	$PROTEO_SCRIPTS
)

INDEX=${#vars_to_check[@]}
for ((i=0;i<$INDEX;i++));
do
  var=${vars_to_check[${i}]}
  if [ ! -e $var ]; then
     echo "HEY ! One of your pathes does not exist: "
     echo "ERROR: $path_to_check[$i] $var "
     exit
  else
     echo "checking "'$'"$path_to_check[$i]  $var  ======> OK "
  fi
done 

echo " "
echo "Add Appris modules into PERL5LIB env...."
echo " " 

# Add Appris modules 
export PERL5LIB=${APPRIS_LIB}:${PERL5LIB}
echo "Add APPRIS_LIB: ${PERL5LIB}"

echo " " 
echo " " 
echo " " 
echo " " 
