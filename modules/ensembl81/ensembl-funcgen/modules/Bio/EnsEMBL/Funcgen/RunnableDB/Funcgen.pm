=pod

=head1 NAME

Bio::EnsEMBL::Funcgen::RunnableDB::Funcgen

=head1 DESCRIPTION

'Funcgen' is a base class for other runnables of the Funcgen Hive Pipeline
It performs common tasks such as connecting to the EFG DB etc...

=cut

package Bio::EnsEMBL::Funcgen::RunnableDB::Funcgen;

use warnings;
use strict;

use Bio::EnsEMBL::Funcgen::Utils::Helper;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Funcgen::InputSet;
use Bio::EnsEMBL::Funcgen::DataSet;
use Bio::EnsEMBL::Funcgen::FeatureSet;
use Bio::EnsEMBL::Funcgen::AnnotatedFeature;
#use Bio::EnsEMBL::Hive::DBSQL::AnalysisDataAdaptor;
use base ('Bio::EnsEMBL::Hive::Process');

use Bio::EnsEMBL::Utils::Exception qw(throw warning stack_trace_dump);
use Data::Dumper;

#This defines a set of parameters based on given parameters from the pipeline:
sub fetch_input {   # nothing to fetch... just the DB parameters...
  my $self = shift @_;

  #An example of debug, in case needed
  #print Dumper $self->param('dnadb');
  if(!$self->param('bin_dir')){ throw "Folder with funcgen binaries bin_dir required"; }
  $self->_bin_dir($self->param('bin_dir'));

  my $dnadb_params = $self->param('dnadb') || throw "No parameters for Core DB";
  my $efgdb_params = $self->param('efgdb') || throw "No parameters for EFG DB";

  #Get efg connection, otherwise fail..
  eval{
	$self->_efgdba(Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor->new
				   (
					   %{ $efgdb_params },
					#let efg dba hanle dnadb
					{
					 -dnadb_name => $dnadb_params->{-dbname},
					 -dnadb_user => $dnadb_params->{-user},
					 -dnadb_host => $dnadb_params->{-host},
					 -dnadb_port => $dnadb_params->{-port},
					 -dnadb_pass => $dnadb_params->{-pass},

					}
				   ));

	#Actually test connections
	$self->_efgdba->dbc->db_handle;
	$self->_efgdba->dnadb->dbc->db_handle;
  };

  if($@) { throw "Error creating the EFG DBAdaptor and/or dna DBAdaptor $@";  }

  #Set some params
  my $cell_type       = $self->param('cell_type')       || throw "No cell_type given";
  my $feature_type    = $self->param('feature_type')    || throw "No feature_type given";
  my $experiment_name = $self->param('experiment_name') || throw "No experiment_name given";
  $self->_experiment_name($experiment_name);
  my $set_name =  $self->param('set_name') || $cell_type."_".$feature_type."_".$experiment_name;
  $self->_set_name($set_name);
  my $group_name = $self->param('group') || 'efg';
  my $species = $self->param('species') || throw "No species defined";
  $self->_species($species);
  my $assembly = $self->param('assembly') || throw "No assembly version given";
  $self->_assembly($assembly);
  my $file_type = $self->param('file_type') || throw "No file type given";
  $self->_file_type($file_type);
  my $work_dir = $self->param('work_dir') || throw "'work_dir' is a required parameter";
  $self->_work_dir($work_dir);



  #Configure DBAdaptors
  my $efgdba = $self->_efgdba();
  #To avoid farm issues...
  $efgdba->dbc->disconnect_when_inactive(1);
  $efgdba->dnadb->dbc->disconnect_when_inactive(1);

  #Fetch & Set object params
  #CellType
  my $cta    = $efgdba->get_CellTypeAdaptor();
  my $ct_obj = $cta->fetch_by_name($cell_type);
  if(!$ct_obj){ throw "Cell type $cell_type does not exist in the database";  }
  $self->_cell_type($ct_obj);

  #FeatureType
  my $fta    = $efgdba->get_FeatureTypeAdaptor();
  my $ft_obj = $fta->fetch_by_name($feature_type);
  if(!$ft_obj){ throw "Feature type $feature_type does not exist in the database";  }
  $self->_feature_type($ft_obj);

  #ExperimentalGroup
  my $ega = $efgdba->get_ExperimentalGroupAdaptor();
  my $eg_obj = $ega->fetch_by_name($group_name);
  if(!$eg_obj){ throw "Experimental Group $group_name does not exist in the database";  }
  $self->_group($eg_obj);


  if($file_type eq 'sam' || $file_type eq 'bam'){
    #Change the directory structure so it will agree with the rest, without the need to do uc()
    my $sam_header = $self->_work_dir()."/sam_header/".$species."/".$species."_";
    $sam_header .= $ct_obj->gender() ? $ct_obj->gender() : 'male';
    #Carefull with naming standards...
    #$sam_header .= "_".$assembly."_unmasked.fa.fai";
    $sam_header .= "_".$assembly."_unmasked.fasta.fai";
    $self->_sam_header($sam_header);
  }

  #Work with conventions here too?? work_dir/output/dbname ??
  my $output_dir = $self->param('output_dir') || throw "'output_dir' is a required parameter";
  $self->_output_dir($output_dir."/".$experiment_name);

  return 1;
}


sub run {
  my $self = shift @_;

  return 1;
}


sub write_output {
  my $self = shift @_;

  return 1;

}


#Private Function to check and create Experiment and Feature/Data sets as needed
#Requires some global parameters that are not set in Funcgen->fetch_input, such as
#'analysis', 'feature_set_name', 'data_set_name' (these could be given as local parameters...)
sub _check_Experiment {
  #Todo make it more generic and accept multiple input_subsets
  #Also maybe pass parameters as hash list...
  my ($self, $analysis, $input_subset, $fset_name) = @_;

  #Global parameters set in Funcgen->fetch_input
  my $efgdba        = $self->_efgdba();
  my $set_name      = $self->_set_name();
  my $group         = $self->_group();
  my $cell_type     = $self->_cell_type();
  my $feature_type  = $self->_feature_type();

  my $iset_name = $set_name;
  my $dset_name = $fset_name;

  # set experiment: Reuse if already exists? (This comes from result sets)
  my $ea  = $efgdba->get_ExperimentAdaptor;
  my $exp = $ea->fetch_by_name($set_name);

  my @date  = (localtime)[5,4,3];
  $date[0] += 1900;
  $date[1]++;

  if (! defined $exp) {

    #Group needs to be set manually, like Cell_Type and Feature_Type
    #Do not create Group on the fly here, as it will cause concurrency issues...
    $exp = Bio::EnsEMBL::Funcgen::Experiment->new (
       -NAME => $set_name,
       -EXPERIMENTAL_GROUP => $group,
       -DATE => join('-', @date),
       -PRIMARY_DESIGN_TYPE => 'binding_site_identification',
       -ADAPTOR => $ea,
      );

    ($exp) =  @{$ea->store($exp)};

  }
  throw("Can't create experiment $set_name ") unless $exp;

  my $isa  = $efgdba->get_InputSetAdaptor();
  my $iset = $isa->fetch_by_name($iset_name);

  if (! defined $iset){

    $iset = Bio::EnsEMBL::Funcgen::InputSet->new (
       -name         => $iset_name,
       -experiment   => $exp,
       -feature_type => $feature_type,
       -cell_type    => $cell_type,
       -vendor       => 'SOLEXA',
       -format       => 'SEQUENCING',
       -feature_class => 'result',
       # Analysis is not being used??
       -analysis     => $analysis,
      );
    warn "Storing new InputSet:\t$iset_name\n";
    ($iset)  = @{$isa->store($iset)};

#$iset->add_new_subset($input_subset);
#$iset->adaptor->store_InputSubsets($iset->get_InputSubsets);
    $input_subset = Bio::EnsEMBL::Funcgen::InputSubset->new (
       -name      => $input_subset,
       -input_set => $iset,
# Analysis is not being used??
#-analysis     => $self->feature_analysis,
      );
  } else {

    #We only expect one subset here (? why??)...
    #shouldn't we be adding the control file also when used?? But this is SWEmbl-specific...
    #And it should be the same file name...
    #Maybe do some file checking  here???
    warn "InputSet already exists:\t$iset_name\n";
    my @issets = @{$iset->get_InputSubsets};

    #if(scalar(@issets) > 1){
    #  throw("InputSet $iset_name has more than one InputSubset:\t".join("\t", (map $_->name, @issets)));
    #} elsif((scalar(@issets) == 1) && ($issets[0]->name ne $self->param('input_file'))){
    #  throw("InputSet $iset_name already has an InputSubset(".$issets[0]->name.") which does not match ".$self->param('input_file'));
    #} elsif(scalar(@issets) == 0){ #we can just add this InputSubset
    #  $iset->add_new_subset($self->input_id);
    #  $iset->adaptor->store_InputSubsets($iset->get_InputSubsets);
    #}

    if(scalar(@issets)==0){
      #we can just add this InputSubset. Add an extra 'input:' as prefix?
#$iset->add_new_subset($input_subset);
#$iset->adaptor->store_InputSubsets($iset->get_InputSubsets);
      $input_subset = Bio::EnsEMBL::Funcgen::InputSubset->new (
         -name         => $input_subset,
         -input_set     => $iset,
# Analysis is not being used??
#-analysis     => $self->feature_analysis,
        );
    } else {
      #warn("Need to uncomment this section!! - it was commented just for testing purposes!!");
      #we just need to check if our file(s) is(are) already here...
      if(!$iset->get_subset_by_name($input_subset)){
      	#throw("InputSet $iset_name has InputSubsets(".join("\t", (map $_->name, @issets)).") which do not match ".$input_subset);
	#warn("InputSet $iset_name has InputSubsets(".join("\t", (map $_->name, @issets)).") which do not match ".$input_subset);
      }
    }
  }

  my $fsa = $efgdba->get_FeatureSetAdaptor();
  my $fset = $fsa->fetch_by_name($fset_name);

  if ( ! defined $fset ) {

    $fset = Bio::EnsEMBL::Funcgen::FeatureSet->new(
       -analysis      => $analysis,
       -feature_type  => $feature_type,
       -cell_type     => $cell_type,
       -name          => $fset_name,
       -feature_class => 'annotated',
#       -experiment_id => $exp->dbID,
       #The adaptor is needed to store!
       -adaptor       => $fsa

      );

    warn "Storing new FeatureSet:\t$fset_name\n";
    ($fset) = @{$fsa->store($fset)};

  }
  else {
    warn "FeatureSet already exists:\t$fset_name\n";

    if(@{$efgdba->get_AnnotatedFeatureAdaptor->fetch_all_by_FeatureSets([$fset])}){
      throw "Feature Set $set_name already contains data. Please rollback before rerunning";
    }

  }

  my $dsa = $efgdba->get_DataSetAdaptor;
  my $dset = $dsa->fetch_by_name($dset_name);

    if ( ! defined $dset ) {
      $dset = Bio::EnsEMBL::Funcgen::DataSet->new (
	 -SUPPORTING_SETS     => [$iset],
	 -FEATURE_SET         => $fset,
	 -DISPLAYABLE         => 1,
	 -NAME                => $dset_name,
	 -SUPPORTING_SET_TYPE => 'input',
	);

      warn "Storing new DataSet:\t$dset_name\n";
      ($dset) = @{$dsa->store($dset)}
    }
  else {

    warn "DataSet already exists:\t$dset_name\n";

    # need to check whether InputSets and supporting_sets are the same and
    # possibly add InputSet to supporting_sets

    my $ssets = $dset->get_supporting_sets();

    my %ssets_dbIDs = ();
    map { $ssets_dbIDs{$_->dbID}='' } (@{$ssets});
    $dset->add_supporting_sets([ $iset ]) if (! exists $ssets_dbIDs{$iset->dbID});

  }

}


#Private Generic getter and setter
sub _getter_setter {
  my ($self, $param_name, $param_value) = @_;
  if(!$param_name){ return undef; }
  if(!$param_value){
    $param_value = $self->param($param_name);
  } else {
    $self->param($param_name, $param_value);
  }
  return $param_value;
}

# Private getter / setters : Maybe do some validation in some cases...

#Private getter / setter to the EFG DB Adaptor
sub _efgdba {
  return $_[0]->_getter_setter('efgdb',$_[1]);
}

#Private getter / setter to the Core DB Adaptor
sub _dnadba {
  return $_[0]->_getter_setter('dnadb',$_[1]);
}

#Private getter / setter to the Cell Type object
sub _cell_type {
  return $_[0]->_getter_setter('cell_type',$_[1]);
}

#Private getter / setter to the Feature Type object
sub _feature_type {
  return $_[0]->_getter_setter('feature_type',$_[1]);
}

#Private getter / setter to the Species name
sub _species {
  return $_[0]->_getter_setter('species',$_[1]);
}

#Private getter / setter to the assembly name
sub _assembly {
  return $_[0]->_getter_setter('assembly',$_[1]);
}

#Private getter / setter to the Analysis object
sub _analysis {
  return $_[0]->_getter_setter('analysis',$_[1]);
}

#Private getter / setter to the Group
sub _group {
  return $_[0]->_getter_setter('group',$_[1]);
}

#Private getter / setter to the Experiment Name (do not mix with the Set Name)
sub _experiment_name {
  return $_[0]->_getter_setter('experiment_name',$_[1]);
}

#Private getter / setter to the Set Name
sub _set_name {
  return $_[0]->_getter_setter('set_name',$_[1]);
}

#Private getter / setter to the file type
sub _file_type {
  return $_[0]->_getter_setter('file_type',$_[1]);
}

#Private getter / setter to the sam header (only set when file type is sam)
sub _sam_header {
  return $_[0]->_getter_setter('sam_header',$_[1]);
}

#Private getter / setter to the work folder
sub _work_dir {
  return $_[0]->_getter_setter('work_dir',$_[1]);
}

#Private getter / setter to the output folder
sub _output_dir {
  return $_[0]->_getter_setter('output_dir',$_[1]);
}

#Private getter / setter to the bin folder
sub _bin_dir {
  return $_[0]->_getter_setter('bin_dir',$_[1]);
}

1;
