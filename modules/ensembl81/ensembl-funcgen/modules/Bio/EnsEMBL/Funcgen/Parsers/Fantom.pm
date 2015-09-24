=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=cut

package Bio::EnsEMBL::Funcgen::Parsers::Fantom;

use strict;
use warnings;
use feature qw(say);
use Bio::EnsEMBL::DBEntry;
use Bio::EnsEMBL::Funcgen::FeatureType;
use Bio::EnsEMBL::Utils::Exception qw( throw );
use Bio::EnsEMBL::Funcgen::Utils::EFGUtils  qw(add_external_db dump_data);


use parent qw( Bio::EnsEMBL::Funcgen::Parsers::BaseExternalParser );

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;



  my $self = $class->SUPER::new(@_, type => 'Fantom');
  # say dump_data($self,1,1);die;

  # Set default feature_type and feature_set config
  $self->{static_config}{feature_types} = {
    (
      'FANTOM'   => {
        -name         => 'FANTOM',
        -class        => 'DNA',
        -description  => 'FANTOM enhancers and TSS',
        -so_accession => 'SO:0000110',
        -so_name      => 'sequence_feature',
        },
      'FANTOM_enhancer_robust'   => {
        -name         => 'FANTOM robust enhancer',
        -class        => 'Enhancer',
        -description  => 'FANTOM enhancers, robust',
        -so_accession => 'SO:0000165',
        -so_name      => 'enhancer',
        },
      'FANTOM_enhancer_permissive'   => {
        -name         => 'FANTOM permissive enhancer',
        -class        => 'Enhancer',
        -description  => 'FANTOM enhancers, permissive',
        -so_accession => 'SO:0000165',
        -so_name      => 'enhancer',
        },
      'FANTOM_tss_strict'   => {
        -name         => 'FANTOM TSS strict',
        -class        => 'Transcription Start Site',
        -description  => 'FANTOM TSS, strict',
        -so_accession => 'SO:0000315',
        -so_name      => 'transcription_start_site',
        },
      'FANTOM_tss_relaxed'   => {
        -name         => 'FANTOM TSS relaxed',
        -class        => 'Transcription Start Site',
        -description  => 'FANTOM TSS, relaxed',
        -so_accession => 'SO:0000315',
        -so_name      => 'transcription_start_site',
        },
    )
  };

  $self->{static_config}{analysis} = {
    FANTOM =>  {
      -logic_name    => 'FANTOM_v5',
      -description   => '<a>FANTOM5</a> functional annotation (http://fantom.gsc.riken.jp/5/)',
      -display_label => 'FANTOM5',
      -displayable   => 1,
    },
  };
# "INSERT INTO analysis (logic_name, created, db_version) values ('FANTOM_v5', NOW(), 5)"


  # This is used as the entry point to store/validate
  # So all of the above needs to be referenced in here
  $self->{static_config}{feature_sets} = {
    'FANTOM' =>{
      #Stored in this order

      # Entries here are flexible
      # Can be omited if defined in feature_set
      # top level analysis/feature_types definition required if no DB defaults available
      # These can be a ref to the whole or subset of the top level analysis/feature_types hash
      # A key with an empty hash or undef(with or without a matching key in the top level analysis/feature_types hash

      #analysis      => $self->{static_config}{analysis},
      feature_types => $self->{static_config}{feature_types},

      # feature_type and analysis values must be string key to top level hash
      # This wont work for feature_types as they are not unique by name!!!!!!
      # This is why we have top level hash where we can define a unique compound key name

      feature_set   =>{
        -feature_type      => 'FANTOM', #feature_types config key name not object
        -display_label     => 'FANTOM predictions',
        -description       => 'FANTOM predictions',
        -analysis          => 'FANTOM_v5', #analysis config key name not object
      },
    }
  };

  # $self->validate_and_store_feature_types;
  $self->validate_and_store_config([keys %{$self->{static_config}{feature_sets}}]);
  $self->set_feature_sets;

  return $self;
}

# RGB for enhancers

my $dispatch = {
  
  enhancers => sub {
    my ($line) = @_;
    my @line = split /\t/, $line;
    my $data = {};
    $line[0] =~ /chr([0-9XYMT]{1,2})/;
    $data->{chr}    = $1;
    $data->{chr}    = 'MT' if($data->{chr} eq 'M');
    $data->{start}  = $line[1] + 1; 
    $data->{end}    = $line[2]; 
    $data->{name}   = $line[3];
    $data->{score}  = $line[4];
    $data->{strand} = _assign_strand($line[5]);
    $data->{rgb}    = $line[8];
    return $data;
  }, # enhancers

  tss => sub{
    my ($line) = @_;
    my @line = split /\t/, $line;
    my $data = {};
    $line[0] =~ /chr([0-9XYMT]{1,2})/;
    $data->{chr}    = $1;
    $data->{chr}   = 'MT' if($data->{chr} eq 'M');
    $data->{start}  = $line[1] + 1; 
    $data->{end}    = $line[2];

    $line[3] =~ /^(p\d*?)@(\w*?),(0\.\d*?)$/;
    $data->{peak}   = $1;
    $data->{gene}   = $2;
    $data->{anTss}  = $3;

    $data->{name}   = $line[3];
    $data->{score}  = $line[4];
    $data->{strand} = _assign_strand($line[5]);
    $data->{rgb}    = $line[8];
    return $data;
  }, #tss
}; # dispatch

sub _assign_strand {
  my ($symbol) = @_;
  my $strand;

  if($symbol eq '.'){
    $strand = 0;
  }
  elsif($symbol eq '+'){
    $strand = 1;
  }
  elsif($symbol eq '-'){
    $strand = -1;
  }
  else{
    throw("Implement strand symbol '$symbol'");
  }
  return $strand;
}

sub parse_and_load {
  my ($self, $files, $old_assembly, $new_assembly) = @_;


  if(scalar @$files < 1){
    throw("No files provided ");
  }
  
  # permissive_enhancers.bed

  # Set release to the release TarBase used to map their miRNA targets 
  my $external_db_name;
  my $external_db_release;
  
  
  # my $ex_Db = add_external_db(
  #   $self->db, 
  #   $external_db_name,
  #   $external_db_release,
  #   'EnsemblGene',
  #   );



  my $extfeat_a   = $self->db->get_ExternalFeatureAdaptor;
  my $dbentry_a   = $self->db->get_DBEntryAdaptor;
  my $feattype_a  = $self->db->get_FeatureTypeAdaptor;
  my $gene_a      = $self->db->dnadb->get_GeneAdaptor;
  my $slice_a     = $self->db->dnadb->get_SliceAdaptor;
  
  

  my $log = {};

  my $fset_config   = $self->{static_config}{feature_sets}{'FANTOM'};
  my $fset          = $fset_config->{feature_set};
  
  $self->rollback_FeatureSet($fset);

  for my $file (@$files){
    $self->log_header("Parsing and loading FANTOM data from:\t$file");
    my $flag_track_line = 0;
    open(my $fh,'<',$file) or throw "Can't access $file";
      my $ft_name;
      my $type; 
      my $feature_type_cfg;

      while (my $line = <$fh>) {
        chomp($line);

        if($flag_track_line == 0){
          $line =~ /track name="(.*?)"/;
          my $track_name = $1;

          if(!$track_name){
            throw("Wrong format. First line is expected to contain 'track name=\"name\"'. Line:\n$line"); 
          }

          if($track_name =~ 'permissive_enhancers'){
            $type = 'enhancers';
            $feature_type_cfg = $self->{static_config}{feature_types}{FANTOM_enhancer_permissive};
            $ft_name = $feature_type_cfg->{name};

          }
          elsif($track_name eq 'robust_enhancers'){
            $type = 'enhancers';
            $feature_type_cfg = $self->{static_config}{feature_types}{FANTOM_enhancer_robust};
            $ft_name = $feature_type_cfg->{name};
          }
          elsif($track_name eq 'Timo TSS predictions'){
            $type = 'tss';
          }
          else{
            throw("Unkown track_name: [$track_name] in line [$line]");
          }
          $flag_track_line = 1;
          next;
        } # first line

      my $data;
      if($type eq 'tss'){
        $data = $dispatch->{tss}->($line);
          
        # relaxed
        if($data->{rgb} eq '30,144,255'){
          $ft_name = $self->{static_config}{feature_types}{FANTOM_tss_relaxed}->{name};
        }
        
        # strict
        elsif($data->{rgb} eq '60,179,113'){
          $ft_name = $self->{static_config}{feature_types}{FANTOM_tss_strict}->{name};
        }
        else{next}
      }
      elsif($type eq 'enhancers'){
        $data = $dispatch->{enhancers}->($line);
          # say dump_data($data,1,1);
      }
      else{
        throw("Unkown type: $type");
      }

      my $slice = $slice_a->fetch_by_region('chromosome', $data->{chr},);
      
      my $feature_type = shift @{$feattype_a->fetch_all_by_name($ft_name)};

      if(! defined($feature_type)){
        throw "No FeatureType for $ft_name";
      }
      
      my $feature = Bio::EnsEMBL::Funcgen::ExternalFeature->new (
       -start         => $data->{start}, 
       -end           => $data->{end},  
       -strand        => $data->{strand},
       -feature_type  => $feature_type,
       -slice         => $slice,
       -display_label => $data->{name}, 
       -feature_set   => $fset,
       );
      # say $line;
      # say $feature->slice->name;
      # say $feature->start;
      # say $feature->end;
      $feature = $self->project_feature($feature, $new_assembly);
      # say $feature->slice->name."\n";
      # say $feature->start;
      # say $feature->end;

      if($feature){
        if($type eq 'enhancers'){
          $feature->{display_label}  = 
          $feature->slice->seq_region_name . ':'.
          $feature->start . '-' .
          $feature->end;
        }
        ($feature) = @{$extfeat_a->store($feature)};
        foreach my $status (qw(DISPLAYABLE MART_DISPLAYABLE)) {
          $fset->adaptor->store_status($status, $fset);
        }
        # say "Stored: " .$feature->dbID;
 # die;

      }

      # my $dbentry = Bio::EnsEMBL::DBEntry->new(
      #  -primary_id             => $gene->stable_id,
      #  -dbname                 => $external_db_name,
      #  -release                => $external_db_release, 
      #  -display_id             => $gene->display_xref->display_id,
      #  -status                 => 'KNOWNXREF',
      #  -db_display_name        => 'EnsemblGene',
      #  -type                   => 'MISC', # this is external_db.type
      #  -info_type              => 'MISC',
      #  -info_text              => 'GENE',
      #  -linkage_annotation     => 'FANTOM ',
      #  # -linkage_annotation     => 'FANTOM '. $track_name ,
      #  -analysis               => $fset->analysis,
      # );
     # $dbentry_a->store($dbentry, $feature->dbID, 'ExternalFeature');#
     #Now set states
    }
    close $fh;
  }
  say dump_data($log,1,1);    

}


1;
