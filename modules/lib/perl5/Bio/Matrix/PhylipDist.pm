# BioPerl module for Bio::Matrix::PhylipDist
#
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Shawn Hoon <shawnh@fugu-sg.org>
#
# Copyright Shawn Hoon
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Matrix::PhylipDist - A Phylip Distance Matrix object 

=head1 SYNOPSIS

  use Bio::Tools::Phylo::Phylip::ProtDist;
  my $dist = Bio::Tools::Phylo::Phylip::ProtDist->new(
    -file=>"protdist.out",
    -program=>"ProtDist");
  #or
   my $dist = Bio::Tools::Phylo::Phylip::ProtDist->new(
    -fh=>"protdist.out",
    -program=>"ProtDist");


  #get specific entries
  my $distance_value = $dist->get_entry('ALPHA','BETA');
  my @columns        = $dist->get_column('ALPHA');
  my @rows           = $dist->get_row('BETA');
  my @diagonal       = $dist->get_diagonal();

  #print the matrix in phylip numerical format
  print $dist->print_matrix;

=head1 DESCRIPTION

Simple object for holding Distance Matrices generated by the following Phylip programs:

1) dnadist
2) protdist
3) restdist

It currently handles parsing of the matrix without the data output option.

    5
Alpha          0.00000  4.23419  3.63330  6.20865  3.45431
Beta           4.23419  0.00000  3.49289  3.36540  4.29179
Gamma          3.63330  3.49289  0.00000  3.68733  5.84929
Delta          6.20865  3.36540  3.68733  0.00000  4.43345
Epsilon        3.45431  4.29179  5.84929  4.43345  0.00000

=head1 FEEDBACK


=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists. Your participation is much appreciated.

  bioperl-l@bioperl.org                  - General discussion
  http://bioperl.org/wiki/Mailing_lists  - About the mailing lists

=head2 Support 

Please direct usage questions or support issues to the mailing list:

I<bioperl-l@bioperl.org>

rather than to the module maintainer directly. Many experienced and 
reponsive experts will be able look at the problem and quickly 
address it. Please include a thorough description of the problem 
with code and data examples if at all possible.

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution.  Bug reports can be submitted via the
web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR - Shawn Hoon

Email shawnh@fugu-sg.org

=head1 CONTRIBUTORS

Jason Stajich, jason-at-bioperl-dot-org

=head1 APPENDIX


The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a "_".

=cut

# Let the code begin...

package Bio::Matrix::PhylipDist;
$Bio::Matrix::PhylipDist::VERSION = '1.7.8';
use strict;


use base qw(Bio::Root::Root Bio::Matrix::MatrixI);

=head2 new

 Title   : new
 Usage   : my $family = Bio::Matrix::PhylipDist->new(-file=>"protdist.out",
                                                     -program=>"protdist");
 Function: Constructor for PhylipDist Object
 Returns : L<Bio::Matrix::PhylipDist>

=cut

sub new {
    my ($class,@args) = @_;
    my $self = $class->SUPER::new(@args);
    my ($matrix,$values, $names,
	$program,$matname,
	$matid) = $self->_rearrange([qw(MATRIX 
					VALUES 
					NAMES 
					PROGRAM
					MATRIX_NAME
					MATRIX_ID
					)],@args);
    
    ($matrix && $values && $names) || 
	$self->throw("Need matrix, values, and names fields all provided!");

    $program && $self->matrix_name($program) if defined $program;
    
    $self->_matrix($matrix) if ref($matrix) =~ /HASH/i;
    $self->_values($values) if ref($values) =~ /ARRAY/i;
    $self->names($names) if ref($names) =~ /ARRAY/i;

    $self->matrix_name($matname) if defined $matname;
    $self->matrix_id  ($matid)   if defined $matid;

    return $self;
}

=head2 get_entry

 Title   : get_entry
 Usage   : $matrix->get_entry();
 Function: returns a particular entry 
 Returns : a float
 Arguments:  string id1, string id2

=cut

sub get_entry {
  my ($self,$row,$column) = @_;
  $row && $column || $self->throw("Need at least 2 ids");
  my %matrix = %{$self->_matrix};
  my @values = @{$self->_values};
  if(ref $matrix{$row}{$column}){
      my ($i,$j) = @{$matrix{$row}{$column}};
      return $values[$i][$j];
  }
  return;

}

=head2 get_row

 Title   : get_row
 Usage   : $matrix->get_row('ALPHA');
 Function: returns a particular row 
 Returns : an array of float
 Arguments:  string id1

=cut

sub get_row {
    my ($self,$row) = @_;
    $row || $self->throw("Need at least a row id");

    my %matrix = %{$self->_matrix};
    my @values = @{$self->_values};
    my @names = @{$self->names};
    $matrix{$row} || return;
    my ($val) = values %{$matrix{$row}};
    my $row_pointer = $val->[0];
    my $index = scalar(@names)-1;
    return @{$values[$row_pointer]}[0..$index];
}

=head2 get_column

 Title   : get_column
 Usage   : $matrix->get_column('ALPHA');
 Function: returns a particular column 
 Returns : an array of floats 
 Arguments:  string id1

=cut

sub get_column {
    my ($self,$column) = @_;
    $column || $self->throw("Need at least a column id");

    my %matrix = %{$self->_matrix};
    my @values = @{$self->_values};
    my @names = @{$self->names}; 
    $matrix{$column} || return ();
    my ($val) = values %{$matrix{$column}};
    my $row_pointer = $val->[0];
    my @ret;
    for(my $i=0; $i < scalar(@names); $i++) {
	push @ret, $values[$i][$row_pointer];
    }
    return @ret;
} 

=head2 get_diagonal

 Title   : get_diagonal
 Usage   : $matrix->get_diagonal();
 Function: returns the diagonal of the matrix
 Returns : an array of float
 Arguments:  string id1

=cut

sub get_diagonal {
  my ($self) = @_;
  my %matrix = %{$self->_matrix};
  my @values = @{$self->_values};
  my @return;
  foreach my $name (@{$self->names}){
    my ($i,$j) = @{$matrix{$name}{$name}};
    push @return,$values[$i][$j];
  }
  return @return;
}

=head2 print_matrix

 Title   : print_matrix
 Usage   : $matrix->print_matrix();
 Function: returns a string of the matrix in phylip format 
 Returns : a string
 Arguments:  

=cut

sub print_matrix {
  my ($self) = @_;
  my @names = @{$self->names};
  my @values = @{$self->_values};
  my %matrix = %{$self->_matrix};
  my $str;
  $str.= (" "x 4). scalar(@names)."\n";
  foreach my $name (@names){
    my $newname = $name. (" " x (15-length($name)));
    if( length($name) >= 15 ) { $newname .= " " }
    $str.=$newname;
    my $count = 0;
    foreach my $n (@names) {
      my ($i,$j) = @{$matrix{$name}{$n}};
      if($count < $#names){
        $str .= $values[$i][$j]. "  ";
      }
      else {
	  if( ! defined $values[$i][$j] ) { 
	      $self->debug("no value for $i,$j cell\n");
	  } else { 
	      $str .= $values[$i][$j];
	  }
      }
      $count++;
    }
    $str.="\n";
  }
  return $str;
}

=head2 _matrix

 Title   : _matrix
 Usage   : $matrix->_matrix();
 Function: get/set for hash reference of the pointers
           to the value matrix 
 Returns : hash reference 
 Arguments: hash reference

=cut

sub _matrix {
  my ($self,$val) = @_;
  if($val){
    $self->{'_matrix'} = $val;
  }
  return $self->{'_matrix'};
}


=head2 names

 Title   : names
 Usage   : $matrix->names();
 Function: get/set for array ref of names of sequences
 Returns : an array reference 
 Arguments: an array reference

=cut

sub names {
  my ($self,$val) = @_;
  if($val){
    $self->{'_names'} = $val;
  }
  return $self->{'_names'};
}

=head2 program

 Title   : program
 Usage   : $matrix->program();
 Function: get/set for the program name generating this 
           matrix
 Returns : string
 Arguments: string

=cut

sub program {
  my ($self) = shift;
  return $self->matrix_name(@_);
}

=head2 _values

 Title   : _values
 Usage   : $matrix->_values();
 Function: get/set for array ref of the matrix containing
           distance values 
 Returns : an array reference 
 Arguments: an array reference

=cut

sub _values {
  my ($self,$val) = @_;
  if($val){
    $self->{'_values'} = $val;
  }
  return $self->{'_values'};
}


=head1 L<Bio::Matrix::MatrixI> implementation


=head2 matrix_id

 Title   : matrix_id
 Usage   : my $id = $matrix->matrix_id
 Function: Get/Set the matrix ID
 Returns : scalar value
 Args    : [optional] new id value to store


=cut

sub matrix_id{
   my $self = shift;
   return $self->{'_matid'} = shift if @_;
   return $self->{'_matid'};

   
}

=head2 matrix_name

 Title   : matrix_name
 Usage   : my $name = $matrix->matrix_name();
 Function: Get/Set the matrix name
 Returns : scalar value
 Args    : [optional] new matrix name value


=cut

sub matrix_name{
   my $self = shift;
   return $self->{'_matname'} = shift if @_;
   return $self->{'_matname'};
}

=head2 column_header

 Title   : column_header
 Usage   : my $name = $matrix->column_header(0)
 Function: Gets the column header for a particular column number
 Returns : string
 Args    : integer


=cut

sub column_header{
    my ($self,$num) = @_;
    my @coln = $self->column_names;
    return $coln[$num];
}


=head2 row_header

 Title   : row_header
 Usage   : my $name = $matrix->row_header(0)
 Function: Gets the row header for a particular row number
 Returns : string
 Args    : integer


=cut

sub row_header{
    my ($self,$num) = @_;
    my @rown = $self->row_names;
   return $rown[$num];
}
=head2 column_num_for_name

 Title   : column_num_for_name
 Usage   : my $num = $matrix->column_num_for_name($name)
 Function: Gets the column number for a particular column name
 Returns : integer
 Args    : string


=cut

sub column_num_for_name{
   my ($self,$name) = @_;
   my $ct = 0;
   foreach my $n ( $self->column_names ) {
       return $ct if $n eq $name;
       $ct++;
   }
   return;
}

=head2 row_num_for_name

 Title   : row_num_for_name
 Usage   : my $num = $matrix->row_num_for_name($name)
 Function: Gets the row number for a particular row name
 Returns : integer
 Args    : string


=cut

sub row_num_for_name{
   my ($self,$name) = @_;
   my $ct = 0;
   foreach my $n ( $self->row_names ) {
       return $ct if $n eq $name;
       $ct++;
   }
}

=head2 num_rows

 Title   : num_rows
 Usage   : my $rowcount = $matrix->num_rows;
 Function: Get the number of rows
 Returns : integer
 Args    : none


=cut

sub num_rows{ return scalar @{shift->names} }

=head2 num_columns

 Title   : num_columns
 Usage   : my $colcount = $matrix->num_columns
 Function: Get the number of columns
 Returns : integer
 Args    : none


=cut

sub num_columns{
   return scalar @{shift->names};
}

=head2 row_names

 Title   : row_names
 Usage   : my @rows = $matrix->row_names
 Function: The names of all the rows
 Returns : array in array context, arrayref in scalar context
 Args    : none


=cut

sub row_names{ return @{shift->names} }

=head2 column_names

 Title   : column_names
 Usage   : my @columns = $matrix->column_names
 Function: The names of all the columns
 Returns : array in array context, arrayref in scalar context
 Args    : none


=cut

sub column_names{ return @{shift->names} }  
1;
