package Exon;
#exon object to store coordinates, omega average and omega per position

use strict;


sub new {
    my($class) = shift;

    my ($self) = { };
    bless $self;
    my ($start, $end, $strand, $om_average, $std, $p_val, $d_val, @omega,@ts) = @_;

    $self->{'start'} = $start;
    $self->{'end'} = $end;
    $self->{'strand'} = $strand;
    $self->{'om_average'} = $om_average;
    $self->{'std'} = $std;
    $self->{'d_value'}=$p_val; 
	$self->{'p_value'}=$d_val;
    $self->{'omega'}  = [];
    @{$self->{'omega'}} = @omega;
    $self ->{'ts'} = [];
    @{$self->{'ts'}} = @ts;

    return $self;
}

sub start{
    my $self = shift;
    if(@_){
	$self->{'start'} = @_;
    }
    return $self->{'start'};
}
sub end{
    my $self = shift;
    if(@_){
	$self->{'end'} = @_;
    }
    return $self->{'end'};
}
sub strand{
    my $self = shift;
    if(@_){
	$self->{'strand'} = @_;
    }
    return $self->{'strand'};
}
sub om_average{
    my $self = shift;
    if(@_){
	$self->{'om_average'} = @_;
    }
    return $self->{'om_average'};
}
sub std{
    my $self = shift;
    if(@_){
	$self->{'std'} = @_;
    }
    return $self->{'std'};
}
sub omega{
    my $self = shift;
    if(@_){
	@{$self->{'omega'}}  = @_;
    }
    return  $self->{'omega'};
}
sub ts{
    my $self = shift;
    if(@_){
	@{$self->{'ts'}}  = @_;
    }
    return  $self->{'ts'};
}

1;
