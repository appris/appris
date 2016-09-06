=head1 CONTACT

  Please email comments or questions to the public INB
  developers list at <inb-tecnico@lists.cnio.es>.

  Questions may also be sent to the developer, 
  Jose Manuel Rodriguez <jmrodriguez@cnio.es>.

=cut

=head1 NAME

APPRIS::DBSQL::DBAdaptor

=head1 SYNOPSIS

  $db = APPRIS::DBSQL::DBAdaptor->new(
    -user   => 'root',
    -dbname => 'pog',
    -host   => 'caldy',
    -driver => 'mysql'
  );

  $gene = $gene_adaptor->fetch_by_stable_id($stable_id);

=head1 DESCRIPTION

Formerly this class provided database connectivity and a means
to retrieve object adaptors.

=head1 METHODS

=cut

package APPRIS::DBSQL::DBAdaptor;

use DBI;
use DBD::mysql;
use Data::Dumper;

{
	#Encapsulated class data

	#___________________________________________________________
	#ATTRIBUTES
	my %_attr_data = #          DEFAULT         ACCESSIBILITY
		(
			driver			=>  ["DBI:mysql",	'read/write'],
			dbh				=>  [undef,			'read/write'],
			dbhost			=>  ['',			undef],
			dbport			=>  ['3306',		undef],
			dbname			=>  ['',			undef],
			dbuser			=>  ['',			undef],
			dbpass			=>  ['',			undef],
		);

	#_____________________________________________________________

	# Classwide default value for a specified object attribute
	sub _default_for {
		my ($self, $attr) = @_;
		$_attr_data{$attr}[0];
	}

	# List of names of all specified object attributes
	sub _standard_keys {
		keys %_attr_data;
	}

	sub driver {
		my ($self, $arg) = @_;
		$self->{driver} = $arg if defined $arg;
		return $self->{driver};
	}
	
	sub dbh {
		my ($self, $arg) = @_;
		$self->{dbh} = $arg if defined $arg;
		return $self->{dbh};
	}

	sub dbhost {
		my ($self, $arg) = @_;
		$self->{dbhost} = $arg if defined $arg;
		return $self->{dbhost};
	}

	sub dbport {
		my ($self, $arg) = @_;
		$self->{dbport} = $arg if defined $arg;
		return $self->{dbport};
	}

	sub dbname {
		my ($self, $arg) = @_;
		$self->{dbname} = $arg if defined $arg;
		return $self->{dbname};
	}

	sub dbuser {
		my ($self, $arg) = @_;
		$self->{dbuser} = $arg if defined $arg;
		return $self->{dbuser};
	}

	sub dbpass {
		my ($self, $arg) = @_;
		$self->{dbpass} = $arg if defined $arg;
		return $self->{dbpass};
	}
	
}

sub new {
	my ($caller, %args) = @_;
	my $caller_is_obj = ref($caller);
	return $caller if $caller_is_obj;
	my $class = $caller_is_obj || $caller;
	my $proxy;
	my ($self) = bless {}, $class;


	foreach my $attrname ( $self->_standard_keys ) {

		if (exists $args{$attrname} && defined $args{$attrname}) {
			$self->{$attrname} = $args{$attrname};
		} else {
			$self->{$attrname} = $self->_default_for($attrname);
		}
	}

	return unless $self->driver;
	my $driver = $self->driver;
	my $dbname = $self->dbname;
	my $dbhost = $self->dbhost;
	my $dbport = $self->dbport;
	my $dbuser = $self->dbuser;
	my $dbpass = $self->dbpass;

	my $dsn = "$driver:dbname=$dbname;host=$dbhost;port=$dbport";
	my $dbh = DBI->connect($dsn, $dbuser, $dbpass, {RaiseError => 1, AutoCommit => 0}) or die "\ncan't connect to database: " . $DBI::errstr . "\n\n";
   	$dbh->{'mysql_enable_utf8'} = 1; #UTF8 enable

	##############################################################
	unless ($dbh) {
		print STDERR "Couldn't connect to the datasource: " . $DBI::errstr ."\n\n";
		return undef;
	}

	$self->dbh($dbh);
	#############################################################

	return undef unless $self->dbh;
	return $self;
}


###################
# PRIVATE methods #
###################
sub _do_query{
	my ($dbh, $statement, @bindvalues) = @_;

	my $sth = $dbh -> prepare($statement);
	if (@bindvalues < 1)
	{
		$sth->execute;
	}
	else
	{
		$sth->execute(@bindvalues);
	}
	# returns an array of hash references
	my $arrayHashRef = $sth->fetchall_arrayref({});
	return $arrayHashRef;
}

sub do_sql_sentence {
	my ($dbh, $statement, @bindvalues) = @_;

	my $sth = $dbh -> prepare($statement);
	if (@bindvalues < 1)
	{
		$sth->execute;
	}
	else
	{
		$sth->execute(@bindvalues);
	}
	if ($dbh->err){
		return 0;
	}
	else{
		$dbh->commit();
		return 1;
	}
}

sub _add_condition {
	my ($statement, @params) = @_;
	my @bindvalues = ();
	my $condition = '';
	$condition = 'where ' if(scalar(@params)>0);

	foreach my $param (@params ) {
		if (($param eq 'and') || ($param eq 'or')) {
			$condition .= $param . " ";
		} else {
			my %pair = %{$param};
			for my $key (keys %pair) {
				if (defined $pair{$key}) {
					$condition .= $key . " = ? ";
					push(@bindvalues, $pair{$key});
				} else {
					$condition .= $key . " IS NULL "
				}
			}
		}
	}

	$statement .= $condition;
	return ($statement, @bindvalues);
}

sub _add_condition2 {
	my ($statement, @params) = @_;
	my @bindvalues = ();
	my $condition = '';
	$condition = 'where ' if(scalar(@params)>0);

	for (my $i=0; $i<scalar(@params);$i++) {
		my $param = $params[$i];
		if (($param eq 'and') || ($param eq 'or')) {
			$condition .= $param . " ";
		} else {
			my %pair = %{$param};
			for my $key (keys %pair) {
				if (defined $pair{$key}) {
					if($params[$i+1] and ($params[$i+1] eq 'or')) {
						$condition .= " ( " . $key . " = ? ";	
					} elsif($params[$i-1] and ($params[$i-1] eq 'or')) {
						$condition .= $key . " = ? )";	
					} else {
						$condition .= $key . " = ? ";	
					}
					push(@bindvalues, $pair{$key});
				} else {
					$condition .= $key . " IS NULL "
				}
			}
		}
	}

	$statement .= $condition;
	return ($statement, @bindvalues);
}

sub _add_condition_like {
	my ($statement, @params) = @_;
	my @bindvalues = ();
	my $condition = '';
	$condition = 'where ' if(scalar(@params)>0);

	foreach my $param (@params ) {
		if (($param eq 'and') || ($param eq 'or')) {
			$condition .= "%" . $param . "% ";
		} else {
			my %pair = %{$param};
			for my $key (keys %pair) {
				if (defined $pair{$key}) {
					$condition .= $key . " like ? ";
					push(@bindvalues, "%" . $pair{$key} . "%");
				} else {
					$condition .= $key . " IS NULL "
				}
			}
		}
	}
	$condition .= 'order by ' if(scalar(@params)>0);
	foreach my $param (@params ) {
		unless (($param eq 'and') and ($param eq 'or')) {
			my %pair = %{$param};
			for my $key (keys %pair) {
				if (defined $pair{$key}) {
					$condition .= $key . ',';
				}
			}
		}
	}
	$condition =~ s/\,$//;

	$statement .= $condition;
	return ($statement, @bindvalues);
}

sub _add_condition_combine2 {
	my ($statement, @params) = @_;
	my @bindvalues = ();
	my $condition = '';
	$condition = 'where x.entity_id=y.entity_id and ' if(scalar(@params)>0);

	foreach my $param (@params ) {
		if (($param eq 'and') || ($param eq 'or')) {
			$condition .= $param . " ";
		} else {
			my %pair = %{$param};
			for my $key (keys %pair) {
				if (defined $pair{$key}) {
					$condition .= $key . " = ? ";
					push(@bindvalues, $pair{$key});
				} else {
					$condition .= $key . " IS NULL "
				}
			}
		}
	}

	$statement .= $condition;
	return ($statement, @bindvalues);
}

sub _add_condition_combine3 {
	my ($statement, @params) = @_;
	my @bindvalues = ();
	my $condition = '';
	$condition = 'where f.firestar_id=r.firestar_id and ' if(scalar(@params)>0);

	my $j = 0;
	for (my $i=0; $i<scalar(@params);$i++) {
		my $param = $params[$i];
		if (($param eq 'and') || ($param eq 'or')) {
			$condition .= $param . " ";
		} else {
			if(ref $param eq 'ARRAY') {
				$j=$j+1;
				if ($j == 1) {
					$condition .= " ( ";	
				}
				foreach my $v1 (@{$param}) {
					my %pair = %{$v1};
					for my $key (keys %pair) {
						if (defined $pair{$key}) {
							if ($j == 1) {
								if($key eq 'trans_start') {
									$condition .= $key . " >= ? ) ";	
								}elsif($key eq 'trans_end') {
									$condition .= " ( ". $key . " <= ? and ";	
								}								
							} else {
								if($key eq 'trans_start') {
									$condition .= $key . " <= ? ) ";	
								}elsif($key eq 'trans_end') {
									$condition .= " ( ". $key . " >= ? and ";	
								}									
							}
							push(@bindvalues, $pair{$key});
						}
					}
				}
				if ($j == 2) {
					$condition .= " ) ";	
				}
			} else {
				my %pair = %{$param};
				for my $key (keys %pair) {
					if (defined $pair{$key}) {
						$condition .= $key . " = ? ";
						push(@bindvalues, $pair{$key});
					} else {
						$condition .= $key . " IS NULL "
					}
				}
			}
		}
	}
	
	$statement .= $condition;
	return ($statement, @bindvalues);
}

sub _add_condition_combine4 {
	my ($statement, @params) = @_;
	my @bindvalues = ();
	my $condition = '';
	$condition = 'where i.inertia_id=o.inertia_id and o.omega_id=r.omega_id and ' if(scalar(@params)>0);

	foreach my $param (@params ) {
		if (($param eq 'and') || ($param eq 'or')) {
			$condition .= $param . " ";
		} else {
			my %pair = %{$param};
			for my $key (keys %pair) {
				if (defined $pair{$key}) {
					$condition .= $key . " = ? ";
					push(@bindvalues, $pair{$key});
				} else {
					$condition .= $key . " IS NULL "
				}
			}
		}
	}

	$statement .= $condition;
	return ($statement, @bindvalues);
}

sub _add_condition_combine5 {
	my ($statement, @params) = @_;
	my @bindvalues = ();
	my $condition = '';
	$condition = 'where m.matador3d_id=a.matador3d_id and ' if(scalar(@params)>0);

	my $j = 0;
	for (my $i=0; $i<scalar(@params);$i++) {
		my $param = $params[$i];
		if (($param eq 'and') || ($param eq 'or')) {
			$condition .= $param . " ";
		} else {
			if(ref $param eq 'ARRAY') {
				$j=$j+1;
				if ($j == 1) {
					$condition .= " ( ";	
				}
				foreach my $v1 (@{$param}) {
					my %pair = %{$v1};
					for my $key (keys %pair) {
						if (defined $pair{$key}) {
							if ($j == 1) {
								if($key eq 'trans_start') {
									$condition .= $key . " >= ? ) ";	
								}elsif($key eq 'trans_end') {
									$condition .= " ( ". $key . " <= ? and ";	
								}								
							} else {
								if($key eq 'trans_start') {
									$condition .= $key . " <= ? ) ";	
								}elsif($key eq 'trans_end') {
									$condition .= " ( ". $key . " >= ? and ";	
								}									
							}
							push(@bindvalues, $pair{$key});
						}
					}
				}
				if ($j == 2) {
					$condition .= " ) ";	
				}
			} else {
				my %pair = %{$param};
				for my $key (keys %pair) {
					if (defined $pair{$key}) {
						$condition .= $key . " = ? ";
						push(@bindvalues, $pair{$key});
					} else {
						$condition .= $key . " IS NULL "
					}
				}
			}
		}
	}

	$statement .= $condition;
	return ($statement, @bindvalues);
}

sub _add_condition_combine6 {
	my ($statement, @params) = @_;
	my @bindvalues = ();
	my $condition = '';
	$condition = 'where m.matador3d_id=a.matador3d_id and a.matador3d_alignments_id=p.matador3d_alignments_id and ' if(scalar(@params)>0);

	foreach my $param (@params ) {
		if (($param eq 'and') || ($param eq 'or')) {
			$condition .= $param . " ";
		} else {
			my %pair = %{$param};
			for my $key (keys %pair) {
				if (defined $pair{$key}) {
					$condition .= $key . " = ? ";
					push(@bindvalues, $pair{$key});
				} else {
					$condition .= $key . " IS NULL "
				}
			}
		}
	}

	$statement .= $condition;
	return ($statement, @bindvalues);
}

sub _add_condition_combine8 {
	my ($statement, @params) = @_;
	my @bindvalues = ();
	my $condition = '';
	$condition = 'where i.inertia_id=o.inertia_id and ' if(scalar(@params)>0);


	foreach my $param (@params ) {
		if (($param eq 'and') || ($param eq 'or')) {
			$condition .= $param . " ";
		} else {
			my %pair = %{$param};
			for my $key (keys %pair) {
				if (defined $pair{$key}) {
					$condition .= $key . " = ? ";
					push(@bindvalues, $pair{$key});
				} else {
					$condition .= $key . " IS NULL "
				}
			}
		}
	}

	$statement .= $condition;
	return ($statement, @bindvalues);
}

sub _add_condition_combine9 {
	my ($statement, @params) = @_;
	my @bindvalues = ();
	my $condition = '';
	$condition = 'where x.entity_id=y.entity_id and ' if(scalar(@params)>0);

	for (my $i=0; $i<scalar(@params);$i++) {
		my $param = $params[$i];
		if (($param eq 'and') || ($param eq 'or')) {
			$condition .= $param . " ";
		} else {
			my %pair = %{$param};
			for my $key (keys %pair) {
				if (defined $pair{$key}) {
					if(($params[$i+1] and $params[$i+1] eq 'or')) {
						$condition .= " ( " . $key . " = ? ";	
					} elsif($params[$i-1] and ($params[$i-1] eq 'or')) {
						$condition .= $key . " = ? )";	
					} else {
						$condition .= $key . " = ? ";	
					}
					push(@bindvalues, $pair{$key});
				} else {
					$condition .= $key . " IS NULL "
				}
			}
		}
	}

	$statement .= $condition;
	return ($statement, @bindvalues);
}

sub _add_condition_combine10 {
	my ($statement, @params) = @_;
	my @bindvalues = ();
	my $condition = '';
	$condition = 'where s.spade_id=a.spade_id and ' if(scalar(@params)>0);

	my $j = 0;
	for (my $i=0; $i<scalar(@params);$i++) {
		my $param = $params[$i];
		if (($param eq 'and') || ($param eq 'or')) {
			$condition .= $param . " ";
		} else {
			if(ref $param eq 'ARRAY') {
				$j=$j+1;
				if ($j == 1) {
					$condition .= " ( ";	
				}
				foreach my $v1 (@{$param}) {
					my %pair = %{$v1};
					for my $key (keys %pair) {
						if (defined $pair{$key}) {
							if ($j == 1) {
								if($key eq 'trans_start') {
									$condition .= $key . " >= ? ) ";	
								}elsif($key eq 'trans_end') {
									$condition .= " ( ". $key . " <= ? and ";	
								}								
							} else {
								if($key eq 'trans_start') {
									$condition .= $key . " <= ? ) ";	
								}elsif($key eq 'trans_end') {
									$condition .= " ( ". $key . " >= ? and ";	
								}									
							}
							push(@bindvalues, $pair{$key});
						}
					}
				}
				if ($j == 2) {
					$condition .= " ) ";	
				}
			} else {
				my %pair = %{$param};
				for my $key (keys %pair) {
					if (defined $pair{$key}) {
						$condition .= $key . " = ? ";
						push(@bindvalues, $pair{$key});
					} else {
						$condition .= $key . " IS NULL "
					}
				}
			}
		}
	}

	$statement .= $condition;
	return ($statement, @bindvalues);
}

sub _add_condition_combine11 {
	my ($statement, @params) = @_;
	my @bindvalues = ();
	my $condition = '';
	$condition = 'where t.thump_id=h.thump_id and ' if(scalar(@params)>0);

	my $j = 0;
	for (my $i=0; $i<scalar(@params);$i++) {
		my $param = $params[$i];
		if (($param eq 'and') || ($param eq 'or')) {
			$condition .= $param . " ";
		} else {
			if(ref $param eq 'ARRAY') {
				$j=$j+1;
				if ($j == 1) {
					$condition .= " ( ";	
				}
				foreach my $v1 (@{$param}) {
					my %pair = %{$v1};
					for my $key (keys %pair) {
						if (defined $pair{$key}) {
							if ($j == 1) {
								if($key eq 'trans_start') {
									$condition .= $key . " >= ? ) ";	
								}elsif($key eq 'trans_end') {
									$condition .= " ( ". $key . " <= ? and ";	
								}								
							} else {
								if($key eq 'trans_start') {
									$condition .= $key . " <= ? ) ";	
								}elsif($key eq 'trans_end') {
									$condition .= " ( ". $key . " >= ? and ";	
								}									
							}
							push(@bindvalues, $pair{$key});
						}
					}
				}
				if ($j == 2) {
					$condition .= " ) ";	
				}
			} else {
				my %pair = %{$param};
				for my $key (keys %pair) {
					if (defined $pair{$key}) {
						$condition .= $key . " = ? ";
						push(@bindvalues, $pair{$key});
					} else {
						$condition .= $key . " IS NULL "
					}
				}
			}
		}
	}

	$statement .= $condition;
	return ($statement, @bindvalues);
}

sub _add_condition_combine12 {
	my ($statement, @params) = @_;
	my @bindvalues = ();
	my $condition = '';
	$condition = 'where x.entity_id=y.entity_id and ' if(scalar(@params)>0);

	for (my $i=0; $i<scalar(@params);$i++) {
		my $param = $params[$i];
		if (($param eq 'and') || ($param eq 'or')) {
			$condition .= $param . " ";
		} else {
			my %pair = %{$param};
			for my $key (keys %pair) {
				if (defined $pair{$key}) {
					if(($params[$i+1] and $params[$i+1] eq 'or')) {
						if($key eq 'start') {
							$condition .= " ( " . $key . " >= ? ";
						}elsif($key eq 'end') {
							$condition .= " ( " . $key . " <= ? ";
						} else {
							$condition .= " ( " . $key . " = ? ";
						}	
					} elsif($params[$i-1] and ($params[$i-1] eq 'or')) {
						if($key eq 'start') {
							$condition .= $key . " >= ? )";
						}elsif($key eq 'end') {
							$condition .= $key . " <= ? )";	
						} else {
							$condition .= $key . " = ? )";	
						}	
					} else {
						if($key eq 'start') {
							$condition .= $key . " >= ? ";	
						}elsif($key eq 'end') {
							$condition .= $key . " <= ? ";	
						} else {
							$condition .= $key . " = ? ";	
						}	
					}
					push(@bindvalues, $pair{$key});
				} else {
					$condition .= $key . " IS NULL "
				}
			}
		}
	}

	$statement .= $condition;
	return ($statement, @bindvalues);
}

sub _add_condition_combine15 {
	my ($statement, @params) = @_;
	my @bindvalues = ();
	my $condition = '';
	$condition = 'where i.inertia_id=r.inertia_id and ' if(scalar(@params)>0);

	foreach my $param (@params ) {
		if (($param eq 'and') || ($param eq 'or')) {
			$condition .= $param . " ";
		} else {
			my %pair = %{$param};
			for my $key (keys %pair) {
				if (defined $pair{$key}) {
					$condition .= $key . " = ? ";
					push(@bindvalues, $pair{$key});
				} else {
					$condition .= $key . " IS NULL "
				}
			}
		}
	}

	$statement .= $condition;
	return ($statement, @bindvalues);
}

sub _add_condition_combine16 {
	my ($statement, @params) = @_;
	my @bindvalues = ();
	my $condition = '';
	$condition = 'where x.entity_id=y.entity_id and ' if(scalar(@params)>0);

	for (my $i=0; $i<scalar(@params);$i++) {
		my $param = $params[$i];
		if (($param eq 'and') || ($param eq 'or')) {
			$condition .= $param . " ";
		} else {
			my %pair = %{$param};
			for my $key (keys %pair) {
				if (defined $pair{$key}) {
					if(($params[$i+1] and $params[$i+1] eq 'or')) {
						if($key eq 'start') {
							$condition .= " ( " . $key . " <= ? ";
						}elsif($key eq 'end') {
							$condition .= " ( " . $key . " >= ? ";
						} else {
							$condition .= " ( " . $key . " = ? ";
						}	
					} elsif($params[$i-1] and ($params[$i-1] eq 'or')) {
						if($key eq 'start') {
							$condition .= $key . " <= ? )";
						}elsif($key eq 'end') {
							$condition .= $key . " >= ? )";	
						} else {
							$condition .= $key . " = ? )";	
						}	
					} else {
						if($key eq 'start') {
							$condition .= $key . " <= ? ";	
						}elsif($key eq 'end') {
							$condition .= $key . " >= ? ";	
						} else {
							$condition .= $key . " = ? ";	
						}	
					}
					push(@bindvalues, $pair{$key});
				} else {
					$condition .= $key . " IS NULL "
				}
			}
		}
	}

	$statement .= $condition;
	return ($statement, @bindvalues);
}

sub _add_condition_combine17 {
	my ($statement, @params) = @_;
	my @bindvalues = ();
	my $condition = '';
	$condition = 'where c.crash_id=r.crash_id and ' if(scalar(@params)>0);

	my $j = 0;
	for (my $i=0; $i<scalar(@params);$i++) {
		my $param = $params[$i];
		if (($param eq 'and') || ($param eq 'or')) {
			$condition .= $param . " ";
		} else {
			if(ref $param eq 'ARRAY') {
				$j=$j+1;
				if ($j == 1) {
					$condition .= " ( ";	
				}
				foreach my $v1 (@{$param}) {
					my %pair = %{$v1};
					for my $key (keys %pair) {
						if (defined $pair{$key}) {
							if ($j == 1) {
								if($key eq 'trans_start') {
									$condition .= $key . " >= ? ) ";	
								}elsif($key eq 'trans_end') {
									$condition .= " ( ". $key . " <= ? and ";	
								}								
							} else {
								if($key eq 'trans_start') {
									$condition .= $key . " <= ? ) ";	
								}elsif($key eq 'trans_end') {
									$condition .= " ( ". $key . " >= ? and ";	
								}									
							}
							push(@bindvalues, $pair{$key});
						}
					}
				}
				if ($j == 2) {
					$condition .= " ) ";	
				}
			} else {
				my %pair = %{$param};
				for my $key (keys %pair) {
					if (defined $pair{$key}) {
						$condition .= $key . " = ? ";
						push(@bindvalues, $pair{$key});
					} else {
						$condition .= $key . " IS NULL "
					}
				}
			}
		}
	}

	$statement .= $condition;
	return ($statement, @bindvalues);
}

sub _add_condition_combine18 {
	my ($statement, @params) = @_;
	my @bindvalues = ();
	my $condition = '';
	$condition = 'where x.entity_id=y.entity_id and ' if(scalar(@params)>0);

	my $j = 0;
	for (my $i=0; $i<scalar(@params);$i++) {
		my $param = $params[$i];
		if (($param eq 'and') || ($param eq 'or')) {
			$condition .= $param . " ";
		} else {
			if(ref $param eq 'ARRAY') {
				$j=$j+1;
				if ($j == 1) {
					$condition .= " ( ";	
				}
				foreach my $v1 (@{$param}) {
					my %pair = %{$v1};
					for my $key (keys %pair) {
						if (defined $pair{$key}) {
							if ($j == 1) {
								if($key eq 'start') {
									$condition .= $key . " >= ? ) ";	
								}elsif($key eq 'end') {
									$condition .= " ( ". $key . " <= ? and ";	
								}								
							} else {
								if($key eq 'start') {
									$condition .= $key . " <= ? ) ";	
								}elsif($key eq 'end') {
									$condition .= " ( ". $key . " >= ? and ";	
								}									
							}
							push(@bindvalues, $pair{$key});
						}
					}
				}
				if ($j == 2) {
					$condition .= " ) ";	
				}
			} else {
				my %pair = %{$param};
				for my $key (keys %pair) {
					if (defined $pair{$key}) {
						$condition .= $key . " = ? ";
						push(@bindvalues, $pair{$key});
					} else {
						$condition .= $key . " IS NULL "
					}
				}
			}
		}
	}

	$statement .= $condition;	
	return ($statement, @bindvalues);
}

sub _add_condition_combine19 {
	my ($statement, @params) = @_;
	my @bindvalues = ();
	my $condition = '';
	$condition = 'where c.proteo_id=r.proteo_id and ' if(scalar(@params)>0);

	my $j = 0;
	for (my $i=0; $i<scalar(@params);$i++) {
		my $param = $params[$i];
		if (($param eq 'and') || ($param eq 'or')) {
			$condition .= $param . " ";
		} else {
			if(ref $param eq 'ARRAY') {
				$j=$j+1;
				if ($j == 1) {
					$condition .= " ( ";	
				}
				foreach my $v1 (@{$param}) {
					my %pair = %{$v1};
					for my $key (keys %pair) {
						if (defined $pair{$key}) {
							if ($j == 1) {
								if($key eq 'trans_start') {
									$condition .= $key . " >= ? ) ";	
								}elsif($key eq 'trans_end') {
									$condition .= " ( ". $key . " <= ? and ";	
								}								
							} else {
								if($key eq 'trans_start') {
									$condition .= $key . " <= ? ) ";	
								}elsif($key eq 'trans_end') {
									$condition .= " ( ". $key . " >= ? and ";	
								}									
							}
							push(@bindvalues, $pair{$key});
						}
					}
				}
				if ($j == 2) {
					$condition .= " ) ";	
				}
			} else {
				my %pair = %{$param};
				for my $key (keys %pair) {
					if (defined $pair{$key}) {
						$condition .= $key . " = ? ";
						push(@bindvalues, $pair{$key});
					} else {
						$condition .= $key . " IS NULL "
					}
				}
			}
		}
	}

	$statement .= $condition;
	return ($statement, @bindvalues);
}

sub _add_condition_update {
	my ($statement, @params) = @_;
	my @bindvalues = ();
	my $condition = '';

	foreach my $param (@params ) {
		if ($param eq ',') {
			$condition .= $param . " ";
		} else {
			my %pair = %{$param};
			for my $key (keys %pair) {
				if (defined $pair{$key}) {
					$condition .= $key . " = ? ";
					push(@bindvalues, $pair{$key});
				} else {
					$condition .= $key . " IS NULL "
				}
			}
		}
	}
	
	$statement .= $condition;
	return ($statement, @bindvalues);
}

##################
# SELECT methods #
##################
sub query_entity {
	my ($self, %args) = @_;
	my $dbh = $self->dbh;

	my @args;
	while (my ($k, $v) = each %args){
		push @args, ({$k => $v}, "and"); # format for the_add_condition subroutine but too bad won't be scalable for "or"
	}
	if (keys(%args)){ pop @args;}  # remove final "and"

        my $statement = "select
		entity_id,
		datasource_id,
		identifier,
		source,
		biotype,
		tsl,
		tag
  
		from entity ";
	my @bindvalues;
	($statement, @bindvalues) = _add_condition($statement, @args);

	my $final = _do_query($dbh, $statement, @bindvalues);
	return $final;
}

sub query_entity_like {
	my ($self, %args) = @_;
	my $dbh = $self->dbh;

	my @args;
	while (my ($k, $v) = each %args){
		push @args, ({$k => $v}, "and"); # format for the_add_condition subroutine but too bad won't be scalable for "or"
	}
	if (keys(%args)){ pop @args;}  # remove final "and"

        my $statement = "select
		entity_id,
		datasource_id,
		identifier,
		source,
		biotype,
		tsl,
		tag

		from entity ";

	my @bindvalues;
	($statement, @bindvalues) = _add_condition_like($statement, @args);

	my $final = _do_query($dbh, $statement, @bindvalues);
	return $final;
}

sub query_datasource {
	my ($self, %args) = @_;
	my $dbh = $self->dbh;

	my @args;
	while (my ($k, $v) = each %args){
		push @args, ({$k => $v}, "and"); # format for the_add_condition subroutine but too bad won't be scalable for "or"
	}
	if (keys(%args)){ pop @args;}  # remove final "and"

        my $statement = "select
		datasource_id,
		name,
		description,
		url
  
		from datasource ";

	my @bindvalues;
	($statement, @bindvalues) = _add_condition($statement, @args);

	my $final = _do_query($dbh, $statement, @bindvalues);
	return $final;
}

sub query_codon {
	my ($self, %args) = @_;
	my $dbh = $self->dbh;

	my @args;
	while (my ($k, $v) = each %args){
		push @args, ({$k => $v}, "and"); # format for the_add_condition subroutine but too bad won't be scalable for "or"
	}
	if (keys(%args)){ pop @args;}  # remove final "and"

        my $statement = "select
		type,
		start,
		end,
		strand,
		phase
  
		from codon ";

	my @bindvalues;
	($statement, @bindvalues) = _add_condition($statement, @args);

	my $final = _do_query($dbh, $statement, @bindvalues);
	return $final;
}

sub query_xref_identify {
	my ($self, %args) = @_;
	my $dbh = $self->dbh;

	my @args;
	while (my ($k, $v) = each %args){
		push @args, ({$k => $v}, "and"); # format for the_add_condition subroutine but too bad won't be scalable for "or"
	}
	if (keys(%args)){ pop @args;}  # remove final "and"

        my $statement = "select
		entity_id,
		datasource_id,
		identifier

		from xref_identify ";

	my @bindvalues;
	($statement, @bindvalues) = _add_condition($statement, @args);

	my $final = _do_query($dbh, $statement, @bindvalues);
	return $final;
}

sub query_xref_identify2 {
	my ($self, %args) = @_;
	my $dbh = $self->dbh;

	my @args;
	while (my ($k, $v) = each %args){
		if(ref $v eq 'ARRAY') {
			foreach my $v1 (@{$v}) {
				push @args, ({$k => $v1}, "or"); # OR 
			}			
		} else {
			push @args, ({$k => $v}, "and"); # AND
		}
	}
	if (keys(%args)){ pop @args;}  # remove final "and"

        my $statement = "select
		entity_id,
		datasource_id,
		identifier

		from xref_identify ";

	my @bindvalues;
	($statement, @bindvalues) = _add_condition2($statement, @args);

	my $final = _do_query($dbh, $statement, @bindvalues);
	return $final;
}

sub query_xref_identify_like {
	my ($self, %args) = @_;
	my $dbh = $self->dbh;

	my @args;
	while (my ($k, $v) = each %args){
		push @args, ({$k => $v}, "and"); # format for the_add_condition subroutine but too bad won't be scalable for "or"
	}
	if (keys(%args)){ pop @args;}  # remove final "and"

        my $statement = "select
		entity_id,
		datasource_id,
		identifier

		from xref_identify ";

	my @bindvalues;
	($statement, @bindvalues) = _add_condition_like($statement, @args);

	my $final = _do_query($dbh, $statement, @bindvalues);
	return $final;
}

sub query_coordinate {
	my ($self, %args) = @_;
	my $dbh = $self->dbh;

	my @args;
	while (my ($k, $v) = each %args){
		push @args, ({$k => $v}, "and"); # format for the_add_condition subroutine but too bad won't be scalable for "or"
	}
	if (keys(%args)){ pop @args;}  # remove final "and"

        my $statement = "select
        entity_id,
		chromosome,
		start,
		end,
		strand
  
		from coordinate ";
		
	my @bindvalues;
	($statement, @bindvalues) = _add_condition($statement, @args);

	my $final = _do_query($dbh, $statement, @bindvalues);
	return $final;
}
  
sub query_entity_coordinate {
	my ($self, %args) = @_;
	my $dbh = $self->dbh;

	my @args;
	while (my ($k, $v) = each %args){
		push @args, ({$k => $v}, "and"); # format for the_add_condition subroutine but too bad won't be scalable for "or"
	}
	if (keys(%args)){ pop @args;}  # remove final "and"

        my $statement = "select
		x.entity_id as entity_id,
		x.datasource_id as datasource_id,
		x.identifier as identifier,
		x.source as source,
		x.biotype as biotype,
		x.tsl as tsl,
		x.tag as tag,
		y.chromosome as chromosome,
		y.start as start,
		y.end as end,
		y.strand as strand
		
  
		from entity as x, coordinate as y ";
	my @bindvalues;
	($statement, @bindvalues) = _add_condition_combine2($statement, @args);

	my $final = _do_query($dbh, $statement, @bindvalues);
	return $final;
}

sub query_entity_coordinate2 {
	my ($self, %args) = @_;
	my $dbh = $self->dbh;

	my @args;
	while (my ($k, $v) = each %args){
		if(ref $v eq 'ARRAY')
		{
			foreach my $v1 (@{$v})
			{
				push @args, ({$k => $v1}, "or"); # OR 
			}
			pop @args; # remove final "or" and add "AND"
			push @args, "and";
		}
		else
		{
			push @args, ({$k => $v}, "and"); # AND
		}
	}
	if (keys(%args)){ pop @args;}  # remove final "and"
		
        my $statement = "select
		x.entity_id as entity_id,
		x.datasource_id as datasource_id,
		x.identifier as identifier,
		x.source as source,
		x.biotype as biotype,
		x.tsl as tsl,
		x.tag as tag,
		y.chromosome as chromosome,
		y.start as start,
		y.end as end,
		y.strand as strand
		
  
		from entity as x, coordinate as y ";

	my @bindvalues;
	($statement, @bindvalues) = _add_condition_combine9($statement, @args);

	my $final = _do_query($dbh, $statement, @bindvalues);
	return $final;
}

sub query_entity_coordinate4 {
	my ($self, %args) = @_;
	my $dbh = $self->dbh;

	my @args;
	while (my ($k, $v) = each %args){
		if(ref $v eq 'ARRAY')
		{
			foreach my $v1 (@{$v})
			{
				push @args, ({$k => $v1}, "or"); # OR 
			}
			pop @args; # remove final "or" and add "AND"
			push @args, "and";
		}
		else
		{
			push @args, ({$k => $v}, "and"); # AND
		}
	}
	if (keys(%args)){ pop @args;}  # remove final "and"
		
        my $statement = "select
		x.entity_id as entity_id,
		x.datasource_id as datasource_id,
		x.identifier as identifier,
		x.source as source,
		x.biotype as biotype,
		x.tsl as tsl,
		x.tag as tag,
		y.chromosome as chromosome,
		y.start as start,
		y.end as end,
		y.strand as strand
		
  
		from entity as x, coordinate as y ";

	my @bindvalues;
	($statement, @bindvalues) = _add_condition_combine16($statement, @args);
	my $final = _do_query($dbh, $statement, @bindvalues);
	return $final;
}

sub query_entity_coordinate5 {
	my ($self, @iargs) = @_;
	my $dbh = $self->dbh;

	my @args;
	foreach my $iarg (@iargs){
		if(ref $iarg eq 'ARRAY') {
			foreach my $v1 (@{$iarg}) {
				my @args1;
				foreach my $k (sort keys %{$v1}){ # important sort the keys!!!
					my $v = $v1->{$k};
					push(@args1, {$k => $v});
					push(@args1, "and");
				}
				pop @args1; # remove final "or" and add "AND"
				push(@args, [@args1]);
				push(@args, "or");				
			}
		} else {
			while (my ($k, $v) = each %{$iarg}){
				push @args, ({$k => $v}, "and"); # AND
			}
		}
	}
	pop @args; # remove final "or" and add "AND"
	
        my $statement = "select
		x.entity_id as entity_id,
		x.datasource_id as datasource_id,
		x.identifier as identifier,
		x.source as source,
		x.biotype as biotype,
		x.tsl as tsl,
		x.tag as tag,
		y.chromosome as chromosome,
		y.start as start,
		y.end as end,
		y.strand as strand
		
  
		from entity as x, coordinate as y ";

	my @bindvalues;
	($statement, @bindvalues) = _add_condition_combine18($statement, @args);
		
	my $final = _do_query($dbh, $statement, @bindvalues);
	return $final;
}

sub query_exon {
	my ($self, %args) = @_;
	my $dbh = $self->dbh;

	my @args;
	while (my ($k, $v) = each %args){
		push @args, ({$k => $v}, "and"); # format for the_add_condition subroutine but too bad won't be scalable for "or"
	}
	if (keys(%args)){ pop @args;}  # remove final "and"

        my $statement = "select
		entity_id,
		exon_id,
		identifier,
		start,
		end,
		strand

		from exon ";

	my @bindvalues;
	($statement, @bindvalues) = _add_condition($statement, @args);

	my $final = _do_query($dbh, $statement, @bindvalues);
	return $final;
}

sub query_cds {
	my ($self, %args) = @_;
	my $dbh = $self->dbh;

	my @args;
	while (my ($k, $v) = each %args){
		push @args, ({$k => $v}, "and"); # format for the_add_condition subroutine but too bad won't be scalable for "or"
	}
	if (keys(%args)){ pop @args;}  # remove final "and"

        my $statement = "select
		entity_id,
		cds_id,
		start,
		end,
		strand,
		phase

		from cds ";

	my @bindvalues;
	($statement, @bindvalues) = _add_condition($statement, @args);

	my $final = _do_query($dbh, $statement, @bindvalues);
	return $final;
}

sub query_type {
	my ($self, %args) = @_;
	my $dbh = $self->dbh;

	my @args;
	while (my ($k, $v) = each %args){
		push @args, ({$k => $v}, "and"); # format for the_add_condition subroutine but too bad won't be scalable for "or"
	}
	if (keys(%args)){ pop @args;}  # remove final "and"

        my $statement = "select
		type_id,
		name,
		description
  
		from type ";

	my @bindvalues;
	($statement, @bindvalues) = _add_condition($statement, @args);

	my $final = _do_query($dbh, $statement, @bindvalues);
	return $final;
}

sub query_sequence {
	my ($self, %args) = @_;
	my $dbh = $self->dbh;

	my @args;
	while (my ($k, $v) = each %args){
		push @args, ({$k => $v}, "and"); # format for the_add_condition subroutine but too bad won't be scalable for "or"
	}
	if (keys(%args)){ pop @args;}  # remove final "and"

        my $statement = "select
        entity_id,
		length,
		sequence
  
		from sequence ";
		
	my @bindvalues;
	($statement, @bindvalues) = _add_condition($statement, @args);

	my $final = _do_query($dbh, $statement, @bindvalues);
	return $final;
}

sub query_alignment {
	my ($self, %args) = @_;
	my $dbh = $self->dbh;

	my @args;
	while (my ($k, $v) = each %args){
		push @args, ({$k => $v}, "and"); # format for the_add_condition subroutine but too bad won't be scalable for "or"
	}
	if (keys(%args)){ pop @args;}  # remove final "and"

        my $statement = "select
        entity_id,
        type_id,
		length,
		num_species,
		alignment
  
		from alignment ";

	my @bindvalues;
	($statement, @bindvalues) = _add_condition($statement, @args);

	my $final = _do_query($dbh, $statement, @bindvalues);
	return $final;
}

sub query_extended_sequence {
	my ($self, %args) = @_;
	my $dbh = $self->dbh;

	my @args;
	while (my ($k, $v) = each %args){
		push @args, ({$k => $v}, "and"); # format for the_add_condition subroutine but too bad won't be scalable for "or"
	}
	if (keys(%args)){ pop @args;}  # remove final "and"

        my $statement = "select
        entity_id,
		length,
		sequence
  
		from extended_sequence ";
		
	my @bindvalues;
	($statement, @bindvalues) = _add_condition($statement, @args);

	my $final = _do_query($dbh, $statement, @bindvalues);
	return $final;
}

sub query_firestar {
	my ($self, %args) = @_;
	my $dbh = $self->dbh;

	my @args;
	while (my ($k, $v) = each %args){
		push @args, ({$k => $v}, "and"); # format for the_add_condition subroutine but too bad won't be scalable for "or"
	}
	if (keys(%args)){ pop @args;}  # remove final "and"

        my $statement = "select
        firestar_id,
        entity_id,
		num_residues,
		result,
		functional_residue
  
		from firestar ";

	my @bindvalues;
	($statement, @bindvalues) = _add_condition($statement, @args);

	my $final = _do_query($dbh, $statement, @bindvalues);
	return $final;
}

sub query_firestar_residues {
	my ($self, %iargs) = @_;
	my $dbh = $self->dbh;
	
	my @args;
	foreach my $k (sort keys %iargs){ # important sort the keys!!!
		my $v = $iargs{$k};
		if ( $k eq 'region' ) {
			my @args1;
			push(@args1, {'trans_end' => $v->{'start'}});
			push(@args1, "and");			
			push(@args1, {'trans_start' => $v->{'end'}});
			push(@args, [@args1]);
			push(@args, "or");
			my @args2;
			push(@args2, {'trans_end' => $v->{'start'}});
			push(@args2, "and");			
			push(@args2, {'trans_start' => $v->{'end'}});
			push(@args, [@args2]);
			push(@args, "or");
		}
		else {
			push @args, ({$k => $v}, "and"); # format for the_add_condition subroutine but too bad won't be scalable for "or"
		}
	}
	if (keys(%iargs)){ pop @args;}  # remove final "and"
		
        my $statement = "select
        f.firestar_id,
        f.num_residues,
        f.functional_residue,
        r.firestar_residues_id,
		r.peptide_position,
		r.domain,
		r.ligands,
		r.trans_start,
		r.trans_end,
		r.trans_strand
   
		from firestar f, firestar_residues r ";

	my @bindvalues;
	($statement, @bindvalues) = _add_condition_combine3($statement, @args);
		
	my $final = _do_query($dbh, $statement, @bindvalues);
	return $final;
}

sub query_spade {
	my ($self, %args) = @_;
	my $dbh = $self->dbh;

	my @args;
	while (my ($k, $v) = each %args){
		push @args, ({$k => $v}, "and"); # format for the_add_condition subroutine but too bad won't be scalable for "or"
	}
	if (keys(%args)){ pop @args;}  # remove final "and"

        my $statement = "select
        spade_id,
        entity_id,
		result,
		num_domains,
		num_possibly_damaged_domains,
		num_damaged_domains,
		num_wrong_domains,
		domain_signal
  
		from spade ";

	my @bindvalues;
	($statement, @bindvalues) = _add_condition($statement, @args);

	my $final = _do_query($dbh, $statement, @bindvalues);
	return $final;
}
 
sub query_spade_alignments {
	my ($self, %iargs) = @_;
	my $dbh = $self->dbh;

	my @args;
	foreach my $k (sort keys %iargs){ # important sort the keys!!!
		my $v = $iargs{$k};
		if ( $k eq 'region' ) {
			my @args1;
			# query where 2 ranges intersect
			#SELECT d.* FROM thedates d WHERE d.jobstart < $myenddate AND d.jobend > $mystartdate
			#push(@args1, {'trans_end' => $v->{'end'}});
			push(@args1, {'trans_end' => $v->{'start'}});
			push(@args1, "and");			
			#push(@args1, {'trans_start' => $v->{'start'}});
			push(@args1, {'trans_start' => $v->{'end'}});
			push(@args, [@args1]);
			push(@args, "or");
			my @args2;
			#push(@args2, {'trans_end' => $v->{'end'}});
			push(@args2, {'trans_end' => $v->{'start'}});
			push(@args2, "and");			
			#push(@args2, {'trans_start' => $v->{'start'}});
			push(@args2, {'trans_start' => $v->{'end'}});
			push(@args, [@args2]);
			push(@args, "or");
		}
		else {
			push @args, ({$k => $v}, "and"); # format for the_add_condition subroutine but too bad won't be scalable for "or"
		}
	}
	if (keys(%iargs)){ pop @args;}  # remove final "and"

        my $statement = "select
        s.spade_id,
		s.result,
		s.num_domains,
		s.num_possibly_damaged_domains,
		s.num_damaged_domains,
		s.num_wrong_domains,
		s.domain_signal,
		a.alignment_start,
		a.alignment_end,
		a.envelope_start,
		a.envelope_end,
		a.hmm_start,
		a.hmm_end,
		a.hmm_length,
		a.hmm_acc,
		a.hmm_name,
		a.hmm_type,
		a.bit_score,
		a.evalue,
		a.significance,
		a.clan,
		a.predicted_active_site_residues,
		a.trans_start,
		a.trans_end,
		a.trans_strand,
		a.type_domain,
		a.external_id,
		a.discarded
		
   
		from spade s, spade_alignments a ";

	my @bindvalues;
	($statement, @bindvalues) = _add_condition_combine10($statement, @args);
	
	my $final = _do_query($dbh, $statement, @bindvalues);
	return $final;
}

sub query_corsair {
	my ($self, %args) = @_;
	my $dbh = $self->dbh;

	my @args;
	while (my ($k, $v) = each %args){
		push @args, ({$k => $v}, "and"); # format for the_add_condition subroutine but too bad won't be scalable for "or"
	}
	if (keys(%args)){ pop @args;}  # remove final "and"

        my $statement = "select
        corsair_id,
        entity_id,
		result,        
		score,
		vertebrate_signal
  
		from corsair ";

	my @bindvalues;
	($statement, @bindvalues) = _add_condition($statement, @args);

	my $final = _do_query($dbh, $statement, @bindvalues);
	return $final;
}

sub query_matador3d {
	my ($self, %args) = @_;
	my $dbh = $self->dbh;

	my @args;
	while (my ($k, $v) = each %args){
		push @args, ({$k => $v}, "and"); # format for the_add_condition subroutine but too bad won't be scalable for "or"
	}
	if (keys(%args)){ pop @args;}  # remove final "and"

        my $statement = "select
        entity_id,
		result,
		score,
		conservation_structure
  
		from matador3d ";
  
	my @bindvalues;
	($statement, @bindvalues) = _add_condition($statement, @args);

	my $final = _do_query($dbh, $statement, @bindvalues);
	return $final;
}

sub query_matador3d_alignments {
	my ($self, %iargs) = @_;
	my $dbh = $self->dbh;

	my @args;
	foreach my $k (sort keys %iargs){ # important sort the keys!!!
		my $v = $iargs{$k};
		if ( $k eq 'region' ) {
			my @args1;
			push(@args1, {'trans_end' => $v->{'start'}});
			push(@args1, "and");			
			push(@args1, {'trans_start' => $v->{'end'}});
			push(@args, [@args1]);
			push(@args, "or");
			my @args2;
			push(@args2, {'trans_end' => $v->{'start'}});
			push(@args2, "and");			
			push(@args2, {'trans_start' => $v->{'end'}});
			push(@args, [@args2]);
			push(@args, "or");
		}
		else {
			push @args, ({$k => $v}, "and"); # format for the_add_condition subroutine but too bad won't be scalable for "or"
		}
	}
	if (keys(%iargs)){ pop @args;}  # remove final "and"

        my $statement = "select
        m.matador3d_id,
        m.result,
        m.score,
        m.conservation_structure,
        a.matador3d_alignments_id,
		a.cds_id,
		a.start,
		a.end,
		a.score as alignment_score,
		a.type,
		a.alignment_start,
		a.alignment_end,
		a.pdb_id,
		a.identity,
		a.external_id,
		a.trans_start,
		a.trans_end,
		a.trans_strand
   
		from matador3d m, matador3d_alignments a ";

	my @bindvalues;
	($statement, @bindvalues) = _add_condition_combine5($statement, @args);

	my $final = _do_query($dbh, $statement, @bindvalues);
	return $final;	
}

sub query_thump {
	my ($self, %args) = @_;
	my $dbh = $self->dbh;

	my @args;
	while (my ($k, $v) = each %args){
		push @args, ({$k => $v}, "and"); # format for the_add_condition subroutine but too bad won't be scalable for "or"
	}
	if (keys(%args)){ pop @args;}  # remove final "and"

        my $statement = "select
        entity_id,
		result,
		num_tmh,
		num_damaged_tmh,
		transmembrane_signal
  
		from thump ";
  
	my @bindvalues;
	($statement, @bindvalues) = _add_condition($statement, @args);

	my $final = _do_query($dbh, $statement, @bindvalues);
	return $final;
}

sub query_thump_helixes {
	my ($self, %iargs) = @_;
	my $dbh = $self->dbh;

	my @args;
	foreach my $k (sort keys %iargs){ # important sort the keys!!!
		my $v = $iargs{$k};
		if ( $k eq 'region' ) {
			my @args1;
			push(@args1, {'trans_end' => $v->{'start'}});
			push(@args1, "and");			
			push(@args1, {'trans_start' => $v->{'end'}});
			push(@args, [@args1]);
			push(@args, "or");
			my @args2;
			push(@args2, {'trans_end' => $v->{'start'}});
			push(@args2, "and");			
			push(@args2, {'trans_start' => $v->{'end'}});
			push(@args, [@args2]);
			push(@args, "or");
		}
		else {
			push @args, ({$k => $v}, "and"); # format for the_add_condition subroutine but too bad won't be scalable for "or"
		}
	}
	if (keys(%iargs)){ pop @args;}  # remove final "and"
	
        my $statement = "select
        t.thump_id,
        t.result,
        t.num_tmh,
        t.num_damaged_tmh,
        t.transmembrane_signal,
        h.thump_helixes_id,
		h.start,
		h.end,
		h.trans_start,
		h.trans_end,
		h.trans_strand,
		h.damaged
   
		from thump t, thump_helixes h ";

	my @bindvalues;
	($statement, @bindvalues) = _add_condition_combine11($statement, @args);

	my $final = _do_query($dbh, $statement, @bindvalues);
	return $final;
}

sub query_slr_type {
	my ($self, %args) = @_;
	my $dbh = $self->dbh;

	my @args;
	while (my ($k, $v) = each %args){
		push @args, ({$k => $v}, "and"); # format for the_add_condition subroutine but too bad won't be scalable for "or"
	}
	if (keys(%args)){ pop @args;}  # remove final "and"

        my $statement = "select
        slr_type_id,
		name,
		description
  
		from slr_type ";
  
	my @bindvalues;
	($statement, @bindvalues) = _add_condition($statement, @args);

	my $final = _do_query($dbh, $statement, @bindvalues);
	return $final;
}

sub query_inertia {
	my ($self, %args) = @_;
	my $dbh = $self->dbh;

	my @args;
	while (my ($k, $v) = each %args){
		push @args, ({$k => $v}, "and"); # format for the_add_condition subroutine but too bad won't be scalable for "or"
	}
	if (keys(%args)){ pop @args;}  # remove final "and"

        my $statement = "select
        inertia_id,
        entity_id,
        result,
		unusual_evolution
  
		from inertia ";
  
	my @bindvalues;
	($statement, @bindvalues) = _add_condition($statement, @args);

	my $final = _do_query($dbh, $statement, @bindvalues);
	return $final;
}

sub query_inertia_residues {
	my ($self, %args) = @_;
	my $dbh = $self->dbh;

	my @args;
	while (my ($k, $v) = each %args){
		push @args, ({$k => $v}, "and"); # format for the_add_condition subroutine but too bad won't be scalable for "or"
	}
	if (keys(%args)){ pop @args;}  # remove final "and"

        my $statement = "select
        i.inertia_id,
        i.entity_id,
        i.result,
        i.unusual_evolution as inertia_unusual_evolution,
		r.inertia_residues_id,
		r.inertia_exon_id,
		r.trans_start,
		r.trans_end,
		r.trans_strand,
		r.unusual_evolution as inertia_residue_unusual_evolution
		  
		from inertia i, inertia_residues r ";
  
	my @bindvalues;
	($statement, @bindvalues) = _add_condition_combine15($statement, @args);

	my $final = _do_query($dbh, $statement, @bindvalues);
	return $final;
}

sub query_omega {
	my ($self, %args) = @_;
	my $dbh = $self->dbh;

	my @args;
	while (my ($k, $v) = each %args){
		push @args, ({$k => $v}, "and"); # format for the_add_condition subroutine but too bad won't be scalable for "or"
	}
	if (keys(%args)){ pop @args;}  # remove final "and"

        my $statement = "select
        o.omega_id,
        o.slr_type_id,
		o.omega_average,
		o.omega_st_desviation,
		o.result,
		o.unusual_evolution
  
		from inertia i, omega o ";
  
	my @bindvalues;
	($statement, @bindvalues) = _add_condition_combine8($statement, @args);

	my $final = _do_query($dbh, $statement, @bindvalues);
	return $final;
}

sub query_omega_residues {
	my ($self, %args) = @_;
	my $dbh = $self->dbh;

	my @args;
	while (my ($k, $v) = each %args){
		push @args, ({$k => $v}, "and"); # format for the_add_condition subroutine but too bad won't be scalable for "or"
	}
	if (keys(%args)){ pop @args;}  # remove final "and"

        my $statement = "select
        i.inertia_id,
        i.entity_id,
        i.unusual_evolution as inertia_unusual_evolution,
        o.omega_id,
        o.slr_type_id,
		o.result,
		o.unusual_evolution,
		r.omega_residues_id,
		r.omega_exon_id,
		r.trans_start,
		r.trans_end,
		r.omega_mean,
		r.st_deviation,
		r.p_value,
		r.difference_value,
		r.unusual_evolution
		  
		from inertia i, omega o, omega_residues r ";
  
	my @bindvalues;
	($statement, @bindvalues) = _add_condition_combine4($statement, @args);

	my $final = _do_query($dbh, $statement, @bindvalues);
	return $final;
}

sub query_crash {
	my ($self, %args) = @_;
	my $dbh = $self->dbh;

	my @args;
	while (my ($k, $v) = each %args){
		push @args, ({$k => $v}, "and"); # format for the_add_condition subroutine but too bad won't be scalable for "or"
	}
	if (keys(%args)){ pop @args;}  # remove final "and"

        my $statement = "select
        crash_id,
        entity_id,
		result,
		sp_score,
		tp_score,
		peptide_signal,
		mitochondrial_signal
  
		from crash ";

	my @bindvalues;
	($statement, @bindvalues) = _add_condition($statement, @args);

	my $final = _do_query($dbh, $statement, @bindvalues);
	return $final;
}

sub query_crash_residues {
	my ($self, %iargs) = @_;
	my $dbh = $self->dbh;

	my @args;
	foreach my $k (sort keys %iargs){ # important sort the keys!!!
		my $v = $iargs{$k};
		if ( $k eq 'region' ) {
			my @args1;
			push(@args1, {'trans_end' => $v->{'start'}});
			push(@args1, "and");			
			push(@args1, {'trans_start' => $v->{'end'}});
			push(@args, [@args1]);
			push(@args, "or");
			my @args2;
			push(@args2, {'trans_end' => $v->{'start'}});
			push(@args2, "and");			
			push(@args2, {'trans_start' => $v->{'end'}});
			push(@args, [@args2]);
			push(@args, "or");
		}
		else {
			push @args, ({$k => $v}, "and"); # format for the_add_condition subroutine but too bad won't be scalable for "or"
		}
	}
	if (keys(%iargs)){ pop @args;}  # remove final "and"
	
        my $statement = "select
        c.crash_id,
		c.sp_score,
		c.tp_score,
		c.peptide_signal,
		c.mitochondrial_signal,
        r.crash_residues_id,
		r.s_mean,
		r.s_prob,
		r.d_score,
		r.c_max,
		r.reliability,
		r.localization,
		r.start,
		r.end,
		r.trans_start,
		r.trans_end,
		r.trans_strand
   
		from crash c, crash_residues r ";

	my @bindvalues;
	($statement, @bindvalues) = _add_condition_combine17($statement, @args);

	my $final = _do_query($dbh, $statement, @bindvalues);
	return $final;
}

sub query_proteo {
	my ($self, %args) = @_;
	my $dbh = $self->dbh;

	my @args;
	while (my ($k, $v) = each %args){
		push @args, ({$k => $v}, "and"); # format for the_add_condition subroutine but too bad won't be scalable for "or"
	}
	if (keys(%args)){ pop @args;}  # remove final "and"

        my $statement = "select
        proteo_id,
        entity_id,
		num_peptides,
		result,
		peptide_evidence
  
		from proteo ";

	my @bindvalues;
	($statement, @bindvalues) = _add_condition($statement, @args);

	my $final = _do_query($dbh, $statement, @bindvalues);
	return $final;
}

sub query_proteo_peptides {
	my ($self, %iargs) = @_;
	my $dbh = $self->dbh;

	my @args;
	foreach my $k (sort keys %iargs){ # important sort the keys!!!
		my $v = $iargs{$k};
		if ( $k eq 'region' ) {
			my @args1;
			push(@args1, {'trans_end' => $v->{'start'}});
			push(@args1, "and");			
			push(@args1, {'trans_start' => $v->{'end'}});
			push(@args, [@args1]);
			push(@args, "or");
			my @args2;
			push(@args2, {'trans_end' => $v->{'start'}});
			push(@args2, "and");			
			push(@args2, {'trans_start' => $v->{'end'}});
			push(@args, [@args2]);
			push(@args, "or");
		}
		else {
			push @args, ({$k => $v}, "and"); # format for the_add_condition subroutine but too bad won't be scalable for "or"
		}
	}
	if (keys(%iargs)){ pop @args;}  # remove final "and"
	
        my $statement = "select
        c.proteo_id,
		c.num_peptides,
		c.peptide_evidence,
        r.proteo_peptides_id,
		r.peptide_id,
		r.sequence,
		r.num_experiments,
		r.experiments,
		r.start,
		r.end,
		r.trans_start,
		r.trans_end,
		r.trans_strand
   
		from proteo c, proteo_peptides r ";

	my @bindvalues;
	($statement, @bindvalues) = _add_condition_combine19($statement, @args);

	my $final = _do_query($dbh, $statement, @bindvalues);
	return $final;
}

sub query_appris {
	my ($self, %args) = @_;
	my $dbh = $self->dbh;

	my @args;
	while (my ($k, $v) = each %args){
		push @args, ({$k => $v}, "and"); # format for the_add_condition subroutine but too bad won't be scalable for "or"
	}
	if (keys(%args)){ pop @args;}  # remove final "and"
	
        my $statement = "select
        entity_id,
		functional_residues_score,
		homologous_structure_score,
		vertebrate_conservation_score,
		domain_score,
		transmembrane_helices_score,
		peptide_score,
		mitochondrial_score,
		unusual_evolution_score,
		peptide_evidence_score,
		principal_isoform_score,
		functional_residues_signal,
		homologous_structure_signal,
		vertebrate_conservation_signal,
		domain_signal,
		transmembrane_helices_signal,
		peptide_signal,
		mitochondrial_signal,
		unusual_evolution_signal,
		peptide_evidence_signal,
		principal_isoform_signal,
		reliability,
		result
		
		from appris ";
  
	my @bindvalues;
	($statement, @bindvalues) = _add_condition($statement, @args);

	my $final = _do_query($dbh, $statement, @bindvalues);
	return $final;
}

##################
# INSERT methods #
##################
sub insert_entity {
	my ($self, %args) = @_;
	my $dbh = $self->dbh;

	eval {
	$dbh->do(q{insert into entity (datasource_id, identifier, source, biotype, tsl, tag) values (?,?,?,?,?,?)},
			undef,(
				$args{'datasource_id'},
				$args{'identifier'},
				$args{'source'},
				$args{'biotype'},
				$args{'tsl'},
				$args{'tag'}	));
	};
	die("\nERROR: $!\n") if $@;
	return $dbh->{'mysql_insertid'};
}  

sub insert_xref_identify {
	my ($self, %args) = @_;
	my $dbh = $self->dbh;

	eval {
	$dbh->do(q{insert into xref_identify (entity_id, datasource_id, identifier) values (?,?,?)},
			undef,(
				$args{'entity_id'},
				$args{'datasource_id'},
				$args{'identifier'} ));
	};
	die("\nERROR: $!\n") if $@;
}

sub insert_coordinate {
	my ($self, %args) = @_;
	my $dbh = $self->dbh;

	eval {
	$dbh->do(q{insert into coordinate (entity_id, chromosome, start, end, strand) values (?,?,?,?,?)},
			undef,(
				$args{'entity_id'},
				$args{'chromosome'},
				$args{'start'},
				$args{'end'},
				$args{'strand'} ));
	};
	die("\nERROR: $!\n") if $@;
}

sub insert_exon {
	my ($self, %args) = @_;
	my $dbh = $self->dbh;

	eval {
	$dbh->do(q{insert into exon (entity_id, exon_id, identifier, start, end, strand) values (?,?,?,?,?,?)},
			undef,(
				$args{'entity_id'},
				$args{'exon_id'},
				$args{'identifier'},
				$args{'start'},
				$args{'end'},
				$args{'strand'} ));
	};
	die("\nERROR: $!\n") if $@;
}

sub insert_codon {
	my ($self, %args) = @_;
	my $dbh = $self->dbh;

	eval {
	$dbh->do(q{insert into codon (entity_id, type, start, end, strand, phase) values (?,?,?,?,?,?)},
			undef,(
				$args{'entity_id'},
				$args{'type'},
				$args{'start'},
				$args{'end'},
				$args{'strand'},
				$args{'phase'} ));
	};
	die("\nERROR: $!\n") if $@;
}

sub insert_cds {
	my ($self, %args) = @_;
	my $dbh = $self->dbh;

	eval {
	$dbh->do(q{insert into cds (entity_id, cds_id, start, end, strand, phase) values (?,?,?,?,?,?)},
			undef,(
				$args{'entity_id'},
				$args{'cds_id'},
				$args{'start'},
				$args{'end'},
				$args{'strand'},
				$args{'phase'} ));
	};
	die("\nERROR: $!\n") if $@;
}

sub insert_sequence {
	my ($self, %args) = @_;
	my $dbh = $self->dbh;

	eval {
	$dbh->do(q{insert into sequence (entity_id, type_id, length, sequence) values (?,?,?,?)},
			undef,(
				$args{'entity_id'},
				$args{'type_id'},
				$args{'length'},
				$args{'sequence'} ));
	};
	die("\nERROR: $!\n") if $@;
}
  
sub insert_alignment {
	my ($self, %args) = @_;
	my $dbh = $self->dbh;

	eval {
	$dbh->do(q{insert into alignment (entity_id, type_id, length, num_species, alignment) values (?,?,?,?,?)},
			undef,(
				$args{'entity_id'},
				$args{'type_id'},
				$args{'length'},
				$args{'num_species'},
				$args{'alignment'} ));
	};
	die("\nERROR: $!\n") if $@;
}

sub insert_extended_sequence {
	my ($self, %args) = @_;
	my $dbh = $self->dbh;

	eval {
	$dbh->do(q{insert into extended_sequence (entity_id, length, sequence, description) values (?,?,?,?)},
			undef,(
				$args{'entity_id'},
				$args{'length'},
				$args{'sequence'},
				$args{'description'} ));
	};
	die("\nERROR: $!\n") if $@;
}

sub insert_signalp {
	my ($self, %args) = @_;
	my $dbh = $self->dbh;

	eval {
	$dbh->do(q{insert into signalp (entity_id, result, s_mean, s_prob, d_score, c_max, score, start, end, trans_start, trans_end, trans_strand, peptide_signal) values (?,?,?,?,?,?,?,?,?,?,?,?,?)},
			undef,(
				$args{'entity_id'},
				$args{'result'},
				$args{'s_mean'},
				$args{'s_prob'},
				$args{'d_score'},
				$args{'c_max'},
				$args{'score'},
				$args{'start'},
				$args{'end'},
				$args{'trans_start'},
				$args{'trans_end'},
				$args{'trans_strand'},
				$args{'peptide_signal'} ));
	};
	die("\nERROR: $!\n") if $@;
}

sub insert_targetp {
	my ($self, %args) = @_;
	my $dbh = $self->dbh;

	eval {
	$dbh->do(q{insert into targetp (entity_id, result, reliability, localization, score, start, end, trans_start, trans_end, trans_strand, mitochondrial_signal) values (?,?,?,?,?,?,?,?,?,?,?)},
			undef,(
				$args{'entity_id'},
				$args{'result'},
				$args{'reliability'},
				$args{'localization'},
				$args{'score'},
				$args{'start'},
				$args{'end'},
				$args{'trans_start'},
				$args{'trans_end'},
				$args{'trans_strand'},
				$args{'mitochondrial_signal'} ));
	};
	die("\nERROR: $!\n") if $@;
}

sub insert_firestar {
	my ($self, %args) = @_;
	my $dbh = $self->dbh;

	eval {
	$dbh->do(q{insert into firestar (entity_id, num_residues, result, functional_residue) values (?,?,?,?)},
			undef,(
				$args{'entity_id'},
				$args{'num_residues'},
				$args{'result'},
				$args{'functional_residue'} ));
	};
	die("\nERROR: $!\n") if $@;
	return $dbh->{'mysql_insertid'};
}

sub insert_firestar_residues {
	my ($self, %args) = @_;
	my $dbh = $self->dbh;

	eval {
	$dbh->do(q{insert into firestar_residues (
										firestar_id,
										peptide_position,
										domain,
										ligands,
										trans_start,
										trans_end,
										trans_strand
									)
									values (?,?,?,?,?,?,?)},
			undef,(
				$args{'firestar_id'},
				$args{'peptide_position'},
				$args{'domain'},
				$args{'ligands'},
				$args{'trans_start'},
				$args{'trans_end'},
				$args{'trans_strand'} ));
	};
	die("\nERROR: $!\n") if $@;
	return $dbh->{'mysql_insertid'};
}

sub insert_matador3d {
	my ($self, %args) = @_;
	my $dbh = $self->dbh;

	eval {
	$dbh->do(q{insert into matador3d (entity_id, result, score, conservation_structure) values (?,?,?,?)},
			undef,(
				$args{'entity_id'},
				$args{'result'},
				$args{'score'},
				$args{'conservation_structure'} ));
	};
	die("\nERROR: $!\n") if $@;
	return $dbh->{'mysql_insertid'};
}

sub insert_matador3d_alignments {
	my ($self, %args) = @_;
	my $dbh = $self->dbh;

	eval {
	$dbh->do(q{insert into matador3d_alignments (matador3d_id, cds_id, start, end, score, type, alignment_start, alignment_end, pdb_id, identity, external_id, trans_start, trans_end, trans_strand) values (?,?,?,?,?,?,?,?,?,?,?,?,?,?)},
			undef,(
				$args{'matador3d_id'},
				$args{'cds_id'},
				$args{'start'},
				$args{'end'},
				$args{'score'},
				$args{'type'},
				$args{'alignment_start'},
				$args{'alignment_end'},
				$args{'pdb_id'},
				$args{'identity'},
				$args{'external_id'},
				$args{'trans_start'},
				$args{'trans_end'},
				$args{'trans_strand'} ));
	};
	die("\nERROR: $!\n") if $@;
	return $dbh->{'mysql_insertid'};
}
  
sub insert_thump {
	my ($self, %args) = @_;
	my $dbh = $self->dbh;

	eval {
	$dbh->do(q{insert into thump (entity_id, result, num_tmh, num_damaged_tmh, transmembrane_signal) values (?,?,?,?,?)},
			undef,(
				$args{'entity_id'},
				$args{'result'},
				$args{'num_tmh'},
				$args{'num_damaged_tmh'},
				$args{'transmembrane_signal'} ));
	};
	die("\nERROR: $!\n") if $@;
	return $dbh->{'mysql_insertid'};	
}

sub insert_thump_helixes {
	my ($self, %args) = @_;
	my $dbh = $self->dbh;

	eval {
	$dbh->do(q{insert into thump_helixes (thump_id, start, end, trans_start, trans_end, trans_strand, damaged) values (?,?,?,?,?,?,?)},
			undef,(
				$args{'thump_id'},
				$args{'start'},
				$args{'end'},
				$args{'trans_start'},
				$args{'trans_end'},
				$args{'trans_strand'},
				$args{'damaged'} ));
	};
	die("\nERROR: $!\n") if $@;
	return $dbh->{'mysql_insertid'};
}

sub insert_spade {
	my ($self, %args) = @_;
	my $dbh = $self->dbh;

	eval {
	$dbh->do(q{insert into spade (entity_id, result, num_domains, num_possibly_damaged_domains, num_damaged_domains, num_wrong_domains, domain_signal) values (?,?,?,?,?,?,?)},
			undef,(
				$args{'entity_id'},
				$args{'result'},
				$args{'num_domains'},
				$args{'num_possibly_damaged_domains'},
				$args{'num_damaged_domains'},
				$args{'num_wrong_domains'},
				$args{'domain_signal'} ));
	};
	die("\nERROR: $!\n") if $@;
	return $dbh->{'mysql_insertid'};	
}

sub insert_spade_alignments {
	my ($self, %args) = @_;
	my $dbh = $self->dbh;

	eval {
	$dbh->do(q{insert into spade_alignments (spade_id, 
												alignment_start, alignment_end, 
												envelope_start, envelope_end, 
												hmm_start, hmm_end, hmm_length, hmm_acc, hmm_name, hmm_type, 
												bit_score, evalue, significance, clan, predicted_active_site_residues,
												trans_start, trans_end, trans_strand, type_domain, external_id, discarded) values (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)},
			undef,(
				$args{'spade_id'},
				$args{'alignment_start'},
				$args{'alignment_end'},
				$args{'envelope_start'},
				$args{'envelope_end'},
				$args{'hmm_start'},
				$args{'hmm_end'},
				$args{'hmm_length'},
				$args{'hmm_acc'},
				$args{'hmm_name'},
				$args{'hmm_type'},
				$args{'bit_score'},
				$args{'evalue'},
				$args{'significance'},
				$args{'clan'},
				$args{'predicted_active_site_residues'},
				$args{'trans_start'},
				$args{'trans_end'},
				$args{'trans_strand'},
				$args{'type_domain'},
				$args{'external_id'},
				$args{'discarded'} ));
	};
	die("\nERROR: $!\n") if $@;
	return $dbh->{'mysql_insertid'};
}

sub insert_corsair {
	my ($self, %args) = @_;
	my $dbh = $self->dbh;

	eval {
	$dbh->do(q{insert into corsair (entity_id, result, score, vertebrate_signal) values (?,?,?,?)},
			undef,(
				$args{'entity_id'},
				$args{'result'},
				$args{'score'},
				$args{'vertebrate_signal'} ));
	};
	die("\nERROR: $!\n") if $@;
	return $dbh->{'mysql_insertid'};	
}

sub insert_corsair_alignments {
	my ($self, %args) = @_;
	my $dbh = $self->dbh;

	eval {
	$dbh->do(q{insert into corsair_alignments (corsair_id, cds_id, start, end, score, maxscore, sp_report, type, trans_start, trans_end, trans_strand) values (?,?,?,?,?,?,?,?,?,?,?)},
			undef,(
				$args{'corsair_id'},
				$args{'cds_id'},
				$args{'start'},
				$args{'end'},
				$args{'score'},
				$args{'maxscore'},
				$args{'sp_report'},
				$args{'type'},
				$args{'trans_start'},
				$args{'trans_end'},
				$args{'trans_strand'} ));
	};
	die("\nERROR: $!\n") if $@;
	return $dbh->{'mysql_insertid'};
}

sub insert_inertia {
	my ($self, %args) = @_;
	my $dbh = $self->dbh;
 
	eval {
	$dbh->do(q{insert into inertia (entity_id, result, unusual_evolution) values (?,?,?)},
			undef,(
				$args{'entity_id'},
				$args{'result'},
				$args{'unusual_evolution'}));
	};
	die("\nERROR: $!\n") if $@;
	return $dbh->{'mysql_insertid'};
}

sub insert_inertia_residues {
	my ($self, %args) = @_;
	my $dbh = $self->dbh;

	eval {
	$dbh->do(q{insert into inertia_residues (inertia_id, inertia_exon_id, trans_start, trans_end, trans_strand, unusual_evolution) values (?,?,?,?,?,?)},
			undef,(
				$args{'inertia_id'},
				$args{'inertia_exon_id'},
				$args{'trans_start'},
				$args{'trans_end'},
				$args{'trans_strand'},
				$args{'unusual_evolution'} ));				
	};
	die("\nERROR: $!\n") if $@;
	return $dbh->{'mysql_insertid'};
}

sub insert_omega {
	my ($self, %args) = @_;
	my $dbh = $self->dbh;
 
	eval {
	$dbh->do(q{insert into omega (inertia_id, slr_type_id, omega_average, omega_st_desviation, result, unusual_evolution) values (?,?,?,?,?,?)},
			undef,(
				$args{'inertia_id'},
				$args{'slr_type_id'},
				$args{'omega_average'},
				$args{'omega_st_desviation'},				
				$args{'result'},
				$args{'unusual_evolution'}));
	};
	die("\nERROR: $!\n") if $@;
	return $dbh->{'mysql_insertid'};
}

sub insert_omega_residues {
	my ($self, %args) = @_;
	my $dbh = $self->dbh;

	eval {
	$dbh->do(q{insert into omega_residues (omega_id, omega_exon_id, trans_start, trans_end, omega_mean, st_deviation, p_value, difference_value, unusual_evolution) values (?,?,?,?,?,?,?,?,?)},
			undef,(
				$args{'omega_id'},
				$args{'omega_exon_id'},
				$args{'trans_start'},
				$args{'trans_end'},
				$args{'omega_mean'},
				$args{'st_deviation'},
				$args{'p_value'},
				$args{'difference_value'},
				$args{'unusual_evolution'} ));				
	};
	die("\nERROR: $!\n") if $@;
	return $dbh->{'mysql_insertid'};
}

sub insert_crash {
	my ($self, %args) = @_;
	my $dbh = $self->dbh;

	eval {
	$dbh->do(q{insert into crash (entity_id, result, sp_score, tp_score, peptide_signal, mitochondrial_signal) values (?,?,?,?,?,?)},
			undef,(
				$args{'entity_id'},
				$args{'result'},
				$args{'sp_score'},
				$args{'tp_score'},
				$args{'peptide_signal'},
				$args{'mitochondrial_signal'} ));
	};
	die("\nERROR: $!\n") if $@;
	return $dbh->{'mysql_insertid'};
}

sub insert_crash_residues {
	my ($self, %args) = @_;
	my $dbh = $self->dbh;

	eval {
	$dbh->do(q{insert into crash_residues (crash_id, s_mean, s_prob, d_score, c_max, reliability, localization, start, end, trans_start, trans_end, trans_strand) values (?,?,?,?,?,?,?,?,?,?,?,?)},
			undef,(
				$args{'crash_id'},
				$args{'s_mean'},
				$args{'s_prob'},
				$args{'d_score'},
				$args{'c_max'},
				$args{'reliability'},
				$args{'localization'},
				$args{'start'},
				$args{'end'},
				$args{'trans_start'},
				$args{'trans_end'},
				$args{'trans_strand'} ));
	};
	die("\nERROR: $!\n") if $@;
	return $dbh->{'mysql_insertid'};
}

sub insert_proteo {
	my ($self, %args) = @_;
	my $dbh = $self->dbh;

	eval {
	$dbh->do(q{insert into proteo (entity_id, num_peptides, result, peptide_evidence) values (?,?,?,?)},
			undef,(
				$args{'entity_id'},
				$args{'num_peptides'},
				$args{'result'},
				$args{'peptide_evidence'} ));
	};
	die("\nERROR: $!\n") if $@;
	return $dbh->{'mysql_insertid'};
}

sub insert_proteo_peptides {
	my ($self, %args) = @_;
	my $dbh = $self->dbh;

	eval {
	$dbh->do(q{insert into proteo_peptides (
										proteo_id,
										peptide_id,
										sequence,
										num_experiments,
										experiments,
										start,
										end,
										trans_start,
										trans_end,
										trans_strand
									)
									values (?,?,?,?,?,?,?,?,?,?)},
			undef,(
				$args{'proteo_id'},
				$args{'peptide_id'},
				$args{'sequence'},
				$args{'num_experiments'},
				$args{'experiments'},
				$args{'start'},
				$args{'end'},
				$args{'trans_start'},
				$args{'trans_end'},
				$args{'trans_strand'} ));
	};
	die("\nERROR: $!\n") if $@;
	return $dbh->{'mysql_insertid'};
}

sub insert_appris {
	my ($self, %args) = @_;
	my $dbh = $self->dbh;
	
	eval {
	$dbh->do(q{insert into appris (
						        entity_id,
								functional_residues_score,
								homologous_structure_score,
								vertebrate_conservation_score,
								domain_score,
								transmembrane_helices_score,
								peptide_score,
								mitochondrial_score,
								unusual_evolution_score,
								peptide_evidence_score,
								principal_isoform_score,
								functional_residues_signal,
								homologous_structure_signal,
								vertebrate_conservation_signal,
								domain_signal,
								transmembrane_helices_signal,
								peptide_signal,
								mitochondrial_signal,
								unusual_evolution_signal,
								peptide_evidence_signal,
								principal_isoform_signal,
								reliability,
								result
								
		) values (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)},
			undef,(
				$args{'entity_id'},
				$args{'functional_residues_score'},
				$args{'homologous_structure_score'},
				$args{'vertebrate_conservation_score'},
				$args{'domain_score'},
				$args{'transmembrane_helices_score'},
				$args{'peptide_score'},
				$args{'mitochondrial_score'},
				$args{'unusual_evolution_score'},
				$args{'peptide_evidence_score'},
				$args{'principal_isoform_score'},				
				$args{'functional_residues_signal'},
				$args{'homologous_structure_signal'},
				$args{'vertebrate_conservation_signal'},
				$args{'domain_signal'},
				$args{'transmembrane_helices_signal'},
				$args{'peptide_signal'},
				$args{'mitochondrial_signal'},
				$args{'unusual_evolution_signal'},
				$args{'peptide_evidence_signal'},
				$args{'principal_isoform_signal'},			
				$args{'reliability'},
				$args{'result'}));
	};
	die("\nERROR: $!\n") if $@;
}

##################
# UPDATE methods #
##################
sub update_entity {
	my ($self, %args) = @_;
	my $dbh = $self->dbh;

	my $id;
	my @args;
	while (my ($k, $v) = each %args) {
		if ($k eq 'entity_id') {
			$id = $args{'entity_id'};
		} else {
			push @args, ({$k => $v}, ","); # format for the_add_condition subroutine but too bad won't be scalable for "or"
		}
	}
	if (keys(%args)){ pop @args;}  # remove final ", "

	my $statement = "update entity set ";
	my @bindvalues;
	($statement, @bindvalues) =_add_condition_update($statement, @args);
	$statement .= "where entity_id = $id ";
 
	my $final = do_sql_sentence($dbh, $statement, @bindvalues);
	return $final;
}

sub update_xref_identify {
	my ($self, %args) = @_;
	my $dbh = $self->dbh;

	my $id;
	my $datasource_id;
	my @args;
	while (my ($k, $v) = each %args) {
		if ($k eq 'entity_id') {
			$id = $args{'entity_id'};
		} elsif ($k eq 'datasource_id') {
			$datasource_id = $args{'datasource_id'};
		} else {
			push @args, ({$k => $v}, ","); # format for the_add_condition subroutine but too bad won't be scalable for "or"
		}
	}
	if (keys(%args)){ pop @args;}  # remove final ", "

	my $statement = "update xref_identify set ";
	my @bindvalues;
	($statement, @bindvalues) =_add_condition_update($statement, @args);
	$statement .= "where entity_id = $id and datasource_id = $datasource_id ";
 
	my $final = do_sql_sentence($dbh, $statement, @bindvalues);
	return $final;
}

sub update_alignment {
	my ($self, %args) = @_;
	my $dbh = $self->dbh;

	my $id;
	my @args;
	while (my ($k, $v) = each %args) {
		if ($k eq 'entity_id') {
			$id = $args{'entity_id'};
		} else {
			push @args, ({$k => $v}, ","); # format for the_add_condition subroutine but too bad won't be scalable for "or"
		}
	}
	if (keys(%args)){ pop @args;}  # remove final ", "

	my $statement = "update alignment set ";
	my @bindvalues;
	($statement, @bindvalues) =_add_condition_update($statement, @args);
	$statement .= "where entity_id = $id ";

	my $final = do_sql_sentence($dbh, $statement, @bindvalues);
	return $final;
}

sub update_exon {
	my ($self, %args) = @_;
	my $dbh = $self->dbh;

	my $id;
	my @args;
	while (my ($k, $v) = each %args) {
		if ($k eq 'entity_id') {
			$id = $args{'entity_id'};
		} else {
			push @args, ({$k => $v}, ","); # format for the_add_condition subroutine but too bad won't be scalable for "or"
		}
	}
	if (keys(%args)){ pop @args;}  # remove final ", "

	my $statement = "update exon set ";
	my @bindvalues;
	($statement, @bindvalues) =_add_condition_update($statement, @args);
	$statement .= "where entity_id = $id ";

	my $final = do_sql_sentence($dbh, $statement, @bindvalues);
	return $final;
}

sub update_exon2 {
	my ($self, %args) = @_;
	my $dbh = $self->dbh;

	my $id;
	my $start;
	my $end;
	my $strand;
	my @args;
	while (my ($k, $v) = each %args) {
		if ($k eq 'entity_id') {
			$id = $args{'entity_id'};
		} elsif ($k eq 'start') {
			$start = $args{'start'};
		} elsif ($k eq 'end') {
			$end = $args{'end'};
		} elsif ($k eq 'strand') {
			$strand = $args{'strand'};
		} else {
			push @args, ({$k => $v}, ","); # format for the_add_condition subroutine but too bad won't be scalable for "or"
		}
	}
	if (keys(%args)){ pop @args;}  # remove final ", "

	my $statement = "update exon set ";
	my @bindvalues;
	($statement, @bindvalues) =_add_condition_update($statement, @args);
	$statement .= "where entity_id = $id and start = $start and end = $end and strand = '$strand' ";

	my $final = do_sql_sentence($dbh, $statement, @bindvalues);
	
	return $final;
}

sub update_signalp {
	my ($self, %args) = @_;
	my $dbh = $self->dbh;

	my $id;
	my @args;
	while (my ($k, $v) = each %args) {
		if ($k eq 'entity_id') {
			$id = $args{'entity_id'};
		} else {
			push @args, ({$k => $v}, ","); # format for the_add_condition subroutine but too bad won't be scalable for "or"
		}
	}
	if (keys(%args)){ pop @args;}  # remove final ", "

	my $statement = "update signalp set ";
	my @bindvalues;
	($statement, @bindvalues) =_add_condition_update($statement, @args);
	$statement .= "where entity_id = $id ";

	my $final = do_sql_sentence($dbh, $statement, @bindvalues);
	return $final;
}

sub update_targetp {
	my ($self, %args) = @_;
	my $dbh = $self->dbh;

	my $id;
	my @args;
	while (my ($k, $v) = each %args) {
		if ($k eq 'entity_id') {
			$id = $args{'entity_id'};
		} else {
			push @args, ({$k => $v}, ","); # format for the_add_condition subroutine but too bad won't be scalable for "or"
		}
	}
	if (keys(%args)){ pop @args;}  # remove final ", "

	my $statement = "update targetp set ";
	my @bindvalues;
	($statement, @bindvalues) =_add_condition_update($statement, @args);
	$statement .= "where entity_id = $id ";

	my $final = do_sql_sentence($dbh, $statement, @bindvalues);
	return $final;
}

sub update_firestar {
	my ($self, %args) = @_;
	my $dbh = $self->dbh;

	my $id;
	my @args;
	while (my ($k, $v) = each %args) {
		if ($k eq 'entity_id') {
			$id = $args{'entity_id'};
		} else {
			push @args, ({$k => $v}, ","); # format for the_add_condition subroutine but too bad won't be scalable for "or"
		}
	}
	if (keys(%args)){ pop @args;}  # remove final ", "

	my $statement = "update firestar set ";
	my @bindvalues;
	($statement, @bindvalues) =_add_condition_update($statement, @args);
	$statement .= "where entity_id = $id ";

	my $final = do_sql_sentence($dbh, $statement, @bindvalues);
	return $final;
}

sub update_matador3d {
	my ($self, %args) = @_;
	my $dbh = $self->dbh;

	my $id;
	my @args;
	while (my ($k, $v) = each %args) {
		if ($k eq 'entity_id') {
			$id = $args{'entity_id'};
		} else {
			push @args, ({$k => $v}, ","); # format for the_add_condition subroutine but too bad won't be scalable for "or"
		}
	}
	if (keys(%args)){ pop @args;}  # remove final ", "

	my $statement = "update matador3d set ";
	my @bindvalues;
	($statement, @bindvalues) =_add_condition_update($statement, @args);
	$statement .= "where entity_id = $id ";

	my $final = do_sql_sentence($dbh, $statement, @bindvalues);
	return $final;
}

sub update_inertia {
	my ($self, %args) = @_;
	my $dbh = $self->dbh;

	my $id;
	my $inertia_id;
	my @args;
	while (my ($k, $v) = each %args) {
		if ($k eq 'entity_id') {
			$id = $args{'entity_id'};
		} elsif ($k eq 'inertia_id') {
			$inertia_id = $args{'inertia_id'};
		} else {
			push @args, ({$k => $v}, ","); # format for the_add_condition subroutine but too bad won't be scalable for "or"
		}
	}
	if (keys(%args)){ pop @args;}  # remove final ", "

	my $statement = "update inertia set ";
	my @bindvalues;
	($statement, @bindvalues) =_add_condition_update($statement, @args);
	$statement .= "where entity_id = $id ";

	my $final = do_sql_sentence($dbh, $statement, @bindvalues);
	return $final;
}

sub update_omega {
	my ($self, %args) = @_;
	my $dbh = $self->dbh;

	my $id;
	my $omega_id = '';
	my @args;
	while (my ($k, $v) = each %args) {
		if ($k eq 'entity_id') {
			$id = $args{'entity_id'};
		} elsif ($k eq 'omega_id') {
			$omega_id = " and omega_id = '".$args{'omega_id'}."' ";
		} else {
			push @args, ({$k => $v}, ","); # format for the_add_condition subroutine but too bad won't be scalable for "or"
		}
	}
	if (keys(%args)){ pop @args;}  # remove final ", "

	my $statement = "update omega set ";
	my @bindvalues;
	($statement, @bindvalues) =_add_condition_update($statement, @args);
	$statement .= "where entity_id = $id $omega_id ";

	my $final = do_sql_sentence($dbh, $statement, @bindvalues);
	return $final;
}

sub update_omega_residues {
	my ($self, %args) = @_;
	my $dbh = $self->dbh;

	my $id;
	my @args;
	while (my ($k, $v) = each %args) {
		if ($k eq 'omega_id') {
			$id = $args{'omega_id'};
		} else {
			push @args, ({$k => $v}, ","); # format for the_add_condition subroutine but too bad won't be scalable for "or"
		}
	}
	if (keys(%args)){ pop @args;}  # remove final ", "

	my $statement = "update omega_residues set ";
	my @bindvalues;
	($statement, @bindvalues) =_add_condition_update($statement, @args);
	$statement .= "where omega_id = $id ";

	my $final = do_sql_sentence($dbh, $statement, @bindvalues);
	return $final;
}

sub update_proteo {
	my ($self, %args) = @_;
	my $dbh = $self->dbh;

	my $id;
	my @args;
	while (my ($k, $v) = each %args) {
		if ($k eq 'entity_id') {
			$id = $args{'entity_id'};
		} else {
			push @args, ({$k => $v}, ","); # format for the_add_condition subroutine but too bad won't be scalable for "or"
		}
	}
	if (keys(%args)){ pop @args;}  # remove final ", "

	my $statement = "update proteo set ";
	my @bindvalues;
	($statement, @bindvalues) =_add_condition_update($statement, @args);
	$statement .= "where entity_id = $id ";

	my $final = do_sql_sentence($dbh, $statement, @bindvalues);
	return $final;
}

sub update_appris {
	my ($self, %args) = @_;
	my $dbh = $self->dbh;

	my $id;
	my @args;
	while (my ($k, $v) = each %args) {
		if ($k eq 'entity_id') {
			$id = $args{'entity_id'};
		} else {
			push @args, ({$k => $v}, ","); # format for the_add_condition subroutine but too bad won't be scalable for "or"
		}
	}
	if (keys(%args)){ pop @args;}  # remove final ", "

	my $statement = "update appris set ";
	my @bindvalues;
	($statement, @bindvalues) =_add_condition_update($statement, @args);
	$statement .= "where entity_id = $id ";

	my $final = do_sql_sentence($dbh, $statement, @bindvalues);
	return $final;
}



sub commit {
	my ($self) = @_;
	my $dbh = $self->dbh;
	$dbh->commit();
}

sub disconnect {
	my ($self) = @_;
	my $dbh = $self->dbh;
	$dbh->disconnect();

}
sub DESTROY {}

1;
