=head1 CONTACT

  Please email comments or questions to the public INB
  developers list at <inb-tecnico@lists.cnio.es>.

  Questions may also be sent to the developer, 
  Jose Manuel Rodriguez <jmrodriguez@cnio.es>.

=cut

=head1 NAME

appris - common process

=head1 SYNOPSIS

=head1 DESCRIPTION


=head1 METHODS

=cut

package appris;

use strict;
use warnings;
use Bio::SeqIO;
use File::Temp;
use Config::IniFiles;
use MIME::Lite;

use APPRIS::Parser qw( parse_gencode parse_infiles );
use APPRIS::Utils::File qw( getTotalStringFromFile );
use APPRIS::Utils::Argument qw( rearrange );
use APPRIS::Utils::Exception qw( info throw warning deprecate );

use Exporter;

use vars qw(@ISA @EXPORT);

@ISA = qw(Exporter);
@EXPORT = qw(
	create_gene_list
	get_gene_list
	retrieve_gene_list
	create_ensembl_input
	reduce_input
	create_gencode_data
	create_indata
	send_email
);


=head2 create_gene_list

  Arg[1]      : (optional) String $text - notification text to present to user
  Example     : # run a code snipped conditionally
                if ($support->user_proceed("Run the next code snipped?")) {
                    # run some code
                }

                # exit if requested by user
                exit unless ($support->user_proceed("Want to continue?"));
  Description : If running interactively, the user is asked if he wants to
                perform a script action. If he doesn't, this section is skipped
                and the script proceeds with the code. When running
                non-interactively, the section is run by default.
  Return type : TRUE to proceed, FALSE to skip.
  Exceptions  : none
  Caller      : general

=cut

sub create_gene_list
{
	my ( $itype, $id,  $list,  $pos,  $gdata, $gtransc, $gtransl, $econf, $ever, $spe) = rearrange( 
		['type', 'id', 'list', 'pos', 'gdata','gtransc','gtransl','econf','ever', 'spe' ],
	@_ );
	my ($report, $data) = (undef, undef);
	
	if ( $itype =~ /datafile/ and defined $gdata )
	{
		my ($data_fh);
		if ( defined $list ) {
			$report = get_gene_list($list);
			$data_fh = reduce_input($gdata, undef, $report);
			if ( UNIVERSAL::isa($data_fh,'File::Temp') ) {
				$gdata = $data_fh->filename;
			}
		}		
		elsif ( defined $pos ) {
			$data_fh = reduce_input($gdata, $pos);
			if ( UNIVERSAL::isa($data_fh,'File::Temp') ) {
				$gdata = $data_fh->filename;
			}
		}
		$data = create_indata($gdata, $gtransc, $gtransl);
						
		# delete tmp file
		if ( defined $data_fh and UNIVERSAL::isa($data_fh,'File::Temp') ) {
			$data_fh->unlink_on_destroy(1);
		}
	}	
	elsif ( $itype =~ /ensembl/ and defined $econf and defined $ever and defined $spe )
	{
		if ( defined $id ) {
			$report->{$id} = 1;
		}
		elsif ( defined $list ) {
			$report = get_gene_list($list);
		}
		elsif ( defined $pos ) {
			$report = create_ensembl_input($pos, $econf, $ever, $spe);			
		}
		else {
			$report = create_ensembl_input(undef, $econf, $ever, $spe);			
		}
	}
	
	File::Temp::cleanup();
	
	return ($report, $data);
}

=head2 get_gene_list

  Arg[1]      : (optional) String $text - notification text to present to user
  Example     : # run a code snipped conditionally
                if ($support->user_proceed("Run the next code snipped?")) {
                    # run some code
                }

                # exit if requested by user
                exit unless ($support->user_proceed("Want to continue?"));
  Description : If running interactively, the user is asked if he wants to
                perform a script action. If he doesn't, this section is skipped
                and the script proceeds with the code. When running
                non-interactively, the section is run by default.
  Return type : TRUE to proceed, FALSE to skip.
  Exceptions  : none
  Caller      : general

=cut

sub get_gene_list($)
{
	my ($file) = @_;
	my ($list);
	
	my ($genes) = getTotalStringFromFile($file);
	foreach my $gene_id (@{$genes}) {
		$gene_id =~ s/\s*//mg;
		$list->{$gene_id} = 1 if ( $gene_id ne '');
	}
	return $list;
} # end get_gene_list

=head2 retrieve_gene_list

  Arg[1]      : (optional) String $text - notification text to present to user
  Example     : # run a code snipped conditionally
                if ($support->user_proceed("Run the next code snipped?")) {
                    # run some code
                }

                # exit if requested by user
                exit unless ($support->user_proceed("Want to continue?"));
  Description : If running interactively, the user is asked if he wants to
                perform a script action. If he doesn't, this section is skipped
                and the script proceeds with the code. When running
                non-interactively, the section is run by default.
  Return type : TRUE to proceed, FALSE to skip.
  Exceptions  : none
  Caller      : general

=cut

sub retrieve_gene_list($;$)
{
	my ($data_file, $position) = @_;
	my ($list);
		
	# customize data for a position
	if ( defined $position ) {
		foreach my $pos ( split(',', $position) ) {
			if ( defined $pos ) {
				eval {
					my ($cmd) =	"awk '{if( (\$1==\"$pos\" || \$1==\"chr$pos\") && \$3==\"gene\") {print \$10}}' $data_file | sed 's/[\"|;]//g' | sed 's/\.[0-9]*\$//'";
					info("** script: $cmd\n");
					my (@aux_list) = `$cmd`;
					foreach my $gene_id (@aux_list) {
						if ( defined $gene_id and ($gene_id ne '') ) {
							$gene_id =~ s/\s*//mg;
							$list->{$gene_id} = 1 if ( $gene_id ne '');
						}
					}
				};
				throw("creating input for $pos") if($@);
			}		
		}		
	}
	return $list;
} # end retrieve_gene_list

=head2 create_ensembl_input

  Arg[1]      : (optional) String $text - notification text to present to user
  Example     : # run a code snipped conditionally
                if ($support->user_proceed("Run the next code snipped?")) {
                    # run some code
                }

                # exit if requested by user
                exit unless ($support->user_proceed("Want to continue?"));
  Description : If running interactively, the user is asked if he wants to
                perform a script action. If he doesn't, this section is skipped
                and the script proceeds with the code. When running
                non-interactively, the section is run by default.
  Return type : TRUE to proceed, FALSE to skip.
  Exceptions  : none
  Caller      : general

=cut

sub create_ensembl_input($$$$)
{
	my ($position, $conf_file, $e_version, $species) = @_;
	my ($gene_list);
	
	# create ensembl registry from default config file
	my ($cfg) = new Config::IniFiles( -file => $conf_file );
	my ($e_core_param) = {
			'-host'       => $cfg->val( 'ENSEMBL_CORE_REGISTRY', 'host'),
			'-user'       => $cfg->val( 'ENSEMBL_CORE_REGISTRY', 'user'),
			'-pass'       => $cfg->val( 'ENSEMBL_CORE_REGISTRY', 'pass'),
			'-verbose'    => $cfg->val( 'ENSEMBL_CORE_REGISTRY', 'verbose'),
			'-db_version' => $e_version,
			'-species'    => $species,
	};
	info(
		"\t-host        => ".$cfg->val( 'ENSEMBL_CORE_REGISTRY', 'host')."\n".
		"\t-user        => ".$cfg->val( 'ENSEMBL_CORE_REGISTRY', 'user')."\n".
		"\t-pass        => ".$cfg->val( 'ENSEMBL_CORE_REGISTRY', 'pass')."\n".
		"\t-verbose     => ".$cfg->val( 'ENSEMBL_CORE_REGISTRY', 'verbose')."\n".
		"\t-db_version  => ".$e_version."\n".										
		"\t-species     => ".$species."\n"
	);
	require Bio::EnsEMBL::Registry;
	my ($registry) = 'Bio::EnsEMBL::Registry';
	eval {
		$registry->load_registry_from_db(%{$e_core_param});
	};
	throw("can not load ensembl registry: $!\n") if $@;
	
	# create slice adaptor
	my ($slice_adaptor) = $registry->get_adaptor( $species, 'Core', 'Slice' );	
	
	# create ensembl input for each given position
	my (@slices);
	if ( defined $position ) {
		@slices = split(',', $position);
	}
	else {
 		foreach my $slice ( @{$slice_adaptor->fetch_all('chromosome')} ) { 			
 			push(@slices, $slice->seq_region_name);
 		}		
	}
	
	# get gene list depending on given position
	foreach my $pos ( @slices ) {
		if ( defined $pos ) {
			$pos =~ s/^chr//;		
			my ($slice) = $slice_adaptor->fetch_by_region( 'chromosome', $pos );
			if ( defined $slice ) {
				my ($genes) = $slice->get_all_Genes();
				if ( defined $genes ) {
					while ( my $gene = shift @{$genes} ) {
						my ($gene_id) = $gene->stable_id;
						my ($gene_eid) = $gene_id;
						#if ( $gene->version ) {
						#	$gene_eid = $gene_id.'.'.$gene->version;
						#}
						$gene_list->{$gene_eid} = $gene->biotype();
					}					
				}
			}
		}		
	}
	
	return $gene_list;
	
} # end create_ensembl_input

=head2 reduce_input

  Arg[1]      : (optional) String $text - notification text to present to user
  Example     : # run a code snipped conditionally
                if ($support->user_proceed("Run the next code snipped?")) {
                    # run some code
                }

                # exit if requested by user
                exit unless ($support->user_proceed("Want to continue?"));
  Description : If running interactively, the user is asked if he wants to
                perform a script action. If he doesn't, this section is skipped
                and the script proceeds with the code. When running
                non-interactively, the section is run by default.
  Return type : TRUE to proceed, FALSE to skip.
  Exceptions  : none
  Caller      : general

=cut

sub reduce_input($;$;$)
{
	my ($data_file, $position, $gene_list) = @_;
	my ($data_tmpfile) = File::Temp->new( UNLINK => 0, SUFFIX => '.appris.dat' );
	my ($data_tmpfilename) = $data_tmpfile->filename;
		
	# customize data for a position
	if ( defined $position ) {
		foreach my $pos ( split(',', $position) ) {
			if ( defined $pos ) {
				eval {
					my ($cmd) =	"awk '{if(\$1==\"$pos\" || \$1==\"chr$pos\") {print \$0}}' $data_file >> $data_tmpfilename";
					info("** script: $cmd\n");
					my (@cmd_out) = `$cmd`;
				};
				throw("creating input for $pos") if($@);
			}		
		}		
	}
	# get customized data for a gene list
	elsif ( defined $gene_list ) {
		my ($g_cond) = '';	
		foreach my $g ( keys(%{$gene_list}) ) {
			$g_cond .= ' $10 ~ /'.$g.'/ ||'; # for rel7 version		
		}
		if ( $g_cond ne '' ) {
	    	$g_cond =~ s/\|\|$//;
	    	$g_cond = '('.$g_cond.')';			
			eval {
				my ($cmd) = "awk '{if(\$9==\"gene_id\" && $g_cond) {print \$0}}' $data_file >> $data_tmpfilename";
				info("** script: $cmd\n");
				my (@cmd_out) = `$cmd`;
			};
			throw("creating gencode data for genelist") if($@);			
		}		
	}	
	
	return $data_tmpfile;
	
} # end reduce_input

=head2 create_gencode_data

  Arg[1]      : (optional) String $text - notification text to present to user
  Example     : # run a code snipped conditionally
                if ($support->user_proceed("Run the next code snipped?")) {
                    # run some code
                }

                # exit if requested by user
                exit unless ($support->user_proceed("Want to continue?"));
  Description : If running interactively, the user is asked if he wants to
                perform a script action. If he doesn't, this section is skipped
                and the script proceeds with the code. When running
                non-interactively, the section is run by default.
  Return type : TRUE to proceed, FALSE to skip.
  Exceptions  : none
  Caller      : general

=cut

sub create_gencode_data($;$;$)
{
	my ($data_file, $transcripts_file, $translations_file) = @_;
	
	my ($data) = parse_gencode($data_file, $transcripts_file, $translations_file);
	unless ( defined $data ) {
		throw("can not create gencode object\n");
	}
	
	return $data;
	
} # end create_gencode_data

=head2 create_ensembl_data

  Arg[1]      : (optional) String $text - notification text to present to user
  Example     : # run a code snipped conditionally
                if ($support->user_proceed("Run the next code snipped?")) {
                    # run some code
                }

                # exit if requested by user
                exit unless ($support->user_proceed("Want to continue?"));
  Description : If running interactively, the user is asked if he wants to
                perform a script action. If he doesn't, this section is skipped
                and the script proceeds with the code. When running
                non-interactively, the section is run by default.
  Return type : TRUE to proceed, FALSE to skip.
  Exceptions  : none
  Caller      : general

=cut

sub create_indata($;$;$)
{
	my ($data_file, $transcripts_file, $translations_file) = @_;
	
	my ($data) = parse_infiles($data_file, $transcripts_file, $translations_file);
	unless ( defined $data ) {
		throw("can not create gencode object\n");
	}
	
	return $data;
	
} # end create_ensembl_data

=head2 send_email

  Arg[1]      : (optional) String $text - notification text to present to user
  Example     : # run a code snipped conditionally
                if ($support->user_proceed("Run the next code snipped?")) {
                    # run some code
                }

                # exit if requested by user
                exit unless ($support->user_proceed("Want to continue?"));
  Description : If running interactively, the user is asked if he wants to
                perform a script action. If he doesn't, this section is skipped
                and the script proceeds with the code. When running
                non-interactively, the section is run by default.
  Return type : TRUE to proceed, FALSE to skip.
  Exceptions  : none
  Caller      : general

=cut

sub send_email($$$$$$)
{
	my ($econf, $email, $wserver, $species, $methods, $outpath) = @_;
	
	# Attach results if webserver
	my ($gzip_file);
	unless ( defined $wserver ) { # local execution
		my (@inpath_n) = split('/', $outpath);		
		if ( scalar(@inpath_n) > 0 ) {
			my ($num) = scalar(@inpath_n);
			$wserver = $inpath_n[$num-1];
		}
	}
	else
	{ # webserver
		# create tmp dir for email sending
		my ($email_tmp_dir) = $outpath."/email";
		eval {
			my ($cmd) =	"mkdir $email_tmp_dir";
			my (@cmd_out) = `$cmd`;
		};
		throw("creating tmp dir $email_tmp_dir") if($@);
		
		# create tmp files from methods
		foreach my $method ( split(',', $methods) ) {
			eval {
				my ($cmd) =	"find $outpath -type f -name '*.$method' ! -size 0 -exec cp -v {} $email_tmp_dir/$wserver.$method \\;";
				my (@cmd_out) = `$cmd`;
			};
			throw("creating tmp files $method in $email_tmp_dir") if($@);
		}
		
		# create gzip file
		my ($gzip_file) = "$email_tmp_dir/$wserver.tar.gz";
		eval {
			my ($cmd) =	"cd $email_tmp_dir && tar -cf - * | gzip -9 > $gzip_file";
			my (@cmd_out) = `$cmd`;
		};
		throw("compressing file $gzip_file in $email_tmp_dir") if($@);
	}
	
	my ($cfg) = new Config::IniFiles( -file => $econf );
	my ($cfg_from) = $cfg->val('APPRIS_EMAIL', 'from');
	my ($cfg_subject) = $cfg->val('APPRIS_EMAIL', 'subject');
	
	# send email
	my ($subject,$content);
	if ( defined $ENV{APPRIS_WSERVER_REPORT_URL} ) {
		$subject = $cfg_subject." Your jobid $wserver has finished";
		$content = "<p>The results for your query $wserver are now available at ";	
		my ($rst_url) = $ENV{APPRIS_WSERVER_REPORT_URL} . '/' . $wserver;
		$content	 .= "<a href='$rst_url'>$rst_url</a></p>";
		$content	 .= "</p>";
		#$content	 .= "<p>Please remember that your results will be stored in our server only for 3 weeks, than will be deleted !!</p>";		
	}
	else {
		$subject = $cfg_subject." Your jobid $species/$wserver has finished";
		$content = "<p>The results for your query $species/$wserver are now available.";	
	}
	_semail($cfg_from, $email, $subject, $content, $gzip_file);
	
	# rm email dir
	#if ( -d $email_tmp_dir ) {
	#	$logger->info("-- delete email dir\n");
	#	eval {
	#		my ($cmd) =	"rm -rf $email_tmp_dir";
	#		$logger->debug("** script: $cmd\n");
	#		my (@cmd_out) = `$cmd`;
	#	};
	#	$logger->error("deleting file\n") if($@);
	#}
					
} # end send_email

sub _semail($$$$;$)
{
	my ($from, $to, $subject, $content, $files) = @_;
	my ($msg) = MIME::Lite->new(
        From     	=> $from,
        To 			=> $to,
        Subject 	=> $subject,
		Type		=> 'multipart/mixed',
	);
	$msg->attach(
		Type    	=> 'text/html',
		Encoding 	=> 'quoted-printable',
		Data    	=> $content,
	);
	if ( defined $files ) {
		foreach my $file (split(';', $files)) {
			$msg->attach(
				Type		 => 'x-gzip',
				Path    	 => $file,
				Disposition	=> 'attachment'
			);		
		}		
	}
 	eval { $msg->send }; die "MIME::Lite->send failed: $@\n" if $@;
 	
} # end _semail


1;
