#!/usr/bin/env perl
use strict;

sub help()
{
        my $usage = qq{
        CheckModulesPerl <directory>
        };
        print STDERR $usage."\n";
}

sub get_version
{
        my (%args)=@_;
        my $RES;
        my $mod=$args{-module};
        my $CMD="perl -M$mod -e 'print \"\$".$mod."::VERSION\"' 2>&1";
        open $RES ,"$CMD |";
       
        while(my $line=<$RES>)
        {
                return $line;
        }
        close($RES);
       
        return -1;
}
sub exists_mod_perl
{
        my (%args)=@_;
        my $RES;
        my $mod=$args{-module};
        my $CMD="perl -e 'use ".$mod."' 2>&1";
        open $RES ,"$CMD |";
        while(my $line=<$RES>)
        {
                return -1;
        }
        close($RES);
       
        return 0;

}

sub search_mod_perl
{
        my (%args)=@_;
        my $RES;
        my %modules;
       
        my $CMD="ls -R ".$args{-base};
        open $RES ,"$CMD |";
        my $path='';
        while(my $line=<$RES>)
        {
                chomp $line;
                if($line ne '')
                {
                        if($line =~ /(.*):$/)
                        {
                                $path=$1;
                        }
                        else
                        {
                                my $file=$path.'/'.$line;
                                if(-f $file)
                                {
                                        $file=quotemeta($file);
                                        my $CMD="grep 'use ' $file";
                                        my $RES;
                                        open $RES ,"$CMD |";
                                        while(my $line=<$RES>)
                                        {
                                                if($line =~ /use\s+(\w+(::\w+)*)(\s+qw\(.+\))?\s*;/)
                                                {
                                                        $modules{$1}=1;
                                                }
                                        }
                                        close($RES);
                                       
                                }
                        }
                }
        }
        close($RES);
        return %modules;
}
my $dir = $ARGV[0];

if(!-e $dir)
{
        print STDERR "ERROR: Dir Not Exists\n";
        help();
        exit;
}


my %modules=search_mod_perl(-base=>$dir);

foreach my $current_module (sort keys %modules)
{
        my $status_msg="Check Module Perl: ".$current_module;
        my $status=exists_mod_perl(-module=>$current_module);
       
        print $status_msg." " x (50-length($status_msg));
       
        if($status==-1)
        {
                print "\t\t\tNot Found\n";
        }
        else
        {
                # print $status."------\n";
                print "\t\t\t".get_version(-module=>$current_module)."v\n";
        }
       
}

__END__
