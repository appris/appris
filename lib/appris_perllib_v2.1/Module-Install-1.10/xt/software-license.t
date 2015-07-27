use strict;
use warnings FATAL => 'all';

use English qw( -no_match_vars );
local $OUTPUT_AUTOFLUSH = 1;

#BEGIN {
#	unless ($ENV{RELEASE_TESTING}) {
#		use Test::More;
#		Test::More::plan(
#			skip_all => 'Author tests, not required for installation.');
#	}
#}
use Test::More;
use Test::Requires {
	'Test::Software::License' => 0.002000,
};

# all_software_license_ok();

# the following is brutal, if it exists it must have a license
all_software_license_ok({ strict => 1 });

done_testing();

__END__

