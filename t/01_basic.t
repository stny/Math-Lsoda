use strict;
use Test::More (tests => 3);

use Math::Lsoda;

my $solver = Math::Lsoda->new(verbose => 1);

isa_ok ($solver, 'Math::Lsoda');
is ($solver->verbose, 1, 'verbose_check');
is ($solver->verbose(0), 0, 'verbose_check');
print STDERR scalar(localtime), "\n";
$solver->register;
$solver->run;

sub eqns {
  my ($t, $x, $y) = @_;
  @$y[0] = 1.0e+4 * @$x[1] * @$x[2] - 0.04 * @$x[0];
  @$y[2] = 3.0e+7 * @$x[1] * @$x[1];
  @$y[1] = -(@$y[0] + @$y[2]);
}
my @x = (1.0, 0.0, 0.0);
my @ans = (59999.96, -120059999.96, 120000000);
#my $f = *open("data.dat", "w");
open( FILEHANDLE, "data.dat");
my $fh = *FILEHANDLE;
Math::Lsoda::solve(\&eqns, \@x, 0.0, 1.0, 0.1, $fh);
close(FILEHANDLE);
#foreach (@x) { print STDERR " $_\n" }
#is(Math::Lsoda::solve(\&eqns, \@x, 0.0, 0.4, 0.1), -3, "error_test");
#foreach (@x) { print STDERR " $_\n" }
