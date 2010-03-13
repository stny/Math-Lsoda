use inc::Module::Install;
name 'Math-Lsoda';
all_from 'lib/Math/Lsoda.pm';

my $mycompiler;
my $mycflags;

require ExtUtils::F77;
if ($^O eq "darwin") {
    ExtUtils::F77->import(qw{Darwin GFortran});
    $mycompiler     = ExtUtils::F77->compiler();
    $mycflags       = "-ffixed-form -std=legacy -O2"; # -fPIC
} else {
    ExtUtils::F77->import();
    $mycflags       = "-ffixed-form -std=legacy -O2 -fPIC";
}

can_cc or die "This module requires a C compiler";
makemaker_args(
    LIBS => ExtUtils::F77->runtime,
#    OBJECT => "Lsoda.o odepack/opkda1.o odepack/opkda2.o odepack/opkdmain.o",
    OBJECT => "Lsoda.o opkda1.o opkda2.o opkdmain.o",
    clean => {
        FILES => q{
            *.stackdump
            *.gcov *.gcda *.gcno
            *.out
            nytprof
            cover_db
        },
    },
);

tests 't/*.t';
author_tests 'xt';

build_requires 'Test::More';
use_test_base;
auto_include;
WriteAll;

sub MY::postamble {
    my $str  = "";
    $str .= ".f.o :\n";
    $str .= "\t$mycompiler $mycflags -c \$<";
    $str;
}
