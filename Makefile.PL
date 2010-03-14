use inc::Module::Install;
name 'Math-Lsoda';
all_from 'lib/Math/Lsoda.pm';

my $compiler;
my $cflags;

require ExtUtils::F77;
if ($^O eq "darwin") {
    ExtUtils::F77->import(qw{Darwin GFortran});
    $compiler     = ExtUtils::F77->compiler();
    $cflags       = "-ffixed-form -std=legacy -O2";
} else {
    ExtUtils::F77->import();
    $cflags       = "-ffixed-form -std=legacy -O2 -fPIC";
}

can_cc or die "This module requires a C compiler";
makemaker_args(
    LIBS => ExtUtils::F77->runtime,
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
    $str .= "\t$compiler $cflags -c \$<";
    $str;
}
