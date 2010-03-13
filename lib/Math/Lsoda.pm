package Math::Lsoda;
use Any::Moose;

our $VERSION = '0.01';
our @ISA;

has verbose => (
    is    => 'rw',
    isa   => 'Int',
    default => 1,
    );

__PACKAGE__->meta->make_immutable;

no Any::Moose;

eval {
    require XSLoader;
    XSLoader::load(__PACKAGE__, $VERSION);
    1;
} or do {
    require DynaLoader;
    push @ISA, 'DynaLoader';
    __PACKAGE__->bootstrap($VERSION);
};

sub register {
    my ($self, @args) = @_;
}
sub run {
    my $self = shift;
}

1;
__END__

=head1 NAME

Math::Lsoda - solve ordinary differential equation systems using lsoda.

=head1 SYNOPSIS

  use Math::Lsoda;

=head1 DESCRIPTION

Math::Lsoda is

=head1 AUTHOR

Naoya Sato E<lt>synclover@gmail.comE<gt>

=head1 SEE ALSO

=head1 LICENSE

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.

=cut
