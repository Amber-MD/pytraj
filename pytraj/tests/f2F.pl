#!/usr/bin/perl
use strict;
use warnings;
use File::Find;
use Getopt::Long;

#script to convert file suffixes of Fortran files and eliminate the
#need for a separate FPP call

my $verbose = 0;
my $dryrun = 0;

GetOptions("-verbose!"=>\$verbose,
           "-dryrun!"=>\$dryrun);

usage() if($#ARGV>=0);
$verbose=1 if $dryrun;

my @dir = ("src", "AmberTools/src");

#recursively go through the two source directories and pass the files to f2F
find(\&f2F,@dir);

#remove the -P FPP flag in configure.  It will be interpreted by the
#compiler and bad things may happend
my $configure = "AmberTools/src/configure";
print "$configure: Modifying FPPFLAGS\n" if ($verbose);
open(CON,"<$configure") || die "ERROR: coud not open '$configure':$!\n";
my @configure = <CON>;
close CON;
for my $line (@configure){
  $line=~s/(fpp.*cpp -traditional)/$1 -P/
    if( !($line=~/ -P[ "]?/));
  $line=~s/(fppflags=)"-P "/$1""/;
}
if(!$dryrun){
  open(CON,">$configure") || die "ERROR: coud not open '$configure':$!\n";
  print CON @configure;
  close CON;
}

exit;

#converts .f suffixes to .F suffixes in both the actual file names and
#in Makefiles, makedepend and depend files. Many files are excluded.
#Mostly those belonging to external projects.
sub f2F {

#
#Exclusions
#

  #exclude external projects
  return if( $File::Find::name =~ /\/(lapack|arpack|odrpack|netcdf|mtkpp|[p]?netcdf|pubpme|mopac6|fftw-(3.3|3.2.2|2.1.5))\//);
  #We need to fix up the BLAS makefile but nothing else
  return if( $File::Find::name =~ /\/blas\// && !(/Makefile/));
  #We need to fix up the BLAS makefile but nothing else
  return if( $File::Find::name =~ /\/arpack\// && !(/Makefile/));

  #stay away from PMEMD and cpptraj
  return if( $File::Find::dir =~ /(pmemd)/);

  #Files that don't change
  return if($File::Find::name =~
     /\/(af-nmr|nab)\//);
  return if($File::Find::name =~
     /\/Rotdif.f$/);

#
#modify source file in special cases
#
  if($File::Find::name =~ /\/sander\/ew_fft.f$/){
    print "$_: Updating filename in include statement to .F90\n" if($verbose);
    open(SRC,"<$_") || die "ERROR: coud not open '$_':$!\n";
    my @source = <SRC>;
    close SRC;

    for my $line (@source){
      $line=~s/^(#include +".*\.)f(")/$1F90$2/;
    }
    if(!$dryrun){
      open(SRC,">$_") || die "ERROR: coud not open '$_':$!\n";
      print SRC @source;
      close SRC;
    }
  }

#
#suffixes
#

  #.f -> .F
  if(/\.f$/ && 
     (($File::Find::name =~ /\/lib\// && !/(random|veclib)/)
     || ($File::Find::name =~ /\/(etc|anal|addles|protonate|mm_pbsa|nmr_aux\/fantasian)\//))
    ){
    (my $newname = $_) =~ s/\.f/.F/;
    my $cmd = "git mv $_ $newname";
    print "$_: $cmd\n" if($verbose);
    if(!$dryrun){
      !system($cmd) || die("ERROR: Git failed:$!\n");
    }
    return;
  }

  #.f90 -> .F90
  if(/\.f90$/ ){
    (my $newname = $_) =~ s/\.f90/.F90/;
    my $cmd = "git mv $_ $newname";
    print "$_: $cmd\n" if($verbose);
    if(!$dryrun){
      !system($cmd) || die("ERROR: Git failed:$!\n");
    }
    return;
  }

  #.f -> .F90
  if(/\.f$/){
    (my $newname = $_) =~ s/\.f/.F90/;
    my $cmd = "git mv $_ $newname";
    print "$_: $cmd\n" if($verbose);
    if(!$dryrun){
      !system($cmd) || die("ERROR: Git failed:$!\n");
    }
  }

#
#Makefiles
#

  #check for makefiles that build Fortran files
  if(/^Makefile$/ || $File::Find::name=~/cpptraj\/src\/Makefile_at$/){
    open(MK,"<$_") || die "ERROR: coud not open '$_':$!\n";
    my @makefile = <MK>;
    close MK;

    #determine if this file handles fortran at all
    #and make a few universal changes:
    #  -get rid of '_*.f's
    #  -change '_$<' to '$<'
    #  -remove leading underscores for explicitly mentioned Fortran files
    my $dotf=0;
    for(my $i=0; $i <=$#makefile; $i++){
      $dotf++ if($makefile[$i] =~ /\.f +/ || $makefile[$i] =~/\.f$/);
      $makefile[$i] =~ s/_\*\.f//g;
      $makefile[$i] =~ s/_(\$<)/$1/g;
      $makefile[$i] =~ s/(\s)_(\S+\.f)/$1$2/g;
      #SPECIAL CASE
      if($File::Find::name =~ /\/pbsa\//){
        if($makefile[$i-1] =~ /\$\(FPP\).*SANDER/){
          $makefile[$i] =~ s/(\$\(FC\))/$1 -DSANDER/;
        }
        if($makefile[$i-1] =~ /\$\(FPP\).*LIBPBSA/){
          $makefile[$i] =~ s/(\$\(FC\))/$1 -DLIBPBSA/;
        }
      }
    }
    return if (!$dotf);
    print "$_: $dotf '.f' lines... " if($verbose);

    #special case: don't actually do anything with fortran
    if($File::Find::name =~ /\/(nab|nmr_aux\/prepare_input)\//){
      print "skipping.\n" if($verbose);
      return;
    }

    #special case: Pure fixed format (and etc/)
    if($File::Find::name =~ /\/(addles|anal|protonate|etc|mm_pbsa|nmr_aux\/fantasian)\//){
      print "fixed format.\n" if($verbose);
      for my $line (@makefile){
        $line =~ s/\.f/.F/g;
      }
    }

    #special case: mixed format
    #
    #Go back and fix this.  We can have a single default rule and no
    #special cases
    #
    if($File::Find::name =~ /\/(lib)\//){
      print "mixed format special case.\n" if($verbose);
      for my $line (@makefile){
        #first treat random.F and veclib.F
        $line =~ s/(random|veclib)\.f/$1.F90/g 
          if ($File::Find::name =~ /\/lib\//);
        #add globs to treat both suffixes
#        $line =~ s/_\*\.f/_*.F _*.F90/g;
        $line =~ s/\.f/.F/g;
      }
    }

    if($File::Find::name =~ /\/(rism|xray)\//){
      print "mixed format special case\n" if($verbose);
      for my $line (@makefile){
        $line =~ s/(rism_nxtsec)\.f/$1.F/g 
          if ($File::Find::name =~ /\/rism\//);
        $line =~ s/(nxtsec)\.f/$1.F/g 
          if ($File::Find::name =~ /\/xray\//);
        #add globs to treat both suffixes
#        $line =~ s/_\*\.f/_*.F _*.F90/g;
        $line =~ s/\.f/.F90/g;
      }
    }

    #cpptraj
    if(/Makefile_at/){
      print "mixed format special case\n" if($verbose);
      for my $line (@makefile){
        next if($line =~ /Rotdif/);
        $line =~ s/\.f/.F90/g;
      }
    }

    #default case: all free format
    if($File::Find::name =~ 
       /\/(nmode|sander|chamber|pbsa|sqm|nmr_aux\/rotdif|ptraj|sff|cpptraj)\//
      && !/Makefile_at/){
      print "free format\n" if($verbose);
      for my $line (@makefile){
        $line =~ s/\.f9/.F9/g;
        $line =~ s/\.f/.F90/g;
      }
    }

    #special case : add .SUFFIXES
    if($File::Find::name =~ /\/pbsa\//){
      splice @makefile,10,0,".SUFFIXES : .F90\n";
    }elsif($File::Find::name =~ /\/nmode\//){
      splice @makefile,22,0,".SUFFIXES : .F90\n";
    }

    #Remove $(FPP) references
    removeFPP(\@makefile);

    if(!$dryrun){
      open(MK,">$_") || die "ERROR: could not open '$_':$!\n";
      print MK @makefile;
      close MK;
    }

  }

#
#MAKEDEPEND
#

  #check for makedepend files that handle Fortran files
  if(/^makedepend$/){
    open(MD,"<$_") || die "ERROR: coud not open '$_':$!\n";
    my @makedepend = <MD>;
    close MD;

    #determine if this file handles fortran at all
    #and removes leading underscores for explicitly mentioned Fortran files
    my $dotf=0;
    for(my $i=0; $i <=$#makedepend; $i++){
      $dotf++ if($makedepend[$i] =~ /\.f[^[:alpha:]]/ || $makedepend[$i] =~/\.f$/);
      $makedepend[$i] =~ s/(\s)_(\S+\.f)/$1$2/g;
    }
    return if (!$dotf);
    print "$_: $dotf .f lines.\n" if($verbose);

    for my $line (@makedepend){

      #Special cases
      specialFPPFlags(\$line);
      $line =~ s/\.f90/.F90/g;
      $line =~ s/\.f/.F90/g;
    }

    #Remove $(FPP) references
    removeFPP(\@makedepend);

    if(!$dryrun){
      open(MD,">$_") || die "ERROR: could not open '$_':$!\n";
      print MD @makedepend;
      close MD;
    }

  }

#
#DEPEND
#

  #check for depend files that handle Fortran files
  if(/^depend$/){
    open(DP,"<$_") || die "ERROR: coud not open '$_':$!\n";
    my @depend = <DP>;
    close DP;

    #determine if this file handles fortran at all
    #and removes leading underscores for explicitly mentioned Fortran files
    my $dotf=0;
    for(my $i=0; $i <=$#depend; $i++){
      $dotf++ if($depend[$i] =~ /\.f[^[:alpha:]]/ || $depend[$i] =~/\.f$/);
      $depend[$i] =~ s/(\s)_(\S+\.f)/$1$2/g;
    }
    return if (!$dotf);
    print "$_: $dotf .f lines.\n" if($verbose);

    for my $line (@depend){
      #Special cases
      specialFPPFlags(\$line);

      $line =~ s/\.f90/.F90/g;
      $line =~ s/\.f([^h])/.F90$1/g;
    }

    #Remove $(FPP) references
    removeFPP(\@depend);


    if(!$dryrun){
      open(DP,">$_") || die "ERROR: could not open '$_':$!\n";
      print DP @depend;
      close DP;
    }

  }
}

#removes explicit FPP calls and adds FPPFLAGS to the FC line if it is
#not already there.  This is intended for makedepend and depend files
# $_[0] : references to an array of lines
sub removeFPP{
  my $line = shift;

  for(my $i=scalar(@$line)-1; $i>=0 ; $i--){
    if(!($$line[$i] =~ /FPPFLAGS/)){
      if($$line[$i] =~/\\\$\(FC\)/){
        $$line[$i] =~ s/(\$\(FC\))/$1 \\\$(FPPFLAGS)/;
      }else{
        $$line[$i] =~ s/(\$\(FC\))/$1 \$(FPPFLAGS)/;
      }
    }
    if($$line[$i] =~/\$\(FPP\)/){
      splice @$line, $i,1;
    }
  }
}


sub specialFPPFlags {
  my $line = shift;
      if($File::Find::name =~ /\/rism\/depend/){
        $$line =~ s/(\$\(FC\))/$1 -DRISM -DRISM_CRDINTERP -DRISM_DX -DRISM_LINPROJ \$(SPEC_FPP)/;
      }elsif($File::Find::name =~ /\/rism\/makedepend/){
        $$line =~ s/(\$\(FC\))/$1 -DRISM -DRISM_CRDINTERP -DRISM_DX -DRISM_LINPROJ \\\$(SPEC_FPP)/;
        $$line =~ s/(rism_nxtsec)\.f/$1.F/g 
          if ($File::Find::name =~ /\/rism\//);
      }
      if($$line =~ /\$\(FC\).*\.SQM\.f/){
        $$line =~ s/(\$\(FC\))/$1 -DSQM/;
        $$line =~ s/\.SQM\.f/.f/;
      }
      if($$line =~ /\$\(FC\).*\.LES\.f/){
        $$line =~ s/(\$\(FC\))/$1 -DLES/;
        $$line =~ s/\.LES\.f/.f/;
      }
      if($$line =~ /\$\(FC\).*\.PUPIL\.f/){
        $$line =~ s/(\$\(FC\))/$1 -DPUPIL_SUPPORT/;
        $$line =~ s/\.PUPIL\.f/.f/;
      }
      if($$line =~ /\$\(FC\).*\.APBS\.f/){
        $$line =~ s/(\$\(FC\))/$1 -DAPBS/;
        $$line =~ s/\.APBS\.f/.f/;
      }
      if($$line =~ /\$\(FC\).*\.RISM\.f/){
        $$line =~ s/(\$\(FC\))/$1 -DRISM -DRISM_CRDINTERP -DRISM_DX -DRISM_LINPROJ -DNO_SANDER_DIVCON/;
        $$line =~ s/\.RISM\.f/.f/;
      }
}

sub usage{
  printf STDERR "USAGE: f2F.pl [-verbose] [-dryrun]\n";
}
