#!perl
#**********************************************************
#*                                                        *
#*     XTD2XYZ - Convert XTD files into XYZ ormat        *
#*                                                        *
#**********************************************************
# Version: 0.1
# Author: Andrea Minoia
# Date: 08/09/2010
#
# Convert MS trajectory xtd file into xYZ trajectory file.
# Backup of files that are about to be overwritten is managed
# by MS. The most recent file is that with higher index number (N)
# The script has to be in the same directory of the
# structure to modify and the user has to update the
# variable $doc (line 31) according to the name of the
# file containing the trajectory.
# The xmol trajectory is stored in trj.txt file and it is not
# possible to rename the file within MS, nor it is possible to
# automatically export it as xyz or car file. You should manage
# the new trajectory manually for further use (e.g. VMD)
#
# Modificator: Sobereva (sobereva@sina.com)
# Date: 2012-May-23
# The range of the frames to be outputted can be altered by line 49 and 51
# Modificator: nanxu (tamas@zju.edu.cn)
# Date: 2022-Jau-03
# Add support for lattice

use strict;
use MaterialsScript qw(:all);

#open the multiframe trajectory structure file or die
my $doc = $Documents{"./Si.xtd"};

if (!$doc) {die "no document";}

my $trajectory = $doc->Trajectory;

if ($trajectory->NumFrames>1) {

    print "Found ".$trajectory->NumFrames." frames in the trajectory\n";
    # Open new xmol trajectory file
    my $xmolFile=Documents->New("trj.txt");
   
    #get atoms in the structure
#    my $atoms = $doc->Atoms;
    my $atoms = $doc->DisplayRange->Atoms;
    my $Natoms=@$atoms;

    # loops over the frames
    my $framebegin=1;
    my $frameend=$trajectory->NumFrames;
    # my $frameend=10;
    for (my $frame=$framebegin; $frame<=$frameend; ++$frame){
        $trajectory->CurrentFrame = $frame;
        #write header xyz
        $xmolFile->Append(sprintf "%i \n", $Natoms);
    my $lattice = $doc->Lattice3D;
    my $a = $lattice->LengthA;
    my $b = $lattice->LengthB;
    my $c = $lattice->Lengthc;    
    my $alpha = $lattice->AngleAlpha /180.0 * 3.1415926;
    my $beta = $lattice->AngleBeta /180.0 * 3.1415926;
    my $gamma = $lattice->AngleGamma /180.0 * 3.1415926;

    my $bc2 = $b**2 + $c**2 - 2*$b*$c*cos($alpha);

    my $h1 = $a;
    my $h2 = $b * cos($gamma);
    my $h3 = $b * sin($gamma);
    my $h4 = $c * cos($beta);
    my $h5 = (($h2 - $h4)**2 + $h3**2 + $c**2 - $h4**2 - $bc2)/(2 * $h3);
    my $h6 = sqrt($c**2 - $h4**2 - $h5**2);
    $xmolFile->Append(sprintf "Lattice=\"%f %f %f %f %f %f %f %f %f\" Properties=species:S:1:pos:R:3 frame=%i pbc=\"T T T\"\n", $h1,0,0,$h2,$h3,0,$h4,$h5,$h6,$frame);
        foreach my $atom (@$atoms) {
            # write atom symbol and x-y-z- coordinates
            $xmolFile->Append(sprintf "%s %f  %f  %f \n",$atom->ElementSymbol, $atom->X, $atom->Y,

$atom->Z);
        }   
    }
    #close trajectory file
    $xmolFile->Close;
}
else {
    print "The " . $doc->Name . " is not a multiframe trajectory file \n";
}

