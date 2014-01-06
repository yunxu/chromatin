#/usr/bin/perl -w
use File::Path;
use File::Basename;
use List::Util qw(min max);
use strict;

my $num_args = $#ARGV+1;
if ($num_args != 6) {
  print "\nUsage: ENm008.stat2.pl sample_size  prefix cutoff pardir celldir\n";
  exit;
}
my $SampleSize = $ARGV[0];
my $Prefix = $ARGV[1];
my $Cutoff = $ARGV[2];
my $MaxLogWeight = $ARGV[3];
my $ParDir = $ARGV[4];
my $CellDir = $ARGV[5];
#my $SampleSize = "640_r6000";
#my $Prefix = "AA";


my $DataDir = "../result/ENm008/chain/".$ParDir."/".$CellDir;

#my $OutDir = "/tmp/".$SampleSize;
# my $OutDir = "/tmp/".$Prefix;
my $OutDir = "../result/analysis/ENm008/chain/".$ParDir."/".$CellDir."/cut".$Cutoff;
#print $OutDir."\n";
if ( ! -d $OutDir){mkpath($OutDir);}

#my $Prefix = "AB";
#my @FileList = <$DataDir/$Prefix/*.pts>;
my @FileList = <$DataDir/$Prefix/*.pts>;
@FileList = @{FileList[0..int($#{FileList})*1/10]};
#my $Count = 0;
my %HoA;
foreach my $LenFile (@FileList){
 print $LenFile."\n";

	my %SegHash;
  open FILE, $LenFile or die "Can not open file ".$LenFile;
  my $LogWeight ;
  my $ExpDiffLogWeight;

  my $SegInd = 1;
  while(<FILE>){
    chomp;
    if (/^#/){
      if (/^# LogWeight= (.*)$/){
        $LogWeight = $1;
        $ExpDiffLogWeight = exp( $LogWeight - $MaxLogWeight);
      }
      next;
    }
    my ($x, $y, $z) = split;
		# print $SegInd;
        # my @ValsArr = ($x, $y, $z, max($r_x,$r_y,$r_z));
		my @ValsArr = ($x, $y, $z, $Cutoff);
		$SegHash{$SegInd} = [@ValsArr];
    $SegInd ++;
  }
  close FILE;

	my @SegArr = sort {$a <=> $b} keys %SegHash;

	my $OutFile = $OutDir."/".$Prefix."_".basename($LenFile,".pts").".con";
	print $OutFile."\n";
	open FILE, ">$OutFile" or die "Can not write file $OutFile\n";
    printf FILE "# ExpDiffLogWeight= %.3e\n", $ExpDiffLogWeight;    
    for (my $i=0; $i<scalar (@SegArr); $i++){
		for (my $j=$i+1; $j<scalar(@SegArr); $j++){
			my $dx = $SegHash{$SegArr[$i]}[0] - $SegHash{$SegArr[$j]}[0];
			my $dy = $SegHash{$SegArr[$i]}[1] - $SegHash{$SegArr[$j]}[1];
			my $dz = $SegHash{$SegArr[$i]}[2] - $SegHash{$SegArr[$j]}[2];
			my $distsquar = sqrt($dx*$dx + $dy*$dy + $dz*$dz);
			my $r1r2 = $Cutoff;
			
			if ($distsquar < $r1r2){
				printf FILE "%s\t%s\t%.3e\n", $SegArr[$i], $SegArr[$j], $ExpDiffLogWeight;
			}
		}
	}
	close FILE;
	
	# last;
}

