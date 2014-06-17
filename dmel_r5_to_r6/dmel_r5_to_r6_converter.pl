#!/usr/local/bin/perl
# standalone_flybase_coord_converter.pl

#$debug= 1;

@ChrMapperStart= ();
@ChrMapperEnd= ();
@ChrMapperDiff= ();
@ChrMapperNewchr= ();
$ChrMapper;
$ChrMapperKey= 'none';
$CurrRelease;
%Mapper= (); 

	my @r5fbids=();
	my @r5notes= (); 
	my @r6fbids=();
	my @r6notes= (); 
	my @errids=();
	my %RelFbids= (); 

	my @ids= ();
	while( @ARGV>0 ) {
		my $coord= shift(@ARGV);
		if( $coord=~/^\s*([^ :]+):([0-9,]+)/ ) { # looks like coord
			push(@ids,$coord);
			}
  	elsif( -f $coord ) { # assume file with coords
			print "ids from file ".$coord." <br />\n" if $debug;
			my @strings= `/bin/cat $coord`;
  		foreach( @strings ){
  	  	push(@ids, split(/[ ;\t\n\r]+/,$_) );
  	  	}
			}
  	}
	foreach my $idd ( @ids ) { $id=~s/[ ,\t]+//g; $id=~s/^Chr([^:]+:)/$1/i; }
	print "No. initial ids == ".@ids." <br />\n" if $debug;
	
	$ChrMapper= setConversionArrays(5);

			@fbids=();
			@notes= (); 
			foreach my $id (@ids) {
				my($converted,$note)= getConvertedCoord(5,$id); 
				push(@fbids,$converted); push(@notes,$note); 
				}

	open(OOO,'>./converted_coordinates.tsv');
					for( my $n= 0; $n<= $#ids; $n++ ) { 
						my $coord1= $ids[$n];
						#my($chr,$start,$stop)= splitRange($coord1);
						my $coord2= $fbids[$n];
						$coord2= $chr.':?..?' unless $coord2=~/^\s*[^ :]+:[^?-]*\d/;
						print $coord1."  ==>  ".$coord2."\n";   
						print OOO $coord1."\t".$coord2."\n";   
						}
	close(OOO);
	
	exit;

#***********************************************************
#
#***********************************************************

sub setConversionArrays
{
	my $version= shift;
		my @sss= ();
		open(FTS,"/bin/cat ./5to6_Dmel_mapfile.out |") || die "Cannot open conversion data file..";
		while( (my $str=<FTS>) ) {
			my($chr,$start,$end,$newchr,$newstart,$newend,$strand,@rest)= split(/[ \t]+/,$str);
			($start,$end)= ($end,$start) if $start>$end;
			($newstart,$newend)= ($newend,$newstart) if $newstart>$newend;
			if( $strand=~/\-/ ) {
				push(@sss, $chr . "\t" . $start . "\t" . $newend . "\t" . ($end-$start) . "\t" . $newchr . "\t" . $strand); }
			else {
				push(@sss, $chr . "\t" . $start . "\t" . $newstart . "\t" . ($end-$start) . "\t" . $newchr . "\n"); }
			}
		close(FTS);
		my $inn= "";
		my $cnt=1;
		open(FTS,'>./t.ttt') if $debug;
		foreach my $str ( sort mapfile @sss ) { # sorting (GFF-like) is critical
			$inn.= $str; 
			print $str."<br>" if $debug && $cnt++<20; 
			print FTS $str if $debug;
			}
		close(FTS) if $debug;

	my %Mapper= ();
	my @strings= split(/\n/,$inn);
	my $oldkey= 'none';
	my $nstrings= 0;
	foreach my $string ( @strings ) { # must be already sorted GFF-like
		my($chr,$start,$newstartornewend,$l,@rest)= split(/[ \t]+/,$string);
		push(@{$Mapper{$version.'-'.$chr.'-start'}},$start); # start of region on old
		push(@{$Mapper{$version.'-'.$chr.'-end'}},$start+$l); # end of region on old
		my $strand= (@rest>1) ? $rest[1] : '+';
		if( $strand=~/\-/ ) { 
			push(@{$Mapper{$version.'-'.$chr.'-diff'}},$start+$newstartornewend); } # what to subtract to convert
		else {
			push(@{$Mapper{$version.'-'.$chr.'-diff'}},$start-$newstartornewend); } # what to subtract to convert
		push(@{$Mapper{$version.'-'.$chr.'-newchr'}},(@rest>0) ? $rest[0] : $chr); # chromosome (might be new)
		$nstrings++;
		}
	print "Configured $nstrings stretches<br>\n" if $debug;
	return(\%Mapper);
}

#***********************************************************
#
#***********************************************************

sub getConvertedCoord
{
	my $version= shift;
	my $coord= shift;
	my($chr,$start,$stop)= splitRange($coord);
	return($chr,' failed to convert') unless defined $chr;
	($start,$stop)= ($stop,$start) if $start>$stop;
	$IsInversion= 0; # !!GLOBAL
	my($newchr,$newchr1,$newstart,$newstop,$shifts,$a,$b,$c)= ('?','?','?','?',"","","","");
	($newchr,$newstart)= convertPosition($version,$chr,$start);
	print "  vars state: $newchr:$newstart..<br>\n" if $debug;
	($newchr1,$newstop)= convertPosition($version,$chr,$stop);
	print "  vars state: $newchr:..$newstop<br>\n" if $debug;
	return("-",' failed to convert: maps to different scaffolds') if $newchr ne $newchr1;
	($c,$shifts,$a,$b)= convertRange($version,$chr,$start,$stop);
	($newstart,$newstop)= ($newstop,$newstart) if $newstart>$newstop;
	print "  vars state: $newchr:$newstart..$newstop<br>\n" if $debug;
	my $converted= $newchr.':'.$newstart.(($SingleCoordInput)?"":"..$newstop");
	my $note= " "; # should not be null
	if( $newstart eq '?' || $newstop eq '?' ) { $note= "in area of non-identity" ; }
	elsif( $shifts==1 ) { $note= "includes 1 area of non-identity"; }
	elsif( $shifts>1 ) { $note= "includes $shifts area(s) of non-identity"; }
	if( $newchr ne $chr ) { $note.= "; " unless $note eq " "; $note.= "different scaffold" ; }
	if( $IsInversion ) { $note.= "; " unless $note eq " "; $note.= "inversion" ; }
	print "  conversion result: $coord => $converted $note<br>\n" if $debug;
	return($converted,$note);
}

#***********************************************************
#
#***********************************************************

sub convertRange
# inside of stretch only!
{
	my($version,$chr,$start,$end)= @_;
	print "Decoding range $chr:$start..$end (v$version)<br>\n" if $debug;
	($start,$end)= ($end,$start) if $start>$end;
	if( $chr ne $ChrMapperKey ) {
		print "Setting arrays for $chr (v$version)<br>\n" if $debug;
		@ChrMapperStart= @{$ChrMapper->{$version.'-'.$chr.'-start'}};
		@ChrMapperEnd= @{$ChrMapper->{$version.'-'.$chr.'-end'}};
		@ChrMapperDiff= @{$ChrMapper->{$version.'-'.$chr.'-diff'}};
		@ChrMapperNewchr= @{$ChrMapper->{$version.'-'.$chr.'-newchr'}};
		$ChrMapperKey= $chr;
		} 
	my $nstretches= $#ChrMapperStart;
	print "Scanning array of $nstretches stretches<br>\n" if $debug;
	my $nbad= 0;
	for( my $n=0; $n<=$nstretches; $n++ ) {
		print "  testing stretch ".$ChrMapperStart[$n]."..".$ChrMapperEnd[$n]."<br>\n" if $debug;
		next if( $ChrMapperStart[$n]>$end );
		next if( $ChrMapperEnd[$n]<$start );
		# if passed, then this stretch of old release somehow overlaps with requested range
		print "    partially fits<br>\n" if $debug;
		if( $start>=$ChrMapperStart[$n] && $end<=$ChrMapperEnd[$n] ) { # fully within stretch
			my $diff= $ChrMapperDiff[$n];
			my $npos1= $pos1 - $diff; if( $npos1<0 ) { $IsInversion= 1; $npos1= -$npos1; }
			my $npos2= $pos2 - $diff; if( $npos2<0 ) { $IsInversion= 1; $npos2= -$npos2; }
			($npos1,$npos2)= ($npos2,$npos1) if $npos1>$npos2;
			return($ChrMapperNewchr[$n],0,$npos1,$npos2); }
		else { $nbad++; }
		}
	return($chr,$nbad-1,'?','?');
}

#***********************************************************
#
#***********************************************************

sub convertPosition
{
	my($version,$chr,$pos)= @_;
	print "Decoding position $chr:$pos (v$version)<br>\n" if $debug;
	if( $chr ne $ChrMapperKey ) {
		print "Setting arrays for $chr (v$version)<br>\n" if $debug;
		@ChrMapperStart= @{$ChrMapper->{$version.'-'.$chr.'-start'}};
		@ChrMapperEnd= @{$ChrMapper->{$version.'-'.$chr.'-end'}};
		@ChrMapperDiff= @{$ChrMapper->{$version.'-'.$chr.'-diff'}};
		@ChrMapperNewchr= @{$ChrMapper->{$version.'-'.$chr.'-newchr'}};
		$ChrMapperKey= $chr;
		} 
	my $nstretches= $#ChrMapperStart;
	for( my $n=0; $n<=$nstretches; $n++ ) {
		print "  testing ".$ChrMapperStart[$n].'-'.$ChrMapperEnd[$n]."<br>\n" if $debug;
		next if( $ChrMapperEnd[$n]<$pos );
		last if( $ChrMapperStart[$n]>$pos );
		my $newpos= $pos - $ChrMapperDiff[$n]; if( $newpos<0 ) { $IsInversion= 1; $newpos= -$newpos; }
		my $newchr= $ChrMapperNewchr[$n];
		print "$newchr:$pos => $chr:$newpos (fully in ".$ChrMapperStart[$n].'-'.$ChrMapperEnd[$n].")<br>\n" if $debug;
		return($newchr,$newpos);
		}
	return($chr,'?');
}

#***********************************************************
#
#***********************************************************

sub splitRange
{
	my $coords= shift;
	print "Splitting $coords<br>\n" if $debug;
	$coords=~s/,//g;
	if( $coords=~/^\s*([^ :]+):([0-9,]+)([\.\-]+)([0-9,]+)\s*$/ ) { # range
		$SingleCoordInput= 0; # !!GLOBAL
		my($chr,$start,$stop)= ($1,$2,$4);
		$chr=~s/^Chr//i;
		#if( $coords=~/[,]/ ) { $start=~s/[,]//; $stop=~s/[,]//; }
		if( $start>$stop ) { my $a= $stop; $stop= $start; $start= $a; }
		print "Split as $chr,$start,$stop<br>\n" if $debug;
		return($chr,$start,$stop);
		}
	elsif( $coords=~/^\s*([^ :]+):([0-9,]+)\s*$/ ) { # single coord
		$SingleCoordInput= 1; # !!GLOBAL
		my($chr,$start,$stop)= ($1,$2,$2);
		$chr=~s/^Chr//i;
		#if( $coords=~/[,]/ ) { $start=~s/[,]//; $stop=~s/[,]//; }
		if( $start>$stop ) { my $a= $stop; $stop= $start; $start= $a; }
		print "  Split as $chr $start $stop<br>\n" if $debug;
		return($chr,$start,$stop);
		}
	return;
}

#***********************************************************
#
#***********************************************************

sub mapfile
{
	my($a1,$a2,@ra)= split(/\t/,$a);
	my($b1,$b2,@rb)= split(/\t/,$b);
	if( $a1 ne $b1 ) { return( $a1 cmp $b1 ); }
	return( $a2 <=> $b2 );
}

