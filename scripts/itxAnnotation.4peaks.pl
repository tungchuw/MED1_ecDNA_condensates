$|++;
use strict;
use Data::Dumper;
#
# purpose: to annotate the ChIP-Seq peak for the interaction
# we will need 3 outputs category
# 1) both anchors with peaks
# 2) 1 anchor with peak
# 3) no anchor with peak
#
# awk -v FS="\t" -v OFS="\t" '{ if ($0!~/^#/) { if (""!=$1 && "chrL"!=$1 && "."!=$14) { print $1,$2,$3,$4,$5,$6,$7,$8,$9 } } else { print $0} }' Ezh2-all.cis.sigf.peaks.xls > Ezh2-all.cis.sigf.withpeaks.xls
#
#modify this block!!!!
my %G_EXECUTABLES=(
	'intersectBed'=>'singularity run --bind /net/nwgc/vol1/nobackup/nocleanup/tungch/:/net/nwgc/vol1/nobackup/nocleanup/tungch/ /net/nwgc/vol1/nobackup/nocleanup/tungch/sifsh/chipseqtools.sif bedtools intersect'
	);
my $G_ReportTranscriptBioType = 0;

my ($itxFile, $annotationFile, $peak1File, $peak2File, $peak3File, $peak4File, $factor) = @ARGV;

my $marker2 = 'ecDNA';
my $marker3 = 'SupEnh';
my $marker4 = 'MED1';
print "Input file: $itxFile\n";
print "Annotation file: $annotationFile\n";
print "Peak file: $peak1File\n";
print "$marker2 region file: $peak2File\n";
print "$marker3 region file: $peak3File\n";
print "$marker4 region file: $peak4File\n";

if (not defined $itxFile ) {
  die "Need significant interaction file at #1\n";
}
if (not defined $annotationFile) {
  die "Need annotation file at #2\n";
}
if (not defined $peak1File ) {
  die "Need Pol2 peak file at #3\n";
}
if (not defined $peak2File ) {
  die "Need ecDNA regions at #4\n";
}
if (not defined $peak3File ) {
  die "Need SuperEnhancer file at #5\n";
}
if (not defined $peak4File ) {
  die "Need MED1 peak at #6\n";
}
if (not defined $factor) {
  die "Please specify RUN name at #7\n";
}

#-------------------------------------------------------------------
#
#my $file = sprintf("%s/%s.%s", $datadir, $factor, $chiasigSuf);
my $outfile = sprintf("%s.annotated_itx.txt", $factor);

my %peaks1 = (
$factor => {factor1=>$factor, file=>sprintf('%s', $peak1File)}
);

my %peaks2 = (
$factor => {factor2=>$factor, file=>sprintf('%s', $peak2File)}
);

my %peaks3 = (
$factor => {factor3=>$factor, file=>sprintf('%s', $peak3File)}
);

my %peaks4 = (
$factor => {factor4=>$factor, file=>sprintf('%s', $peak4File)}
);

my $anchorsBed = sprintf("%s.bedpe.anchors.bed", $factor);
my @commands = ('awk', '-v', 'FS="\t"', '-v', 'OFS="\t"', '\'{if ($0!~/^#/) { print $1,$2,$3,"AL_"NR"\n"$4,$5,$6,"AR_"NR} }\'', $itxFile, '>', $anchorsBed);
my $command = join(' ', @commands);
print STDERR "Generating anchor file $anchorsBed\n";
0==system($command) || die "system @commands failed: $?\n";
print STDERR " done.\n";

open INFILE, $itxFile || die "Fail to open $itxFile\n$!\n";
open OUTFILE, ">$outfile" || die "Fail to open $outfile\n$!\n";
print STDERR "Reading $itxFile..";
my @rows = ();
while (<INFILE>) {
	if (/^#/) {
		print  ".";
	} else {
		chomp ();
		my @bits = split(/\t/);
		push @bits, $bits[2]-$bits[1];
		push @bits, $bits[5]-$bits[4];
		push @bits, int((($bits[5]+$bits[4])-($bits[2]+$bits[1]))/2);
		push @bits, $bits[0].':'.$bits[1].'-'.$bits[5];
		$rows[$.] = { 'interaction'=>\@bits, 'annotations'=>{} };
	}
}
close INFILE;
print STDERR " done.\n";

# inject additional comments
#print OUTFILE "# Legends (FinalCode): E=End; TES, G=Gene body; Exon + Intron, I=Intergenic, P=Promotor +/-2.5Kb, prefix 's'=shared common gene(s) between left and right anchor\n";
#print OUTFILE "# Legends (FullCode): S=overlap TSS +/-2.5Kb, E=Overlap exon(s), I=ovlerap intron(s), H=overlap TES+/-2.5Kb\n";


my @headers = ('chrL', 'startL', 'endL', 'chrR', 'startR', 'endR', 'ID', 'Nread', 'type', 'spanAL', 'spanAR', 'interactSpan', 'genome_coord');

#####
# let's get the peak1 binding information
if (1==1) {
	push @headers, 'Peak';
	
	my $anchorsIntersectBed1 = sprintf("%s.cis.sigf.anchors.peak1.intersect.bed", $factor);
	my @commands = ($G_EXECUTABLES{'intersectBed'}, '-wao', '-a', $anchorsBed, '-b', $peaks1{$factor}->{file}, '>', $anchorsIntersectBed1);
	my $command = join(' ', @commands);
	print STDERR "Checking overlap for $peak1File --> $anchorsIntersectBed1..";
	0==system($command) || die "system @commands failed: $?\n";
	print STDERR " done.\n";
	
	# work on the intersect results
	open INFILE, $anchorsIntersectBed1 || die "Fail to open $anchorsIntersectBed1\n$!\n";
	print STDERR "Processing peak $factor overlaps..";
	while (<INFILE>) {
		chomp();
		my @bits = split(/\t/);
		my $chr = shift @bits;
		my $start = shift @bits;
		my $end = shift @bits;
		my $name = shift @bits;
		my ($anchorType, $lineno) = split(/\_/, $name);
		$rows[$lineno]->{factors1}->{$factor} = {'AL'=>[], 'AR'=>[]} if (!exists $rows[$lineno]->{factors1}->{$factor});
		if ('.' ne $bits[3]) {
			my ($peak1) = $bits[3] =~ /(peak_\d+$)/;
			push @{$rows[$lineno]->{factors1}->{$factor}->{$anchorType}}, $peak1;
		}
	}
	close INFILE;
	print STDERR " done.\n";
}
# END - let's get the peak1 binding information
#####

#####
# let's get the peak2 binding information
if (1==1) {
	push @headers, sprintf("%s", $marker2);
	
	my $anchorsIntersectBed2 = sprintf("%s.cis.sigf.anchors.peak2.intersect.bed", $factor);
	my @commands = ($G_EXECUTABLES{'intersectBed'}, '-wao', '-a', $anchorsBed, '-b', $peaks2{$factor}->{file}, '>', $anchorsIntersectBed2);
	my $command = join(' ', @commands);
	print STDERR "Checking overlap for $peak2File --> $anchorsIntersectBed2..";
	0==system($command) || die "system @commands failed: $?\n";
	print STDERR " done.\n";
	
	
	# work on the intersect results
	open INFILE, $anchorsIntersectBed2 || die "Fail to open $anchorsIntersectBed2\n$!\n";
	print STDERR "Processing $marker2 $factor overlaps..";
	
	while (<INFILE>) {
		chomp();
		my @bits = split(/\t/);
		my $chr = shift @bits;
		my $start = shift @bits;
		my $end = shift @bits;
		my $name = shift @bits;
		my ($anchorType, $lineno) = split(/\_/, $name);
		$rows[$lineno]->{factors2}->{$factor} = {'AL'=>[], 'AR'=>[]} if (!exists $rows[$lineno]->{factors2}->{$factor});
		if ('.' ne $bits[3]) {
			my ($peak2) = $bits[3] =~ /(peak_\d+$)/;
			push @{$rows[$lineno]->{factors2}->{$factor}->{$anchorType}}, $peak2;
		}
	}
	close INFILE;
	print STDERR " done.\n";
}
# END - let's get the peak2 binding information
#####

#####
# let's get the peak3 binding information
if (1==1) {
	push @headers, sprintf("%s", $marker3);
	
	my $anchorsIntersectBed3 = sprintf("%s.cis.sigf.anchors.peak3.intersect.bed", $factor);
	my @commands = ($G_EXECUTABLES{'intersectBed'}, '-wao', '-a', $anchorsBed, '-b', $peaks3{$factor}->{file}, '>', $anchorsIntersectBed3);
	my $command = join(' ', @commands);
	print STDERR "Checking overlap for $peak3File --> $anchorsIntersectBed3..";
	0==system($command) || die "system @commands failed: $?\n";
	print STDERR " done.\n";
	
	
	# work on the intersect results
	open INFILE, $anchorsIntersectBed3 || die "Fail to open $anchorsIntersectBed3\n$!\n";
	print STDERR "Processing $marker3 $factor overlaps..";
	
	while (<INFILE>) {
		chomp();
		my @bits = split(/\t/);
		my $chr = shift @bits;
		my $start = shift @bits;
		my $end = shift @bits;
		my $name = shift @bits;
		my ($anchorType, $lineno) = split(/\_/, $name);
		$rows[$lineno]->{factors3}->{$factor} = {'AL'=>[], 'AR'=>[]} if (!exists $rows[$lineno]->{factors3}->{$factor});
		if ('.' ne $bits[3]) {
			my ($peak3) = $bits[3] =~ /(peak_\d+$)/;
			push @{$rows[$lineno]->{factors3}->{$factor}->{$anchorType}}, $peak3;
		}
	}
	close INFILE;
	print STDERR " done.\n";
}
# END - let's get the peak3 binding information
#####
#
# let's get the peak4 binding information
if (1==1) {
	push @headers, sprintf("%s", $marker4);
	
	my $anchorsIntersectBed4 = sprintf("%s.cis.sigf.anchors.peak4.intersect.bed", $factor);
	my @commands = ($G_EXECUTABLES{'intersectBed'}, '-wao', '-a', $anchorsBed, '-b', $peaks4{$factor}->{file}, '>', $anchorsIntersectBed4);
	my $command = join(' ', @commands);
	print STDERR "Checking overlap for $peak4File --> $anchorsIntersectBed4..";
	0==system($command) || die "system @commands failed: $?\n";
	print STDERR " done.\n";
	
	
	# work on the intersect results
	open INFILE, $anchorsIntersectBed4 || die "Fail to open $anchorsIntersectBed4\n$!\n";
	print STDERR "Processing $marker4 $factor overlaps..";
	
	while (<INFILE>) {
		chomp();
		my @bits = split(/\t/);
		my $chr = shift @bits;
		my $start = shift @bits;
		my $end = shift @bits;
		my $name = shift @bits;
		my ($anchorType, $lineno) = split(/\_/, $name);
		$rows[$lineno]->{factors4}->{$factor} = {'AL'=>[], 'AR'=>[]} if (!exists $rows[$lineno]->{factors4}->{$factor});
		if ('.' ne $bits[3]) {
			my ($peak4) = $bits[3] =~ /(peak_\d+$)/;
			push @{$rows[$lineno]->{factors4}->{$factor}->{$anchorType}}, $peak4;
		}
	}
	close INFILE;
	print STDERR " done.\n";
}
# END - let's get the peak3 binding information
#####
#####


# let's process
# intersectBed -wao -a Ezh2-all.cis.sigf.anchors.bed -b mm10-2015-11-02-UCSCGenes.piece.sorted.bed > Ezh2-all.cis.sigf.anchors.gene.bed
my $anchorsAnnotationBed = sprintf("%s.bedpe.anchors.gene.bed", $factor);
my @commands = ($G_EXECUTABLES{'intersectBed'}, '-wao', '-a', $anchorsBed, '-b', $annotationFile, '>', $anchorsAnnotationBed);
my $command = join(' ', @commands);
print STDERR "Checking overlap with gene annotation --> $anchorsAnnotationBed..";
0==system($command) || die "system @commands failed: $?\n";
print STDERR " done.\n";

# this is for gene annotation
my %pieceCodeToCategory = ('TSS'=>'TSS', 'EX'=>'Exon', 'IT'=>'Intron', 'TES'=>'TES');
my @pieceCategories = ('TSS', 'Exon', 'Intron', 'TES');
push @headers, 'FinalCode', 'GeneL', 'GeneR';
# WCH: 2018-09-07 Harianto wants gene's biotype to be included to remove separate steps
push @headers, 'LGeneType', 'RGeneType', 'FinalCodeL', 'FinalCodeR';
# WCH: 2018-09-07 As gene and transcript biotype can disagree, we are reporting transcript level information
push @headers, 'TranscriptL', 'TranscriptR', 'LTranscriptType', 'RTranscriptType' if (0!=$G_ReportTranscriptBioType);
push @headers, 'FullCodeL', 'FullCodeR';
push @headers, 'NO_common', 'common';
foreach my $pieceCategory (@pieceCategories) {
	push @headers, 'NO_'.$pieceCategory.'_L', $pieceCategory.'_L';
	push @headers, 'NO_'.$pieceCategory.'_R', $pieceCategory.'_R';
}
my %G_pieceCategory2Code = ('TSS'=>'S', 'Exon'=>'E', 'Intron'=>'I', 'TES'=>'H');
my %G_finalCodes = ();
# WCH: this will not consider TES as gene rather than TES
%G_finalCodes = (
'GI'=>'GI', 'IG'=>'GI'
,'PI'=>'PI', 'IP'=>'PI'
,'PG'=>'PG', 'GP'=>'PG'
,'II'=>'II', 'PP'=>'PP', 'GG'=>'GG'
,'sPP'=>'sPP','sGG'=>'sGG','sPG'=>'sPG','sGP'=>'sPG'
,'PE'=>'PG', 'EP'=>'PG'
,'GE'=>'GG', 'EG'=>'GG'
,'EI'=>'GI', 'IE'=>'GI'
,'EE'=>'GG'
,'sEE'=>'sGG'
,'sEG'=>'sGG','sGE'=>'sGG'
,'sPE'=>'sPG','sEP'=>'sPG'
);

# work on the intersect results
open INFILE, $anchorsAnnotationBed || die "Fail to open $$anchorsAnnotationBed\n$!\n";
print STDERR "Processing gene annotation overlaps..";

# chunk processing
#print OUTFILE "\n";
print OUTFILE join("\t", @headers), "\n";
my $currLineNo = -1;
my %interaction = ();
my %peak1_bindingCounts = ('.'=>0, 'L'=>0, 'R'=>0, 'LR'=>0);
my %peak2_bindingCounts = ('.'=>0, 'L'=>0, 'R'=>0, 'LR'=>0);
my %peak3_bindingCounts = ('.'=>0, 'L'=>0, 'R'=>0, 'LR'=>0);
my %peak4_bindingCounts = ('.'=>0, 'L'=>0, 'R'=>0, 'LR'=>0);
grep { $interaction{$_} = {'AL'=>[],'AR'=>[]} } @pieceCategories;
while (<INFILE>) {
	chomp();
	my @bits = split(/\t/);
	my $chr = shift @bits;
	my $start = shift @bits;
	my $end = shift @bits;
	my $name = shift @bits;
	my ($anchorType, $lineno) = split(/\_/, $name);
	
	if ($lineno!=$currLineNo) {
		if (-1!=$currLineNo) {
			# trigger to report if we are done accummulating for an interaction
			
			# report the peak1 binding information
			my @peak1_bindings = ();
			my $peak1_bindingsRef = $rows[$currLineNo]->{factors1}->{$factor};
			push @peak1_bindings, 'L' if (scalar(@{$peak1_bindingsRef->{'AL'}}));
			push @peak1_bindings, 'R' if (scalar(@{$peak1_bindingsRef->{'AR'}}));
			my $peak1_binding = (0==scalar(@peak1_bindings)) ? '.' : join("", @peak1_bindings);
			$peak1_bindingCounts{$peak1_binding}++;
				
			# report the peak2 binding information
			my @peak2_bindings = ();
			my $peak2_bindingsRef = $rows[$currLineNo]->{factors2}->{$factor};
			push @peak2_bindings, 'L' if (scalar(@{$peak2_bindingsRef->{'AL'}}));
			push @peak2_bindings, 'R' if (scalar(@{$peak2_bindingsRef->{'AR'}}));
			my $peak2_binding = (0==scalar(@peak2_bindings)) ? '.' : join("", @peak2_bindings);
			$peak2_bindingCounts{$peak2_binding}++;
			
			# report the peak3 binding information
			my @peak3_bindings = ();
			my $peak3_bindingsRef = $rows[$currLineNo]->{factors3}->{$factor};
			push @peak3_bindings, 'L' if (scalar(@{$peak3_bindingsRef->{'AL'}}));
			push @peak3_bindings, 'R' if (scalar(@{$peak3_bindingsRef->{'AR'}}));
			my $peak3_binding = (0==scalar(@peak3_bindings)) ? '.' : join("", @peak3_bindings);
			$peak3_bindingCounts{$peak3_binding}++;
			
			# report the peak4 binding information
			my @peak4_bindings = ();
			my $peak4_bindingsRef = $rows[$currLineNo]->{factors4}->{$factor};
			push @peak4_bindings, 'L' if (scalar(@{$peak4_bindingsRef->{'AL'}}));
			push @peak4_bindings, 'R' if (scalar(@{$peak4_bindingsRef->{'AR'}}));
			my $peak4_binding = (0==scalar(@peak4_bindings)) ? '.' : join("", @peak4_bindings);
			$peak4_bindingCounts{$peak4_binding}++;
			
				print OUTFILE join("\t", @{$rows[$currLineNo]->{interaction}});
				print OUTFILE "\t", $peak1_binding;
				print OUTFILE "\t", $peak2_binding;
				print OUTFILE "\t", $peak3_binding;
				print OUTFILE "\t", $peak4_binding;

				# report the conclusion
				my @fullcodesAL = ();
				my @fullcodesAR = ();
				my %countAL = ();
				my %countAR = ();
				my @genesLBiotypes = ();
				my @genesRBiotypes = ();
				foreach my $pieceCategory (@pieceCategories) {
					$countAL{$pieceCategory} = scalar(@{$interaction{$pieceCategory}->{'AL'}});
					if ($countAL{$pieceCategory}>0) {
						push @fullcodesAL, $G_pieceCategory2Code{$pieceCategory};
						foreach my $name (@{$interaction{$pieceCategory}->{'AL'}}) {
							my @nameBits = split(/,/, $name);
							push @genesLBiotypes, {gene=>$nameBits[0], gene_type=>$nameBits[3], transcript=>$nameBits[1], transcript_type=>$nameBits[4]};
						}
					}
					$countAR{$pieceCategory} = scalar(@{$interaction{$pieceCategory}->{'AR'}});
					if ($countAR{$pieceCategory}>0) {
						push @fullcodesAR, $G_pieceCategory2Code{$pieceCategory};
						foreach my $name (@{$interaction{$pieceCategory}->{'AR'}}) {
							my @nameBits = split(/,/, $name);
							push @genesRBiotypes, {gene=>$nameBits[0], gene_type=>$nameBits[3], transcript=>$nameBits[1], transcript_type=>$nameBits[4]};
						}
					}
				}
				my $fullcodeAL = '.';
				my $fullcodeAR = '.';
				my $finalcodeAL = '.';
				my $finalcodeAR = '.';
				my $finalcodeL = 'I';
				my $finalcodeR = 'I';
				if (scalar(@fullcodesAL)>0) {
					$fullcodeAL = join('', @fullcodesAL);
					# decide on the final code
					$finalcodeL = prioritizeCode(\%countAL);
				}
				if (scalar(@fullcodesAR)>0) {
					$fullcodeAR = join('', @fullcodesAR);
					# decide on the final code
					$finalcodeR = prioritizeCode(\%countAR);
				}
				my $finalcode = $finalcodeL.$finalcodeR;
				
				my @sPP = ();
				my %genesL = (); grep { $genesL{$_->{gene}}++; } @genesLBiotypes;
				my %genesR = (); grep { $genesR{$_->{gene}}++; } @genesRBiotypes;
				if ('I' ne $finalcodeL && 'I' ne $finalcodeR) {
					my @commons = ();
					foreach my $gene (keys %genesL) {
						if (exists $genesR{$gene}) {
							push @commons, $gene;
						}
					}
					my $numCommonds = scalar(@commons);
					if ($numCommonds>0) {
						$finalcode = 's'.$finalcode;
						push @sPP, $numCommonds, join(";", @commons);
					} else {
						push @sPP, '.', '.';
					}
					
				} else {
					push @sPP, '.', '.';
				}
				my @genesL = (); my $genesLBiotypes = geneBioTypes(\@genesLBiotypes, \@genesL);
				my @genesR = (); my $genesRBiotypes = geneBioTypes(\@genesRBiotypes, \@genesR);
				my @transcriptsL = (); my $transcriptsLBiotypes = ''; transcriptBioTypes(\@genesLBiotypes, \@transcriptsL) if (0!=$G_ReportTranscriptBioType);
				my @transcriptsR = (); my $transcriptsRBiotypes = ''; transcriptBioTypes(\@genesRBiotypes, \@transcriptsR) if (0!=$G_ReportTranscriptBioType);
				my @cols = ($G_finalCodes{$finalcode}, 0==scalar(@genesL) ? '.' : join(";", @genesL), 0==scalar(@genesR) ? '.' : join(";", @genesR));
				# WCH: 2018-09-07 add biotype columns
				push @cols, ('' eq $genesLBiotypes) ? '.' : $genesLBiotypes, ('' eq $genesRBiotypes) ? '.' : $genesRBiotypes;
				push @cols, $finalcodeL, $finalcodeR;
				if (0!=$G_ReportTranscriptBioType) {
					# WCH: 2018-09-07 add biotype columns
					push @cols, (0==scalar(@transcriptsL)) ? '.' : join(";", @transcriptsL), (0==scalar(@transcriptsR)) ? '.' : join(";", @transcriptsR);
					push @cols, ('' eq $transcriptsLBiotypes) ? '.' : $transcriptsLBiotypes, ('' eq $transcriptsRBiotypes) ? '.' : $transcriptsRBiotypes;
				}
				push @cols, $fullcodeAL, $fullcodeAR;
				print OUTFILE "\t", join("\t", @cols);
				print OUTFILE "\t", join("\t", @sPP);
				
				# report the individual piece-wise information
				foreach my $pieceCategory (@pieceCategories) {
					# WCH: 2018-09-07 remove the suffix gene and transcript biotype to reduce storage
					if ($countAL{$pieceCategory}>0) {
						print OUTFILE "\t", $countAL{$pieceCategory};
						my @elements = (); trimBioType(\@{$interaction{$pieceCategory}->{'AL'}}, \@elements);
						print OUTFILE "\t", join(";", @elements);
					} else {
						print OUTFILE "\t.\t.";
					}
					if ($countAR{$pieceCategory}>0) {
						print OUTFILE "\t", $countAR{$pieceCategory};
						my @elements = (); trimBioType(\@{$interaction{$pieceCategory}->{'AR'}}, \@elements);
						print OUTFILE "\t", join(";", @elements);
					} else {
						print OUTFILE "\t.\t.";
					}
				}
				
				print OUTFILE "\n";
			}
		%interaction = ();
		grep { $interaction{$_} = {'AL'=>[],'AR'=>[]} } @pieceCategories;
		$currLineNo = $lineno;
	}
	
	if ('.' ne $bits[3]) {
		#TSS:Xkr4:uc007aeu.1:1/3
		#EX:Xkr4:uc007aeu.1:1/3
		#IT:Sox17:uc007aey.1:1/1
		#TES:Sox17:uc033fhy.1:4/4
		my @nameBits = split(/:/, $bits[3]);
		my $pieceCategory = shift @nameBits;
		push @{$interaction{$pieceCodeToCategory{$pieceCategory}}->{$anchorType}}, join(',', @nameBits);
	}
	
}
close INFILE;
if (-1!=$currLineNo) {
	# trigger to report if we are done accummulating for an interaction
	
	my @peak1_bindings = ();
	my $peak1_bindingsRef = $rows[$currLineNo]->{factors1}->{$factor};
	push @peak1_bindings, 'L' if (scalar(@{$peak1_bindingsRef->{'AL'}}));
	push @peak1_bindings, 'R' if (scalar(@{$peak1_bindingsRef->{'AR'}}));
	my $peak1_binding = (0==scalar(@peak1_bindings)) ? '.' : join("", @peak1_bindings);
	$peak1_bindingCounts{$peak1_binding}++;

	my @peak2_bindings = ();
	my $peak2_bindingsRef = $rows[$currLineNo]->{factors2}->{$factor};
	push @peak2_bindings, 'L' if (scalar(@{$peak2_bindingsRef->{'AL'}}));
	push @peak2_bindings, 'R' if (scalar(@{$peak2_bindingsRef->{'AR'}}));
	my $peak2_binding = (0==scalar(@peak2_bindings)) ? '.' : join("", @peak2_bindings);
	$peak2_bindingCounts{$peak2_binding}++;

	my @peak3_bindings = ();
	my $peak3_bindingsRef = $rows[$currLineNo]->{factors3}->{$factor};
	push @peak3_bindings, 'L' if (scalar(@{$peak3_bindingsRef->{'AL'}}));
	push @peak3_bindings, 'R' if (scalar(@{$peak3_bindingsRef->{'AR'}}));
	my $peak3_binding = (0==scalar(@peak3_bindings)) ? '.' : join("", @peak3_bindings);
	$peak3_bindingCounts{$peak3_binding}++;

	my @peak4_bindings = ();
	my $peak4_bindingsRef = $rows[$currLineNo]->{factors4}->{$factor};
	push @peak4_bindings, 'L' if (scalar(@{$peak4_bindingsRef->{'AL'}}));
	push @peak4_bindings, 'R' if (scalar(@{$peak4_bindingsRef->{'AR'}}));
	my $peak4_binding = (0==scalar(@peak4_bindings)) ? '.' : join("", @peak4_bindings);
	$peak4_bindingCounts{$peak4_binding}++;

		print OUTFILE join("\t", @{$rows[$currLineNo]->{interaction}});
		print OUTFILE "\t", $peak1_binding;
		print OUTFILE "\t", $peak2_binding;
		print OUTFILE "\t", $peak3_binding;
		print OUTFILE "\t", $peak4_binding;
		
		# report the conclusion
		my @fullcodesAL = ();
		my @fullcodesAR = ();
		my %countAL = ();
		my %countAR = ();
		my @genesLBiotypes = ();
		my @genesRBiotypes = ();
		foreach my $pieceCategory (@pieceCategories) {
			$countAL{$pieceCategory} = scalar(@{$interaction{$pieceCategory}->{'AL'}});
			if ($countAL{$pieceCategory}>0) {
				push @fullcodesAL, $G_pieceCategory2Code{$pieceCategory};
				foreach my $name (@{$interaction{$pieceCategory}->{'AL'}}) {
					my @nameBits = split(/,/, $name);
					push @genesLBiotypes, {gene=>$nameBits[0], gene_type=>$nameBits[3], transcript=>$nameBits[1], transcript_type=>$nameBits[4]};
				}
			}
			$countAR{$pieceCategory} = scalar(@{$interaction{$pieceCategory}->{'AR'}});
			if ($countAR{$pieceCategory}>0) {
				push @fullcodesAR, $G_pieceCategory2Code{$pieceCategory};
				foreach my $name (@{$interaction{$pieceCategory}->{'AR'}}) {
					my @nameBits = split(/,/, $name);
					push @genesRBiotypes, {gene=>$nameBits[0], gene_type=>$nameBits[3], transcript=>$nameBits[1], transcript_type=>$nameBits[4]};
				}
			}
		}
		my $fullcodeAL = '.';
		my $fullcodeAR = '.';
		my $finalcodeAL = '.';
		my $finalcodeAR = '.';
		my $finalcodeL = 'I';
		my $finalcodeR = 'I';
		if (scalar(@fullcodesAL)>0) {
			$fullcodeAL = join('', @fullcodesAL);
			# decide on the final code
			$finalcodeL = prioritizeCode(\%countAL);
		}
		if (scalar(@fullcodesAR)>0) {
			$fullcodeAR = join('', @fullcodesAR);
			# decide on the final code
			$finalcodeR = prioritizeCode(\%countAR);
		}
		my $finalcode = $finalcodeL.$finalcodeR;
		
		my @sPP = ();
		my %genesL = (); grep { $genesL{$_->{gene}}++; } @genesLBiotypes;
		my %genesR = (); grep { $genesR{$_->{gene}}++; } @genesRBiotypes;
		if ('I' ne $finalcodeL && 'I' ne $finalcodeR) {
			# TODO: test if within the same promoter/gene
			
			#my %tss_ALs = ();
			#grep { $tss_ALs{$_} = $_ } @{$interaction{'TSS'}->{'AL'}};
			#my @commons = ();
			#grep { push @commons, $_ if (exists $tss_ALs{$_}); } @{$interaction{'TSS'}->{'AR'}};
			#my $numCommonds = scalar(@commons);
			#if ($numCommonds>0) {
			#	$finalcode = 'sPP';
			#	push @sPP, $numCommonds, join(";", @commons);
			#} else {
			#	push @sPP, '.', '.';
			#}
			
			my @commons = ();
			foreach my $gene (keys %genesL) {
				if (exists $genesR{$gene}) {
					push @commons, $gene;
				}
			}
			my $numCommonds = scalar(@commons);
			if ($numCommonds>0) {
				$finalcode = 's'.$finalcode;
				push @sPP, $numCommonds, join(";", @commons);
			} else {
				push @sPP, '.', '.';
			}
			
		} else {
			push @sPP, '.', '.';
		}
		my @genesL = (); my $genesLBiotypes = geneBioTypes(\@genesLBiotypes, \@genesL);
		my @genesR = (); my $genesRBiotypes = geneBioTypes(\@genesRBiotypes, \@genesR);
		my @transcriptsL = (); my $transcriptsLBiotypes = ''; transcriptBioTypes(\@genesLBiotypes, \@transcriptsL) if (0!=$G_ReportTranscriptBioType);
		my @transcriptsR = (); my $transcriptsRBiotypes = ''; transcriptBioTypes(\@genesRBiotypes, \@transcriptsR) if (0!=$G_ReportTranscriptBioType);
		my @cols = ($G_finalCodes{$finalcode}, 0==scalar(@genesL) ? '.' : join(";", @genesL), 0==scalar(@genesR) ? '.' : join(";", @genesR));
		# WCH: 2018-09-07 add biotype columns
		push @cols, ('' eq $genesLBiotypes) ? '.' : $genesLBiotypes, ('' eq $genesRBiotypes) ? '.' : $genesRBiotypes;
		push @cols, $finalcodeL, $finalcodeR;
		if (0!=$G_ReportTranscriptBioType) {
			# WCH: 2018-09-07 add biotype columns
			push @cols, (0==scalar(@transcriptsL)) ? '.' : join(";", @transcriptsL), (0==scalar(@transcriptsR)) ? '.' : join(";", @transcriptsR);
			push @cols, ('' eq $transcriptsLBiotypes) ? '.' : $transcriptsLBiotypes, ('' eq $transcriptsRBiotypes) ? '.' : $transcriptsRBiotypes;
		}
		push @cols, $fullcodeAL, $fullcodeAR;
		print OUTFILE "\t", join("\t", @cols);
		print OUTFILE "\t", join("\t", @sPP);
		
		# report the individual piece-wise information
		foreach my $pieceCategory (@pieceCategories) {
			# TODO WCH: 2018-09-07 remove the suffix gene and transcript biotype to reduce storage
			if ($countAL{$pieceCategory}>0) {
				print OUTFILE "\t", $countAL{$pieceCategory};
				my @elements = (); trimBioType(\@{$interaction{$pieceCategory}->{'AL'}}, \@elements);
				print OUTFILE "\t", join(";", @elements);
			} else {
				print OUTFILE "\t.\t.";
			}
			if ($countAR{$pieceCategory}>0) {
				print OUTFILE "\t", $countAR{$pieceCategory};
				my @elements = (); trimBioType(\@{$interaction{$pieceCategory}->{'AR'}}, \@elements);
				print OUTFILE "\t", join(";", @elements);
			} else {
				print OUTFILE "\t.\t.";
			}
		}
		
		print OUTFILE "\n";
	}

print STDERR " done.\n";
close OUTFILE;
print STDERR "Peak_Binding\tCount\n";
foreach my $peak1_binding ('.', 'L', 'R', 'LR') {
	print STDERR $peak1_binding, "\t", $peak1_bindingCounts{$peak1_binding}, "\n";
}
print STDERR "ecDNA\tCount\n";
foreach my $peak2_binding ('.', 'L', 'R', 'LR') {
	print STDERR $peak2_binding, "\t", $peak2_bindingCounts{$peak2_binding}, "\n";
}
print STDERR "SuperEnhancer\tCount\n";
foreach my $peak3_binding ('.', 'L', 'R', 'LR') {
	print STDERR $peak3_binding, "\t", $peak3_bindingCounts{$peak3_binding}, "\n";
}

print STDERR "MED1_binding\tCount\n";
foreach my $peak4_binding ('.', 'L', 'R', 'LR') {
	print STDERR $peak4_binding, "\t", $peak4_bindingCounts{$peak4_binding}, "\n";
}

print STDERR "Done.\n";
exit 0;

sub prioritizeCode {
	my ($countsRef) = @_;
	
	if ($countsRef->{'TSS'}>0) {
		return 'P';
	} elsif ($countsRef->{'TES'}>0) {
		return 'E';
	} elsif ($countsRef->{'Exon'}>0 || $countsRef->{'Intron'}>0
		#|| $countsRef->{'TES'}>0 # TODO: do we want to include this?
		) {
		return 'G';
	} else {
		return 'I';
	}
}

sub trimBioType {
	my ($featuresRef, $trimsRef) = @_;
	@{$trimsRef} = ();
	foreach my $feature (@{$featuresRef}) {
		my @bits = split(/,/, $feature);
		push @{$trimsRef}, join(',', @bits[0..$#bits-2]);
	}
}

sub getElementAndBioTypes {
	my ($elementTag, $elementTypeTag, $biotypesRef, $uniqElementsRef) = @_;

	@{$uniqElementsRef} = ();
	my @orders = sort { return ($a->{$elementTag} cmp $b->{$elementTag} || $a->{$elementTypeTag} cmp $b->{$elementTypeTag})} @{$biotypesRef};
	my @biotypes = ();
	my $prevKey = '';
	foreach my $biotypeRef (@orders) {
		my $key = $biotypeRef->{$elementTag} . ':' . $biotypeRef->{$elementTypeTag};
		if ($prevKey ne $key) {
			push @{$uniqElementsRef}, $biotypeRef->{$elementTag};
			push @biotypes, $biotypeRef->{$elementTypeTag};
			$prevKey = $key;
		}
	}
	return join(';', @biotypes);
}

sub geneBioTypes {
	my ($biotypesRef, $uniqElementsRef) = @_;
	return getElementAndBioTypes('gene', 'gene_type', $biotypesRef, $uniqElementsRef);
}

sub transcriptBioTypes {
	my ($biotypesRef, $uniqElementsRef) = @_;
	return getElementAndBioTypes('transcript', 'transcript_type', $biotypesRef, $uniqElementsRef);
}
