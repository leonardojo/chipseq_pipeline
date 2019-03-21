#!/usr/bin/perl
use strict;use warnings;
use List::Util qw[min max];


## This script aims to get replicated ChIP Datas (Processed BAM files) and their corresponding Input Controls (Processed BAM files) and perform:
## 1) Cross-Correlation Analysis;
## 2) Call Peaks in each independent Biorep with macs2 and Call Peaks for concatenated files -> Oracle File;
## 3) IDR analysis to find reproducible peaks in the oracle file(IDR =< 0.01);
## 4) Select the peaks in the oracle that are reproducible based on IDR;
## 5) Assigning Peaks to genes;
## 6) Create a chipseqsummary file;
## GENERAL Usage: perl chipseqpipeline.ljo_version.pl chip.br1.dedup.bam chip.br2.dedup.bam input.br1.dedup.bam input.br2.dedup.bam;
## OUPUT: 4 folders:
##	$samplename/CC/ -> tab and plot files for CC analysis
##  $samplename/macs2/ -> all peaks identified and summit files for BR1, BR2 and ORACLE
##  $samplename/IDR/ -> IDR output file with peaks in the oracle that are reproducible in BR1 and BR2 + Plot for IDR
##  $samplename/OUTPUT/ -> Final peak files (reproducible peaks, summits) + List of Bound genes from reproducible peaks + ChIP-Seq Summary file that contains the summary for the analysis

my $start_run = time();
my $bamfile_br1=$ARGV[0];
my $bamfile_br2=$ARGV[1];
my $bamfile_Input1=$ARGV[2];
my $bamfile_Input2=$ARGV[3];

## File names (Should BE STAGE-TF-PEPTIDE-BR#.DEDUP.BAM)
my $bamfile_br1_short=substr($bamfile_br1,0,-10);
my $bamfile_br2_short=substr($bamfile_br2,0,-10);
my $bamfile_Input1_short=substr($bamfile_Input1,0,-10);
my $bamfile_Input2_short=substr($bamfile_Input2,0,-10);
my $samplename=substr($bamfile_br1,0,-18); #LEC1.BR1.dedup.bam -> LEC1

system("rm -r $samplename; mkdir $samplename");
system("mkdir $samplename/macs2");
system("mkdir $samplename/IDR");
system("mkdir $samplename/CC");
system("mkdir $samplename/OUTPUT");

## 1-Cross-Correlation -- OUTPUT: TAB FILE CONTAINING CC Graph, and Filename.CC.txt containing the NSC and RSC values
system("Rscript /home/jpelletier/phantompeakqualtools/trunk/run_spp_nodups.R -c=$bamfile_br1 -savp -odir=$samplename/CC/ -out=$samplename/CC/$bamfile_br1_short.CC.tab");
system("Rscript /home/jpelletier/phantompeakqualtools/trunk/run_spp_nodups.R -c=$bamfile_br2 -savp -odir=$samplename/CC/ -out=$samplename/CC/$bamfile_br2_short.CC.tab");
system("head $samplename/CC/*.CC.tab > $samplename/CC/$samplename.CC.summary.txt"); 

## 2-Call Peaks with macs2
system("macs2 callpeak -t $bamfile_br1 -c $bamfile_Input1 $bamfile_Input2 -g 7.88e8 -n $bamfile_br1_short-catInp --outdir $samplename/macs2/ -p 1e-1 --keep-dup all");
system("macs2 callpeak -t $bamfile_br2 -c $bamfile_Input1 $bamfile_Input2 -g 7.88e8 -n $bamfile_br2_short-catInp --outdir $samplename/macs2/ -p 1e-1 --keep-dup all");
system("macs2 callpeak -t $bamfile_br1 $bamfile_br2  -c $bamfile_Input1 $bamfile_Input2 -g 7.88e8 -n $samplename-catIP-catInp --outdir $samplename/macs2/ -p 1e-1 --keep-dup all");

## 3-IDR Analysis
my $FILE_NAME_1="$samplename/macs2/$bamfile_br1_short-catInp_peaks.narrowPeak";
my $FILE_NAME_2="$samplename/macs2/$bamfile_br2_short-catInp_peaks.narrowPeak";
my $SORTED_NAME_1="$samplename/macs2/SortedPeak-$bamfile_br1_short-catInp_peaks.narrowPeak";
my $SORTED_NAME_2="$samplename/macs2/SortedPeak-$bamfile_br2_short-catInp_peaks.narrowPeak";
my $ORACLE_NAME="$samplename/macs2/$samplename-catIP-catInp_peaks.narrowPeak";
my $cnt1; open(FH,$FILE_NAME_1) or die "$!"; $cnt1++ while <FH>; close FH;
my $cnt2; open(FH,$FILE_NAME_2) or die "$!"; $cnt2++ while <FH>; close FH;
my $npeaks=min($cnt1,$cnt2);
system("sort -k 8nr,8nr $FILE_NAME_1 | head -n $npeaks > $samplename/macs2/SortedPeak-$bamfile_br1_short-catInp_peaks.narrowPeak; done");
system("sort -k 8nr,8nr $FILE_NAME_2 | head -n $npeaks > $samplename/macs2/SortedPeak-$bamfile_br2_short-catInp_peaks.narrowPeak; done");
my $cnt3; open(FH,$ORACLE_NAME) or die "$!"; $cnt3++ while <FH>; close FH;
system("idr --samples $SORTED_NAME_1 $SORTED_NAME_2 --peak-list $ORACLE_NAME --idr-threshold 0.01 --plot --output-file $samplename/IDR/$samplename-BR1_vs_BR2.bed >$samplename/IDR/$samplename-BR1_vs_BR2.summary.txt");

## 4-Assessing the REPRODUCIBLE peaks and summits in the ORACLE file
my $IDR="$samplename/IDR/$samplename-BR1_vs_BR2.bed";
my $summit="$samplename/macs2/$samplename-catIP-catInp_summits.bed";
my $ORACLE_REPRODUCIBLE_PEAKS="$samplename/OUTPUT/reproducible-$samplename-catIP-catInp_peaks.narrowPeak";
system("bedtools intersect -a $ORACLE_NAME -b $IDR -wa -nonamecheck | sort -k 8nr,8nr > $samplename/OUTPUT/reproducible-$samplename-catIP-catInp_peaks.narrowPeak");
system("bedtools intersect -a $summit -b $IDR -wa -nonamecheck | sort -k 8nr,8nr > $samplename/OUTPUT/reproducible-$samplename-catIP-catInp_summits.bed");

## 5-Assigning Peaks to Genes
system("bedtools window -l 1000 -r 0 -sw -a /home/ljo/Gmax2BED/Gmax_275_Wm82.a2.v1.TSS.bed -b $ORACLE_REPRODUCIBLE_PEAKS > $samplename/OUTPUT/boundgenes-reproducible-$samplename-catIP-catInp.bed"); 

## 6 - Final Summary file
my $end_run=time();
open (FH,"$samplename/CC/$bamfile_br1_short.CC.tab") or die "something went wrong $!";
        my $firstline1=<FH>;
        my @spl1 = split('\t',$firstline1);
close FH;
my $NSC_BR1= $spl1[8];
my $RSC_BR1= $spl1[9];

open (FH,"$samplename/CC/$bamfile_br2_short.CC.tab") or die "something went wrong $!";
        my $firstline2=<FH>;
        my @spl2 = split('\t',$firstline2);
close FH;
my $NSC_BR2= $spl2[8];
my $RSC_BR2= $spl2[9];

my $idrcnt; open(FH,$IDR) or die "$!"; $idrcnt++ while <FH>; close FH;
my $boundgenes; open(FH,"$samplename/OUTPUT/boundgenes-reproducible-$samplename-catIP-catInp.bed") or die "$!"; $boundgenes++ while <FH>; close FH;
my $macs2_br1="macs2 callpeak -t $bamfile_br1 -c $bamfile_Input1 $bamfile_Input2 -g 7.88e8 -n $bamfile_br1_short-catInp --outdir $samplename/macs2/ -p 1e-1 --keep-dup all";
my $macs2_br2="macs2 callpeak -t $bamfile_br2 -c $bamfile_Input1 $bamfile_Input2 -g 7.88e8 -n $bamfile_br1_short-catInp --outdir $samplename/macs2/ -p 1e-1 --keep-dup all";
my $macs2_oracle="macs2 callpeak -t $bamfile_br1 $bamfile_br2  -c $bamfile_Input1 $bamfile_Input2 -g 7.88e8 -n $samplename-catIP-catInp --outdir $samplename/macs2/ -p 1e-1 --keep-dup all";
my $end_run=time();
my $run_time = $end_run - $start_run;
my $SUMMARY="$samplename/OUTPUT/$samplename.chipsummary.txt" ;
open(my $fh, '>', "$samplename/OUTPUT/$samplename.ChIPseq.summary.txt") or die "Could not open file $!";
print $fh "$samplename ChIP-Seq Summary\nBiological Replicate Bam File 1 = $bamfile_br1\nBiological Replicate Bam File 2 = $bamfile_br2\nInput Bam File 1 = $bamfile_Input1\nInput Bam File 2 = $bamfile_Input2\nCross-Correlation NSC and RSC Biorep1 = $NSC_BR1 $RSC_BR1\nCross-Correlation NSC and RSC Biorep2 = $NSC_BR2 $RSC_BR2\nMacs2 command Biorep1 = $macs2_br1\nMacs2 command Biorep1 = $macs2_br2\nMacs2 command Biorep1 = $macs2_oracle\nNumber of peaks in Biorep1 after macs2 = $cnt1\nNumber of peaks in Biorep2 after macs2 = $cnt2\nAdjusting the total number of peaks to each biorep to = $npeaks\nNumber of peaks in Oracle after macs2 = $cnt3\nPeaks in Oracle are reproducible after IDR Analysis = $idrcnt\nNumber of boundgenes = $boundgenes\nTime to run the pipeline = $run_time\n";
close $fh;

print "done\n";





