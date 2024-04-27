#!/usr/bin/perl -w
use strict;

while (<>) {
    chomp;
    next if /^$/;
    my $info = $_;
    my ($path, $file)=('','');
    #    ($path, $file)=($1, $2) if $info=~/(.*)\/(.*)/;
    $file=$info;
    my $from = "/home/qiushi.li/periwinkle/job8-ovl1/4-polish";
    my $to = "/home/qiushi.li/periwinkle/gcpp/gcpp_files";
    
    `mkdir $to/$file && cp $from/quiver-run/$file/uow-00/aln-$file.bam $from/quiver-run/$file/uow-00/aln-$file.bam.bai $from/quiver-run/$file/uow-00/aln-$file.bam.pbi $from/quiver-split/refs/$file/ref.fasta $from/quiver-split/refs/$file/ref.fasta.fai $from/quiver-split/refs/$file/ctg_type $to/$file`;
	}
