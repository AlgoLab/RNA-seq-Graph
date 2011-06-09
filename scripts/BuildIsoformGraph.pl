#!/usr/bin/perl -w

####
#
# Script che produce il grafo delle isoforme (STDOUT) dal file GTF dei trascritti (STDIN) dal file
#
####

use Time::HiRes;
use strict;
use Getopt::Long;

our $debug = 0;
my $contract_path;
my $wdir;
my $gtfFile;
my $gen_file;
my $outIsoformFile;
my $outIsoformFileRegions;
my $outIsoformFileGDL;

GetOptions (
			'contract-path=i' => \$contract_path,
			'working-dir=s' => \$wdir,
            'gtf=s' => \$gtfFile,
            'genomic=s' => \$gen_file, 
            'isoforms=s' => \$outIsoformFile,
            'isoforms-regions=s' => \$outIsoformFileRegions,
            'isoforms-gdl=s' => \$outIsoformFileGDL,            
            'debug' => \$debug,
           );
my $usage="Usage: perl BuildIsoformGraph.pl [options]
 --contract-path= <integer> (mandatory) type of path contraction
 							(0=regions; 1=chunks; 2=RNA-read-like chunks)
 								(default: 1)
 --working-dir=<file>		 (optional) path of the working directory
 								(default: ./)
 --gtf= <file>          	(optional) name of the input GTF file
 								(default: exons.gtf)
 --genomic= <file>           (optional) name of the genomic sequence file (5'3')
 								(default: genomic.txt)
 --isoforms= <file>          (optional) name of the output Isoform Graph file
 								(default: isoform-graph.txt)
 --isoforms-regions= <file>  (optional) name of the output Isoform Graph file with regions
 								(default: isoform-graph-regions.txt)
 --isoforms-gdl= <file>  	(optional) name of the output Isoform Graph file in GDL format
 								(default: isoform-graph.gdl)
 --debug                     emits debugging information
";

die $usage unless (defined $contract_path and $contract_path >= 0 and $contract_path <= 2);

$wdir="./" if (not defined $wdir or $wdir eq '' or $wdir eq ".");
print "Working directory: ", $wdir, "\n";
chop $wdir if(substr($wdir, length($wdir)-1, 1) eq "/");

$gtfFile="exons.gtf" if (not defined $gtfFile or $gtfFile eq '');
print "GTF file: ", $gtfFile, "\n";

$gen_file="genomic.txt" if (not defined $gen_file or $gen_file eq '');
print "Genomic sequence file: ", $gen_file, "\n";

$outIsoformFile="isoform-graph.txt" if (not defined $outIsoformFile or $outIsoformFile eq '');
print "Output isoform graph file: ", $outIsoformFile, "\n";

$outIsoformFileRegions="isoform-graph-regions.txt" if (not defined $outIsoformFileRegions or $outIsoformFileRegions eq '');
print "Output isoform graph file with regions: ", $outIsoformFileRegions, "\n";

$outIsoformFileGDL="isoform-graph.gdl" if (not defined $outIsoformFileGDL or $outIsoformFileGDL eq '');
print "Output GDL isoform graph file: ", $outIsoformFileGDL, "\n";

print "Path contraction: ", $contract_path, "\n";

$gtfFile=$wdir."/".$gtfFile;
$gen_file=$wdir."/".$gen_file;
$outIsoformFile=$wdir."/".$outIsoformFile;
$outIsoformFileRegions=$wdir."/".$outIsoformFileRegions;
$outIsoformFileGDL=$wdir."/".$outIsoformFileGDL;

for my $f (($gen_file, $gtfFile)){
    die "File $f does not exists\n" unless (-f $f);
}

open OUT, ">$outIsoformFile" or die "Could not open $outIsoformFile: $!\n";
open OUTREG, ">$outIsoformFileRegions" or die "Could not open $outIsoformFileRegions: $!\n";
open OUTGDL, ">$outIsoformFileGDL" or die "Could not open $outIsoformFileGDL: $!\n";

open GEN, "<", $gen_file or die "Could not read $gen_file: $!\n";
my $header=<GEN>;
die "Genomic file $gen_file not in FASTA format!" unless($header =~ m/^>/);
my $gen_seq="";
while(<GEN>){
	chomp;
	$gen_seq.=$_;
}
close GEN;
my $gen_length=length($gen_seq);

open GTF, $gtfFile or die "Could not open $gtfFile: $!\n";

my $strand=1;
my %exon_hash=();	#key=exon_start, value=hash_ref1(key=exon_end, value=hash_ref2(key=transcript_ID, value=1))
while(<GTF>){
	chomp;
	my @record_list=split(/\s/, $_);
	
	if($record_list[2] eq "exon"){
		my $exon_start=$record_list[3];
		my $exon_end=$record_list[4];
		
		#Inversione per gene -
		if($record_list[6] eq "-"){
			$strand=-1;
			$exon_start=$gen_length-$record_list[4]+1;
			$exon_end=$gen_length-$record_list[3]+1;
		}
		else{
			$exon_start=$record_list[3];
			$exon_end=$record_list[4];
		}	
		
		my %exon_hash1;
		if(exists $exon_hash{$exon_start}){
			#print "\tEsiste start\n";
			%exon_hash1=%{$exon_hash{$exon_start}};
			my %exon_hash2;
			if(exists $exon_hash1{$exon_end}){
				#print "\tEsiste end\n";
				%exon_hash2=%{$exon_hash1{$exon_end}};
			}
			else{
				#print "\tNON Esiste end\n";
				%exon_hash2=();
			}
			$exon_hash2{$record_list[9]}=1;
			$exon_hash1{$exon_end}=\%exon_hash2;
			
		}
		else{
			#print "\tNON Esiste start\n";
			%exon_hash1=();
			my %exon_hash2=();
			$exon_hash2{$record_list[9]}=1;
			$exon_hash1{$exon_end}=\%exon_hash2;
		}
		$exon_hash{$exon_start}=\%exon_hash1;
	}
}

close GTF;

if($strand != 1){
	$gen_seq=reverse_complement($gen_seq);
}

my %B_hash=();	#key=position, value=0|1|2 (0:leftend, 1:rightend, 2:both)
foreach my $key(keys %exon_hash){
	if(exists $B_hash{$key}){
		if($B_hash{$key} == 1){
			$B_hash{$key}=2;
		}
	}
	else{
		$B_hash{$key}=0;
	}
	my %exon_hash1=%{$exon_hash{$key}};
	foreach my $key1(keys %exon_hash1){
		if(exists $B_hash{$key1}){
			if($B_hash{$key1} == 0){
				$B_hash{$key1}=2;
			}
		}
		else{
			$B_hash{$key1}=1;
		}
		
		#my %exon_hash2=%{$exon_hash1{$key1}};
		#print $key, "-", $key1, "\n";
		#my @keys_id=keys %exon_hash2;
		#print "\t".join(",", @keys_id), "\n";
	}
}

my @B_list=();	#[genomic_position, 0|1]	0:leftend, 1:rightend
foreach my $key(sort {$a <=> $b} (keys %B_hash)){
	if($B_hash{$key}==2){
		push @B_list, [$key, 0];
		push @B_list, [$key, 1];
	}
	else{
		push @B_list, [$key, $B_hash{$key}];
	}
	#print "POS=".$key." FLAG=".$B_hash{$key}."\n";
}

my @G_list=();	#[list_ref:list of regions, 0|1, 0|1, trs_hash_ref ,0|1:contracted, list_of_out_ref, list_of_in_ref]	list of no-variation regions
for(my $i=0; $i < $#B_list; $i++){
	my @pos1=@{$B_list[$i]};
	my @pos2=@{$B_list[$i+1]};
	#print $pos[0]."-".$pos[1]."\n";

	my @region=();
	my $add_region=0;
	#Region of type start-* ==> no-variation region
	if($pos1[1] == 0){
		$add_region=1;
		my @chunks=();
		push @chunks, $pos1[0];
		if($pos2[1] == 0){
			push @chunks, $pos2[0]-1;
		}
		else{
			push @chunks, $pos2[0];
		}
		push @region, \@chunks;
		push @region, $pos1[1], $pos2[1];
	}
	else{
		#Region of type end-end ==> no-variation region
		if($pos1[1] == $pos2[1]){
			$add_region=1;
			my @chunks=();
			push @chunks, $pos1[0]+1;
			push @chunks, $pos2[0];
			push @region, \@chunks;
			push @region, $pos1[1], $pos2[1];
		}
		else{
			#The region is not null
			if($pos2[0]-$pos1[0]-2 >= 0){
				my @exon_starts=sort {$a <=> $b} keys %exon_hash;
				my $i=0;
				my $found=0;
				while($i <= $#exon_starts && ($exon_starts[$i] <= $pos1[0] && $found == 0)){
					my %exon_hash1=%{$exon_hash{$exon_starts[$i]}};
					my @exon_ends=sort {$b <=> $a} keys %exon_hash1;
					if($exon_ends[0] >= $pos2[0]){
						$found=1;	
					}
					$i++;
				}
				if($found == 1){
					$add_region=1;
					my @chunks=();
					push @chunks, $pos1[0]+1;
					push @chunks, $pos2[0]-1;
					push @region, \@chunks;
					push @region, $pos1[1], $pos2[1];
				}
			}
		}
	}
	if($add_region == 1){
		my %region_trs=();
		my @exon_starts=sort {$a <=> $b} keys %exon_hash;
		my $i=0;
		while($i <= $#exon_starts && $exon_starts[$i] <= ${$region[0]}[0]){
			my %exon_hash1=%{$exon_hash{$exon_starts[$i]}};
			my @exon_ends=sort {$b <=> $a} keys %exon_hash1;
			my $j=0;
			while($j <= $#exon_ends && ($exon_ends[$j] >= ${$region[0]}[1])){
				my %exon_hash2=%{$exon_hash1{$exon_ends[$j]}};
				my @ids=keys %exon_hash2;
				
				#push @region_trs, @ids;
				foreach my $id(@ids){
					$region_trs{$id}=1;
				}
				
				$j++;
			}
			$i++;
		}
		push @region, \%region_trs;	
		push @region, 0;
		push @region, [];
		push @region, [];
		push @G_list, \@region;
	}
}

my $adj_matx;
for(my $i=0; $i <= $#G_list; $i++){
	for(my $j=0; $j <= $#G_list; $j++){
		$adj_matx->[$i][$j]=0;
	}
}

for(my $i=0; $i <= $#G_list; $i++){
	my @region=@{$G_list[$i]};
	
	#my @chunks=@{$region[0]};
	#print $chunks[0]."-".$chunks[1]."-".$region[1]."-".$region[2]."\n";
	
	if($region[1] == $region[2]){
		if($region[1] == 0){
			$adj_matx->[$i][$i+1]=1;
		}
		else{
			$adj_matx->[$i-1][$i]=1;
		}
	}
	else{
		if($region[1] == 1){
			$adj_matx->[$i][$i+1]=1;
			$adj_matx->[$i-1][$i]=1;
		}
	}
}

for(my $i=0; $i < $#G_list; $i++){
	my @region=@{$G_list[$i]};

	my @chunks=@{$region[0]};
	
	#print "Node ".$i.": ".$chunks[0]."-".$chunks[1]."-".$region[1]."-".$region[2]."\n";

	my %hash_region_trs=%{$region[3]};

	#print "\t".join(",", keys %hash_region_trs)."\n";

	if($region[1] == 0 && $region[2] == 0){
		if($adj_matx->[$i][$i+1] != 1){
			die "Region ".$chunks[0]."-".$chunks[1]." must be linked to the next one!\n";
		}
	}
	else{
		if($region[1] == 1 && $region[2] == 0){
			if($adj_matx->[$i][$i+1] != 1 || $adj_matx->[$i-1][$i] != 1){
				die "Region ".$chunks[0]."-".$chunks[1]." must be linked to the next one!\n";
			}
		}
		else{
			if($region[1] == 1 && $region[2] == 1){
				if($adj_matx->[$i-1][$i] != 1){
					die "Region ".$chunks[0]."-".$chunks[1]." must be linked to the previous one!\n";
				}
			}

			for(my $j=$i+1; $j <= $#G_list; $j++){
				my @next_region=@{$G_list[$j]};
				my %next_region_trs=%{$next_region[3]};
				my $count_inters=0;
				foreach my $tr(keys %next_region_trs){
					if(exists $hash_region_trs{$tr}){
						if($hash_region_trs{$tr} == 1){
							$hash_region_trs{$tr}=0;
							$count_inters++;
						}
					}
				}				

				if($adj_matx->[$i][$j] == 1){
					if($j > $i+1 || $count_inters == 0){
						die "Problem 1!\n";
					}
				}
				else{
					if($count_inters > 0){
						$adj_matx->[$i][$j]=1;
					}
				}
			}
		}		
	}
}

#Add the adjacency list
for(my $i=0; $i <= $#G_list; $i++){
	my @adj_node=();
	push @adj_node, $G_list[$i];
	
	for(my $j=0; $j <= $#G_list; $j++){
		if($adj_matx->[$i][$j] == 1){
			push @{${$G_list[$i]}[5]}, $G_list[$j];
			push @{${$G_list[$j]}[6]}, $G_list[$i];
		}	
	}
}


#Path contraction
my $i=0;
while($contract_path >= 1 && $i <= $#G_list){
	my $region_ref=$G_list[$i];
	#Non-contracted node
	if(${$region_ref}[4] == 0){
		#Region of type *-end	
		if($contract_path == 2 || ${$region_ref}[2] == 1){
			#print "Region ".${${$region_ref}[0]}[0]."-".${${$region_ref}[0]}[1]."\n";
			my $stop=0;
			while($stop == 0){
				my @outs=@{${$region_ref}[5]};
				#Just one adj
				if($#outs == 0){
					my $out_ref=$outs[0];
					if(${$out_ref}[4] == 1){
						die "Problem in contraction flag!\n";
					}
					#Indegree = 1
					if(scalar(@{${$out_ref}[6]}) == 1){
						#Type start-*
						if($contract_path == 2 || ${$out_ref}[1] == 0){
							#print "\tto region ".${${$out_ref}[0]}[0]."-".${${$out_ref}[0]}[1]."\n";
														
							my $link_nodes=1;
							
							#Check if the transcript sets are the same
							if($contract_path == 1){
								my %tr_first=%{${$region_ref}[3]};
								my %tr_second=%{${$out_ref}[3]};
							
								my @keys_first=keys %tr_first;
								foreach my $key(@keys_first){
									if(!exists $tr_second{$key}){
										$link_nodes=0;
									}
								}
								if($link_nodes == 1){
									my @keys_second=keys %tr_second;
									foreach my $key(@keys_second){
										if(!exists $tr_first{$key}){
											$link_nodes=0;
										}
									}
								}
							}
							
							if($link_nodes == 1){													
								${$out_ref}[4]=1; #Contract node
								${$G_list[$i]}[5]=${$out_ref}[5];
								push @{${$G_list[$i]}[0]}, @{${$out_ref}[0]};
							
								if($contract_path == 1 && ${$out_ref}[2] == 0){
									$stop=1;	
								}
								else{
									$region_ref=$out_ref;	
								}
							}
							else{
								$stop=1;	
							}
						}
						else{
							$stop=1;	
						}
					}
					else{
						$stop=1;	
					}
				}
				else{
					$stop=1;	
				}
			}
			
		}
	}
	$i++;
}

#Node labelling
my $k=0;
my $label=1;
my $label_contr=1;
while($k <= $#G_list){
	if(${$G_list[$k]}[4] == 0){
		push @{$G_list[$k]}, $label;
		$label++;
	}
	else{
		push @{$G_list[$k]}, "Contracted ".$label_contr;
		$label_contr++;
	}
	$k++;
}

print OUTGDL "graph: \{\n";
print OUTGDL "node.shape\t: circle\n";
print OUTGDL "node.color\t: blue\n";
print OUTGDL "node.height\t: 80\n";
print OUTGDL "node.width\t: 80\n";
foreach my $node_ref(@G_list){
	my @region=@{$node_ref};

	#Non-contracted node
	if($region[4] == 0){
		my @chunks=@{$region[0]};
		my $title="";
		my $i=0;
		my $dimension=0;
		my $chunk_seq="";
		while($i <= $#chunks){
			$title.=$chunks[$i]."-".$chunks[$i+1].";";
			$dimension=$dimension+$chunks[$i+1]-$chunks[$i]+1;
			$chunk_seq=$chunk_seq.substr($gen_seq, $chunks[$i]-1, $chunks[$i+1]-$chunks[$i]+1);
			$i+=2;	
		}
		print OUT "node\#".$region[7]." ".$chunk_seq."\n";
		print OUTREG "node\#".$region[7]." ".$title."\n";
		
		#print OUTGDL "\tnode: \{\n\t\ttitle: "."\"".$title."\""."\n\t\}\n";
		print OUTGDL "\tnode: \{\n\t\ttitle: "."\"".$region[7]."\"\n";
		print OUTGDL "\t\tlabel: "."\"".$region[7]." - ".$dimension."\""."\n\t\}\n";
	}
}

my $edge_index=1;
foreach my $node_ref(@G_list){
	my @region=@{$node_ref};

	#Non-contracted node
	if($region[4] == 0){
		my @chunks=@{$region[0]};
		my $title="";
		my $i=0;
		while($i <= $#chunks){
			$title.=$chunks[$i]."-".$chunks[$i+1].",";
			$i+=2;	
		}
		my @outs=@{$region[5]};
		foreach my $adj_ref(@outs){
			my @dest_region=@{$adj_ref};
			
			my @dest_chunks=@{$dest_region[0]};
			my $dest_title="";
			my $j=0;
			while($j <= $#dest_chunks){
				$dest_title.=$dest_chunks[$j]."-".$dest_chunks[$j+1].",";
				$j+=2;	
			}

			print OUT "edge\#".$edge_index." ".$region[7].";".$dest_region[7]."\n";
			print OUTREG "edge\#".$edge_index." ".$region[7].";".$dest_region[7]."\n";

			#print OUTGDL "\tedge: \{\n\t\tsource: "."\"".$title."\""."\n";
			#print OUTGDL "\t\ttarget: "."\"".$dest_title."\""."\n\t\}\n";
			print OUTGDL "\tedge: \{\n\t\tsource: "."\"".$region[7]."\""."\n";
			print OUTGDL "\t\ttarget: "."\"".$dest_region[7]."\""."\n\t\}\n";
			
			$edge_index++;
		}
	}
}
print OUTGDL "\}\n";

close OUT;
close OUTREG;
close OUTGDL;

exit;

sub reverse_complement{
	my $gen_to_be_rev=shift;
	
	my $rev_seq=reverse $gen_to_be_rev;
	my @car_list=split(//, $rev_seq);
	my $new_seq="";
	foreach my $car(@car_list){
		$new_seq=$new_seq.complement($car);
	}
		
	return $new_seq;
}

sub complement{
	my $car=shift;
	
	my $new_car;
	if($car eq "a"){
		$new_car="t";
	}elsif($car eq "A"){
		$new_car="T";
	}elsif($car eq "t"){
		$new_car="a";
	}elsif($car eq "T"){
		$new_car="A";
	}elsif($car eq "c"){
		$new_car="g";
	}elsif($car eq "C"){
		$new_car="G";
	}elsif($car eq "g"){
		$new_car="c";
	}elsif($car eq "G"){
		$new_car="C";
	}else{
		$new_car=$car;
	}
	
	return $new_car;
	
}
