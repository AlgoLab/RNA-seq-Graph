#!/usr/bin/perl -w

#Script per arricchire un set di RNA-seq con copertura bassa (lunghezza pari e costante per tutto il set)

use strict;
use POSIX;
use Getopt::Long;

our $debug = 0;
my $wdir;
my $rnafile;
my $out_file;
my $junction_file;
my $log_file;

GetOptions (
	    'working-dir=s' => \$wdir,
            'rna-seq=s' => \$rnafile,
            'output=s' => \$out_file,
	    'junctions=s' => \$junction_file,
	    'log=s' => \$log_file,
            'debug' => \$debug,
           );
my $usage="Usage: perl enrich-rna-set.pl [options]
 --working-dir=<file>		 (optional) path of the working directory
 								(default: ./)
 --rna-seq= <file>          (mandatory) name of the read file
 								
  --output=<file>			 (optional) name of the enriched read output file
 								(default: RNA-seqs-enriched.fa)
  --junctions=<file>		(optional) name of the output file of junctions
 								(default: junctions.txt)
  --log=<file>		(optional) name of the log file
 								(default: junctions.log)

 --debug                     emits debugging information
";

die $usage unless (defined $rnafile and $rnafile ne '');

$wdir="./" if (not defined $wdir or $wdir eq '' or $wdir eq ".");
print STDERR "Working directory: ", $wdir, "\n";
chop $wdir if(substr($wdir, length($wdir)-1, 1) eq "/");
            
print STDERR "RNA-seq file: ", $rnafile, "\n";

$out_file="RNA-seqs-enriched.fa" if (not defined $out_file or $out_file eq '');
print STDERR "Output read file: ", $out_file, "\n";

$junction_file="junctions.txt" if (not defined $junction_file or $junction_file eq '');
print STDERR "Output junction file: ", $junction_file, "\n";

$log_file="junctions.log" if (not defined $log_file or $log_file eq '');
print STDERR "Log file: ", $log_file, "\n";

my $dirh;
opendir ($dirh, "$wdir") or die "Cannot open $wdir!\n";
my @lists=readdir($dirh);
foreach my $dir(@lists){
		
	if($dir =~ m/out_(.+)/){

		my $gene=$1;
		print STDERR "$gene\n";
		
		my $current_dir=$wdir."/".$dir;
		#print STDERR "$current_dir\n";
		
		my $current_rnafile=$current_dir."/".$rnafile;
		my $current_out_file=$current_dir."/".$out_file;
		my $current_junction_file=$current_dir."/".$junction_file;
		my $current_log_file=$current_dir."/".$log_file;

		open(IN, "<$current_rnafile") or die "Cannot open $current_rnafile!\n";

		my $read="";
		my $header="";
		my %headers=();

		#key=index; value=[read_sequence, [@left_junction_sequences], is_left_junction_init, is_included_in_left_junction, [@right_junction_sequences], is_right_junction_init, is_included_in_right_junction, left_starting_read_index, right_starting_read_index, left_shared, right_shared]
		#Le left junctions sono viste da sinistra (cioé il taglio è nella metà destra, e riguardano una giunzione tra un blocco e i suoi uscenti), mentre le right sono viste da destra (cioé il taglio è nella metà sinistra, e riguardano una giunzione tra un blocco e i suoi entranti). Per le giunzioni da destra, le sequenze in questa struttura sono reverse
		#left_starting_read_index è la lista (CAMBIARE IN HASH) degli indici dei left-most reads delle giunzioni sinistre in cui il read è incluso, analogo per right_starting_read_index.
		#left_shared, right_shared sono le lunghezze della regione comune (prefisso e suffisso)
		my %reads=();

		my $index=1;
		my %presence=();

		while(<IN>){
			chomp;
			if($_ =~ /^>/){
				#print STDERR $_, "\n";
				if($header ne ""){
					if(!exists $presence{$read}){
						$reads{$index}=[$read, [], 0, 0, [], 0, 0, [], [], -1, -1];
						$presence{$read}=1;
						$headers{$index}=$header;
						$index++;
					}				
				}
				$header=$_;
				$read="";
				
			}else{
				$read=$read.$_;
			}
		}
		if(!exists $presence{$read}){
			$reads{$index}=[$read, [], 0, 0, [], 0, 0, [], [], -1, -1];
			$presence{$read}=1;
			$headers{$index}=$header;
		}

		close IN;

		open(OUT, ">$current_out_file") or die "Cannot create $current_out_file!\n";
		open(OUT_JUN, ">$current_junction_file") or die "Cannot create $current_junction_file!\n";
		open(LOG, ">$current_log_file") or die "Cannot create $current_log_file!\n";

		my %left_hash=();
		my %right_hash=();
		foreach my $key(keys %reads){
			my $read=${$reads{$key}}[0];
			my $left_half=substr($read, 0, length($read)/2);
			my $right_half=substr($read, length($read)/2, length($read));
			if(exists $left_hash{$left_half}){
				push @{$left_hash{$left_half}}, $key;
			}else{
				$left_hash{$left_half}=[$key];
			}
			if(exists $right_hash{$right_half}){
				push @{$right_hash{$right_half}}, $key;
			}else{
				$right_hash{$right_half}=[$key];
			}			
		}

		print LOG "\n**Costruzione delle giunzioni viste da sinistra\n\n";
		#Per ogni input read cerco l'eventuale taglio nella metà destra (usando l'hash sinistra)
		my @indices=keys %reads;
		foreach my $in(@indices){
			my $sequence=${$reads{$in}}[0];
			my $l=length($sequence);
			print LOG "...read $in:\n", ${$reads{$in}}[0], "\n";
			#Se non è stato incluso in una giunzione
			if(${$reads{$in}}[3] == 0){
				print LOG "\t...inizio di una possibile giunzione\n";		
				my @extensions=($sequence);
				my $i=0;
				my $stop=0;
				my @all_ids;
				while($i <= $l/2-1 && $stop == 0){
					my $hash_seq=substr($sequence, $i, $l/2);
					if(exists $left_hash{$hash_seq}){
						my @ids=@{$left_hash{$hash_seq}};
						push @all_ids, @ids;
				
						my $cfr_str=join("-", @all_ids);							
						
						foreach my $ext_in(@ids){
							if($ext_in != $in){
								print LOG "*"x$i, ${$reads{$ext_in}}[0], "\tread $ext_in\n";
								#Se ha già dato origine a qualche giunzione
								if(${$reads{$ext_in}}[2] == 1){
									print LOG "\t\tesiste già una giunzione:\n";
									my @junctions=@{${$reads{$ext_in}}[1]};
									my $j=0;
									while($j <= $#junctions){
										$junctions[$j]=("*" x $i).$junctions[$j];

										print LOG "\t\t", $junctions[$j], "\n";

										$j++;
									}
									print LOG "\n";
									push @extensions, @junctions;
								}else{
									print LOG "\t\tda includere:\n";
									push @extensions, ("*" x $i).${$reads{$ext_in}}[0];
								}
							}
						}
					}
					$i++;
				}

				if(scalar(@extensions) > 1){
					my $junctions_ref=search_cut(\@extensions, length($extensions[0])/2, \*LOG);
					my @junctions=@{${$junctions_ref}[0]};
					my $shared_length=${$junctions_ref}[1];
								
					if(scalar(@junctions) > 0){
						
						${$reads{$in}}[1]=\@junctions;
						${$reads{$in}}[2]=1;
						push @{${$reads{$in}}[7]}, $in;
						${$reads{$in}}[9]=$shared_length;
											
						print LOG "...read $in fornisce la giunzione:\n";
						print LOG join("\n", @junctions), "\n";

						my $j=0;
						while($j <= $#all_ids){
							if($all_ids[$j] != $in){
								print LOG "\tread ", $all_ids[$j] ,"\n";
								if(${$reads{$all_ids[$j]}}[2] == 1){
									${$reads{$all_ids[$j]}}[1]=[];
									${$reads{$all_ids[$j]}}[2]=0;
									${$reads{$all_ids[$j]}}[9]=-1;
									print LOG "\t\t(non è più origine di giunzione)\n";
								}
								push @{${$reads{$all_ids[$j]}}[7]}, $in;
								${$reads{$all_ids[$j]}}[3]=1;
								print LOG "\t\tincluso nella giunzione\n";
							}
							$j++;
						}						

						print LOG "***************************************\n";
					}else{
						print LOG "\t...NON fornisce una giunzione\n";
					}
				}else{
					print LOG "\t...NON ha read che lo estendono\n";
					print LOG "***************************************\n";
				}
			}
			else{
				print LOG "\t...è già incluso in una giunzione\n";
				print LOG "***************************************\n";
			}
		}

		print LOG "Giunzioni viste da sinistra al momento della costruzione\n\n";
		foreach my $in(keys %reads){
			if(${$reads{$in}}[2] == 1){
				print LOG "\t..read $in: ", ${$reads{$in}}[0], " shared: ", ${$reads{$in}}[9], "\n";		
				print LOG join("\n", @{${$reads{$in}}[1]}), "\n\n";
			}
		}
		print LOG "***************************************\n";

		print LOG "Costruzione giunzioni viste da destra\n";
		#Per ogni input read cerco l'eventuale taglio nella metà sinistra (usando l'hash destra)
		@indices=keys %reads;
		foreach my $in(@indices){
			my $sequence=${$reads{$in}}[0];
			my $l=length($sequence);
			print LOG "Read numero $in:\n", ${$reads{$in}}[0], "\n";
			#Se non è stato incluso in una giunzione
			if(${$reads{$in}}[6] == 0){
				print LOG "\t...inizio di una possibile giunzione (si considera il reverse)\n";
				my $reverse_seq=reverse($sequence);
				print LOG $reverse_seq, "\n";
				my @extensions=($reverse_seq);
				my $i=0;
				my $stop=0;
				my @all_ids;
				while($i <= $l/2-1 && $stop == 0){
					my $hash_seq_rev=substr($reverse_seq, $i, $l/2);
					my $hash_seq=reverse($hash_seq_rev);

					if(exists $right_hash{$hash_seq}){
						my @ids=@{$right_hash{$hash_seq}};
						push @all_ids, @ids;

						my $cfr_str=join("-", @all_ids);					
						
						foreach my $ext_in(@ids){
							if($ext_in != $in){
								my $rev_read=reverse(${$reads{$ext_in}}[0]);
								print LOG "*"x$i, $rev_read, "\tread $ext_in\n";
								#Se ha già dato origine a qualche giunzione
								if(${$reads{$ext_in}}[5] == 1){
									print LOG "\t\tesiste già una giunzione:\n";
									my @junctions=@{${$reads{$ext_in}}[4]};
									my $j=0;
									while($j <= $#junctions){
										$junctions[$j]=("*" x $i).$junctions[$j];

										print LOG $junctions[$j], "\n";

										$j++;
									}
									print LOG "\n";
									push @extensions, @junctions;
								}else{
									print LOG "\t\tda includere:\n";
									my $rev_seq=reverse(${$reads{$ext_in}}[0]);
									push @extensions, ("*" x $i).$rev_seq;
								}
							}
						}
					}
					$i++;
				}

				if(scalar(@extensions) > 1){
					my $junctions_ref=search_cut(\@extensions, length($extensions[0])/2, \*LOG);
					my @junctions=@{${$junctions_ref}[0]};
					my $shared_length=${$junctions_ref}[1];
			
					if(scalar(@junctions) > 0){
						
						${$reads{$in}}[4]=\@junctions;
						${$reads{$in}}[5]=1;
						push @{${$reads{$in}}[8]}, $in;
						${$reads{$in}}[10]=$shared_length;

						print LOG "Read numero $in corrisponde alla giunzione dx (sul reverse):\n";
						print LOG join("\n", @junctions), "\n";

						my $j=0;
						while($j <= $#all_ids){
							if($all_ids[$j] != $in){
								print LOG "\tRead numero ", $all_ids[$j] ,"\n";
								if(${$reads{$all_ids[$j]}}[5] == 1){
									${$reads{$all_ids[$j]}}[4]=[];
									${$reads{$all_ids[$j]}}[5]=0;
									${$reads{$all_ids[$j]}}[10]=-1;
									print LOG "\t\tnon è più origine della giunzione dx\n";
								}
								push @{${$reads{$all_ids[$j]}}[8]}, $in;
								${$reads{$all_ids[$j]}}[6]=1;
								print LOG "\t\tincluso nella giunzione dx\n";
							}
							$j++;
						}
						
						print LOG "***************************************\n";
					}else{
						print LOG "no taglio\n";
					}
				}else{
					print LOG "no estensione\n";
					print LOG "***************************************\n";
				}
			}
			else{
				print LOG "già incluso in giunzione dx!\n";
				print LOG "***************************************\n";
			}
		}

		print LOG "Giunzioni viste da destra al momento della costruzione\n\n";
		foreach my $in(keys %reads){
			if(${$reads{$in}}[5] == 1){
				print LOG "Giunzione read numero $in: ", ${$reads{$in}}[0], " shared ", ${$reads{$in}}[10], "\n";
				my $max=0;
				foreach my $g(@{${$reads{$in}}[4]}){
					if($max < length($g)){
						$max=length($g);
					}					
				}
				foreach my $g(@{${$reads{$in}}[4]}){
					my $rev=reverse($g);
					print LOG " "x($max-length($g)), $rev, "\n";
				}
				print LOG "\n";
			}
		}
		print LOG "***************************************\n";
		
		#Per giunzioni da sinistra:
			#Estensione a sinistra del taglio (dalla parte del prefisso condiviso) fino ad almeno 64bp prima 				#del taglio (se possibile). Aggiornamento di left_shared (left_starting_read_index rimane il 				#leftmost read spliced). L'estensione avviene solo con read che non sono già stati inclusi in 				#alcuna giunzione. Scopo: facilitare il completamento delle catene ai bordi includendo anche 				#unspliced mancanti.			
		foreach my $in(keys %reads){
			if(${$reads{$in}}[2] == 1){
				#print LOG "Giunzione sx read numero $in: ", ${$reads{$in}}[0], " shared: ", ${$reads{$in}}[9], "bp\n";		
				my $shared_region_length=${$reads{$in}}[9];
				my $l=length(${$reads{$in}}[0]);

				if($shared_region_length >= $l/2){
					my @junctions=@{${$reads{$in}}[1]};
					my $some_junction=$junctions[0];
					my $shared_region=substr($some_junction, 0, $shared_region_length);
					#print LOG "Shared region ", $shared_region, "\n";	
					my $i=$shared_region_length-1;
					my $stop=0;
					while($i >= $l/2-1 && $stop == 0){
						my $hash_seq=substr($shared_region, $i-$l/2+1, $l/2);
						
						#print LOG $hash_seq, "\n";
						
						if(exists $right_hash{$hash_seq}){
							my @ids=@{$right_hash{$hash_seq}};
							#Se il read è unico (RIVEDERE!!!)
							if(scalar(@ids) == 1){
								#Se ancora non è incluso in giunzioni
								if(${$reads{$ids[0]}}[2] == 0 && ${$reads{$ids[0]}}[3] == 0){									
									#print LOG "Read: ", ${$reads{$ids[0]}}[0], "\n";
									#Se esiste overlapping
									my $min=($i-$l/2+1 <= $l/2)?($i-$l/2+1):($l/2);
									my $check1=substr($shared_region, $i-$l/2-$min+1, $min);
									my $check2=substr(${$reads{$ids[0]}}[0], $l/2-1-$min+1, $min);
									if($check1 eq $check2){
										my $add_length=$l/2-$min;
										if($add_length > 0){
											$stop=1;
											${$reads{$in}}[9]=${$reads{$in}}[9]+$add_length;
											${$reads{$ids[0]}}[3]=1;
											push @{${$reads{$ids[0]}}[7]}, $in;
											my $add_string=substr(${$reads{$ids[0]}}[0], 0, $add_length);
											my $k=0;
											while($k <= $#junctions){
												${${$reads{$in}}[1]}[$k]=$add_string.${${$reads{$in}}[1]}[$k];
												$k++;
											}
										}
									}								
								}
							}
						}
				
						$i--;
					}
				}

				
			}
		}

		print LOG "Giunzioni viste da sinistra dopo estensione verso sinistra\n\n";
		foreach my $in(keys %reads){
			if(${$reads{$in}}[2] == 1){
				print LOG "Giunzione read numero $in: ", ${$reads{$in}}[0], " shared: ", ${$reads{$in}}[9], "\n";		
				print LOG join("\n", @{${$reads{$in}}[1]}), "\n\n";
			}
		}
		print LOG "***************************************\n";

		#ESTENSIONE VERSO DESTRA
		foreach my $in(keys %reads){
			if(${$reads{$in}}[2] == 1){
				#print LOG "Giunzione sx read numero $in: ", ${$reads{$in}}[0], " shared: ", ${$reads{$in}}[9], "bp\n";		
				my $l=length(${$reads{$in}}[0]);
				my @junctions=@{${$reads{$in}}[1]};
				my $i=0;
				while($i <= $#junctions){
					my $junction=$junctions[$i];
					my $after_cut=length($junction)-${$reads{$in}}[9];
					#print LOG "\tgiunzione $i: ", $junction, " after cut: $after_cut\n";
					my $enlarge=1;
					while($enlarge == 1){
						my $k=length($junction)-$l+1;
						my $stop=0;
						while($after_cut-length($junction)+$k <= 0 && $stop == 0){
							my $hash_seq=substr($junction, $k, $l/2);
							#print LOG "\t\t$hash_seq\n";
							if(exists $left_hash{$hash_seq}){
								my @ids=@{$left_hash{$hash_seq}};
								#Se il read è unico (RIVEDERE!!!)
								if(scalar(@ids) == 1){
									#print LOG "\t\tRead: ", ${$reads{$ids[0]}}[0], "\n";
									#Se ancora non è incluso in giunzioni
									if(${$reads{$ids[0]}}[2] == 0 && ${$reads{$ids[0]}}[3] == 0){									
										my $check1=substr($junction, $k, length($junction)-$k);
										my $check2=substr(${$reads{$ids[0]}}[0], 0, length($junction)-$k);
										#print LOG $check1, " ", $check2, "\n";
										if($check1 eq $check2){
											$stop=1;
											${$reads{$ids[0]}}[3]=1;
											push @{${$reads{$ids[0]}}[7]}, $in;
											my $ext=$l-length($junction)+$k;
											$after_cut+=$ext;
											$junction=$junction.substr(${$reads{$ids[0]}}[0], length(${$reads{$ids[0]}}[0])-$ext, $ext);
											#print LOG "New after cut $after_cut\n";
											#print LOG "New junction $junction\n";
											
										}
									}
								}
							}
							$k++;
						}
						if($stop == 0){
							$enlarge=0;
						}						
					}

					#Aggiorno la giunzione $i-esima nella struttura con $junction
					${${$reads{$in}}[1]}[$i]=$junction;
					$i++;
				}
			}
			
		}

		print LOG "Giunzioni viste da sinistra dopo estensione verso destra\n\n";
		foreach my $in(keys %reads){
			if(${$reads{$in}}[2] == 1){
				print LOG "Giunzione read numero $in: ", ${$reads{$in}}[0], " shared: ", ${$reads{$in}}[9], "\n";		
				print LOG join("\n", @{${$reads{$in}}[1]}), "\n\n";
			}
		}
		print LOG "***************************************\n";

		#Merging tra giunzioni viste da sinistra
		my @id_keys=keys %reads;
		my $z=0;
		my @merge_flags=();
		while($z <= $#id_keys){
			if(${$reads{$id_keys[$z]}}[2] == 1){
				push @merge_flags, $id_keys[$z];
			}
			$z++;
		}
		$z=0;
		while($z <= $#merge_flags){
			#Se la giunzione non è stata ancora fusa (cioé contiene ancora l'indice del read iniziale)
			if($merge_flags[$z] != -1){
				my $l=length(${$reads{$merge_flags[$z]}}[0]);
				my @junctions1=@{${$reads{$merge_flags[$z]}}[1]};
				#print STDERR "RAFFA: ", $merge_flags[$z], "\n", join("\n", @junctions1), "\n";
				my $left_shared1=${$reads{$merge_flags[$z]}}[9];
				my $shared_region1=substr($junctions1[0], 0, $left_shared1);
				my $y=$z+1;
				while($y <= $#merge_flags){
					if($merge_flags[$y] != -1){
						my @junctions2=@{${$reads{$merge_flags[$y]}}[1]};
						my $left_shared2=${$reads{$merge_flags[$y]}}[9];
						my $shared_region2=substr($junctions2[0], 0, $left_shared2);
						#print STDERR "RAFFA\n", join("\n", @junctions2), "\n";

						my $merge=0;
						my $k=0;
						my $left_ext=0;
						while($k <= $#junctions1 && $merge == 0){
							#print STDERR "\t", $shared_region2, " in ", $junctions1[$k], "\n";
							if($junctions1[$k] =~ m/$shared_region2/){
								$left_ext=length($`);
								
								#$merge=1;
								if($left_ext < $left_shared1){
									$merge=1;
								}
							}
							$k++;
						}

						$k=0;
						while($k <= $#junctions2 && $merge == 0){
							#print STDERR "\t", $shared_region1, " in ", $junctions2[$k], "\n";
							if($junctions2[$k] =~ m/$shared_region1/){
								
								$left_ext=length($`);
								#Caso gene EEF1A1
								#$merge=2;
								if($left_ext < $left_shared2){
									$merge=2;
								}
							}
							$k++;
						}

						if($merge != 0){
							my @extensions;
							if($merge == 1){
								@extensions=(@junctions1);
								foreach my $j2(@junctions2){
									push @extensions, ("*" x $left_ext).$j2;
								}
							}else{
								@extensions=(@junctions2);
								foreach my $j1(@junctions1){
									push @extensions, ("*" x $left_ext).$j1;
								}

							}
							#print STDERR "NEW\n", join("\n", @extensions), "\n";
							my $junctions_ref=search_cut(\@extensions, $l/2, \*LOG);
							my @merge_junctions=@{${$junctions_ref}[0]};
							my $merge_shared_length=${$junctions_ref}[1];
							if(scalar(@merge_junctions) > 0){
								#print STDERR "CUT ", $merge_flags[$y], "!\n";
								${$reads{$merge_flags[$y]}}[1]=[];
								${$reads{$merge_flags[$y]}}[2]=0;
								${$reads{$merge_flags[$y]}}[9]=-1;
								@junctions1=@merge_junctions;
								$left_shared1=$merge_shared_length;
								#print STDERR "\n", join("\n", @junctions1), "\n";
								$merge_flags[$y]=-1;
							}
						}
					}
					$y++;
				}
				#print STDERR "\nHERE\n", join("\n", @junctions1), "\n";
				${$reads{$merge_flags[$z]}}[1]=\@junctions1;
				${$reads{$merge_flags[$z]}}[9]=$left_shared1;
				$merge_flags[$z]=-1;
			}
			
			$z++;
		}

		print LOG "Giunzioni viste da sinistra dopo merging\n\n";
		foreach my $in(keys %reads){
			if(${$reads{$in}}[2] == 1){
				print LOG "Giunzione read numero $in: ", ${$reads{$in}}[0], " shared: ", ${$reads{$in}}[9], "\n";		
				print LOG join("\n", @{${$reads{$in}}[1]}), "\n\n";
			}
		}
		print LOG "***************************************\n";	

		#Per giunzioni da destra:
			#Estensione a destra del taglio (dalla parte del suffisso condiviso) fino ad almeno 64bp dopo 				#il taglio (se possibile). Aggiornamento di right_shared (right_starting_read_index rimane il 				#rightmost read spliced). L'estensione avviene solo con read che non sono già stati inclusi in 				#alcuna giunzione. Scopo: facilitare il completamento delle catene ai bordi includendo anche 				#unspliced mancanti.
		foreach my $in(keys %reads){
			if(${$reads{$in}}[5] == 1){
				#print LOG "Giunzione dx read numero $in: ", ${$reads{$in}}[0], " shared: ", ${$reads{$in}}[10], "bp\n";	
				my $shared_region_length=${$reads{$in}}[10];
				my $l=length(${$reads{$in}}[0]);

				if($shared_region_length >= $l/2){
					my @junctions=@{${$reads{$in}}[4]};
					my $some_junction=$junctions[0];
					my $shared_region=substr($some_junction, 0, $shared_region_length);
					my $rev=reverse($shared_region);
					#print LOG "Shared region ", $rev, "\n";
					my $i=0;
					my $stop=0;
					while($i <= $shared_region_length-$l/2 && $stop == 0){
						my $hash_seq=substr($rev, $i, $l/2);

						#print LOG $hash_seq, "\n";
						if(exists $left_hash{$hash_seq}){
							my @ids=@{$left_hash{$hash_seq}};			
							#Se il read è unico (RIVEDERE!!!)
							if(scalar(@ids) == 1){
								#Se ancora non è incluso in giunzioni
								if(${$reads{$ids[0]}}[5] == 0 && ${$reads{$ids[0]}}[6] == 0){									
									#print LOG "Read: ", ${$reads{$ids[0]}}[0], "\n";
									#Se esiste overlapping
									my $min=($shared_region_length-$i-$l/2 <= $l/2)?($shared_region_length-$i-$l/2):($l/2);
									my $check1=substr($rev, $i+$l/2, $min);
									my $check2=substr(${$reads{$ids[0]}}[0], $l/2, $min);									
									
									if($check1 eq $check2){
										my $add_length=$l/2-$min;
										if($add_length > 0){
											$stop=1;
											${$reads{$in}}[10]=${$reads{$in}}[10]+$add_length;
											${$reads{$ids[0]}}[6]=1;
											push @{${$reads{$ids[0]}}[8]}, $in;
											my $add_string=substr(${$reads{$ids[0]}}[0], $l-$add_length, $add_length);
											my $add_rev=reverse($add_string);
											my $k=0;
											while($k <= $#junctions){
												${${$reads{$in}}[4]}[$k]=$add_rev.${${$reads{$in}}[4]}[$k];
												$k++;
											}
										}
																		
									}								
								}
							}
						}
						
						$i++;
					}					
				}				
			}
		}
		
		print LOG "Giunzioni viste da destra dopo estensione verso destra\n\n";
		foreach my $in(keys %reads){
			if(${$reads{$in}}[5] == 1){
				print LOG "Giunzione read numero $in: ", ${$reads{$in}}[0], " shared ", ${$reads{$in}}[10], "\n";
				my $max=0;
				foreach my $g(@{${$reads{$in}}[4]}){
					if($max < length($g)){
						$max=length($g);
					}					
				}
				foreach my $g(@{${$reads{$in}}[4]}){
					my $rev=reverse($g);
					print LOG " "x($max-length($g)), $rev, "\n";
				}
				print LOG "\n";
			}
		}
		print LOG "***************************************\n";

		#ESTENSIONE VERSO SINISTRA
		foreach my $in(keys %reads){
			if(${$reads{$in}}[5] == 1){
				#print LOG "Giunzione dx read numero $in: ", ${$reads{$in}}[0], " shared: ", ${$reads{$in}}[10], "bp\n";		
				my $l=length(${$reads{$in}}[0]);
				my @junctions=@{${$reads{$in}}[4]};
				my $i=0;
				while($i <= $#junctions){
					my $junction=$junctions[$i];
					my $after_cut=length($junction)-${$reads{$in}}[10];
					#print LOG "\tgiunzione (su reverse) $i: ", $junction, " after cut: $after_cut\n";
					my $enlarge=1;
					while($enlarge == 1){
						my $k=length($junction)-$l+1;
						my $stop=0;
						while($after_cut-length($junction)+$k <= 0 && $stop == 0){
							my $hash_seq=substr($junction, $k, $l/2);
							my $rev=reverse($hash_seq);
							#print LOG "\t\t$rev\n";						
							if(exists $right_hash{$rev}){
								my @ids=@{$right_hash{$rev}};
								#Se il read è unico (RIVEDERE!!!)
								if(scalar(@ids) == 1){
									#print LOG "\t\tRead: ", ${$reads{$ids[0]}}[0], "\n";
									#Se ancora non è incluso in giunzioni
									if(${$reads{$ids[0]}}[5] == 0 && ${$reads{$ids[0]}}[6] == 0){									
										my $check1=substr($junction, $k, length($junction)-$k);
										my $rev_read=reverse(${$reads{$ids[0]}}[0]);
										my $check2=substr($rev_read, 0, length($junction)-$k);
										#print LOG $check1, " ", $check2, "\n";
										if($check1 eq $check2){
											$stop=1;
											${$reads{$ids[0]}}[6]=1;
											push @{${$reads{$ids[0]}}[8]}, $in;
											my $ext=$l-length($junction)+$k;
											$after_cut+=$ext;
											$junction=$junction.substr($rev_read, length(${$reads{$ids[0]}}[0])-$ext, $ext);
											#print LOG "New after cut $after_cut\n";
											#print LOG "New junction $junction\n";
											
										}
									}
								}
							}
							$k++;
						}
						if($stop == 0){
							$enlarge=0;
						}						
					}

					#Aggiorno la giunzione $i-esima nella struttura con $junction
					${${$reads{$in}}[4]}[$i]=$junction;
					$i++;
				}
			}
			
		}

		print LOG "Giunzioni viste da destra dopo estensione verso sinistra\n\n";
		foreach my $in(keys %reads){
			if(${$reads{$in}}[5] == 1){
				print LOG "Giunzione read numero $in: ", ${$reads{$in}}[0], " shared ", ${$reads{$in}}[10], "\n";
				my $max=0;
				foreach my $g(@{${$reads{$in}}[4]}){
					if($max < length($g)){
						$max=length($g);
					}					
				}
				foreach my $g(@{${$reads{$in}}[4]}){
					my $rev=reverse($g);
					print LOG " "x($max-length($g)), $rev, "\n";
				}
				print LOG "\n";
			}
		}
		print LOG "***************************************\n";

		#Merging tra giunzioni viste da destra
		@id_keys=keys %reads;
		$z=0;
		@merge_flags=();
		while($z <= $#id_keys){
			if(${$reads{$id_keys[$z]}}[5] == 1){
				push @merge_flags, $id_keys[$z];
			}
			$z++;
		}
		$z=0;
		while($z <= $#merge_flags){
			#Se la giunzione non è stata ancora fusa (cioé contiene ancora l'indice del read finale)
			if($merge_flags[$z] != -1){
				my $l=length(${$reads{$merge_flags[$z]}}[0]);
				my @junctions1=@{${$reads{$merge_flags[$z]}}[4]};
				#print STDERR "RAFFA: ", $merge_flags[$z], "\n", join("\n", @junctions1), "\n";
				my $right_shared1=${$reads{$merge_flags[$z]}}[10];
				my $shared_region1=substr($junctions1[0], 0, $right_shared1);
				my $y=$z+1;
				while($y <= $#merge_flags){
					if($merge_flags[$y] != -1){
						my @junctions2=@{${$reads{$merge_flags[$y]}}[4]};
						my $right_shared2=${$reads{$merge_flags[$y]}}[10];
						my $shared_region2=substr($junctions2[0], 0, $right_shared2);
						#print STDERR "RAFFA\n", join("\n", @junctions2), "\n";

						my $merge=0;
						my $k=0;
						my $right_ext=0;
						while($k <= $#junctions1 && $merge == 0){
							#print STDERR "\t", $shared_region2, " in ", $junctions1[$k], "\n";
							if($junctions1[$k] =~ m/$shared_region2/){
							
								$right_ext=length($`);
								#$merge=1;
								if($right_ext < $right_shared1){
									$merge=1;
								}
							}
							$k++;
						}

						$k=0;
						while($k <= $#junctions2 && $merge == 0){
							#print STDERR "\t", $shared_region1, " in ", $junctions2[$k], "\n";
							if($junctions2[$k] =~ m/$shared_region1/){
								
								$right_ext=length($`);
								#$merge=1;
								if($right_ext < $right_shared2){
									$merge=2;
								}
							}
							$k++;
						}

						if($merge != 0){
							my @extensions;
							if($merge == 1){
								@extensions=(@junctions1);
								foreach my $j2(@junctions2){
									push @extensions, ("*" x $right_ext).$j2;
								}
							}else{
								@extensions=(@junctions2);
								foreach my $j1(@junctions1){
									push @extensions, ("*" x $right_ext).$j1;
								}

							}
							#print STDERR "NEW\n", join("\n", @extensions), "\n";
							my $junctions_ref=search_cut(\@extensions, $l/2, \*LOG);
							my @merge_junctions=@{${$junctions_ref}[0]};
							my $merge_shared_length=${$junctions_ref}[1];
							if(scalar(@merge_junctions) > 0){
								#print STDERR "CUT ", $merge_flags[$y], "!\n";
								${$reads{$merge_flags[$y]}}[4]=[];
								${$reads{$merge_flags[$y]}}[5]=0;
								${$reads{$merge_flags[$y]}}[10]=-1;
								@junctions1=@merge_junctions;
								$right_shared1=$merge_shared_length;
								#print STDERR "\n", join("\n", @junctions1), "\n";
								$merge_flags[$y]=-1;
							}
						}
					}
					$y++;
				}
				#print STDERR "\nHERE\n", join("\n", @junctions1), "\n";
				${$reads{$merge_flags[$z]}}[4]=\@junctions1;
				${$reads{$merge_flags[$z]}}[10]=$right_shared1;
				$merge_flags[$z]=-1;
			}
			
			$z++;
		}

		print LOG "Giunzioni viste da destra dopo merging\n\n";
		foreach my $in(keys %reads){
			if(${$reads{$in}}[5] == 1){
				print LOG "Giunzione read numero $in: ", ${$reads{$in}}[0], " shared ", ${$reads{$in}}[10], "\n";
				my $max=0;
				foreach my $g(@{${$reads{$in}}[4]}){
					if($max < length($g)){
						$max=length($g);
					}					
				}
				foreach my $g(@{${$reads{$in}}[4]}){
					my $rev=reverse($g);
					print LOG " "x($max-length($g)), $rev, "\n";
				}
				print LOG "\n";
			}
		}
		print LOG "***************************************\n";

		print LOG "Controllo delle ripetizioni che generano false giunzioni";
		my %left_check_false=();
		my %right_check_false=();

		my $min_repeat_length=length(${$reads{1}}[0])/2;	#32bp
		#my $min_repeat_length=length(${$reads{1}}[0])/4;	#16bp

		#Controllo delle false giunzioni dovute a ripetizioni di almeno 32bp
		foreach my $in(keys %reads){
			if(${$reads{$in}}[2] == 1){
				my @l_junctions=@{${$reads{$in}}[1]};
				my $left_shared=${$reads{$in}}[9];
				die "Error3!\n" if(($left_shared-$min_repeat_length) < 0);
				my $cfr_seq=substr($l_junctions[0], $left_shared-$min_repeat_length, $min_repeat_length);
				if(exists $left_check_false{$cfr_seq}){
					push @{$left_check_false{$cfr_seq}}, $in;
				}
				else{
					$left_check_false{$cfr_seq}=[$in];
				}
			}

			if(${$reads{$in}}[5] == 1){
				my @r_junctions=@{${$reads{$in}}[4]};
				my $right_shared=${$reads{$in}}[10];
				die "Error4!\n" if(($right_shared-$min_repeat_length) < 0);
				my $cfr_seq=substr($r_junctions[0], $right_shared-$min_repeat_length, $min_repeat_length);
				if(exists $right_check_false{$cfr_seq}){
					push @{$right_check_false{$cfr_seq}}, $in;
				}
				else{
					$right_check_false{$cfr_seq}=[$in];
				}
			}
		}

		foreach my $repeat(keys %left_check_false){
			my @ids=@{$left_check_false{$repeat}};
			if(scalar(@ids) > 1){
				print LOG "Candidate repeat (from left) of ", length($repeat),"bp: $repeat\n";
				foreach my $id(@ids){
					die "Error4!\n" if(${$reads{$id}}[2] != 1);
					my @junctions=@{${$reads{$id}}[1]};
					print LOG "Junction:\n", join("\n", @junctions), "\n";
				}		
				print LOG "********************************************\n";

				my $i=0;
				my @check_flags=();
				foreach my $id(@ids){
					push @check_flags, 0;
				}
				my $cluster_check=1;
				my %cluster_check_hash=();
				while($i <= $#ids){
					my @junctions1=@{${$reads{$ids[$i]}}[1]};
					if($check_flags[$i] == 0){
						
					 	my $j=$i+1;
					 	my $check=0;
						my $left_shared1=${$reads{$ids[$i]}}[9];	
						my %suffices1=();
						foreach my $j(@junctions1){
							$suffices1{substr($j, $left_shared1, length($j)-$left_shared1)}=1;
						}
						die "Error5!\n" if(scalar(@junctions1) != scalar(keys %suffices1));
						print LOG "SUFF1 ", join("\n", keys %suffices1), "\n";

						my @check_list=($ids[$i]);
					 	while($j <= $#ids){
							my @junctions2=@{${$reads{$ids[$j]}}[1]};
							if(scalar(@junctions1) == scalar(@junctions2)){
								my $left_shared2=${$reads{$ids[$j]}}[9];	
								my %suffices2=();
								foreach my $j(@junctions2){
									$suffices2{substr($j, $left_shared2, length($j)-$left_shared2)}=1;
								}
								die "Error6!\n" if(scalar(@junctions2) != scalar(keys %suffices2));
								print LOG "SUFF2 ", join("\n", keys %suffices2), "\n";

								my $count_equals=0;

								my @keys1=keys %suffices1;
								my @keys2=keys %suffices2;
								my $p=0;
								while($p <= $#keys1){
									my $q=0;
									my $ok=0;
									while($q <= $#keys2 && $ok == 0){
										if($keys1[$p] =~ m/^$keys2[$q]/ || $keys2[$q] =~ m/^$keys1[$p]/){
											$ok=1;
											$count_equals++;
										}
										$q++;			
									}
									$p++;						
								}
								print LOG "Count $count_equals\n";

								if($count_equals == scalar(@junctions2)){
									$check=1;
									push @check_list, $ids[$j];
									$check_flags[$j]=$cluster_check;
								}							
							}							

							$j++;
					 	}
						
						if($check == 1){
							$check_flags[$i]=$cluster_check;
							$cluster_check_hash{$cluster_check}=\@check_list;
							$cluster_check++;
						}
					}
					$i++;
				}
				print LOG "False junctions:\n";
				foreach my $cluster(keys %cluster_check_hash){
					my @check_list=@{$cluster_check_hash{$cluster}};
					print LOG "Cluster to be checked:\n";
					my @shared_regions=();
					foreach my $id(@check_list){
						die "Error6!\n" if(${$reads{$id}}[2] != 1);
						my @junctions=@{${$reads{$id}}[1]};
						print LOG "Junction:\n", join("\n", @junctions), "\n";
						my $left_shared=${$reads{$id}}[9];
						push @shared_regions, substr($junctions[0], 0, $left_shared);
					}
					print LOG "Shared:\n", join("\n", @shared_regions), "\n";
					
					my $min_limit=length($shared_regions[0]);
					my $k=1;
					while($k <= $#shared_regions){
						if($min_limit > length($shared_regions[$k])){
							$min_limit=length($shared_regions[$k]);
						}
						$k++;
					}
					$k=$min_repeat_length;
					my $stop=0;
					while($k < $min_limit && $stop == 0){
						my $p=1;
						while($p <= $#shared_regions && $stop == 0){
							if(substr($shared_regions[0], length($shared_regions[0])-$k-1, 1) ne substr($shared_regions[$p], length($shared_regions[$p])-$k-1, 1)){
								$stop=1;
							}
							$p++;
						}					

						$k++;
					}
					my $repeat_length=$k-1;
					print LOG "Lunghezza max della ripetizione: $repeat_length\n";

					$k=0;
					my @valid_junctions_ref=();
					while($k <= $#check_list){
						my @junctions=@{${$reads{$check_list[$k]}}[1]};
						my $left_shared=${$reads{$check_list[$k]}}[9];
						my @valid_junctions=();

						foreach my $junction(@junctions){
							print LOG "Junction: $junction\n";

							my $valid=0;

							my $l=length(${$reads{$check_list[$k]}}[0]);
							my $min_range=($l-2 < $left_shared)?($l-2):($left_shared);
							my $max_range=($repeat_length > ($l-length($junction)+$left_shared))?($repeat_length):($l-length($junction)+$left_shared);
							#print LOG "MIN: $min_range MAX: $max_range\n";
							
							if($min_range >= $max_range && $valid == 0){
								my $p=$min_range;
								
								while($p >= $max_range && $valid == 0){
									my $validation_str=substr($junction, $left_shared-$p-1, $l);
									#print LOG "Candidate validation read: $validation_str\n";
									if(exists $presence{$validation_str}){
										#print LOG "\tis a read\n";
										my $q=0;
										my $found=0;
										while($q <= $#check_list && $found==0){
											if($q != $k){
												my @cfr_junctions=@{${$reads{$check_list[$q]}}[1]};
												my $w=0;
												while($w <= $#cfr_junctions && $found == 0){
													my $cfr_junction=$cfr_junctions[$w];
													if($cfr_junction =~ m/$validation_str/){
														$found=1;

													}
													$w++;
												}
											}
											$q++;
										}
										if($found == 0){
											$valid=1;
										}
									}

									$p--;
								}
								

							}
							
							if($valid == 1){
								print LOG "...IS VALID!\n";
								push @valid_junctions, $junction;
							}
							else{
								print LOG "...IS NOT VALID!\n";
							}
						}

						push @valid_junctions_ref, \@valid_junctions;

						print LOG "********************************************\n";
						
						$k++;
					}
					$k=0;
					while($k <= $#check_list){
						${$reads{$check_list[$k]}}[1]=$valid_junctions_ref[$k];
						$k++;
					}

					print LOG "********************************************\n";
				}
			}			
		}

		print LOG "Giunzioni viste da sinistra dopo validazione\n\n";
		print OUT_JUN "Giunzioni viste da sinistra\n\n";
		foreach my $in(keys %reads){
			if(${$reads{$in}}[2] == 1){
				print LOG "Giunzione read numero $in: ", ${$reads{$in}}[0], " shared: ", ${$reads{$in}}[9], "\n";		
				print OUT_JUN ">>>\n";		
				print LOG join("\n", @{${$reads{$in}}[1]}), "\n\n";
				print OUT_JUN join("\n", @{${$reads{$in}}[1]}), "\n\n";
			}
		}
		print LOG "***************************************\n";
		print OUT_JUN "***************************************\n";

		foreach my $repeat(keys %right_check_false){
			my @ids=@{$right_check_false{$repeat}};
			if(scalar(@ids) > 1){
				my $rev=reverse($repeat);
				print LOG "Candidate repeat (from right) of ", length($rev),"bp: $rev\n";
				foreach my $id(@ids){
					die "Error4!\n" if(${$reads{$id}}[5] != 1);
					my @junctions=@{${$reads{$id}}[4]};
					print LOG "Junction (on reverse):\n", join("\n", @junctions), "\n";
				}		
				print LOG "********************************************\n";

				my $i=0;
				my @check_flags=();
				foreach my $id(@ids){
					push @check_flags, 0;
				}
				my $cluster_check=1;
				my %cluster_check_hash=();
				while($i <= $#ids){
					my @junctions1=@{${$reads{$ids[$i]}}[4]};
					if($check_flags[$i] == 0){
						
					 	my $j=$i+1;
					 	my $check=0;
						my $right_shared1=${$reads{$ids[$i]}}[10];	
						my %suffices1=();
						foreach my $j(@junctions1){
							$suffices1{substr($j, $right_shared1, length($j)-$right_shared1)}=1;
						}
						die "Error5!\n" if(scalar(@junctions1) != scalar(keys %suffices1));
						print LOG "SUFF1 ", join("\n", keys %suffices1), "\n";

						my @check_list=($ids[$i]);
					 	while($j <= $#ids){
							my @junctions2=@{${$reads{$ids[$j]}}[4]};
							if(scalar(@junctions1) == scalar(@junctions2)){
								my $right_shared2=${$reads{$ids[$j]}}[10];	
								my %suffices2=();
								foreach my $j(@junctions2){
									$suffices2{substr($j, $right_shared2, length($j)-$right_shared2)}=1;
								}
								die "Error6!\n" if(scalar(@junctions2) != scalar(keys %suffices2));
								print LOG "SUFF2 ", join("\n", keys %suffices2), "\n";

								my $count_equals=0;

								my @keys1=keys %suffices1;
								my @keys2=keys %suffices2;
								my $p=0;
								while($p <= $#keys1){
									my $q=0;
									my $ok=0;
									while($q <= $#keys2 && $ok == 0){
										if($keys1[$p] =~ m/^$keys2[$q]/ || $keys2[$q] =~ m/^$keys1[$p]/){
											$ok=1;
											$count_equals++;
										}
										$q++;			
									}
									$p++;						
								}
								print LOG "Count $count_equals\n";

								if($count_equals == scalar(@junctions2)){
									$check=1;
									push @check_list, $ids[$j];
									$check_flags[$j]=$cluster_check;
								}							
							}							

							$j++;
					 	}
						
						if($check == 1){
							$check_flags[$i]=$cluster_check;
							$cluster_check_hash{$cluster_check}=\@check_list;
							$cluster_check++;
						}
					}
					$i++;
				}
				print LOG "False junctions:\n";
				foreach my $cluster(keys %cluster_check_hash){
					my @check_list=@{$cluster_check_hash{$cluster}};
					print LOG "Cluster to be checked:\n";
					my @shared_regions=();
					foreach my $id(@check_list){
						die "Error6!\n" if(${$reads{$id}}[5] != 1);
						my @junctions=@{${$reads{$id}}[4]};
						print LOG "Junction (on reverse):\n", join("\n", @junctions), "\n";
						my $right_shared=${$reads{$id}}[10];
						push @shared_regions, substr($junctions[0], 0, $right_shared);
					}
					print LOG "Shared:\n", join("\n", @shared_regions), "\n";
					
					my $min_limit=length($shared_regions[0]);
					my $k=1;
					while($k <= $#shared_regions){
						if($min_limit > length($shared_regions[$k])){
							$min_limit=length($shared_regions[$k]);
						}
						$k++;
					}
					$k=$min_repeat_length;
					my $stop=0;
					while($k < $min_limit && $stop == 0){
						my $p=1;
						while($p <= $#shared_regions && $stop == 0){
							if(substr($shared_regions[0], length($shared_regions[0])-$k-1, 1) ne substr($shared_regions[$p], length($shared_regions[$p])-$k-1, 1)){
								$stop=1;
							}
							$p++;
						}					

						$k++;
					}
					my $repeat_length=$k-1;
					print LOG "Lunghezza max della ripetizione: $repeat_length\n";

					$k=0;
					my @valid_junctions_ref=();
					while($k <= $#check_list){
						my @junctions=@{${$reads{$check_list[$k]}}[4]};
						my $right_shared=${$reads{$check_list[$k]}}[10];
						my @valid_junctions=();

						foreach my $junction(@junctions){
							print LOG "Junction (on reverse): $junction\n";

							my $valid=0;

							my $l=length(${$reads{$check_list[$k]}}[0]);
							my $min_range=($l-2 < $right_shared)?($l-2):($right_shared);
							my $max_range=($repeat_length > ($l-length($junction)+$right_shared))?($repeat_length):($l-length($junction)+$right_shared);
							#print LOG "MIN: $min_range MAX: $max_range\n";
							
							if($min_range >= $max_range && $valid == 0){
								my $p=$min_range;
								
								while($p >= $max_range && $valid == 0){
									my $validation_str=substr($junction, $right_shared-$p-1, $l);
									my $rev_validation_str=reverse($validation_str);

									#print LOG "Candidate validation read (on forward): $rev_validation_str\n";
									if(exists $presence{$rev_validation_str}){
										#print LOG "\tis a read\n";
										my $q=0;
										my $found=0;
										while($q <= $#check_list && $found==0){
											if($q != $k){
												my @cfr_junctions=@{${$reads{$check_list[$q]}}[4]};
												my $w=0;
												while($w <= $#cfr_junctions && $found == 0){
													my $cfr_junction=$cfr_junctions[$w];
													if($cfr_junction =~ m/$validation_str/){
														$found=1;

													}
													$w++;
												}
											}
											$q++;
										}
										if($found == 0){
											$valid=1;
										}
									}

									$p--;
								}
								

							}
							
							if($valid == 1){
								print LOG "...IS VALID!\n";
								push @valid_junctions, $junction;
							}
							else{
								print LOG "...IS NOT VALID!\n";
							}
						}

						push @valid_junctions_ref, \@valid_junctions;

						print LOG "********************************************\n";
						
						$k++;
					}
					$k=0;
					while($k <= $#check_list){
						${$reads{$check_list[$k]}}[4]=$valid_junctions_ref[$k];
						$k++;
					}

					print LOG "********************************************\n";
				}
			}			
		}

		print LOG "Giunzioni viste da destra dopo validazione\n\n";
		print OUT_JUN "Giunzioni viste da destra\n\n";
		foreach my $in(keys %reads){
			if(${$reads{$in}}[5] == 1){
				print LOG "Giunzione read numero $in: ", ${$reads{$in}}[0], " shared ", ${$reads{$in}}[10], "\n";
				print OUT_JUN ">>>\n";

				my $max=0;
				foreach my $g(@{${$reads{$in}}[4]}){
					if($max < length($g)){
						$max=length($g);
					}					
				}
				foreach my $g(@{${$reads{$in}}[4]}){
					my $rev=reverse($g);
					print LOG " "x($max-length($g)), $rev, "\n";
					print OUT_JUN " "x($max-length($g)), $rev, "\n";
				}
				print LOG "\n";
				print OUT_JUN "\n";
			}
		}
		print LOG "***************************************\n";
		print OUT_JUN "***************************************\n";

		#ARRICCHIMENTO (SOLO CON GIUNZIONI SINGOLARMENTE PRESE)
		#Arricchisco il set di input
		foreach my $in(keys %reads){
			my $l=length(${$reads{$in}}[0]);
						
			if(${$reads{$in}}[2] == 1){
				#print LOG "I: ", $in, "\n";
				
				my @l_junctions=@{${$reads{$in}}[1]};

				#Per controllo di giunzione...
				my @suffices=();
				my $left_shared=${$reads{$in}}[9];			

				foreach my $junction(@l_junctions){
					#print LOG "Junction raffa: $junction\n";

					push @suffices, substr($junction, $left_shared, length($junction)-$left_shared);
				
					my $i;
					for($i=0; $i < length($junction)-$l; $i++){
						my $add_read=substr($junction, $i, $l);
						
						if(!exists $presence{$add_read}){
							$presence{$add_read}=1;
						}
					}
				}

				#Controllo di giunzione
				my $p=0;
				while($p <= $#suffices){
					my $q=$p+1;
					while($q <= $#suffices){
						if($suffices[$p]=~ m/^$suffices[$q]/ || $suffices[$q]=~ m/^$suffices[$p]/){
							warn "\tError1!\n";
						}
						$q++;
					}
					$p++;
				}
			
			}

			if(${$reads{$in}}[5] == 1){
				#print LOG "I: ", $in, "\n";
				my @r_junctions=@{${$reads{$in}}[4]};

				my @rev_prefices=();
				my $right_shared=${$reads{$in}}[10];

				foreach my $junction(@r_junctions){
					my $rev=reverse($junction);
					#print LOG "Junction raffa: $rev\n";

					push @rev_prefices, substr($junction, $right_shared, length($junction)-$right_shared);
					my $i;
					for($i=0; $i < length($rev)-$l; $i++){
						my $add_read=substr($rev, $i, $l);
						
						if(!exists $presence{$add_read}){
							$presence{$add_read}=1;
						}
					}
				}

				#Controllo di giunzione
				my $p=0;
				while($p <= $#rev_prefices){
					my $q=$p+1;
					while($q <= $#rev_prefices){
						if($rev_prefices[$p]=~ m/^$rev_prefices[$q]/ || $rev_prefices[$q]=~ m/^$rev_prefices[$p]/){
							warn "\tError2!\n";
						}
						$q++;
					}
					$p++;
				}
			}
		}

		my @r_list=keys %presence;
		my $incr=1;
		foreach my $r(@r_list){
			print OUT ">read$incr\n";
			print OUT $r, "\n";
			$incr++;
		}

		#NON METTERE IF(1) PERCHE' NON FUNZIONEREBBE A CAUSA DEL FATTO CHE ${$reads{$ids[$j]}}[8] E' ORA UNA LISTA!!!
		if(0){
		#Il matching tra giunzioni viste da sx e da dx potrebbe servire per evidenziare ripetizioni
		print OUT_JUN "***************************************\n";
		print OUT_JUN "Matching giunzioni\n\n";

		#Per tenere traccia delle giunzioni da sinistra che sono accoppiate con una giunzione da destra (e viceversa) -> key=read_index, value=\%hash, %hash -> key=junction_index (indice relativo alle liste @{${$reads{$in}}[1]} e @{${$reads{$in}}[4]}), value=[match_read_index, match_junction_index, left_before_align, left_after_align, right_before_align, right_after_align]
		my %left_matches=();	#Da giunzioni sinistre a destre
		my %right_matches=();	#Viceversa

		#RICORDARSI CHE LE GIUNZIONI DA DESTRA SONO MEMORIZZATE IN REVERSE!!!

		foreach my $in(keys %reads){
			if(${$reads{$in}}[2] == 1){
				#print LOG "Giunzione sinistra relativa a read numero $in: ", ${$reads{$in}}[0], "\n";
				my $junction_index=0;
				foreach my $seq (@{${$reads{$in}}[1]}){
					#print LOG "\tSEQ: $seq\n";
					my $junction_length=length($seq);
					my $read_length=length(${$reads{$in}}[0]);
					#print LOG "\tL: $junction_length l $read_length\n";

					my $i=$junction_length-1;
					my $stop=0;
					while($i >= $read_length-1 && $stop == 0){
						my $hash_seq=substr($seq, $i-$read_length/2+1, $read_length/2);
						#print LOG "hash $hash_seq\n";

						my @ids=();
						push @ids, @{$left_hash{$hash_seq}} if(exists $left_hash{$hash_seq});
						push @ids, @{$right_hash{$hash_seq}} if(exists $right_hash{$hash_seq});

						my $j=0;
						while($j <= $#ids && $stop == 0){
						   if(${$reads{$ids[$j]}}[5] == 1 || ${reads{$ids[$j]}}[6] == 1){
							#print LOG "Read ", ${$reads{$ids[$j]}}[0], "\n";

							#DA CORREGGERE PERCHE' ${$reads{$ids[$j]}}[8] E' UNA LISTA!!!!
							my $ref_in=${$reads{$ids[$j]}}[8];
							#Il ciclo dovrebbe risolvere il fatto che questo campo viene aggiornato solo per il read che termina una giunzione destra e che viene inclusa in una giunzione che si estende verso destra (non viene aggiornato per tutti gli altri che mantengono l'originale)
							while(${$reads{$ref_in}}[5] == 0){
								$ref_in=${$reads{$ref_in}}[8];
							}

							my @match_junctions=@{${$reads{$ref_in}}[4]};
							my $match_junction_index=0;
							while($match_junction_index <= $#match_junctions && $stop == 0){
								my $match_j=$match_junctions[$match_junction_index];
								my $rev=reverse($match_j);
								#print LOG "Junction: ", $rev, "\n";

								$seq =~ m/$hash_seq/;
								my $prefix_seq=$`;
								my $suffix_seq=$';
								while($rev =~ m/$hash_seq/g && $stop == 0){
									 
								 my $prefix_rev=$`;
								 my $suffix_rev=$';
								 #print LOG "Before: $prefix_seq\n";
								 #print LOG "Before: $prefix_rev\n";
								 #print LOG "After: $suffix_seq\n";
								 #print LOG "After: $suffix_rev\n";

								 my $min1=(length($prefix_seq) < length($prefix_rev))?(length($prefix_seq)):(length($prefix_rev));
								 my $min2=(length($suffix_seq) < length($suffix_rev))?(length($suffix_seq)):(length($suffix_rev));
								 my $prefix_seq_check=substr($prefix_seq, length($prefix_seq)-$min1, $min1);
								 my $prefix_rev_check=substr($prefix_rev, length($prefix_rev)-$min1, $min1);
								 my $suffix_seq_check=substr($suffix_seq, 0, $min2);
								 my $suffix_rev_check=substr($suffix_rev, 0, $min2);

								 #print LOG $prefix_seq_check, " VS ", $prefix_rev_check, "\n";
								 #print LOG $suffix_seq_check, " VS ", $suffix_rev_check, "\n";
								 if($prefix_seq_check eq $prefix_rev_check && $suffix_seq_check eq $suffix_rev_check){
									$stop=1;
									#print LOG "($in,$junction_index) ($ref_in,$match_junction_index)\n";
									if(exists $left_matches{$in}){
										${$left_matches{$in}}{$junction_index}=[$ref_in, $match_junction_index, length($prefix_seq)-$min1, length($suffix_seq)-$min2, length($prefix_rev)-$min1, length($suffix_rev)-$min2];
									}else{
										my %temp_hash=($junction_index, [$ref_in, $match_junction_index, length($prefix_seq)-$min1, length($suffix_seq)-$min2, length($prefix_rev)-$min1, length($suffix_rev)-$min2]);
										$left_matches{$in}=\%temp_hash;
#{$junction_index, [$ref_in, $match_junction_index, length($prefix_seq)-$min1, length($suffix_seq)-$min2, length($prefix_rev)-$min1, length($suffix_rev)-$min2]};
									}

									if(exists $right_matches{$ref_in}){
										${$right_matches{$ref_in}}{$match_junction_index}=[$in, $junction_index, length($prefix_seq)-$min1, length($suffix_seq)-$min2, length($prefix_rev)-$min1, length($suffix_rev)-$min2];
									}else{
										$right_matches{$ref_in}={$match_junction_index, [$in, $junction_index, length($prefix_seq)-$min1, length($suffix_seq)-$min2, length($prefix_rev)-$min1, length($suffix_rev)-$min2]};
									}
								 }
								}
									
								$match_junction_index++;
							}								
						   }
						   $j++;
						}
			
						$i--;
					}
					$junction_index++;

				}				
				
				#print LOG "Shared prefix of ", ${$reads{$in}}[9], "bp\n";	
				#print LOG "***************************************\n";			
			}
		}

		foreach my $read_index(keys %left_matches){
			foreach my $junction_index(keys %{$left_matches{$read_index}}){
				print OUT_JUN "Giunzione sx: ", ${${$reads{$read_index}}[1]}[$junction_index], "\n";
				print OUT_JUN "\ttaglio: ", ${$reads{$read_index}}[9], "\n";
				my @list=@{${$left_matches{$read_index}}{$junction_index}};
				my $rev=reverse(${${$reads{$list[0]}}[4]}[$list[1]]);
				print OUT_JUN "\tmatch to giunzione dx: ", $rev, "\n";
				print OUT_JUN "\t\ttaglio: ", ${$reads{$list[0]}}[10], "\n";
				print OUT_JUN "Before common region on sx: ", $list[2], "\n";
				print OUT_JUN "After common region on sx: ", $list[3], "\n";
				print OUT_JUN "Before common region on dx: ", $list[4], "\n";
				print OUT_JUN "After common region on dx: ", $list[5], "\n";				
			}
		}

		foreach my $read_index(keys %right_matches){
			foreach my $junction_index(keys %{$right_matches{$read_index}}){
				my $rev=reverse(${${$reads{$read_index}}[4]}[$junction_index]);
				print OUT_JUN "Giunzione dx: ", $rev, "\n";
				print OUT_JUN "\ttaglio: ", ${$reads{$read_index}}[10], "\n";
				my @list=@{${$right_matches{$read_index}}{$junction_index}};
				print OUT_JUN "\tmatch to giunzione sx: ", ${${$reads{$list[0]}}[1]}[$list[1]], "\n";
				print OUT_JUN "\t\ttaglio: ", ${$reads{$list[0]}}[9], "\n";
				print OUT_JUN "Before common region on sx: ", $list[2], "\n";
				print OUT_JUN "After common region on sx: ", $list[3], "\n";
				print OUT_JUN "Before common region on dx: ", $list[4], "\n";
				print OUT_JUN "After common region on dx: ", $list[5], "\n";				
			}
		}
		}#end if(0)

		close OUT;
		close OUT_JUN;
		close LOG;
	}

}

exit;

sub search_cut{
	my $ref_reads=shift;
	my $start=shift;
	local *FH=shift;

	my @reads=@{$ref_reads};

	print FH "LAST ", join("\n", @reads), "\n";

	#my $i=length($reads[0])/2;
	my $i=$start;
	my $stop=0;
	while($i < length($reads[0]) && $stop == 0){
		my $j=1;
		my $col_char=substr($reads[0], $i, 1);
		while($j <= $#reads && $stop == 0){
			if(substr($reads[$j], $i, 1) ne $col_char && substr($reads[$j], $i, 1) ne "*"){
				$stop=1;
			}
			$j++;
		}
		$i++;
	}

	if($stop == 1){
		$i--;
		my $shared_prefix=substr($reads[0], 0, $i);
		my @junction_suffices=();
		my $j=$#reads;
		while($j >= 0){
			my $s=substr($reads[$j], $i, length($reads[$j])-$i);
			my $k=0;
			my $stop=0;
			while($k <= $#junction_suffices && $stop == 0){
				#my $pref=substr($junction_suffices[$k], 0, length($s));
				print FH $junction_suffices[$k], "\n";
				if($s =~ m/^$junction_suffices[$k]/){
					$stop=1;
					$junction_suffices[$k]=$s;
				}else{
					if($junction_suffices[$k] =~ m/^$s/){
						$stop=1;
					}
				}

				#if($pref eq $s){
				#	$stop=1;
				#}

				$k++;
			}
			if($stop == 0){
				push @junction_suffices, $s;
			}
			
			$j--;
		}
		$j=0;
		while($j <= $#junction_suffices){
			$junction_suffices[$j]=$shared_prefix.$junction_suffices[$j];
			$j++;
		}

		return [\@junction_suffices, length($shared_prefix)];
	}
	else{
		return [[], -1];
	}
}

sub reverse{
	my $sequence=shift;

	my $reverse="";
	
	my $i=0;
	while($i < length($sequence)){
		$reverse=substr($sequence, $i, 1).$reverse;
		$i++;
	}

	return $reverse;
}


